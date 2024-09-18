#' @title Process and Summarise Data data for T1
#' @description This function processes raw data to filter, transform, and summarize it for table generation.
#' It also adds medication flags, computes baseline values, and calculates changes in various features.
#' Finally, it generates a summary table and exports it to HTML and PDF.
#' @param data_raw_path The file path to the raw data file (e.g., "data/data_raw.Rdata").
#' @param baseline_change_path The file path to the baseline change data file (e.g., "data/data_baseline_change.Rdata").
#' @param output_html The file path for the output HTML table (e.g., "out/tab1.html").
#' @param output_pdf The file path for the output PDF table (e.g., "out/tab1.pdf").
#' @return A gt table object.
#' @export
process_and_summarise_data_T1 <- function(data_raw_path, baseline_change_path, output_html, output_pdf) {
  
  # Load raw data
  cat("Loading raw data...\n")
  load(data_raw_path)
  
  # Convert colData to data frame and filter rows
  cat("Filtering and transforming raw data...\n")
  df <- as.data.frame(colData(data)) |> 
    dplyr::filter(visitId == "M0" & interventionId != "S") |> 
    dplyr::mutate(
      complabel = ifelse(is.na(dropout_date), paste0(compliance, "<br>compliance"), "dropout<br>"),
      complabel = factor(complabel, levels = c("high<br>compliance", "medium<br>compliance", "low<br>compliance", "dropout<br>")),
      
      smoking_history = case_when(
        smoking_ever == "yes" ~ "Former smoker",
        smoking_ever == "no" ~ "Never smoker",
        smoking_ever == "unknown" ~ "unknown"
      ),
      smoking_history = factor(smoking_history, levels = c("Never smoker", "Former smoker", "unknown")),
      
      diet = case_when(
        diet == 'normal' ~ 'normal',
        diet != 'unknown' ~ 'vegetarian/vegan/pescetarian',
        diet == 'unknown' ~ 'unknown'
      ),
      diet = factor(diet, levels = c("normal", 'vegetarian/vegan/pescetarian', 'unknown')),
      
      activity = case_when(
        is.na(intactcurr) ~ "unknown",
        intactcurr %in% c("0") ~ "no intense activity",
        intactcurr %in% c("30 min", "1 h", "1.5 h", "2 h") ~ "30-149 min",
        TRUE ~ "≥150 min"
      ),
      activity = factor(activity, levels = c("no intense activity", "30-149 min", "≥150 min", "unknown"))
    ) |> 
    dplyr::select(subjectId, age_at_consent, bmi_at_consent, mpstatrs, complabel, diet, smoking_history, activity, medtype, med_cur)
  
  # Process medications
  cat("Processing medication data...\n")
  df_sep <- df |> 
    dplyr::select(subjectId, medtype, med_cur) |> 
    tidyr::separate_rows(medtype, sep = ";") |> 
    dplyr::mutate(medtype = case_when(
      medtype %in% c("antidepressant", "antihypertensive") ~ medtype,
      is.na(medtype) ~ NA_character_,
      TRUE ~ "other"
    ))
  
  medications_with_flags <- df_sep |> 
    dplyr::mutate(
      antidepressant = ifelse(med_cur == 'yes' & medtype == "antidepressant", 1, 0),
      antihypertensive = ifelse(med_cur == 'yes' & medtype == "antihypertensive", 1, 0),
      other = ifelse(med_cur == 'yes' & !(medtype %in% c("antidepressant", "antihypertensive")), 1, 0)
    )
  
  med_summary <- medications_with_flags |> 
    dplyr::group_by(subjectId) |> 
    dplyr::summarise(
      antidepressant = coalesce(max(antidepressant), 0),
      antihypertensive = coalesce(max(antihypertensive), 0),
      other = coalesce(max(other), 0)
    ) |> 
    dplyr::mutate(across(c("antidepressant", "antihypertensive", "other"), ~ ifelse(. == 1, "yes", "no"))) |> 
    dplyr::ungroup()
  
  df <- df |> 
    dplyr::left_join(med_summary)
  
  # Compute baseline values
  cat("Computing baseline values...\n")
  baseline <- as.data.frame(wideFormat(data[c("sysbp", "diabp", "cholesterol", "triglycerides", "hba1c", "glucose", "hemoglobin"), , c('Functional sports exam', 'Blood haemogram')], colDataCols = c("visitId", "subjectId"))) |> 
    dplyr::filter(visitId == "M0") 
  
  df <- df |> 
    dplyr::left_join(baseline, by = 'subjectId') |> 
    dplyr::select(-c(primary, visitId)) |> 
    dplyr::mutate(Functional.sports.exam_diabp = as.numeric(Functional.sports.exam_diabp)) 
  
  # Load baseline change data and compute changes
  cat("Loading baseline change data and computing changes...\n")
  load(baseline_change_path)
  
  change <- as.data.frame(wideFormat(data[c("sysbp", "diabp", "cholesterol", "triglycerides", "hba1c", "glucose", "hemoglobin", "bmi"), , c('Functional sports exam', 'Blood haemogram', 'Body composition')], colDataCols = c("visitId", "subjectId"))) |> 
    dplyr::filter(visitId == "M6") |> 
    dplyr::rename_with(~ paste0("delta_", .), -c("primary", "visitId", "subjectId"))
  
  df <- df |> 
    dplyr::left_join(change, by = 'subjectId') |> 
    dplyr::select(-c(primary, visitId, medtype)) |> 
    dplyr::mutate(delta_Functional.sports.exam_diabp = as.numeric(delta_Functional.sports.exam_diabp)) 
  
  # Process BMI changes for M2 and M4
  cat("Processing BMI changes for M2 and M4...\n")
  change <- as.data.frame(wideFormat(data["bmi", , 'Body composition'], colDataCols = c("visitId", "subjectId"))) |> 
    dplyr::filter(visitId %in% c("M2", "M4")) |> 
    tidyr::pivot_wider(id_cols = c("subjectId"), names_from = visitId, values_from = Body.composition_bmi) |> 
    dplyr::rename(M2_bmi = M2, M4_bmi = M4)
  
  df <- df |> 
    dplyr::left_join(change, by = 'subjectId') |> 
    dplyr::select(-c(subjectId))
  
  # Generate summary table
  cat("Generating summary table...\n")
  labels <- list(
    age_at_consent = '<b>Age at consent</b>',
    bmi_at_consent = gt::html('<b>BMI at consent</b> (kg/m²)'),
    M2_bmi = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;∆ M2'),
    M4_bmi = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;∆ M4'),
    delta_Body.composition_bmi = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;∆ M6'),
    smoking_history = gt::html('<b>Smoking history</b>'),
    diet = gt::html('<b>Dietary preference</b>'),
    activity = gt::html('<b>Weekly intense activity (min)</b>'),
    mpstatrs = gt::html('<b>Postmenopausal</b>'),
    med_cur = gt::html('<b>Regular medication use</b>'),
    antihypertensive = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;Antihypertensive medication'),
    antidepressant = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;Antidepressant medication'),
    other = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;Other medication'),
    Functional.sports.exam_sysbp = gt::html('<b>Systolic blood pressure</b> (mmHg)'),
    delta_Functional.sports.exam_sysbp = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;∆ M6'),
    Functional.sports.exam_diabp = "<b>Diastolic blood pressure</b> (mmHg)",
    delta_Functional.sports.exam_diabp = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;∆ M6'),
    Functional.sports.exam_cholesterol = gt::html('<b>Total cholesterol</b> (mg/dL)'),
    delta_Functional.sports.exam_cholesterol = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;∆ M6'),
    Functional.sports.exam_triglycerides = gt::html('<b>Triglycerides</b> (mg/dL)'),
    delta_Functional.sports.exam_triglycerides = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;∆ M6'),
    Functional.sports.exam_hba1c = gt::html('<b>HbA1c</b> (%)'),
    delta_Functional.sports.exam_hba1c = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;∆ M6'),
    Functional.sports.exam_glucose = gt::html('<b>Fasting glucose</b> (mg/dL)'),
    delta_Functional.sports.exam_glucose = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;∆ M6'),
    Blood.haemogram_hemoglobin = gt::html('<b>Haemoglobin</b> (g/dL)'),
    delta_Blood.haemogram_hemoglobin = gt::html('&nbsp;&nbsp;&nbsp;&nbsp;∆ M6')
  )
  
  t1 <- df |> 
    dplyr::select(complabel, dplyr::any_of(names(labels))) |> 
    gtsummary::tbl_summary(
      by = 'complabel',
      type = list(Functional.sports.exam_diabp ~ "continuous", delta_Functional.sports.exam_diabp ~ "continuous"),
      label = labels,
      digits = gtsummary::all_continuous() ~ 1,
      missing = "no",
      statistic = list(
        gtsummary::all_continuous() ~ "{mean} ({sd})",
        Functional.sports.exam_triglycerides ~ "{median} ({IQR})",
        delta_Functional.sports.exam_triglycerides ~ "{median} ({IQR})"
      )
    ) |> 
    gtsummary::add_overall() |> 
    gtsummary::modify_header(update = gtsummary::all_stat_cols() ~ "<b>{level}</b><br>n = {n}") |> 
    gtsummary::modify_header(update = stat_0 ~ "<b>{level}</b><br><br>n = {n}")
  
  # Convert to gt and customize the table appearance
  cat("Customizing table appearance...\n")
  t2 <- t1 |> 
    gtsummary::as_gt() |> 
    gt::fmt_markdown(columns = c("label")) |>
    gt::text_transform(
      locations = gt::cells_body(),
      fn = function(x) {
        gsub("NA \\(NA\\)", "--", x)
      }
    ) |> 
    gt::tab_options(
      table.width = gt::px(680),
      table.font.size = 12,
      column_labels.font.size = 13,
      table.font.names = "Helvetica"
    ) |> 
    gt::opt_stylize(style = 3) |> 
    gtExtras::gt_add_divider(columns = "stat_0", sides = "right", style = 'solid', weight = gt::px(1)) |> 
    gt::cols_align(align = 'right', columns = dplyr::contains("stat_"))
  
  # Save the table to HTML and PDF
  cat("Saving the table to HTML and PDF...\n")
  gt::gtsave(t2, file = output_html)
  pagedown::chrome_print(output_html, output = output_pdf)
  
  cat("Process completed successfully!\n")
  
  return(list(table = t2,
              data = df))
}
