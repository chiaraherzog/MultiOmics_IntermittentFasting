# Export populations
library(here)
load(here("src/populations_names_annotated.Rdata"))
populations_exp <- populations |> 
  dplyr::mutate(variable = name,
                name = `second name`,
                staining = `fixable viability dye and antibodies α-`) |> 
  dplyr::select(type, name, staining, targets, antibodies, variable, main_analysis)

writexl::write_xlsx(populations_exp,  path = "out/ImmuneStaining.xlsx")
