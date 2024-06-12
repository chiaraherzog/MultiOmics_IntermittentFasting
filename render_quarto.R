quarto::quarto_render('index.qmd')
input <- list.files(".", recursive = F)[grepl(".qmd$", list.files(".", recursive = F))]

for (i in input[!input %in% c('index.qmd', 'helper-variance-decomposition.qmd')]){
  quarto::quarto_render(i)
}


quarto::quarto_render("ImmuneAge_all.qmd")

