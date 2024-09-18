quarto::quarto_render('index.qmd')
input <- list.files(".", recursive = F)[grepl(".qmd$", list.files(".", recursive = F))]

for (i in input[!input %in% c('index.qmd')]){
  quarto::quarto_render(i)
}
