load("out/out_lmm_factor.Rdata")
load("out/outPerm.R")

observed <- out_lmm$`Basic model with BMI`
random <- outPerm
variable = observed$x[grepl("bmi$", observed$x)]



source("src/getPermutationComparisonPval.R")


corr <- getPermutationComparisonPval(observed, random, variables = unique(observed$x)[1:300])
