# convert matrices into numeric
convert_numeric_matrix <- function(x){
  
  out <- vector("list", length(x))
  
  out <- sapply(1:length(x), function(i){
   
   tmp <- x[[i]] |> 
     dplyr::mutate(across(c(M0:M6), ~ as.numeric(stringr::str_trim(stringr::str_split(., "\\(", simplify = T)[,1])))) |> 
     dplyr::mutate(across(c("M2", "M4", "M6"), ~ (./M0)*100))
 
   out[[i]] <- tmp
   
   })
  
}
