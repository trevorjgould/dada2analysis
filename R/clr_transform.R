#' clr + imputation function
#'
#' This function runs clr transform
#' 
#' @param x table
#'
#' @export
#' 
## Borrowed from Gabe
## Note, this assumes taxa as rows and samples as columns
clr_transform <- function(x){
  clr.taxa <- x
  clr.taxa = t(clr.taxa); eps = 0.5
  clr.taxa = clr.taxa*(1 -rowSums(clr.taxa==0)*eps/rowSums(clr.taxa))
  clr.taxa[clr.taxa==0]=eps
  clr.taxa = sweep(clr.taxa,1,rowSums(clr.taxa),'/');
  ls = log(clr.taxa)
  clr.taxa = t(ls - rowMeans(ls))
  clr.taxa = clr.taxa[,!is.nan(colSums(clr.taxa))]
  return(clr.taxa)
}
