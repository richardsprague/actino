# actinoUtils.R

#' @title Dataframe at a specific rank
#' @description turn a uBIome data frame into a matrix at a specific rank
#' @param df a well-formed dataframe created by join_all_ubiome_files_full
#' @param rank tax_rank (default = "genus")
#' @importFrom stats na.omit
#' @export
dataframe_at_rank <- function(df, rank="genus" ){
  z <- lapply(df[,1],function(x) {
    l = tax_rank_of_full_taxa(x)
    return(c(l[[1]],l[[2]]))
  })
  names(z) = NULL
  z.df = t(as.data.frame(z))
  row.names(z.df) = NULL
  new.df <- na.omit(data.frame(z.df, df[,-1]))
  colnames(new.df)[1:2] = c("tax_rank","tax_name")
  return(new.df[new.df$tax_rank==rank,][-1])

}

#' @title Matrix version of a dataframe at a specific rank
#' @description turn a uBIome data frame into a matrix at a specific rank
#' @param df a well-formed dataframe created by join_all_ubiome_files_full
#' @param rank tax_rank (default = "genus")
#' @export
matrix_at_rank <- function(df, rank = "genus"){
  mf <- dataframe_at_rank(df,rank)
  m <- as.matrix(mf[,-1])
  rownames(m) <- mf[,1]
  m
}
