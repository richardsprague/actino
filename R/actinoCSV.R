# actinoCSV.R
# CSV to Phyloseq functions

# Requirements:
#
# Experiment dataframe with the following exact column names (additional columns are okay but will be ignored):
#   "ssr": sequencing revision (e.g. 42578).  Your column probably repeats many SSRs.
#   "tax_name": (e.g. "Bifidobacterium").
#   "count" : reads for that taxa at that SSR (e.g. 2137)
#
# Mapping file: a dataframe that contains:
#   "ssr": exact same SSRs as experiment dataframe above
#   attributes columns: as many as you like.  (e.g. "yogurteater", "geo", "kefireater", "gender")

#' @title Turn a dataframe into a format that Phyloseq likes
#' @description Convenience function to make a dataframe where sample names are in columns, not rows
#' @param DF data frame with columns: SSR, tax_name, count, (optional) tax_rank
#' @importFrom tidyr %>% spread
#' @importFrom dplyr select count
#' @export
samples_as_column_dataframe <- function (DF) {
  # Convenience function to make a dataframe where sample names are in columns, not rows
  #
  # Args:
  #   DF: a dataframe with the following columns:
  #     SSR
  #     tax_name
  #     count
  #     (optional): tax_rank
  #
  # Returns:
  #   A dataframe where
  #     first col: tax_name
  #     other cols : SSRs (aka sample names)
  #     rows : taxa names (e.g. Bifidobacterium, etc.)

  cols <- names(DF)
  SSR <- match("ssr",cols)
  t = DF %>% select(tax_name, SSR, count) %>% spread(ssr,count)
  t[is.na(t)]<-0 # replace NA with 0
  return(t)
}


#' @title convert a uBiome experiment dataframe and mapping file to a valid phyloseq object
#' @description Make a valid Phyloseq object out of a uBiome dataframe and a map file
#' @param experiment.df dataframe representation of a uBiome experiment. Rows are unique ssrs.
#' @param mapfile a dataframe with column 'ssr'
#' @param rank (not used currently -- just assumes everything is always just genus)
#' @return a Phyloseq object
#' @importFrom phyloseq phyloseq parse_taxonomy_qiime build_tax_table otu_table sample_data tax_table
#' @examples simple.df <- data.frame(tax_name=c("taxa1"),"1234"=c(2352))
#' simple.map <- data.frame(ssr=c(1234))
#' @export
experiment_to_phyloseq <- function(experiment.df, mapfile, rank="genus"){
  # convert a uBiome experiment dataframe and mapping file to a valid phyloseq object
  #
  # Args:
  #   experiment.df:
  #   mapfile:
  #   rank: (not used currently -- just assumes everything is always just genus)
  #
  e <- experiment.df %>% filter(tax_rank == rank) %>% samples_as_column_dataframe()
  ssrs <- unique(experiment.df$ssr)
  e.map <- mapfile[match(ssrs,mapfile$ssr),]
  row.names(e.map)<-e.map$ssr

  qiime_tax_names<-sapply(e[1],function (x) paste("g__",x,sep=""))  # assumes every taxa is genus level
  e.taxtable<-build_tax_table(lapply(qiime_tax_names,parse_taxonomy_qiime))
  dimnames(e.taxtable)[[1]]<-qiime_tax_names


  e.matrix<-as.matrix(e[,2:ncol(e)]) # eliminate final col, which is a repeat of the taxonomy info
  colnames(e.matrix)<-ssrs
  rownames(e.matrix)<-qiime_tax_names


  e.ps <-phyloseq(otu_table(e.matrix,taxa_are_rows=TRUE),
                  sample_data(e.map),
                  tax_table(e.taxtable))
  return(e.ps)
}


