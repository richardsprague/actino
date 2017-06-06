# actinoJSON.R
# handy functions for reading JSON files

# read uBiome JSON files and convert to tables that can be used by BIOM

# converts every json file in directory to a CSV file

# you can convert a single JSON file to CSV like this:
# convert_json_files_to_csv(pattern="sprague-nose-1410.json")


## Dependencies
# the following packages must be loaded before using this function
#install.packages("jsonlite")

MAPFILE_ATTRIBUTES = c("SSR", "Username", "Site", "Date", "Geo", "Label", "Notes")

#' @title returns the abbreviated character form of the character vector "tax_rank"
#' @description Needed by functions that convert tax_names to QIIME (or other) format.
#' @param tax_rank a character vector (e.g. "genus")
#' @param inverse if true, return the full tax_rank when given the abbreviated version
#' @export
tax_abbrev<-function(tax_rank, inverse=FALSE){
  # returns the abbreviated character form of the character vector "tax_rank"
  # Needed by functions that convert tax_names to QIIME (or other) format.
  # Args:
  #   tax_rank: a character vector (e.g. "genus")
  # Returns:
  #   A one character abbreviation of tax_rank (e.g. "g")
  taxes<-unlist(strsplit("root superkingdom genus species order family phylum class superphylum subclass suborder no_rank species_group subgenus subphylum subkingdom", split=" "))
  tax_abbrev<-unlist(strsplit("r k g s o f p c l 1 r n h i 3 2", split=" "))
  tax<-data.frame(taxes,tax_abbrev)
  ifelse(inverse!=FALSE,
          as.character(tax[tax$tax_abbrev==tax_rank,1]),
          as.character(tax[tax$taxes==tax_rank,2]))
}

full_taxa<- function(df,taxaRow){
  # taxaRow: list
  # # given a row taxaRow in df, return a character vector showing the full set of tax_ranks
  # taxaRow = as.data.frame(t(taxaRow))
  # parent = as.numeric(taxaRow[match("parent",names(df))])#taxaRow[2])
  # taxName = taxaRow[,match("tax_name",names(df))] #taxaRow[5]
  # taxRank = taxaRow[,match("tax_rank",names(df)),] # taxaRow[6]
  parent = as.numeric(taxaRow["parent"])
  taxName = as.character(unlist(taxaRow["tax_name"]))
  taxRank = as.character(unlist(taxaRow["tax_rank"]))



  if(is.na(taxName)){browser(text="Null Taxname")} # some kind of error occurred.
  if(taxName == "root")
    return("Root")
  else
    return(
      paste(full_taxa(df,as.list(df[df$taxon == parent,])),#t(df[df$taxon==parent,])),
            ";",
            tax_abbrev(as.character(taxRank)),"__",
            as.character(taxName),
            sep=""
      )
    )

}

#' @title find tax_rank for a given full taxa
#' @description opposite function for full_taxa(). Return a list [tax_name, tax_rank]
#' @param parse_taxon string representation for qiime full taxonomy for this taxon
#' @export
tax_rank_of_full_taxa <- function(parse_taxon){
# typical input: "Root;n__cellular organisms;k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Alcaligenaceae"
# correct output: list(family,"Alcaligenaceae")

  s = strsplit(parse_taxon,";")
  t = s[[1]][length(s[[1]])]  # last element of s
  if(t=="Root")
    return(list("Root","Root"))
  u = strsplit(t,"__")
  return(list(tax_abbrev(u[[1]][1],inverse=TRUE),u[[1]][2]))

  if(parse_taxon == "Root;n__cellular organisms;k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Burkholderiales;f__Alcaligenaceae")
    return(list("family","Alcaligenaceae"))
   parse_taxon
}

#' @title Convert to CSV all uBiome JSON files in a directory.
#' @description Convenient way to automatically convert all JSON files in a directory to CSV.
#' @param pattern a regular expression specifying which files to look for (default is anything JSON)
#' @param directory the pathname of the directory to search (default is current working directory)
#' @importFrom jsonlite fromJSON
#' @export
convert_json_files_to_csv <- function(pattern= "[[:alnum:]].json", directory=getwd()){
  json_files_in_directory <- list.files(directory, full.names =TRUE,pattern=pattern)
  for (i in json_files_in_directory){
    json_version <- fromJSON(i)  # returns as a list, so we must convert to a data frame
    asFrame <- json_version$ubiome_bacteriacounts # do.call("rbind",lapply(json_version$ubiome_bacteriacounts,as.data.frame))
    fname <-strsplit(i,split=".json")
    csvName <- paste0(fname[[1]],".csv")
    cat("converting...",i," to ",csvName,"\n")
    utils::write.csv(asFrame,file=csvName)
  }
}

#' @title Just JSON files in a directory
#' @description returns a vector of just the json files in a directory d
#' @param d pathname to a directory
#' @export
just_json_files_in<-function(d){
  #json_files_in_directory <-
  list.files(d, full.names = TRUE, pattern="[[:alnum:]].json") #todo: prevent tilde at end of line.
}



#' @title  returns a big dataframe joining all tax_names from a list of JSON files
#' @description returns a big dataframe joining all tax_names from a list of JSON files
#' @param flist list of file pathnames, each pointing to a uBiome JSON file
#' @param tax_rank taxonomic rank to read. Default is "genus"
#' @param count.normalized use the count_norm field instead of raw read numbers (default is FALSE)
#' @return data frame where first col is all taxa names, other cols are samples
#' @importFrom dplyr full_join
#' @export
join_all_ubiome_files <- function(flist,tax_rank="genus",count.normalized=FALSE){
  f.all <- read_ubiome_json(flist[1],rank=tax_rank,count.normalized = count.normalized)
  #f.all <- rbind(f.all, read_ubiome_json(flist[1], rank = "root", count.normalized = count.normalized))
  for(f in flist[-1]){
    f.part<-read_ubiome_json(f,rank=tax_rank,count.normalized)
    # test: # f.part<-read_ubiome_json_full(f,count.normalized)
    if(nrow(f.part)>0 & nrow(f.all)>0)  f.all<-full_join(f.all,f.part)
    else if (nrow(f.part)>0 & nrow(f.all)==0) f.all<-f.part
  }
  f.all
}

#' @title Read a uBiome JSON file into a dataframe
#' @description Reads JSON into a dataframe.
#' @param fname file name as character vector
#' @param count.normalized use the count_norm field instead of raw read numbers (default is FALSE)
#' @importFrom dplyr filter
#' @importFrom jsonlite fromJSON
#' @return data frame representation of the sample
read_ubiome_json_full<-function(fname,count.normalized=FALSE){


  # convenience wrapper that returns a character name based on sample date and the notes field
  # expects dateVal to be a date of type character
  name_for_sample <- function(dateVal,notesVal,ssr){
    paste(as.character(as.Date(dateVal)),"$",ssr,"$",substring(notesVal,6,10))
  }
  j<-jsonlite::fromJSON(fname)

  rj<-j[["ubiome_bacteriacounts"]]
  #p<-apply(rj,1,function (x) full_taxa(rj,as.data.frame(t(x))))
  p<-apply(rj,1,function (x) full_taxa(rj,as.list(x)))
  rj$tax_name<-p  # new column stores the whole tax path for each organism
  r<-data.frame(tax_name=rj$tax_name,
                reads=if(count.normalized) rj$count_norm else rj$count,
                stringsAsFactors = FALSE)

  names(r)<-c("tax_name",name_for_sample(j$sampling_time,j$notes,j$sequencing_revision))

  return(r)

}

#' @title Read a uBiome JSON file into a dataframe, returning only the taxa from a single rank.
#' @description Reads JSON into a dataframe.
#' @param fname file name as character vector
#' @param rank tax_rank as character vector
#' @param count.normalized use the count_norm field instead of raw read numbers (default is FALSE)
#' @importFrom dplyr filter
#' @importFrom jsonlite fromJSON
#' @return data frame representation of the sample
#returns a dataframe of just the genus and the count (or count_norm if count.normalized is TRUE)
read_ubiome_json<-function(fname,rank="genus",count.normalized=FALSE){

  # convenience wrapper that returns a character name based on sample date and the notes field
  # expects dateVal to be a date of type character
  name_for_sample <- function(dateVal,notesVal,ssr){
    paste(as.character(as.Date(dateVal)),"$",ssr,"$",substring(notesVal,6,10))
  }
  j<-jsonlite::fromJSON(fname)


  rj<-filter(j[["ubiome_bacteriacounts"]],tax_rank == rank)
  #p<-apply(rj,1,function (x) full_taxa(rj,as.list(x)))

  r<- data.frame(rj$tax_name,rj$tax_rank,reads=if(count.normalized) rj$count_norm else rj$count,stringsAsFactors = FALSE)
 # r<- data.frame(p,rj$tax_rank,reads=if(count.normalized) rj$count_norm else rj$count,stringsAsFactors = FALSE)

  names(r)<-c("tax_name","tax_rank",name_for_sample(j$sampling_time,j$notes,j$sequencing_revision))
  return(r)

}


#' @title  returns a big dataframe joining all tax_names from a list of JSON files
#' @description returns a big dataframe joining all tax_names from a list of JSON files
#' @param flist list of file pathnames, each pointing to a uBiome JSON file
#' @param count.normalized use the count_norm field instead of raw read numbers (default is FALSE)
#' @param site character vector
#' @return data frame where first col is all taxa names, other cols are samples
#' @importFrom dplyr full_join
#' @export
join_all_ubiome_files_full <- function(flist,count.normalized=FALSE,site="gut"){
  f.all<-read_ubiome_json_full(flist[1],count.normalized )
  for(f in flist[-1]){
    f.part<-read_ubiome_json_full(f,count.normalized)
    if(nrow(f.part)>0 & nrow(f.all)>0)  f.all<-full_join(f.all,f.part)
    else if (nrow(f.part)>0 & nrow(f.all)==0) f.all<-f.part
  }
  f.all
}


#' return a valid Phyloseq object created from the JSON files in flist and an Excel formatted mapfile
#'
#' @param flist character vector of file names
#' @param mapfile XLSX filename that contains mapping info for the same SSRS in flist
#' @param rank taxonomical rank (generally "genus" or "family" or "species")
#' @param count.normalized use the count_norm field instead of raw read numbers (default is FALSE)
#' @return valid Phyloseq object
#' @importFrom readxl read_excel
#' @importFrom phyloseq taxa_names sample_names build_tax_table parse_taxonomy_default tax_table
#' @export
phyloseq_from_JSON_at_rank <- function(flist, mapfile, rank="genus", count.normalized = FALSE) {
  # return a valid Phyloseq object created from the JSON files in flist and an Excel formatted mapfile
  # Args:
  #   flist: character vector of file names.
  #   mapfile: XLSX filename that contains mapping info for the same SSRS in flist
  # Returns:  valid Phyloseq object

  f.all <- join_all_ubiome_files(flist,tax_rank = rank, count.normalized)
  f.all[is.na(f.all)] <- 0
  f.mat <- f.all[, c(-1,-2)] %>% matrix()
  ssrs<-sapply(strsplit(names(f.all)[c(-1,-2)],"\\$"),function(x) as.numeric(x[2]))

  #names(f.mat)[2:(length(f.mat)-1)] <- ssrs
  map <- readxl::read_excel(mapfile)
  map.data <-
    map[match(ssrs, map$SSR), which(colnames(map) %in% MAPFILE_ATTRIBUTES)]

  map.data.ps <- phyloseq::sample_data(map.data)
  phyloseq::sample_names(map.data.ps) <- map.data$SSR
  # phyloseq::sample_data(map.data.ps)$Date <- sapply(lapply(sample_data(map.data.ps)[['Date']],
  #                                                           as.Date,
  #                                                           origin = "1899-12-30"), # not sure why this origin gives correct answers
  #                                                    as.character)
  #row.names(map.data) <- map.data$SSR # todo: maybe delete this line?  not necessary to set rownames?
  f.taxtable <-
    phyloseq::build_tax_table(lapply(f.all$tax_name, phyloseq::parse_taxonomy_default))
  colnames(f.taxtable) <- c(Hmisc::upFirst(rank))
  dimnames(f.taxtable)[[1]] <- f.all$tax_name
  f.biome <- as.matrix(f.all[,c(-1,-2)])
  colnames(f.biome) <- ssrs
  rownames(f.biome) <- f.all$tax_name

  return(phyloseq::phyloseq(
    phyloseq::otu_table(f.biome, taxa_are_rows = TRUE),
    map.data.ps, #phyloseq::sample_data(map.data),
    phyloseq::tax_table(f.taxtable)
  ))

}

#' @title Make Phyloseq object from JSON files
#' @description  return a valid Phyloseq object created from the JSON files in flist and an Excel formatted mapfile
#' @param flist character vector of file names
#' @param count.normalized use the count_norm field instead of raw read numbers (default is FALSE)
#' @param mapfile XLSX filename that contains mapping info for the same SSRS in flist
#' @return valid Phyloseq object
#' @importFrom phyloseq sample_names taxa_names build_tax_table parse_taxonomy_qiime tax_table otu_table sample_data
#' @importFrom readxl read_excel
#' @export
phyloseq_from_JSON <- function(flist, mapfile, count.normalized = FALSE) {
  # return a valid Phyloseq object created from the JSON files in flist and an Excel formatted mapfile
  # Args:
  #   flist: character vector of file names.
  #   mapfile: XLSX filename that contains mapping info for the same SSRS in flist
  # Returns:  valid Phyloseq object

  f.all <- join_all_ubiome_files_full(flist, count.normalized)
  f.all[is.na(f.all)] <- 0
  f.mat <- f.all[, -1] %>% matrix()
  ssrs<-sapply(strsplit(names(f.all)[-1],"\\$"),function(x) as.numeric(x[2]))

  map <- readxl::read_excel(mapfile)
  map.data <-
    map[match(ssrs, map$SSR), which(colnames(map) %in% MAPFILE_ATTRIBUTES)]

  map.data.ps <- phyloseq::sample_data(map.data)
  phyloseq::sample_names(map.data.ps) <- map.data$SSR

  #row.names(map.data) <- map.data$SSR # todo: maybe delete this line?  not necessary to set rownames?
  f.taxtable <-
    phyloseq::build_tax_table(lapply(f.all$tax_name, phyloseq::parse_taxonomy_qiime))

  #dimnames(f.taxtable)[[1]] <- f.all$tax_name
  phyloseq::taxa_names(f.taxtable) <- f.all$tax_name
  f.biome <- as.matrix(f.all[,-1])
  colnames(f.biome) <- ssrs
  rownames(f.biome) <- f.all$tax_name

  return(phyloseq::phyloseq(
    phyloseq::otu_table(f.biome, taxa_are_rows = TRUE),
    map.data.ps, # phyloseq::sample_data(map.data),
    phyloseq::tax_table(f.taxtable)
  ))
}
