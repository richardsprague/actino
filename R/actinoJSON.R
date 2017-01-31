# actinoJSON.R
# handy functions for reading JSON files

# read uBiome JSON files and convert to tables that can be used by BIOM

# converts every json file in directory to a CSV file

# you can convert a single JSON file to CSV like this:
# convert_json_files_to_csv(pattern="sprague-nose-1410.json")


## Dependencies
# the following packages must be loaded before using this function
#install.packages("jsonlite")


tax_abbrev<-function(tax_rank){
  # returns the abbreviated character form of the character vector "tax_rank"
  # Needed by functions that convert tax_names to QIIME (or other) format.
  # Args:
  #   tax_rank: a character vector (e.g. "genus")
  # Returns:
  #   A one character abbreviation of tax_rank (e.g. "g")
  taxes<-unlist(strsplit("root superkingdom genus species order family phylum class superphylum subclass suborder no_rank species_group", split=" "))
  tax_abbrev<-unlist(strsplit("r k g s o f p c l k r n g", split=" "))
  tax<-data.frame(taxes,tax_abbrev)
  as.character(tax[tax$taxes==tax_rank,2])
}

full_taxa<- function(df,taxaRow){
  # # given a row taxaRow in df, return a character vector showing the full set of tax_ranks

  parent = as.numeric(taxaRow[match("parent",names(df))])#taxaRow[2])
  taxName = taxaRow[match("tax_name",names(df))] #taxaRow[5]
  taxRank = taxaRow[match("tax_rank",names(df))] # taxaRow[6]
  if(is.na(taxName)){browser(text="Null Taxname")} # some kind of error occurred.
  if(taxName == "root")
    return("Root")
  else
    return(
      paste(full_taxa(df,as.character(df[df$taxon==parent,])),
            ";",
            tax_abbrev(taxRank),"__",
            taxName,
            sep=""
      )
    )

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
  list.files(d, full.names = TRUE, pattern="[[:alnum:]].json")
}



#returns a big dataframe joining all tax_names from a list of JSON files
join_all_ubiome_files <- function(flist,tax_rank="genus",count.normalized=FALSE){
  f.all<-read_ubiome_json(flist[1],rank=tax_rank,count.normalized)
  for(f in flist[-1]){
    f.part<-read_ubiome_json(f,rank=tax_rank,count.normalized)
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
#' @return data frame representation of the sample
read_ubiome_json_full<-function(fname,count.normalized=FALSE){


  # convenience wrapper that returns a character name based on sample date and the notes field
  # expects dateVal to be a date of type character
  name_for_sample <- function(dateVal,notesVal,ssr){
    paste(as.character(as.Date(dateVal)),"$",ssr,"$",substring(notesVal,6,10))
  }
  j<-fromJSON(fname)

  rj<-j[["ubiome_bacteriacounts"]]
  p<-apply(rj,1,function (x) full_taxa(rj,x))
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
#' @return data frame representation of the sample
#returns a dataframe of just the genus and the count (or count_norm if count.normalized is TRUE)
read_ubiome_json<-function(fname,rank="genus",count.normalized=FALSE){

  # convenience wrapper that returns a character name based on sample date and the notes field
  # expects dateVal to be a date of type character
  name_for_sample <- function(dateVal,notesVal,ssr){
    paste(as.character(as.Date(dateVal)),"$",ssr,"$",substring(notesVal,6,10))
  }
  j<-fromJSON(fname)


  rj<-filter(j[["ubiome_bacteriacounts"]],tax_rank==rank)
  r<- data.frame(rj$tax_name,rj$tax_rank,reads=if(count.normalized) rj$count_norm else rj$count,stringsAsFactors = FALSE)

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
  f.all<-read_ubiome_json_full(flist[1],count.normalized,site )
  for(f in flist[-1]){
    f.part<-read_ubiome_json_full(f,count.normalized,site)
    if(nrow(f.part)>0 & nrow(f.all)>0)  f.all<-full_join(f.all,f.part)
    else if (nrow(f.part)>0 & nrow(f.all)==0) f.all<-f.part
  }
  f.all
}



