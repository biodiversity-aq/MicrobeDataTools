#==============================================================
# The MicrobeDataTools package
#       MicrobeDataTools is a collection data management tools for microbial 'omics datasets.
#       They allow to download, structure, quqlity-controll and standardize microbial datasets
#==============================================================
# Author Maxime Sweetlove
# lisence CC 4.0
# Part of the POLA3R website (successor or mARS.biodiversity.aq)
# version 1.0 (2019-09-20)
# file encdong UTF-8
#
#==============================================================
# Tools for general use
#==============================================================
# internal tools
#==============================================================
#' collect warning messages troughout a running function
#' @author Maxime Sweetlove ccBY 4.0 2019
#' @description collect warning messages troughout a running function
#' @usage multi.warnings(message_text, warningmessages)
#' @param message_text a character string. The error message to be added
#' @param warningmessages a vector with one ore more character strings. The previous error messages.
#' @return a character string with unique warning messages
multi.warnings <- function(message_text, warningmessages){
  if(! message_text %in% warningmessages){
    warningmessages <- c(warningmessages, message_text)
  }
  return(warningmessages)
}

#' Convert a taxon name to an NCBI taxon ID
#' @author Maxime Sweetlove ccBY 4.0 2019
#' @description converts taxon names of common taxa (superkingdom and phylum level) to it's NCBI taxID using an internal library. For taxa not in the internal library, please see https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi
#' @param taxon character or character vector. The taxon names to be converted to NCBI tax IDs.
#' @param fill.unknown character. The string to return when a taxon was not found in the list with common taxa. Default is NCBI:txid12908 (the ID for "unknown sequences"), other options include NA or ""
#' @details This function makes use of a limited internal library, not all the taxIDs are present
#' @return a character or character vector with the matching taxon IDs, with fill.unknown for those not found.
commonTax.to.NCBI.TaxID<-function(taxon, fill.unknown="12908"){

  #note: ID 12908 stands for "unclassified sequences"
  taxIDs<-c()
  failed_taxa<-c()
  for(tx in taxon){
    if(tolower(tx) %in% rownames(TaxIDLib)){
      if(tolower(tx) %in% c("eukaryotes", "eukarya", "eukaryote")){
        tx<-"eukaryota"
      }
      taxIDs <- c(taxIDs, as.character(TaxIDLib[tolower(tx),]$NCBItaxID))
    }else{
      failed_taxa <- tx
      taxIDs <- c(taxIDs, fill.unknown)
    }
  }
  if(length(failed_taxa)>0){
    messagex <- paste("The following taxa were not found under the most common taxa: ",
                      paste(unique(failed_taxa), collapse=", "), "\n an NCBI taxon ID for these taxa can be found at https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi",
                      sep="")
    warning(messagex)
  }
  return(taxIDs)
}

#' convert a MIxS term into it's ENA variant
#' @author Maxime Sweetlove ccBY 4.0 2019
#' @description get the ENA variant a MIxS term
#' @param variable character a MIxS term.
#' @return a character
get.ENAName <- function(variable){
  ENAName <- as.character(TermsLib[TermsLib$name==variable,]$name_variant_ENA)
  return(ENAName)
}

#' find a column in a dataframe
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @description find a column in a dataframe where it is a priori unkown wether the column is present or not.
#' @param dataset data.frame. A dataframe where to find a column based on the columnames
#' @param TermsList a character vector with the columnnames to look for. They will be handled in the order given, and the search will stop at the first match
#' @details An internal function to find information in incomming data of which the format and content is unknown
#' @return a character vactor with the content of the matching column (or the first match), or an empty vector
find.dataset <- function(dataset, TermsList=c()){

  ctu<-TRUE
  i=1
  colMatch<-c("")
  while(ctu){
    if(TermsList[i] %in% colnames(dataset)){
      colMatch <- dataset[,TermsList[i]]
      ctu<-FALSE
    }else{
      i=i+1
    }
    if(i>length(TermsList)){
      ctu<-FALSE
    }
  }
  return(colMatch)
}

#==============================================================
# tools for standardizing data content
#==============================================================

#' convert a coordinate in degrees to decimal
#' @author Maxime Sweetlove CC-BY 4.0 2019
#' @description Turns a latutude or longitude value in a degrees-minutes-seconds (DMS) format into a decimal value
#' @usage degree.to.decimal(val)
#' @param val a character string. A single latitude or longitude value to be transformed. Can include non-numeric character like the degree symbol, N-S-E-W wind directions,...
#' @details N-S-W-E wind directions as well as degrees, minutes and seconds are recognized and turned into a numeric decimal coordinate value
#' @return a numeric value
#' @export
coordinate.to.decimal<-function(val){
  val<-as.character(val)

  # remove illegal characters:
  val <- gsub(" ", "", val)[[1]]
  val <- gsub("~", "", val, fixed=TRUE)[[1]]
  val <- gsub("%", "", val, fixed=TRUE)[[1]]
  val <- gsub("$", "", val, fixed=TRUE)[[1]]
  val <- gsub("#", "", val, fixed=TRUE)[[1]]
  val <- gsub("*", "", val, fixed=TRUE)[[1]]
  val <- gsub("??", "", val, fixed=TRUE)[[1]]
  val <- gsub("?", "", val, fixed=TRUE)[[1]]
  val <- gsub("!", "", val, fixed=TRUE)[[1]]
  val <- gsub("ca.", "", val, fixed=TRUE)[[1]]
  Encoding(val)<-"UTF-8"

  if(grepl("S|W", val)){
    s=-1
    val <- gsub("S|W", "", val)[[1]]
  }else{
    s=1
    val <- gsub("N|E", "", val)[[1]]
  }

  # identify values:
  # replace all the weird reincarnations of the latitude symbols with something understandable
  # need to use fixed=T, so can't do it in a single line...
  # DEGREE == DDD
  val <- gsub('??', "DDD", val, fixed=TRUE)[[1]]

  val <- gsub('\302\260', "DDD", val, fixed=TRUE)[[1]]
  val <- gsub("\\241", "DDD", val, fixed=TRUE)[[1]]
  #val <- gsub("\241", "DDD", val, fixed=TRUE)[[1]] #single slash gave error...
  val <- gsub("\u00b0", "DDD", val, fixed=TRUE)[[1]]
  val <- gsub("\u00c2", "DDD", val, fixed=TRUE)[[1]]
  val <- gsub("\u00c2\u00b0", "DDD", val, fixed=TRUE)[[1]]

  val <- gsub("<U+00B0>", "DDD", val, fixed=TRUE)[[1]]
  val <- gsub("<U+00C2>", "DDD", val, fixed=TRUE)[[1]]
  val <- gsub("<c2><b0>>", "DDD", val, fixed=TRUE)[[1]]

  # MINUTE == MMM
  val <- gsub('\'', "MMM", val, fixed=TRUE)[[1]]
  val <- gsub("'", "MMM", val, fixed=TRUE)[[1]]

  val <- gsub('\342\200\262', "MMM", val, fixed=TRUE)[[1]]
  val <- gsub('\u00b4', "MMM", val, fixed=TRUE)[[1]]

  val <- gsub('<U+00B4>', "MMM", val, fixed=TRUE)[[1]]

  # SECOND == SSS
  val <- gsub('\"', "SSS", val, fixed=TRUE)[[1]]

  val <- gsub('\342\200\263', "SSS", val, fixed=TRUE)[[1]]

  # now interpret the values
  if(grepl("DDD", val)){
    degrees <- as.numeric(strsplit(val, "DDD")[[1]][1])
    degrees2 <- strsplit(val, "DDD")[[1]][2]
    if(grepl("MMM", degrees2)){
      minutes <- as.numeric(strsplit(degrees2, "MMM")[[1]][1])
      minutes2 <- strsplit(degrees2, "MMM")[[1]][2]
      if(grepl('SSS', minutes2)){
        seconds <- as.numeric(strsplit(minutes2, 'SSS')[[1]][1])
      }else{
        seconds <- 0
      }
    } else{
      minutes <- seconds <- 0
    }
    decimal <- s*(degrees + minutes/60 + seconds/3600)
  }else{
    decimal <- s*(as.numeric(val))
  }
  return(decimal)
}


#' get a rectagular bounding box from a set of coordinates
#' @author Maxime Sweetlove CC-BY 4.0 2019
#' @description finds the (rectangular) bounding box given decimal latitudes and longitudes
#' @usage get.boundingBox(latitudes, longitudes)
#' @param latitudes numeric vector. one ore more decimal latitude values
#' @param longitudes numeric vector. one ore more decimal longitude values
#' @details Documenting data with metadata on an Integrated Publishing Toolkit (IPT) typically requires a bounding box to simplify the goegraphic content of the dataset. This function is made to quickly calculate this.
#' @return a character string with the boundingbox printed to the console
#' @export
get.boundingBox<-function(latitudes, longitudes){

  latitudes<- as.numeric(latitudes)
  longitudes<- as.numeric(longitudes)

  S<-paste("South = ", min(latitudes))
  N<-paste("North = ", max(latitudes))

  E<-paste("East = ", max(longitudes))
  W<-paste("West = ", min(longitudes))

  cat(paste(S, "\n", N, '\n\n', W, "\n", E, "\n", sep=""))
  }

#' Get the sequence length of the first sequence in a fastq or fastq.gz file
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @description get the insert size (i.e. number of bases) of a fastq or fastq.gz file
#' @param file_path the full file path to the file
#' @details submitting sequence to a database of the Nucleotide Sequence Database Consortium typically requires to give the insert size of the sequences. This function gets the length of the first sequence in a file, assuming all sequences are of the same length.
#' @return a numeric value
#' @export
get.insertSize <- function(file_path=NA){
  con <- file(file_path,"r")
  first_seq <- readLines(con,n=2)
  close(con)
  insertSize <- paste("the insert size of the sequences is:   ", nchar(first_seq[[2]]), sep="")
  return(insertSize)
}

#==============================================================
# tools for getting help
#==============================================================
#' Find the MIxS or DarwinCore standard term and definition of a variable
#' @author Maxime Sweetlove ccBY 4.0 2019
#' @description retrieve the MIxS or DarwinCore standard term and definition of a variable.
#' @param term a character string. The variable to look for among the MIxS and DarwinCore vocabularies
#' @details Standerdizing microbial sequence data, metadata and environmental data can be quite difficult given the plethora of standard terms already in existance. This function returns a definition of any term that is used on the POLAAAR portal at biodiversity.aq.
#' @return chracater string printed to the console. The best matching terms and their the definitions.
#' @export
term.definition <- function(term){
  if(term %in% TermsLib$name){
    def_out <- as.character(TermsLib[TermsLib$name==term,]$definition)
    out_message <- paste(term, "\t\n", def_out)
  }else{
    potential_match <- c()
    def_out <- c()
    for(ls in TermsSyn){
      if(grepl(term, ls)){
        lst <- names(TermsSyn[ls[1]])
        potential_match <- c(potential_match, lst)
        def_out <- c(def_out, as.character(TermsLib[TermsLib$name==lst,]$definition))
      }
    }

    if(length(potential_match)>1){
      out_message <- paste("Multiple matches found for \"" , term, "\"\n\n",
                           "\t", paste(c(1:length(potential_match)),
                                       rep(". ",length(potential_match)),
                                       potential_match,
                                       rep(" :\n",length(potential_match)),
                                       rep("\"",length(potential_match)),
                                       def_out,
                                       rep("\"",length(potential_match)),
                                       collapse="\n\t", sep=""), sep="")
    }else if(length(potential_match)==1){
      out_message <- paste("Best match for \"" , term, "\"\n\n",
                                 potential_match, "\n", def_out, sep="")
    }else{
      out_message <-  paste("Could not find any matches for \"" , term, "\"\n",sep="")
    }
  }
  cat(out_message)

}



