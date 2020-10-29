#==============================================================
# Author Maxime Sweetlove
# lisence CC 4.0
# Part of the POLA3R website (successor or mARS.biodiversity.aq)
# version 1.0 (2019-09-20)
# file encdong UTF-8
#
#==============================================================

#--------------------------------------------------------------
# classes for metadata formatted as MIxS
#--------------------------------------------------------------
# The Minimum Information on any Sequence (MIxS) standard was designed by the Genomics Standards Consortium to archive metadata and environmental data that is associated with sequence data.
# the MIxS.metadata class is made to work with this type of data in R

#' an S4 class to document data formated in the MIxS standard
#' @slot data a data.frame with the meta- and/or environmental data, formatted with samples/events as rows, and variables/MIxS terms as columns.
#' @slot section a character, a vector with the same length and order as number of columns of the dataframe in the data slot, containing a MIxS section/cathegory for each term
#' @slot units a character. a vector with the same length and order as number of columns of the dataframe in the data slot, containing a unit for each variable (use "alphanumeric" if there is no unit)
#' @slot type a character. Either "versatile", meaning the object is a loose following of the MIxS standard, or "strict.MIxS", meaning the object strictly follows all the MIxS rules.
#' @slot env_package a character. The MIxS environmental package that is appropriate for this dataset. Can also be not_specified or multiple_packages.
#' @slot QC boolean to indicate that the data been Quality controlled (TRUE) or not (FALSE) 
setClass("MIxS.metadata", slots=list(data="data.frame", 
                                     section="character", 
                                     units="character", 
                                     type="character", 
                                     env_package="character", 
                                     QC="logical" 
)
)

# methods for MIxS.metadata
#' the show method for MIxS.metadata
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @method show MIxS.metadata
#' @export
setMethod("show",
          "MIxS.metadata",
          function(object) {
            N <- nrow(object@data)
            C <- ncol(object@data)
            message(paste("a MIxS.metadata class object with", as.character(N), "samples and",
                      as.character(C), "variables\n"))
            if(object@QC){
              message("\tThe quality of this object has been checked.\n")
            }else{
              message("\tThe quality of this object has NOT been checked.\n")
            }
          }
)

#' validate a MIxS.metadata class object
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @description validates if a MIxS.metadata object is implemented correctly
#' @usage check.valid.MIxS.metadata(d)
#' @param d a MIxS.metadata object
#' @return a boolean, TRUE when valid, FALSE when not
#' @export
check.valid.MIxS.metadata <- function(d){
  valid <- TRUE
  #class must be MIxS.metadata
  if(!class(d)=="MIxS.metadata"){
    valid <- FALSE
  }else{
    dcol <- ncol(d@data)
    drow <- nrow(d@data)
    #The meta must have at least 1 sample (row) and 1 variable (column)
    if(dcol<=0 & drow<=0){
      valid <- FALSE
    }
    # each variable must have an associated unit
    if(length(d@units) < dcol |
       length(d@section) < dcol |
       length(d@units) != length(d@section)){
      valid <- FALSE
    }
    # only 2 types allowed: strict MIxS or a versatile form that not strictly follows the MIxS rules
    if(!d@type %in% c("versatile", "strict.MIxS")){
      valid <- FALSE
    }
    # the environmental package must be valid
    if(!d@env_package %in% c("air", "built_environment", "host_associated", "human_associated",
                             "human_gut", "human_oral", "human_skin", "human_vaginal",
                             "microbial_mat_biofilm", "miscellaneous_natural_or_artificial_environment",
                             "plant_associated", "sediment", "soil", "wastewater_sludge", "water",
                             "not_specified", "multiple_packages")){
      valid <- FALSE
    }
  }
  return(valid)
}

#--------------------------------------------------------------
# classes for data formatted as DarwinCore
#--------------------------------------------------------------
# DarwinCore a data standard for ecologcal and biodiversity data that was develloped by TDWG. It differs from MIxS in that is centred around observations of species, and not sequence or environmental data.
# flavour 1: the DarwinCore EventCore
# In the EventCore, data is structured as (hierarchical and nested) events (the core), that can be associated with different types of information (e.g. environmental measurements (emof extension), species occurrences (occurrence extension) or sequence data)

#' an S4 class to document data formated in the DarwinCore event core format
#' @slot core a data.frame. The DwC event core, formatted with samples/events as rows, and variables/DwCTerms as columns.
#' @slot occurrence a data.frame. The occurrence extension, formatted with occurrences as rows and variables/DwCTerms as columns. Must contain "eventID", "occurrenceID" and "basisOfRecord"
#' @slot emof a data.frame. The DwC extended Measurement or Fact (eMoF) extention. This can be used to store environmental data using the MIxS vocabulary.
#' @slot EML.url a character. the url to the Ecological Metadata Language (EML) file from the IPT
#' @slot QC boolean to indicate that the data been Quality controlled (TRUE) or not (FALSE) 
setClass("DwC.event", slots=list(core="data.frame",
                                occurrence="data.frame", 
                                emof="data.frame",
                                EML.url="character",
                                QC="logical"
)
)

# flavour 2: the DarwinCore OccurrenceCore
# In the OccurrenceCore, data is structured as occurrences of species (core), that can be associated with environmental measurements (emof extension)

#' an S4 class to document data formated in the DarwinCore occurrence core format
#' @slot core a data.frame. The DwC occurrence core, formatted with occurrences as rows and variables/DwCTerms as columns. Must contain "eventID", "occurrenceID" and "basisOfRecord"
#' @slot emof a data.frame. The DwC extended Measurement or Fact (eMoF) extention. This can be used to store environmental data using the MIxS vocabulary.
#' @slot EML.url a character. the url to the Ecological Metadata Language (EML) file from the IPT
#' @slot QC boolean to indicate that the data been Quality controlled (TRUE) or not (FALSE) 
setClass("DwC.occurrence", slots=list(core="data.frame", #a dataframe with the DwC event core (events are rows, variables are columns)
                                      emof="data.frame", #a dataframe with the DwC eMoF extention (=the environmental data)
                                      EML.url="character", #the url to the EML file from the IPT
                                      QC="logical" #Has the data been Quality comtrolled (TRUE/FALSE)
)
)

# methods for DarwinCore class objects
#' the show method for DwC.event
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @method show DwC.event
#' @export
setMethod("show",
          "DwC.event",
          function(object) {
            N <- nrow(object@core)
            if(nrow(object@emof)>0 & nrow(object@occurrence)>0){
              E <- "and an occurrence and eMoF extensions"
            } else if(nrow(object@emof)>0 & nrow(object@occurrence)==0){
              E <- "and an eMoF extension"
            }else if(nrow(object@emof)==0 & nrow(object@occurrence)>0){
              E <- "and an occurrence extension"
            }else{ E <- ""}
            message(paste("a  DwC.event class  object with a core of ", as.character(N),
                      " events", as.character(E), ".\n", sep=""))
            if(EML.url!=""){
              message(paste("\tThe metadata can be found at", EML.url,
                        "\n", sep=""))
            }
            if(object@QC){
              message("\tThe quality of this object has been checked.\n")
            }else{
              message("\tThe quality of this object has NOT been checked.\n")
            }
          }
)

#' the show method for DwC.occurrence
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @method show DwC.occurrence
#' @export
setMethod("show",
          "DwC.occurrence",
          function(object) {
            N <- nrow(object@core)
            if(nrow(object@emof)>0){
              E <- "and an eMoF extension"
            }else{ E <- ""}
            message(paste("a  DwC.occurrence  class object with a core of ", as.character(N),
                      " occurrences", as.character(E), ".\n", sep=""))
            if(EML.url!=""){
              message(paste("\tThe metadata can be found at", EML.url,
                        "\n", sep=""))
            }
            if(object@QC){
              message("\tThe quality of this object has been checked.\n")
            }else{
              message("\tThe quality of this object has NOT been checked.\n")
            }
          }
)


#' validate a DarwinCore class object
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @description validates if a DwC.event or DwC.occurrence object is implemented correctly
#' @usage check.valid.metadata.DwC(d)
#' @param d a DwC.event or DwC.occurrence object
#' @return a boolean, TRUE when valid, FALSE when not
#' @export
check.valid.metadata.DwC <- function(d){
  valid <- TRUE
  #class must be DwC.event or DwC.occurrence
  if(!class(d) %in% c("DwC.event", "DwC.occurrence")){
    valid <- FALSE
  }else{
    dcol <- ncol(d@core)
    drow <- nrow(d@core)
    #The core must have at least 1 sample (row) and 1 variable (column)
    if(dcol<=0 & drow<=0){
      valid <- FALSE
    }
    #type must be occurence or event
    if(class(d) == "DwC.occurrence" &
       !"occurrenceID" %in% colnames(d@core) |
       !"basisOfRecord" %in% colnames(d@core)){
      valid <- FALSE
    }else if(class(d) == "DwC.event" &
             !"eventID" %in% colnames(d@core)){
      valid <- FALSE
    }
  }
  return(valid)
}



#' Data Output: write a MIxS.metadata class object to a CSV file
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @description write a MIxS.metadata class object to a CSV file
#' @usage write.MIxS(x, file="")
#' @param x the MIxS.metadata class object to be written.
#' @param file either a character string naming a file or a connection open for writing. "" indicates output to the console.
#' @export
write.MIxS <- function(x, file = ""){
  if(!check.valid.MIxS.metadata(x)){
    stop("Input needs to be a MIxS.metadata class object.")
  }
  out.data <- rbind(x@section, x@units, x@data)
  rownames(out.data)[1]<-"section"
  rownames(out.data)[2]<-"units"
  write.csv(out.data, file=file, quote=FALSE, 
            row.names = TRUE, na = "")
  
}
