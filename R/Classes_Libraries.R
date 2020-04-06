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

#--------------------------------------------------------------
# metadata.MIxS
#--------------------------------------------------------------
# metadata.MIxS Class object
#--------------------------------------------------------------
setClass("metadata.MIxS", slots=list(data="data.frame", #a dataframe with the metadata (samples are rows, variables are columns)
                                     section="character", #a character vector containing a MIxS section for each term (variable)
                                     units="character", #a character vector containing a unit for each variable
                                     type="character", #versatile (loose following of MIxS) or strict.MIxS (following all the MIxS rules)
                                     env_package="character", #the environmental package. Can be not_specified or multiple_packages.
                                     QC="logical" #Has the data been Quality comtrolled (TRUE/FALSE)
)
)

#--------------------------------------------------------------
# metadata.MIxS Methods
#--------------------------------------------------------------
setMethod("show",
          "metadata.MIxS",
          function(object) {
            N <- nrow(object@data)
            C <- ncol(object@data)
            cat(paste("a metadata.MIxS class object with", as.character(N), "samples and",
                      as.character(C), "variables\n"))
            if(object@QC){
              cat("\tThe quality of this object has been checked.\n")
            }else{
              cat("\tThe quality of this object has NOT been checked.\n")
            }
          }
)

check.valid.metadata.MIxS <- function(d){
  valid <- TRUE
  #class must be metadata.MIxS
  if(!class(d)=="metadata.MIxS"){
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
# metadata.DwC
#--------------------------------------------------------------
# metadata.DwC Class Object
#--------------------------------------------------------------
setClass("DwC.event", slots=list(core="data.frame", #a dataframe with the DwC event core (events are rows, variables are columns)
                                occurrence="data.frame", #a dataframe with the occurrence extension
                                emof="data.frame", #a dataframe with the DwC eMoF extention (=the environmental data)
                                EML.url="character", #the url to the EML file from the IPT
                                QC="logical" #Has the data been Quality comtrolled (TRUE/FALSE)
)
)

setClass("DwC.occurrence", slots=list(core="data.frame", #a dataframe with the DwC event core (events are rows, variables are columns)
                                      emof="data.frame", #a dataframe with the DwC eMoF extention (=the environmental data)
                                      EML.url="character", #the url to the EML file from the IPT
                                      QC="logical" #Has the data been Quality comtrolled (TRUE/FALSE)
)
)

#--------------------------------------------------------------
# metadata.DwC Methods
#--------------------------------------------------------------
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
            cat(paste("a  DwC.event class  object with a core of ", as.character(N),
                      " events", as.character(E), ".\n", sep=""))
            if(EML.url!=""){
              cat(paste("\tThe metadata can be found at", EML.url,
                        "\n", sep=""))
            }
            if(object@QC){
              cat("\tThe quality of this object has been checked.\n")
            }else{
              cat("\tThe quality of this object has NOT been checked.\n")
            }
          }
)
setMethod("show",
          "DwC.occurrence",
          function(object) {
            N <- nrow(object@core)
            if(nrow(object@emof)>0){
              E <- "and an eMoF extension"
            }else{ E <- ""}
            cat(paste("a  DwC.occurrence  class object with a core of ", as.character(N),
                      " occurrences", as.character(E), ".\n", sep=""))
            if(EML.url!=""){
              cat(paste("\tThe metadata can be found at", EML.url,
                        "\n", sep=""))
            }
            if(object@QC){
              cat("\tThe quality of this object has been checked.\n")
            }else{
              cat("\tThe quality of this object has NOT been checked.\n")
            }
          }
)



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
