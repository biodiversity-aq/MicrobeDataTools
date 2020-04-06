#==============================================================
# The MicrobeDataTools package
#       MicrobeDataTools is a collection data management tools for microbial 'omics datasets.
#       They allow to download, structure, quqlity-controll and standardize microbial datasets
#==============================================================
# Author Maxime Sweetlove
# lisence CC 4.0
# Part of the POLA3R website (successor or mARS.biodiversity.aq)
# version 1.0 (2020-01-28)
# file encdong UTF-8
#
#==============================================================
# data Quality Controll (QC) for DarwinCore
#==============================================================

#' Quality Controll (QC) for DarwinCore data
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @description performs a general QC on a darwinCore file. Desired output is either and event (core) file, occurrence (core or extension) or extended Measurement or Fact (eMoF; extension) file. Not an Ecological Metadata Language (EML) file.
#' @usage dataQC.DwC_general(DwC.data = NA, DwC.type = "event", ask.input = TRUE, complete.data=TRUE)
#' @param DwC.data data.frame. A dataframe with the data, structured with DarwinCore terms
#' @param DwC.type character. The type of DarwinCore of the output, either event or occurrence. If event, the output will have an event core with possible occurrence and eMoF extensions. If occurrence, the output will have an occurrence core with possibly an eMoF extension. Default event, if the parameter Event is not NA, out.type will be fixed as event.
#' @param ask.input logical. If TRUE, console input will be requested to the user when a problem occurs. Default TRUE
#' @param complete.data logical. If TRUE, datathat has not been provided, but can be completed automatically will be added to the DarwinCore data.frame. For instance, footprintWKT can be generated from the coordinates, or higher taxonomic level names (like kingdom, phylum,...) can looked up with the species name. Defaut TRUE.
#' @details DarwinCore is the biodiversity data standard develloped by TDWG, and is used by the Global Biodiversity Information Facility (GBIF). This function performs a basic and user-supervised quality control. This includes cheking all variables terms adhere to the DarwinCore vocabulary, listing other variables in an eMoF file or in the dynamicProperties field, and checking for obvious errors in the content of the data (typos, different NA values,...)
#' @return a dataframe that is formatted as either: an event or occurrence
#' @export
dataQC.DwC_general<-function(DwC.data = NA, DwC.type = "event", ask.input = TRUE,
                             complete.data=TRUE){
  #requires stringr and worrms
  warningmessages<-c()

  # check input
  if(tolower(DwC.type)=="event"){
    DwC.type<-"event"
    DwcType <- "Event"
    DwCLib<-"DwC_Event"
  }else if(tolower(DwC.type)=="occurrence"){
    DwC.type<-"occurrence"
    DwcType <- "Occurrence"
    DwCLib<-"DwC_Occurrence"
  }else if(tolower(DwC.type)=="emof"){
    DwC.type<-"emof"
    DwcType <- "eMoF"
    DwCLib<-"DwC_eMoF"
  }else{
    stop("invalid input for DwC.type")
  }

  # check columnnames
  termsQC <- dataQC.TermsCheck(observed=colnames(DwC.data),
                               exp.standard = "DwC", exp.section = DwC.type,
                               fuzzy.match = TRUE, out.type = "full")

  # run over fuzzy matches with possible solution
  if(length(termsQC$terms_wrongWithSolution)>0 & ask.input){
    for(tm in names(termsQC$terms_wrongWithSolution)){
      #note: in terms_wrongWithSolution: names= observed term, value is most likely match from standard
      tm_obs <- tm
      tm_match <- unname((termsQC$terms_wrongWithSolution)[tm])
      cat(paste("The column name \"", tm_obs,"\" is not a DarwinCore ", DwC.type, " term...\n", sep=""))
      cat(paste("\tdid you mean \"", tm_match,"\"? (y/n)\n", sep=""))
      doNext <- tolower(readline())
      if(doNext %in% c("y", "yes")){
        colnames(DwC.data)[colnames(DwC.data)==tm_obs] <- tm_match
      }else if(!doNext %in% c("n", "no")){
        stop("invalid input, expected n or y")
      }
    }
  }

  # run over non-matched terms
  if(length(termsQC$terms_notFound)>0 & ask.input){
    for(tm in termsQC$terms_notFound){
      cat(paste("The column name \"", tm,"\" is not a DarwinCore ", DwC.type, " term...\nPlease choose what to do:\n", sep=""))

      if(DwcType=="Event"){
        cat(paste("\t1) drop the term\n\t2) add to dynamicProperties\n\t3) add to eventRemarks\n\t4) add to fieldNotes\n", sep=""))
      }else if(DwcType=="Occurrence"){
        cat(paste("\t1) drop the term\n\t2) add to dynamicProperties\n\t3) add to eventRemarks\n\t4) add to fieldNotes",
                  "\n\t5) add to identificationRemarks\n\t6) add to taxonRemarks\n\t7) add to measurementRemarks\n\t8) add to occurrenceRemarks\n", sep=""))
      }
      doNext <- tolower(readline())
      if(doNext == 1){
        # drop the term
        DwC.data <- DwC.data[,!colnames(DwC.data) %in% tm]
      }else if(doNext == 2){
        #add to dynamicProperties
        dynProp <- as.character(DwC.data[,tm])
        dynProp <- sapply(dynProp, function(x){
          if(x!="" & !is.na(x)){
            x<-paste("{\"", tm, "\":", x, "}", sep="")
          }else{x}
        })
        dynProp<-unname(dynProp)
        if("dynamicProperties" %in% colnames(DwC.data)){
          dynProp <- paste(DwC.data$dynamicProperties, dynProp, sep=", ")
          dynProp <- sapply(dynProp, function(x){
            gsub(", $", "", x, fixed=FALSE)
          })
        }
        DwC.data$dynamicProperties <- dynProp
      }else if(doNext %in% 3:8){
        if(doNext == 3){
          #add to eventRemark
          doName <-"eventRemark"
        }else if(doNext == 4){
          #add to fieldNotes
          doName <-"fieldNotes"
        }else if(doNext == 5){
          #add to identificationRemarks
          doName <-"identificationRemarks"
        }else if(doNext == 6){
          #add to taxonRemarks
          doName <-"taxonRemarks"
        }else if(doNext == 7){
          #add to measurementRemarks
          doName <-"measurementRemarks"
        }else if(doNext == 8){
          #add to occurrenceRemarks
          doName <-"occurrenceRemarks"
        }
        remarks <- as.character(DwC.data[,tm])
        remarks <- sapply(remarks, function(x){
          if(x!="" & !is.na(x)){
            x<-paste(tm, ":", x, sep="")
          }else{x}
        })
        remarks<-unname(remarks)

        if(doName %in% colnames(DwC.data)){
          remarks <- paste(DwC.data[,doName], remarks, sep=", ")
          remarks <- sapply(remarks, function(x){
            gsub(", $", "", x, fixed=FALSE)
          })
        }
        DwC.data[,doName] <- remarks
      }else{
        if(DwcType=="Event"){
          stop("invalid input, expected 1, 2, 3 or 4")
        }else if(DwcType=="Occurrence"){
          stop("invalid input, expected 1, 2, 3, ... or 8")
        }
      }
    }
  }

  # check presence of required terms
  req_terms <- as.character(TermsLib[TermsLib[,DwCLib]==2,]$name)
  if(!all(req_terms %in% colnames(DwC.data))){
    terms_short <- setdiff(req_terms, colnames(DwC.data))
    warningmessages <- multi.warnings(paste("Please provide data for the following required DwC terms that were not found:\n\t", paste(terms_short, collapse=","), sep=""), warningmessages)
    for(tm in terms_short){
      DwC.data[,tm]<-""
    }
    ## user input for basisOfRecord
    if("basisOfRecord" %in% terms_short & ask.input){
      cat("No basisOfRecord found. Please specify if the data is:\n\t1)HumanObservation\n\t2)MachineObservation\n\t3)LivingSpecimen\n\t4)PreservedSpecimen\n\t5)FossilSpecimen")
      doNext <- tolower(readline())
      if(doNext == 1){
        DwC.data$basisOfRecord <- rep("HumanObservation", nrow(DwC.data))
      }else if(doNext == 2){
        DwC.data$basisOfRecord <- rep("MachineObservation", nrow(DwC.data))
      }else if(doNext == 3){
        DwC.data$basisOfRecord <- rep("LivingSpecimen", nrow(DwC.data))
      }else if(doNext == 4){
        DwC.data$basisOfRecord <- rep("PreservedSpecimen", nrow(DwC.data))
      }else if(doNext == 5){
        DwC.data$basisOfRecord <- rep("FossilSpecimen", nrow(DwC.data))
      }else{
        stop("invalid input, expected 1, 2, 3, 4 or 5")
      }
    }
    ## user input for eventID
    if("eventID" %in% terms_short & ask.input){
      cat("No eventID found.\n\tPlease type an eventID prefix, which will be used to generate a unique ID per event")
      doNext <- readline()
      if(doNext == ""){
        prefix<-"event_"
      }else{
        prefix<-paste(doNext, "_", sep="")
      }
      DwC.data$eventID <- paste("prefix",
                                stringr::str_pad(1:nrow(DwC.data), nchar(nrow(DwC.data)), pad = "0"),
                                sep="")
    }
    ## user input for occurrenceID
    if("occurrenceID" %in% terms_short & ask.input){
      cat("No occurrenceID found.\n\tPlease type an occurrenceID prefix, which will be used to generate a unique ID per event")
      doNext <- readline()
      if(doNext == ""){
        prefix<-"occ_"
      }else{
        prefix<-paste(doNext, "_", sep="")
      }
      DwC.data$occurrenceID <- paste("prefix",
                                     stringr::str_pad(1:nrow(DwC.data), nchar(nrow(DwC.data)), pad = "0"),
                                     sep="")
    }
  }

  # some other checks
  if(!"footprintWKT" %in% colnames(DwC.data) &
     complete.data &
     DwC.type == "Event" &
     "decimalLatitude" %in% colnames(DwC.data) &
     "decimalLongitude" %in% colnames(DwC.data)){
    DwC.data$footprintWKT <- dataQC.generate.footprintWKT(DwC.data)
    warningmessages <- multi.warnings("added footprintWKT", warningmessages)
  }
  #if("dynamicProperties" %in% colnames(DwC.data)){
  #  structure dynamicProperties
  #}

  # check latitude-longitude


  # dealing with the collection date, and putting it in the YYYY-MM-DD format
  if(!"eventDate" %in% colnames(DwC.data)){
    if("year" %in% colnames(DwC.data) && !grepl("/", DwC.data$year)){
      DwC.data$eventDate <- DwC.data$year
      if("month" %in% colnames(DwC.data)){
        for(i in nrow(DwC.data)){
          if(DwC.data[i,]$eventDate!="" &
             !is.na(DwC.data[i,]$eventDate)){
            if(DwC.data[i,]$month!="" &
               !is.na(DwC.data[i,]$month)){
              DwC.data[i,]$eventDate <- paste(DwC.data[i,]$eventDate, DwC.data[i,]$month, sep="-")
            }
          }else{
            DwC.data[i,]$eventDate <- ""
          }
        }
        if("day" %in% colnames(DwC.data)){
          for(i in nrow(DwC.data)){
            if(DwC.data[i,]$eventDate!=""){
              if(DwC.data[i,]$day!="" &
                 !is.na(DwC.data[i,]$day)){
                DwC.data[i,]$eventDate <- paste(DwC.data[i,]$eventDate, DwC.data[i,]$day, sep="-")
              }
            }
          }
        }
      }
    }
  }

  if("eventDate" %in% colnames(DwC.data)){
    QCDate <- dataQC.dateCheck(DwC.data, "eventDate")
  }
  warningmessages<-c(QCDate$warningmessages, warningmessages)
  if(length(QCDate$values)==nrow(DwC.data)){
    DwC.data$eventDate <- QCDate$values
  }

  #check taxa and look for scientificNameID
  if(DwC.type == "occurrence" & complete.data){
    # try to get the species names
    species <- dataQC.TaxonListFromData(DwC.data)

    species <- dataQC.taxaNames(species)

    # make a small table with the unique taxa to collect all the info from WORMS or GBIF
    cat("Completing the taxonomic information (this might take a while)...\n")
    taxid_key <- dataQC.completeTaxaNamesFromRegistery(species$scientificName)

    # now use the taxid_key table to complete the data in DwC.data
    for(term in setdiff(colnames(taxid_key), c("aphID"))){
      DwC.data[,term]  <- as.character(species)
      DwC.data[,term] <- unname(unlist(sapply(as.character(DwC.data[,term]),
                                              FUN = function(x){
                                                gsub(x,taxid_key[taxid_key$taxname==x,][,term],x)
                                              })))
    }
    warningmessages <- multi.warnings("added scientificNameID and additional species data", warningmessages)
  }

  # retrun the input file with the executed adjustments
  return(DwC.data)
}


#' format dataframes into a DarwinCore object
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @description Takes one or more dataframes or a DarwinCore Archive and performs a basic Quality Controll. (see details)
#' @usage DataQC.DwC(Event=NA, Occurrence=NA, eMoF=NA, EML.url=NA, out.type="event", ask.input=TRUE))
#' @param Event data.frame. A dataframe with the event data, structured with DarwinCore terms
#' @param Occurrence data.frame. A dataframe with the occurrence data, structured with DarwinCore terms
#' @param eMoF data.frame. A dataframe with the extended MeasurementOrFact (eMoF) data, structured with DarwinCore terms
#' @param out.type character. The type of DarwinCore of the output, either event or occurrence. If event, the output will have an event core with possible occurrence and eMoF extensions. If occurrence, the output will have an occurrence core with possibly an eMoF extension. Default event, if the parameter Event is not NA, out.type will be fixed as event.
#' @param ask.input logical. If TRUE, console input will be requested to the user when a problem occurs. Default TRUE
#' @param EML.url character. A URL to the Ecological Metadata Language (EML) file.
#' @details DarwinCore is the biodiversity data standard develloped by TDWG, and is used by the Global Biodiversity Information Facility (GBIF). This function performs a basic and user-supervised quality control. This includes cheking all variables terms adhere to the DarwinCore vocabulary, listing other variables in an eMoF file or in the dynamicProperties field, and checking for obvious errors in the content of the data (typos, different NA values,...). All dataframes are combined into a DwC.event or DwC.occurrence object
#' @return a DwC.event (event core with occurrence and/or eMoF extensions) or DwC.occurence (occurrence core with eMoF extension) object
#' @export
dataQC.DwC <- function(Event=NA, Occurrence=NA, eMoF=NA, EML.url=NA,
                       out.type="event", ask.input=TRUE){

  warningmessages<-c()

  #check out.format and  out.type
  if(!tolower(out.type) %in% c("event", "occurrence")){
    stop("invalid input for out.type")
  }

  # 1. checking the Event data.
  if(!all(is.na(Event))){
    # fix out.type as event, regardless of user input
    out.type <- "event"
    has.event <- TRUE
    eventQC <- dataQC.DwC_general(DwC.data = Event, DwC.type = "event", ask.input = ask.input)
  }else{
    has.event <- FALSE
    eventQC<-data.frame
  }
  # 2. checking the Occurrence data.
  if(!all(is.na(Occurrence))){
    has.occurrence <- TRUE
    occurrenceQC <- dataQC.DwC_general(DwC.data = Occurrence, DwC.type = "occurrence", ask.input = ask.input)
    if(out.type == "event" & !"eventID" %in% colnames(occurrenceQC)){
      stop("eventID missing from the occurrence extension. Please correct.")
    }
    if(!has.event){
      out.type <- "occurrence"
    }
  }else{
    has.occurrence <- FALSE
    occurrenceQC<-data.frame
  }
  # 3. checking the Extended Measurement or Fact data.
  if(!all(is.na(eMoF))){
    has.emof <- TRUE
    emofQC <- dataQC.DwC_general(DwC.data = eMoF, DwC.type = "emof", ask.input = ask.input)
  }else{
    has.emof <- FALSE
    emofQC<-data.frame()
  }

  # 4. check EML (if present)
  if(!is.na(EML.url)){
    if(!RCurl::url.exists(EML.url)){
      stop(paste("The following EML page does not seem to exist:\n\t", EML.url, sep=""))
    }
  }else{
    EML.url<-""
  }

  # 5. check IDs if multiple files provided
  if(tolower(out.type) == "event" & has.event){
    u_ev_ev <- unique(eventQC$eventID)
    if(has.occurrence){
      u_ev_occ <- unique(occurrenceQC$eventID)
      if(length(setdiff(u_ev_ev, u_ev_occ))>0){
        warningmessages <- multi.warnings("there are eventIDs that do not occur in the occurrence file", warningmessages)
      }else if(length(setdiff(u_ev_occ, u_ev_ev))>0){
        ev_diff <- setdiff(u_ev_occ, u_ev_ev)
        stop(paste("The following eventIDs in the occurrence table were not fount in the event core:\n\t", paste(ev_diff, collapse=", "), sep=""))
      }
    }
    if(has.emof){
      u_ev_emof <- unique(emofQC$eventID)
      if(length(setdiff(u_ev_emof, u_ev_ev))>0){
        ev_diff <- setdiff(u_ev_emof, u_ev_ev)
        stop(paste("The following eventIDs in the eMoF table were not fount in the event core:\n\t", paste(ev_diff, collapse=", "), sep=""))
      }
    }
  }else if(tolower(out.type) == "occurrence" & has.occurrence & has.emof){
    u_ev_occ <- unique(occurrenceQC$occurrence)
    if(length(setdiff(u_ev_emof, u_ev_occ))>0){
      ev_diff <- setdiff(u_ev_emof, u_ev_occ)
      stop(paste("The following occurrenceIDs in the eMoF table were not fount in the occurrence core:\n\t", paste(ev_diff, collapse=", "), sep=""))
    }
  }

  # 6. finalizing and formating the output
  if(tolower(out.type) == "event" & has.event){
    DwC_out <- new("DwC.event",
                   core = eventQC,
                   occurrence = occurrenceQC,
                   emof = emofQC,
                   EML.url=EML.url,
                   QC=TRUE)
  }else if(tolower(out.type) == "occurrence" & has.occurrence){
    DwC_out <- new("DwC.occurrence",
                   core = occurrenceQC,
                   emof = emofQC,
                   EML.url=EML.url,
                   QC=TRUE)
  }else{
    stop("out.type conflicts with provided data")
  }
  warning(warningmessages)
  return(DwC_out)

}


#==============================================================
# data Quality Controll (QC) for MIxS
#==============================================================

#' format dataframes into a MIxS object
#' @author Maxime Sweetlove ccBY 4.0 2019
#' @description takes a dataframe with contextual data and metadata from a sequencing dataset and performs a basis Quality Controll. (see details)
#' @usage DataQC.MIxS(metadata = NA, ask.input=TRUE, add_to = NA))
#' @param metadata data.frame. The raw metadata downloaded from INSDC to be cleaned up. Rows are samples, columns variables. Units can be listed in the first or second row, and will be automatically detected if the row names of this row includes the word "units". Different units per sample are not allowed.
#' @param ask.input logical. If TRUE, console input will be requested to the user when a problem occurs (process runs user-supervised). Default TRUE
#' @param add_to a metadata.MIxS object. An already present dataset with quality-comtrolled metadata. must be formatted as metadata.MIxS to ensure the correct input format of the data.
#' @details Any sequencing project typically has important additional data associated with it. This goes from laboratory protocols, sequencing platform settings or environmental measurements. Thisfunction was develloped to sort through these metadata (provided in a dataframe), and perform a basic quality controll, correcting the most common mistakes, like incorrectly formatting the geographic coordinates, formatting dates, typos or variants of variable names, etc. To do this, the function makes use of a build-in dictionary of (MIxS) terms and their synonyms (that is: spelling errors, writing differences, true synonyms,...). Note that it is possible some terms are not recognized. In that case contact the author of the package to update the dictionary in the upcomming version.
#' @seealso get.BioProject.metadata.INSDC, get.sample.attributes.INSDC
#' @return a metadata.MIxS object that is compatible with the MIxS standard
#' @export
dataQC.MIxS <- function(metadata = NA, ask.input=TRUE, add_to = NA){
  warningmessages<-c()

  # 0. pre-process input
  # 0.1. check input data
  if(!is.data.frame(metadata)){
    stop("The input must be a dataframe, with samples as rows/ variables as columns")
  }
  if(!is.na(add_to) && !check.valid.metadata.MIxS(add_to)){
    stop("The input for the add_to argument must be a metadata.MIxS object to ensure correct merging of the datasets")
  }


  # 0.2 formatting
  # remove empty columns
  metadata <- metadata[,colSums(is.na(metadata)) < nrow(metadata)]
  # clean up columnames
  colnames(metadata) <- gsub("[\\.]+", "_", colnames(metadata)) # replace double dots with underscore
  colnames(metadata) <- gsub("_$", "", colnames(metadata))  # remove trailing underscore
  metadata_origcolNames <- metadata
  colnames(metadata) <- tolower(colnames(metadata)) # all to lowercase

  # check columnnames, and correct errors
  termsQC <- dataQC.TermsCheck(observed=colnames(metadata),
                               exp.standard = "MIxS", exp.section = NA,
                               fuzzy.match = FALSE, out.type = "full")
  for(tQC in names(termsQC$terms_wrongWithSolution)){
    colnames(metadata)[colnames(metadata)==tQC] <- termsQC$terms_wrongWithSolution[tQC]
  }

  # 0.3 check if data is in one-header table or if there are additional MiMARKS header lines
  # for additional MIxS headers: "environmental package", "units template" => only units of importance
  pre_def_units <- FALSE
  if(grepl("unit", tolower(row.names(metadata)[1]))){
    units <- data.frame(var_name=c(colnames(metadata)), unit= unlist(metadata[1,]), stringsAsFactors = FALSE)
    metadata <- data.frame(metadata[-1,], stringsAsFactors = FALSE)
    metadata_origcolNames <- data.frame(metadata_origcolNames[-1,], stringsAsFactors = FALSE)
    pre_def_units <- TRUE
    warningmessages <- multi.warnings("the units were taken from the row \"units template\" in the input data", warningmessages)
  } else if(grepl("unit", tolower(row.names(metadata)[2]))){
    units <- data.frame(var_name=c(colnames(metadata)), unit= unlist(metadata[2,]), stringsAsFactors = FALSE)
    metadata <- data.frame(metadata[-c(1,2),], stringsAsFactors = FALSE)
    metadata_origcolNames <- data.frame(metadata_origcolNames[-c(1,2),], stringsAsFactors = FALSE)
    pre_def_units <- TRUE
    warningmessages <- multi.warnings("the units were taken from the row \"units template\" in the input data", warningmessages)
  }


  # 0.4 prepare output data:
  # make an empty output file to fill along the way
  New_metadata <- data.frame(row.names=rownames(metadata))

  # 1. looking for the original sample name
  metadataNames <- dataQC.findNames(dataset = metadata_origcolNames, ask.input=ask.input)
  New_metadata$original_name <- (metadataNames$Names)$original_name
  warningmessages <- multi.warnings(metadataNames$warningmessages, warningmessages)
  tryCatch({
    rownames(New_metadata) <- (metadataNames$Names)$original_name
  },
  error=function(x){
    warningmessages <- multi.warnings("duplicate or missing original sample names", warningmessages)
  })
  if(! all((metadataNames$Names)$INSDC_SampleID=="")){New_metadata$INSDC_SampleID <- (metadataNames$Names)$INSDC_SampleID}
  if(! all((metadataNames$Names)$eventID=="")){New_metadata$eventID <- (metadataNames$Names)$eventID}
  if(! all((metadataNames$Names)$parentEventID=="")){New_metadata$parentEventID <- (metadataNames$Names)$parentEventID}
  if(! all((metadataNames$Names)$occurrenceID=="")){New_metadata$occurrenceID <- (metadataNames$Names)$occurrenceID}


  # 2. some basic info from insdc
  TermsSyn_insdc<-TermsSyn[as.character(TermsLib[TermsLib$name_origin=="INSDC",]$name)]
  for(item in names(TermsSyn_insdc)){
    item_shared <- intersect(TermsSyn_insdc[item][[1]], colnames(metadata))
    if(length(item_shared)==1){
      New_metadata[,item] <- metadata[,item_shared]
    }
  }

  # 3. dealing with latitude-longitude, and it's many possible formats...
  TermsSyn_latlon<-TermsSyn[as.character(TermsLib[TermsLib$name=="lat_lon",]$name)]
  TermsSyn_lat<-TermsSyn[as.character(TermsLib[TermsLib$name=="decimalLatitude",]$name)]
  TermsSyn_lon<-TermsSyn[as.character(TermsLib[TermsLib$name=="decimalLongitude",]$name)]

  QClatlon <- dataQC.LatitudeLongitudeCheck(metadata,
                                            latlon.colnames=list(TermsSyn_latlon[[1]],
                                                                 TermsSyn_lat[[1]],
                                                                 TermsSyn_lon[[1]]))
  warningmessages<-c(QClatlon$warningmessages, warningmessages)
  New_metadata$lat_lon <- QClatlon$values
  New_metadata$decimalLatitude <- sapply(New_metadata$lat_lon, function(x){strsplit(x, " ")[[1]][1]})
  New_metadata$decimalLongitude <- sapply(New_metadata$lat_lon, function(x){strsplit(x, " ")[[1]][2]})

  New_metadata[is.na(New_metadata)] <- ""

  #change the units for the coordinates
  if(pre_def_units){
    for(unitx in c("lat_lon", "decimalLatitude", "decimalLongitude")){
      if(unitx %in% units$var_name){
        units[units$var_name==unitx,]$unit <- as.character(TermsLib[TermsLib$name==unitx,]$expected_unit)
      }else{
        units <- rbind(units, data.frame(var_name=unitx, unit=as.character(TermsLib[TermsLib$name==unitx,]$expected_unit)))
        rownames(units)[nrow(units)] <- unitx
      }
    }
  }


  # 4. dealing with the collection date, and putting it in the YYYY-MM-DD format
  TermsSyn_date<-TermsSyn[as.character(TermsLib[TermsLib$name=="collection_date",]$name)]
  QCDate <- dataQC.dateCheck(metadata, TermsSyn_date[[1]])
  warningmessages<-c(QCDate$warningmessages, warningmessages)
  if(length(QCDate$values)==nrow(metadata)){
    New_metadata$collection_date <- QCDate$values
  }

  # 5. the core MIxS terms
  TermsSyn_MIxS <- TermsSyn[as.character(TermsLib[TermsLib$MIxS_core>0,]$name)]
  TermsSyn_MIxS <- TermsSyn_MIxS[!names(TermsSyn_MIxS) %in% c("lat_lon", "collection_date")]
  for(item in names(TermsSyn_MIxS)){
    item_shared <- intersect(TermsSyn_MIxS[item][[1]], colnames(metadata))
    if(length(item_shared)==1){
      New_metadata[,item] <- metadata[,item_shared]
    } else if(length(item_shared)>1){
      New_metadata[,item] <- metadata[,item_shared[[1]]]
    }
  }

  # 6. the MIxS package terms
  # 6.1 find the best package/ask user if no package was specified
  env_package<-""
  if(!"env_package" %in% colnames(metadata)){
    if(ask.input){
      cat("No env_package was specified.\nPlease specify what to do next:\n1) Make an educated guess based on the data\n2) Ask user for the package name\n3) stop executing\n(type 1, 2 or 3)\n")
      doNext <- readline()
      if(doNext==1){
        env_package <- dataQC.guess.env_package.from.data(metadata)
        warningmessages<-c(warningmessages, env_package$warningmessages)
        env_package <- env_package$values
      }else if(doNext==2){
        cat("Please provide a single MIxS environmental package\nThe choices are: air, built_environment, host_associated, human_associated,human_gut,\nhuman_oral, human_skin, human_vaginal,microbial_mat_biofilm,\nmiscellaneous_natural_or_artificial_environment,\nplant_associated, soil, sediment, wastewater_sludge, water\n")
        env_package <- readline()
        if(!env_package %in% colnames(TermsLib)){
          stop("incorrect environmental package provided. Be sure to use underscores and lowercase letters")
        } else{
          env_package <- rep(env_package, nrow(metadata))
        }
      }else if(doNext==3){
        stop("you chose to interrupt execution.")
      }else{
        stop("incorrect input. Interrupted execution.")
      }
    }else{
      env_package <- dataQC.guess.env_package.from.data(metadata)
      warningmessages<-c(warningmessages, env_package$warningmessages)
      env_package <- env_package$values
    }
  }else{# 6.2 check if package is valid or if the ENA checklist number needs to be converted to a package
    New_metadata$env_package <- metadata$env_package
    if(length(setdiff(unique(New_metadata$env_package), ENA_checklistAccession$env_package))>0 |
       length(setdiff(unique(New_metadata$env_package), ENA_checklistAccession$ena_package))>0){
      #no correct package name, check if it is one of the ENA accession numbers
      if(sum(grepl("ERC", unique(New_metadata$env_package)))==length(unique(New_metadata$env_package))){
        #convert all the ENA checklist accession numbers to MIxS packages
        ENA_checklistAccession$ena_checklist_accession
        for(i in 1:nrow(New_metadata)){
          pk <- ENA_checklistAccession[ENA_checklistAccession$ena_checklist_accession == New_metadata[i,]$env_package,]$env_package
          if(length(pk)>0 && ENA_checklistAccession[ENA_checklistAccession$env_package==pk,]$MIxS==TRUE){
            New_metadata[i,]$env_package <- pk
          } else{
            New_metadata[i,]$env_package <- ""
            warningmessages<-multi.warnings("For some samples, the given environmental package did not correspond to a MIxS package", warningmessages)
          }
        }
      }else{
        env_package <- dataQC.guess.env_package.from.data(metadata)
        warningmessages<-c(warningmessages, env_package$warningmessages)
        env_package <- env_package$values
        warningmessages<-multi.warnings("It seems the way the environmatla package is written is not allowed, better correct", warningmessages)
      }
    }
  }

  if(all(nchar(env_package)==0)){
    warningmessages<-multi.warnings("No env_package could be inferred", warningmessages)
  } else{
    New_metadata$env_package <- env_package
  }

  TermsSyn_MIxSpackage<-TermsSyn[as.character(TermsLib[TermsLib$name_origin=="MIxS" & TermsLib$MIxS_core==0,]$name)]
  for(item in names(TermsSyn_MIxSpackage)){
    item_shared <- intersect(TermsSyn_MIxSpackage[item][[1]], colnames(metadata))
    if(length(item_shared)==1){
      New_metadata[,item] <- metadata[,item_shared]
    } else if(length(item_shared)>1){
      New_metadata[,item] <- metadata[,item_shared[[1]]]
    }
  }

  # 7. any other additionalinformation
  # 7.1 already registered terms
  TermsSyn_add <- TermsSyn[as.character(TermsLib[TermsLib$name_origin %in% c("DwC", "miscellaneous"),]$name)]
  TermsSyn_add <- TermsSyn_add[!names(TermsSyn_add) %in% c("decimalLatitude", "decimalLongitude")]
  for(item in names(TermsSyn_add)){
    item_shared <- intersect(TermsSyn_add[item][[1]], colnames(metadata))
    if(length(item_shared)==1){
      New_metadata[,item] <- metadata[,item_shared]
    }
  }
  # 7.2 novel terms
  unknown_terms <- setdiff(colnames(metadata), unlist(TermsSyn))
  if(length(unknown_terms) > 5){
    # if there are too much novel terms, don't go over all of them
    if(ask.input){
      t<-paste(unknown_terms, collapse=", ")
      cat(paste("The following unknown variables were encountered:\n",t," \n\tAdd all to the QC'd data? (y/n)\n", sep=" "))
      ctu <- readline()
      if(ctu %in% c("y", "Y", "yes", "YES", "Yes")){
        for(t in unknown_terms){
          New_metadata[,t] <- metadata[,t]
        }
      }
    }else{
      warningmessages<-multi.warnings("Some unknown variables present in the data", warningmessages)
      for(t in unknown_terms){
        New_metadata[,t] <- metadata[,t]
      }
    }
  }else{
    for(t in unknown_terms){
      if(ask.input){
        cat(paste("The following unknown variable was encountered:",t," \n\tAdd to the QC'd data? (y/n)\n", sep=" "))
        ctu <- readline()
        if(ctu %in% c("y", "Y", "yes", "YES", "Yes")){
          New_metadata[,t] <- metadata[,t]
        }
      }else{
        New_metadata[,t] <- metadata[,t]
        warningmessages<-multi.warnings("Some unknown variables present in the data", warningmessages)
      }
    }
  }

  # 8. some additional quality controll
  if("investigation_type" %in% colnames(New_metadata)){
    New_metadata$investigation_type <- sapply(New_metadata$investigation_type, FUN=function(x){
      if(x %in% c("AMPLICON", "amplicon", "metabarcode")){x<-"mimarks-survey"
      }else if(x=="WGS"){x<-"metagenome"}
      return(x)})
  }else{
    if(ask.input){
      cat("No investigation_type was found...\n\tPlease provide an investigation_type. Common ones include mimarks-survey or metagenome. Type n to ignore.\n")
      invtype <- readline()
      if(! invtype %in% c("n", "N")){
        New_metadata$investigation_type <- rep(invtype, nrow(New_metadata))
      }
    }
  }
  if("target_gene" %in% colnames(New_metadata)){
    New_metadata$target_gene <- sapply(New_metadata$target_gene, function(x){
      if(grepl("16S", x)){x<-"16S ssu rRNA"
      }else if(grepl("18S", toupper(x))){x<-"18S ssu rRNA"
      }else if(grepl("ITS", toupper(x))){x<-"ITS"
      }else if(grepl("COI", toupper(x))){x<-"COI"
      }
      return(x)})
  }
  if("specific_host" %in% colnames(New_metadata)){
    host_val <- setdiff(unique(New_metadata$specific_host), c(NA, "NA", "-", "not applicable"))
    if(identical(host_val, character(0))){
      New_metadata <- New_metadata[,!colnames(New_metadata) %in% "specific_host"]
    }
  }

  #if("eventID" %in% colnames(New_metadata)){
  #  if(length(unique(New_metadata$eventID))==nrow(New_metadata)){
  #    cat("There are the same number of events as there are samples...\n\tkeep $eventID? (y/n)\n")
  #    ctu <- readline()
  #    if(ctu %in% c("n", "N", "no", "NO", "No")){
  #      New_metadata <- New_metadata[,!colnames(New_metadata) %in% "eventID"]
  #    }
  #  }
  #}

  # QC on length/depth/size measurements
  # if multiple units found: do nothing
  # if a string found that is not one of the expected units: leave it and work with the rest
  size_related <- c("depth", "elev", "alt_elev", "filter_size", "tot_depth_water_col")
  items_shared <- intersect(size_related, colnames(New_metadata))
  QC_units<-c()
  if(length(items_shared)>0){
    alternative_units<-data.frame(name=items_shared, unit_full=rep(NA, length(items_shared)))
    possible_units <- c("nanometer", "micrometer", "milimeter", "centimeter", "decimeter", "meter", "kilometer")
    for(name_unit in items_shared){
      vals <- as.vector(as.character(New_metadata[,colnames(New_metadata) %in% name_unit]))
      names(vals)<-row.names(New_metadata)
      #extract the units
      if(pre_def_units){
        val_units <- as.character(units[units$var_name==name_unit,]$unit)
        val_units <- rep(val_units,nrow(New_metadata))
      }else{
        val_units <- unlist(lapply(vals, function(v){gsub("[0-9]|\\.|,|-", "", v)} ))
      }
      val_units <- gsub(" ", "", tolower(val_units))
      val_units <- gsub("^nm$", "nanometer", tolower(val_units), fixed=FALSE)
      val_units <- gsub("^um$", "micrometer", tolower(val_units), fixed=FALSE)
      val_units <- gsub("^mm$", "milimeter", tolower(val_units), fixed=FALSE)
      val_units <- gsub("^cm$", "centimeter", tolower(val_units), fixed=FALSE)
      val_units <- gsub("^dm$", "decimeter", tolower(val_units), fixed=FALSE)
      val_units <- gsub("^km$", "kilometer", tolower(val_units), fixed=FALSE)
      val_units <- gsub("^m$", "meter", tolower(val_units))
      val_units <- intersect(unique(val_units), possible_units) #danger in this step: discards any text that is not an expected unit

      vals<-unlist(lapply(vals, function(v){gsub(",", ".", v)} ))
      if(length(val_units)==1){
        for(v in 1:length(vals)){
          if(grepl(val_units[1], vals[v])){
            vals[v]<-gsub(val_units[1], "", vals[v])
          }
        }
        QC_units[name_unit]<-val_units[1]
        New_metadata[,colnames(New_metadata)==name_unit] <- c(vals[rownames(New_metadata)])
      }

    }
  }

  # 9. finalizing and formatting the output
  # 9.1 the units
  New_metadata_units <- c()
  if(!pre_def_units){
    # there were no pre defined units in a additional header line
    for(i in 1:ncol(New_metadata)){
      unit_i <- as.character(TermsLib[TermsLib$name %in% colnames(New_metadata)[i],]$expected_unit)
      if(length(unit_i)>0){
        New_metadata_units[colnames(New_metadata)[i]] <- unit_i
      } else{
        # if unit completely unkown, put alphanumeric
        New_metadata_units[colnames(New_metadata)[i]] <- "alphanumeric"
      }
    }
  }else{
    # there were pre defined units in a additional header line
    # user defined units obviously have priority over assumed standard units in the TermsLib file
    for(i in 1:ncol(New_metadata)){
      if(colnames(New_metadata)[i] %in% units$var_name){
        predef_unit <- as.character(units[units$var_name==colnames(New_metadata)[i],]$unit)
        if(is.na(predef_unit) | predef_unit=="" | predef_unit=="NA"){
          predef_unit <- "alphanumeric"
        }
        New_metadata_units[colnames(New_metadata)[i]] <- predef_unit
      } else{
        # if unit completely unkown, put alphanumeric
        New_metadata_units[colnames(New_metadata)[i]] <- "alphanumeric"
      }
    }
  }
  # QC_units
  if(length(QC_units)>0){
    for(u in 1:length(QC_units)){
      New_metadata_units[names(QC_units)[u]] <- QC_units[u]
    }
  }


  # 9.2 concatenate all non-MIxS terms in the misc_param term
  #MIxS_terms <- setdiff(as.character(TermsLib[TermsLib$official_MIxS==TRUE,]$name), "misc_param")
  #misc_units <- New_metadata_units[which(!colnames(New_metadata) %in% MIxS_terms)]
  #New_metadata_units <- New_metadata_units[which(colnames(New_metadata) %in% MIxS_terms)]
  #misc_metadata <- New_metadata[,!colnames(New_metadata) %in% MIxS_terms, drop=FALSE]
  #New_metadata <- New_metadata[,colnames(New_metadata) %in% MIxS_terms, drop=FALSE]
  #if(ncol(misc_metadata)>0){
  #  for(cl in 1:ncol(misc_metadata)){
  #    if(!grepl("alphanumeric", misc_units[cl])){
  #      un <- paste("(", misc_units[cl], ")", sep="")
  #    }else{
  #      un<-""
  #    }
  #    misc_metadata[,cl]<-paste(colnames(misc_metadata)[cl], ":", misc_metadata[,cl], un, sep="")
  #  }
  #  New_metadata$misc_param<-apply(misc_metadata, 1, function(x) paste(x, collapse = ";"))
  #  New_metadata_units["misc_param"] <- "alphanumeric"
  #}

  # 9.3 the section
  New_metadata_section <- c()
  for(i in 1:ncol(New_metadata)){
    if(colnames(New_metadata)[i] %in%TermsLib$name){
      New_metadata_section[colnames(New_metadata)[i]] <- as.character(TermsLib[TermsLib$name %in% colnames(New_metadata)[i],]$MIxS_section)
    } else{
      New_metadata_section[colnames(New_metadata)[i]] <- "miscellaneous"
    }
  }

  # 9.4 env_package
  if("env_package" %in% colnames(New_metadata)){
    env_package <- unique(New_metadata$env_package)
    if(length(env_package)>1){
      env_package <- "multiple_packages"
    } else if(is.na(env_package) | is.null(env_package)){
      env_package <- "not_specified"
    }
  } else{
    env_package <- "not_specified"
  }

  # 9.5 warning messages
  if(length(warningmessages)>0){
    for(i in 1:length(warningmessages)){
      warningmessages[i] <- paste(i, warningmessages[i], sep=". ")
    }
    warningmessages <- c("Please consider the following warning messages carefully before proceding:", warningmessages)
    warning(paste(warningmessages, collapse='\n'))
  }

  # 9.6 convert to the right output format (data.frame or metadata.MIxS)
  New_metadata <- new("metadata.MIxS",
                      data   = New_metadata,
                      section = New_metadata_section,
                      units      = New_metadata_units,
                      env_package = env_package,
                      type = "versatile",
                      QC = TRUE
  )
  if(!is.na(add_to)){
    New_metadata <- combine.data(New_metadata, add_to)
  }
  return(New_metadata)
}
