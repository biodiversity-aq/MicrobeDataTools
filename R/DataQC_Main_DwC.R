#==============================================================
# Author Maxime Sweetlove
# lisence CC 4.0
# Part of the POLA3R website (successor or mARS.biodiversity.aq)
# version 1.0 (2020-01-28)
# file encdong UTF-8
#
# assumptions:
#    - use NA for missing values
#
#==============================================================
# data Quality Controll (QC) for DarwinCore
#==============================================================

#' format dataframes into a DarwinCore object
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @family standardization functions
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
    EML.url<-NA
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

