#==============================================================
# Author Maxime Sweetlove
# lisence CC 4.0
# Part of the POLA3R website (successor or mARS.biodiversity.aq)
# version 1.0 (2020-01-28)
# file encdong UTF-8
#
#==============================================================
# data Quality Controll (QC) function
#==============================================================
#' find and standardize dates in a dataframe
#' @author Maxime Sweetlove CC-BY 4.0 2019
#' @family quality control functions
#' @description looks in the columns of a dataset for a column with dates and transforms them to the YYYY-MM-DD format.
#' @usage dataQC.dateCheck(dataset, date.colnames)
#' @param dataset dataframe. The dataset where the date column should be found
#' @param date.colnames character vector. a list of potential names for the column with the date. e.g. c("date", "Date", "collection date")
#' @details The date column is found based on a user-provided list of possible names to look for (data.colnames argument). If a columnname is found that corresponds to a term in the list, the dates will be convered to the YYYY-MM-DD format, if the original format can be recognized.
#' @return a list of length 2, with "$values" a vactor with the same number of rows as the dataset argument containing the corrected date values, and "$warningmessages" a vector with potential warning messages as character strings.
#' @export
dataQC.dateCheck <- function(dataset, date.colnames=c("date", "Date", "collection_date")){

  warningmessages<-c()
  NAvals <- c("NA", "ND", "unnkown", "", NA, NULL, "na")
  xdate <- intersect(date.colnames, colnames(dataset))
  if(length(xdate)>1){
    warningmessages <- multi.warnings(paste("multiple date columns found, only executed QC on", xdate[1]), warningmessages)
    xdate <- xdate[1]
  }

  if(length(xdate)==1){
    date_values <- as.character(dataset[,xdate])
    # try to format date in YYYY-MM-DD format
    for(i in 1:length(date_values)){
      if(!is.na(date_values[i]) && !(gsub("/", "", date_values[i]) %in% NAvals) && !(gsub("-", "", date_values[i]) %in% NAvals)){
        date_values[i] <- gsub("/", "-", date_values[i], fixed=TRUE)
        date_values[i] <- gsub(" ", "", date_values[i], fixed=TRUE)
        date_values[i] <- gsub(".", "-", date_values[i], fixed=TRUE)
        date_values[i] <- gsub("_", "-", date_values[i], fixed=TRUE)
        if(grepl("[A-Za-z]", date_values[i])){
          date_values[i] <- gsub("jan.+-|jan-", "01-", tolower(date_values[i]))
          date_values[i] <- gsub("feb.+-|feb-", "02-", tolower(date_values[i]))
          date_values[i] <- gsub("mar.+-|mar-", "03-", tolower(date_values[i]))
          date_values[i] <- gsub("apr.+-|apr-", "04-", tolower(date_values[i]))
          date_values[i] <- gsub("may.+-|may-", "05-", tolower(date_values[i]))
          date_values[i] <- gsub("jun.+-|jun-", "06-", tolower(date_values[i]))
          date_values[i] <- gsub("jul.+-|jul-", "07-", tolower(date_values[i]))
          date_values[i] <- gsub("aug.+-|aug-", "08-", tolower(date_values[i]))
          date_values[i] <- gsub("sep.+-|sep-", "09-", tolower(date_values[i]))
          date_values[i] <- gsub("oct.+-|oct-", "10-", tolower(date_values[i]))
          date_values[i] <- gsub("nov.+-|nov-", "11-", tolower(date_values[i]))
          date_values[i] <- gsub("dec.+-|dec-", "12-", tolower(date_values[i]))
          #how it was written before: date_values[i] <- gsub("\\sdec.+\\s|\\sdec\\s", "-12-", tolower(date_values[i]))
        }
        date_values[i] <- strsplit(date_values[i], "T")[[1]][1]
        date_split <- strsplit(date_values[i], "-")

        if(length(date_split[[1]])==3 && nchar(date_split[[1]][3])==4){ #assume DD-MM-YYYY
          year <- as.character(date_split[[1]][3])
          mnth <- as.character(sprintf("%02d", as.numeric(date_split[[1]][2])))
          day <- as.character(sprintf("%02d", as.numeric(date_split[[1]][1])))
          date_values[i]<-paste(year, mnth, day, sep="-")
        }else if(length(date_split[[1]])==3 && nchar(date_split[[1]][1])==4){ #assume YYYY-MM-DD
          year <- as.character(date_split[[1]][1])
          mnth <- as.character(sprintf("%02d", as.numeric(date_split[[1]][2])))
          day <- as.character(sprintf("%02d", as.numeric(date_split[[1]][3])))
          date_values[i]<-paste(year, mnth, day, sep="-")
        }else if(length(date_split[[1]])==2 && nchar(date_split[[1]][2])==4){ #assume MM-YYYY
          year <- as.character(date_split[[1]][2])
          mnth <- as.character(sprintf("%02d", as.numeric(date_split[[1]][1])))
          date_values[i]<-paste(year, mnth, sep="-")
        }else if(length(date_split[[1]])==2 && nchar(date_split[[1]][1])==4){ #assume YYYY-MM
          year <- as.character(date_split[[1]][1])
          mnth <- as.character(sprintf("%02d", as.numeric(date_split[[1]][2])))
          date_values[i]<-paste(year, mnth, sep="-")
        }else if(!(length(date_split)==1 && nchar(date_split)==4)){ #assume YYYY, leave value as is.
          warningmessages <- multi.warnings("check the formats of the collection dates", warningmessages)
        }
      }else{
        warningmessages <- multi.warnings("some collection dates are unknown", warningmessages)
        date_values[i] <- NA
      }
    }
  }else{
    warningmessages <- multi.warnings("no collection dates found", warningmessages)
    date_values <- rep(NA, nrow(dataset))
  }
  return(list(values=date_values, warningmessages=warningmessages))
}

#' find and standardize geographic coordinates in a dataframe
#' @author Maxime Sweetlove CC-BY 4.0 2019
#' @family quality control functions
#' @description looks in the columns of a dataset (dataframe) for a column with coordinates and transforms them a standardized decimal format (see details).
#' @usage dataQC.LatitudeLongitudeCheck(dataset, latlon.colnames)
#' @param dataset dataframe. The dataset where the date column should be found
#' @param latlon.colnames a list of length 3 with character vectors. Three vectors of potential names for the columns with the latitude-longidude values. The first vector in the list are names if latitude and longitude were to be in the same column (e.g. the MIxS lat_lon format), the second and third are for when latitude and longitude, respectively, are in seperate columns. Example: list(c("lat_lon"), c("latitude"), c("longitude"))
#' @details The date column is found based on a user-provided list of possible names to look for (latlon.colnames argument). First, a single column is searched where latitude and longitude are noted in a single field, if this returns no result, latitude and longitude are looked for in seperate fields. When found, the coordinates are transformed to decimals and returned as a single field, with values separated by a single space That is: (X Y), with X a numeric decimal latitude value and Y a numeric decimal longitude value.
#' @return a list of length 2, with "$values" a vector of same length as the number of rows in the dataset argument, and "$warningmessages" a vector with potential warning messages as character strings.
#' @export
dataQC.LatitudeLongitudeCheck <- function(dataset, latlon.colnames=list(c("lat_lon"),c("latitude"), c("longitude"))){

  warningmessages<-c()
  latlonName <- intersect(latlon.colnames[[1]], colnames(dataset))
  latName <- intersect(latlon.colnames[[2]], colnames(dataset))
  lonName <- intersect(latlon.colnames[[3]], colnames(dataset))
  lat_space_lon_output<-c()

  NAvals <- c("NA", "ND", "-", "/", "?", "unnkown", "")

  # check format lat_lon: 1 field
  if(length(latlonName)>=1){ #case 1: latitude and longitude are in the same field
    if(length(latlonName)>1){ #more than one column detected
      latlonName <- latlonName[1]
      warningmessages <- multi.warnings("multiple lat_lon columns found, just used the first one", warningmessages)
    }
    latlon_values <- dataset[,latlonName]
    for(i in 1:length(latlon_values)){
      if(!is.na(latlon_values[i]) && ! latlon_values[i] %in% NAvals){
        #try to find the right separator
        #first look for a tab
        latlon_val <- gsub('\t', " ", latlon_values[i])
        #ten split on for spaces
        latlon_val <- gsub("^\\s+|\\s+$", "", latlon_val) #remove leading and trailing spaces
        latlon_val <- strsplit(latlon_val," ")[[1]] #split on space
        if(length(latlon_val)==1){ #it's not separated by a space or a tab
          if(grepl(",", latlon_val)){ #comma can be decimal separator
            if(lengths(regmatches(latlon_val, gregexpr(",", latlon_val)))>1){
              # more than one comma: assume it is a decimal separator
              latlon_val <- gsub(",", ".", latlon_val)
            }else{
              #if there is just one comma, assume it is the value separator
              latlon_val <- gsub(",", " ", latlon_val)
            }
          }
          latlon_val <- gsub(",", " ", latlon_val)
          latlon_val <- gsub(";", " ", latlon_val)
          latlon_val <- gsub(":", " ", latlon_val)
          latlon_val <- gsub("|", " ", latlon_val)
          latlon_val <- gsub("\\s+", " ", latlon_val) #collapse multiple spaces
          latlon_val <- gsub("^\\s+|\\s+$", "", latlon_val) #remove leading and trailing spaces
          latlon_val <- strsplit(latlon_val," ")[[1]]
        }
        #weed out the NAs
        if(length(latlon_val)==1){#is the lat_lon is still not split, it must be NA
          if(is.na(latlon_val) | latlon_val %in% NAvals){
            lat_space_lon_output<-c(lat_space_lon_output, NA)
            warningmessages <- multi.warnings("there are unknowns in the coordinates", warningmessages)
          }else{
            warningmessages <- multi.warnings("for some coordinates, the format could not be recognized", warningmessages)
            lat_space_lon_output<-c(lat_space_lon_output, latlon_val)
          }
        }else{
          #re-formatting the lat_lon values
          if(length(latlon_val)==2){
            lat<-latlon_val[1]
            lon<-latlon_val[2]
          }else if(length(latlon_val)==4){
            lat<-paste(latlon_val[1], latlon_val[2], sep="")
            lon<-paste(latlon_val[3], latlon_val[4], sep="")
          }

          #converto to decimal format
          lat <- coordinate.to.decimal(lat)
          lon <- coordinate.to.decimal(lon)
          lat_space_lon_output<-c(lat_space_lon_output, paste(lat, lon))
        }
      }else{
        lat_space_lon_output<-c(lat_space_lon_output, NA)
        warningmessages <- multi.warnings("there are unknowns in the coordinates", warningmessages)
      }
    }
    # check format lat_lon: 2 fields
  }else if(length(latName)>=1 & length(lonName)>=1){ # case2: latitude and longitude are be in seperate fields
    if(length(latName)>1 | length(lonName)>1){
      latName <- latName[1]
      lonName <- lonName[1]
      warningmessages <- multi.warnings("multiple latitude or longitude columns found, just used the first ones", warningmessages)
    }
    lat_values <- dataset[,latName]
    lon_values <- dataset[,lonName]
    if(length(lat_values)==length(lon_values)){
      for(i in 1:length(lat_values)){
        lat <- lat_values[i]
        lon <- lon_values[i]

        if(is.na(lat) | lat %in% NAvals | is.na(lon) | lon %in% NAvals){
          warningmessages <- multi.warnings("there are NAs in the coordinates", warningmessages)
          lat_space_lon_output<-c(lat_space_lon_output, NA)
        }else{
          lat <- coordinate.to.decimal(lat)
          lon <- coordinate.to.decimal(lon)
          lat_space_lon_output<-c(lat_space_lon_output, paste(lat, lon))
        }
      }
    }else{
      warningmessages <- multi.warnings("the length of the latitude and longitude columns do not match, process not executed", warningmessages)
      lat_space_lon_output<-c(lat_space_lon_output, NA)
    }
  }else{ #case 3: no latitude and longitude data found
    warningmessages <- multi.warnings("no latitude longitudes found", warningmessages)
    lat_space_lon_output<-NA
  }
  return(list(values=lat_space_lon_output, warningmessages=warningmessages))
}

#' make an educated guess to the MIxS environmental package of a dataset
#' @author Maxime Sweetlove CC-BY 4.0 2019
#' @family quality control functions
#' @description looks in the columns of a dataset for clues to what the most appropriate MIxS environmental could be for each sample.
#' @usage dataQC.guess.env_package.from.data(dataset, pckge.colnames=c("env_package", "ScientificName"))
#' @param dataset dataframe. The dataset for which the MIxS environmental package should be found or guessed
#' @param pckge.colnames a character vector. A vector with the potential names for the column where the environmental package can be found. Place the terms in order of decreasing likeliness.
#' @details The "Minimul Information on any (x) Sequence" (MIxS) standard requires an appropriate environmental package with MIxS terms to be selected to document the data. Some data and nucleotide archives enforce their users to select such a package. This function is made to automatically either find the package in a dataset, of guess it based on the data that is present.
#' @return a list of length 2, with "$values" a vector of same length as the number of rows in the dataset argument, and "$warningmessages" a vector with potential warning messages as character strings.
#' @export
dataQC.guess.env_package.from.data <- function(dataset, pckge.colnames=c("env_package", "ScientificName")){

  warningmessages <- c()
  warningmessages<-multi.warnings("the env_package was not specified. An educated guess was made, but should be checked", warningmessages)

  #look if a collumn can be found
  pck <- intersect(pckge.colnames, colnames(dataset))
  if(length(pck)>0){
    env_package<-dataset[,pck[1]]
  } else{
    env_package<-NULL
    warningmessages<-multi.warnings("No env_package could be infered", warningmessages)
  }

  # try to convert the input to a correct package name (necessary when the package is not given directly and needs to be infered)
  if(!is.null(env_package)){
    env_package <- as.character(env_package)
    for(pk in 1:length(env_package)){
      pk_val <- tolower(env_package[pk])
      if(grepl("water", pk_val)|grepl("sea", pk_val)|grepl("ocean", pk_val)|
         grepl("pond", pk_val)|grepl("river", pk_val)|grepl("lake", pk_val)|
         grepl("aquatic", pk_val)|grepl("lagustrine", pk_val)|grepl("marine", pk_val)){
        env_package[pk]<-"water"
      }else if(grepl("soil", pk_val)|grepl("earth", pk_val)|grepl("sand", pk_val)|
               grepl("loam", pk_val)|grepl("clay", pk_val)|grepl("silt", pk_val)|
               grepl("peat", pk_val)|grepl("chalk", pk_val)|grepl(".[a-z]sol", pk_val)){
        env_package[pk]<-"soil"
      }else if(grepl("built_environment", pk_val)|grepl("built-environment", pk_val)|
               grepl("building", pk_val)|grepl("concrete", pk_val)|grepl("brick", pk_val)|
               grepl("cement", pk_val)|grepl("pavement", pk_val)|grepl("house", pk_val)|
               grepl("lobby", pk_val)|grepl("room", pk_val)){
        env_package[pk]<-"built_environment"
      }else if(grepl("air", pk_val)|grepl("wind", pk_val)){
        env_package[pk]<-"air"
      }else if(grepl("sediment", pk_val)|grepl("floor", pk_val)|grepl("bottom", pk_val)){
        env_package[pk]<-"sediment"
      }else if(grepl("microbial_mat_biofilm", pk_val)|grepl("microbial mat/biofilm", pk_val)|
               grepl("biofilm", pk_val)|grepl("microbial", pk_val)|grepl("mat", pk_val)){
        env_package[pk]<-"microbial_mat_biofilm"
      }else if(grepl("human_associated", pk_val)|grepl("human-associated", pk_val)){
        env_package[pk]<-"human_associated"
      }else if(grepl("human_gut", pk_val)|grepl("human-gut", pk_val)){
        env_package[pk]<-"human_gut"
      }else if(grepl("human_oral", pk_val)|grepl("human-oral", pk_val)){
        env_package[pk]<-"human_oral"
      }else if(grepl("human_skin", pk_val)|grepl("human-skin", pk_val)){
        env_package[pk]<-"human_skin"
      }else if(grepl("human_vaginal", pk_val)|grepl("human-vaginal", pk_val)){
        env_package[pk]<-"human_vaginal"
      }else if(grepl("wastewater_sludge", pk_val)|grepl("wastewater/sludge", pk_val)){
        env_package[pk]<-"wastewater_sludge"
      }else if(grepl("plant_associated", pk_val)|grepl("plant-associated", pk_val)){
        env_package[pk]<-"plant_associated"
      }else if(grepl("host_associated", pk_val)|grepl("host-associated", pk_val)|
               grepl("host", pk_val)|grepl("tissue", pk_val)|grepl("gut", pk_val)|
               grepl("skin", pk_val)|grepl("stomach", pk_val)|grepl("mouth", pk_val)|
               grepl("anus", pk_val)|grepl("ear", pk_val)|grepl("vagina", pk_val)|
               grepl("lungs", pk_val)|grepl("genital", pk_val)|grepl("penis", pk_val)|
               grepl("saliva", pk_val)|grepl("urine", pk_val)|grepl("feaces", pk_val)|
               grepl("oesophagus", pk_val)){
        env_package[pk]<-"host_associated"
      }else if(grepl("miscellaneous_natural_or_artificial_environment", pk_val)|
               grepl("miscellaneous natural or artificial environment", pk_val)|
               grepl("miscellaneous", pk_val)){
        env_package[pk]<-"miscellaneous_natural_or_artificial_environment"
      }else{
        env_package[pk]<-"miscellaneous_natural_or_artificial_environment"
        warningmessages<-multi.warnings("for some samples no env_package could be infered", warningmessages)
      }
    }
  }

  return(list(values=env_package, warningmessages=warningmessages))
}

#' make variable names in a dataset to comply to a standrad
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @family quality control functions
#' @description checks a set of terms (e.g. columnnames) to a Standard, and flags inconsistencies, gives solutions if possible.
#' @usage dataQC.TermsCheck(observed=NA, exp.standard="MIxS", exp.section=NA, fuzzy.match=TRUE,out.type="full")
#' @param observed character vector. The terms to be checked
#' @param exp.standard character. The expected standard to which he terms should comply. Either MIxS (Minimum Information on any x sequence), DwC (DarwinCore), INSDC (International Nucleotide Sequence Database Consortium).
#' @param exp.section character. Optionally an specific section standard where the terms should come from. When exp.standard is MIxS, the allowed sections are: core,  air, built_environment, host_associated, human_associated, human_gut, human_oral, human_skin, human_vaginal, microbial_mat_biofilm, miscellaneous_natural_or_artificial_environment, plant_associated, sediment, soil, wastewater_sludge, water. When exp.standard is DwC, the allowed sections are: event, occurence, emof
#' @param fuzzy.match logical. If TRUE, fuzzy matching will be done when no corresponding term is found in the Standard
#' @param out.type character. The type of the output. Either "full" (the output is a list of three lists: terms_OK=the correct terms, terms_wrongWithSolution = wrong terms with a proposed solution, and terms_notFound = terms that had no match), "logical" (output is logical vector for exact matches or not) or "best_match" (output returns a vector with the best matching terms). default full
#' @details For interoperability of data and data archiving, variable names in datasets need to comply to vocabulary standards. This function compares a list of existing terms to the MIxS or DwC standard, and tries to fiend the best matches.
#' @return depending on out.type, either a boolean, or a list of length 3, with "$terms_OK" (terms that comply to the standard), "$terms_wrongWithSolution" (terms that do not comply to the standard but have a close match), and "$terms_notFound" (terms that do not comply to the standard, and that not match any term in it)
#' @export
dataQC.TermsCheck <- function(observed=NA, exp.standard="MIxS", exp.section=NA, fuzzy.match=TRUE,
                        out.type="full"){
  #requires RecordLinkage
  MIxS_sections <- c("core",  "air", "built_environment", "host_associated",
                     "human_associated", "human_gut", "human_oral", "human_skin",
                     "human_vaginal", "microbial_mat_biofilm", "sediment", "soil",
                     "miscellaneous_natural_or_artificial_environment", "plant_associated",
                     "wastewater_sludge", "water")
  DwC_sections <- c("event", "occurrence", "emof")

  #chage terms to lowercase: lowers number of possible differences
  names(observed) <- observed
  observed <- tolower(observed)

  # 1. check input
  if(!tolower(out.type) %in% c("full", "logical", "best_match")){
    stop("out.type should be either:\n\t \"full\": output will be a list of three lists, with:\n\t\t1) the correct terms\n\t\t2) wrong terms with a solution\n\t\t3) terms that had no match\n\t\"logical\": the output will be a TRUE/FALSE vector\n\t\"best_match\": the output will be the best matching term")
  }

  if(!tolower(exp.standard) %in% c("mixs", "dwc", "insdc")){
    stop("expected standard should be either \"MIxS\", \"DwC\" or \"INSDC\".")
  }
  if(tolower(exp.standard) == "mixs"){
    exp.standard<-"MIxS"
    if(!(is.na(exp.section)) && !tolower(exp.section) %in% MIxS_sections || length(exp.section)>1){
      stop(paste("expected section for MIxS should be NA or ONE of the following:\n\t", paste(MIxS_sections, collapse="\n\t"), sep=""))
    }
  }
  if(tolower(exp.standard) == "dwc"){
    exp.standard<-"DwC"
    if(!(is.na(exp.section)) && !tolower(exp.section) %in% DwC_sections || length(exp.section)>1){
      stop(paste("expected section for DwC should be NA or ONE of the following:\n\t", paste(DwC_sections, collapse="\n\t"), sep=""))
    }
  }
  if(tolower(exp.standard) == "insdc"){
    exp.standard<-"INSDC"
    exp.section <- NA
  }

  # 2. get the right synonym list
  if(is.na(exp.section)){
    synonymList <- sapply(as.character(TermsLib[TermsLib$name_origin %in% exp.standard,]$synonyms), function(x){strsplit(x, ";")})
    names(synonymList) <- TermsLib[TermsLib$name_origin %in% exp.standard,]$name
  }else{
    exp.section <- gsub("^core$", "MIxS_core", exp.section, fixed=FALSE)
    exp.section <- gsub("^event$", "DwC_Event", exp.section, fixed=FALSE)
    exp.section <- gsub("^occurrence$", "DwC_Occurrence", exp.section, fixed=FALSE)
    exp.section <- gsub("^emof$", "DwC_eMoF", exp.section, fixed=FALSE)
    synonymList <- sapply(as.character(TermsLib[TermsLib[,exp.section] > 0,]$synonyms), function(x){strsplit(x, ";")})
    names(synonymList) <- TermsLib[TermsLib[,exp.section] > 0,]$name
  }

  # 3. check terms
  terms_OK <- c()
  terms_wrongWithSolution <- c()
  terms_notFound <- c()
  for(obs_term in observed){
    obs_termName <- names(observed[grepl(paste("\\b", obs_term, "\\b", sep=""), observed, fixed=FALSE)])
    # first look for complete matches
    term_match <- synonymList[grepl(paste("\\b", obs_term, "\\b", sep=""), synonymList, fixed=FALSE)]
    if(length(term_match)==1){
      # one complete match
      if(obs_termName==names(term_match)){
        # terms is completely OK
        terms_OK <- c(terms_OK, obs_termName)
      }else{
        # term is one of the synonyms
        solutionNames <- names(terms_wrongWithSolution) #the names are the obsrved terms
        terms_wrongWithSolution <- c(terms_wrongWithSolution, names(term_match)) #the values are the matching term from the standard
        names(terms_wrongWithSolution) <- c(solutionNames, obs_termName)
      }
    }else{
      # no match found: check for spelling errors with a fuzzy match
      if(fuzzy.match){
        distances <- RecordLinkage::jarowinkler(obs_term, unlist(synonymList)) # caluculate all word distances
        if(max(distances)>0.85){ #treshold step, not bothering with worse matches
          dist_mtch <- match(max(distances),distances) #get index of best match
          if(length(dist_mtch)==1){ #only one match can be best, otherwise need human input
            dist_mtch <- unlist(synonymList)[dist_mtch] #get term that matched best
            term_match <- synonymList[grepl(paste("\\b", dist_mtch, "\\b", sep=""), synonymList, fixed=FALSE)]
            # return best solution
            solutionNames <- names(terms_wrongWithSolution)
            terms_wrongWithSolution <- c(terms_wrongWithSolution, names(term_match))
            names(terms_wrongWithSolution) <- c(solutionNames, obs_termName)
          }else{
            terms_notFound <- c(terms_notFound, obs_termName)
          }
        }else{
          terms_notFound <- c(terms_notFound, obs_termName)
        }
      }else{
        terms_notFound <- c(terms_notFound, obs_termName)
      }
    }
  }

  if(out.type=="full"){
    out<-list(terms_OK = terms_OK,
              terms_wrongWithSolution = terms_wrongWithSolution,
              terms_notFound = terms_notFound)
  }else if(out.type=="logical"){
    out <- sapply(observed, function(w){
      if(w %in% terms_OK){w<-TRUE}else{w<-FALSE}
    })
  }else if(out.type=="best_match"){
    out <- sapply(observed, function(w){
      if(w %in% terms_wrongWithSolution){w<-names(terms_wrongWithSolution[terms_wrongWithSolution %in% w])
      }else if(w %in% terms_notFound){w<-NA
      }else{w}})
  }


  return(out)
}

#' generate a footprintWKT from coordinates
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @family quality control functions
#' @description generate a footprintWKT column for a dataset (preferably DarwinCore eventCore). Lowest level events must be point locations.
#' @param dataset data.frame. A data.frame with Event data, formatted as a DarwinCore EventCore file. Must include decimalLatitude, decimalLongitude, eventID and parentEventID.
#' @param NA.val character. What to fill in when there is no data. Default NA
#' @usage dataQC.generate.footprintWKT(dataset, NA.val=NA)
#' @details This function will format the geographic coordinates of an event (in the fields decimalLatitude and decimalLongitude) into the WellKnownText format. For events that are at the lowest level of the hierarchy, it assumes a single point location. Higher level events are desribed by polygons based on the coordinates of the child events.
#' @return a vector with a footprintWKT value for each row in the dataset
#' @export
dataQC.generate.footprintWKT <- function(dataset, NA.val=NA){
  #requires Orcs
  #format POINT(lon lat)
  if(!"decimalLatitude" %in% colnames(dataset) |
     !"decimalLongitude" %in% colnames(dataset) |
     !"eventID" %in% colnames(dataset)){
    stop("The data must include the following columns:\n\tdecimalLatitude, decimalLongitude, eventID")
  }
  footprintWKT_vec <- c()
  if("parentEventID" %in% colnames(dataset)){
    for(evi in 1:nrow(dataset)){
      if(!dataset[evi,]$eventID %in% dataset$parentEventID){
        # in this case the event is at the lowest level, assume it is a point
        lat<-dataset[evi,]$decimalLatitude
        lon<-dataset[evi,]$decimalLongitude
        if(!is.na(lat) & !is.na(lon) & lat!="" & lon!=""){
          footprintWKT_vec <- c(footprintWKT_vec, paste("POINT(", as.character(lon), " ", as.character(lat), ")", sep=""))
        } else{
          footprintWKT_vec <- c(footprintWKT_vec, NA.val)
        }
      }else{
        # this event is the parent of one or more child events. we need to find all the coordinates of the lowest events under this parent
        child_ev <- dataset[dataset$parentEventID == dataset[evi,]$eventID,]$eventID #the direct descendants
        doneList<-c()
        n_added <- length(child_ev)
        while(! n_added == 0){ #loop through the events to add all the other descendants
          n_added = 0
          for(d in child_ev){
            if(!d %in% doneList){
              d2 <- dataset[dataset$parentEventID == d,]$eventID
              doneList <- c(doneList,d)
              if(length(d2)>0){
                child_ev<- c(child_ev, d2)
                n_added <- n_added +1
              }
            }
          }
        }
        child_ev <- setdiff(child_ev, unique(dataset$parentEventID))

        lat<-dataset[dataset$eventID %in% child_ev,]$decimalLatitude
        lon<-dataset[dataset$eventID %in% child_ev,]$decimalLongitude
        if(!all(is.na(lat)) & !all(is.na(lon)) & all(lat!="") & all(lon!="")){
          coords <- as.matrix(data.frame(x=as.numeric(lon), y=as.numeric(lat)))
          polygx <- Orcs::coords2Polygons(coords, hole=FALSE, ID=ev)
          footprintWKT_vec <- c(footprintWKT_vec, writeWKT(polygx))
        } else{
          footprintWKT_vec <- c(footprintWKT_vec, NA.val)
        }
      }
    }
  }else{
    footprintWKT_vec <- paste("POINT(", as.character(dataset$decimalLongitude),
                              " ", as.character(dataset$decimalLatitude), ")",
                              sep="")
  }

  footprintWKT_vec <- gsub("POINT(NA NA)", NA.val, footprintWKT_vec, fixed=TRUE)
  footprintWKT_vec <- gsub("POINT( )", NA.val, footprintWKT_vec, fixed=TRUE)
  return(footprintWKT_vec)

}

#' find the samplenames in a dataset
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @family quality control functions
#' @description find the sample names in a given (meta-)dataset where at least an attempt has been made to standardize the data following MIxS or DarwinCore
#' @param dataset data.frame. The data.frame where to look for the sample names
#' @param ask.input logical. If TRUE, console input will be requested to the user when a problem occurs. Default TRUE
#' @usage dataQC.findNames(dataset, ask.input=TRUE)
#' @return a list of length 3 with: "$Names" a named vector with the most likely sample names, "$Names.column" the column name where the sample names were found, "$warningmessages" a vector with warning messages
#' @details It is often not clear where the sample names are in a dataset. This function makes an educated guess, based on rownames or tags that are often used to indicate sample names. If ask.input, then the process happens user-supervised.
#' @export
dataQC.findNames <- function(dataset = NA, ask.input=TRUE){
  orig_names <- data.frame(matrix(data=NA, ncol=5, nrow=nrow(dataset)))
  colnames(orig_names)<-c("original_names", "eventID", "parentEventID", "occurrenceID", "INSDC_SampleID")
  rownames(orig_names) <- rownames(dataset)
  warningmessages <- ""
  item_shared <- intersect(TermsSyn["original_name"][[1]], tolower(colnames(dataset)))
  # loop over the potential names in a fixed order
  if(length(item_shared)>=1){
    if(ask.input){
      for(si in 1:length(item_shared)){
        likely_sampNames <- item_shared[si] #order of item_shared is in decreasing likelyness
        likely_sampNames_orig <- colnames(dataset)[tolower(colnames(dataset))==likely_sampNames]
        message(paste("The original sample names could be in the \"",likely_sampNames_orig,"\" column.\n", 
                      "\tThe first five names in this column are:\n\t\t", 
                      paste(dataset[1:5,tolower(colnames(dataset))==likely_sampNames], collapse="\n\t\t"), 
                      " ...\n\tDoes this seems correct? (y/n)\n",
                      sep=""))
        doNext <- readline()
        if(tolower(doNext) %in% c("y", "yes")){
          orig_names$original_names <- as.character(dataset[,tolower(colnames(dataset))==likely_sampNames])
          warningmessages <- paste("assumed the \"",likely_sampNames_orig,"\" column contained the original sample names", sep="")
          break
        }else if(!tolower(doNext) %in% c("n", "no")){
          stop("incorrect input... only yes or no allowed.")
        }
      }
    }else{
      likely_sampNames <- item_shared[1]
      likely_sampNames_orig <- colnames(dataset)[tolower(colnames(dataset))==likely_sampNames]
      orig_names$original_names <- as.character(dataset[,tolower(colnames(dataset))==likely_sampNames])
      warningmessages <- paste("assumed the \"",likely_sampNames_orig,"\" column contained the original sample names", sep="")
    }
  }else{
    likely_sampNames_orig <- NA
    if(.row_names_info(dataset)>0){#rownames provided by the user might be the original names
      if(ask.input){
        message(paste("No original sample names found...\n\tthe rownames can be used as an alternative.\n", 
                      "\tThe first five names are:\n\t\t", 
                      paste(row.names(dataset)[1:5], collapse="\n\t\t"), 
                      " ...\n\tuse rownames as sample names? (y/n)\n",
                      sep=""))
        doNext <- readline()
        if(tolower(doNext) %in% c("y", "yes")){
          orig_names$original_names <- row.names(dataset)
          warningmessages <- "no original sample names found, used the rownames instead"
        }else if(tolower(doNext) %in% c("n","no")){
          message("you chose no.\n\t Can you give the columnname with the original sample names instead? Leave blank and hit enter to ignore\n")
          doNext2 <- readline()
          if(doNext2!=""){
            if (doNext2 %in% colnames(dataset)){
              orig_names$original_name <- dataset[,doNext2]
            }else{
              stop(paste("could not find the column", doNext2, "in the colnames of the dataset provided..."))
            }
          }
        } else{
          stop("incorrect input... only yes or no allowed.")
        }
      }
    }else{
      warningmessages <- "no original sample names found"
    }
  }

  #find and return some other common sample IDs
  if("eventID" %in% colnames(dataset)){
    orig_names$eventID <- as.character(dataset[,colnames(dataset)=="eventID"])
  }
  if("parentEventID" %in% colnames(dataset)){
    orig_names$parentEventID <- as.character(dataset[,colnames(dataset)=="parentEventID"])
  }
  if("INSDC_SampleID" %in% colnames(dataset)){
    orig_names$INSDC_SampleID <- as.character(dataset[,colnames(dataset)=="INSDC_SampleID"])
  }
  if("occurrenceID" %in% colnames(dataset)){
    orig_names$occurrenceID <- as.character(dataset[,colnames(dataset)=="occurrenceID"])
  }

  return(list(Names=orig_names, Names.column=likely_sampNames_orig, warningmessages = warningmessages))
}

#' extract the finest resolution taxonomy from a dataset
#' @author Maxime Sweetlove CC-BY 4.0 2020
#' @family quality control functions
#' @description tries to find taxonomic names for samples (rows) in a dataset (data.frame),
#' @usage find.sampleTaxon(dataset)
#' @param dataset a data.frame. The dataset with samples as rows, and taxonomy information in the columns. using the MIxS or DarwinCore taxonomy terms, the taxonomy information will be extracted
#' @return a vector with the highest level taxonomic name found, with genus and species epithet separated by a space.
#' @export
dataQC.TaxonListFromData <- function(dataset){
  taxaTerms<-c("specificEpithet", "subgenus", "genus", "family",
               "order","class","phylum","kingdom", "domain")
  rankNames<-intersect(colnames(dataset), taxaTerms)
  # 1. try subspecf_gen_lin (for MIxS data)
  if("subspecf_gen_lin" %in% colnames(dataset)){
    taxaNames<-dataset$subspecf_gen_lin
  }else if("scientificName" %in% colnames(dataset)){
    # 2. try scientificName (for DwC data)
    taxaNames<-dataset$scientificName
  }else{
    # 3. look for the taxonomy columns and try to get a name frome there
    if("genus" %in% rankNames){
      if("specificEpithet" %in% rankNames){
        taxaNames <- paste(dataset$genus, dataset$specificEpithet, sep=" ")
      }else{
        taxaNames <- DwC.data$genus
      }
    }else{
      taxaNames <- rep(NA, nrow(dataset))
    }
  }
  #standardize all unknown values to be NA
  taxaNames <- as.character(taxaNames)
  taxaNames[taxaNames=="NA"]<-NA
  taxaNames[taxaNames==""]<-NA

  for(tx in 1:length(taxaNames)){
    if(is.na(taxaNames[tx])){
      # means both genus and species names were unknown
      tx_i=1
      while(tx_i<(length(rankNames)+1)){
        tx_term <- rankNames[tx_i]
        taxName <- dataset[tx,tx_term]
        if(!is.na(taxName) & !(taxName %in% c("", "NA"))){
          taxaNames[tx] <- taxName
          tx_i<-length(rankNames)+100
        }
        tx_i <- tx_i+1
      }
    }
  }

  taxaNames <- gsub("\\s$", "", taxaNames, fixed=FALSE) #remove trailing spaces
  return(taxaNames)
}

#' perform a quality control onf taxonomic names
#' @author Maxime Sweetlove CC-BY 4.0 2020
#' @family quality control functions
#' @description checks a list of taxonomic names, and perform a basic quality control
#' @usage dataQC.taxaNames(taxaNames)
#' @param taxaNames character vector. a list of scientific taxonomic names
#' @details looks for trailing spaces or common typos in taxonomic names.
#' @return a vector with the checked taxon names
#' @export
dataQC.taxaNames <- function(taxaNames){
  taxaNamesQC <- data.frame(matrix(data=NA, ncol=3, nrow=length(taxaNames)), stringsAsFactors = FALSE)
  colnames(taxaNamesQC) <- c("speciesLevelName", "scientificName", "identificationQualifier")

  for(tx in 1:length(taxaNames)){
    taxon <- taxaNames[tx]
    taxon <- gsub("\\s$", "", taxon, fixed=FALSE) #remove trailing spaces
    taxon <- gsub("Eukaryotes", "Eukaryota", taxon, fixed=TRUE)
    taxon <- gsub("Eukarya", "Eukaryota", taxon, fixed=TRUE)
    taxon_split <- strsplit(taxon, " ")[[1]]

    if(length(taxon_split)>1){
      taxaNamesQC[tx,]$speciesLevelName <- taxon
      if(taxon_split[2]=="cf."){
        taxaNamesQC[tx,]$scientificName <- tx_split[1]
        taxaNamesQC[tx,]$identificationQualifier <- paste(tx_split[2], tx_split[3])
      }else if(taxon_split[2]=="sp." | taxon_split[2]=="sp"){
        taxaNamesQC[tx,]$scientificName <- tx_split[1]
        taxaNamesQC[tx,]$identificationQualifier <- "sp."
      }else if(taxon_split[2]=="spp." | taxon_split[2]=="spp"){
        taxaNamesQC[tx,]$scientificName <- tx_split[1]
        taxaNamesQC[tx,]$identificationQualifier <- "spp."
      }else if(taxon_split[2]=="gr."){
        taxaNamesQC[tx,]$scientificName <- tx_split[1]
        taxaNamesQC[tx,]$identificationQualifier <- paste(taxon_split[3], "group")
      }else{
        taxaNamesQC[tx,]$scientificName <- taxon
      }
    }else{
      taxaNamesQC[tx,]$speciesLevelName <- paste(taxon, "sp.")
      taxaNamesQC[tx,]$scientificName <- taxon
      taxaNamesQC[tx,]$identificationQualifier <- "sp."
    }
  }
  return(taxaNamesQC)
}

#' Complete a list of taxonomic names with information form an authorotive taxonomic backbone
#' @author Maxime Sweetlove CC-BY 4.0 2020
#' @family quality control functions
#' @description complete a list of taxonomic names by looking-up missing information on an accepted taxonomic registery
#' @param taxaNames a character vector. A list with the taxonomic names to look for
#' @param taxBackbone a character string. The taxonomic backbone to querry. Either "worms" or "gbif"
#' @usage dataQC.completeTaxaNamesFromRegistery(taxaNames, taxBackbone="gbif")
#' @details using the API-client connection to the World Registry of Marine Species (WORMS), additional taxonomic information can be added to an existing list of taxa
#' @return a dataframe with the following fields:scientificName, scientificNameID, aphID, kingdom, phylum, class, order, family, genus, specificEpithet, scientificNameAuthorship, namePublishedInYear
#' @export
dataQC.completeTaxaNamesFromRegistery <- function(taxaNames, taxBackbone="worms"){
  #requires worrms and rgbif

  taxid_key <- data.frame(scientificName = unique(taxaNames), scientificNameID=NA,
                          aphID=NA, kingdom=NA, phylum=NA, class=NA, order=NA, family=NA,
                          genus=NA, specificEpithet=NA, scientificNameAuthorship=NA,
                          namePublishedInYear=NA)

  if(tolower(taxBackbone) == "worms"){
    for(nc in 1:nrow(taxid_key)){
      taxon<-as.character(taxid_key[nc,]$scientificName)
      if(!taxon %in% c("", "NA", NA)){
        taxid <- tryCatch({
          tx <- worrms::wm_name2id(taxon)
        }, error = function(e){
          tx <- NA
        }
        )

        if(!is.na(taxid)){
          taxnum <- taxid
          taxid<-paste("urn:lsid:marinespecies.org:taxname:", taxid, sep="")
          taxdata <- data.frame(worrms::wm_record(taxnum))

          #fill the taxid_key table
          taxid_key[nc,]$aphID <- taxnum
          taxid_key[nc,]$scientificNameID <- taxid
          taxid_key[nc,]$kingdom<-taxdata$kingdom
          taxid_key[nc,]$phylum<-taxdata$phylum
          taxid_key[nc,]$class<-taxdata$class
          taxid_key[nc,]$order<-taxdata$order
          taxid_key[nc,]$family<-taxdata$family
          taxid_key[nc,]$genus<-taxdata$genus
          if(!is.na(taxdata$scientificname)){
            taxid_key[nc,]$specificEpithet<-strsplit(taxdata$scientificname, " ")[[1]][2]
          }else{
            taxid_key[nc,]$specificEpithet <- NA
          }

          if(!is.na(taxdata$authority)){
            authoryear<-strsplit(taxdata$authority, ", ")[[1]]
            taxid_key[nc,]$scientificNameAuthorship<-gsub("\\(", "", authoryear[1], fixed=FALSE)
            taxid_key[nc,]$namePublishedInYear<-gsub("\\)", "", authoryear[2], fixed=FALSE)
          }else{
            taxid_key[nc,]$scientificNameAuthorship <- NA
            taxid_key[nc,]$namePublishedInYear <- NA
          }

        }else{
          taxid_key[nc,]$scientificNameID <- taxid
          taxid_key[nc,]$aphID <- taxid
        }
      }
    }
  }else if(tolower(taxBackbone) == "gbif"){
    for(nc in 1:nrow(taxid_key)){
      taxon<-as.character(taxid_key[nc,]$scientificName)
      if(!taxon %in% c("", "NA", NA)){
        taxid <- tryCatch({
          tx <- rgbif::name_backbone(name="taxon")
        }, error = function(e){
          tx <- NA
        }
        )
      }
      ### to be finished
    }

  }else{
    stop("invalid taxBackbone argument")
  }

  return(taxid_key)
}


#' check a dataset for an event-structure
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @family quality control functions
#' @description checks if an eventID is present in a (QC'd) dataset, and generates one if not. Does the same of parentEventID if check.parentEventID is TRUE, and check the hierarchical relationships between eventID and parentEventID
#' @param dataset data.frame. The dataset for which the event structure should be checked
#' @param eventID.col character. The column where the names of the events are given. Default eventID. If NA, and event.prefix must be provided to be able to create unique eventIDs
#' @param parentEventID.col character. The column where the names of the parentEvents are given. If NA, parentEvents are not considered. Default NA
#' @param project.col character. The column where the names of the projects are given. Projects are high-level parentEvents that group a large number of events. If NA, projects are not considered. If there is just one project, please consider the project parameter. Default NA. This parameter overrides the project parameter. Only effective if complete.hierarchy is TRUE.
#' @param project character. The project that groups all samples. This parameter is overriden by the the project.col parameter if not NA. Only effective if complete.hierarchy is TRUE.
#' @param event.prefix character. A prefix to feauture in the eventIDs if these have to be created from scratch. This parameter overrides the eventID.col parameter.
#' @param complete.hierarchy logical. If TRUE, the event-parentEvent hierarchy structure will be completed upt to a root (project). Any parentEvents that are not listed among the eventIDs will also be created. Default FALSE. if parentEventID.col was NA, new parentEventIDs will be created.
#' @usage dataQC.eventStructure(dataset, eventID.col = "eventID", parentEventID.col = NA, project.col = NA, project = NA, event.prefix = NA, complete.hierarchy=FALSE)
#' @details An event structure is a hierarchical grouping of (sub-) samples (events, that is, something that occurs at a place and time) into higher parentEvents, like expeditions, projects,... This interlinked structure is for instance used in the DarwinCore Event standard. The algorithm here looks for specific columnames that may give an indication to the event structure (e.g. a column with sampe names, projects,...). If disered, the user can complete the event structure by rooting it into an over-arching project.
#' @return a dataframe with 2 columns: eventID and parentEventID
#' @export
dataQC.eventStructure <- function(dataset, eventID.col = "eventID", parentEventID.col = NA,
                                  project.col = NA, project = NA, event.prefix = NA,
                                  complete.hierarchy=FALSE){

  #requires stringr
  eventTable <- data.frame(original_name = rownames(dataset), eventID=rep(NA, nrow(dataset)),
                           parentEventID=rep(NA, nrow(dataset)), type=NA,
                           stringsAsFactors = FALSE)
  rownames(eventTable)<-rownames(dataset)

  # 0. checking the input and error handling
  if(complete.hierarchy & is.na(parentEventID.col) & is.na(project) & is.na(project.col)){
    # complete.hierarchy makes no sense if there is no project or parentEvent information
    complete.hierarchy<-FALSE
  }

  if(is.na(eventID.col)){
    if(is.na(event.prefix)){
      stop("eventID.col or an event.prefix must be provided")
    }else{
      eventID.col <- "eventID"
      dataset$eventID <- paste(event.prefix, ":N",
                               stringr::str_pad(1:nrow(dataset), nchar(nrow(dataset)), pad = "0"), sep="")
    }
  }else if(!eventID.col %in% colnames(dataset)){
    stop("eventID.col not found")
  }else if(length(unique(dataset[,eventID.col])) < nrow(dataset)){
    stop("all events must be unique")
  }

  if(!is.na(project.col)){
    if(!project.col %in% colnames(dataset)){
      stop("project.col not found")
    }else{
      eventTable$project <- as.character(dataset[,colnames(dataset)==project.col])
    }
  }else if(!is.na(project)){
    project.col<-"project"
    eventTable$project <- rep(project, nrow(dataset))
  }

  ## 1. get eventID
  eventTable$eventID <- as.character(dataset[,colnames(dataset)==eventID.col])
  # 1.1 QC the eventIDs: they have to be unique
  if(length(unique(eventTable$eventID)) != nrow(dataset)){
    eventTable$eventID <- make.names(eventTable$eventID, unique=TRUE)
  }

  # 2. get parentEventID
  if(!is.na(parentEventID.col)){
    if(!parentEventID.col %in% colnames(dataset)){
      stop("parentEventID.col not found")
    }else{
      eventTable$parentEventID <- as.character(dataset[,colnames(dataset)==parentEventID.col])
    }
  }

  # make all parentEvents that are NA ""
  if(any(is.na(eventTable$parentEventID))){
    eventTable$parentEventID[is.na(eventTable$parentEventID)]<-""
  }

  # complete the type column
  eventTable$type[!eventTable$eventID %in% eventTable$parentEventID]<-"event"
  eventTable$type[eventTable$eventID %in% eventTable$parentEventID]<-"parentEvent"

  # 3. complete hierarchy
  if(complete.hierarchy){
    eventTable$original <- TRUE
    # 3.1 make sure there is a parentEvent column to work with
    if(is.na(parentEventID.col)){
      # create a parentEventID column with the provided project
      eventTable$parentEventID <- eventTable$project
    }

    # 3.2 check if the parentEventIDs are in the eventID column
    # if not: complete
    parentEvents <- setdiff(unique(eventTable$parentEventID), "")
    if(!all(parentEvents %in% eventTable$eventID)){
      # not all parentEvents are in the event column
      for(pi in 1:length(parentEvents)){
        if(!parentEvents[pi] %in% eventTable$eventID & !is.na(parentEvents[pi])){
          if(!is.na(project.col) & !is.na(parentEventID.col)){
            # add the (true) parentEvents and their projects
            prj <- eventTable[eventTable$parentEventID==parentEvents[pi],]$project
            eventTable <- rbind(data.frame(original_name=NA, eventID=parentEvents[pi],
                                           parentEventID=prj, project=prj, type="parentEvent",
                                           original=FALSE, stringsAsFactors = FALSE),
                                eventTable)
          }else if(!is.na(project.col) & is.na(parentEventID.col)){
            # add the parentEvents, but the are the same as the projects
            eventTable <- rbind(data.frame(original_name=NA, eventID=parentEvents[pi],
                                           parentEventID=NA, project=parentEvents[pi], type="project",
                                           original=FALSE, stringsAsFactors = FALSE),
                                eventTable)
          }else{
            # add parentEvents, there are no projects
            eventTable <- rbind(data.frame(original_name="", eventID=parentEvents[pi],
                                           parentEventID="", type="parentEvent",
                                           original=FALSE, stringsAsFactors = FALSE),
                                eventTable)
          }
          
          rownames(eventTable)[1]<-parentEvents[pi]
        }
      }
    }

    # 3.3 if there is a project, put this as the higest event.
    # else=> don't change anything about the parentEvents anymore.
    if(!is.na(project.col)){
      projects<-as.character(unique(eventTable$project))
      if(!all(is.na(projects))){
        if(length(projects)==1){
          if(!projects %in% unique(eventTable$parentEventID)){
            eventTable[eventTable$parentEventID=="",]$parentEventID <- projects
          }else{
            for(pri in 1:nrow(eventTable)){
              if(eventTable[pri,]$parentEventID=="" & eventTable[pri,]$eventID!=projects){
                eventTable[pri,]$parentEventID <- projects
              }
            }
          }
        }
        if(!all(projects %in% eventTable$eventID)){
          # add the projects as events
          for(pi in 1:length(projects)){
            if(!projects[pi] %in% eventTable$eventID){
              eventTable <- rbind(data.frame(original_name=NA, eventID=projects[pi],
                                             parentEventID=NA, project="", type="project",
                                             original=FALSE, stringsAsFactors = FALSE),
                                  eventTable)
              rownames(eventTable)[1]<-projects[pi]
            }
          }
        }else{
          # make sure projects as events have no parentEvent
          eventTable[eventTable$eventID %in% projects,]$parentEventID <- NA
        }
      }
    }
  }
  # 4. finalize
  return(eventTable)
}

