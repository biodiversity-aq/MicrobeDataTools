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
# Tools for downloading data
#==============================================================

#' Fetch BioProject metadata
#' @author Maxime Sweetlove ccBY 4.0 2019
#' @description downloads a minimal set of sample metadata of all samples ("Runs") within a BioProject (argument BioPrjct) from the International Nucleotide Sequence Database Consortium (INSDC) databases.
#' @usage get.BioProject.metadata.INSDC(BioPrjct, just.names=FALSE)
#' @param BioPrjct  A chracter string. A single Bioproject ID.
#' @param just.names  Boolean. If TRUE, only the INSDC sample names are returned, else all data are returned. default FALSE
#' @details BioProjects combine all biological nucleotide sequence data related to a single initiative, originating from a single organization. With each sample ("Run") within a BioProject, there is additional data associated that is crucial to the correct interpretation of the nucleotide sequence data, but is not automatically downloaded along with it. The get.BioProject.metadata.INSDC function will fetch the most basic metadata of a BioProject from the INSDC repositories to complete the nucleotide sequence dataset, using E-utils API function of NCBI. These basic metadata typically include Run number, relsease date, load date, spots, bases, av_MB and download path. Note that the data returned by get.BioProject.metadata.INSDC does not include all the metadata associated with a BioProject. Other information, like coordinates, sampling dates or environmental measurements may also be available, but require the user to register at NCBI and request a personal API-key (this is required by NCBI to acces their data since 2017). The complete set of additional data can be downloaded using the get.sample.attributes.INSDC function, given a user-specified API-key. see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
#' @seealso get.sample.attributes.INSDC
#' @return if get.BioProject.metadata.INSDC(just.names=FALSE) (default) a data.frame with n rows an m columns is returned, n being the number of samples in the BioProject, and m being the number of variables found. If get.BioProject.metadata.INSDC(just.names=TRUE), a character vector of length n is returned with the sample numbers ("Run numbers", SRR numbers)
#' @references Sayers, E. (2009) The E-utilities In-Depth: Parameters, Syntax and More, https://www.ncbi.nlm.nih.gov/books/NBK25499/
#' @example metadata_PRJNA369175 <- get.BioProject.metadata.INSDC(BioPrjct="PRJNA369175", just.names=FALSE)
#' @export
get.BioProject.metadata.INSDC <- function(BioPrjct=NA, just.names=FALSE){
  if(! is.character(BioPrjct) | c(NULL,NA) %in% BioPrjct | length(BioPrjct)>1){
    stop("incorrect BioProject provided.
         Expected input is a single Bioproject ID as a character string.")
  }

  sra_url <- paste("http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=", BioPrjct, sep="")
  tmpFile <- tempfile()
  download.file(sra_url, destfile = tmpFile, "wget", quiet = TRUE)
  if(readLines(tmpFile, n=1)==""){
    stop(paste("No metadata found for ", BioPrjct, ". Could not download the project metadata.", sep=""))
  } else{
    RawMetadata <- read.csv(tmpFile, header=TRUE)
  }
  file.remove(tmpFile)
  if(just.names){
    return(c(as.character(RawMetadata$Run)))
  } else{
    return(RawMetadata)
  }
}

#' Downloads all sequence sample attributes
#' @author Maxime Sweetlove ccBY 4.0 2019
#' @description Downloads all sample attributes (that is, additional environmental or other associated data) from INSDC. Note: requires a user-specified API-key to acces the INSDC databases. see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ to generate an API-key
#' @param sampleID a list. a list of one ore more SRA sample IDs (that is "Run" numbers). This argument can be left blank if input is provided for the BioPrjct argument (see further)
#' @param apiKey a character string. A personal API-key to the acces the NCBI databases, and required to use the Entrez Programming Utilities (E-utilities). An API-key (API stands for application programming interface) is a unique identifier used to authenticate a user. A personal API-key can easily be generated at requested at https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
#' @param BioPrjct a character string or a vector with character strings. Providing the associated BioProject numbers helps to give more understandable error messages. Alternativey, if the sampleID argument is empty, all sample IDs from the given BioProject numbers will be given to the sampleID argument.
#' @details Each sequence data sample ("Run") typically has additional measurements or metadata associated with it. However, these are difficult to find, and cannot be downloaded together with the nucleotide sequence data. The get.sample.attributes.INSDC function will fetch the all the metadata of the samples given to the sampleID (or BioPrjct) argument from the INSDC. To do this, get.sample.attributes.INSDC will use the Entrez Programming Utilities (E-utilities) from NCBI, and will access the databases with the API-key provided by the user. see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/ for more info on getting an API-key.
#' @return a dataframe with the data found
#' @export
get.sample.attributes.INSDC <- function(sampleID=NA, apiKey=NA, BioPrjct=NA){
  #requires xml2 to parse the xml file, reshape2 long to wide format

  if(is.na(apiKey)){stop("No apiKey provided.
                         A personal API key for the Entrez Programming Utilities (E-utilities) should be requested at the NCBI website.")
  }

  # check if param sampleID is present. if not, check for BioPrjct
  if(class(sampleID)!="character" | length(c(sampleID))==0 |
     c(NA, NULL) %in% sampleID | c("") %in% sampleID){
    if(class(BioPrjct)!="character" | length(c(BioPrjct))==0 |
       c(NA, NULL) %in% BioPrjct | c("") %in% BioPrjct){
      stop("No sampleID or BioPrjct provided,
           or one of the sampleIDs/BioPrjct was NA, NULL or \"\"")
    }else{
      sampleID<-c()
      for(BP in BioPrjct){
        sampleID <- c(sampleID, get.BioProject.metadata.INSDC(BP, just.names=TRUE))
      }
    }
  }

  # see how many requests (sample ID's) will need to be handled
  # SRA can't handle large requests (lets say >50 at a time), so the IDs will be split into chunks of no mare than 50
  sampleID_setlist <- split(sampleID, ceiling(seq_along(sampleID)/50))
  env_metadata_full <- data.frame()
  for(sampleID_set in sampleID_setlist){
    sampleID_setSring <- paste(c(unlist(sampleID_set)), collapse="+")
    request_url <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=",
                         sampleID_setSring, "&retmax=99999&api_key=", apiKey, sep="")
    metaXML <- tryCatch(read_xml(request_url[1]),
                        error = function(x){return(NULL)}
    )
    if(is.null(metaXML)){
      warning("No additional data found for one or more samples: better check this out.")
      env_metadata <- NULL
      next
    } else{
      # parse the xml file
      recs <- xml_find_all(metaXML, "//EXPERIMENT_PACKAGE/SAMPLE")
      attr_sample <- trimws(xml_attr(recs, "accession"))
      attr_name <- xml_find_all(metaXML, "//EXPERIMENT_PACKAGE/RUN_SET/RUN")
      attr_name <- trimws(xml_attr(attr_name, "accession"))
      attr_vals <- xml_find_all(metaXML, "//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/*[self::VALUE]")
      attr_vals <- trimws(xml_text(attr_vals))
      attr_cols <- xml_find_all(metaXML, "//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/*[self::TAG]")
      attr_cols <- xml_text(attr_cols)
      attr_BioPrj <- xml_find_all(metaXML, "//STUDY/IDENTIFIERS/EXTERNAL_ID")
      attr_BioPrj <- xml_text(attr_BioPrj)

      attr_LibNname <- xml_find_all(metaXML, "//EXPERIMENT_PACKAGE/EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/*[self::LIBRARY_NAME]")
      attr_LibNname <- trimws(xml_text(attr_LibNname))

      n_samples <- as_list(xml_find_all(metaXML, "//SAMPLE_ATTRIBUTES"))
      n_samples <- lapply(n_samples, function(x){length(x)})
      n_samples <- unlist(n_samples)
      namesVec <- c()
      for(ni in 1:length(attr_name)){
        namesVec <- c(namesVec, rep(attr_name[ni], n_samples[ni]))
      }
      env_metadata <- data.frame(attr_name=namesVec,
                                 attr_col=attr_cols, value=attr_vals)

      env_metadata <- reshape2::dcast(env_metadata, attr_name ~ attr_cols, value.var="value")
      env_metadata$BioProject <- attr_BioPrj
      env_metadata$SRA_sample <- attr_sample

      env_metadata$LibraryName <- attr_LibNname

    }
    if(nrow(env_metadata_full)==0){
      env_metadata_full <- env_metadata
    }else{
      env_metadata_full <- combine.data(env_metadata_full, env_metadata,
                                              fill=NA, variables.as.cols=TRUE)
    }
  }
  row.names(env_metadata_full)<-env_metadata_full$attr_name
  colnames(env_metadata_full)<-gsub(" ", "_", colnames(env_metadata_full))
  colnames(env_metadata_full)<-gsub("-", "_", colnames(env_metadata_full))
  colnames(env_metadata_full)<-gsub("/", "_", colnames(env_metadata_full))
  colnames(env_metadata_full)<-gsub("\\(", "", colnames(env_metadata_full))
  colnames(env_metadata_full)<-gsub("\\)", "", colnames(env_metadata_full))
  return(env_metadata_full)
  }



#' Download open sequence data to your computer.
#' @author Maxime Sweetlove ccBY 4.0 2019
#' @description Downloads (high throughput) nucleotide sequence datasets that are deposited on the International Nucleotide Sequence Database Consortium (INSDC, e.g. SRA, ENA, GenBank,...). Any possible metadata and environmental data is also downloaded.
#' @usage download.sequences.INSDC(BioPrj = c(), destination.path = NA, apiKey=NA, unzip = FALSE, keep.metadata = TRUE, download.sequences = TRUE)
#' @param BioPrj a list with character strings. A list of one or more BioProject numbers to be downloaded. Required argument.
#' @param destination.path a character string. The path to the directory where all the downloaded sequence data needs to go
#' @param apiKey a character string. Only required if download.sequences.INSDC(keep.metadata=TRUE). A personal API-key to the acces the NCBI databases, and required to use the Entrez Programming Utilities (E-utilities). An API-key (API stands for application programming interface) is a unique identifier used to authenticate a user. You can easily generate an API-key: see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
#' @param unzip boolean. If TRUE, the all *.fastq.gz files in the destination.path will unzipped. Default FALSE
#' @param keep.metadata boolean. If TRUE, the downloaded metadata can be saved to a file (Console), if FALSE it is discarded. Default TRUE
#' @param download.sequences boolean. If TRUE, the sequences will be downloaded to the destination.path. If FALSE, no sequences are downloaded. Default TRUE
#' @details download.sequences.INSDC will write the sequence data (*.fastq.gz files) to a destiation path (e.g. a designated file), the metadata (if it should be kept) is written to the Console and should be caught in an R-varaiable that can be later written to a csv file by the user. Point of entry to INSDC is the SRA database from NCBI.
#' @return the sequence data are written to the destination.path, the metadata is returned as a data.frame to the console.
#' @export
download.sequences.INSDC <- function(BioPrj = c(), destination.path = NA, apiKey=NA,
                                     unzip = FALSE, keep.metadata = TRUE, download.sequences = TRUE){

  ### 1. some initial checks before continuing
  # 1.1. check destination.path
  if(is.na(destination.path)){destination.path <- getwd()}
  # 1.2. check BioPrj
  if(length(BioPrj) == 0){
    stop("No BioProject number was given (required).
         Valid input: a list of one or more BioProject numbers")
  }
  if(keep.metadata){
    metadata_all<-data.frame()
    # 1.3. check apiKey (only needed if keep.metadata is TRUE)
    if(is.na(apiKey)){stop("No apiKey provided.
                                         A personal API key for the Entrez Programming Utilities (E-utilities) should be requested at the NCBI website.\n")
    }
    cat("Notice!\nthe metadata will be retruned to the Console\nIf you did assign the output of this function to an R-object (using \"<-\"): better abort and restart now\n\n")
  }

  # 2. start going though the BioProjects
  cat("Getting the metadata ...\n")
  downloads_failed<-0
  for(BP in BioPrj){
    cat(paste("Processing BioProject ",BP,"...\n", sep=""))
    faultyBioPrj<-FALSE
    ### 2.1. Downloading the metadata from SRA to get the run numbers (=samples) of the BioProject
    ###    These run numbers are the required input to download the sequence data, as sequence data cannot be redectly downloaded via the BioProject number
    faultyBioPrj <- tryCatch({
      RawMetadata <- get.BioProject.metadata.INSDC(BP)
      cat(paste("\t",as.character(nrow(RawMetadata))," samples (Runs)...\n", sep=""))
    }, error = function(e) {
      return(TRUE)
    })

    if(! is.null(faultyBioPrj) && faultyBioPrj){
      downloads_failed <- downloads_failed+1
      cat(paste("\t No data found for ", BP,". Could not download the data.\n", sep=""))
    } else{
      faultyBioPrj<-FALSE
    }




    ### 2.2 add the FTP path to download the sequences
    ftp_url <- paste("https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=", BP,
                     "&result=read_run&fields=run_accession,fastq_ftp&download=txt", sep="")
    tmpFile <- tempfile()
    download.file(ftp_url, destfile = tmpFile, method="auto", quiet = TRUE)
    ftps <- read.table(tmpFile, header=TRUE)
    ftps <- ftps[ftps$run_accession[order(RawMetadata$Run)],]
    RawMetadata$ftp <- ftps$fastq_ftp
    file.remove(tmpFile)

    ### 2.2. if keep.metadata==TRUE, the raw metadata must be cleaned up, and if applicable added to the other metadata
    if(keep.metadata & !faultyBioPrj){
      cat("\tprocessing the metadata ...\n")
      sample_numbers<-as.character(RawMetadata$Run)
      env_data <- get.sample.attributes.INSDC(sampleID = sample_numbers,
                                              apiKey = apiKey, BioPrjct = BP)
      rownames(RawMetadata) <- RawMetadata$Run
      RawMetadata <- RawMetadata[order(row.names(RawMetadata)),]
      rownames(env_data) <- env_data$attr_name
      env_data <- env_data[order(row.names(env_data)),]
      MetaData <- cbind(RawMetadata, env_data)

      if(nrow(metadata_all)>0){
        metadata_all <- process.metadata(metadata=MetaData, add_to = metadata_all, strict.MIxS = FALSE, ask.input = FALSE)
      } else{
        metadata_all <- process.metadata(metadata=MetaData, strict.MIxS = FALSE, ask.input = FALSE)
      }
    }

    ### 2.3. Downloading the sequence data from SRA
    #      This requires a sqlite database of what is on SRA. This database just needs to be downloaded once, but can be very big.
    #      No way of getting te data without the SRA sqlite database...
    #      Also using the sample names gathered in (1) instead of the BioProject number
    if(download.sequences & !faultyBioPrj){
      cat("\tDownloading the sequence data ...\n")
      ftps <- as.character(RawMetadata$ftp)
      ftps <- sapply(ftps, function(x){strsplit(x, ";")})
      ftps <- unname(unlist(ftps))
      for(ftp_url in ftps){
        destfile <- strsplit(ftp_url, "/")[[1]]
        destfile <- destfile[length(destfile)]
        download.file(ftp_url, destfile = file.path(destination.path, destfile),method="auto")
      }
    }else{
      downloads_failed <- downloads_failed+1
    }
  }

  ### 2.4. unzipping the data, if required by the user
  if(unzip == TRUE & !faultyBioPrj){
    cat("Unzipping files ...\n")
    fileList<- list.files(path=getwd(), pattern='*.fastq.gz')
    for(file in fileList){
      fileName= strsplit(file, ".gz")[[1]]
      system2(command="gunzip", args=paste(" -c ", file," > ", fileName, sep=''), stdout=TRUE)
      system2(command="rm", args=file)
    }
  }

  ### 3. finalizing
  if(downloads_failed < length(BioPrj)){
    outputmessage <- paste("Finished processing\nThe files are in ",destination.path,"\n", sep="")
    if(keep.metadata){
      outputmessage <- paste(outputmessage, "Please consider the warning messages for a manual quality controll of the metadata\nand environmental data\n")
    }
  } else{
    outputmessage <- paste("Finished processing\nNo files were downloaded\n", sep="")
  }
  cat(outputmessage)
  if(keep.metadata){return(metadata_all)}
}





