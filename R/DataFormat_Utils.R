#==============================================================
# Author Maxime Sweetlove
# lisence CC 4.0
# Part of the POLA3R website (successor or mARS.biodiversity.aq)
# version 1.0 (2019-09-20)
# file encdong UTF-8
#
#==============================================================

#' convert a DarwinCore extended Measurement or Fact (eMoF) file to a wide tabular format
#' @family formating functions
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @description converts a DarwinCore extended Measurement Or Fact (eMOF) file, which is has a "long" file format into a wide tabular format
#' @usage eMoF.to.wideTable(dataset)
#' @param dataset a dataframe of the eMoF file
#' @details Long formated data if great for data arciving, but is difficult to use in day-to-day statistical analyses. This function extracts the data from an eMoF file and puts it in a wide sample x variable table
#' @return a list of length 3: "$data" the data in a wide formt, "$units" the units, and "$method" the methods
#' @export
eMoF.to.wideTable <- function(dataset){
  ## converts an extended Measurement of Fact table to a regular wide table
  #requires tidyr
  # check input
  if(!all(c("measurementType", "measurementValue", "measurementUnit") %in% colnames(dataset))){
    stop("Invalid DwC eMoF file.")
  }

  if("eventID" %in% colnames(dataset)){
    dataLong<-dataset[,colnames(dataset) %in% c("eventID", "measurementType", "measurementValue")]
    dataWide <- tidyr::spread(dataLong, eventID, measurementValue)
  }else if("occurrenceID" %in% colnames(dataLong)){
    dataLong<-dataset[,colnames(dataset) %in% c("occurrenceID", "measurementType", "measurementValue")]
    dataWide <- tidyr::spread(dataLong, occurrenceID, measurementValue)
  }else{
    stop("Invalid DwC eMoF file.")
  }

  dataWide_rows <- colnames(dataWide)[-1]
  dataWide<-data.frame(dataWide)
  rownames(dataWide) <- dataWide$measurementType
  dataWide <- data.frame(t(dataWide[,!colnames(dataWide) %in% "measurementType"]))
  rownames(dataWide) <- dataWide_rows

  dataUnits<-dataset[,colnames(dataset) %in% c("measurementType", "measurementUnit")]
  dataUnits <- unique(dataUnits)
  Units_list<-c()
  Units_list[as.character(dataUnits$measurementType)]<-as.character(dataUnits$measurementUnit)

  Method_list<-NA
  if("measurementMethod" %in% colnames(dataset)){
    Method_list<-c()
    dataMethods<-dataset[,colnames(dataset) %in% c("measurementType", "measurementMethod")]
    dataMethods <- unique(dataMethods)
    Method_list[as.character(dataMethods$measurementType)]<-as.character(dataMethods$measurementMethod)
  }

  # not yet included:
  #measurementAccuracy
  #measurementDeterminedDate
  #measurementDeterminedBy
  #measurementRemarks
  #measurementID
  #measurementTypeID
  #measurementValueID
  #measurementUnitID

  return(list(data=dataWide, units=Units_list, method=Method_list))

}

#' convert dataframe to a DarwinCore extended Measurement or Fact (eMoF) file
#' @family formating functions
#' @author Maxime Sweetlove ccBY 4.0 2020
#' @description converts a dataframe to a DarwinCore extended Measurement Or Fact (eMOF) file
#' @usage wideTable.to.eMoF(dataset)
#' @param metadata.object a MIxS.metadata object
#' @param variables a character vector. a list of the variables that need to be included in the eMoF
#' @details extended Measurement or Fact (eMoF) as a DarwinCore extension to standardize environmental or other additional data in a computer readable fashon. This standard structures data into a long format (a column with sample name, variable name and value). This function converts more commonly used wide format tables (that is, structured like a matrix, e.g. samples as rows and variables as columns) into the correct long format that complies to eMoF
#' @return a data.frame formatted as an extended Measurement or Fact table
#' @export
wideTable.to.eMoF <- function(metadata.object, variables=NA){
  ## converts a regular wide table to an extended Measurement of Fact table 
  #requires tidyr
  
  # check input
  if(!check.valid.MIxS.metadata(metadata.object)){
    stop("invalid input")
  }
  bad.vars<-setdiff(variables, colnames(metadata.object@data))
  if(length(bad.vars)>0){
    stop(paste("The following variables were not present in the dataset: ", paste(bad.vars, collapse=", "), sep=""))
  }
  
  #format the data
  if("eventID" %in% colnames(metadata.object@data)){
    dataset <- metadata.object@data[,c("eventID", variables)]
  }else if("occurrenceID" %in% colnames(metadata.object@data)){
    dataset <- metadata.object@data[,c("occurrenceID", variables)]
  }else if("original_name" %in% colnames(metadata.object@data)){
    dataset <- metadata.object@data[,c("original_name", variables)]
    colnames(dataset)[1]<-"eventID"
  }else{
    stop("no valid samplenames found (eventID or occurrenceID)")
  }
  
  data_long <- suppressWarnings(gather(dataset, measurementType, measurementValue, variables, factor_key=TRUE))
  
  # measurementUnit
  if(length(metadata.object@units)>0){
    data_long$measurementUnit <- data_long$measurementType
    data_long$measurementUnit <- unname(sapply(data_long$measurementUnit, function(x){
      gsub(x,metadata.object@units[names(metadata.object@units)==x],x)
    }))
  }
  
  # not yet included:
  # measurementMethod
  #measurementAccuracy
  #measurementDeterminedDate
  #measurementDeterminedBy
  #measurementRemarks
  #measurementID
  #measurementTypeID
  #measurementValueID
  #measurementUnitID
  
  return(data_long)
  
}

#' merge dataframes
#' @author Maxime Sweetlove ccBY 4.0 2019
#' @family formating functions
#' @description combine.data.frame merges two dataframes, completing the rows and columns that are not shared by the dataframes.
#' @usage combine.data.frame(df1, df2, fill=NA, merge.cols=TRUE, original_rowName.col=TRUE, merge.rows="df1")
#' @param df1 a dataframe
#' @param df2 a dataframe
#' @param fill character or NA. A value to put into the cells that have no data. Default NA
#' @param merge.cols logical. Merge colomns with common name. Default TRUE
#' @param original_rowName.col logical. if TRUE, add a column with the original rownames. Default TRUE
#' @param merge.rows character. If the data.frames df1 and df2 have common rows, the values for these rows of this dataset must be kept in favor of the other. If NA, both rows will be kept. Either "df1", "df2" or NA. Default "df1"
#' @details Columns with matching names can or cannot be merged, rows are automatically bound (a wrapper of rbind), not merged. Any missing data as a result of the non-matchng columns will be filled by the fill argument.
#' @return a data.frame
#' @export
combine.data.frame <- function(df1, df2, fill=NA, 
                               merge.cols=TRUE, original_rowName.col=TRUE, merge.rows="df1"){
  if(!class(df1) %in% c("data.frame", "matrix") |
     !class(df2) %in% c("data.frame", "matrix") |
     ncol(df1) == 0 | nrow(df1) == 0 |
     ncol(df2) == 0 | nrow(df2) == 0){
    stop("invalid input, must be a dataframe of matrix with at leaste 1 column and row.")
  }
  if(!length(merge.rows)==1 && !merge.rows %in% c("df1", "df2", NA)){
    stop("incorrect value for parameter \"merge.rows\". Must be either \"df1\", \"df2\" or NA")
  }
  for(nc in 1:ncol(df1)){ #merging columns with factors is always trouble, so convert to character
    if(class(df1[,nc])=="factor"){
      df1[,nc] <- as.character(df1[,nc])
    }
  }
  for(nc in 1:ncol(df2)){ #idem
    if(class(df2[,nc])=="factor"){
      df2[,nc] <- as.character(df2[,nc])
    }
  }

  shared_cols <- intersect(colnames(df1), colnames(df2))
  if(merge.cols){
    df11 <- df1[,!colnames(df1) %in% shared_cols, drop=FALSE]
    df22 <- df2[,!colnames(df2) %in% shared_cols, drop=FALSE]
  }else{
    df22 <- df2
    df11 <- df1
    if(length(shared_cols)>0){
      colnames(df22)[colnames(df22) %in% shared_cols] <- paste(colnames(df22)[colnames(df22) %in% shared_cols], "_2", sep="")
      colnames(df2) <- colnames(df22)
      shared_cols <- c()
    }
  }
  
  shared_rows <- intersect(rownames(df1), rownames(df2))
  if(length(shared_rows)>0){
    if(merge.rows=="df1"){
      df2<-df2[!rownames(df2) %in% shared_rows,]
      df22<-df22[!rownames(df22) %in% shared_rows,]
    }else if(merge.rows=="df2"){
      df1<-df1[!rownames(df1) %in% shared_rows,]
      df11<-df11[!rownames(df11) %in% shared_rows,]
    }
  }

  if(merge.cols){
    if(length(shared_rows)>0 & original_rowName.col){
      df1$original_rowName <- row.names(df1)
      df2$original_rowName <- row.names(df2)
      shared_cols<-c(shared_cols, "original_rowName")
    }
    if(!ncol(df22)==0 && !nrow(df22)==0){ #if ncol(df22)==0, then df2 has no other colnames than df1
      df_UpRight <- data.frame(matrix(nrow = nrow(df1), ncol = ncol(df22), data=fill))
      colnames(df_UpRight) <- colnames(df22)
      rownames(df_UpRight) <- rownames(df1)
      df1_b <- cbind(df1, df_UpRight)
    }else{
      df1_b <- df1
    }
    if(!ncol(df11)==0 & !nrow(df11)==0){
      df11 <- df1[,!colnames(df1) %in% shared_cols, drop=FALSE]
      df_DownLeft <- data.frame(matrix(nrow = nrow(df2), ncol = ncol(df11), data=fill))
      rownames(df_DownLeft) <- rownames(df2)
      colnames(df_DownLeft) <- colnames(df11)
      df2_b <- cbind(df_DownLeft, df2)
    }else{
      df2_b <- df2
    }

    df_out <- data.frame(rbind(df1_b, df2_b))
    colnames(df_out)<-colnames(df1_b)

    df_out[c((nrow(df1_b)+1):(nrow(df_out))),shared_cols]<-df2[,shared_cols, drop=FALSE]
  }else{ #not merging columns
    if(length(setdiff(rownames(df1), shared_rows))>0){
      df_x <- data.frame(matrix(nrow = nrow(df1)-length(shared_rows), ncol = ncol(df22), data=fill))
      colnames(df_x) <- colnames(df22)
      rownames(df_x) <- setdiff(rownames(df1), shared_rows)
      df_notshared_1 <- cbind(df1[setdiff(rownames(df1), shared_rows),], df_x)
    }else{
      df_notshared_1 <- data.frame()
    }
    if(length(setdiff(rownames(df2), shared_rows))>0){
      df_x <- data.frame(matrix(nrow = nrow(df2)-length(shared_rows), ncol = ncol(df11), data=fill))
      colnames(df_x) <- colnames(df11)
      rownames(df_x) <- setdiff(rownames(df2), shared_rows)
      df_notshared_1 <- cbind(df_x, df2[setdiff(rownames(df2), shared_rows),])
    }else{
      df_notshared_2 <- data.frame()
    }

    df_notshared <- rbind(df_notshared_1, df_notshared_2)
    df_shared <- cbind(df1[shared_rows,], df2[shared_rows,])
    df_out <- data.frame(rbind(df_shared, df_notshared))
  }

  return(df_out)
}

#' merge MIxS.metadata objects
#' @author Maxime Sweetlove ccBY 4.0 2019
#' @family formating functions
#' @description combine.data combines two MIxS.metadata objects into one, merging all common variables (columns as default).
#' @usage combine.data(d1, d2, fill=NA, variables.as.cols=TRUE)
#' @param d1 a MIxS.metadata object
#' @param d2 a MIxS.metadata object
#' @param fill character or NA. A value to put into the cells that have no data.
#' @param variables.as.cols boolean. If TRUE, the input data is assumed to have rows as samples and variables/parameters as columns. If FALSE the data is formatted the other was around. default TRUE
#' @details Variables with matching names are merged, if variables.as.cols=TRUE this means columns are merged, if FALSE rows are merged. Any missing data as a result of the non-matchng variable-name will be filled by the fill argument.
#' @return a MIxS.metadata object
#' @export
combine.data <- function(d1, d2, fill=NA, variables.as.cols=TRUE){
  if(class(d1) %in% c("data.frame", "matrix") &
     class(d2) %in% c("data.frame", "matrix")){
    if(ncol(d1) == 0 | nrow(d1) == 0 |
       ncol(d2) == 0 | nrow(d2) == 0){
      stop("invalid input, each dataset must at leaste have 1 column and row.")
    }
    class_out <- "data.frame"
    if(variables.as.cols){
      df1<-d1
      df2<-d2
    } else{
      df1<-data.frame(t(d1))
      df2<-data.frame(t(d2))
    }
  } else if(check.valid.MIxS.metadata(d1) & check.valid.MIxS.metadata(d2)){
    class_out <- "MIxS.metadata"
    variables.as.cols<-TRUE
    df1<-d1@data
    df2<-d2@data
  } else{
    stop("invalid input, both datasets must be a dataframe, matrix or MIxS.metadata class object.
         possible causes of this error include:
         - a dataframe/matrix or MIxS.metadata@data object had 0 columns or rows.
         - the two input dataset differed in their object class.
         - the MIxS.metadata object was invalid.")
  }

  for(nc in 1:ncol(df1)){ #merging columns with factors is always trouble, so convert to character
    if(class(df1[,nc])=="factor"){
      df1[,nc] <- as.character(df1[,nc])
    }
  }
  for(nc in 1:ncol(df2)){ #idem
    if(class(df2[,nc])=="factor"){
      df2[,nc] <- as.character(df2[,nc])
    }
  }
  shared_cols <- intersect(colnames(df1), colnames(df2))
  df11 <- df1[,!colnames(df1) %in% shared_cols, drop=FALSE]
  df22 <- df2[,!colnames(df2) %in% shared_cols, drop=FALSE]

  shared_rows <- intersect(rownames(df1), rownames(df2))

  if(length(shared_rows)>0 & original_rowName.col){ #new column with original rownames
    df1$original_rowName <- row.names(df1)
    df2$original_rowName <- row.names(df2)
    shared_cols<-c(shared_cols, "original_rowName")
  }

  if(!ncol(df22)==0){ #if ncol(df22)==0, then df2 has no other colnames than df1
    df_UpRight <- data.frame(matrix(nrow = nrow(df1), ncol = ncol(df22), data=fill))
    colnames(df_UpRight) <- colnames(df22)
    rownames(df_UpRight) <- rownames(df1)
    df1_b <- cbind(df1, df_UpRight)
  }else{
    df1_b <- df1
  }

  if(!ncol(df11)==0){
    df11 <- df1[,!colnames(df1) %in% shared_cols, drop=FALSE]
    df_DownLeft <- data.frame(matrix(nrow = nrow(df2), ncol = ncol(df11), data=fill))
    rownames(df_DownLeft) <- rownames(df2)
    colnames(df_DownLeft) <- colnames(df11)
    df2_b <- cbind(df_DownLeft, df2)
  }else{
    df2_b <- df2
  }

  df_out <- data.frame(rbind(df1_b, df2_b))
  colnames(df_out)<-colnames(df1_b)
  df_out[c((nrow(df1_b)+1):(nrow(df_out))),shared_cols]<-df2[,shared_cols, drop=FALSE]

  if(!variables.as.cols){
    df_out <- data.frame(t(df_out))
  }

  if(class_out=="MIxS.metadata"){
    df_units <- c()
    df_section <- c()
    for(cname in colnames(df_out)){
      df_u1 <- d1@units[cname]
      df_u2 <- d2@units[cname]
      if(!is.na(df_u1) & !is.na(df_u2)){
        if(df_u1 != df_u2){
          stop("some variable units do not match. Correct and try again")
        }
      }
      if(is.na(df_u1)){
        df_u1 <- df_u2
      }
      df_units[cname] <- as.character(df_u1)

      df_s1 <- d1@section[cname]
      if(is.na(df_s1)){
        df_s1 <- d2@section[cname]
      }
      df_section[cname] <- as.character(df_s1)
    }

    df_t1<-d1@type
    df_t2<-d2@type
    if(df_t1=="flexible.MIxS" | df_t2=="flexible.MIxS"){
      df_t<-"flexible.MIxS"
    }else{
      df_t<-"strict.MIxS"
    }

    df_env1<-d1@env_package
    df_env2<-d2@env_package
    if(df_env1==df_env2){
      df_env<-df_env1
    }else{
      df_env<-"multiple_packages"
    }

    df_out <- new("MIxS.metadata",
                  data = df_out,
                  units = df_units,
                  section = df_section,
                  type = df_t,
                  env_package = df_env
    )

  }

  return(df_out)
  }

#' convert a wide table to hierarchical table
#' @author Maxime Sweetlove CC-BY 4.0 2019
#' @family formating functions
#' @description turns a regular wide table into a hierarchical recursive table.
#' @param dataset a data.frame. The wide table to be transformed, with columns as hierarchical variables. T
#' @param col_hierarchy a data.frame with two columns named "child" and "parent". The hierarchical relations between the columns formatted in a child-parent table listing the column names. Use "root" of NA for columns with nu perent.
#' @usage wideTab.to.hierarchicalTab(dataset, col_hierarchy)
#' @details Columns in the input data with different levels of a hierarchy are poolled into one column with values and a second column indicating the parent of the value.
#' @return a data.frame with three columns: one for the value (cell content in the origical data.frame), it's rank (the original column name), and parent (cell content of the column that is one up in the hierarchy)
#' @export
wideTab.to.hierarchicalTab <- function(dataset, col_hierarchy){
  if(colnames(col_hierarchy) != c("child", "parent")){
    stop("invalid input of the col_hierarchy argument.\nMust be a data.frame with two columns, names \"child\" and \"parent\".")
  }
  if(nrow(col_hierarchy) != ncol(dataset)){
    stop("invalid input of the col_hierarchy argument.\n Every column in the dataTab input must be represented as a row in the col_hierarchy data.frame.")
  }
  recs_out <- data.frame(matrix(nrow=0, ncol=3))
  colnames(recs_out) <- c("value", "rank", "parent")

  dataset <- dataset[!duplicated(dataset),]

  for(i in 1:nrow(dataset)){
    for(j in 1:ncol(dataset)){
      val <- as.character(dataset[i,j])
      rank <- colnames(dataset[,j, drop=FALSE])
      rankParent <- as.character(col_hierarchy[col_hierarchy$child == rank,]$parent)
      if(rankParent %in% colnames(dataset)){
        valParent <- dataset[i, rankParent]
      } else{
        valParent <- NA
      }
      recs_out <- rbind(recs_out,data.frame(value=val, rank=rank, parent=valParent))
    }
  }
  return(recs_out)
}



