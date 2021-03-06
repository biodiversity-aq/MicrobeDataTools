% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Format_sequenceData_ENA.R
\name{renameSequenceFiles}
\alias{renameSequenceFiles}
\title{renames sequence files in a directory}
\usage{
renameSequenceFiles(name_df, colNameOld="OldName", colNameNew="NewName", file.dir=NULL, paired=TRUE, seq.file.extension=".fastq.gz", ask.input=TRUE, pairedEnd.extension=c("_1", "_2"))
}
\arguments{
\item{name_df}{a data.frame. The dataframe that lists the original sample same (without file extension) alongside the new file names (also without their file extension)}

\item{colNameOld}{a character string. The name of the column with the original file names}

\item{colNameNew}{a character string. The name of the column with the new file names}

\item{file.dir}{a character string. The path to the directory where the sequence files are stored}

\item{paired}{boolean. wether or not the sequence files are paired-end (forward _1, reverse_2) or single-end}

\item{seq.file.extension}{a character string. The file-extension of the sequence files}

\item{ask.input}{boolean. Will give a warning message before continuing.}

\item{pairedEnd.extension}{a character vector of length 2. If the data is paired-end data, specify the forward (first element of te vector) and reverse (second) extension tags here. Default is c("_1", "_2")}
}
\value{
the number of files changed.
}
\description{
This function changes the filenames in the given directory according to the information provided in name_df
}
\details{
fastq files from sequencing facilities often come with long and complex file names that were automatically generated by the sequencer machine and no longer resemble the original name of the sample. This function is part of tools that help to get file names and easily convert them back into the original file names. It makes use of a table where each sequence file name is linked to a new name desired by the user. This can be generated by the FileNames.to.Table function.
}
\seealso{
Other data archiving functions: 
\code{\link{FileNames.to.Table}()},
\code{\link{commonTax.to.NCBI.TaxID}()},
\code{\link{get.ENAName}()},
\code{\link{prep.metadata.ENA}()},
\code{\link{sync.metadata.sequenceFiles}()}
}
\author{
Maxime Sweetlove ccBY 4.0 2020
}
\concept{data archiving functions}
