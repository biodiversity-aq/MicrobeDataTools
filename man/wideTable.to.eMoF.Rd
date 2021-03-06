% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataFormat_Utils.R
\name{wideTable.to.eMoF}
\alias{wideTable.to.eMoF}
\title{convert dataframe to a DarwinCore extended Measurement or Fact (eMoF) file}
\usage{
wideTable.to.eMoF(dataset)
}
\arguments{
\item{metadata.object}{a MIxS.metadata object}

\item{variables}{a character vector. a list of the variables that need to be included in the eMoF}
}
\value{
a data.frame formatted as an extended Measurement or Fact table
}
\description{
converts a dataframe to a DarwinCore extended Measurement Or Fact (eMOF) file
}
\details{
extended Measurement or Fact (eMoF) as a DarwinCore extension to standardize environmental or other additional data in a computer readable fashon. This standard structures data into a long format (a column with sample name, variable name and value). This function converts more commonly used wide format tables (that is, structured like a matrix, e.g. samples as rows and variables as columns) into the correct long format that complies to eMoF
}
\seealso{
Other formating functions: 
\code{\link{combine.data.frame}()},
\code{\link{combine.data}()},
\code{\link{eMoF.to.wideTable}()},
\code{\link{wideTab.to.hierarchicalTab}()}
}
\author{
Maxime Sweetlove ccBY 4.0 2020
}
\concept{formating functions}
