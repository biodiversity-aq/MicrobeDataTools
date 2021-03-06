% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataQC_Utils.R
\name{dataQC.dateCheck}
\alias{dataQC.dateCheck}
\title{find and standardize dates in a dataframe}
\usage{
dataQC.dateCheck(dataset, date.colnames)
}
\arguments{
\item{dataset}{dataframe. The dataset where the date column should be found}

\item{date.colnames}{character vector. a list of potential names for the column with the date. e.g. c("date", "Date", "collection date")}
}
\value{
a list of length 2, with "$values" a vactor with the same number of rows as the dataset argument containing the corrected date values, and "$warningmessages" a vector with potential warning messages as character strings.
}
\description{
looks in the columns of a dataset for a column with dates and transforms them to the YYYY-MM-DD format.
}
\details{
The date column is found based on a user-provided list of possible names to look for (data.colnames argument). If a columnname is found that corresponds to a term in the list, the dates will be convered to the YYYY-MM-DD format, if the original format can be recognized.
}
\seealso{
Other quality control functions: 
\code{\link{dataQC.LatitudeLongitudeCheck}()},
\code{\link{dataQC.TaxonListFromData}()},
\code{\link{dataQC.TermsCheck}()},
\code{\link{dataQC.completeTaxaNamesFromRegistery}()},
\code{\link{dataQC.eventStructure}()},
\code{\link{dataQC.findNames}()},
\code{\link{dataQC.generate.footprintWKT}()},
\code{\link{dataQC.guess.env_package.from.data}()},
\code{\link{dataQC.taxaNames}()}
}
\author{
Maxime Sweetlove CC-BY 4.0 2019
}
\concept{quality control functions}
