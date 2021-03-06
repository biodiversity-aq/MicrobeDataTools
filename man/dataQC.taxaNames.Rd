% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataQC_Utils.R
\name{dataQC.taxaNames}
\alias{dataQC.taxaNames}
\title{perform a quality control onf taxonomic names}
\usage{
dataQC.taxaNames(taxaNames)
}
\arguments{
\item{taxaNames}{character vector. a list of scientific taxonomic names}
}
\value{
a vector with the checked taxon names
}
\description{
checks a list of taxonomic names, and perform a basic quality control
}
\details{
looks for trailing spaces or common typos in taxonomic names.
}
\seealso{
Other quality control functions: 
\code{\link{dataQC.LatitudeLongitudeCheck}()},
\code{\link{dataQC.TaxonListFromData}()},
\code{\link{dataQC.TermsCheck}()},
\code{\link{dataQC.completeTaxaNamesFromRegistery}()},
\code{\link{dataQC.dateCheck}()},
\code{\link{dataQC.eventStructure}()},
\code{\link{dataQC.findNames}()},
\code{\link{dataQC.generate.footprintWKT}()},
\code{\link{dataQC.guess.env_package.from.data}()}
}
\author{
Maxime Sweetlove CC-BY 4.0 2020
}
\concept{quality control functions}
