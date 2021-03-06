% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataQC_Utils.R
\name{dataQC.completeTaxaNamesFromRegistery}
\alias{dataQC.completeTaxaNamesFromRegistery}
\title{Complete a list of taxonomic names with information form an authorotive taxonomic backbone}
\usage{
dataQC.completeTaxaNamesFromRegistery(taxaNames, taxBackbone="gbif")
}
\arguments{
\item{taxaNames}{a character vector. A list with the taxonomic names to look for}

\item{taxBackbone}{a character string. The taxonomic backbone to querry. Either "worms" or "gbif"}
}
\value{
a dataframe with the following fields:scientificName, scientificNameID, aphID, kingdom, phylum, class, order, family, genus, specificEpithet, scientificNameAuthorship, namePublishedInYear
}
\description{
complete a list of taxonomic names by looking-up missing information on an accepted taxonomic registery
}
\details{
using the API-client connection to the World Registry of Marine Species (WORMS), additional taxonomic information can be added to an existing list of taxa
}
\seealso{
Other quality control functions: 
\code{\link{dataQC.LatitudeLongitudeCheck}()},
\code{\link{dataQC.TaxonListFromData}()},
\code{\link{dataQC.TermsCheck}()},
\code{\link{dataQC.dateCheck}()},
\code{\link{dataQC.eventStructure}()},
\code{\link{dataQC.findNames}()},
\code{\link{dataQC.generate.footprintWKT}()},
\code{\link{dataQC.guess.env_package.from.data}()},
\code{\link{dataQC.taxaNames}()}
}
\author{
Maxime Sweetlove CC-BY 4.0 2020
}
\concept{quality control functions}
