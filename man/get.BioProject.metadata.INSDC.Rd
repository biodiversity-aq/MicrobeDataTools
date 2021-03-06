% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Get_SequenceData_INSDC.R
\name{get.BioProject.metadata.INSDC}
\alias{get.BioProject.metadata.INSDC}
\title{Fetch BioProject metadata}
\usage{
get.BioProject.metadata.INSDC(BioPrjct, just.names=FALSE)
}
\arguments{
\item{BioPrjct}{A chracter string. A single Bioproject ID, e.g. "PRJNA369175"}

\item{just.names}{Boolean. If TRUE, only the INSDC sample names are returned, else all data are returned. default FALSE}
}
\value{
if get.BioProject.metadata.INSDC(just.names=FALSE) (default) a data.frame with n rows an m columns is returned, n being the number of samples in the BioProject, and m being the number of variables found. If get.BioProject.metadata.INSDC(just.names=TRUE), a character vector of length n is returned with the sample numbers ("Run numbers", SRR numbers)
}
\description{
downloads a minimal set of sample metadata of all samples ("Runs") within a BioProject (argument BioPrjct) from the International Nucleotide Sequence Database Consortium (INSDC) databases.
}
\details{
BioProjects combine all biological nucleotide sequence data related to a single initiative, originating from a single organization. With each sample ("Run") within a BioProject, there is additional data associated that is crucial to the correct interpretation of the nucleotide sequence data, but is not automatically downloaded along with it. The get.BioProject.metadata.INSDC function will fetch the most basic metadata of a BioProject from the INSDC repositories to complete the nucleotide sequence dataset, using E-utils API function of NCBI. These basic metadata typically include Run number, relsease date, load date, spots, bases, av_MB and download path. Note that the data returned by get.BioProject.metadata.INSDC does not include all the metadata associated with a BioProject. Other information, like coordinates, sampling dates or environmental measurements may also be available, but require the user to register at NCBI and request a personal API-key (this is required by NCBI to acces their data since 2017). The complete set of additional data can be downloaded using the get.sample.attributes.INSDC function, given a user-specified API-key. see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
}
\references{
Sayers, E. (2009) The E-utilities In-Depth: Parameters, Syntax and More, https://www.ncbi.nlm.nih.gov/books/NBK25499/
}
\seealso{
get.sample.attributes.INSDC

Other downloading data functions: 
\code{\link{download.sequences.INSDC}()},
\code{\link{get.sample.attributes.INSDC}()}
}
\author{
Maxime Sweetlove ccBY 4.0 2019
}
\concept{downloading data functions}
