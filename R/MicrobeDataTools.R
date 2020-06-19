#==============================================================
# The MicrobeDataTools package
#       MicrobeDataTools is a collection data management tools for microbial 'omics datasets.
#       They allow to download, structure, quqlity-controll and standardize microbial datasets
#==============================================================
# Author Maxime Sweetlove
# lisence CC 4.0
# Part of the POLA3R website (successor or mARS.biodiversity.aq)
# version 1.0 (2020-01-28)
# file encdong UTF-8
#
#==============================================================

#' MicrobeDataTools: A package with tools to format and standardize microbial 'omics datasets.
#' 
#' This package provides 5 cathegoties of tools, these are:
#' formating functions, standardization functions, quality control functions, data archiving functions and downloading data functions
#' In addition, there are also different libraries with terms of the MIxS and DarwinCore standards, term variants, synonyms and translations.
#' 
#' @section Classes:
#' MIxS.metadata
#' DwC.event
#' DwC.occurrence
#' 
#' @section Libraries:
#' TaxIDLib
#' TermsLib
#' TermsSyn
#' TermsSyn_DwC
#' MarsLib
#' ENA_allowed_terms
#' 
#' @section Functions:
#' check.valid.metadata.DwC
#' check.valid.metadata.MIxS  
#' combine.data
#' combine.data.frame
#' commonTax.to.NCBI.TaxID
#' coordinate.to.decimal
#' dataQC.completeTaxaNamesFromRegistery
#' dataQC.dateCheck
#' dataQC.DwC
#' dataQC.DwC_general
#' dataQC.eventStructure
#' dataQC.findNames
#' dataQC.generate.footprintWKT
#' dataQC.guess.env_package.from.data
#' dataQC.LatitudeLongitudeCheck
#' dataQC.MIxS
#' dataQC.taxaNames
#' dataQC.TaxonListFromData
#' dataQC.TermsCheck
#' download.sequences.INSDC
#' eMoF.to.wideTable
#' ENA_checklistAccession
#' ENA_geoloc
#' ENA_instrument
#' ENA_select
#' ENA_strat
#' FileNames.to.Table
#' find.dataset
#' get.BioProject.metadata.INSDC
#' get.boundingBox
#' get.ENAName
#' get.insertSize
#' get.sample.attributes.INSDC
#' multi.warnings
#' prep.metadata.ENA
#' rename.sequenceFiles
#' sync.metadata.sequenceFiles
#' term.definition
#' wideTab.to.hierarchicalTab
#' wideTable.to.eMoF
#'
#' @docType package
#' @name MicrobeDataTools
NULL

#' function to navigate through package
#' @export
MicrobeDataTools.help <- function(){
  message(paste(ls("package:MicrobeDataTools"), collapse="\n"))
}

#' Find the MIxS or DarwinCore standard term and definition of a variable
#' @author Maxime Sweetlove ccBY 4.0 2019
#' @description retrieve the MIxS or DarwinCore standard term and definition of a variable.
#' @param term a character string. The variable to look for among the MIxS and DarwinCore vocabularies
#' @details Standerdizing microbial sequence data, metadata and environmental data can be quite difficult given the plethora of standard terms already in existance. This function returns a definition of any term that is used on the POLAAAR portal at biodiversity.aq.
#' @return chracater string printed to the console. The best matching terms and their the definitions.
#' @export
term.definition <- function(term){
  if(term %in% TermsLib$name){
    def_out <- as.character(TermsLib[TermsLib$name==term,]$definition)
    out_message <- paste(term, "\t\n", def_out)
  }else{
    potential_match <- c()
    def_out <- c()
    for(ls in TermsSyn){
      if(grepl(term, ls)){
        lst <- names(TermsSyn[ls[1]])
        potential_match <- c(potential_match, lst)
        def_out <- c(def_out, as.character(TermsLib[TermsLib$name==lst,]$definition))
      }
    }
    
    if(length(potential_match)>1){
      out_message <- paste("Multiple matches found for \"" , term, "\"\n\n",
                           "\t", paste(c(1:length(potential_match)),
                                       rep(". ",length(potential_match)),
                                       potential_match,
                                       rep(" :\n",length(potential_match)),
                                       rep("\"",length(potential_match)),
                                       def_out,
                                       rep("\"",length(potential_match)),
                                       collapse="\n\t", sep=""), sep="")
    }else if(length(potential_match)==1){
      out_message <- paste("Best match for \"" , term, "\"\n\n",
                           potential_match, "\n", def_out, sep="")
    }else{
      out_message <-  paste("Could not find any matches for \"" , term, "\"\n",sep="")
    }
  }
  message(out_message)
  
}
