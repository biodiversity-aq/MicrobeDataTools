# README

## package name
working name: MicrobeDataTools

## Badges


## Package details
### discription
Since the advent of high throughput DNA sequencing techniques, 'omics datasets have become oen of the most commonly applied methods so study questions in ecology (e.g. eDNA), microbiology (e.g. amplicon and shotgun metagenome sequencing) or taxonomy (gene and genome sequencing). While these 'omics datasets are centred around DNA or RNA sequences, there are usually metadata and environmental measurement data associated with it that are vital for the correct use and interpretation of the sequence data. While standards have been develloped to deal with meta- and environmental data, such as the Minimum Information about any (x) sequence (MIxS) or DarwinCore (DwC), the implementaltion of these standards, which are optimized for data storage and computer redibility, is often challenging, significantly complicating data archiving or the re-use of data by scientists.

This package contains tools to download, format and archive genomic metadata and environmental data, with the aim to encourage long-term interoperability of data and facilitate its re-use. 

More specifically, this package focuses on the MIxS/MIMARKS standard, which is most often used in molecula biology, and the DarwinCore standard, which has become the most commonly used stadard in ecology, as well as allowing interchangebility between both standards (e.g. in case of ecological research using molecular techniques). The tools that are provided help to download data from one of the International Nucleotide Sequence Database Consortium (INSDC), including sequences (fastq-fasta format) as well as metadata and environmental data (structured using the MIxS standard). The data content can be checked for errors and for violations to the data standards using various data quality control function, including the format of date-time and geographic coordinates. User-provided data can be standardized with MIxS or DwC in a supervised way, and can be readied for upload and archiving on one of the INSDC databases. Finally, there are several tools for common manipulations in ecological molecular data, such as changing between wide and long formats,...

Link to vignetes:
  - Downloading data
  - Data quality controll
  - Data standardisation
  - upload data
  - data manipulation
  
### Background: Data standards and flavours
The MIxS and DwC data standards are implemented through a regularly updated internal library of terms, including the defintion, commonly encountered synonyms, and translation into other standards or flavours of the standard.

In theory, each standard contains a glossary of terms (in other contexts these might be called properties, elements, fields, columns or attributes), as well as some rules on it's implementation. In practice, however, variants can exist on the implenetation of the standard rules (i.e. different flavours), and non-official terms can be used when offical acceptance of a novel concept by the complete community is pending.

#### MIxS 
The MIxS has been developed by the Genomic Standards Consortium (GSC) for reporting information on nucleotide sequences (Yilmaz et al., 2011). The standard is built up of a core of common descriptors, that can be supplemented with packages of environment-specific information components.

#### DarwinCore
The Darwin Core standard was intended to facilitate the sharing of information about biological diversity, and is primarily based on taxa and their occurrences. Recently, a new structure has also been accepted, which centres around an event that can be linked to sub-events and/or occurrences. For molecular biodiversity data, this event structure allows environmental nucleotide sequence data about taxa that cannot be directly observed to fit into the DarwinCore standard.

#### flavours and implementation of data standards
In some instances, the implementation of a datastandard between platforms. For instance differences in the use of capitals, underscores, or additional "non-official" terms can complicate interoperability between platforms. In this package some of these flavours have been taken into account, including the ENA version of MIxS as well as other commonly encountered variants (read: typos) of MIxS terms.

For the internal workings of this package, we also needed to adapt the MIxS standard. For isntance, geographic coordinates are associated with the lat_lon term, but are also seperately stored in decimalLatitude and decimalLongitude due to the tendancy of some table editing software to automatically sum the values in the lat_lon column. Also because most datasets require terms from different environmental packages, all MIxS terms can be used without being confined to a single package. Finaly, to allow interoperability with the DarwinCore format, and better capture complex sampling designs, samples are structured as events (eventID) and parent events (parentEventID), analogougly to the eventCore variant of DarwinCore.

### Disclaimer
This package (version 1.0.0 2020) was specifically develloped for microbial ecological data from polar and alpine regions. This type of datasets are usually assembled from environmental samples (e.g. soil, water, tissue of a host,...) from which species identification is done using DNA sequencing methods (e.g. high throughput amplicon sequencing of the 16S rDNA gene). Often, environmental measurements on the sample or the original natural matrix are also taken to look at correlates of changes in comunity structure and composition. 
Thus, at present, the target dataset of this package have to following characteristics: 1) They follow an event structure (sample->species+other measurements); 2) they have a core or of biological information (DNA sequences, cell counts, individual species counts,...); 3) they have important associated metadata (lab protocols, DNA sequencing methods, geographic sample location,...) and 4) they can have additional environmental measurements (pH, conductivity, weather variables, chemical ions, soil composition,...)
When users have ecological datasets that somehow can't be fit into this structure, or when the package does not perform optimally on a specific dataset (e.g. not recognizing a term synonym), they are advised to mke a GIT issue or contact the package author.Upcomming versions of the package can then try to accomodate more different types of data, depending on user feedback.

## Installation instructions
  - make API key for ENA
  - install with devtools

<<<<<<< HEAD
=======
 > library(devtools)
 > devtools::install(_path_to_package_)
 

>>>>>>> ae4cff7e1564ee318f00a61a9c58db1692254f5a
## Brief demonstration usage

## Citation information
Please cite as:
Sweetlove M., Gan Y.M., Van de Putte A., 2020, MicrobeDataTools: an R package to download, format and standardize genomic meta- and environmental data.

<<<<<<< HEAD
=======


>>>>>>> ae4cff7e1564ee318f00a61a9c58db1692254f5a
