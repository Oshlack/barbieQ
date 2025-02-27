#' Monkey HSPC Cell Barcoding Data
#'
#' A subset of data from a study on monkey HSPC cell expansion using barcoding technique.
#' Unique barcodes were initially integrated into hematopoietic stem and progenitor cells (HSPCs)  
#' and subsequently passed to daughter cells.
#' After a defined period of expansion, offspring cells were collected and sorted into various cell types. 
#' Barcode counts within different cell types were used to interpret the patterns of HSPC differentiation.
#' 
#' The original monkey HSPC barcoding data were published in the following paper:
#'   [Wu, Chuanfeng, et al. "Clonal expansion and compartmentalized maintenance of rhesus macaque NK cell subsets." Science Immunology (2018)](http://dx.doi.org/10.1126/sciimmunol.aat9781)
#' 
#' However, the barcode count data were analysed by the original authors using the [barcodetrackR](10.18129/B9.bioc.barcodetrackR) package 
#' and made available through the compatible [barcodetrackRData](https://github.com/dunbarlabNIH/barcodetrackRData) repository on GitHub.
#' 
#' Here, we directly source the raw data from `barcodetrackRData`.
#' 
#' This dataset includes data from the "ZG66" monkey only. Additional datasets are available at the source link.
#'
#' @format ## `monkeyHSPC`
#' A list containing:
#' - A barcode count matrix with 16,603 rows and 62 columns,
#' - A a data frame of sample metadata.
#' - other components associated with this `barbieQ` package.
#' 
#' \describe{
#'   \item{Monkey}{Individual monkey ID}
#'   \item{Months}{Time (in months) after cell expansion}
#'   \item{Celltype}{Collected cell type}
#'   ...
#' }
#' @source <https://github.com/dunbarlabNIH/barcodetrackRData/tree/main/datasets/WuC_etal>
"monkeyHSPC"
#' The following code is used to process the raw data into the package exported data

## read raw barcode count matrix.
count_path <- system.file("extdata/monkey_ZG66.txt", package = "barbieQ")
## data is separated by "\t", so use read.table()
count <- read.table(count_path, header = TRUE, sep = "\t", row.names = 1)

## rownames are barcode sequence, basically consist of ATGCs: "^[ATGC]+$"
# lapply(rownames(count), function(x) grepl("^[ATGC]+$", x)) %>% unlist() %>% table()

## read metadata of samples.
metadata_path <- system.file("extdata/monkey_ZG66_metadata.txt", package = "barbieQ")
metadata <- read.table(metadata_path, header = 1, row.names = 1)

## update the variable names following camelCase style.
colnames(metadata)[grepl("Cell_type", colnames(metadata))] <- "Celltype"
# print(colnames(metadata))

## Timepoint and Months variables in the `metadata` are simply repetition, delete one
metadata$Timepoint <- NULL

## create an object to store the barcode count matrix and sample metadata.
monkeyHSPC <- barbieQ::createBarbieQ(object = count, target = metadata)

## save the package exported data.
usethis::use_data(monkeyHSPC, overwrite = TRUE)

