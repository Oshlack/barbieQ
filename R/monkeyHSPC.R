#' Monkey HSPC Cell Barcoding Data
#'
#' A subset of data from a study on monkey HSPC cell expansion using barcoding technique.
#' Unique barcodes were initially integrated into hematopoietic stem and progenitor cells (HSPCs)  
#' and subsequently passed to daughter cells.
#' After a defined period of cell expansion, progeny cells were collected and sorted into various cell types as individual samples. 
#' Barcode counts within different samples were used to interpret the patterns of HSPC differentiation.
#' 
#' The original monkey HSPC barcoding data were published in the following paper:
#'   [Wu, Chuanfeng, et al. 'Clonal expansion and compartmentalized maintenance of rhesus macaque NK cell subsets.' Science Immunology (2018)](http://dx.doi.org/10.1126/sciimmunol.aat9781)
#' 
#' However, the barcode count data were analysed by the original authors using the [barcodetrackR](10.18129/B9.bioc.barcodetrackR) package 
#' and made available through the compatible [barcodetrackRData](https://github.com/dunbarlabNIH/barcodetrackRData) repository on GitHub.
#' 
#' Here, we directly source the raw data from `barcodetrackRData`.
#' 
#' This dataset includes data from the 'ZG66' monkey individual only. Additional datasets are available at the source link.
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
