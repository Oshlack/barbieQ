#' Plot Barcode contributions as mean Barcode proportion across samples
#'
#' `plotBarcodeProportion()` visualizes the overall proportion of Barcodes
#'  within groups of significance, as determined by [testBarcodeBias], for each sample.
#'
#' @param barbieQ A `SummarizedExperiment` object created by the [createBarbieQ] function,
#'  updated with Barcode test results by calling the [testBarcodeBias] function.
#'
#' @return A `ggplot` S3 class object displaying Barcode contributions in a
#'  bar plot.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom circlize colorRamp2
#' @import data.table
#' @importFrom magrittr %>% 
#' @importFrom dplyr  group_by
#' @importFrom dplyr  summarise
#' @importFrom dplyr  mutate
#' @importFrom tidyr  pivot_longer
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors  metadata
#'
#' @examples
#' ## sample conditions and color palettes
#' sampleConditions <- data.frame(
#'   Treat = factor(rep(c('ctrl', 'drug'), each = 6)),
#'   Time = rep(rep(seq_len(2), each = 3), 2)
#' )
#' conditionColor <- list(
#'   Treat = c(ctrl = '#999999', drug = '#112233'),
#'   Time = c('1' = '#778899', '2' = '#998877')
#' )
#' ## Barcode count data
#' nbarcodes <- 50
#' nsamples <- 12
#' barcodeCount <- abs(matrix(10, nbarcodes, nsamples))
#' barcodeCount[seq(21, 50), ] <- 0.0001
#' rownames(barcodeCount) <- paste0('Barcode', seq_len(nbarcodes))
#' ## create a `barbieQ` object
#' myBarbieQ <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
#' plotSignifBarcodeProportion(myBarbieQ)
plotSignifBarcodeProportion <- function(barbieQ) {
  ## extract testing results and information
  statsDf <- SummarizedExperiment::rowData(barbieQ)$testingBarcode
  design <- S4Vectors::metadata(statsDf)$design
  contrastGroups <- S4Vectors::metadata(statsDf)$contrastGroups
  method <- S4Vectors::metadata(statsDf)$method
  
  
  ## extract color code
  colorCode <- S4Vectors::metadata(barbieQ)$factorColors$testingBarcode
  
  ## extract sample groups
  if(all(contrastGroups %in% colnames(design))) {
    ## when contrast is factor: two levels
    Group <- design[, contrastGroups]
    GroupVec <- Group[,contrastGroups["levelHigh"]] - Group[,contrastGroups["levelLow"]]
    GroupVec[GroupVec == -1] <- contrastGroups["levelLow"]
    GroupVec[GroupVec == 1] <- contrastGroups["levelHigh"]
  } else {
    ## when contrast is numeric
    if(method == "diffOcc") {
      contrastVariables <- S4Vectors::metadata(statsDf)$contrastVariables
      primaryVar <- contrastVariables$primaryVar
      GroupVec <- design[,primaryVar]
    } else if(method == "diffProp") {
      contrasts <- S4Vectors::metadata(statsDf)$contrasts
      primaryVar <- colnames(contrasts)[1]
      GroupVec <- design[,primaryVar]
    } 
  }
  
  ## extract barcode groups
  proportion <- SummarizedExperiment::assays(barbieQ)$proportion
  
  propGroupDf <- data.frame(tendencyTo = statsDf$tendencyTo, proportion)
  
  groupSums <- propGroupDf %>%
    dplyr::group_by(tendencyTo) %>%
    dplyr::summarise(across(everything(), \(x) sum(x, na.rm = TRUE)))
  
  longData <- tidyr::pivot_longer(
    groupSums, 
    cols = colnames(groupSums)[-1], 
    names_to = "sample", 
    values_to = "proportion")

  ## assign sample groups based on contrast
  longData <- longData %>% 
    dplyr::mutate(sampleGroup =  GroupVec[match(sample, names(GroupVec))])
  ## convert sample and sampleGroup into factors
  if(all(contrastGroups %in% colnames(design))) {
    ## when contrast is factor: two levels
    longData$sampleGroup <- factor(longData$sampleGroup, levels = contrastGroups)
  }
  longData$sample <- factor(longData$sample, levels = unique(longData$sample[order(longData$sampleGroup)]))
  
  p <- ggplot(longData, 
         aes(x = sample, y = proportion, fill = tendencyTo)) +
    geom_bar(stat = "identity", alpha = 0.7) +  
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, size = 4),
          aspect.ratio = 1,
          panel.grid.major.x = element_blank()) + 
    labs(x = "Sample", y = "Overall Barcode proportion within signif groups") + 
    scale_fill_manual(values = colorCode) 
  
  ## add sample annotation when group is factor
  if(all(contrastGroups %in% colnames(design))) {
    p <- p +
      geom_tile(
        data = longData,
        aes(x = sample, y = -0.02, fill = sampleGroup),
        height = 0.02, width = 1
      ) +
      geom_text(
        data = longData %>% 
          ## convert 'sample' to numeric for positions
          mutate(id = as.numeric(sample)) %>%  
          group_by(sampleGroup) %>% 
          ## calculate the mean of the x positions
          summarize(mid_point = mean(id), .groups = "drop"),  
        aes(x = mid_point, y = -0.05, label = sampleGroup),
        ## disable inheritance of global aesthetics
        size = 2, vjust = 0.5, color = "black", inherit.aes = FALSE
      )
  }
  
  return(p)
}