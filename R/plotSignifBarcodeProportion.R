#' Plot Barcode contributions as mean Barcode proportion across samples
#'
#' `plotBarcodeProportion()` visualizes the overall proportion of Barcodes
#'  within groups of significance, as determined by [testBarcodeSignif], for each sample.
#'
#' @param barbieQ A `SummarizedExperiment` object created by the [createBarbieQ] function,
#'  updated with Barcode test results by calling the [testBarcodeSignif] function.
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
#' Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
#' Treat <- factor(rep(c('ctrl', 'drug'), each = 6))
#' Time <- rep(rep(seq_len(2), each = 3), 2)
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
#' testBB <- testBarcodeSignif(barbieQ, sampleGroup = 'Treat')
#' plotSignifBarcodeProportion(testBB)
plotSignifBarcodeProportion <- function(barbieQ) {
  ## extract testing results and information
  statsDf <- SummarizedExperiment::rowData(barbieQ)$testingBarcode
  contrastGroups <- S4Vectors::metadata(statsDf)$contrastGroups
  method <- S4Vectors::metadata(statsDf)$method
  ## extract design based on tests
  if(method == "diffProp") {
    design <- S4Vectors::metadata(statsDf)$design
  } else if (method == "diffOcc") {
    design <- S4Vectors::metadata(statsDf)$pseudo.design
  }
  
  ## extract color code
  colorCode <- S4Vectors::metadata(barbieQ)$factorColors$testingBarcode
  
  splitGroupHigh <- ""
  splitGroupLow <- ""
  ## extract relevant conditions
  if(all(names(contrastGroups) %in% c("levelLow", "levelHigh"))) {
    splitGroupHigh <- strsplit(contrastGroups["levelHigh"], " \\+ ")[[1]]
    splitGroupLow <- strsplit(contrastGroups["levelLow"], " \\+ ")[[1]]
  }
  ## extract sample groups
  if(all(splitGroupHigh %in% colnames(design)) ||
     all(splitGroupLow %in% colnames(design))) {
    ## when contrast is factor: two levels
    GroupHigh <- design[, splitGroupHigh, drop = FALSE] |> rowSums()
    GroupLow <- design[, splitGroupLow, drop = FALSE] |> rowSums()
    GroupVec <- rep("others", ncol(barbieQ))
    names(GroupVec) <- rownames(design)
    GroupVec[GroupHigh == 1] <- contrastGroups["levelHigh"]
    GroupVec[GroupLow == 1] <- contrastGroups["levelLow"]
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
  if(all(splitGroupHigh %in% colnames(design)) ||
     all(splitGroupLow %in% colnames(design))) {
    ## when contrast is factor: two levels
    longData$sampleGroup <- factor(longData$sampleGroup, levels = c(contrastGroups, "others"))
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
  if(all(splitGroupHigh %in% colnames(design)) ||
     all(splitGroupLow %in% colnames(design))) {
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