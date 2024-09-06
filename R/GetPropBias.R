library(limma)
library(WeightIt)

getPropBiasGroup <- function(Barbie, mydesign = NULL, group_vector = NULL, myblock = NULL, mytargets = NULL){

  # # targets
  # mytargets <- data.frame(Barbie$metadata, category = group_vector)
  #
  # # model design matrix
  # mydesign <- mytargets %>%
  #   with(model.matrix(~0 + category + treat + mouse + tissue))

  # make design matrix full rank
  mydesign2 <- WeightIt::make_full_rank(mydesign, with.intercept = FALSE) # delete column with nested factors

  # check if mydesign matrix is full rank or not.
  is_full_rank2 <- qr(mydesign2)$rank == min(dim(mydesign2))
  print(paste("Design matrix full rank?", is_full_rank2))

  # set up groups for correction of sample duplication
  # myblock <- mytargets %>%
  #   with(paste(treat, mouse, tissue, category))

  # arcsin squreroot transformation for data
  mydata <- Barbie$prop %>% sqrt() %>% asin()

  # perform duplication correction
  dup <- limma::duplicateCorrelation(object = mydata,
                                     design = mydesign2,
                                     block = myblock)
  print(paste("Consensus correlation: ", dup$consensus.correlation))


  # set up contrast
  mycontrasts <- mytargets %>%
    with(limma::makeContrasts(Bias = groupgroup1 - groupgroup2,
                              levels = colnames(mydesign2)))

  # fit limma model
  myfit <- limma::lmFit(object = mydata,
                        design = mydesign2,
                        block = myblock,
                        correlation = dup$consensus.correlation)

  myfit2 <- limma::contrasts.fit(fit = myfit, contrasts = mycontrasts)
  myfit3 <- limma::eBayes(myfit2)

  # apply test
  myresults <- limma::decideTests(myfit3)

  Group <- dplyr::recode(as.vector(myresults), "1" = "group1", "-1" = "group2", "0" = "Unbiased")
  names(Group) <- rownames(myresults)

  # Extract the moderated p-values
  moderated_pvalues <- myfit3$p.value %>% as.vector()

  Barbie$Bias_Prop <- data.frame(group = Group,
                                 pvalue = moderated_pvalues)

  # # Extract adjusted p-values using topTable
  top_bars <- limma::topTable(myfit3, number = Inf)  # Set number to Inf to get all barcodes
  adjusted_pvalues <- top_bars[rownames(Barbie$Bias_Prop),]$adj.P.Val

  Barbie$Bias_Prop <- data.frame(group = Group,
                                 pvalue = adjusted_pvalues)

  if(is.null(Barbie$color_panel$bias_group)) {
    Barbie$color_panel$bias_group <- c("group1" = "#33AAFF", "group2" = "#FF5959", "Unbiased" = "#FFC000")
  }

  return(Barbie)
}
