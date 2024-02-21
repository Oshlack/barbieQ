## ----dummy, include=FALSE-----------------------------------------------------
# dummy data
data <- data.frame(
  METRIC = paste('var', 1:10),
  VAL = runif(10),
  stringsAsFactors = FALSE
)

## ----results = 'asis', include=FALSE------------------------------------------
# automatically generate <style> block of css for each tab
tab_number <- 1:nrow(data)
tab_color <- "#CCCCCC"
css <- sprintf(".color-tabs>.nav-pills>li:nth-child(%d){background:%s;}", tab_number, tab_color)
css <- paste(css, collapse = "\n")
css <- paste("<style>", css, "</style>")
cat(css)

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = "#>")

## ----load library, results='hide'---------------------------------------------
library(here)
library(magrittr)
library(knitr)
library(grid)
library(tidyverse)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggbreak)
library(patchwork)
library(plotly)
library(ComplexHeatmap)
library(eulerr)
library(colorRamp2)
library(igraph)
library(htmlwidgets)
library(speckle)
library(limma)
library(edgeR)
library(pracma)
library(WeightIt) # make full rank

source(here::here("R", "Barbie_object.R")) # create Barbie object
source(here::here("R", "Pair_Correlation.R")) # predict correlating barcodes
source(here::here("R", "Pareto_contribution.R")) # visualize barcode contribution
source(here::here("R", "Sankey_contribution.R"))
source(here::here("R", "Bar_contribution.R"))
source(here::here("R", "Get_OccBias.R")) # test barcode occurrence bias
source(here::here("R", "Get_PropBias0.R")) # test barcode proportion bias
source(here::here("R", "CountHP0.R")) # plot complexheatmap


## ----read data----------------------------------------------------------------
counts <- read.csv(here::here("inst", "extdata", "barcodes_example.csv"), row.names = 1) 
# make sure "counts" is a numeric count table with barcodes in rows and samples in columns

targets <- read.table(here::here("inst", "extdata", "targets_example.txt"))
# make sure "targets" is a table saving experiment designs with samples in rows and factors/conditions in columns

colors <- readRDS(here::here("inst", "extdata", "colors_example.RDS"))
# "colors" is a list of color palette you designed for different conditions of samples.

## ----create barbie------------------------------------------------------------
example_bb <- createBarbie(counts = counts, 
                           metadata = targets, 
                           color_panel = colors) # no need to pass color_panel if you don't have it.

## ----pre-filter, include=FALSE------------------------------------------------
# trim rows
flag <- rowSums(example_bb$presence) >= 2
example_bb <- trimRow(Barbie = example_bb, 
                      keep_rows = flag)

# trim unwanted samples
example_bb <- trimObjectByMetadata(Barbie = example_bb, 
                                   condition = "time", 
                                   specified = c("week4","week8")) # or use function: trimColumn

# YOUR BARBIE <- trimColumn(Barbie = YOUR BARBIE,
#                           keep_columns = YOUR ARRAY)

# select top contributing barcodes
example_bb <- getTopBar(Barbie = example_bb)

# trim rows by "is_top"
example_bb <- trimRow(Barbie = example_bb, 
                      keep_rows = example_bb$is_top)

## ----collapse correlating barcodes, fig.width=6, fig.height=4.5---------------
# group.vec is meant to pass the "ground truth" of known correlating barcodes from the same clones
plot_bar_pair_cor <- plotBarcodePairCor_Count(count.mx = example_bb$assay, # this is your barcode count table (make sure barcodes in rows)
                                              # group.vec = ADD YOUR OWN GROUP OR LEAVE IT AS NULL, # this is a vector of known barcode groups
                                              plot = "mean") # options for plotting, including "mean", "max", "histogram", and "table"

# plot the pairwise correlation of barcodes
plot_bar_pair_cor 
# if you prefer interactive plot, go for this:
# ggplotly(plot_bar_pair_cor)

# save the results of pairwise correlation
result_bar_pair_cor <- plotBarcodePairCor_Count(count.mx = example_bb$assay, 
                                                # group.vec = ADD YOUR OWN GROUP OR LEAVE IT AS NULL, 
                                                plot = "table")

# select barcode pairs that are previously undefined, highly correlating, and making big contribution.
req1 <- result_bar_pair_cor$Cor > 0.95 # pair correlation > 0.95
req2 <- result_bar_pair_cor$Mean > 2^14 # pair mean count > 2^14
req3 <- result_bar_pair_cor$group == "undefined barcode pair"

proposed_bar_pair <- result_bar_pair_cor[req1 & req2 & req3,]
kable(head(proposed_bar_pair))

proposed_bar_pair_list <- lapply(proposed_bar_pair$barID, function(x) strsplit(x, split = "\\.") %>% unlist())

barcode_clone_ref <- createBarcodeCloneRef(new_pairlist = proposed_bar_pair_list)

# If you have an old barcode-clone reference list, go for this:
# updateBarcodeCloneRef(old_reflist = OLD LIST,
#                       new_pairlist = proposed_bar_pair_list)

## ----recreate barbie----------------------------------------------------------

example_bb <- readRDS(here::here("inst", "extdata", "example_bb.RDS"))

collapsed_bb <- CollapseRow(Barbie = example_bb,
                            group_array = example_bb$clone_group)

# re-create new Barbie Object, because CPM needs to be recalculated.
example_bb <- createBarbie(counts = collapsed_bb$assay,
                           metadata = collapsed_bb$metadata,
                           color_panel = collapsed_bb$color_panel)


## ----filter, fig.width=5, fig.height=4----------------------------------------
# trim rows
flag <- rowSums(example_bb$presence) >= 2
example_bb <- trimRow(Barbie = example_bb, 
                      keep_rows = flag)

# trim unwanted samples
# example_bb <- trimObjectByMetadata(Barbie = example_bb, 
#                                    condition = "time", 
#                                    specified = c("week4","week8")) # or use function: trimColumn

# YOUR BARBIE <- trimColumn(Barbie = YOUR BARBIE,
#                           keep_columns = YOUR ARRAY)

# Hi Enid, if you let specified = "week4" as the following code, you will get the same results as I got for week4 data.
week4_bb <- trimObjectByMetadata(
  Barbie = example_bb,
  condition = "time",
  specified = "week4"
  )

# select top contributing barcodes
week4_bb <- getTopBar(Barbie = week4_bb)

# plot barcode contribution
PlotCircularContribution(Barbie = week4_bb)

PlotTotalContribution(Barbie = week4_bb)

# trim rows by "is_top"
week4_top <- trimRow(Barbie = week4_bb, 
                      keep_rows = week4_bb$is_top)

# plot top barcode contribution
PlotBarContribution(Barbie = week4_top)

## ----collapse samples---------------------------------------------------------
# create a group vector that identify samples of different conditions. 
# samples of PCR reps are in the same group in this vector.
vector_PCRrep <- week4_top$metadata %>% 
  with(
    paste(mouse, tissue, celltype, sep = ".")
  )

# collapse the columns in Barbie Object by the group vector.
# count / prop / PCM data takes the average, presence data takes the max.
week4_top_coll <- CollapseColumn(Barbie = week4_top, group_array = vector_PCRrep)

## ----contingency table bias test----------------------------------------------
# trim samples as you need. 
week4_test <- trimObjectByMetadata(Barbie = week4_top_coll, 
                                    condition = "treat", specified = "IT", 
                                    keep = FALSE) # exclude "IT" samples
week4_test <- trimObjectByMetadata(Barbie = week4_test, 
                                    condition = "lineage", specified = "immature", 
                                    keep = FALSE) # exclude "immature" celltype from samples

# customize sample groups as you need.
Vector_customized <- week4_test$metadata$lineage
Vector_customized[Vector_customized %in% c("Bcell", "Tcell")] <- "Lymphoid"# group samples

# Generate tables and create the list
c_tables <- GetContingencyTable(week4_test, Vector_customized = Vector_customized)
#print the first 5 contingency table
lapply(c_tables[1:5], function(x) {knitr::kable(x)})

## ----apply bias test----------------------------------------------------------
# Apply test, and get Bias group
week4_test <- GetOccBiasGroup(Barbie = week4_test, contingency_table_ls = c_tables)

kableExtra::kable(week4_test$Bias_Occ)

## ----visualize occ bias test, fig.width=6, fig.height=4-----------------------
# Visualize Bias group of barcodes
PlotBiasVsRank(Barbie = week4_test, passing_data = "contribution")

PlotBiasVsRank(Barbie = week4_test, passing_data = "output")

PlotBiasVsRank(Barbie = week4_test, passing_data = "rank")

PlotBiasVsRank(Barbie = week4_test, passing_data = "rank_var")

PlotCpmHP_0(Barbie = week4_test, show_bias = TRUE, Vector_customized = Vector_customized, split_by = Vector_customized) 

PlotPreHP_0(Barbie = week4_test, show_bias = TRUE, Vector_customized = Vector_customized, split_by = Vector_customized) 

## ----design matrix for porp bias test-----------------------------------------
# trim samples as you need. 
week4_test <- trimObjectByMetadata(Barbie = week4_top, 
                                    condition = "treat", specified = "IT", 
                                    keep = FALSE) # exclude "IT" samples
week4_test <- trimObjectByMetadata(Barbie = week4_test, 
                                    condition = "lineage", specified = "immature", 
                                    keep = FALSE) # exclude "immature" celltype from samples

# customize sample groups as you need for comparison in the following test.
Vector_customized <- week4_test$metadata$lineage
Vector_customized[Vector_customized %in% c("Bcell", "Tcell")] <- "Lymphoid"# group samples

# make your own targets and design matrix
  # targets
  mytargets <- data.frame(week4_test$metadata, category = Vector_customized)

  # model design matrix
  mydesign <- mytargets %>%
    with(model.matrix(~0 + category + treat + mouse + tissue))

# apply test for Bias_Prop
week4_test <- getPropBiasGroup(Barbie = week4_test, mydesign = mydesign)

kableExtra::kable(week4_test$Bias_Prop)

## ----visualize prop bias test, fig.width=6, fig.height=4----------------------
# Visualize Bias group of barcodes
PlotBiasVsRank(Barbie = week4_test, passing_data = "contribution", 
               bias_group = week4_test$Bias_Prop$group,
               bias_pvalue = week4_test$Bias_Prop$pvalue)

PlotBiasVsRank(Barbie = week4_test, passing_data = "output", 
               bias_group = week4_test$Bias_Prop$group,
               bias_pvalue = week4_test$Bias_Prop$pvalue)

PlotBiasVsRank(Barbie = week4_test, passing_data = "rank", 
               bias_group = week4_test$Bias_Prop$group,
               bias_pvalue = week4_test$Bias_Prop$pvalue)

PlotBiasVsRank(Barbie = week4_test, passing_data = "rank_var", 
               bias_group = week4_test$Bias_Prop$group,
               bias_pvalue = week4_test$Bias_Prop$pvalue)

PlotCpmHP_0(Barbie = week4_test, show_bias = TRUE, 
            Vector_customized = Vector_customized, split_by = Vector_customized, 
            bias_group = week4_test$Bias_Prop$group) 

PlotPreHP_0(Barbie = week4_test, show_bias = TRUE, 
            Vector_customized = Vector_customized, split_by = Vector_customized,
            bias_group = week4_test$Bias_Prop$group) 

