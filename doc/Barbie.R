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

source(here::here("R", "Barbie_object.R")) # create Barbie object
source(here::here("R", "Pair_Correlation.R")) # predict correlating barcodes
source(here::here("R", "Pareto_contribution.R")) # visualize barcode contribution
source(here::here("R", "Sankey_contribution.R"))
source(here::here("R", "Bar_contribution.R"))
source(here::here("R", "Get_Potential.R")) # test barcode potential (bias)
source(here::here("R", "CountHP.R")) # plot complexheatmap


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
example_bb <- trimObjectByMetadata(Barbie = example_bb, 
                                   condition = "time", 
                                   specified = c("week4","week8")) # or use function: trimColumn

# YOUR BARBIE <- trimColumn(Barbie = YOUR BARBIE,
#                           keep_columns = YOUR ARRAY)

# select top contributing barcodes
example_bb <- getTopBar(Barbie = example_bb)

# plot barcode contribution
PlotCircularContribution(Barbie = example_bb)

PlotTotalContribution(Barbie = example_bb)

# trim rows by "is_top"
example_top <- trimRow(Barbie = example_bb, 
                      keep_rows = example_bb$is_top)

# plot top barcode contribution
PlotBarContribution(Barbie = example_top)

## ----contingency table bias test----------------------------------------------
# trim samples as you need. 
example_test <- trimObjectByMetadata(Barbie = example_top, 
                                    condition = "treat", specified = "IT", 
                                    keep = FALSE) # exclude "IT" samples
example_test <- trimObjectByMetadata(Barbie = example_test, 
                                    condition = "lineage", specified = "immature", 
                                    keep = FALSE) # exclude "immature" celltype from samples

# customize sample groups as you need.
Vector_customized <- example_test$metadata$lineage
Vector_customized[Vector_customized %in% c("Bcell", "Tcell")] <- "Lymphoid"# group samples

# Generate tables and create the list
c_tables <- GetContingencyTable(example_test, Vector_customized = Vector_customized)
#print the first 5 contingency table
lapply(c_tables[1:5], function(x) {knitr::kable(x)})

## ----apply bias test----------------------------------------------------------
# Apply test, and get Bias group
example_test <- GetFisherBiasGroup(Barbie = example_test, contingency_table_ls = c_tables)

## ----visualize bias test, fig.width=5, fig.height=4---------------------------
# Visualize Bias group of barcodes
PlotBiasVsRank(Barbie = example_test, passing_data = "output")

PlotBiasVsRank(Barbie = example_test, passing_data = "rank")

PlotBiasVsRank(Barbie = example_test, passing_data = "rank_var")

PlotCpmHP_0(Barbie = example_test, show_bias = TRUE, Vector_customized = Vector_customized) 

PlotPreHP_0(Barbie = example_test, show_bias = TRUE, Vector_customized = Vector_customized) 

