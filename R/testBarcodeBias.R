#' Test Barcodes bias across sample groups
#'
#' @param Barbie a Barbie object created by createBarbie()
#' @param method a string, options: "diffProp" and "diffOcc", choosing test method
#' @param targets a matrix or data.frame storing sample conditions
#' @param groupBy a string, or a vector of sample conditions indicating the primary effector to be tested
#' @param contrastLevels a charactor vector indicating the levels in specified groupBy
#' @param designFormula a formula to compute the designMatrix
#' @param designMatrix a numaric matrix to standardize 'targets'
#' @param block a vector of indicating sample duplicates
#' @param regularization a string, options: "firth" and "none", choosing regularization method in "diffOcc" test
#'
#' @return an updated Barbie object storing test results
#' @export
#'
#' @import limma
#' @import magrittr
#' @import logistf
#'
#' @examples
#' Block <- c(1,1,2,3,3,4,1,1,2,3,3,4)
#' Treat <- factor(rep(1:2, each=6))
#' Time <- rep(rep(1:2, each=3), 2)
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- matrix(rnorm(nbarcodes*nsamples), nbarcodes, nsamples) %>% abs()
#' rownames(y) <- paste0("Barcode", 1:nbarcodes)
#' Barbie <- Barbie::createBarbie(count, data.frame(Treat=Treat, Time=Time))

# myblock <- mytargets %>%
#   with(paste(treat, mouse, tissue, group))
testBarcodeBias <- function(Barbie, method="diffProp",
                            targets=NULL, groupBy=NULL, contrastLevels=NULL,
                            designFormula=NULL, designMatrix=NULL,
                            block=NULL, regularization="firth"){
  ## check Barbie dimensions
  if(!checkBarbieDimensions(Barbie))
    stop("Barbie components are not in right format or dimensions don't match.
         use createBarbie() and other functions in the 'Barbie' package to modify the object and don't do it manually.")
  ## check method: confirm method is chosen from "diffProp" and "diffOcc"
  method <- match.arg(method, c("diffProp", "diffOcc"))
  ## check targets: if 'targets' is not specified, assign Barbie$metadata (could still be NULL)
  if(is.null(targets)) targets <- Barbie$metadata
  ## case when targets is specified or provided by Barbie$metadata in right format
  if(is.vector(targets) || is.factor(targets)) {
    targets <- data.frame(V1=targets)
    } else if(is.matrix(targets) || is.data.frame(targets)) {
      if(nrow(targets) != ncol(Barbie$assay))
        stop("the row dimension of 'targets' doesn't match the column dimension (sample size) of 'Barbie$assay'.")
      } else {
    ## case when targets is still NULL or not in right format
    ## if group is a vector or factor of correct length, add it to targets
    if(is.vector(groupBy) || is.factor(groupBy)) {
      ## check groupBy length
      if(length(groupBy) != ncol(Barbie$assay)) {
        stop("the length of 'groupBy' doesn't match the column dimention (sample size) of 'Barbie$assay'.")
      } else {
        mytargets <- data.frame(groupBy=groupBy)
        pointer <- which(colnames(mytargets) == "groupBy")
        message("adding groupBy to targets.")}
    } else {
      stop("target not properly specified; Barbie$metadata not provided; groupBy not properly specified.
           at least one of them is needed in right format.")}
      }
  ## now targets should be a matrix or data.frame already
  ## check groupBy: if 'groupBy' is a specified effector name, extract the entire vector
  if(is.character(groupBy)) {
    if(groupBy %in% colnames(targets)) {
      groupBy <- targets[,groupBy]
      mytargets <- targets
      pointer <- which(colnames(mytargets) == groupBy)
      message("found", groupBy, "as an effector in targets or Barbie$metadata.")
    } else {stop("the groupBy specified is a charactor value,
                 but it's not an effector name found in targets or Barbie$metadata.
                 make sure you spell it correctly.")}
  } else if(is.vector(groupBy) || is.factor(groupBy)) {
    if(length(groupBy) != ncol(Barbie$assay)) {
      stop("the length of 'groupBy' doesn't match the column dimention (sample size) of 'targets' or'Barbie$assay'.")
    } else {
      mytargets <- data.frame(groupBy=groupBy, targets)
      pointer <- which(colnames(mytargets) == "groupBy")
      message("binding 'groupBy' to 'targets'.")
    }
  } else {
    groupBy <- rep(1, ncol(Barbie$assay)) %>% factor()
    mytargets <- data.frame(groupBy=groupBy, targets)
    pointer <- which(colnames(mytargets) == "groupBy")
    message("no properly specified 'groupBy'. setting samples by homogenenous group.")
  }

  ## confirm all effectors (columns) in 'mytargets' are factor() or numeric()
  ## convert columns that are neither factor nor numeric into factor
  nonFac <- sapply(mytargets, function(x) !(is.factor(x) | is.numeric(x)))
  for(col in seq(nonFac)[nonFac]) {
    mytargets[,col] <- factor(mytargets[,col])
  }

  ## import 'design' using tidy evaluation
  ## if designFormula not specified, taking all effectors from 'mytargets' into account
  if (is.null(designFormula)) {
    designFormula <- paste("~0", paste0("+ ", colnames(mytargets), collapse = " ")) %>%
      as.formula()
  }
  ## check designFormula format
  if (!inherits(designFormula, "formula")) {
    stop("The 'designFormula' argument must be a valid formula.")
  } else {
    ## check if all variables in designFormula are present in 'mytargets'
    missingTerms <- setdiff(all.vars(designFormula), colnames(mytargets))
    if (length(missingTerms) > 0) {
      stop("The following variables in the 'designFormula' are missing from 'targets' or 'Barbie$metadata':",
           paste(missingTerms, collapse = ", "))
      }
    }

  ## if designMatrix not specified, generate it by designFormula
  if (is.null(designMatrix)) {
    designMatrix <- model.matrix(designFormula, data=mytargets)
  } else {
    ## check designMatrix format and dimension
    if (is.matrix(designMatrix) || is.data.frame(designMatrix)) {
      if (nrow(designMatrix) != nrow(mytargets)) {
        stop("row dimension of 'designMatrix' doean't match row dimension of 'targets' or 'Barbie$metadata'.")
      }
    } else {
      stop("'designMatrix' should always be a matrix.
           use fucntion model.matrix() to generate a 'designMatrix'.")
    }
    }

  ## make designMatrix full rank by deleting columns of nested effectors, ie. linearly related vectors
  ## compute QR decomposition of the designMatrix
  q <- qr(designMatrix)
  keep <- rep(TRUE, ncol(designMatrix))
  ## select the indices in the pivot vector after the rank of the matrix
  ## the columns of matrix that are linearly dependent (those that do not contribute to the rank)
  keep[q$pivot[-seq(q$rank)]] <- FALSE
  designMatrix <- designMatrix[,keep, drop=FALSE]
  ## message the users if any linearly related vectors are deleted
  if(any(!keep))
    message(sum(!keep), " nested effector(s) in the designMatrix were deleted because designMatrix must be full rank.")

  ## check block groups if it's specified
  if(!(is.null(block))) {
    if(length(block) != ncol(Barbie$assay))
      stop("the length of 'block' doesn't match the row dimention (sample size) of specified 'targets' or 'Barbie$metadata'.")
  }

  ## 'pointer' indicates which column relates to 'groupBy': either a imported 'groupBy' column or a column name specified by 'groupBy' like 'Treat'
  ## case when groupBy column is factor
  if(is.factor(mytargets[, pointer])) {
    ## make group contrast, extract levels if not specified
    if(is.null(contrastLevels)) {
      contrastLevels <- levels(mytargets[, pointer])
    } else if (is.vector(contrastLevels)){
      missingLevels <- setdiff(contrastLevels, levels(mytargets[, pointer]))
      if(length(missingLevels) > 0) {
        stop("'contrastLevels' constains levels: ", missingLevels, " missing from the 'groupBy' column in 'targets' or 'Barbie$metadata'.")
      }
    } else {
      stop("'contrastLevels' argument should be a vector indicating levels in the 'groupBy' column in 'targets' or 'Barbie$metadata'.")
    }
    ## now 'contrastLevels' should be a vector indicating levels in the 'groupBy' column
    ## 'contrastLevels' has one, two, or several levels.
    if(length(contrastLevels) == 2L) {
      ## create contrast for the first two levels of 'groupBy'
      groupTitle <- colnames(mytargets)[pointer]
      contrastFormula <- paste0(groupTitle, contrastLevels[2], " - ", groupTitle, contrastLevels[1])
      ## generate contrast for designMatrix
      mycontrasts <- limma::makeContrasts(contrasts = contrastFormula,
                                          levels = colnames(designMatrix))
    } else {
      stop("'contrastLevels' has ", length(contrastLevels), " level(s), it should contain two levels for the 'groupBy' factor in 'targets' or 'Barbie$metadata'.")
    }

  } else if(is.numeric(mytargets[, pointer])) {
    ## case when groupBy column is numeric
    ## generate contrast for designMatrix
    mycontrasts <- limma::makeContrasts(contrasts = colnames(mytargets)[pointer],
                                        levels = colnames(designMatrix))
  }

  ## dispatch test functions based on the specified method
  ## default setting is "diffProp"
  if(method == "diffProp") {
    Barbie <- testDiffProp(Barbie = Barbie, mycontrasts = mycontrasts,
                           contrastLevels = contrastLevels, designMatrix = designMatrix,
                           block = block)
  } else if(method == "diffOcc") {
    ## logistic regression, default regularization is "firth"
    Barbie <- testDiffOcc(Barbie, regularization="firth",
                          mycontrasts = mycontrasts, contrastLevels = contrastLevels,
                          designMatrix = designMatrix)
  } else {stop("please choose test method from 'diffProp' or 'diffOcc'.")}

  ## assign colors for the test results
  if(is.null(Barbie$factorColors[[colnames(mytargets)[pointer]]])) {
    Barbie$factorColors[[colnames(mytargets)[pointer]]] <- setNames(
      c("#33AAFF", "#FF5959", "#FFC000"),
      c(contrastLevels[1], contrastLevels[2], "n.s."))
  }

  ## visualize test results by Heatmap and dotplots

  return(Barbie)
}
