#' Differential occurrence test on Barcodes across sample groups
#'
#' @param occurrence The `occurrence` assay in `SummarizedExperiment::assays(barbieQ)`.
#' @param regularization A string specifying the regularization method when
#'  testing 'diffOcc'. Options: 'firth' and 'none'. Defaults to 'firth'.
#' @param mycontrasts A numeric vector generated by [limma::makeContrasts]
#' @param designMatrix A numeric matrix standardizing `targets`, generated by
#'  the [stats::model.matrix] function. Defaults to be generated by
#'  `designFormula`.
#'
#' @importFrom logistf logistf
#' @importFrom magrittr %>%
#' @importFrom stats as.formula
#' @importFrom stats model.matrix
#' @importFrom stats glm
#' @importFrom stats p.adjust
#'
#' @noRd
#'
#' @return A list containing: 
#'  * a`data.frame` of statistical test results for each Barcode;
#'  * a `matrix` indicating the design of samples;
#'  * a `vector` indicating the primary, intercept, and additional variables.
#'
#' @examples \donttest{
#' Treat <- factor(rep(c('ctrl', 'drug'), each = 6))
#' Time <- rep(rep(seq_len(2), each = 3), 2)
#' nbarcodes <- 50
#' nsamples <- 12
#' count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
#' rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
#' barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
#' barbieQ <- testDiffOcc(barbieQ,
#'   regularization = 'firth',
#'   mycontrasts = mycontrasts, contrastLevels = contrastLevels,
#'   designMatrix = designMatrix
#' )
#' barbieQ:::testDiffOcc(
#'   occurrence = SummarizedExperiment::assays(barbieQ)$occurrence,
#'   mycontrasts = setNames(c(-1,1,0), c('Treatctrl', 'Treatdrug', 'Time')),
#'   designMatrix = model.matrix(~ 0 + Treat + Time)
#' )
#' }
testDiffOcc <- function(occurrence, regularization = "firth", mycontrasts = NULL, designMatrix = NULL) {
    ## check 'regularization'
    regularization <- match.arg(regularization, c("firth", "none"))

    ## extract binary data
    mydata <- occurrence + 1 - 1

    ## for now, this function only supports simple contrasts
    if (!all(mycontrasts %in% c(-1, 0, 1))) {
        stop("testDiffOcc doesn't support complex contrasts yet; please customize the design and contrast by re-grouping your samples.")
    }

    ## regenerate designMatrix based on the mycontrasts specified identify variables
    ## based on contrast rules
    mycontrasts <- data.frame(mycontrasts)
    rownames(mycontrasts) <- colnames(designMatrix)
    interceptVar <- rownames(mycontrasts)[mycontrasts == -1]
    primaryVar <- rownames(mycontrasts)[mycontrasts == 1]
    additionalVar <- rownames(mycontrasts)[mycontrasts == 0]
    ## recreate the formula add intercept if existing
    formulaStr <- paste(primaryVar)
    ## add additional variable if existing
    formulaStr <- paste(c(formulaStr, additionalVar), collapse = " + ")

    ## recreate design without the intercept (interceptVar serves as reference)
    formulaStr <- paste0("~ ", formulaStr)
    formula <- stats::as.formula(formulaStr)
    design <- stats::model.matrix(formula, data = data.frame(designMatrix))

    ## make designMatrix full rank by deleting columns of nested effectors, ie.
    ## linearly related vectors compute QR decomposition of the designMatrix
    q <- qr(design)
    keep <- rep(TRUE, ncol(design))
    ## select the indices in the pivot vector after the rank of the matrix the columns
    ## of matrix that are linearly dependent (those that do not contribute to the
    ## rank)
    keep[q$pivot[-seq(q$rank)]] <- FALSE
    design <- design[, keep, drop = FALSE]

    ## case when 'none' regularization, fit classic model
    if (regularization == "none") {
        results <- lapply(seq_len(nrow(occurrence)), function(i) {
            stats::glm(occurrence[i, ] ~ design - 1, family = "binomial") %>%
                summary()
        })
        ## extract stats
        tag <- paste0("design", primaryVar)
        P.Value <- lapply(results, function(x) {
            x$coefficients[tag, "Pr(>|z|)"]
        }) %>%
            unlist()
        logOR <- lapply(results, function(x) {
            x$coefficients[tag, "Estimate"]
        }) %>%
            unlist()
        z <- lapply(results, function(x) {
            x$coefficients[tag, "z value"]
        }) %>%
            unlist()
        ## extract stats result sheet
        statMat <- data.frame(P.Value = P.Value, logOR = logOR, z = z)
        rownames(statMat) <- rownames(occurrence)
    } else if (regularization == "firth") {
        ## case when 'firth' regularization, fit penalized logistic model
        results <- lapply(seq_len(nrow(occurrence)), function(i) {
            invisible(logistf::logistf(occurrence[i, ] ~ design - 1))
        })
        tag <- paste0("design", primaryVar)
        ## extract stats
        P.Value <- lapply(results, function(x) {
            x$prob[tag]
        }) %>%
            unlist()
        logOR <- lapply(results, function(x) {
            x$coefficients[tag]
        }) %>%
            unlist()
        lowerCI <- lapply(results, function(x) {
            x$ci.lower[tag]
        }) %>%
            unlist()
        upperCI <- lapply(results, function(x) {
            x$ci.upper[tag]
        }) %>%
            unlist()
        ## extract stats result sheet
        statMat <- data.frame(logOR = logOR, P.Value = P.Value, lowerCI = lowerCI, upperCI = upperCI)
        rownames(statMat) <- rownames(occurrence)
    }
    ## compute adjusted P.Values
    adj.P.Value <- stats::p.adjust(statMat$P.Value, method = "BH")
    ## decide direction
    direction <- ifelse(adj.P.Value >= 0.05, 0, ifelse(statMat$logOR > 0, 1, -1))

    statsDf <- data.frame(statMat, adj.P.Val = adj.P.Value, direction = direction)

    return(list(statsDf = statsDf, design = design, variables = list(interceptVar = interceptVar,
        primaryVar = primaryVar, additionalVar = additionalVar)))
}
