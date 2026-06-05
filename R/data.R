#' Myeloma gene-expression dataset (CTA benchmark)
#'
#' A data frame with 256 observations and 19 variables, formatted for use
#' with \code{\link{cta_fit}} and \code{\link{oda_fit}}.  Derived from the
#' publicly available myeloma gene-expression dataset (GEO accession GSE4581),
#' as distributed in the \pkg{survminer} package.
#'
#' This dataset is used throughout the oda documentation and vignettes to
#' illustrate weighted CTA, MINDENOM constraints, LOO STABLE validation, and
#' missing-code handling.  Reference CTA.exe golden outputs for MINDENOM = 1,
#' 30, and 56 are used as regression anchors.
#'
#' @format A data frame with 256 rows and 19 columns:
#' \describe{
#'   \item{V1}{Survival event indicator (0 = censored, 1 = event).
#'     Used as the class variable \code{y} in CTA/ODA.}
#'   \item{V2}{Case weight (observation time in months).
#'     Use as \code{w} in \code{cta_fit}; rows with V2 == 0 should be excluded.}
#'   \item{V3}{CCND1 gene expression.}
#'   \item{V4}{CRIM1 gene expression.}
#'   \item{V5}{DEPDC1 gene expression.}
#'   \item{V6}{IRF4 gene expression.}
#'   \item{V7}{TP53 expression / mutation burden.}
#'   \item{V8}{WHSC1 gene expression.}
#'   \item{V9}{Molecular group: Cyclin D-1 (binary).}
#'   \item{V10}{Molecular group: Cyclin D-2 (binary).}
#'   \item{V11}{Molecular group: Hyperdiploid (binary).}
#'   \item{V12}{Molecular group: Low bone disease (binary).}
#'   \item{V13}{Molecular group: MAF (binary).}
#'   \item{V14}{Molecular group: MMSET (binary).}
#'   \item{V15}{Molecular group: Proliferation (binary).}
#'   \item{V16}{Chr1q21 status: 2 copies (binary).}
#'   \item{V17}{Chr1q21 status: 3 copies (binary).}
#'   \item{V18}{Chr1q21 status: 4+ copies (binary).}
#'   \item{V19}{Chr1q21 status: NA-coded (binary).
#'     Missing values are coded as -9 (\code{miss_codes = -9}).}
#' }
#'
#' @details
#' Use \code{miss_codes = -9} and \code{w = myeloma$V2} when calling
#' \code{cta_fit}.  With \code{mindenom = 1}, the enumerated CTA tree roots
#' at V14 with a V15 child (OVERALL ESS = 26.32\%, WEIGHTED ESS = 27.69\%).
#' With \code{mindenom = 30}, the selected tree is a V17 stump
#' (WEIGHTED ESS = 16.51\%).  With \code{mindenom = 56}, no admissible
#' tree exists.
#'
#' @source
#' Derived from the \code{myeloma} dataset in the \pkg{survminer} package.
#' Original data: NCBI GEO accession GSE4581.  No PHI; no institutional data.
#' See \code{tests/testthat/fixtures/myeloma/README.md} in the source tree.
#'
#' @name myeloma
#' @docType data
#' @keywords datasets
NULL

#' CTA demonstration dataset
#'
#' A simulated data frame with 200 observations and 6 variables, designed to
#' illustrate Classification Tree Analysis with \code{\link{cta_fit}}.
#' This is the dataset used in the CTA.exe demonstration program.
#'
#' The CTA.exe golden output for MINDENOM = 1 selects V2 as root
#' (cut = 4.5, ESS = 52.63\%).  MINDENOM = 8 requires \code{mc_iter = 25000}
#' for parity.
#'
#' @format A data frame with 200 rows and 6 columns:
#' \describe{
#'   \item{V1}{Class label (integer; 1 or 2).}
#'   \item{V2}{Ordered attribute (root in MINDENOM = 1 solution).}
#'   \item{V3}{Ordered attribute.}
#'   \item{V4}{Binary attribute (0/1).}
#'   \item{V5}{Ordered attribute.}
#'   \item{V6}{Ordered attribute.}
#' }
#'
#' @details
#' Simulated dataset; no real subjects or PHI.  Used as the primary
#' introductory CTA example in the oda package vignettes and in the
#' CTA.exe demonstration program (\code{CTA_DEMO.pgm}).
#'
#' @name cta_demo
#' @docType data
#' @keywords datasets
NULL
