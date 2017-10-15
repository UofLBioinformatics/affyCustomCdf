#' Create custom  CDFs for Affymetrix GeneChips based on region, gene and
#' transcript annotations.
#'
#' \tabular{ll}{
#' Package: \tab affyCustomCdf\cr

#' Type: \tab Package\cr
#' Version: \tab 0.1\cr
#' Date: \tab 2017-06-01\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' Long description of the tool
#'
#' \code{\link{CreateAffyCustomCdf}} creates the file based on inputs
#'
#' @name affyCustomCdf_package
#' @aliases affyCustomCdf
#' @docType package
#' @title Create Affymetrix GeneChip custom CDF
#' @author Ernur Saka \email{ernur.saka@louisville.edu}
#' @references
#' \url{http://biorxiv.org/content/early/2017/04/11/126573}
#' @keywords Affymetrix Custom Cdf
#' @examples
#' createAffyCustomCdf("Rat230_2.cdf","rat230Probes.txt",
#' "Rattus_norvegicus.Rnor_6.0.85.gtf")
