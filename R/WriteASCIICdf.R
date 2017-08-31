# TODO: Add comment
#
# Author: Ernur
###############################################################################
#' @include affyCustomCdfClassDef.R

setGeneric("WriteASCIICdf", function(.object,reportFile)
  standardGeneric("WriteASCIICdf"))

setMethod("WriteASCIICdf", signature("affyCustomCdf"),
          definition=function(.object,reportFile) {

  x = unlist(lapply(strsplit(.object@modifiedTable$values,":"), "[" , 1))
  y = unlist(lapply(strsplit(.object@modifiedTable$values,":"), "[" , 2))


  res = .Call("CreateCustomFile",.object@orginalCdfName,
              .object@newCDFName, reportFile, .object@modifiedTable$keys,
              as.integer(x),as.integer(y),.object@unitType,
              .object@controlProbeSetNumber,.object@probeLength,
              .object@minProbeSetNumber,PACKAGE = "affyCustomCdf")

  # # Sanity check
  # if (is.null(res)) {
  #   stop("Failed to read BPMAP file: ", filename);
  # }

  return(.object)

          }
)
