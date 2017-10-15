# This call the ASCII cdf file writer written in
# C++
# Author: Ernur Saka
###############################################################################
#' @include affyCustomCdfClassDef.R

setGeneric("WriteASCIICdf", function(.object,reportFile,missiggProbesFile)
  standardGeneric("WriteASCIICdf"))

setMethod("WriteASCIICdf", signature("affyCustomCdf"),
          definition=function(.object,reportFile,missiggProbesFile) {

  if(length(.object@modifiedTable) != 0) {
      x = unlist(lapply(strsplit(.object@modifiedTable$values,":"), "[" , 1))
      y = unlist(lapply(strsplit(.object@modifiedTable$values,":"), "[" , 2))
  } else{
    x = unlist(list())
    y = unlist(list())
  }


  res = .Call("CreateCustomFile",.object@orginalCdfName,
              .object@newCDFName, reportFile, .object@modifiedTable$keys,
              as.integer(x),as.integer(y),.object@unitType,
              .object@controlProbeSetNumber,.object@probeLength,
              .object@minProbeSetNumber,missiggProbesFile,
              PACKAGE = "affyCustomCdf")

  # # Sanity check
   # if (is.null(res)) {
   #   #if (!res) {
   #   stop("Failed to write Custom CDF");
   # }

  return(.object)

          }
)
