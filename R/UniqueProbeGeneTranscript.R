# TODO: Add comment
#
# Author: Ernur
###############################################################################
#' @include affyCustomCdfClassDef.R
setGeneric("UniqueProbeGeneTranscript", function(.object)
  standardGeneric("UniqueProbeGeneTranscript"))

setMethod("UniqueProbeGeneTranscript", signature("affyCustomCdf"),
          definition=function(.object) {

  data.table::setkeyv(.object@total,cols = "keys")

  .object@modifiedTable = data.table()
  probesIn = split(seq(nrow(.object@total)),.object@total$values)
  indexes =  unlist(lapply(probesIn, `[[`, 1),use.names = FALSE)

  .object@modifiedTable = data.table(keys=.object@total$keys[indexes],
                                     values=.object@total$values[indexes])

  data.table::setkeyv(.object@modifiedTable,cols = "keys")

  return(.object)

}
)
