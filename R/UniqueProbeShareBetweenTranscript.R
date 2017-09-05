# TODO: Add comment
#
# Author: Ernur
###############################################################################
#' @include affyCustomCdfClassDef.R
#'
setGeneric("UniqueProbeShareBetweenTranscript", function(.object)
  standardGeneric("UniqueProbeShareBetweenTranscript"))

setMethod("UniqueProbeShareBetweenTranscript", signature("affyCustomCdf"),
          definition=function(.object) {

    setkey(.object@total,keys)

    .object@modifiedTable = data.table()

    #Splited based on probes splited index is ordered based on probes
    print(1)
    probes = split(seq(nrow(.object@total)),.object@total$values)
    print(2)
    indexes =  unlist(lapply(probes, `[[`, 1),use.names = FALSE)
    #This part gives the gene_id probe combination. It will be used to add
    #same probe to another transcript of the same gene. Since it will be ordeed
    #based on transcript id.
    #alfabaticly smaller transcript id's gene id will be used to select which
    #trnascripts will be selected
    print(3)
    keys = paste(.object@total$values[indexes],.object@total$geneId[indexes],
                 sep='_')
    print(4)
    temporary = data.frame(tmpkeys = paste(.object@total$values,
                                           .object@total$geneId,sep='_'),
                           keys = .object@total$keys,
                           values=.object@total$values)
    print(5)
    splitTemporary = split(seq(nrow(temporary)),temporary$tmpkeys)
    print(6)
    newIndexes =  unlist(splitTemporary[as.character(keys)],use.names = FALSE)
    print(7)
    .object@modifiedTable = data.table(tmpkeys=temporary$tmpkeys[newIndexes],
                                       keys = temporary$keys[newIndexes],
                                       values=temporary$values[newIndexes])

    print(8)
    setkey(.object@modifiedTable,keys)
    print(9)
    return(.object)

}
)

