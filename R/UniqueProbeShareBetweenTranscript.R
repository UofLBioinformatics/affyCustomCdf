# It eliminates the probe share between different annotations at transcript
# level but it lets sharing probes between transcripts.
#
# Author: Ernur Saka
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
    probes = split(seq(nrow(.object@total)),.object@total$values)
    indexes =  unlist(lapply(probes, `[[`, 1),use.names = FALSE)

    #This part gives the gene_id probe combination. It will be used to add
    #same probe to another transcript of the same gene. Since it will be ordeed
    #based on transcript id.
    #alfabaticly smaller transcript id's gene id will be used to select which
    #trnascripts will be selected
    keys = paste(.object@total$values[indexes],.object@total$geneId[indexes],
                 sep='_')
    temporary = data.frame(tmpkeys = paste(.object@total$values,
                                           .object@total$geneId,sep='_'),
                           keys = .object@total$keys,
                           values=.object@total$values,stringsAsFactors=FALSE)

    splitTemporary = split(seq(nrow(temporary)),temporary$tmpkeys)

    newIndexes =  unlist(splitTemporary[as.character(keys)],use.names = FALSE)

    .object@modifiedTable = data.table(tmpkeys=temporary$tmpkeys[newIndexes],
                                       keys = temporary$keys[newIndexes],
                                       values=temporary$values[newIndexes])

    setkey(.object@modifiedTable,keys)

    return(.object)

}
)

