# It eliminates the probe share between different annotations at region
#level.
#
# Author: Ernur Saka
###############################################################################

setGeneric("UniqueProbeRegion", function(.object)
  standardGeneric("UniqueProbeRegion"))

setMethod("UniqueProbeRegion", signature("affyCustomCdf"),
          definition=function(.object) {

  fiveUtrTotals = .object@total[grep("FIVE|_5_)",.object@total$keys)]
  threeUtrTotals = .object@total[grep("THREE|_3_)",.object@total$keys)]
  utrTotals = .object@total[grep("UTR",.object@total$keys)]

  exonTotals = .object@total[grep("EXON",.object@total$keys)]

  CDSTotals = .object@total[grep("CDS",.object@total$keys)]

  #Remove shared probes between regions
  threeLeft = threeUtrTotals[!(threeUtrTotals$values  %in%
                                fiveUtrTotals$values)]

  utrLeft = utrTotals[!(utrTotals$values  %in% fiveUtrTotals$values)]
  utrLeft = utrLeft[!(utrLeft$values  %in% threeLeft$values)]

  exonLeft = exonTotals[!(exonTotals$values  %in% fiveUtrTotals$values)]
  exonLeft = exonLeft[!(exonLeft$values  %in% threeLeft$values)]
  exonLeft = exonLeft[!(exonLeft$values  %in% utrLeft$values)]

  CDSLeft = CDSTotals[!(CDSTotals$values  %in% fiveUtrTotals$values)]
  CDSLeft = CDSLeft[!(CDSLeft$values  %in% threeLeft$values)]
  CDSLeft = CDSLeft[!(CDSLeft$values  %in% utrLeft$values)]
  CDSLeft = CDSLeft[!(CDSLeft$values  %in% exonLeft$values)]


  #modifiedTable = data.table(keys=character(),values=character())

  #split is used to not share probes between genes of same region
  if(dim(fiveUtrTotals)[1] > 0){
    probes = split(seq(nrow(fiveUtrTotals)),fiveUtrTotals$values)
    firstIndexes = unlist(lapply(probes, `[[`, 1),use.names = FALSE)
    .object@modifiedTable = rbindlist(list(.object@modifiedTable,
                                  as.list(data.table(
                                  keys=fiveUtrTotals$keys[firstIndexes],
                                  values=fiveUtrTotals$values[firstIndexes]))))
  }

  if(dim(threeLeft)[1] > 0){
    probes = split(seq(nrow(threeLeft)),threeLeft$values)
    firstIndexes = unlist(lapply(probes, `[[`, 1),use.names = FALSE)
    .object@modifiedTable = rbindlist(list(.object@modifiedTable,
                                  as.list(data.table(
                                  keys=threeLeft$keys[firstIndexes],
                                  values=threeLeft$values[firstIndexes]))))
  }

  if(dim(utrLeft)[1] > 0){
    probes = split(seq(nrow(utrLeft)),utrLeft$values)
    firstIndexes = unlist(lapply(probes, `[[`, 1),use.names = FALSE)
    .object@modifiedTable = rbindlist(list(.object@modifiedTable,
                                  as.list(data.table(
                                  keys=utrLeft$keys[firstIndexes],
                                  values=utrLeft$values[firstIndexes]))))
  }

  if(dim(exonLeft)[1] > 0){
    probes = split(seq(nrow(exonLeft)),exonLeft$values)
    firstIndexes = unlist(lapply(probes, `[[`, 1),use.names = FALSE)
    .object@modifiedTable = rbindlist(list(.object@modifiedTable,
                                  as.list(data.table(
                                  keys=exonLeft$keys[firstIndexes],
                                  values=exonLeft$values[firstIndexes]))))
  }

  if(dim(CDSLeft)[1] > 0){
    probes = split(seq(nrow(CDSLeft)),CDSLeft$values)
    firstIndexes = unlist(lapply(probes, `[[`, 1),use.names = FALSE)
    .object@modifiedTable = rbindlist(list(.object@modifiedTable,
                                  as.list(data.table(
                                  keys=CDSLeft$keys[firstIndexes],
                                  values=CDSLeft$values[firstIndexes]))))
  }

  data.table::setkeyv(.object@modifiedTable,cols = "keys")
  return(.object)
  }
)

