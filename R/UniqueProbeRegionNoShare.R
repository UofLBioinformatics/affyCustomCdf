# It eliminates the probe share between different annotations at region
# level.
#
# Author: Ernur Saka
###############################################################################

#' @include affyCustomCdfClassDef.R
setGeneric("UniqueProbeRegionNoShare", function(.object,...)
  standardGeneric("UniqueProbeRegionNoShare"))

setMethod("UniqueProbeRegionNoShare", signature("affyCustomCdf"),
          definition=function(.object) {

  data.table::setkeyv(.object@total,cols = "keys")

  fiveUtrTotals = .object@total[grep("FIVE|_5_)",.object@total$keys)]
  threeUtrTotals = .object@total[grep("THREE|_3_)",.object@total$keys)]
  utrTotals = .object@total[grep("UTR",.object@total$keys)]

  exonTotals = .object@total[grep("EXON",.object@total$keys)]

  CDSTotals = .object@total[grep("CDS",.object@total$keys)]

  #Unique to 5-Utr
  fiveLeft = fiveUtrTotals[!(fiveUtrTotals$values  %in% threeUtrTotals$values)]
  fiveLeft = fiveLeft[!(fiveLeft$values  %in% utrTotals$values)]
  fiveLeft = fiveLeft[!(fiveLeft$values  %in% exonTotals$values)]
  fiveLeft = fiveLeft[!(fiveLeft$values  %in% CDSTotals$values)]

  #Unique to 3-Utr
  threeLeft = threeUtrTotals[!(threeUtrTotals$values  %in%
                                fiveUtrTotals$values)]
  threeLeft = threeLeft[!(threeLeft$values  %in% utrTotals$values)]
  threeLeft = threeLeft[!(threeLeft$values  %in% exonTotals$values)]
  threeLeft = threeLeft[!(threeLeft$values  %in% CDSTotals$values)]

  #Unique to UTR
  utrLeft = utrTotals[!(utrTotals$values  %in% fiveUtrTotals$values)]
  utrLeft = utrLeft[!(utrLeft$values  %in% threeUtrTotals$values)]
  utrLeft = utrLeft[!(utrLeft$values  %in% utrTotals$values)]
  utrLeft = utrLeft[!(utrLeft$values  %in% CDSTotals$values)]

  #Unique to EXON
  exonLeft = exonTotals[!(exonTotals$values  %in% fiveUtrTotals$values)]
  exonLeft = exonLeft[!(exonLeft$values  %in% threeUtrTotals$values)]
  exonLeft = exonLeft[!(exonLeft$values  %in% exonTotals$values)]
  exonLeft = exonLeft[!(exonLeft$values  %in% CDSTotals$values)]

  #Unique to CDS
  CDSLeft = CDSTotals[!(CDSTotals$values  %in% fiveUtrTotals$values)]
  CDSLeft = CDSLeft[!(CDSLeft$values  %in% threeUtrTotals$values)]
  CDSLeft = CDSLeft[!(CDSLeft$values  %in% utrTotals$values)]
  CDSLeft = CDSLeft[!(CDSLeft$values  %in% exonTotals$values)]

  .object@modifiedTable = data.table(keys=character(),values=character())

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

  #make probes unique in each region such as one region of two different genes
  #can't share probe.

  data.table::setkeyv(.object@modifiedTable,cols = "keys")

  return(.object)
}
)
