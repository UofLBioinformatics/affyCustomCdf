#' @include affyCustomCdfClassDef.R
#'
setGeneric("TranscriptNotSameD", function(.object,geneStructure,probes)
  standardGeneric("TranscriptNotSameD"))

setMethod("TranscriptNotSameD", signature("affyCustomCdf"),
          definition=function(.object,geneStructure,probes) {

  chr = split(geneStructure,geneStructure$seqname)
  chrN = length(chr)

  .object@total = data.table()
  #For each chromosome
  for(c in 1:chrN ){

    chrName =names(chr)[c]

    if(dim(stats::na.omit(probes[chrName]))[1] != 0){

      Schr = data.frame(chr[c])
      colnames(Schr) =  c("seqname", "geneId","feature","start","end",
                          "strand","attributes")

      nameToKey = paste(Schr$attributes,Schr$strand,sep="_")
      geneId = Schr$geneId[S4Vectors::queryHits(hits)]
      R = IRanges(Schr$start,Schr$end)

      query <- IRanges(probes[chrName]$start, probes[chrName]$start+
                         (.object@probeLength)-1)
      hits = findOverlaps(query, R,type="within",select="all")

      value = probes[chrName]$names[S4Vectors::queryHits(hits)]
      key = nameToKey[S4Vectors::subjectHits(hits)]
      geneId = Schr$geneId[S4Vectors::subjectHits(hits)]
      allKeys = data.table(keys=as.vector(key),values=value,
                           geneId = as.vector(geneId))

      .object@total = rbind(.object@total,allKeys)

    }#end of if(dim(na.omit(probes[chrName]))[1]

  }#end of for(c in 1:chrN)

  data.table::setkeyv(.object@total,cols = "keys")
  return(.object)

}
)
