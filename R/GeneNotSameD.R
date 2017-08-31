# TODO: Add comment
#
# Author: Ernur
###############################################################################

setGeneric("GeneNotSameD", function(.object,geneStructure,probes)
  standardGeneric("GeneNotSameD"))

setMethod("GeneNotSameD", signature("affyCustomCdf"),
          definition=function(.object,geneStructure,probes) {

  chr = split(geneStructure,geneStructure$seqname)
  chrN = length(chr)

  data.table::setkeyv(probes,cols = "chromosome")
  .object@total = data.table()

  #For each chromosome
  for(c in 1:chrN ){

    chrName =names(chr)[c]

    if(dim(stats::na.omit(probes[chrName]))[1] != 0){

      Schr = data.frame(chr[c])
      colnames(Schr) = c("seqname", "feature","start","end","strand",
                         "attributes")

      regionSplit = split(Schr,Schr$attributes)

      start = list()
      end = list()
      name = list()
      length = list()

      l = length(regionSplit)
      for(i in 1:l){

        R = IRanges(regionSplit[[i]]$start,regionSplit[[i]]$end)
        R = reduce(R)
        start =  c(start,R@start)
        end = c(end, (R@start+ R@width - 1))

        nameToKey = paste(regionSplit[[i]]$attributes[[1]],
                          regionSplit[[i]]$strand[[1]],sep="_")
        name = c(name,nameToKey)
        length = c(length,length(R))
      }

      v1 = rep(as.character(name),as.integer(length))
      allKeysNew = S4Vectors::Rle(v1)
      rangesNew = IRanges(as.integer(start),as.integer(end))

      query <- IRanges(probes[chrName]$start,
                       probes[chrName]$start+.object@probeLength-1)
      hits = findOverlaps(query, rangesNew,type="within",select="all")

      value = probes[chrName]$names[S4Vectors::queryHits(hits)]
      key = allKeysNew[S4Vectors::subjectHits(hits)]
      allKeys = data.table(keys=as.vector(key),values=value)

      .object@total = rbind(.object@total,allKeys)

    }#end of if(dim(na.omit(probes[chrName]))[1]

  }#end of for(c in 1:chrN)

  data.table::setkeyv(.object@total,cols = "keys")
  return(.object)

}
)
