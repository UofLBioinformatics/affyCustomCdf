# TODO: Add comment
#
# Author: Ernur
###############################################################################
#' @include affyCustomCdfClassDef.R
NULL
setGeneric("GeneDirected", function(.object,geneStructure,probes)
  standardGeneric("GeneDirected"))

setMethod("GeneDirected", signature("affyCustomCdf"),
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
      colnames(Schr) = c("seqname", "feature","start","end","strand","attributes")

      directedGenes =  split(Schr,Schr$strand)
      directedProbes = split(probes[chrName],probes[chrName]$direction)
      #for genes on both strands
      for(d in 1:length(directedGenes)){

        if(.object@SD == 1){
          curDirect = directedGenes[[d]]$strand[1]
        } else{
          curDirect = directedGenes[[((d)%%2) +1]]$strand[1]
        }

        #Are there probes on that strand?
        if(directedProbes[curDirect] != "NULL"){

          regionSplit = split(directedGenes[[d]],directedGenes[[d]]$attributes)

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

          query <- IRanges(directedProbes[[curDirect]]$start,
                           directedProbes[[curDirect]]$start+
                             (.object@probeLength)-1)
          hits = findOverlaps(query, rangesNew,type="within",select="all")

          value = directedProbes[[curDirect]]$names[S4Vectors::queryHits(hits)]
          key = allKeysNew[S4Vectors::subjectHits(hits)]
          allKeys = data.table(keys=as.vector(key),values=value)

          .object@total = rbind(.object@total,allKeys)

        }#end of if(dim(na.omit(probes[chrName]))[1]

      }#end of for(d in 1:length(directedGenes))
    }

  }#end of for(c in 1:chrN)

  data.table::setkeyv(.object@total,cols = "keys")
  return(.object)

}
)
