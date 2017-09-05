# TODO: Add comment
#
# Author: Ernur
###############################################################################
#' @include affyCustomCdfClassDef.R
#'
setGeneric("TranscriptDirected", function(.object,geneStructure,probes)
  standardGeneric("TranscriptDirected"))

setMethod("TranscriptDirected", signature("affyCustomCdf"),
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
      colnames(Schr) = c("seqname", "geneId","feature","start","end",
                         "strand","attributes")

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

          nameToKey = paste(directedGenes[[d]]$attributes,
                            directedGenes[[d]]$strand,sep="_")
          geneId = directedGenes[[d]]$geneId[S4Vectors::queryHits(hits)]
          R = IRanges(directedGenes[[d]]$start,directedGenes[[d]]$end)

          query <- IRanges(directedProbes[[curDirect]]$start,
                           directedProbes[[curDirect]]$start +
                             (.object@probeLength)-1)
          hits = findOverlaps(query, R,type="within",select="all")

          value = directedProbes[[curDirect]]$names[S4Vectors::queryHits(hits)]
          key = nameToKey[S4Vectors::subjectHits(hits)]
          geneId = directedGenes[[d]]$geneId[S4Vectors::subjectHits(hits)]
          allKeys = data.table(keys=as.vector(key),values=value,
                               geneId = as.vector(geneId))

          .object@total = rbind(.object@total,allKeys)
          #total = rbind(total,allKeys)
        }

      }

    }#end of if(dim(na.omit(probes[chrName]))[1]

  }#end of for(c in 1:chrN)

  data.table::setkeyv(.object@total,cols = "keys")
  return(.object)

}
)
