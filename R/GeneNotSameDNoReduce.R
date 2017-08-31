# TODO: Add comment
#
# Author: Ernur
###############################################################################

setGeneric("GeneNotSameDNoReduce", function(.object,geneStructure,probes)
  standardGeneric("GeneNotSameDNoReduce"))

setMethod("GeneNotSameDNoReduce", signature("affyCustomCdf"),
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

                #Are there probes on that strand?
                if(dim(stats::na.omit(probes[chrName]))[1] != 0){

                  allKeys = paste(Schr$attributes,Schr$strand,sep="_")
                  ranges = IRanges(as.integer(Schr$start),as.integer(Schr$end))

                  query <- IRanges(probes[chrName]$start,
                                   probes[chrName]$start+.object@probeLength-1)
                  hits = findOverlaps(query, ranges,type="within",select="all")

                  value = probes[chrName]$names[S4Vectors::queryHits(hits)]
                  key = allKeys[S4Vectors::subjectHits(hits)]
                  allKeys = data.table(keys=as.vector(key),values=value)

                  .object@total = rbind(.object@total,allKeys)

                }#end of if(dim(na.omit(probes[chrName]))[1]

              }#end of if(dim(na.omit(probes[chrName]))[1]

            }#end of for(c in 1:chrN)

            data.table::setkeyv(.object@total,cols = "keys")
            return(.object)

          }
)
