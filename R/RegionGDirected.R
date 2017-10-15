# Anotates probes based on the suplied gtf file at region (UTR-Exon) level
# Direction of probes and annotations is taken into account.
# Intervals of annotations are combined.
#
# Author: Ernur Saka
###############################################################################
setGeneric("RegionGDirected", function(.object,geneStructure,probes)
  standardGeneric("RegionGDirected"))

setMethod("RegionGDirected", signature("affyCustomCdf"),
          definition=function(.object,geneStructure,probes) {

  chr = split(geneStructure,geneStructure$seqname)
  chrN = length(chr)

  # #For each chromosome
  for(c in 1:chrN ){

    chrName =names(chr)[c]

    if(dim(stats::na.omit(probes[chrName]))[1] != 0){

      Schr = data.frame(chr[c])
      colnames(Schr) = c("seqname", "feature","start","end","strand",
                         "attributes")

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

          regionGenes =  split(directedGenes[[d]],directedGenes[[d]]$attributes)
          start = list()
          end = list()
          name = list()
          number = list()

          for( r in 1:length(regionGenes)){

            nameToKey = paste(regionGenes[[r]]$attributes[1],
                              regionGenes[[r]]$strand[1],sep="_")
            R = IRanges(regionGenes[[r]]$start,regionGenes[[r]]$end)
            R = reduce(R)

            start = c(start, R@start)
            end = c(end, (R@start+ R@width - 1))

            #ranges = append(ranges,R)
            #addition = Rle(nameToKey, length(R))
            name = c(name,nameToKey)
            number = c(number,length(R))
            #allKeys = append(allKeys,addition)

          }

          v1 = rep(as.character(name),as.integer(number))
          allKeys =S4Vectors::Rle(v1)
          ranges = IRanges(as.integer(start),as.integer(end))

          query <- IRanges(directedProbes[[curDirect]]$start,
                           directedProbes[[curDirect]]$start +
                           .object@probeLength-1)
          hits = findOverlaps(query, ranges,type="within",select="all")

          value = directedProbes[[curDirect]]$names[S4Vectors::queryHits(hits)]
          key = allKeys[S4Vectors::subjectHits(hits)]
          all = data.table(keys=as.vector(key),values=value)

          .object@total = rbind(.object@total,all)

        }#end of if(directedProbes[d] != "NULL"

      }#end of directions loop

    }#end of if(dim(na.omit(probes[chrName]))[1]

  }#end of for(c in 1:chrN)

  data.table::setkeyv(.object@total,cols = "keys")
  return(.object)
}
)

