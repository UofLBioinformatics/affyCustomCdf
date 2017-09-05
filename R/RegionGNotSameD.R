setGeneric("RegionGNotSameD", function(.object,geneStructure,probes)
  standardGeneric("RegionGNotSameD"))

setMethod("RegionGNotSameD", signature("affyCustomCdf"),
          definition=function(.object,geneStructure,probes) {

  chr = split(geneStructure,geneStructure$seqname)
  chrN = length(chr)

  data.table::setkeyv(probes,cols = "chromosome")

  #For each chromosome
  for(c in 1:chrN ){

    chrName = names(chr)[c]

    if(dim(stats::na.omit(probes[chrName]))[1] != 0){

      Schr = data.frame(chr[c])
      colnames(Schr) = c("seqname", "feature","start","end","strand","attributes")

      allKeys = S4Vectors::Rle()
      ranges = IRanges()
      st = list()
      en = list()
      name = list()
      number = list()

      regionGenes =  split(Schr,Schr$attributes)

      l = length(regionGenes)

      #ptm = proc.time()
      for( r in 1:l){

        nameToKey = paste(regionGenes[[r]]$attributes[1],
                          regionGenes[[r]]$strand[1],sep="_")
        #R = IRanges(regionGenes[[r]]$start,regionGenes[[r]]$end,names=rep(nameToKey,length(regionGenes[[r]]$start)))
        R = IRanges(regionGenes[[r]]$start,regionGenes[[r]]$end)
        R = reduce(R)

        st =  c(st,R@start)
        en = c(en, (R@start+ R@width - 1))

        #ranges = append(ranges,R)
        #addition = Rle(nameToKey, length(R))
        name = c(name,nameToKey)
        number = c(number,length(R))
        #allKeys = append(allKeys,addition)

      }

      v1 = rep(as.character(name),as.integer(number))
      allKeys = S4Vectors::Rle(v1)
      ranges = IRanges(as.integer(st),as.integer(en))

      #proc.time() - ptm

      #	  ptm = proc.time()
      #	  for( r in 1:l){
      #
      #		  nameToKey = paste(regionGenes[[r]]$attributes[1],regionGenes[[r]]$strand[1],sep="_")
      #		  R = IRanges(regionGenes[[r]]$start,regionGenes[[r]]$end)
      #		  R = reduce(R)
      #		  ranges = append(ranges,R)
      #		  addition = Rle(nameToKey, length(R))
      #		  allKeys = append(allKeys,addition)
      #
      #	  }
      #	  proc.time() - ptm

      query <- IRanges(probes[chrName]$start, probes[chrName]$start +
                         .object@probeLength-1)
      hits = findOverlaps(query, ranges,type="within",select="all")

      value = probes[chrName]$names[S4Vectors::queryHits(hits)]
      key = allKeys[S4Vectors::subjectHits(hits)]
      all = data.table(keys=as.vector(key),values=value)

      .object@total = rbind(.object@total,all)

    }
    #}end of if(dim(na.omit(probes[chrName]))[1]
  }

  data.table::setkeyv(.object@total,cols = "keys")
  return(.object)

  }
)
