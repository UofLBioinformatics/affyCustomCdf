# TODO: Add comment
#
# Author: Ernur
###############################################################################

#TRERE IS A TRY LINE MUST BE CHECKES FOR THE SPEED
setGeneric("ReadGeneStructure", function(.object)
  standardGeneric("ReadGeneStructure"))

setMethod("ReadGeneStructure", signature("affyCustomCdf"),
          definition=function(.object) {

            gtf = read_delim(.object@gtfFileName,delim="\t",
                             escape_backslash = FALSE,
                             escape_double = FALSE,
                             col_names = c("seqname", "source", "feature",
                                           "start", "end", "score", "strand",
                                           "frame","attributes"),
                             col_types = "ccciicccc",
                             locale = default_locale(), na = c("", "NA"),
                             comment = "#", skip = 0,
                             n_max = -1, progress = FALSE)

            if(.object@type != "transcript")
              gtf$source=NULL

            gtf$frame=NULL
            gtf$score=NULL

            gtf$feature = toupper(gtf$feature)

            #Changed in order to make gene gtf similar to old one
            #Get rows from gtf based on the gtf and feature .object@type
            #if(.object@type == "gene"){
            #  gtf = subset(gtf,gtf$feature == "gene")
            #}else
            if(.object@type == "gene" | .object@type == "regionG" |
               .object@type == "regionT"){

              gtf = subset(gtf, grepl("*EXON*", gtf$feature)
                           | grepl("*CDS*", gtf$feature)
                           | grepl("*UTR*", gtf$feature))

            }else if(.object@type == "transcript"){ #Left it as else if since we may add the exon based gtf option
              gtf = subset(gtf,grepl("*TRANSCRIPT*", gtf$feature))
            }

            res = lapply(1:length(as.character(gtf$attributes)),function(i)
              strsplit(as.character(gtf$attributes)[i],";"))
            #X = unlist(lapply(strsplit(gtf$attributes,";"), "[[" , 2))
            #res = noquote(res)


            ptm = proc.time()
            id = vector("list", length(res))
            geneId = vector("list", length(res))

            if(.object@type == "regionG" || .object@type == "gene"){

              id = lapply(1:length(res),function(i) res[[i]][[1]][1])
              id = lapply(1:length(id),function(i)
                str_replace(id[[i]],"gene_id ",""))
            }else if(.object@type == "regionT" | .object@type == "transcript"){

              geneId = lapply(1:length(res),function(i) res[[i]][[1]][1])
              geneId = lapply(1:length(geneId),function(i)
                str_replace(geneId[[i]],"gene_id ",""))
              id = lapply(1:length(res),function(i) res[[i]][[1]][3])
              id = lapply(1:length(id),function(i)
                str_replace(res[[i]][[1]][3]," transcript_id ",""))

            }

            ptm = proc.time()
            #Add region to the id
            if(.object@type == "regionG" || .object@type == "regionT"){
              #id = lapply(1:length(id),function(i) paste(id[[i]],gtf$feature[[i]],sep = "_"))
              id = paste(as.character(id),as.character(gtf$feature),sep = "_")
            }


            geneId = lapply(1:length(geneId),function(i)
              str_replace_all(geneId[[i]],"\"",""))

            id = lapply(1:length(id),function(i) str_replace_all(id[[i]],"\"",""))

            if(.object@type == "transcript"){
              gtf$source = as.character(geneId)
              #remove(geneId)
            }

            gtf$attributes = as.character(id)


            #remove(id)

            return(gtf)
          }
)
