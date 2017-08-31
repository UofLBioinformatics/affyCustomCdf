setGeneric("ReadProbes", function(.object) standardGeneric("ReadProbes"))

setMethod("ReadProbes", signature("affyCustomCdf"),
          definition=function(.object) {

  progress = FALSE

  probes = utils::read.csv2(.object@probeAlignmentFile, header = FALSE, sep = "\t", quote = "\"",
                     dec = ",", fill = TRUE, col.names = c("X","Y","sense","chromosome", "start"),
                     comment.char = "",stringsAsFactors = FALSE)

  probeTable = data.table(names=paste(as.numeric(probes$X),
                                      as.numeric(probes$Y),sep=":"),
                          chromosome=unlist(probes$chromosome),
                          start=unlist(probes$start),
                          direction=unlist(probes$sense))

  remove(probes)

  data.table::setkeyv(probeTable,cols = "chromosome")

  return(probeTable)

          }
)
