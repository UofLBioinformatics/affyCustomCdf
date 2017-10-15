#' An S4 class to represent the Affymetrix Probe alignmnets to an annotation
#'  An S4 class is used to create a custom CDF sturucture
#'
#' 'Long descrition of the function will be here!!!!!!!!!!!!!!!!!
#' @slot total is a data table which keeps the aligment results of probes to the
#'  annotations
#' @slot modifiedTable is a data tabe which  contains the remainder of total
#' after probe shaing procedues applied.
#' @slot orginalCdfName the name of the original cdf file.
#' @slot probeAlignmentFile name of the file that contains the genome aligment
#' results of probes.
#' @slot gtfFileName gtf file name that is being used as a annotation source.
#' @slot newCDFName The name of the created custom CDF
#' @slot controlProbeSetNumber The number of control probe set located in the
#' custom CDF
#' @slot probeLength The number of probe bases.
#' @slot minProbeSetNumber Desired minimum number of probes in a probe set.
#' @slot unitType The type of the probe set units
#' @slot SD Direction relation between annotation and probe. Values 0 no
#' directional relation checked
#' @slot type type of the created custom CDF. It can be regionG, gene or
#' transcript.
#' @slot MM It show wheather the original CDF has mismatched probes or not.
#' @name affyCustomCdfClass
#' @author Ernur Saka
setOldClass(c('data.frame'))
setOldClass(c('data.table', 'data.frame'))

setClass("affyCustomCdf",
                slots = c(total = 'data.table', modifiedTable = 'data.table',
                orginalCdfName = "character", probeAlignmentFile = "character",
                gtfFileName = "character", newCDFName = "character",
                controlProbeSetNumber = "numeric", probeLength = "numeric",
                minProbeSetNumber = "numeric", unitType = "numeric",
                SD = "numeric",type = "character")
)

setMethod("initialize","affyCustomCdf",
            function(.Object, orginalCdfName = "x",probeAlignmentFile = "x",
                    gtfFileName = "x", newCDFName="x",
                    controlProbeSetNumber = as.integer(0),
                    probeLength = as.integer(25),
                    minProbeSetNumber = as.integer(1) ,
                    unitType = as.integer(3), SD = as.integer(0),
                    type = "regionG"){

            .Object@orginalCdfName = orginalCdfName
            .Object@probeAlignmentFile = probeAlignmentFile
            .Object@gtfFileName = gtfFileName
            .Object@newCDFName = newCDFName
            .Object@controlProbeSetNumber = controlProbeSetNumber
            .Object@probeLength = probeLength
            .Object@minProbeSetNumber = minProbeSetNumber
            .Object@unitType = unitType
            .Object@SD = SD
            .Object@type = type
            .Object
        })

    .affyCustomCdf.valid <- function(object){
        if(object@probeLength < 0){
            return("length mismatch")
        }
        if(object@minProbeSetNumber < 0){
            return("length mismatch")
        }
        return(TRUE)
     }

setValidity("affyCustomCdf", .affyCustomCdf.valid)
