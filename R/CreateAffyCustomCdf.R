#' A function to create Affymetrix GeneChip custom Chip Description File (CDF)
#' createAffyCustomCdf function creates custom CDF for Affymetrix GeneChips based
#' on the supplied annotations of interest in a General/Gene Transfer Format(GTF)
#'
#' createAffyCustomCdf function creates custom CDFs via removing nonspecific
#' probes, updating probe target mapping based on the supplied genome information
#' and grouping probes into region (exon, UTR), gene and transcript levels.
#' User has to supply annotations as a GTF file (Ensemble format), probe start
#' points in a tab separated text file and the original CDF file of a specific
#' Affymetrix GeneChip. It also takes several other parameters includes annotation
#' level (region, gene, transcript), annotation direction, number of probes per
#' probe set. When default parameters are being used, the produced custom CDF
#' will consists of region based probe sets with minimum number of three probes
#' per probe set. Probe sets will not contain original control probe sets and
#' probes aligned to junctions will not be included in the probe sets. The number
#' of probes, will be assumed 25 and alignment will be performed same directed
#' probes and annotations.
#'
#' @param originalCdfName String The original CDF file of the selected Affymetrix
#' GeneChip technology. It can be obtained from Affymetrix NetAffx web site.
#' i.e. Rat230_2.cdf file for Affymetrix GeneChip Rat 230
#' @param probeAlignmentFile String The tab separated probe alignment file.
#' The file has five columns, in order: X Location of a Probe, Y location
#' of a probe, Sense/antisense of a probe (-/+), chromosome name, Start Point.
#' i.e. 335	337	-	8	22161426. Mapping of probes onto a genome cab be performed
#' via Bowtie and tab separated file of exact columns, can be produced via R
#' scripts.
#' @param gtfFileName String The General/Gene Transfer Format(GTF) file name.
#' @param newCDFName String Name of the created custom CDF. Default name is
#' orginalCdfName_min_minProbeSetNumber_D_SD.cdf
#' @param reportFile String Name of the created report file name. Default name
#' is report_orginalCdfNameorginalCdfName.txt
#' @param controlProbeSetNumber Number The number of control probe sets in
#' the original CDF. Control probe set names usually starts with AFFX prescript.
#' If the number of control probe set is given, they will be included in the
#' custom CDF without changes otherwise they will not be included. Default value
#' is 0 (not included).
#' @param minProbeSetNumber number The minimum number of probes in a probe set.
#' Probe sets with less than the minimum number will not be included in the
#' custom CDF. Default value is 3.
#' @param probeLength Number Number of bases in a probe. This information can be
#' obtained via checking a sequence of a probe belongs to an Affymetrix GeneChip.
#' Default value is 25.
#' @param SD Number Sense/antisense relationship between probes and annotations.
#' Possible values are 0 (No direction), 1 (Same direction) and 0 (opposite
#' direction). When SD is 0, direction will not be considered during
#' probes to annotations mapping. When SD is 1, mapping is performed between
#' same directions such as sense probes are being mapped to sense annotations.
#' When SD is 0, mapping is performed between opposite directions such as sense
#' probes are being mapped to anti-sense annotations. Default values is 1 (Same
#' direction).
#' @param type String The type of the created custom CDF. Options are regionG,
#' gene, transcript. When regionG is selected, probe sets are designed to target
#' a specific region (exon, UTR) of a gene. When gene is selected, probe sets
#' are designed to target genes. When transcript is selected, probe sets are
#' designed to target a specific transcript of a gene. Default value is regionG.
#' @param transcriptShare Bool It defines whether to allow probe share
#' between transcripts of a gene or not. It is only used when transcript is
#' selected as the type. Possible values are FALSE for not allowing probe share and
#' TRUE for allowing probe share. Default value is FALSE.
#' @param junction Bool It defines whether we are adding probes to
#' probe sets that map onto a junction of a specific region of a gene or not. It is
#' used when regionG or gene option is selected. Options are FALSE for no
#' junction and TRUE to allow junctions. Default value is FALSE.
#' @param uniqueProbe Bool It defines whether to use probes that map
#' onto more than one annotations or not. It is an option for user to examine unique
#' probes to one specific annotation. It can only be used with regionG option.
#' Possible values are FALSE for not unique, TRUE for unique. Default value is
#' FALSE.
#' @return If all inputs are proper, createAffyCustomCdf returns TRUE and
#' produces two output files. One of the files is the ASCII custom CDF file and
#' the other one is the report file.
#' @examples
#' \dontrun{
#' createAffyCustomCdf("Rat230_2.cdf","rat230HG38.txt","Rattus_norvegicus.Rnor_6.0.85.gtf")
#' }
#' @import stringr
#' @import data.table
#' @import readr
#' @importFrom IRanges IRanges reduce findOverlaps
#' @importFrom methods new
#' @importFrom stats na.omit
#' @importFrom S4Vectors queryHits subjectHits
#' @useDynLib affyCustomCdf, .registration = TRUE
#' @export createAffyCustomCdf
#' @include affyCustomCdfClassDef.R
#' @name CreateAffyCustomCdf
NULL

createAffyCustomCdf = function(orginalCdfName, probeAlignmentFile, gtfFileName,
                               newCDFName, reportFile, controlProbeSetNumber,
                               minProbeSetNumber = 1, probeLength=25,SD=0,
                               type="regionG",transcriptShare = FALSE,
                               junction = FALSE, uniqueProbe=FALSE){

  orginalCdfName  = file.path(dirname(orginalCdfName), basename(orginalCdfName))
  if (!file.exists(orginalCdfName)){
    stop("Original Cdf File not found: ", orginalCdfName)
  }

  probeAlignmentFile  = file.path(dirname(probeAlignmentFile),
                                  basename(probeAlignmentFile))
  if (!file.exists(probeAlignmentFile))
    stop("Probe File not found: ", probeAlignmentFile)

  gtfFileName  = file.path(dirname(gtfFileName), basename(gtfFileName))
  if (!file.exists(gtfFileName))
    stop("Gtf File not found: ", gtfFileName)

  if (missing(controlProbeSetNumber) )
    controlProbeSetNumber = 0

  if (missing(minProbeSetNumber))
    minProbeSetNumber = 3

  unitType = 3

  if (missing(probeLength))
    probeLength = 25

  if (missing(SD))
    SD=0 # 1

  if (missing(type))
    type="regionG"

  #!!!!!!!!!!!!!DOES IT NEED TO BE A CLASS MEMBER FOR NOT IT IS NOT????!!!!!!
  if (missing(junction))
    junction = FALSE

  #!!!!!!!!!!!!!DOES IT NEED TO BE A CLASS MEMBER FOR NOT IT IS NOT????!!!!!!
  if (missing(transcriptShare))
    transcriptShare = FALSE

  #!!!!!!!!!!!!!DOES IT NEED TO BE A CLASS MEMBER FOR NOT IT IS NOT????!!!!!!
  if (missing(probeShare))
    probeShare=1

  if(missing(newCDFName)){

    newCDFName = str_replace(toupper(basename(orginalCdfName)),".CDF","")
    newCDFName = paste(newCDFName,type,"min",minProbeSetNumber,"D",SD,sep = "_")
    newCDFName = paste(newCDFName,".cdf",sep = "")
  }
  if(missing(reportFile)){

    reportFile = str_replace(toupper(basename(orginalCdfName)),".CDF","")
    reportFile = paste("report",reportFile,type,sep = "_")
    reportFile = paste(reportFile,".txt",sep = "")
  }

  #Write inputs into report file
  fileConn<-file(reportFile,"wt")
  cat(file=fileConn, "Original Cdf Name: ", orginalCdfName, "\n")
  cat(file=fileConn, "Probe file name:   ", probeAlignmentFile, "\n")
  cat(file=fileConn, "Probe file name:   ", probeAlignmentFile, "\n")
  cat(file=fileConn, "Gtf File Name: ", gtfFileName, "\n")
  cat(file=fileConn, "Control probe set number:   ", controlProbeSetNumber, "\n")
  cat(file=fileConn, "Unit type:   ", unitType, "\n")
  cat(file=fileConn, "Probe length: ", probeLength, "\n")
  cat(file=fileConn, "Direction between annotatin and probe:   ", SD, "\n")
  cat(file=fileConn, "Type of the produced CDF: ", type, "\n")
  cat(file=fileConn, "Junction:   ", junction, "\n")
  if(type == "trascript")
    cat(file=fileConn, "Transcript share: ", transcriptShare, "\n")
  cat(file=fileConn, "Unique probe:   ", probeShare, "\n")
  close(fileConn)

  newCdf = methods::new("affyCustomCdf",orginalCdfName = orginalCdfName,
               probeAlignmentFile = probeAlignmentFile,
               gtfFileName = gtfFileName, newCDFName = newCDFName,
               controlProbeSetNumber = controlProbeSetNumber,
               probeLength=probeLength, minProbeSetNumber = minProbeSetNumber,
               unitType = unitType, SD = SD, type = type)

  if(is.null(newCdf))
    stop("Could not create the new cdf object ")

  probes = ReadProbes(newCdf)

  gtf = ReadGeneStructure(newCdf)

  if (type == "regionG"){

    if (SD == 0) {
      if(junction == FALSE)
        newCdf = RegionGNotSameDNoReduce(newCdf,gtf,probes)
      else
        newCdf = RegionGNotSameD(newCdf,gtf,probes)
      #RegionGNotSameDNoReduece(gtf,probes,25,"RatGeneND5-02Fast.Rda")
    } else {
      if (junction == FALSE)
        newCdf = RegionGDirectedNoReduce(newCdf,gtf,probes)
      else
        newCdf = RegionGDirected(newCdf,gtf,probes)
    }

    if (uniqueProbe == FALSE)
      newCdf = UniqueProbeRegion(newCdf)
    else
      newCdf = UniqueProbeRegionNoShare(newCdf)

  } else if(type == "gene") {

     if (SD ==  0){
       if (junction == FALSE)
         #!!!!!!NEW CODE WRITTTEN BUT NOT CHECKED IN THE RESULT CDF FILE COMPARE FILES
         newCdf = GeneNotSameDNoReduce(newCdf,gtf,probes)
       else {
         newCdf = GeneNotSameD(newCdf,gtf,probes)
       }
     } else {
       if(junction == FALSE)
         newCdf = GeneDirected(newCdf,gtf,probes)
       else {
         newCdf = GeneDirectedNoReduce(newCdf,gtf,probes)
       }
     }

    newCdf = UniqueProbeGeneTranscript(newCdf)

  } else if(type == "transcript"){

    if (SD == 0)
       newCdf = TranscriptNotSameD(newCdf, gtf,probes)
    else
       newCdf = TranscriptDirected(newCdf,gtf,probes)

    if (transcriptShare == FALSE)
      newCdf = UniqueProbeGeneTranscript(newCdf)
    else
      newCdf = UniqueProbeShareBetweenTranscript(newCdf)
  }

  fileConn<-file(reportFile,open="at")
  cat(file=fileConn, "Number of annotation detected by probes before probe
      sharing eliminated", length(unique(newCdf@total$keys)),"\n")
  cat(file=fileConn, "Number of probe sets without control probes",
      length(unique(newCdf@modifiedTable$keys)), "\n")
  cat(file=fileConn, "Number of probes aligned to an annotation   ",
      length(newCdf@modifiedTable$values), "\n")


  newCdf = WriteASCIICdf(newCdf,reportFile)
  #newCdf = WriteASCIICdf(newCdf)

  return(newCdf)

}
