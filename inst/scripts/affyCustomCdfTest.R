affyCustomCdfTest = function(){

  require("affyCustomCdf")

  pathLib = system.file(package="affyCustomCdf")
  pathLib = file.path(pathLib, "scripts")
  cdf = file.path(pathLib, "Rat230_2Sort.cdf")
  probeMap = file.path(pathLib, "rat230HG38Short2.txt")
  gtf = file.path(pathLib, "Rattus_norvegicus.Rnor_6.0.85Sort.gtf")
  createAffyCustomCdf(cdf,probeMap,gtf,minProbeSetNumber=1)
  createAffyCustomCdf(cdf,probeMap,gtf)

  #Files for the following examples are supplied in the affyCustomCdfFull github
  #branch.
  #https://github.com/UofLBioinformatics/affyCustomCDF/tree/affyCustomCdfFull
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=1, SD=0)
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=3, SD=0)
  createAffyCustomCdf("HG-U133_Plus_2.cdf","HGU133Plus2-human38DNAAligned.txt",
                      "Homo_sapiens.GRCh38.85.gtf",
                      controlProbeSetNumber = 62,
                      minProbeSetNumber=1, SD=0)
  createAffyCustomCdf("HG-U133_Plus_2.cdf","HGU133Plus2-human38DNAAligned.txt",
                      "Homo_sapiens.GRCh38.85.gtf",
                      controlProbeSetNumber = 62,
                      minProbeSetNumber=1, SD=1)
  createAffyCustomCdf("Mouse430_2.cdf","Mouse430MouseGRCm3ChrAligned.txt",
                      "Mus_musculus.GRCm38.85.gtf",
                      controlProbeSetNumber = 64,
                      minProbeSetNumber=1, SD=0)
  createAffyCustomCdf("MoGene-1_0-st-v1.r3.cdf","mouseGene1Aligned.txt",
                      "Mus_musculus.GRCm38.85.gtf",
                      controlProbeSetNumber = 6568,
                      SD=0)

  #Creation of different types of custom CDFs for Rat 230
  #TYPE regionG
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=1, SD=0,junction=FALSE,
                      uniqueProbe=FALSE,type="regionG")

  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=1, SD=0,junction=FALSE,
                      uniqueProbe=TRUE,type="regionG")

  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=1, SD=0,junction=TRUE,
                      uniqueProbe=FALSE, type="regionG")
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=1, SD=0,junction=TRUE,
                      uniqueProbe=TRUE, type="regionG")

  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=1, SD=1,junction=FALSE,
                      uniqueProbe=FALSE, type="regionG")
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=1, SD=1,junction=FALSE,
                      uniqueProbe=TRUE, type="regionG")

  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=1, SD=1,junction=TRUE,
                      uniqueProbe=FALSE, type="regionG")
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=1, SD=1,junction=TRUE,
                      uniqueProbe=TRUE, type="regionG")

  #TYPE gene
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=3, SD=0, junction=TRUE, type="gene")
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=3, SD=0, junction=FALSE, type="gene")
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=3, SD=1, junction=TRUE, type="gene")
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=3, SD=1, junction=FALSE, type="gene")
  #TYPE transcript
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=3, SD=0,transcriptShare= FALSE,
                      type="transcript")
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=3, SD=0,transcriptShare= TRUE,
                      type="transcript")
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=3, SD=1,transcriptShare= FALSE,
                      type="transcript")
  createAffyCustomCdf("Rat230_2.cdf","Rat230-2Rn6ChrAligned.txt",
                      "Rattus_norvegicus.Rnor_6.0.85.gtf",
                      controlProbeSetNumber = 57,
                      minProbeSetNumber=3, SD=1,transcriptShare= TRUE,
                      type="transcript")


  #setwd("C:\\Users\\Ernur\\Documents\\MicroArrayProject\\Packages\\affyCustomCdf\\inst\\scripts")

}
