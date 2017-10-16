# TODO: It demonstrates how to obtain the proper column arrangement for probe
# mapping text file. The original files are produced via Bowtie. For more
# information please check introduction to affyCustomCdf file.
#
# Author: Ernur Saka
###############################################################################

ProbeFileClear = function(){

	require(stringr)
	require (data.table)

	#For HUMAN File Remove probe names leave only X and Y
	fileName = "human.txt"
	probes = read.csv2(fileName, header = FALSE, sep = "\t", quote = "\"",
			dec = ",", fill = TRUE, col.names = c("names","sense","chromosome",
			                                      "start","empty"),comment.char = "",
			stringsAsFactors = FALSE)

	probeNames = unlist(lapply(strsplit(probes[,1],";"), "[[" , 1))
	X = unlist(lapply(strsplit(probeNames,":"), "[[" , 4))
	Y = unlist(lapply(strsplit(probeNames,":"), "[[" , 5))

	probeTable = data.table(X = as.numeric(X), Y = as.numeric(Y),
	                        direction=unlist(probes$sense),
	                        chromosome=unlist(probes$chromosome),
	                        start=unlist(probes$start))

	write.table(probeTable, "HGU133Plus2-human38DNAAligned.txt", quote = FALSE, row.names = FALSE,
	            col.names = FALSE, sep="\t")
	remove(probeTable,X,Y,probes)

	#For RAT File Remove probe names leave only X and Y
	fileName = "rat.txt"
	probes = read.csv2(fileName, header = FALSE, sep = "\t", quote = "\"",
			dec = ",", fill = TRUE, col.names = c("names","sense","chromosome",
			                                      "start","empty"),comment.char = "",
			stringsAsFactors = FALSE)

	X = unlist(lapply(strsplit(probes[,1],";"), "[[" , 2))
	Y = unlist(lapply(strsplit(probes[,1],";"), "[[" , 3))

	probeTable = data.table(X = as.numeric(X), Y = as.numeric(Y),
	                        direction=unlist(probes$sense),
	                        chromosome=unlist(probes$chromosome),
	                        start=unlist(probes$start))

	write.table(probeTable, "Rat230-2Rn6ChrAligned.txt", quote = FALSE, row.names = FALSE,
	            col.names = FALSE, sep="\t")
	remove(probeTable,X,Y,probes)

	#For MOUSE GENE File Remove probe names leave only X and Y
	fileName = "mouseGene1.txt"
	probes = read.csv2(fileName, header = FALSE, sep = "\t", quote = "\"",
			dec = ",", fill = TRUE, col.names = c("names","sense","chromosome",
			                                      "start","empty"),comment.char = "",
			stringsAsFactors = FALSE)

	probeNames = unlist(lapply(strsplit(probes[,1],";"), "[[" , 2))
	X = unlist(lapply(strsplit(probeNames,":"), "[[" , 1))
	Y = unlist(lapply(strsplit(probeNames,":"), "[[" , 2))

	probeTable = data.table(X = as.numeric(X), Y = as.numeric(Y),
	                        direction=unlist(probes$sense),
	                        chromosome=unlist(probes$chromosome),
	                        start=unlist(probes$start))

	write.table(probeTable, "mouseGene1Aligned.txt", quote = FALSE, row.names = FALSE,
	            col.names = FALSE, sep="\t")
	remove(probeTable,X,Y,probes)

	#For MOUSE 430 File Remove probe names leave only X and Y
	fileName = "mouse430.txt"
	probes = read.csv2(fileName, header = FALSE, sep = "\t", quote = "\"",
			dec = ",", fill = TRUE, col.names = c("names","sense","chromosome",
			                                      "start","empty"),comment.char = "",
			stringsAsFactors = FALSE)

	probeNames = unlist(lapply(strsplit(probes[,1],";"), "[[" , 1))
	X = unlist(lapply(strsplit(probes[,1],":"), "[[" , 4))
	Y = unlist(lapply(strsplit(probeNames,":"), "[[" , 5))

	probeTable = data.table(X = as.numeric(X), Y = as.numeric(Y),
	                        direction=unlist(probes$sense),
	                        chromosome=unlist(probes$chromosome),
	                        start=unlist(probes$start))

	write.table(probeTable, "Mouse430MouseGRCm3ChrAligned.txt", quote = FALSE, row.names = FALSE,
	            col.names = FALSE, sep="\t")
	remove(probeTable,X,Y,probes)

}
