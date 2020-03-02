#!/usr/bin/env Rscript

###if not installed already, install following libraries!
###load required libraries


gnLimit<-opt[2] ###1 or 0 if you want to limit to a particular organism
includePhagosome<-opt[3] ###1 or 0 if you want to include the proteins found in the physical phagosome in the analysis, regardless of whether or not they are found in one group or another.

if(gnLimit=="Y"){
	limitOrg<-readLines("AdvancedMode/orgGroups/limitOrgs.txt")
	limitOrg<-limitOrg[limitOrg != ""]
	limitname<-paste(limitOrg,collapse=".")
	pcaname<-paste("modelOUTPUT/advancedModeOUTPUT/PCAplots/Gr1.vs.Gr2.limit",limitname,"PCAplot.pdf",sep=".")
	GOheat.name<-paste("modelOUTPUT/advancedModeOUTPUT/heatmaps/Gr1.vs.Gr2.limit",limitname,"GO_heatmap.pdf",sep=".")
	geneHeat.name<-paste("modelOUTPUT/advancedModeOUTPUT/heatmaps/Gr1.vs.Gr2.limit",limitname,"gene_heatmap.pdf",sep=".")
	predplot<-paste("modelOUTPUT/advancedModeOUTPUT/Gr1.vs.Gr2.limit",limitname,"predictionPlot.pdf",sep=".")
	predTable<-paste("modelOUTPUT/advancedModeOUTPUT/tables/predictionsDataFrame.limit",limitname,"txt",sep=".")
	genes.name<-paste("modelOUTPUT/advancedModeOUTPUT/geneAndGOsets/Group1-significantGeneSet.limit",limitname,"txt",sep=".")
	GOs.name<-paste("modelOUTPUT/advancedModeOUTPUT/geneAndGOsets/Group1-predictiveGOcategories.limit",limitname,"txt",sep=".")
}else{
	pcaname<-"modelOUTPUT/advancedModeOUTPUT/PCAplots/Gr1.vs.Gr2.PCAplot.pdf"
	GOheat.name<-"modelOUTPUT/advancedModeOUTPUT/heatmaps/Gr1.vs.Gr2.GO_heatmap.pdf"
	geneHeat.name<-"modelOUTPUT/advancedModeOUTPUT/heatmaps/Gr1.vs.Gr2.gene_heatmap.pdf"
	predplot<-"modelOUTPUT/advancedModeOUTPUT/Gr1.vs.Gr2.predictionPlot.pdf"
	predTable<-"modelOUTPUT/advancedModeOUTPUT/tables/predictionsDataFrame.txt"
	genes.name<-"modelOUTPUT/advancedModeOUTPUT/geneAndGOsets/Group1-significantGeneSet.txt"
	GOs.name<-"modelOUTPUT/advancedModeOUTPUT/geneAndGOsets/Group1-predictiveGOcategories.txt"
}

###training species groupings used in paper
###Group1 (eg. phagocytes)
orgsNamesGr1<-readLines("AdvancedMode/orgGroups/Group1.txt")
orgsNamesGr1<-orgsNamesGr1[orgsNamesGr1 != ""]
###Group2 (eg. non-phagocytes)
orgsNamesGr2<-readLines("AdvancedMode/orgGroups/Group2.txt")
orgsNamesGr2<-orgsNamesGr2[orgsNamesGr2 != ""]

######################read in training set presence absence matrix######################################################################
presAbsAll.all<-read.table("dataFiles/presAbsMatrix.txt")
colnames(presAbsAll.all)<-gsub("\\.\\.","\\. ",colnames(presAbsAll.all))
colnames(presAbsAll.all)<-gsub("\\.sp"," sp",colnames(presAbsAll.all))

######################read in annotations######################################################################
annot<-read.table("dataFiles/annotations/phag_nonphag-aVa-mcl.hmms_UniIDs-table.txt")
UIDs.phnph<-strsplit(as.character(annot[,2]),"\\|")
UIDs.phnph<-do.call(rbind,UIDs.phnph)
annot[,2]<-UIDs.phnph[,2]

######################read in test genome significant HMMs######################################################################
# test if there is at least one argument: if not, return an error
if (length(list.files("TestGenomes", pattern="*.txt"))>0) {
print("Reading in test genome lists")
testgenomes<-list.files("TestGenomes", pattern="*.txt", full.names=TRUE)
for (i in 1:length(testgenomes)){
sigmodname<-basename(testgenomes[i])
modname<-gsub("_sigModels.txt","",sigmodname)
#print(modname)
genomehmms<-readLines(testgenomes[i])
genomehmmsAnnot<-annot[match(genomehmms,annot[,1]),2]
genomehmmsMod<-rownames(presAbsAll.all) %in% genomehmmsAnnot
genomehmmsMod<-genomehmmsMod+0
presAbsAll.all[[modname]]<-genomehmmsMod
}

}else{
	print("No test genomes!")
	exit();
}

if(opt[5]=="Y"){
	orgsRemove<-readLines("AdvancedMode/orgGroups/removeOrgs.txt")
	orgsRemove<-orgsRemove[orgsRemove != ""]
	removeOrgs<-vector(length=length(orgsRemove))
	if(length(na.omit(match(orgsRemove,orgsNamesGr1)))>0 || length(na.omit(match(orgsRemove,orgsNamesGr2)))>0 ){
		
		print("Cannot remove organisms that are contained in the training set! Either keep the organism in the data or remove it from the training sets")
		print("training set Group 1")
		print(orgsNamesGr1)
		print("training set Group 2")
		print(orgsNamesGr2)
		print("species you are trying to remove")
		print(orgsRemove)
		stop()
	}
	for(i in 1:length(orgsRemove)){
		if(length(na.omit(match(orgsRemove[i],colnames(presAbsAll.all))))>0){
			removeOrgs[i]<-match(orgsRemove[i],colnames(presAbsAll.all))
		}
	}
	if(sum(removeOrgs+0)>0){presAbsAll.all<-presAbsAll.all[,-removeOrgs]}
}

print("All species in this analysis")
print(colnames(presAbsAll.all))

###############read in genes plus weights for proteins shared among phagocytes and phagosome proteins########################################################
print("Reading in weights matrix")

geneWeights<-read.table("dataFiles/geneWeightsMatrix.txt",row.names=1)

##############find genes enriched in one group compared to another###################################################################################

###if limiting search to proteins present in only one genome (or the union of several genomes), make sure that the limiting genome is not part of the training sets.
if(gnLimit=="Y"){
	
	if(length(na.omit(match(limitOrg,orgsNamesGr1)))>0){
		orgsNamesGr1<-orgsNamesGr1[-match(limitOrg,orgsNamesGr1)]
	}
	if(length(na.omit(match(limitOrg,orgsNamesGr2)))>0){
		orgsNamesGr2<-orgsNamesGr2[-match(limitOrg,orgsNamesGr2)]
	}

}


if( !exists("GOtermList")){
GOtermList <- as.list(GOTERM)
}

testGr1<-match(orgsNamesGr1,colnames(presAbsAll.all))
if(length(na.omit(testGr1)) != length(orgsNamesGr1)){
	stop("Not all of the genomes in Group 1 are represented in the presence/absence matrix. Check the names!")
}
testGr2<-match(orgsNamesGr2,colnames(presAbsAll.all))
if(length(na.omit(testGr2)) != length(orgsNamesGr2)){
	stop("Not all of the genomes in Group 2 are represented in the presence/absence matrix. Check the names!")
}

source("scripts/predictGroups.r")


