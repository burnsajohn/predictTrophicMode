#!/usr/bin/env Rscript

###if not installed already, install following libraries!
###load required libraries




###training species groupings used in paper
###phagocytes
orgsNamesPh<-c("Aca. castellani","Acy. subglobosum","Bo. saltans","Di. discoideum","Dr. melanogaster","Mu. musculus","Pa. tetraurelia","Te. thermophila","Th. trahens","F. alba","Re. filosa")
###phagomixotrophs
orgsNamesPM<-c("Cym. tetramitiformis","Bi. natans","Chr. tobin")
###plants
orgsNamesPL<-c("Ar. thaliana","O. sativa","Se. moellendorffii","Br. distachyon","Mi. guttatus","Ph. patens","Pi. abies")
###fungi
orgsNamesFN<-c("Ba. dendrobatidis","Sa. cerevisiae","N. californiae","Sc. pombe","Al. macrogynus","Rh. irregularis","Pu. sorghi","Con. coronatus")
###green algae
orgsNamesGA<-c("Chl. reinhardtii","Chlorella sp.","V. carteri")
###red algae
orgsNamesRA<-c("Cya. merolae")


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

testgenomes<-list.files("TestGenomes", pattern="*.txt", full.names=TRUE)
for (i in 1:length(testgenomes)){
sigmodname<-basename(testgenomes[i])
modname<-gsub("_sigModels.txt","",sigmodname)
print(modname)
genomehmms<-readLines(testgenomes[i])
genomehmmsAnnot<-annot[match(genomehmms,annot[,1]),2]
genomehmmsMod<-rownames(presAbsAll.all) %in% genomehmmsAnnot
genomehmmsMod<-genomehmmsMod+0
presAbsAll.all[[modname]]<-genomehmmsMod
}

}else{
	stop("No test genomes!")
}


###############read in genes plus weights for proteins shared among phagocytes and phagosome proteins########################################################
print("Reading in weights matrix")

geneWeights<-read.table("dataFiles/geneWeightsMatrix.txt",row.names=1)

#weightedCategoryScores<-read.table("datasets/weightedCategoryScores.txt",row.names=1)

##############read in list with GO categories and associated genes###################################################################################

PHAG_GOlist<-list.load("dataFiles/PHAG_GOlist.RData")
PHOTO_GOlist<-list.load("dataFiles/PHOTO_GOlist.RData")
PROTOTROPH_GOlist<-list.load("dataFiles/PROTOTROPH_GOlist.RData")
ENT_PHAG_GOlist<-list.load("dataFiles/entamoebids_PHAG_GOlist.RData")
ROZ_PHAG_GOlist<-list.load("dataFiles/r.allomyces_PHAG_GOlist.RData")

#PHAG_GOlist<-list.load("dataFiles/PHAG_GOlist.RData")
#PHAG_GOlist<-list.load("dataFiles/PHAG_GOlist.RData")


preds.list<-list(PHAG_GOlist,PHOTO_GOlist,PROTOTROPH_GOlist,ENT_PHAG_GOlist,ROZ_PHAG_GOlist)


####calculate max category scores
weightedCategoryScores<-list()
for(i in 1:length(preds.list)){
	weightedCategoryScores[[i]]<-data.frame()
	for(j in 1:length(preds.list[[i]])){
		weightedCategoryScores[[i]][names(preds.list[[i]])[j],1]<-sum(1/geneWeights[preds.list[[i]][[j]],1])
	}
}

####################################################
#score GO categories for each organism using the weighted gene scores (multiplicity in GO categories). Initialize outputs
print("Scoring GO categories")
orgsNames2<-c(colnames(presAbsAll.all))

clustscore.ph<-matrix(nrow=length(PHAG_GOlist), ncol=4)
clustscore.photo<-matrix(nrow=length(PHOTO_GOlist), ncol=4)
clustscore.proto<-matrix(nrow=length(PROTOTROPH_GOlist), ncol=4)
clustscore.ent<-matrix(nrow=length(ENT_PHAG_GOlist), ncol=4)
clustscore.roz<-matrix(nrow=length(ROZ_PHAG_GOlist), ncol=4)


clustscore.list<-list(clustscore.ph,clustscore.photo,clustscore.proto,clustscore.ent,clustscore.roz)

clust.filt.ph<-list()
clust.filt.photo<-list()
clust.filt.proto<-list()
clust.filt.ent<-list()
clust.filt.roz<-list()

clust.filt.list<-list(clust.filt.ph,clust.filt.photo,clust.filt.proto,clust.filt.ent,clust.filt.roz)

clustorgscores.ph<-list()
clustorgscores.photo<-list()
clustorgscores.proto<-list()
clustorgscores.ent<-list()
clustorgscores.roz<-list()

clustorgscores.list<-list(clustorgscores.ph,clustorgscores.photo,clustorgscores.proto,clustorgscores.ent,clustorgscores.roz)


###this loop assigns a weighted and normalized category score for each organism in each category.
for (pred in 1:length(preds.list)){
pb <- txtProgressBar(min = 0, max = length(weightedCategoryScores[[pred]][,1]), style = 3)
	for (i in 1:length(preds.list[[pred]])){
		mysums<-list()
		clustgenes<-preds.list[[pred]][[i]]
		clustgenescT<-presAbsAll.all[clustgenes,]
		for (j in 1:length(clustgenescT)){
		mysums[[j]]<-sum(1/(na.omit(geneWeights[match(clustgenes[which(clustgenescT[[j]]>0)],rownames(geneWeights)),1])))/sum(weightedCategoryScores[[pred]][names(preds.list[[pred]])[i],1])#,1/-na.omit(geneWeights[match(clustgenes[which(clustgenescT[[j]]==0)],rownames(geneWeights)),1])
	}
	names(mysums)<-names(clustgenescT)
	clustorgscores.list[[pred]][[i]]<-mysums
	clustscore.list[[pred]][i,1]<-i
	clustscore.list[[pred]][i,2]<-length(clustgenescT)
	setTxtProgressBar(pb, i)
	}
}

##convert to matrix and assign category names
clusts.ph<-matrix(nrow=length(clustorgscores.list[[1]]),ncol=length(orgsNames2))
clusts.photo<-matrix(nrow=length(clustorgscores.list[[2]]),ncol=length(orgsNames2))
clusts.proto<-matrix(nrow=length(clustorgscores.list[[3]]),ncol=length(orgsNames2))
clusts.ent<-matrix(nrow=length(clustorgscores.list[[4]]),ncol=length(orgsNames2))
clusts.roz<-matrix(nrow=length(clustorgscores.list[[5]]),ncol=length(orgsNames2))

clusts.list<-list(clusts.ph,clusts.photo,clusts.proto,clusts.ent,clusts.roz)

hmdf.list<-list()
hmdfc.list<-list()
clustRowsdf.list<-list()

for(pred in 1:length(clusts.list)){
for (i in 1:length(clusts.list[[pred]][,1])){
	clusts.list[[pred]][i,1:length(clusts.list[[pred]][1,])]<-as.numeric(clustorgscores.list[[pred]][[i]][1:length(clustorgscores.list[[pred]][[i]])])
}
colnames(clusts.list[[pred]])<-orgsNames2



clustnames<-names(preds.list[[pred]])

clusts.list[[pred]]<-as.data.frame(clusts.list[[pred]])
rownames(clusts.list[[pred]])<-clustnames
clusts.list[[pred]]<-as.matrix(clusts.list[[pred]])

hmdf.list[[pred]]<-clusts.list[[pred]]
hmdfc.list[[pred]]<-hmdf.list[[pred]]
hmdfc.list[[pred]]<-hmdfc.list[[pred]][!duplicated(hmdfc.list[[pred]]),]
#rownames(hmdfc.list[[pred]])<-make.names(rownames(hmdfc.list[[pred]]), unique=TRUE)
clustRowsdf.list[[pred]]<-hmdfc.list[[pred]]

}

####clustRowsdf is now the working data matrix with weighted summary scores for each category in each organism

###make the heatmap




for(pred in 1:length(clustRowsdf.list)){

clustRowsdfplot<-clustRowsdf.list[[pred]]


PL<-vector(length=length(orgsNamesPL))
for(i in 1:length(orgsNamesPL)){
PL[i]<-grep(orgsNamesPL[i],colnames(clustRowsdfplot))
}

PH<-vector(length=length(orgsNamesPh))
for(i in 1:length(orgsNamesPh)){
PH[i]<-grep(orgsNamesPh[i],colnames(clustRowsdfplot))
}

if(length(orgsNamesPM)>0){
	PM<-vector(length=length(orgsNamesPM))
	for(i in 1:length(orgsNamesPM)){
	PM[i]<-grep(orgsNamesPM[i],colnames(clustRowsdfplot))
	}
}else{PM<-numeric()}


FNall<-vector(length=length(orgsNamesFN))
for(i in 1:length(orgsNamesFN)){
FNall[i]<-grep(orgsNamesFN[i],colnames(clustRowsdfplot))
}

GA<-vector(length=length(orgsNamesGA))
for(i in 1:length(orgsNamesGA)){
GA[i]<-grep(orgsNamesGA[i],colnames(clustRowsdfplot))
}

RA<-vector(length=length(orgsNamesRA))
for(i in 1:length(orgsNamesRA)){
RA[i]<-grep(orgsNamesRA[i],colnames(clustRowsdfplot))
}

d <- dist(t(clustRowsdfplot),method="manhattan") # distance matrix
fitC <- hclust(d, method="ward.D2")

###plot unscaled data
plotHeat<-clustRowsdfplot[order(-(rowSums(clustRowsdfplot[,c(PH,PM)])-rowSums(clustRowsdfplot[,c(PL,FNall,GA,RA)]))),]
plotHeat.genes<-presAbsAll.all[unique(unlist(preds.list[[pred]])),]

d.genes <- dist(t(plotHeat.genes),method="manhattan") # distance matrix
fitC.genes <- hclust(d.genes, method="ward.D2")


if(pred==1){GOheat.name<-"modelOUTPUT/defaultModelOUTPUT/heatmaps/Phagocyte_GOcategory_heatmap.pdf";geneHeat.name<-"modelOUTPUT/defaultModelOUTPUT/heatmaps/Phagocyte_genes_heatmap.pdf"}
if(pred==2){GOheat.name<-"modelOUTPUT/defaultModelOUTPUT/heatmaps/Photosynthesis_GOcategory_heatmap.pdf";geneHeat.name<-"modelOUTPUT/defaultModelOUTPUT/heatmaps/Photosynthesis_genes_heatmap.pdf"}
if(pred==3){GOheat.name<-"modelOUTPUT/defaultModelOUTPUT/heatmaps/Prototrophy_GOcategory_heatmap.pdf";geneHeat.name<-"modelOUTPUT/defaultModelOUTPUT/heatmaps/Prototrophy_genes_heatmap.pdf"}
if(pred==4){GOheat.name<-"modelOUTPUT/defaultModelOUTPUT/heatmaps/Entamoebid.phagocyte_GOcategory_heatmap.pdf";geneHeat.name<-"modelOUTPUT/defaultModelOUTPUT/heatmaps/Entamoebid.phagocyte_genes_heatmap.pdf"}
if(pred==5){GOheat.name<-"modelOUTPUT/defaultModelOUTPUT/heatmaps/R.allomyces.phagocyte_GOcategory_heatmap.pdf";geneHeat.name<-"modelOUTPUT/defaultModelOUTPUT/heatmaps/R.allomyces_genes_heatmap.pdf"}


plotHeat<-as.matrix(plotHeat)
my_palette <- colorRampPalette(brewer.pal(5, "Spectral"))(n = 1000)
pdf(GOheat.name, width = 55, height = 35 )
heatmap.2(plotHeat, Rowv=F,Colv=as.dendrogram(fitC),trace="none",density.info="none",dendrogram="column",col=my_palette,srtCol=45,margins=c(20,60),keysize=0.5,symkey=F,symbreaks=F, cexRow=3, cexCol=3)
dev.off()


plotHeat.genes<-as.matrix(plotHeat.genes)
my_palette <- colorRampPalette(brewer.pal(5, "Spectral"))(n = 1000)
pdf(geneHeat.name, width = 55, height = 35 )
heatmap.2(plotHeat.genes, Rowv=T,Colv=as.dendrogram(fitC.genes),trace="none",density.info="none",dendrogram="column",col=my_palette,srtCol=45,margins=c(20,60),keysize=0.5,symkey=F,symbreaks=F, cexRow=3, cexCol=3)
dev.off()

}

#######################################################################################################################################################################
######################principal components analysis##################################################################################################################
#######################################################################################################################################################################


for(pred in 1:length(clustRowsdf.list)){

clustRowsdfplot<-clustRowsdf.list[[pred]]

PL<-vector(length=length(orgsNamesPL))
for(i in 1:length(orgsNamesPL)){
PL[i]<-grep(orgsNamesPL[i],colnames(clustRowsdfplot))
}

PH<-vector(length=length(orgsNamesPh))
for(i in 1:length(orgsNamesPh)){
PH[i]<-grep(orgsNamesPh[i],colnames(clustRowsdfplot))
}

if(length(orgsNamesPM)>0){
	PM<-vector(length=length(orgsNamesPM))
	for(i in 1:length(orgsNamesPM)){
	PM[i]<-grep(orgsNamesPM[i],colnames(clustRowsdfplot))
	}
}else{PM<-numeric()}


FNall<-vector(length=length(orgsNamesFN))
for(i in 1:length(orgsNamesFN)){
FNall[i]<-grep(orgsNamesFN[i],colnames(clustRowsdfplot))
}

GA<-vector(length=length(orgsNamesGA))
for(i in 1:length(orgsNamesGA)){
GA[i]<-grep(orgsNamesGA[i],colnames(clustRowsdfplot))
}

RA<-vector(length=length(orgsNamesRA))
for(i in 1:length(orgsNamesRA)){
RA[i]<-grep(orgsNamesRA[i],colnames(clustRowsdfplot))
}

if(pred==1 || pred ==4){
learndat<-t(clustRowsdfplot[,c(PH,PM,PL,GA,FNall,RA)])
learndat2<-cbind(c(rep("a",length(c(PH,PM))),rep("b",length(c(PL,GA,FNall,RA)))),learndat)
learndat2<-as.data.frame(learndat2,stringsAsFactors=F)
colnames(learndat2)[1]<-"Category"
for(i in 2:length(learndat2[1,])){learndat2[,i]<-as.numeric(as.character(learndat2[,i]))}
colnames(learndat2)[1]<-"Category"
for(i in 2:length(learndat2[1,])){learndat2[,i]<-as.numeric(as.character(learndat2[,i]))}
}

if(pred==5){

###C.tobin does not cluster with Rozellid phagocytes, stick with training based on cluster from heatmap.
orgsNamesPhRoz<-c("Aca. castellani","Acy. subglobosum","Bo. saltans","Di. discoideum","Dr. melanogaster","Mu. musculus","Th. trahens","F. alba","Re. filosa","Bi. natans")
PHroz<-vector(length=length(orgsNamesPhRoz))
for(i in 1:length(orgsNamesPhRoz)){
PHroz[i]<-grep(orgsNamesPhRoz[i],colnames(clustRowsdfplot))
}
learndat<-t(clustRowsdfplot[,c(PHroz,PL,GA,FNall,RA)])
learndat2<-cbind(c(rep("a",length(c(PHroz))),rep("b",length(c(PL,GA,FNall,RA)))),learndat)
learndat2<-as.data.frame(learndat2,stringsAsFactors=F)
colnames(learndat2)[1]<-"Category"
for(i in 2:length(learndat2[1,])){learndat2[,i]<-as.numeric(as.character(learndat2[,i]))}
colnames(learndat2)[1]<-"Category"
for(i in 2:length(learndat2[1,])){learndat2[,i]<-as.numeric(as.character(learndat2[,i]))}
}

if(pred==2){
learndat<-t(clustRowsdfplot[,c(PH,FNall,PM,PL,GA,RA)])
learndat2<-cbind(c(rep("a",length(c(PH,FNall))),rep("b",length(c(PM,PL,GA,RA)))),learndat)
learndat2<-as.data.frame(learndat2,stringsAsFactors=F)
colnames(learndat2)[1]<-"Category"
for(i in 2:length(learndat2[1,])){learndat2[,i]<-as.numeric(as.character(learndat2[,i]))}
colnames(learndat2)[1]<-"Category"
for(i in 2:length(learndat2[1,])){learndat2[,i]<-as.numeric(as.character(learndat2[,i]))}
}

if(pred==3){
learndat<-t(clustRowsdfplot[,c(PH,FNall,PM,PL,GA,RA)])
learndat2<-cbind(c(rep("a",length(c(PH))),rep("b",length(c(FNall,PM,PL,GA,RA)))),learndat)
learndat2<-as.data.frame(learndat2,stringsAsFactors=F)
colnames(learndat2)[1]<-"Category"
for(i in 2:length(learndat2[1,])){learndat2[,i]<-as.numeric(as.character(learndat2[,i]))}
colnames(learndat2)[1]<-"Category"
for(i in 2:length(learndat2[1,])){learndat2[,i]<-as.numeric(as.character(learndat2[,i]))}
}

###principal compenents on training data
source("scripts/pcaPlots.r")

###plot full predictions in mollweide projection plot
source("scripts/plotMollweide.r")

}
