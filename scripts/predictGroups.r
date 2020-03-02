

###limit training set to only genes present in one or several genomes
if(gnLimit=="Y"){
	limit.gn<-match(limitOrg,colnames(presAbsAll.all))
	if(length(na.omit(limit.gn)) != length(limitOrg)){
		print("Error: A limit genome is missing from the data! (Or maybe the names don't match perfectly (and case sensitive)... no fuzzy matching here). Here are all species in the presence/absense matrix:")
		print(colnames(presAbsAll.all))
		print("Here are the limit species:")
		print(limitOrg)
		stop()
	}
	if(length(limit.gn)>1){
		presAbsAll.limit<-presAbsAll.all[rowSums(presAbsAll.all[,limit.gn])>0,]
	}else{presAbsAll.limit<-presAbsAll.all[presAbsAll.all[,limit.gn]>0,]}
	clustRowsdf.genes<-presAbsAll.limit
}else{clustRowsdf.genes<-presAbsAll.all}

orgsNames<-c(colnames(presAbsAll.all))

#########Find genes with significantly different pattern of presence/absence between groups#########################################
Gr1.gn<-vector(length=length(orgsNamesGr1))
for(i in 1:length(orgsNamesGr1)){
Gr1.gn[i]<-grep(orgsNamesGr1[i],colnames(clustRowsdf.genes))
}

Gr2.gn<-vector(length=length(orgsNamesGr2))
for(i in 1:length(orgsNamesGr2)){
Gr2.gn[i]<-grep(orgsNamesGr2[i],colnames(clustRowsdf.genes))
}

###function for running proportion test
getPredictors.propTest<-function(gr1,gr2,datamat){
	gr1.tot<-length(unlist(gr1))
	gr2.tot<-length(unlist(gr2))
	ptst.out<-apply(datamat,1,function(propt){
			gr1.obs<-sum(propt[unlist(gr1)])
			gr2.obs<-sum(propt[unlist(gr2)])
			sigs<-prop.test(c(gr1.obs,gr2.obs),c(gr1.tot,gr2.tot),alternative="g")
		}
	)
	preddf<-data.frame(nrow=length(ptst.out),ncol=2)
	for(i in 1:length(ptst.out)){
		preddf[i,1]<-names(ptst.out)[i]
		preddf[i,2]<-ptst.out[[i]]$p.value
	}
	preddf
}

gr1<-list(Gr1.gn)
gr2<-list(Gr2.gn)	
sigGenes<-getPredictors.propTest(gr1,gr2,clustRowsdf.genes)
sigGenes<-sigGenes[order(sigGenes[,2]),]
sigGenes[1:10,]

alpha.gn<-genepval
sigGenes<-na.omit(sigGenes[sigGenes[,2]<alpha.gn,1])

sigGenes<-unique(sigGenes)

clustRowsdfplot.genes<-clustRowsdf.genes[match(sigGenes,rownames(clustRowsdf.genes)),]
dim(clustRowsdfplot.genes)



###can include proteins found physically associated with the phagosome. They do not tend to be differentially present between groups of eukaryotes.
if(includePhagosome=="Y"){
	phagosomeHMMs<-readLines("dataFiles/annotations/MousePhagosome_sigModels.txt")
	phagAnnot<-annot[annot[,1] %in% phagosomeHMMs,2]
	phagpresabs<-presAbsAll.all[unique(phagAnnot),]
	if(gnLimit=="Y"){
		sigGenes<-c(sigGenes,rownames(phagpresabs[phagpresabs[,limit.gn]==1,]))
		clustRowsdfplot.genes<-clustRowsdf.genes[match(sigGenes,rownames(clustRowsdf.genes)),]
	}else{
		sigGenes<-c(sigGenes,rownames(phagpresabs))
		clustRowsdfplot.genes<-clustRowsdf.genes[match(sigGenes,rownames(clustRowsdf.genes)),]
	}
}


####Group genes differently represented between groups by GO category
geneID2GO <- readMappings(file = "dataFiles/annotations/phag-nonphag-hmmAnnot.goa",IDsep=";|,")
str(head(geneID2GO))

myInterestingGenes <- sigGenes

geneNames<-names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))

names(geneList) <- geneNames

GOdataBP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisherBP.weight <- runTest(GOdataBP, algorithm = "weight", statistic = "fisher")
allResBP<- GenTable(GOdataBP, weight = resultFisherBP.weight, orderBy = "weight", topNodes = length(attributes(attributes(GOdataBP)$graph)$nodes)) 


####function to pull proteins associated with each GO term after GO enrichment analysis.
getGOgenes<-function(allResTable,updown,GOdat){
	GOtermGenes<-list()
	for (i in 1:length(allResTable[,1])){
		#print(i)
		cmd<-paste("ls(attributes(attributes(attributes(GOdat)$graph)$nodeData)$data$`",allResTable[i,1],"`$genes)",sep="")
		GOgenes<-eval(parse(text=cmd))
		updowngenes<-intersect(GOgenes,updown)
		#updowngenes<-GOgenes
		GOtermGenes[[i]]<-updowngenes
	}
	names(GOtermGenes)<-allResTable[,1]
	gogenesDF<-ldply(GOtermGenes,cbind)
	colnames(gogenesDF)<-c("GO.ID","UNIgene")
	return(gogenesDF)
}

allResBPsig<-allResBP[as.numeric(allResBP[,"weight"])<GOpval,]

sigGOgenesClust<-getGOgenes(allResBPsig,myInterestingGenes,GOdataBP)
sigGOgenesClust[,1]<-as.factor(sigGOgenesClust[,1])

###convert GO clusters to a list
sigGeneGOlist<-list()
for(i in 1:length(allResBPsig[,1])){
	sigGeneGOlist[[i]]<-as.character(sigGOgenesClust[grep(allResBPsig[i,1],sigGOgenesClust[,1]),2])
	
}

allResBPsig<-allResBPsig[!duplicated(sigGeneGOlist),]
sigGeneGOlist<-sigGeneGOlist[!duplicated(sigGeneGOlist)]

###function to get names for each GO category
getGOnames<-function(allResTable){
	geneGOnames<-vector()
	for (i in 1:length(allResTable[,1])){
		cmd<-paste("attributes(GOtermList$`",allResTable[i,1],"`)$Term",sep="")
		GOnames<-eval(parse(text=cmd))
		geneGOnames[i]<-paste(allResTable[i,1],GOnames,sep="_")
	}
	return(geneGOnames)
}

names(sigGeneGOlist)<-make.names(getGOnames(allResBPsig),unique=T)

###add category for genes with no GO annotation########3
noannot<-character()
for(i in 1:length(geneID2GO)){
	if(length(geneID2GO[[i]])<1){
		geneID2GO[[i]]<-"GO:9999999"
		noannot<-c(noannot,names(geneID2GO[i]))
	}
}
noannot.sig<-intersect(noannot,sigGenes)
if(length(noannot.sig)>0){
noannotRes<-c("GO:9999999","noGOannotation",length(noannot),length(noannot.sig),"<0.05","<0.05")
noannotClust<-cbind(rep("GO:9999999",length(noannot.sig)),noannot.sig)
colnames(noannotClust)<-c("GO.ID","UNIgene")

sigGOgenesClust<-getGOgenes(allResBPsig,myInterestingGenes,GOdataBP)
sigGOgenesClust<-rbind(sigGOgenesClust,noannotClust)
sigGOgenesClust[,1]<-as.factor(sigGOgenesClust[,1])

allResBPsig<-rbind(allResBPsig,noannotRes)
sigGeneGOlist[[length(allResBPsig[,1])]]<-noannot.sig
names(sigGeneGOlist)[length(allResBPsig[,1])]<-"GO.9999999_noGOannot"
}
############################



###find how often each protein appears in different GO categories
gene.promiscuity<-ldply(geneID2GO,cbind)
geneWeights<-table(gene.promiscuity[,1])




########score GO categories per organism using weighted scoring scheme
clustscore.all<-matrix(nrow=length(sigGeneGOlist), ncol=2)
clust.filt.all<-list()
clustorgscores.all<-list()
pb <- txtProgressBar(min = 0, max = length(sigGeneGOlist), style = 3)
for (i in 1:length(sigGeneGOlist)){
mysums<-list()
clustgenes<-sigGeneGOlist[[i]]
clustgenescT<-presAbsAll.all[match(clustgenes,rownames(presAbsAll.all)),]

for (j in 1:length(clustgenescT)){
	mysums[[j]]<-sum(1/as.numeric(geneWeights[na.omit(match(clustgenes[which(clustgenescT[[j]]>0)],names(geneWeights)))]))/sum(1/as.numeric(geneWeights[na.omit(match(clustgenes,names(geneWeights)))]))
}
names(mysums)<-names(clustgenescT)

clustorgscores.all[[i]]<-mysums
clustscore.all[i,1]<-i
clustscore.all[i,2]<-length(clustgenescT)

setTxtProgressBar(pb, i)
}

clusts.all<-matrix(nrow=length(clustorgscores.all),ncol=length(orgsNames))

for (i in 1:length(clusts.all[,1])){
	clusts.all[i,1:length(clusts.all[1,])]<-as.numeric(clustorgscores.all[[i]][1:length(clustorgscores.all[[i]])])
}

colnames(clusts.all)<-orgsNames

clustnames.all<-names(sigGeneGOlist)
for(i in 1:length(sigGeneGOlist)){
	clustnames.all[i]<-paste(names(sigGeneGOlist)[i],"_",length(sigGeneGOlist[[i]]),sep="")
}

clustnames.all<-replace(clustnames.all,is.na(clustnames.all),"NoName")

clusts.all<-as.data.frame(clusts.all)
rownames(clusts.all)<-make.names(clustnames.all,unique=T)
clusts.all<-as.matrix(clusts.all)

hmdf.all<-clusts.all

genesFunctList<-sigGeneGOlist
names(genesFunctList)<-make.names(names(genesFunctList), unique=TRUE)
hmdfc<-hmdf.all
rownames(hmdfc)<-make.names(rownames(hmdfc), unique=TRUE)


###Find set of predictive GO categories (each gene is predictive in a way, but they may have been recombined in a way that makes the functional category non-predictive over all!)
clustRowsdf<-hmdfc

dim(clustRowsdf)
###get columns for each group
Gr1.GO<-vector(length=length(orgsNamesGr1))
for(i in 1:length(orgsNamesGr1)){
Gr1.GO[i]<-grep(orgsNamesGr1[i],colnames(clustRowsdf))
}

Gr2.GO<-vector(length=length(orgsNamesGr2))
for(i in 1:length(orgsNamesGr2)){
Gr2.GO[i]<-grep(orgsNamesGr2[i],colnames(clustRowsdf))
}

###reorder data frame by training and test genomes
clustRowsdf.trn<-clustRowsdf[,c(Gr1.GO,Gr2.GO)]
clustRowsdf.tst<-clustRowsdf[,-c(Gr1.GO,Gr2.GO)]

clustRowsdf<-cbind(clustRowsdf.trn,clustRowsdf.tst)

###get columns for each group in re-ordered data frame
Gr1.GO<-vector(length=length(orgsNamesGr1))
for(i in 1:length(orgsNamesGr1)){
Gr1.GO[i]<-grep(orgsNamesGr1[i],colnames(clustRowsdf))
}

Gr2.GO<-vector(length=length(orgsNamesGr2))
for(i in 1:length(orgsNamesGr2)){
Gr2.GO[i]<-grep(orgsNamesGr2[i],colnames(clustRowsdf))
}


###make training matrix
learndat<-t(clustRowsdf[,c(Gr1.GO,Gr2.GO)])
learndat2<-cbind(c(rep("a",length(c(Gr1.GO))),rep("b",length(c(Gr2.GO)))),learndat)
learndat2<-as.data.frame(learndat2,stringsAsFactors=F)
colnames(learndat2)[1]<-"Category"
for(i in 2:length(learndat2[1,])){learndat2[,i]<-as.numeric(as.character(learndat2[,i]))}

####function for plotting Boruta output (from XXX)
make.plots <- function(b, num,
                       true.var = NA,
                       main = paste("Boruta feature selection for test", num)) {
    write.text <- function(b, true.var) {
        if ( !is.na(true.var) ) {
            text(1, max(attStats(b)$meanZ), pos = 4,
                 labels = paste("True vars are V.1-V.",
                     true.var, sep = ""))        
        }
    }
    plot(b, main = main, las = 3, xlab = "")
    write.text(b, true.var)
    png(paste(run.name, num, "png", sep = "."), width = 8, height = 8,
        units = "cm", res = 300, pointsize = 4)
    plot(b, main = main, lwd = 0.5, las = 3, xlab = "")
    write.text(b, true.var)
    dev.off()
}

###run the Boruta algorithm **10000 times here.
btst <- Boruta(learndat2[,2:length(learndat2[1,])],factor(learndat2[,1]), doTrace = 2, maxRuns=as.numeric(opt[4]))
print(btst)
#plot(btst)

###check out the leftover categories... but don't use this!
btst2<-TentativeRoughFix(btst, averageOver = Inf)
print(btst2)

###keep the predictive categories
BorutaCats<-c(names(btst$finalDecision[btst$finalDecision=="Confirmed"]))

###in case you decide to use those fishy categories
#BorutaCats<-c(names(btst2$finalDecision[btst2$finalDecision=="Confirmed"]))

clustRowsdfplot<-clustRowsdf[rownames(clustRowsdf) %in% BorutaCats,]

##create test sets
tsttheseNames<-c(colnames(clustRowsdf.tst))
tsttheseRnum<-c(match(colnames(clustRowsdf.tst),colnames(clustRowsdfplot)))
tsttheseRnum<-na.omit(tsttheseRnum)
testSet<-t(clustRowsdfplot[,c(tsttheseRnum)])

###make the gene and GO category heatmaps
bcats.GOs<-do.call(rbind,strsplit(BorutaCats,"_"))[,1]
geneListMatches <- unique(grep(paste(bcats.GOs,collapse="|"),names(sigGeneGOlist), value=TRUE))
bcats.genes<-unique(unlist(sigGeneGOlist[geneListMatches]))
plotHeat.genes<-presAbsAll.all[bcats.genes,]
plotHeat<-clustRowsdfplot
plotHeat<-na.omit(plotHeat)
plotHeat.genes<-na.omit(plotHeat.genes)

###heirarchical clustering using the manhattan metric of columns (organisms)
d <- dist(t(plotHeat),method="manhattan") # distance matrix
fitC <- hclust(d, method="ward.D2")
d.genes <- dist(t(plotHeat.genes),method="manhattan") # distance matrix
fitC.genes <- hclust(d.genes, method="ward.D2")

###order by max category scores across organisms in Gr1
plotHeat<-clustRowsdfplot[order(-(rowSums(clustRowsdfplot[,c(Gr1.GO)]))),]
plotHeat.genes<-plotHeat.genes[order(-(rowSums(plotHeat.genes[,c(Gr1.GO)]))),]

###plot heatmap for GO categories matrix
plotHeat<-as.matrix(plotHeat)
my_palette <- colorRampPalette(brewer.pal(5, "Spectral"))(n = 1000)
pdf(GOheat.name, width = 55, height = 35 )
heatmap.2(plotHeat, Rowv=F,Colv=as.dendrogram(fitC),trace="none",density.info="none",dendrogram="column",col=my_palette,srtCol=45,margins=c(20,60),keysize=0.5,symkey=F,symbreaks=F, cexRow=3, cexCol=3)
dev.off()

###plot heatmap for proteins/genes matrix
plotHeat.genes<-as.matrix(plotHeat.genes)
my_palette <- colorRampPalette(brewer.pal(5, "Spectral"))(n = 1000)
pdf(geneHeat.name, width = 55, height = 35 )
heatmap.2(plotHeat.genes, Rowv=T,Colv=as.dendrogram(fitC.genes),trace="none",density.info="none",dendrogram="column",col=my_palette,srtCol=45,margins=c(20,60),keysize=0.5,symkey=F,symbreaks=F, cexRow=3, cexCol=3)
dev.off()


BorutaCats2<-strsplit(BorutaCats,"_")
sigGeneGOlist2<-list()
for(i in 1:length(BorutaCats)){
	sigGeneGOlist2[[i]]<-as.character(sigGOgenesClust[grep(BorutaCats2[[i]][1],sigGOgenesClust[,1]),2])
	
}
names(sigGeneGOlist2)<-BorutaCats2

Gr1.gns<-unique(as.character(unlist(sigGeneGOlist2)))
Gr1.bcats<-BorutaCats

writeLines(Gr1.gns,genes.name)
writeLines(Gr1.bcats,GOs.name)

###########################################################################################################################################################
####make pca plot

###make training matrix
learndat<-t(clustRowsdfplot[,c(Gr1.GO,Gr2.GO)])
learndat2<-cbind(c(rep("a",length(c(Gr1.GO))),rep("b",length(c(Gr2.GO)))),learndat)
learndat2<-as.data.frame(learndat2,stringsAsFactors=F)
colnames(learndat2)[1]<-"Category"
for(i in 2:length(learndat2[1,])){learndat2[,i]<-as.numeric(as.character(learndat2[,i]))}

###principal compenents on training data
prclust<-learndat2
prclust<-prclust[,-1]
prcompFeedModel<-prcomp(prclust,scale=F)
bestcols<-brewer.pal(9,"Set1")

if(length(clustRowsdfplot[1,])>dim(learndat2)[1]){
prclust2<-t(clustRowsdfplot[,c((dim(learndat2)[1]+1):length(clustRowsdfplot[1,]))])
prcompFeedModel2<-(scale(prclust2, prcompFeedModel$center, prcompFeedModel$scale) %*% prcompFeedModel$rotation)
allpca2<-as.data.frame(prcompFeedModel2)
allpca2$groups<-rep("Test",length(allpca2[,1]))
mylabels2 <- rownames(allpca2)
}
           
allpca<-(scale(prclust, prcompFeedModel$center, prcompFeedModel$scale)%*%prcompFeedModel$rotation)
allpca<-as.data.frame(allpca)

cate2<-c(rep("Group 1",length(c(Gr1.GO))),rep("Group 2",length(c(Gr2.GO))))
allpca$groups<-cate2
if(exists("allpca2")){
allpca<-rbind(allpca,allpca2)
}else{
allpca<-rbind(allpca)
}

mylabels<-rownames(allpca)
mylabels[1:length(cate2)]<-""


PC1var = paste("PC1 (",100*round(as.numeric(summary(prcompFeedModel)$importance[,1][2]),2),"%)",sep="")
PC2var = paste("PC2 (",100*round(as.numeric(summary(prcompFeedModel)$importance[,2][2]),2),"%)",sep="")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


if(length(grep("Test",allpca$groups))>0){
allpca$groups <- factor(allpca$groups, levels=c("Group 1","Group 2","Test"), labels=c("Group 1","Group 2","Test genome"))
pcaplot <- ggplot(data=allpca)
pcaplot <- pcaplot + geom_point(aes(x=PC1,y=PC2,shape=factor(groups,labels=c("Group 1","Group 2","Test genome")),fill=factor(groups,labels=c("Group 1","Group 2","Test genome"))),size=3) 
pcaplot <- pcaplot + scale_shape_manual(name = "Groups",labels=c("Group 1","Group 2","Test genome"), values = c("Group 1"=21,"Group 2"=24,"Test genome"=23))  
pcaplot <- pcaplot + scale_fill_manual(name = "Groups",labels=c("Group 1","Group 2","Test genome"), values = c("Group 1"=cbPalette[1],"Group 2"=cbPalette[4], "Test genome"=cbPalette[3]))
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Group 1",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[1],lty=2) 
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Group 2",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[4],lty=2) 
pcaplot <- pcaplot + geom_vline(xintercept = 0,lty=2)+ geom_hline(yintercept = 0,lty=2) 
pcaplot <- pcaplot + labs(x=PC1var, y=PC2var) + guides(shape = guide_legend(override.aes = list(size=3))) 
pcaplot <- pcaplot + geom_text_repel(data=allpca,aes(x=PC1,y=PC2,label=mylabels),fontface="italic",point.padding = unit(0.8, "lines"),box.padding = unit(0.4, "lines"),nudge_y =0.2, colour="black")
pcaplot <- pcaplot + theme_minimal()
}else{
allpca$groups <- factor(allpca$groups, levels=c("Group 1","Group 2"), labels=c("Group 1","Group 2"))
pcaplot <- ggplot(data=allpca)
pcaplot <- pcaplot + geom_point(aes(x=PC1,y=PC2,shape=factor(groups,labels=c("Group 1","Group 2")),fill=factor(groups,labels=c("Group 1","Group 2"))),size=3) 
pcaplot <- pcaplot + scale_shape_manual(name = "Groups",labels=c("Group 1","Group 2"), values = c("Group 1"=21,"Group 2"=24))  
pcaplot <- pcaplot + scale_fill_manual(name = "Groups",labels=c("Group 1","Group 2"), values = c("Group 1"=cbPalette[1],"Group 2"=cbPalette[4]))
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Group 1",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[1],lty=2) 
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Group 2",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[4],lty=2) 
pcaplot <- pcaplot + geom_vline(xintercept = 0,lty=2)+ geom_hline(yintercept = 0,lty=2) 
pcaplot <- pcaplot + labs(x=PC1var, y=PC2var) + guides(shape = guide_legend(override.aes = list(size=3))) 
pcaplot <- pcaplot + theme_minimal()
}
#print(pcaplot)

ggsave(pcaname, plot = pcaplot, device = "pdf",width=15,height=10)


learndat3<-learndat2
learndat3$Category<-as.character(learndat3$Category)
learndat3$Category[learndat2$Category=="a"]<-1
learndat3$Category[learndat2$Category=="b"]<-0
learndat3$Category<-factor(learndat3$Category)

if(exists("allpca2")){

testDat<-as.data.frame(t(clustRowsdfplot[,-c(Gr1.GO,Gr2.GO)]))

if(opt[8]=="pnn"){
pnn <- learn(learndat3)
pnn <- smooth(pnn,sigma=1)
predictions<-apply(testDat,1,function(x){
	guess(pnn,as.numeric(x))$probabilities})
var_hat=matrix(nrow=length(testDat[,1]),ncol=2)
for(i in 1:length(predictions[1,])){
var_hat[i,]<-predictions[,i]
}
rownames(var_hat)<-rownames(testDat)
colnames(var_hat)<-c(0,1)
}

if(opt[8]=="rf"){
rfmodel<-randomForest(as.factor(Category)~.,data=learndat3,ntree=10000) 
var_hat<-predict(rfmodel, testDat,type="prob")
}

}



rfpreds<-matrix(nrow=length(testDat[,1]),ncol=3)
rfpreds[,1]<-rep(rownames(testDat),1)
rfpreds[1:length(var_hat[,2]),2]<-"Gr1 prediction"
rfpreds[1:length(var_hat[,2]),3]<-var_hat[,2]
rfpreds<-as.data.frame(rfpreds)
colnames(rfpreds)<-c("Species","prediction","probability")
rfpreds[,3]<-as.numeric(as.character(rfpreds[,3]))

rfpreds[,2]  <- factor(rfpreds[,2] , levels=c("Gr1 prediction"))

p1 <- ggplot(rfpreds, aes(Species,probability,fill=prediction)) + geom_col(position = position_dodge()) +coord_flip() + guides(fill = guide_legend(reverse = TRUE))
p1<- p1 + geom_abline(intercept=0.5,slope=0,col="red",lty=2) +
          xlab("species") +
          ylab("Prediction Probability") +
          ggtitle("Gr1 prediction") +
		  scale_fill_grey(start = 0.2, end = 0.8) + theme_classic()
#plot(p1)

ggsave(predplot, plot = p1, device = "pdf")

write.table(rfpreds,predTable,row.names=F,quote=F)

