

if(pred==1){

pcaname<-"modelOUTPUT/defaultModelOUTPUT/PCAplots/PhagocyteGOs-PCAplot.pdf"

prclust.ph<-learndat2
prclust.ph<-prclust.ph[,-1]
prcompFeedModel<-prcomp(prclust.ph,scale=F)
bestcols<-brewer.pal(9,"Set1")

if(length(clustRowsdfplot[1,])>dim(learndat2)[1]){
prclust2<-t(clustRowsdfplot[,c((dim(learndat2)[1]+1):length(clustRowsdfplot[1,]))])
rownames(prclust2)<-colnames(clustRowsdfplot)[c((dim(learndat2)[1]+1):length(clustRowsdfplot[1,]))]
##project new genome data to pca space defined by training data
prcompFeedModel2<-(scale(prclust2, prcompFeedModel$center, prcompFeedModel$scale) %*% prcompFeedModel$rotation)
##
allpca2<-as.data.frame(prcompFeedModel2)
allpca2$groups<-rep("Test",length(allpca2[,1]))
mylabels2 <- rownames(allpca2)

pcaPM2<-c(grep("Cym. tetramitiformis",rownames(allpca2)),grep("Bi. natans",rownames(allpca2)),grep("Chr. tobin",rownames(allpca2)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM2<-c(pcaPM2,grep(orgsNamesPM[i],rownames(allpca2)))
	}
}
if(length(pcaPM2)>0){allpca2<-allpca2[-pcaPM2,]}
}
           
pcaPM<-c(grep("Cym. tetramitiformis",colnames(clustRowsdfplot)),grep("Bi. natans",colnames(clustRowsdfplot)),grep("Chr. tobin",colnames(clustRowsdfplot)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM<-c(pcaPM,grep(orgsNamesPM[i],colnames(clustRowsdfplot)))
	}
}
pcaPM<-unique(pcaPM)
prclust.pm<-t(clustRowsdfplot[,pcaPM])
#prclust.pm<-prclust.pm[,na.omit(colnames(prclust))]
prcompFeedModel.pm<-(scale(prclust.pm, prcompFeedModel$center, prcompFeedModel$scale) %*% prcompFeedModel$rotation)
allpca.pm<-as.data.frame(prcompFeedModel.pm)
allpca.pm$groups<-rep("Phagomixotroph",length(allpca.pm[,1]))
mylabels.pm <- rownames(allpca.pm)




allpca<-(scale(prclust.ph, prcompFeedModel$center, prcompFeedModel$scale)%*%prcompFeedModel$rotation)
allpca<-as.data.frame(allpca)
pcaPM3<-c(grep("Cym. tetramitiformis",rownames(allpca)),grep("Bi. natans",rownames(allpca)),grep("Chr. tobin",rownames(allpca)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM3<-c(pcaPM3,grep(orgsNamesPM[i],rownames(allpca)))
	}
}
if(length(pcaPM3>0)){allpca<-allpca[-pcaPM3,]}

cate2<-c(rep("Phagocytes",length(c(PH))),rep("Non-phagocytes",length(c(PL,GA,FNall,RA))))
allpca$groups<-cate2
if(exists("allpca2")){
allpca<-rbind(allpca,allpca.pm,allpca2)
}else{
allpca<-rbind(allpca,allpca.pm)
}

mylabels<-rownames(allpca)
mylabels[1:(length(cate2)+length(allpca.pm[,1]))]<-""


PC1var = paste("PC1 (",100*round(as.numeric(summary(prcompFeedModel)$importance[,1][2]),2),"%)",sep="")
PC2var = paste("PC2 (",100*round(as.numeric(summary(prcompFeedModel)$importance[,2][2]),2),"%)",sep="")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


if(length(grep("Test",allpca$groups))>0){
allpca$groups <- factor(allpca$groups, levels=c("Phagocytes","Phagomixotroph","Non-phagocytes","Test"), labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome"))
pcaplot <- ggplot(data=allpca)
pcaplot <- pcaplot + geom_point(aes(x=PC1,y=PC2,shape=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome")),fill=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome"))),size=5) 
pcaplot <- pcaplot + scale_shape_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome"), values = c("Phagocytotic eukaryote"=21,"Phagomixotrophic eukaryote"=22, "Non-phagocytotic eukaryote"=24,  "Test genome"=23))  
pcaplot <- pcaplot + scale_fill_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome"), values = c("Phagocytotic eukaryote"=cbPalette[1],"Phagomixotrophic eukaryote"=cbPalette[2],"Non-phagocytotic eukaryote"=cbPalette[4],  "Test genome"=cbPalette[3]))
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[1],lty=2) 
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Non-phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[4],lty=2) 
pcaplot <- pcaplot + geom_vline(xintercept = 0,lty=2)+ geom_hline(yintercept = 0,lty=2) 
pcaplot <- pcaplot + labs(x=PC1var, y=PC2var) + guides(shape = guide_legend(override.aes = list(size=5))) 
#pcaplot <- pcaplot + geom_text_repel(data=allpca2,aes(x=PC1,y=PC2,label=rownames(allpca2)),fontface="italic",point.padding = unit(0.8, "lines"),box.padding = unit(0.4, "lines"),nudge_y =0.2, colour="black") 
pcaplot <- pcaplot + geom_text_repel(data=allpca,aes(x=PC1,y=PC2,label=mylabels),fontface="italic",point.padding = unit(0.8, "lines"),box.padding = unit(0.4, "lines"),nudge_y =0.2, colour="black")
pcaplot <- pcaplot + theme_minimal() + theme(legend.text=element_text(size=12))
}else{
allpca$groups <- factor(allpca$groups, levels=c("Phagocytes","Phagomixotroph","Non-phagocytes"), labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote"))
pcaplot <- ggplot(data=allpca)
pcaplot <- pcaplot + geom_point(aes(x=PC1,y=PC2,shape=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote")),fill=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote"))),size=5) 
pcaplot <- pcaplot + scale_shape_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote"), values = c("Phagocytotic eukaryote"=21,"Phagomixotrophic eukaryote"=22, "Non-phagocytotic eukaryote"=24))  
pcaplot <- pcaplot + scale_fill_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote"), values = c("Phagocytotic eukaryote"=cbPalette[1],"Phagomixotrophic eukaryote"=cbPalette[2],"Non-phagocytotic eukaryote"=cbPalette[4]))
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[1],lty=2) 
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Non-phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[4],lty=2) 
pcaplot <- pcaplot + geom_vline(xintercept = 0,lty=2)+ geom_hline(yintercept = 0,lty=2) 
pcaplot <- pcaplot + labs(x=PC1var, y=PC2var) + guides(shape = guide_legend(override.aes = list(size=5))) 
pcaplot <- pcaplot + theme_minimal() + theme(legend.text=element_text(size=12))
}
#print(pcaplot)

ggsave(pcaname, plot = pcaplot, device = "pdf",width=15,height=10)


learndat3<-learndat2
learndat3$Category<-as.character(learndat3$Category)
learndat3$Category[learndat2$Category=="a"]<-1
learndat3$Category[learndat2$Category=="b"]<-0
learndat3$Category<-factor(learndat3$Category)

if(exists("allpca2")){

testDat.ph<-as.data.frame(t(clustRowsdfplot[,-c(PH,PM,PL,GA,FNall,RA)]))
rownames(testDat.ph)<-colnames(clustRowsdfplot)[c((dim(learndat2)[1]+1):length(clustRowsdfplot[1,]))]


if(opt[8]=="pnn"){
pnn <- learn(learndat3)
pnn <- smooth(pnn,sigma=1)
predictions<-apply(testDat.ph,1,function(x){
	guess(pnn,as.numeric(x))$probabilities})
var_hat.ph=matrix(nrow=length(testDat.ph[,1]),ncol=2)
for(i in 1:length(predictions[1,])){
var_hat.ph[i,]<-predictions[,i]
}
rownames(var_hat.ph)<-rownames(testDat.ph)
colnames(var_hat.ph)<-c(0,1)
}

if(opt[8]=="rf"){
rfmodel.ph<-randomForest(as.factor(Category)~.,data=learndat3,ntree=10000) 
var_hat.ph<-predict(rfmodel.ph, testDat.ph,type="prob")
}


}
}


if(pred==2){

pcaname<-"modelOUTPUT/defaultModelOUTPUT/PCAplots/PhotosynthesisGOs-PCAplot.pdf"




learndat<-t(clustRowsdfplot[,c(PH,FNall,PM,PL,GA,RA)])
learndat2<-cbind(c(rep("a",length(c(PH,FNall))),rep("b",length(c(PM,PL,GA,RA)))),learndat)
learndat2<-as.data.frame(learndat2,stringsAsFactors=F)
colnames(learndat2)[1]<-"Category"
for(i in 2:length(learndat2[1,])){learndat2[,i]<-as.numeric(as.character(learndat2[,i]))}
colnames(learndat2)[1]<-"Category"
for(i in 2:length(learndat2[1,])){learndat2[,i]<-as.numeric(as.character(learndat2[,i]))}


###principal compenents on training data
prclust.ph<-learndat2
prclust.ph<-prclust.ph[,-1]
prcompFeedModel<-prcomp(prclust.ph,scale=F)
bestcols<-brewer.pal(9,"Set1")

if(length(clustRowsdfplot[1,])>dim(learndat2)[1]){
prclust2<-t(clustRowsdfplot[,c((dim(learndat2)[1]+1):length(clustRowsdfplot[1,]))])
rownames(prclust2)<-colnames(clustRowsdfplot)[c((dim(learndat2)[1]+1):length(clustRowsdfplot[1,]))]
prcompFeedModel2<-(scale(prclust2, prcompFeedModel$center, prcompFeedModel$scale) %*% prcompFeedModel$rotation)
allpca2<-as.data.frame(prcompFeedModel2)
allpca2$groups<-rep("Test",length(allpca2[,1]))
mylabels2 <- rownames(allpca2)
pcaPM2<-c(grep("Cym. tetramitiformis",rownames(allpca2)),grep("Bi. natans",rownames(allpca2)),grep("Chr. tobin",rownames(allpca2)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM2<-c(pcaPM2,grep(orgsNamesPM[i],rownames(allpca2)))
	}
}
if(length(pcaPM2)>0){allpca2<-allpca2[-pcaPM2,]}
}        
		  
pcaPM<-c(grep("Cym. tetramitiformis",colnames(clustRowsdfplot)),grep("Bi. natans",colnames(clustRowsdfplot)),grep("Chr. tobin",colnames(clustRowsdfplot)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM<-c(pcaPM,grep(orgsNamesPM[i],colnames(clustRowsdfplot)))
	}
}
pcaPM<-unique(pcaPM)
prclust.pm<-t(clustRowsdfplot[,pcaPM])
#prclust.pm<-prclust.pm[,na.omit(colnames(prclust))]
prcompFeedModel.pm<-(scale(prclust.pm, prcompFeedModel$center, prcompFeedModel$scale) %*% prcompFeedModel$rotation)
allpca.pm<-as.data.frame(prcompFeedModel.pm)
allpca.pm$groups<-rep("Phagomixotroph",length(allpca.pm[,1]))
mylabels.pm <- rownames(allpca.pm)




allpca<-(scale(prclust.ph, prcompFeedModel$center, prcompFeedModel$scale)%*%prcompFeedModel$rotation)
allpca<-as.data.frame(allpca)
pcaPM3<-c(grep("Cym. tetramitiformis",rownames(allpca)),grep("Bi. natans",rownames(allpca)),grep("Chr. tobin",rownames(allpca)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM3<-c(pcaPM3,grep(orgsNamesPM[i],rownames(allpca)))
	}
}
if(length(pcaPM3>0)){allpca<-allpca[-pcaPM3,]}

cate2<-c(rep("No photosynthesis",length(c(PH,FNall))),rep("Photosynthesis",length(c(PL,GA,RA))))
allpca$groups<-cate2
if(exists("allpca2")){
allpca<-rbind(allpca,allpca.pm,allpca2)
}else{
allpca<-rbind(allpca,allpca.pm)
}
mylabels<-rownames(allpca)
mylabels[1:(length(cate2)+3)]<-""


if(length(grep("Test",allpca$groups))>0){
allpca$groups <- factor(allpca$groups, levels=c("No photosynthesis","Phagomixotroph","Photosynthesis","Test"), labels=c("Non-photosynthetic eukaryote","Phagomixotrophic eukaryote","Photosynthetic eukaryote", "Test genome"))
pcaplot <- ggplot(data=allpca)
pcaplot <- pcaplot + geom_point(aes(x=PC1,y=PC2,shape=factor(groups,labels=c("Non-photosynthetic eukaryote","Phagomixotrophic eukaryote","Photosynthetic eukaryote", "Test genome")),fill=factor(groups,labels=c("Non-photosynthetic eukaryote","Phagomixotrophic eukaryote","Photosynthetic eukaryote", "Test genome"))),size=5) 
pcaplot <- pcaplot + scale_shape_manual(name = "Groups",labels=c("Non-photosynthetic eukaryote","Phagomixotrophic eukaryote","Photosynthetic eukaryote", "Test genome"), values = c("Non-photosynthetic eukaryote"=21,"Phagomixotrophic eukaryote"=22, "Photosynthetic eukaryote"=24,  "Test genome"=23))  
pcaplot <- pcaplot + scale_fill_manual(name = "Groups",labels=c("Non-photosynthetic eukaryote","Phagomixotrophic eukaryote","Photosynthetic eukaryote", "Test genome"), values = c("Non-photosynthetic eukaryote"=cbPalette[1],"Phagomixotrophic eukaryote"=cbPalette[2],"Photosynthetic eukaryote"=cbPalette[4],  "Test genome"=cbPalette[3]))
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Non-photosynthetic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[1],lty=2) 
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Photosynthetic eukaryote" | allpca$groups=="Phagomixotrophic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[4],lty=2) 
pcaplot <- pcaplot + geom_vline(xintercept = 0,lty=2)+ geom_hline(yintercept = 0,lty=2) 
pcaplot <- pcaplot + labs(x=PC1var, y=PC2var) + guides(shape = guide_legend(override.aes = list(size=5))) 
pcaplot <- pcaplot + geom_text_repel(data=allpca,aes(x=PC1,y=PC2,label=mylabels),fontface="italic",point.padding = unit(0.8, "lines"),box.padding = unit(0.4, "lines"),nudge_y =0.2, colour="black")
pcaplot <- pcaplot + theme_minimal() + theme(legend.text=element_text(size=12))
}else{
allpca$groups <- factor(allpca$groups, levels=c("No photosynthesis","Phagomixotroph","Photosynthesis"), labels=c("Non-photosynthetic eukaryote","Phagomixotrophic eukaryote","Photosynthetic eukaryote"))
pcaplot <- ggplot(data=allpca)
pcaplot <- pcaplot + geom_point(aes(x=PC1,y=PC2,shape=factor(groups,labels=c("Non-photosynthetic eukaryote","Phagomixotrophic eukaryote","Photosynthetic eukaryote")),fill=factor(groups,labels=c("Non-photosynthetic eukaryote","Phagomixotrophic eukaryote","Photosynthetic eukaryote"))),size=5) 
pcaplot <- pcaplot + scale_shape_manual(name = "Groups",labels=c("Non-photosynthetic eukaryote","Phagomixotrophic eukaryote","Photosynthetic eukaryote"), values = c("Non-photosynthetic eukaryote"=21,"Phagomixotrophic eukaryote"=22, "Photosynthetic eukaryote"=24))  
pcaplot <- pcaplot + scale_fill_manual(name = "Groups",labels=c("Non-photosynthetic eukaryote","Phagomixotrophic eukaryote","Photosynthetic eukaryote"), values = c("Non-photosynthetic eukaryote"=cbPalette[1],"Phagomixotrophic eukaryote"=cbPalette[2],"Photosynthetic eukaryote"=cbPalette[4]))
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Non-photosynthetic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[1],lty=2) 
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Photosynthetic eukaryote" | allpca$groups=="Phagomixotrophic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[4],lty=2) 
pcaplot <- pcaplot + geom_vline(xintercept = 0,lty=2)+ geom_hline(yintercept = 0,lty=2) 
pcaplot <- pcaplot + labs(x=PC1var, y=PC2var) + guides(shape = guide_legend(override.aes = list(size=5))) 
pcaplot <- pcaplot + theme_minimal() + theme(legend.text=element_text(size=12))
}
#print(pcaplot)

ggsave(pcaname, plot = pcaplot, device = "pdf",width=15,height=10)


learndat3<-learndat2
learndat3$Category<-as.character(learndat3$Category)
learndat3$Category[learndat2$Category=="a"]<-0
learndat3$Category[learndat2$Category=="b"]<-1
learndat3$Category<-factor(learndat3$Category)

if(length(orgsNamesPM)>0){
	PM<-vector(length=length(orgsNamesPM))
	for(i in 1:length(orgsNamesPM)){
	PM[i]<-grep(orgsNamesPM[i],colnames(clustRowsdfplot))
	}
}else{PM<-numeric()}

if(exists("allpca2")){

testDat.photo<-as.data.frame(t(clustRowsdfplot[,-c(PH,PM,PL,GA,FNall,RA)]))

if(opt[8]=="pnn"){
pnn <- learn(learndat3)
pnn <- smooth(pnn,sigma=1)
predictions<-apply(testDat.photo,1,function(x){
	guess(pnn,as.numeric(x))$probabilities})
var_hat.photo=matrix(nrow=length(testDat.photo[,1]),ncol=2)
for(i in 1:length(predictions[1,])){
var_hat.photo[i,]<-predictions[,i]
}
rownames(var_hat.photo)<-rownames(testDat.photo)
colnames(var_hat.photo)<-c(0,1)
}

if(opt[8]=="rf"){
rfmodel.photo<-randomForest(as.factor(Category)~.,data=learndat3,ntree=10000) 
var_hat.photo<-predict(rfmodel.photo, testDat.photo,type="prob")
}

}
}


if(pred==3){

pcaname<-"modelOUTPUT/defaultModelOUTPUT/PCAplots/PrototrophGOs-PCAplot.pdf"

###principal compenents on training data
prclust.ph<-learndat2
prclust.ph<-prclust.ph[,-1]
prcompFeedModel<-prcomp(prclust.ph,scale=F)
bestcols<-brewer.pal(9,"Set1")

if(length(clustRowsdfplot[1,])>dim(learndat2)[1]){
prclust2<-t(clustRowsdfplot[,c((dim(learndat2)[1]+1):length(clustRowsdfplot[1,]))])
rownames(prclust2)<-colnames(clustRowsdfplot)[c((dim(learndat2)[1]+1):length(clustRowsdfplot[1,]))]
prcompFeedModel2<-(scale(prclust2, prcompFeedModel$center, prcompFeedModel$scale) %*% prcompFeedModel$rotation)
allpca2<-as.data.frame(prcompFeedModel2)
allpca2$groups<-rep("Test",length(allpca2[,1]))
mylabels2 <- rownames(allpca2)
pcaPM2<-c(grep("Cym. tetramitiformis",rownames(allpca2)),grep("Bi. natans",rownames(allpca2)),grep("Chr. tobin",rownames(allpca2)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM2<-c(pcaPM2,grep(orgsNamesPM[i],rownames(allpca2)))
	}
}
if(length(pcaPM2)>0){allpca2<-allpca2[-pcaPM2,]}
}
            
pcaPM<-c(grep("Cym. tetramitiformis",colnames(clustRowsdfplot)),grep("Bi. natans",colnames(clustRowsdfplot)),grep("Chr. tobin",colnames(clustRowsdfplot)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM<-c(pcaPM,grep(orgsNamesPM[i],colnames(clustRowsdfplot)))
	}
}
pcaPM<-unique(pcaPM)
prclust.pm<-t(clustRowsdfplot[,pcaPM])
#prclust.pm<-prclust.pm[,na.omit(colnames(prclust))]
prcompFeedModel.pm<-(scale(prclust.pm, prcompFeedModel$center, prcompFeedModel$scale) %*% prcompFeedModel$rotation)
allpca.pm<-as.data.frame(prcompFeedModel.pm)
allpca.pm$groups<-rep("Phagomixotroph",length(allpca.pm[,1]))
mylabels.pm <- rownames(allpca.pm)




allpca<-(scale(prclust.ph, prcompFeedModel$center, prcompFeedModel$scale)%*%prcompFeedModel$rotation)
allpca<-as.data.frame(allpca)
pcaPM3<-c(grep("Cym. tetramitiformis",rownames(allpca)),grep("Bi. natans",rownames(allpca)),grep("Chr. tobin",rownames(allpca)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM3<-c(pcaPM3,grep(orgsNamesPM[i],rownames(allpca)))
	}
}
if(length(pcaPM3>0)){allpca<-allpca[-pcaPM3,]}

cate2<-c(rep("Phagocytes",length(c(PH))),rep("Prototrophs",length(c(PL,GA,FNall,RA))))
allpca$groups<-cate2
if(exists("allpca2")){
allpca<-rbind(allpca,allpca.pm,allpca2)
}else{
allpca<-rbind(allpca,allpca.pm)
}
mylabels<-rownames(allpca)
mylabels[1:(length(cate2)+3)]<-""



PC1var = paste("PC1 (",100*round(as.numeric(summary(prcompFeedModel)$importance[,1][2]),2),"%)",sep="")
PC2var = paste("PC2 (",100*round(as.numeric(summary(prcompFeedModel)$importance[,2][2]),2),"%)",sep="")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

if(length(grep("Test",allpca$groups))>0){
allpca$groups <- factor(allpca$groups, levels=c("Phagocytes","Phagomixotroph","Prototrophs","Test"), labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Prototrophic eukaryote", "Test genome"))
pcaplot <- ggplot(data=allpca)
pcaplot <- pcaplot + geom_point(aes(x=PC1,y=PC2,shape=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Prototrophic eukaryote", "Test genome")),fill=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Prototrophic eukaryote", "Test genome"))),size=5) 
pcaplot <- pcaplot + scale_shape_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Prototrophic eukaryote", "Test genome"), values = c("Phagocytotic eukaryote"=21,"Phagomixotrophic eukaryote"=22, "Prototrophic eukaryote"=24,  "Test genome"=23))  
pcaplot <- pcaplot + scale_fill_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Prototrophic eukaryote", "Test genome"), values = c("Phagocytotic eukaryote"=cbPalette[1],"Phagomixotrophic eukaryote"=cbPalette[2],"Prototrophic eukaryote"=cbPalette[4],  "Test genome"=cbPalette[3]))
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[1],lty=2) 
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Prototrophic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[4],lty=2) 
pcaplot <- pcaplot + geom_vline(xintercept = 0,lty=2)+ geom_hline(yintercept = 0,lty=2) 
pcaplot <- pcaplot + labs(x=PC1var, y=PC2var) + guides(shape = guide_legend(override.aes = list(size=5))) 
#pcaplot <- pcaplot + geom_text_repel(data=allpca2,aes(x=PC1,y=PC2,label=rownames(allpca2)),fontface="italic",point.padding = unit(0.8, "lines"),box.padding = unit(0.4, "lines"),nudge_y =0.2, colour="black") 
pcaplot <- pcaplot + geom_text_repel(data=allpca,aes(x=PC1,y=PC2,label=mylabels),fontface="italic",point.padding = unit(0.8, "lines"),box.padding = unit(0.4, "lines"),nudge_y =0.2, colour="black")
pcaplot <- pcaplot + theme_minimal() + theme(legend.text=element_text(size=12))
}else{
allpca$groups <- factor(allpca$groups, levels=c("Phagocytes","Phagomixotroph","Non-phagocytes"), labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Prototrophic eukaryote"))
pcaplot <- ggplot(data=allpca)
pcaplot <- pcaplot + geom_point(aes(x=PC1,y=PC2,shape=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Prototrophic eukaryote")),fill=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Prototrophic eukaryote"))),size=5) 
pcaplot <- pcaplot + scale_shape_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Prototrophic eukaryote"), values = c("Phagocytotic eukaryote"=21,"Phagomixotrophic eukaryote"=22, "Prototrophic eukaryote"=24))  
pcaplot <- pcaplot + scale_fill_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Prototrophic eukaryote"), values = c("Phagocytotic eukaryote"=cbPalette[1],"Phagomixotrophic eukaryote"=cbPalette[2],"Prototrophic eukaryote"=cbPalette[4]))
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[1],lty=2) 
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Prototrophic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[4],lty=2) 
pcaplot <- pcaplot + geom_vline(xintercept = 0,lty=2)+ geom_hline(yintercept = 0,lty=2) 
pcaplot <- pcaplot + labs(x=PC1var, y=PC2var) + guides(shape = guide_legend(override.aes = list(size=5))) 
pcaplot <- pcaplot + theme_minimal() + theme(legend.text=element_text(size=12))
}
#print(pcaplot)

ggsave(pcaname, plot = pcaplot, device = "pdf",width=15,height=10)


learndat3<-learndat2
learndat3$Category<-as.character(learndat3$Category)
learndat3$Category[learndat2$Category=="a"]<-0
learndat3$Category[learndat2$Category=="b"]<-1
learndat3$Category<-factor(learndat3$Category)

if(length(orgsNamesPM)>0){
	PM<-vector(length=length(orgsNamesPM))
	for(i in 1:length(orgsNamesPM)){
	PM[i]<-grep(orgsNamesPM[i],colnames(clustRowsdfplot))
	}
}else{PM<-numeric()}

if(exists("allpca2")){

testDat.nph<-as.data.frame(t(clustRowsdfplot[,-c(PH,PM,PL,GA,FNall,RA)]))

if(opt[8]=="pnn"){
pnn <- learn(learndat3)
pnn <- smooth(pnn,sigma=1)
predictions<-apply(testDat.nph,1,function(x){
	guess(pnn,as.numeric(x))$probabilities})
var_hat.nph=matrix(nrow=length(testDat.nph[,1]),ncol=2)
for(i in 1:length(predictions[1,])){
var_hat.nph[i,]<-predictions[,i]
}
rownames(var_hat.nph)<-rownames(testDat.nph)
colnames(var_hat.nph)<-c(0,1)
}

if(opt[8]=="rf"){
rfmodel.nph<-randomForest(as.factor(Category)~.,data=learndat3,ntree=10000) 
var_hat.nph<-predict(rfmodel.nph, testDat.nph,type="prob")
}

}

}


if(pred==4){

pcaname<-"modelOUTPUT/defaultModelOUTPUT/PCAplots/Entamoebid.PhagocyteGOs-PCAplot.pdf"

prclust.ph<-learndat2
prclust.ph<-prclust.ph[,-1]
prcompFeedModel<-prcomp(prclust.ph,scale=F)
bestcols<-brewer.pal(9,"Set1")

if(length(clustRowsdfplot[1,])>dim(learndat2)[1]){
prclust2<-t(clustRowsdfplot[,c((dim(learndat2)[1]+1):length(clustRowsdfplot[1,]))])
rownames(prclust2)<-colnames(clustRowsdfplot)[c((dim(learndat2)[1]+1):length(clustRowsdfplot[1,]))]
prcompFeedModel2<-(scale(prclust2, prcompFeedModel$center, prcompFeedModel$scale) %*% prcompFeedModel$rotation)
allpca2<-as.data.frame(prcompFeedModel2)
allpca2$groups<-rep("Test",length(allpca2[,1]))
mylabels2 <- rownames(allpca2)

pcaPM2<-c(grep("Cym. tetramitiformis",rownames(allpca2)),grep("Bi. natans",rownames(allpca2)),grep("Chr. tobin",rownames(allpca2)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM2<-c(pcaPM2,grep(orgsNamesPM[i],rownames(allpca2)))
	}
}
if(length(pcaPM2)>0){allpca2<-allpca2[-pcaPM2,]}
}
           
pcaPM<-c(grep("Cym. tetramitiformis",colnames(clustRowsdfplot)),grep("Bi. natans",colnames(clustRowsdfplot)),grep("Chr. tobin",colnames(clustRowsdfplot)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM<-c(pcaPM,grep(orgsNamesPM[i],colnames(clustRowsdfplot)))
	}
}
pcaPM<-unique(pcaPM)
prclust.pm<-t(clustRowsdfplot[,pcaPM])
#prclust.pm<-prclust.pm[,na.omit(colnames(prclust))]
prcompFeedModel.pm<-(scale(prclust.pm, prcompFeedModel$center, prcompFeedModel$scale) %*% prcompFeedModel$rotation)
allpca.pm<-as.data.frame(prcompFeedModel.pm)
allpca.pm$groups<-rep("Phagomixotroph",length(allpca.pm[,1]))
mylabels.pm <- rownames(allpca.pm)




allpca<-(scale(prclust.ph, prcompFeedModel$center, prcompFeedModel$scale)%*%prcompFeedModel$rotation)
allpca<-as.data.frame(allpca)
pcaPM3<-c(grep("Cym. tetramitiformis",rownames(allpca)),grep("Bi. natans",rownames(allpca)),grep("Chr. tobin",rownames(allpca)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM3<-c(pcaPM3,grep(orgsNamesPM[i],rownames(allpca)))
	}
}
if(length(pcaPM3>0)){allpca<-allpca[-pcaPM3,]}

cate2<-c(rep("Phagocytes",length(c(PH))),rep("Non-phagocytes",length(c(PL,GA,FNall,RA))))
allpca$groups<-cate2
if(exists("allpca2")){
allpca<-rbind(allpca,allpca.pm,allpca2)
}else{
allpca<-rbind(allpca,allpca.pm)
}

mylabels<-rownames(allpca)
mylabels[1:(length(cate2)+length(allpca.pm[,1]))]<-""


PC1var = paste("PC1 (",100*round(as.numeric(summary(prcompFeedModel)$importance[,1][2]),2),"%)",sep="")
PC2var = paste("PC2 (",100*round(as.numeric(summary(prcompFeedModel)$importance[,2][2]),2),"%)",sep="")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


if(length(grep("Test",allpca$groups))>0){
allpca$groups <- factor(allpca$groups, levels=c("Phagocytes","Phagomixotroph","Non-phagocytes","Test"), labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome"))
pcaplot <- ggplot(data=allpca)
pcaplot <- pcaplot + geom_point(aes(x=PC1,y=PC2,shape=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome")),fill=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome"))),size=5) 
pcaplot <- pcaplot + scale_shape_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome"), values = c("Phagocytotic eukaryote"=21,"Phagomixotrophic eukaryote"=22, "Non-phagocytotic eukaryote"=24,  "Test genome"=23))  
pcaplot <- pcaplot + scale_fill_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome"), values = c("Phagocytotic eukaryote"=cbPalette[1],"Phagomixotrophic eukaryote"=cbPalette[2],"Non-phagocytotic eukaryote"=cbPalette[4],  "Test genome"=cbPalette[3]))
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[1],lty=2) 
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Non-phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[4],lty=2) 
pcaplot <- pcaplot + geom_vline(xintercept = 0,lty=2)+ geom_hline(yintercept = 0,lty=2) 
pcaplot <- pcaplot + labs(x=PC1var, y=PC2var) + guides(shape = guide_legend(override.aes = list(size=5))) 
#pcaplot <- pcaplot + geom_text_repel(data=allpca2,aes(x=PC1,y=PC2,label=rownames(allpca2)),fontface="italic",point.padding = unit(0.8, "lines"),box.padding = unit(0.4, "lines"),nudge_y =0.2, colour="black") 
pcaplot <- pcaplot + geom_text_repel(data=allpca,aes(x=PC1,y=PC2,label=mylabels),fontface="italic",point.padding = unit(0.8, "lines"),box.padding = unit(0.4, "lines"),nudge_y =0.2, colour="black")
pcaplot <- pcaplot + theme_minimal() + theme(legend.text=element_text(size=12))
}else{
allpca$groups <- factor(allpca$groups, levels=c("Phagocytes","Phagomixotroph","Non-phagocytes"), labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote"))
pcaplot <- ggplot(data=allpca)
pcaplot <- pcaplot + geom_point(aes(x=PC1,y=PC2,shape=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote")),fill=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote"))),size=5) 
pcaplot <- pcaplot + scale_shape_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote"), values = c("Phagocytotic eukaryote"=21,"Phagomixotrophic eukaryote"=22, "Non-phagocytotic eukaryote"=24))  
pcaplot <- pcaplot + scale_fill_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote"), values = c("Phagocytotic eukaryote"=cbPalette[1],"Phagomixotrophic eukaryote"=cbPalette[2],"Non-phagocytotic eukaryote"=cbPalette[4]))
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[1],lty=2) 
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Non-phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[4],lty=2) 
pcaplot <- pcaplot + geom_vline(xintercept = 0,lty=2)+ geom_hline(yintercept = 0,lty=2) 
pcaplot <- pcaplot + labs(x=PC1var, y=PC2var) + guides(shape = guide_legend(override.aes = list(size=5))) 
pcaplot <- pcaplot + theme_minimal() + theme(legend.text=element_text(size=12))
}
#print(pcaplot)

ggsave(pcaname, plot = pcaplot, device = "pdf",width=15,height=10)


learndat3<-learndat2
learndat3$Category<-as.character(learndat3$Category)
learndat3$Category[learndat2$Category=="a"]<-1
learndat3$Category[learndat2$Category=="b"]<-0
learndat3$Category<-factor(learndat3$Category)

if(exists("allpca2")){

testDat.ent.ph<-as.data.frame(t(clustRowsdfplot[,-c(PH,PM,PL,GA,FNall,RA)]))

if(opt[8]=="pnn"){
pnn <- learn(learndat3)
pnn <- smooth(pnn,sigma=1)
predictions<-apply(testDat.ent.ph,1,function(x){
	guess(pnn,as.numeric(x))$probabilities})
var_hat.ent.ph=matrix(nrow=length(testDat.ent.ph[,1]),ncol=2)
for(i in 1:length(predictions[1,])){
var_hat.ent.ph[i,]<-predictions[,i]
}
rownames(var_hat.ent.ph)<-rownames(testDat.ent.ph)
colnames(var_hat.ent.ph)<-c(0,1)
}

if(opt[8]=="rf"){
rfmodel.ent.ph<-randomForest(as.factor(Category)~.,data=learndat3,ntree=10000) 
var_hat.ent.ph<-predict(rfmodel.ent.ph, testDat.ent.ph,type="prob")
}

}
}


if(pred==5){

pcaname<-"modelOUTPUT/defaultModelOUTPUT/PCAplots/R.allomyces.PhagocyteGOs-PCAplot.pdf"

if(opt[8]=="pnn"){
predplot<-"modelOUTPUT/defaultModelOUTPUT/predictions-pnn.pdf"
}
if(opt[8]=="rf"){
predplot<-"modelOUTPUT/defaultModelOUTPUT/predictions-rf.pdf"
}

prclust.ph<-t(clustRowsdfplot[,c(PH,PM,PL,GA,FNall,RA)])
prcompFeedModel<-prcomp(prclust.ph,scale=F)
bestcols<-brewer.pal(9,"Set1")

if(length(clustRowsdfplot[1,])>dim(prclust.ph)[1]){
prclust2<-t(clustRowsdfplot[,c((dim(prclust.ph)[1]+1):length(clustRowsdfplot[1,]))])
rownames(prclust2)<-colnames(clustRowsdfplot)[c((dim(prclust.ph)[1]+1):length(clustRowsdfplot[1,]))]
prcompFeedModel2<-(scale(prclust2, prcompFeedModel$center, prcompFeedModel$scale) %*% prcompFeedModel$rotation)
allpca2<-as.data.frame(prcompFeedModel2)
allpca2$groups<-rep("Test",length(allpca2[,1]))
mylabels2 <- rownames(allpca2)

pcaPM2<-c(grep("Cym. tetramitiformis",rownames(allpca2)),grep("Bi. natans",rownames(allpca2)),grep("Chr. tobin",rownames(allpca2)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM2<-c(pcaPM2,grep(orgsNamesPM[i],rownames(allpca2)))
	}
}
if(length(pcaPM2)>0){allpca2<-allpca2[-pcaPM2,]}
}
           
pcaPM<-c(grep("Cym. tetramitiformis",colnames(clustRowsdfplot)),grep("Bi. natans",colnames(clustRowsdfplot)),grep("Chr. tobin",colnames(clustRowsdfplot)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM<-c(pcaPM,grep(orgsNamesPM[i],colnames(clustRowsdfplot)))
	}
}
pcaPM<-unique(pcaPM)
prclust.pm<-t(clustRowsdfplot[,pcaPM])
#prclust.pm<-prclust.pm[,na.omit(colnames(prclust))]
prcompFeedModel.pm<-(scale(prclust.pm, prcompFeedModel$center, prcompFeedModel$scale) %*% prcompFeedModel$rotation)
allpca.pm<-as.data.frame(prcompFeedModel.pm)
allpca.pm$groups<-rep("Phagomixotroph",length(allpca.pm[,1]))
mylabels.pm <- rownames(allpca.pm)




allpca<-(scale(prclust.ph, prcompFeedModel$center, prcompFeedModel$scale)%*%prcompFeedModel$rotation)
allpca<-as.data.frame(allpca)
pcaPM3<-c(grep("Cym. tetramitiformis",rownames(allpca)),grep("Bi. natans",rownames(allpca)),grep("Chr. tobin",rownames(allpca)))
if(length(orgsNamesPM)>0){
	for(i in 1:length(orgsNamesPM)){
		pcaPM3<-c(pcaPM3,grep(orgsNamesPM[i],rownames(allpca)))
	}
}
if(length(pcaPM3>0)){allpca<-allpca[-pcaPM3,]}

cate2<-c(rep("Phagocytes",length(c(PH))),rep("Non-phagocytes",length(c(PL,GA,FNall,RA))))
allpca$groups<-cate2
if(exists("allpca2")){
allpca<-rbind(allpca,allpca.pm,allpca2)
}else{
allpca<-rbind(allpca,allpca.pm)
}

mylabels<-rownames(allpca)
mylabels[1:(length(cate2)+length(allpca.pm[,1]))]<-""


PC1var = paste("PC1 (",100*round(as.numeric(summary(prcompFeedModel)$importance[,1][2]),2),"%)",sep="")
PC2var = paste("PC2 (",100*round(as.numeric(summary(prcompFeedModel)$importance[,2][2]),2),"%)",sep="")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


if(length(grep("Test",allpca$groups))>0){
allpca$groups <- factor(allpca$groups, levels=c("Phagocytes","Phagomixotroph","Non-phagocytes","Test"), labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome"))
pcaplot <- ggplot(data=allpca)
pcaplot <- pcaplot + geom_point(aes(x=PC1,y=PC2,shape=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome")),fill=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome"))),size=5) 
pcaplot <- pcaplot + scale_shape_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome"), values = c("Phagocytotic eukaryote"=21,"Phagomixotrophic eukaryote"=22, "Non-phagocytotic eukaryote"=24,  "Test genome"=23))  
pcaplot <- pcaplot + scale_fill_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote", "Test genome"), values = c("Phagocytotic eukaryote"=cbPalette[1],"Phagomixotrophic eukaryote"=cbPalette[2],"Non-phagocytotic eukaryote"=cbPalette[4],  "Test genome"=cbPalette[3]))
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[1],lty=2) 
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Non-phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[4],lty=2) 
pcaplot <- pcaplot + geom_vline(xintercept = 0,lty=2)+ geom_hline(yintercept = 0,lty=2) 
pcaplot <- pcaplot + labs(x=PC1var, y=PC2var) + guides(shape = guide_legend(override.aes = list(size=5))) 
#pcaplot <- pcaplot + geom_text_repel(data=allpca2,aes(x=PC1,y=PC2,label=rownames(allpca2)),fontface="italic",point.padding = unit(0.8, "lines"),box.padding = unit(0.4, "lines"),nudge_y =0.2, colour="black") 
pcaplot <- pcaplot + geom_text_repel(data=allpca,aes(x=PC1,y=PC2,label=mylabels),fontface="italic",point.padding = unit(0.8, "lines"),box.padding = unit(0.4, "lines"),nudge_y =0.2, colour="black")
pcaplot <- pcaplot + theme_minimal() + theme(legend.text=element_text(size=12))
}else{
allpca$groups <- factor(allpca$groups, levels=c("Phagocytes","Phagomixotroph","Non-phagocytes"), labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote"))
pcaplot <- ggplot(data=allpca)
pcaplot <- pcaplot + geom_point(aes(x=PC1,y=PC2,shape=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote")),fill=factor(groups,labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote"))),size=5) 
pcaplot <- pcaplot + scale_shape_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote"), values = c("Phagocytotic eukaryote"=21,"Phagomixotrophic eukaryote"=22, "Non-phagocytotic eukaryote"=24))  
pcaplot <- pcaplot + scale_fill_manual(name = "Groups",labels=c("Phagocytotic eukaryote","Phagomixotrophic eukaryote","Non-phagocytotic eukaryote"), values = c("Phagocytotic eukaryote"=cbPalette[1],"Phagomixotrophic eukaryote"=cbPalette[2],"Non-phagocytotic eukaryote"=cbPalette[4]))
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[1],lty=2) 
pcaplot <- pcaplot + stat_ellipse(data=allpca[allpca$groups=="Non-phagocytotic eukaryote",],aes(x=PC1,y=PC2), geom = "path",level = 0.95, segments = 51,colour=cbPalette[4],lty=2) 
pcaplot <- pcaplot + geom_vline(xintercept = 0,lty=2)+ geom_hline(yintercept = 0,lty=2) 
pcaplot <- pcaplot + labs(x=PC1var, y=PC2var) + guides(shape = guide_legend(override.aes = list(size=5))) 
pcaplot <- pcaplot + theme_minimal() + theme(legend.text=element_text(size=12))
}
#print(pcaplot)

ggsave(pcaname, plot = pcaplot, device = "pdf",width=15,height=10)


learndat3<-learndat2
learndat3$Category<-as.character(learndat3$Category)
learndat3$Category[learndat2$Category=="a"]<-1
learndat3$Category[learndat2$Category=="b"]<-0
learndat3$Category<-factor(learndat3$Category)

if(exists("allpca2")){

testDat.roz.ph<-as.data.frame(t(clustRowsdfplot[,-c(PH,PM,PL,GA,FNall,RA)]))

if(opt[8]=="pnn"){
pnn <- learn(learndat3)
pnn <- smooth(pnn,sigma=1)
predictions<-apply(testDat.roz.ph,1,function(x){
	guess(pnn,as.numeric(x))$probabilities})
var_hat.roz.ph=matrix(nrow=length(testDat.roz.ph[,1]),ncol=2)
for(i in 1:length(predictions[1,])){
var_hat.roz.ph[i,]<-predictions[,i]
}
rownames(var_hat.roz.ph)<-rownames(testDat.roz.ph)
colnames(var_hat.roz.ph)<-c(0,1)
}

if(opt[8]=="rf"){
rfmodel.roz.ph<-randomForest(as.factor(Category)~.,data=learndat3,ntree=10000) 
var_hat.roz.ph<-predict(rfmodel.roz.ph, testDat.roz.ph,type="prob")
}

mylen<-length(testDat.ph[,1])
rfpreds<-matrix(nrow=3*mylen,ncol=3)
rfpreds[,1]<-rep(rownames(testDat.ph),3)
rfpreds[1:mylen,2]<-"Phagocytosis prediction"
#rfpreds[(length(var_hat.ph[,2])+1):(2*length(var_hat.ph[,2])),2]<-"Phagocyte-entamoebid prediction"
#rfpreds[(2*length(var_hat.ph[,2])+1):(3*length(var_hat.ph[,2])),2]<-"Phagocyte-rozellid prediction"
rfpreds[(mylen+1):(2*mylen),2]<-"Prototrophy prediction"
rfpreds[(2*mylen+1):(3*mylen),2]<-"Photosynthesis prediction"

###get max phagocytosis prediction
phagpreds<-cbind(var_hat.ph[,2],var_hat.ent.ph[,2],var_hat.roz.ph[,2])
for (i in 1:mylen){
	colmax<-which(phagpreds[i,]==max(phagpreds[i,]))
	thispred<-max(phagpreds[i,])
	if(phagpreds[i,1]<0.5 & thispred>0.5){
		reportpred<-thispred
		if(colmax==2 | (phagpreds[i,2]>0.5 & phagpreds[i,3]>0.5)){
			rfpreds[c(i,(mylen+i),(2*mylen+i)),1]<-paste(rownames(testDat.ph)[i],"*",sep="")
		}
		if(colmax==3){
			rfpreds[c(i,(mylen+i),(2*mylen+i)),1]<-paste(rownames(testDat.ph)[i],"**",sep="")
		}
	}else{reportpred<-phagpreds[i,1]}
	rfpreds[i,3]<-reportpred
}

#rfpreds[1:length(var_hat.ph[,2]),3]<-var_hat.ph[,2]
#rfpreds[(length(var_hat.ph[,2])+1):(2*length(var_hat.ph[,2])),3]<-var_hat.ent.ph[,2]
#rfpreds[(2*length(var_hat.ph[,2])+1):(3*length(var_hat.ph[,2])),3]<-var_hat.roz.ph[,2]
rfpreds[(mylen+1):(2*mylen),3]<-var_hat.nph[,2]
rfpreds[(2*mylen+1):(3*mylen),3]<-var_hat.photo[,2]
rfpreds<-as.data.frame(rfpreds)
colnames(rfpreds)<-c("Species","prediction","probability")
rfpreds[,3]<-as.numeric(as.character(rfpreds[,3]))

#rfpreds$prediction  <- factor(rfpreds$prediction , levels=c("Prototrophy prediction","Phagocyte-generalist prediction","Phagocyte-entamoebid prediction","Phagocyte-rozellid prediction","Photosynthesis prediction"))
rfpreds$prediction  <- factor(rfpreds$prediction , levels=c("Prototrophy prediction","Phagocytosis prediction","Photosynthesis prediction"))
rfpreds$Species <- factor(rfpreds$Species,levels=rev(levels(rfpreds$Species)))

cols<-c("Prototrophy prediction"="gray","Phagocytosis prediction"="orange4","Photosynthesis prediction"="seagreen")

p1 <- ggplot(rfpreds, aes(Species,probability,fill=prediction)) + geom_col(position = position_dodge()) +coord_flip() + guides(fill = guide_legend(reverse = TRUE))
p1<- p1 + geom_abline(intercept=0.5,slope=0,col="red",lty=2) +
          xlab("species") +
          ylab("Prediction Probability") +
          ggtitle("Phagocytosis, prototrophy, and photosynthesis predictions") +
		  scale_fill_manual(values=cols) + theme_classic()
#plot(p1)

ggsave(predplot, plot = p1, device = "pdf")

write.table(rfpreds,"modelOUTPUT/defaultModelOUTPUT/tables/predictionsDataFrame.txt",row.names=F,quote=F,sep="\t")

predTable<-matrix(nrow=length(testDat.ph[,1]),ncol=6)
predTable[,1]<-rownames(testDat.ph)
predTable[,2]<-var_hat.ph[,2]
predTable[,3]<-var_hat.ent.ph[,2]
predTable[,4]<-var_hat.roz.ph[,2]
predTable[,5]<-var_hat.nph[,2]
predTable[,6]<-var_hat.photo[,2]
colnames(predTable)<-c("Genome","Phagocyte-generalist prediction","Phagocyte-entamoebid prediction","Phagocyte-rozellid prediction","Prototrophy prediction","Photosynthesis prediction")

write.table(predTable,"modelOUTPUT/defaultModelOUTPUT/tables/predictionsDataTable.txt",row.names=F,quote=F,sep="\t",col.name=T)

}
}

