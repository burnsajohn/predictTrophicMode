
pdf("modelOUTPUT/defaultModelOUTPUT/Predictions_MollweideProj.pdf",useDingbats=FALSE,width=7,height=7)
source("scripts/trophicMode_mollweideProject.r")
x<-read.table("modelOUTPUT/defaultModelOUTPUT/tables/predictionsDataTable.txt",sep="\t",header=T)
y<-x[,c(1,2,6,5)] ###column order for predictions table is: 1) column 1: species name, 2) column 2: phagocyte generalist prediction, 3) column 6: photosynthesis prediction, 4) column 5: prototroph prediction
y<-as.data.frame(y,stringsAsFactors=F)
y[,2:4] = map_df(y[,2:4], as.numeric)
y$nopred <- 1-(rowSums(y[2:4]))/3
y.df<-y[2:5]
colnames(y.df)<-c("u","s","m","l")
tcs.pred<-tcspace(y.df)
preds<-mollweidePoints(tcs.pred)
points(preds[[1]]~preds[[2]],pch=20)
addTextLabels(preds[[2]],preds[[1]], labels=seq(1:length(preds[[1]])))
dev.off()
pdf("modelOUTPUT/defaultModelOUTPUT/Legend_MollweideProj.pdf",useDingbats=FALSE,width=7,height=7)
plot.new()
legend("center",legend=paste(1:length(x[,1]),x[,1],sep=". "))
dev.off()
