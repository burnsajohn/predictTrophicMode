#!/usr/bin Rscript

print("Initializing R")

suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(matrixStats))
suppressMessages(library(RColorBrewer))
suppressMessages(library(randomForest))
suppressMessages(library(rlist))
suppressMessages(library(optparse))
suppressMessages(library(pnn))
suppressMessages(library(basicPlotteR))
suppressMessages(library(pavo))
suppressMessages(library(concaveman))
suppressMessages(library(mapproj))
suppressMessages(library(purrr))

option_list = list(
  make_option(c("-m", "--mode"), type="character", default="manuscript", 
              help="use in \"manuscript\" (fast) or \"advanced\" (slower) mode? Manuscript mode runs the model as outlined in the paper. Advanced mode allows comparisons using the 14k HMMs discussed in the paper on any two groups of organisms. [default= %default]"),
  make_option(c("-l", "--limit"), type="character", default="N", 
              help="limit to specific genome/set of genomes Y or N? [default= %default]"),
  make_option(c("-p", "--phagosome"), type="character", default="N", 
              help="include proteins from phagosome Y or N? [default= %default]"),
  make_option(c("-b", "--boruta"), type="numeric", default=10000, 
              help="how many times to iterate Boruta algorithm? For testing, 100 is a good number [default= %default]"),
  make_option(c("-r", "--remove_organisms"), type="character", default="N", 
              help="remove some or all default organisms? Y or N ***Note: cannot remove default organisms if they are used in the training set!!! default= %default]"),
  make_option(c("--gene_pval"), type="numeric", default=0.05, 
              help="pvalue threshold for difference between group1 and group2 for each gene [default= %default]"),
  make_option(c("--GO_pval"), type="numeric", default=0.2, 
              help="pvalue for GO term enrichment [default= %default]"),
  make_option(c("--machine"), type="character", default="pnn", 
              help="which type of learning algorithm should be used to predict outcomes? pnn or rf [default= %default]")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#opt<-c("advanced","Y","N",1000,"N",0.15,0.2,"pnn")

genepval<-as.numeric(opt[6])
GOpval<-as.numeric(opt[7])

if(opt[1]=="manuscript"){
	print("Done Initializing")
	source("scripts/predictTrophicMode_Tool.r")

}else if(opt[1]=="advanced"){
	
	suppressMessages(library(topGO))
	suppressMessages(library(plyr))
	suppressMessages(library(Boruta))
	print("Done Initializing")
	print(genepval)
	print(GOpval)
	source("scripts/predictTrophicMode_Flex.r")
	
}


