# predictTrophicMode

Trophic Mode Prediction Tool, (c) John A. Burns, 2017

This tool runs in R, using the rscript command line tool, which should be platform independent.

##################################################################################################
Dependencies:
To run this tool in default mode, matching the analysis done in the corresponding manuscript, you must install the following packages in your R distribution:

gplots
ggplot2
ggrepel
matrixStats
RColorBrewer
randomForest
rlist
optparse
pnn
basicPlotteR
pavo
concaveman
mapproj
purrr

To run this tool in advanced mode, with additional functionality for comparing between any two groups and limiting the analysis to the proteins in your favorite genome/s, you must install the following additional packages in your R distribution:

topGO
plyr
Boruta

########################################################################################################
#
#Usage:
#
#For either mode, the first step is to search the ~14k HMMs against the proteins from your genome or transcriptome. The HMMs included with the model (the file with all 14k HMMs is here: "TrophicModePredictionTool/HMMs/phag_nonphag-allVall-any3diverse.hmmCAT.hmm") will only work with proteins. The script expects a list of models that have a significant hit (evalue<=1e-5) to some protein in your genome/transcriptomes of interest. To get such a list, run the following commands on a compter/server with HMMER3 installed.
#
#define files, names
	
	species=[your species here]
	filename=[name of your protein file here]
	HMMs=[path to]/phag_nonphag-allVall-any3diverse.hmmCAT.hmm
#
#run hmmsearch
#
	hmmsearch --tblout $species.x.phag_nonphag-allVall-any3diverse.hmmsearchOUT-tbl.txt --cpu 8 $HMMs $filename
#
#filter for significant hits (evalue <= 1e-5) *evalue threshold used to generate the model was <=1e-5, must use same threshold here. If you would like to experiment with changing that threshold, you must rebuild the presense absence matrix by re-thresholding all of the default organisms. To do so, run hmmsearch against all of those genomes, or contact the author for the raw hmmsearch output files. 
#
	
#
#the following commands can be used to get the list of signifcant hits in the correct format for the prediction tool (in BASH):
#
	sigfile=${species}_sigHits.txt
	sigModel=${sigfile//_sigHits.txt/_sigModels.txt}
	echo $sigfile
	echo $sigModel
	grep -v "^#" $species.x.phag_nonphag-allVall-any3diverse.hmmsearchOUT-tbl.txt | awk '$5<=1e-5 && $8<=1e-4' > $sigfile
	awk '{print $3}' $sigfile | sort -u > $sigModel
#
#the file with the format "species_sigModels.txt" should contain a list of HMM names that are significant hits to at least 1 protein in your target genome. 
#
#############################################################################################################
#
#Place the file/s containing the list of significant HMMs into the directory "TrophicModePredictionTool/TestGenomes"
#
#In Windows (tested in Windows PowerShell) or Linux, cd to the directory "TrophicModePredictionTool" 
#
#Default mode, testing against the models outlined in the corresponding manuscript:
#
#run the script: runModel.r
#
#In Windows PowerShell it might be run like so (from within TrophicModePredictionTool):
#
	'C:\Program Files\R\R-3.3.1\bin\Rscript.exe' .\runModel.r
#
#In Linux, it is much the same: [path to]Rscript runModel.r
#
#If the model is run exactly as above, it will finish in about 30 seconds and will output a series of plots and predictions into the directory "TrophicModePredictionTool/modelOUTPUT/defaultModelOUTPUT" and its subdirectories. The outputs are organized by type of analysis.
#
#
#*NOTE* on Mollweide projection: Legend is in the file: "modelOUTPUT/defaultModelOUTPUT/Legend_MollweideProj.pdf"
#Shaded region legend is in the file: "modelOUTPUT/defaultModelOUTPUT/MollweideShading.pdf"
#To make a publication-ready figure, combine the Mollweide projection plot with the information from the two legend files.
#
#############################################################################################################
#
#Advanced mode:
#
#Advanced mode allows the user control over the groupings, and the ability to limit the model to proteins present in a single or set of genomes. 
#
#To change the groupings, edit the files: "AdvancedMode/orgGroups/Group1.txt" and/or "Group2.txt" ***Note: The species names in the files Group1.txt and Group2.txt must exactly match (case sensitive) the species names given to the model.
#
#The advanced model will look for proteins and GO categories relatively enriched in Group1 and depleted in Group2.
#
#The model comes with a presence/absence matrix containing set of default genomes, to remove any of them from the analysis entirely, include them in the file: AdvancedMode/orgGroups/removeOrgs.txt ***Note, if they are removed from the analysis, they also need to be removed from the training set files. ***Note2: must have exact matches here as well.
#
#To limit the model to proteins present in a single genome or set of genomes (For example if it is a small genome, but you are certain it has the character you are interested in), include those species names in the file: "AdvancedMode/orgGroups/limitOrgs.txt". Those species will automatically be removed from the training sets.
#
#############################################################################################################################################################
#
#to see a list of options, use:
#
	Rscript runModel.r --help
#
#That will print the following:
Options:
        -m MODE, --mode=MODE
                use in "manuscript" (fast) or "advanced" (slower) mode? Manuscript mode runs the model as outlined in the paper. Advanced mode allows comparisons using the 14k HMMs discussed in the paper on any two groups of organisms. [default= manuscript]

        -l LIMIT, --limit=LIMIT
                limit to specific genome/set of genomes Y or N? [default= N]

        -p PHAGOSOME, --phagosome=PHAGOSOME
                include proteins from phagosome Y or N? [default= N]

        -b BORUTA, --boruta=BORUTA
                how many times to iterate Boruta algorithm? For testing, 100 is a good number [default= 10000]

        -r REMOVE_ORGANISMS, --remove_organisms=REMOVE_ORGANISMS
                remove some or all default organisms? Y or N ***Note: cannot remove default organisms if they are used in the training set!!! default= N]

        --gene_pval=GENE_PVAL
                pvalue threshold for difference between group1 and group2 for each gene [default= 0.05]

        --GO_pval=GO_PVAL
                pvalue for GO term enrichment [default= 0.2]

        -h, --help
                Show this help message and exit

###############################################################################
#
#Running advanced mode:
#
#use the flags listed above to change the different options. None are required. If the default files are unaltered and no flags are given other than "-m advanced", the script will regenerate the model from the manuscript. If the files are altered, or any flags are used, the model will change. The "-p" flag will force the model to include proteins identified as part of the physical phagosome, regardless of whether or not they are differentially found between group 1 and group 2. Advanced mode takes between 3-10 minutes to run. ***Note, in Windows PowerShell, I've noticed that the script sometimes gets paused somehow and it helps to occasionally hit enter if it looks like it is stuck.
#
#As in "manuscript" (default) mode, the output for advanced mode is organized by type of analysis and appears in the directory "TrophicModePredictionTool/modelOUTPUT/advancedModeOUTPUT" and its subdirectories.
#
#ENJOY!
#################################################################################

