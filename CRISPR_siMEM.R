suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("blme"))
suppressPackageStartupMessages(library("doMC"))
suppressPackageStartupMessages(library("genefilter"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("locfit"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("preprocessCore"))
suppressPackageStartupMessages(library("reshape"))
source("/Users/rkumar/mac_icr/work/colt2/mem_code/data_format_lib.R")
source("/Users/rkumar/mac_icr/work/colt2/mem_code/model_lib.R")
source("/Users/rkumar/mac_icr/work/colt2/mem_code/simem_lib.R")



# this data can contain all zero counts for a gene - bad
exprsFile<-"CRISPR_expression_set.txt"
exprs<-as.matrix(read.table(exprsFile,header=TRUE,sep="\t",row.names=1,as.is=TRUE))



# get gene names for each row
exprs_genes <- gsub("_.+","",rownames(exprs),perl=TRUE)
unique_genes <- unique(exprs_genes)

# sum the counts for all guides of each gene
gene_sum_counts <- NULL
for(gene in unique_genes){
	gene_sum_counts <- c(
		gene_sum_counts,
		sum(
			exprs[which(exprs_genes == gene),]
			)
		)
}
names(gene_sum_counts) <- unique_genes

# exprs_clean contains all genes with > 0 counts
exprs_clean <- exprs[-(which(
	exprs_genes %in% names(gene_sum_counts)[which(gene_sum_counts == 0)]
	)),]

exprs_clean <- log10(exprs_clean+0.5)

pData="CRISPR_pdata.txt"
pdata<-read.table(pData,header=TRUE,sep="\t",row.names=1,as.is=TRUE)


fData="CRISPR_fdata.txt"
fdata<-read.table(fData,header=TRUE,sep="\t",row.names=1,as.is=TRUE)

# need to remove genes from fdata that were also removed from exprs_clean
# Note that this assumes rows (and order) are equal in exprs and fdata
fdata_clean <- fdata[-(which(
	exprs_genes %in% names(gene_sum_counts)[which(gene_sum_counts == 0)]
	)),]

featureData<-AnnotatedDataFrame(data=fdata_clean)

phenoData<-AnnotatedDataFrame(data=pdata)
crispr_screens<-ExpressionSet(exprs_clean,featureData=featureData,phenoData=phenoData)
hp=read.delim("guide_annotation.txt",header=T,as.is=T,check.names=F)
hpWeights=hp[,c("trcn_id","gene_id","weight")]
subtypes=read.delim("class.txt",header=T,as.is=T,check.names=F)
status=subtypes[,c("cell_line","drug")]
status$drug=ifelse(status$drug=="1","high","low")
#genesOfInterest=c(1, 2, 3)
results=simem(
	screens=crispr_screens,
	#geneIds=genesOfInterest,
	covariate="drug",
	reagentWeights=hpWeights,
	annotationsPerCellLine=status,
	inverseVarianceWeights=FALSE,
	signalProbWeights=FALSE,
	analyzeReagents=FALSE,
	covariateFactorOrder=c("low","high"),
	parallelNodes=1
	)

# Writing output
result <- results$gene
write.csv(result,file="result")
