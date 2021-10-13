library(gplots)
library(ggplot2)
library(dplyr)


load("SAMPLES43.RData")

#############################################
load("hisat2.cufflk.OI.n42.RData")
colnames(fpkmMat)<- gsub("_[0-9].*","",basename(colnames(fpkmMat)))

load("htseq-count.gtf.bulk_RNAseq_2021.n43.RData")
colnames(dat.counts)<- gsub("_[0-9].*.gtf.counts.gz","",basename(colnames(dat.counts)))

table(colnames(dat.counts)==colnames(fpkmMat))
table(gsub("-.*","",colnames(dat.counts))==SAMPLES43$RNAsample)
table(colnames(dat.counts)==SAMPLES43$RNAsample)

#############################################

SAMPLES<-readr::read_tsv("SAMPLES-30.txt")
IDS<-as.character(SAMPLES$sampleID)
ind43<-match(gsub("[-_].*","",colnames(dat.counts)),IDS)

Sillence<-readr::read_tsv("P:/HKU-SZH//Sillence_subtyping/Summary.of.Sillence.and.Geno.Sep29.2021.tsv")
ind30b<-match(gsub(" .*$","",SAMPLES$PatientID),Sillence$HospitalIDs)
Sillence30<-Sillence[ind30b,]
IDS<-as.character(SAMPLES$sampleID)
ind43<-match(gsub("[-_].*","",colnames(dat.counts)),IDS)
Sillence43<-Sillence30[ind43,]
names(Sillence43)[5]<-paste0(names(Sillence43)[5],"2")

################################################
################################################
load("F:/OneDrive - The University Of Hong Kong/resources/TF.CD.receptors.human.mouse.RData")

library(DESeq2)
source("step6.help.funcs.R")

indNZ<-which(rowSums(dat.counts>0)>0)
ANNOT<-annot.gtf[indNZ,]
##########
	indBat2Col<-which(SAMPLES43$Batch==2 & (grepl("COL|pending",SAMPLES43$Mutation)|SAMPLES43$Index%in%c(14,15, 29, 30))&!SAMPLES43$Index%in%c(20) &colSums(dat.counts)>1e7)
	SAMPLES9<-SAMPLES43[indBat2Col,]
	FPKM9<-fpkmMat[,indBat2Col]

	SAMPLES43[indBat2Col,]
	Sillence43[indBat2Col,]
	SAMPLES43[["COLvsCtrl"]]<-"COL1A1"
	SAMPLES43[["COLvsCtrl"]][grepl("Healthy",SAMPLES43$Mutation)]<-"Ctrl"
	table(SAMPLES43[["COLvsCtrl"]][indBat2Col])

	DDS_batch12<-DESEQ(indBat2Col,"COLvsCtrl")
	resultsNames(DDS_batch12) # lists the coefficients
	res.Bat12<- lfcShrink(DDS_batch12, coef="COLvsCtrl_Ctrl_vs_COL1A1", type="apeglm")
	rownames(res.Bat12)<-ANNOT[,4]
##########
##########
save(DDS_batch12, res.Bat12, file="DEGs.batch2.2Col1.vs.2Ctrls.RData")

pdf("Volcano.plots.of.Batch2.BulkHumanOI.Col1-vs-Ctrl.Oct-9.n9.pdf",width=16,height=12)
	
	DatDEG<-data.frame(res.Bat12,ANNOT)%>%filter(padj<0.01 & abs(log2FoldChange)>2)%>%
		mutate(Direction=factor(c("higher in Col1","","Higher in Ctrl")[sign(log2FoldChange)+2]))
	LG<-split(DatDEG$gene_name,DatDEG$Direction)
	boxplot(FPKM9[genes$gene_short_name%in%LG[[1]],], outline=F)
	boxplot(FPKM9[genes$gene_short_name%in%LG[[2]],], outline=F)

	barplot(FPKM9[genes$gene_short_name=="COL10A1",])
	barplot(FPKM9[genes$gene_short_name=="IHH",])

	DatDEGtf<-DatDEG[DatDEG[,"gene_name"]%in%TF.human,]

	ggplot(as.data.frame(res.Bat12)%>%filter(abs(log2FoldChange)>0.1),
		aes(x=-log2FoldChange,y=-log10(pvalue))) +
		geom_point(color="grey") +
		geom_point(data=DatDEG,aes(x=-log2FoldChange,y=-log10(pvalue),color=Direction),size=2) +
		geom_text(data=DatDEG,
			aes(x=-log2FoldChange,y=-log10(pvalue),label=gene_name))+
		ggtitle("Contrl vs Col1")

	

dev.off()



