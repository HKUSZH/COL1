library(dplyr)
library("xlsx")

Sillence<-readr::read_tsv("../../Sillence_subtyping/Summary.of.Sillence.and.Geno.Sep29.2021.tsv")
#SillenceCOL<-Sillence%>%filter(grepl("COL", AffectGenes))

NovelMut<-c("c.G415T", "c.G770T", "c.805-2A>G", "c.G923T", "c.1051delG", "c.G1517C", "c.G1408T", "c.1298_1299insT", "c.C1145A", "c.G4198C", "c.C3870A", "c.C3468G")

SillenceCOL<-Sillence%>%filter(AffectGenes%in%c("COL1A1", "COL1A2", "COL1A1, COL1A2"))
table(NovelMut%in%SillenceCOL$MutantAllele)
setdiff(NovelMut,SillenceCOL$MutantAllele)

SillenceCOL1<-SillenceCOL[grep("COL1A1",SillenceCOL$AffectGenes),]
SillenceCOL2<-SillenceCOL[grep("COL1A2",SillenceCOL$AffectGenes),]

COL1A1exons<-read.table("COL1A1-exons.txt")[-1,] 
COL1A1exons$V2<- COL1A1exons$V2-118
COL1A1exons$V3<- COL1A1exons$V3-118
COL1A1exons$V2[1]<-0
COL1A1exons$V3[nrow(COL1A1exons)]<-4500
COL1A1exons$V4<- -0.5
COL1A1exons$V5<- 0.5
COL1A1exons$xpos<-1
COL1A1exons$newPos<-1
COL1A1exons$height<-1

COL1A2exons<-read.table("COL1A2-exons.txt")
COL1A2exons$V2<- COL1A2exons$V2-137
COL1A2exons$V3<- COL1A2exons$V3-137
COL1A2exons$V2[1]<-0
COL1A2exons$V3[nrow(COL1A2exons)]<-4000
COL1A2exons$V4<- -0.5
COL1A2exons$V5<- 0.5
COL1A2exons$xpos<-1
COL1A2exons$newPos<-1
COL1A2exons$height<-1

KOLvec<-hue_pal()(length(table(SillenceCOL$MutationType)))
names(KOLvec)<-names(table(SillenceCOL$MutationType))
####################################
dataPOS_COL1<-data.frame(xpos=as.integer(gsub("[a-zA-Z]+$","",
			gsub("^[a-zA-Z]+","",gsub("[[:punct:]].*","",gsub("c.","",
				gsub("^.*\\(","",
					gsub("\\).*$","",SillenceCOL1$MutantAllele))))))),
	SillenceCOL1)%>%arrange(xpos)
dataPOS_COL1$height<-1
dataPOS_COL1$height[dataPOS_COL1$Novel=="YES"]<- -1
dataPOS_COL1$newPos<-seq(range(dataPOS_COL1$xpos)[1],range(dataPOS_COL1$xpos)[2],len=nrow(dataPOS_COL1))
dataPOS_COL1$height2<-dataPOS_COL1$height*4

as.numeric(dataPOS_COL1[,1])


dataPOS_COL2<-data.frame(xpos=as.integer(gsub("[a-zA-Z]+$","",
			gsub("^[a-zA-Z]+","",gsub("[[:punct:]].*","",gsub("c.","",
				gsub("\\).*$","",
					gsub("^.*\\(","",SillenceCOL2$MutantAllele))))))),
	SillenceCOL2)%>%arrange(xpos)
dataPOS_COL2$height<-1
dataPOS_COL2$height[dataPOS_COL2$Novel=="YES"]<- -1
dataPOS_COL2$newPos<-seq(range(dataPOS_COL2$xpos)[1],range(dataPOS_COL2$xpos)[2],len=nrow(dataPOS_COL2))
dataPOS_COL2$height2<-dataPOS_COL2$height*4

as.numeric(dataPOS_COL2[,1])

####################################
library(ggpubr)

pdf("Mutation.Pos.Chart.pdf", height=5, width=14)
	datCol1<-dataPOS_COL1%>%mutate(TEXT=paste0(dataPOS_COL1$MutantAllele, "; ",dataPOS_COL1$MutantAA))
	ggplot(datCol1, 
			aes(x=xpos, y=height)) +
		geom_segment(data=datCol1,  aes(x=xpos, xend=xpos, y=height*0.5, yend=height, colour=MutationType)) +
		geom_segment(data=datCol1,  aes(x=xpos, xend=newPos, y=height, yend=height2, colour=MutationType)) +
		geom_rect(data=COL1A1exons, mapping=aes(xmin=V2, xmax=V3, ymin=V4, ymax=V5), color="black", alpha=0.5) +
		geom_point(data=datCol1, aes(x=newPos,y=height2,color=MutationType), size=3) +
		coord_cartesian(clip = 'off') +
		geom_text(data=datCol1, aes(x=newPos, y=height2, label=TEXT),
				hjust = 0, vjust = 0.5, angle=90) +
		theme(plot.margin = margin(6, 2, 2, 2, "cm")) + xlim(c(0,4500)) +
		ggtitle("COL1A1")+
		scale_colour_manual( values=KOLvec, breaks=names(KOLvec))

	datCol2<-dataPOS_COL2%>%mutate(TEXT=paste0(dataPOS_COL2$MutantAllele, "; ",dataPOS_COL2$MutantAA))
	ggplot(datCol2, 
			aes(x=xpos, y=height)) +
		geom_segment(data=datCol2, aes(x=xpos, xend=xpos, y=height*0.5, yend=height, colour=MutationType)) +
		geom_segment(data=datCol2, aes(x=xpos, xend=newPos, y=height, yend=height2, colour=MutationType)) +
		geom_rect(data=COL1A2exons, mapping=aes(xmin=V2, xmax=V3, ymin=V4, ymax=V5), color="black", alpha=0.5) +
		geom_point(data=datCol2, aes(x=newPos,y=height2,color=MutationType), size=3) +
		coord_cartesian(clip = 'off') +
		geom_text(data=datCol2,aes(x=newPos, y=height2, label=TEXT),
				hjust = 0, vjust = 0.5, angle=90) +
		theme(plot.margin = margin(6, 2, 2, 2, "cm")) + xlim(c(0,4000)) +
		ggtitle("COL1A2")+
		scale_colour_manual( values=KOLvec, breaks=names(KOLvec))


dev.off()




