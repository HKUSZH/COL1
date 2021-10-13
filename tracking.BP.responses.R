library(dplyr)
library(ggplot2)
library(ggbeeswarm)

load("Hip_BMD2.RData")
load("Spine_BMD2.RData")
##########################################
Sillence<-readr::read_tsv("../../Sillence_subtyping/Summary.of.Sillence.and.Geno.Sep29.2021.tsv")
#SillenceCOL<-Sillence%>%filter(grepl("COL", AffectGenes))

SillenceCOL<-Sillence%>%filter(AffectGenes%in%c("COL1A1", "COL1A2", "COL1A1, COL1A2"))

##########################################


Spine_BMD2[["ScanDay"]]<-as.Date(Spine_BMD2$SCAN_DATE,"%d/%m/%Y")

Spine_ctrl<-Spine_BMD2%>%filter(OIcat2=="ctrl")
Spine_ctrl[["LABEL"]]<-"ctrl"


Spine_ctrlM<-Spine_BMD2%>%filter(OIcat2=="ctrl")%>%filter(SEX=="M")
Spine_ctrlF<-Spine_BMD2%>%filter(OIcat2=="ctrl")%>%filter(SEX=="F")

fitCtrl_bmd_M<-loess(Spine_ctrlM$TOT_BMD~Spine_ctrlM$ScanAGE)
fitCtrl_bmd_F<-loess(Spine_ctrlF$TOT_BMD~Spine_ctrlF$ScanAGE)

fitCtrl_height_M<-loess(Spine_ctrlM$HEIGHT_76~Spine_ctrlM$ScanAGE)
fitCtrl_height_F<-loess(Spine_ctrlF$HEIGHT_76~Spine_ctrlF$ScanAGE)

predF<-predict(fitCtrl_F, seq(min(Spine_ctrlF$ScanAGE), max(Spine_ctrlF$ScanAGE), len=1e4))
predF<-predict(fitCtrl_F, Spine_ctrl$ScanAGE)

ggplot(bind_cols(Spine_ctrl, predF), aes(ScanAGE, TOT_BMD, group=SEX)) +
	geom_point(aes(color=SEX)) +
	geom_point(aes(ScanAGE,predF,color=SEX), size=4) +
	geom_smooth(aes(color=SEX))
##########################################
load("../../REDCAP/REDCAP.2021.sep.24.RData")

Spine_OIseq<-Spine_BMD2%>%filter(OIcat2!="otherOI"&OIcat2!="ctrl")
Sillence_Spine<-Sillence[match(Spine_OIseq$IDENTIFIER1, Sillence$HospitalIDs),]

for(i in 1:nrow(Sillence_Spine)){
	ID<-strsplit(Sillence_Spine$RedCap_ID[i],", ")[[1]]
	his1<-nurseBase$BPtreat[match(ID, nurseBase$record_id)]
	x1<-setdiff(unlist(strsplit(his1, " ") ),"")
	x2<-strsplit(gsub("=$","",x1), "=")
	mat1<-sapply(x2,function(x){
		if(length(grep("[a-zA-Z]", x[1]))>0)
			x[1]<-as.character(as.Date(x[1], "%d-%b-%Y") )
		c(i,Sillence_Spine$HospitalIDs[i], Sillence_Spine$DOB[i], ID[1], ID[2], x[1:2])
		})
	if(i==1){
		BPTREAT<-t(mat1)
	}else{
		BPTREAT<-rbind(BPTREAT, t(mat1))
	}
}
colnames(BPTREAT)<-c("spineID","patID", "DOB", "red1","red2","BPdate","BPdrug")
BPTREAT<-as.tbl(as.data.frame(BPTREAT))

table(sapply(strsplit(BPTREAT$BPdate,"-"),length))

BPTREAT[sapply(strsplit(BPTREAT$BPdate,"-"),length)==1, ]
BPTREAT[sapply(strsplit(BPTREAT$BPdate,"-"),length)==2, ]
ind2<-which(sapply(strsplit(BPTREAT$BPdate,"-"),length)==2)
BPTREAT$BPdate[ind2]<-paste0(gsub("-$","",BPTREAT$BPdate[ind2]),"-1")

BPTREAT$spineID<-as.integer(BPTREAT$spineID)
##########################################

Spine_COL<-bind_cols(spineINDEX=seq(nrow(Spine_OIseq)), Spine_OIseq, Sillence_Spine)%>%arrange(IDENTIFIER1, ScanAGE) 
Spine_COL[["LABEL"]]<-paste0("S",gsub(", .*","",Spine_COL$RedCap_ID))
Spine_COL[["isMis"]]<-"nonMis"
Spine_COL[["isMis"]][Spine_COL$MutationType=="missense"]<-"missense"


BPTREAT[["BPage"]]<-as.numeric(difftime(as.Date(BPTREAT$BPdate, "%Y-%m-%d"), as.Date(Spine_COL$BIRTHDATE[BPTREAT$spineID],"%d/%m/%Y"), unit="weeks"))/52.25
BPTREAT[["BPage2"]]<-as.numeric(difftime(as.Date(BPTREAT$BPdate, "%Y-%m-%d"), as.Date(BPTREAT$DOB,"%d/%m/%Y"), unit="weeks"))/52.25

plot(BPTREAT$BPage, BPTREAT$BPage2)

Spine_COL[["BPtreatedSixMoAgo"]]<-"no"
Spine_COL[["BPtreatedOneYrAgo"]]<-"no"
Spine_COL[["BPtreatedTwoYrAgo"]]<-"no"
Spine_COL[["BPtreated"]]<-"no"


for(i in 1:nrow(Spine_COL)){
	iBPTREAT<-BPTREAT[which(BPTREAT$patID%in%Spine_COL$IDENTIFIER1[i]&!is.na(BPTREAT$BPdrug)),,drop=F]
	table(iBPTREAT$BPage>=(Spine_COL$ScanAGE[i]-0.5), iBPTREAT$BPage <= (Spine_COL$ScanAGE[i]-1/12))
	table(iBPTREAT$BPage>=(Spine_COL$ScanAGE[i]-1), iBPTREAT$BPage <= (Spine_COL$ScanAGE[i]-1/12))
	indTreated0.5<-which(iBPTREAT$BPage>=(Spine_COL$ScanAGE[i]-0.5)& iBPTREAT$BPage <= (Spine_COL$ScanAGE[i]-1/12))
	indTreated1<-which(iBPTREAT$BPage>=(Spine_COL$ScanAGE[i]-1)& iBPTREAT$BPage <= (Spine_COL$ScanAGE[i]-1/12))
	indTreated2<-which(iBPTREAT$BPage>=(Spine_COL$ScanAGE[i]-2)& iBPTREAT$BPage <= (Spine_COL$ScanAGE[i]-1/12))
	if(length(indTreated0.5)>0)
		Spine_COL[["BPtreatedSixMoAgo"]][i]<-"yes"
	if(length(indTreated1)>0)
		Spine_COL[["BPtreatedOneYrAgo"]][i]<-"yes"
	if(length(indTreated2)>0)
		Spine_COL[["BPtreatedTwoYrAgo"]][i]<-"yes"
	if(nrow(iBPTREAT)>0)
		Spine_COL[["BPtreated"]][i]<-"yes"
}

table(Spine_COL[["BPtreatedSixMoAgo"]])
table(Spine_COL[["BPtreatedOneYrAgo"]])
table(Spine_COL[["BPtreatedTwoYrAgo"]])
table(Spine_COL[["BPtreated"]])

table(Spine_COL[["BPtreatedSixMoAgo"]],	Spine_COL[["BPtreatedOneYrAgo"]])
table(Spine_COL[["BPtreatedSixMoAgo"]],	Spine_COL[["BPtreatedTwoYrAgo"]])
##########################################
ggplot(Spine_COL, aes(ScanAGE, TOT_BMD, group=IDENTIFIER1)) +
	geom_line(aes(color=BPtreatedTwoYrAgo), size=2) +
	geom_point(aes(color=BPtreatedTwoYrAgo), size=4) 

ggplot(Spine_COL, aes(ScanAGE, TOT_BMD, group=IDENTIFIER1)) +
	geom_line(aes(color=BPtreated), size=2) +
	geom_point(aes(color=BPtreated), size=4) 


##########################################

Spine_COL_neverTreated<-Spine_COL%>% filter(BPtreatedTwoYrAgo=="no")
fitNoBP<-loess(TOT_BMD ~ ScanAGE, data=Spine_COL_neverTreated)
fitNoBP_Male<-loess(TOT_BMD ~ ScanAGE, data=Spine_COL_neverTreated %>% filter(SEX=="M"))
fitNoBP_Female<-loess(TOT_BMD ~ ScanAGE, data=Spine_COL_neverTreated %>% filter(SEX=="F"))

atan(1)*180/pi

Spine_COL_BP6mo<-Spine_COL%>% filter(BPtreatedSixMoAgo=="yes")
Spine_COL_BP1yr<-Spine_COL%>% filter(BPtreatedOneYrAgo=="yes")
Spine_COL_BP2yr<-Spine_COL%>% filter(BPtreatedTwoYrAgo=="yes")

Spine_COL_NoBP6mo<-Spine_COL%>% filter(BPtreatedSixMoAgo=="no")
Spine_COL_NoBP1yr<-Spine_COL%>% filter(BPtreatedOneYrAgo=="no")
Spine_COL_NoBP2yr<-Spine_COL%>% filter(BPtreatedTwoYrAgo=="no")



source("step6.help.funcs.R")
Spine_COL_BMD<-GetAngleDiff(Spine_COL, fitCtrl_bmd_F,fitCtrl_bmd_M, "TOT_BMD")
Spine_COL_HEIGHT<-GetAngleDiff(Spine_COL, fitCtrl_height_F,fitCtrl_height_M, "HEIGHT_76")

ggplot(Spine_COL_BMD%>%mutate(SEX_BP=paste0(SEX,"__", BPtreatedTwoYrAgo,"__",OIcat2)), 
		aes(ScanAGE, AngleDiff, group=SEX_BP)) + 
	geom_point(aes(color=SEX_BP, shape=BPtreatedTwoYrAgo), size=4)


boxplot(AngleDiff~paste0(FiveYearAgeGroupSex,"@",OIcat), data=Spine_COL_BMD, las=2)

pdf("AngleDiff.Spine.pdf")
	Spine_COL_BMD6 <- Spine_COL_BMD%>%filter(BPtreatedSixMoAgo=="yes")
	Spine_COL_BMD12 <- Spine_COL_BMD%>%filter(BPtreatedOneYrAgo=="yes")
	Spine_COL_BMD24 <- Spine_COL_BMD%>%filter(BPtreatedTwoYrAgo=="yes")
	Spine_COL_BMDever <- Spine_COL_BMD%>%filter(BPtreated=="yes")
	Spine_COL_BMDnever <- Spine_COL_BMD%>%filter(BPtreated=="no")

	GET.WILCOX<-function(OBJ){
		sapply(split(OBJ$AngleDiff, OBJ$FiveYearAgeGroupSex),function(x){
			y<-x[!is.na(x)]
			if(length(y)>2)
				wilcox.test(y)$p.value
			else
				NA
		})
	}

	B6<-GET.WILCOX(Spine_COL_BMD6)
	B12<-GET.WILCOX(Spine_COL_BMD12)
	B24<-GET.WILCOX(Spine_COL_BMD24)
	Bever<-GET.WILCOX(Spine_COL_BMDever )
	Bnever<-GET.WILCOX(Spine_COL_BMDnever )

	tab6<-table(Spine_COL_BMD$OIcat2, Spine_COL_BMD$BPtreatedSixMoAgo)
	tab12<-table(Spine_COL_BMD$OIcat2, Spine_COL_BMD$BPtreatedOneYrAgo)
	tab24<-table(Spine_COL_BMD$OIcat2, Spine_COL_BMD$BPtreatedTwoYrAgo)
	tabEver<-table(Spine_COL_BMD$OIcat2, Spine_COL_BMD$BPtreated)

	boxplot(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_BMD6, las=2, outline=F)
	abline(h=0,lty=2, lwd=3, col=3)
	boxplot(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_BMD12, las=2, outline=F)
	abline(h=0,lty=2, lwd=3, col=3)
	boxplot(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_BMD24, las=2, outline=F)
	abline(h=0,lty=2, lwd=3, col=3)

	boxplot(list(Spine_COL_BMDnever$AngleDiff, Spine_COL_BMDever$AngleDiff))


##############
	bpp1<-boxplot(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_BMDever, las=2, outline=F, ylim=c(-80,70), col="#F3766E")
	text(seq_along(bpp1$names), bpp1$stats[4,]+2,paste0(signif(Bever[bpp1$names],3),""), xpd=T)
	abline(h=0,lty=2, lwd=3, col=3)

	bpp1<-boxplot(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_BMDnever , las=2, outline=F, ylim=c(-80,70))
	text(seq_along(bpp1$names), bpp1$stats[4,],signif(Bnever[bpp1$names],3), xpd=T)
	abline(h=0,lty=2, lwd=3, col=3)
##############
	summary(lm(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_BMD6))
	summary(lm(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_BMD12))
	summary(lm(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_BMD24))
	summary(lm(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_BMDever ))

	t.test((Spine_COL_BMD6)$AngleDiff)
	t.test((Spine_COL_BMD12)$AngleDiff)
	t.test((Spine_COL_BMD24)$AngleDiff)
	t.test((Spine_COL_BMDever)$AngleDiff)



	tabBP<-rbind(tab6, tab12, tab24, tabEver)[,c(2,1)]
	barplot(t(tabBP), las=2)

	H6<-GET.WILCOX(Spine_COL_HEIGHT6)
	H12<-GET.WILCOX(Spine_COL_HEIGHT12)
	H24<-GET.WILCOX(Spine_COL_HEIGHT24)
	Hever<-GET.WILCOX(Spine_COL_BMDever)
	Hnever<-GET.WILCOX(Spine_COL_BMDnever)

	Spine_COL_HEIGHT6<-Spine_COL_HEIGHT%>%filter(BPtreatedSixMoAgo=="yes")
	Spine_COL_HEIGHT12<-Spine_COL_HEIGHT%>%filter(BPtreatedOneYrAgo=="yes")
	Spine_COL_HEIGHT24<-Spine_COL_HEIGHT%>%filter(BPtreatedTwoYrAgo=="yes")
	Spine_COL_HEIGHTever<-Spine_COL_HEIGHT%>%filter(BPtreated=="yes")
	Spine_COL_HEIGHTnever<-Spine_COL_HEIGHT%>%filter(BPtreated=="no")

	boxplot(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_HEIGHT6, las=2, outline=F)
	abline(h=0,lty=2, lwd=3, col=3)
	boxplot(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_HEIGHT12, las=2, outline=F)
	abline(h=0,lty=2, lwd=3, col=3)
	boxplot(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_HEIGHT24, las=2, outline=F)
	abline(h=0,lty=2, lwd=3, col=3)

##############
##############
	boxplot(list(Spine_COL_HEIGHTnever$AngleDiff, Spine_COL_HEIGHTever$AngleDiff))

	bpp1<-boxplot(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_HEIGHTever, las=2, outline=F, col="#F3766E")
	text(seq_along(bpp1$names), bpp1$stats[4,]+2,paste0(signif(Hever[bpp1$names],3),""), xpd=T)
	abline(h=0,lty=2, lwd=3, col=3)

	bpp1<-boxplot(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_BMDnever , las=2, outline=F, ylim=c(-80,70))
	text(seq_along(bpp1$names), bpp1$stats[4,],signif(Hnever[bpp1$names],3), xpd=T)
	abline(h=0,lty=2, lwd=3, col=3)
##############
	summary(lm(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_HEIGHT6))
	summary(lm(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_HEIGHT12))
	summary(lm(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_HEIGHT24))

	t.test((Spine_COL_HEIGHT6)$AngleDiff)
	t.test((Spine_COL_HEIGHT12)$AngleDiff)
	t.test((Spine_COL_HEIGHT24)$AngleDiff)




dev.off()

tabBP<-table(Spine_COL_BMDneverever$IDENTIFIER1, Spine_COL_BMDneverever$BPtreated)

##############
##############
library(ggridges)

	WILCOX<-function(OBJ, VAR1, VAR2){
		var1<-unique(OBJ[[VAR1]])
		var2<-unique(OBJ[[VAR2]])

		MATres<-matrix(NA, 6, length(var1)*length(var2))
		print(dim(MATres))
		matnames<-apply(expand.grid(var2, var1),1, paste0, collapse="_")
		#print(var1)
		#print(var2)
		colnames(MATres)<-matnames
		rownames(MATres)<-seq(0,5)
		#print(head(MATres))
		for(k in 0:5){
			for(i in 1:length(var1)){
				for(j in 1:length(var2)){
					cat(var2[j],"_", var1[i], ";  ")
					ij<-(i-1)*length(var1)+j
					z<-OBJ$AngleDiff[which(OBJ[[VAR1]]==var1[i] & OBJ[[VAR2]]==var2[j] & OBJ$FiveYearAgeGroup==k)]
					if(length(z)>2)
						MATres[(k+1),ij]<-wilcox.test(z)$p.value
				}
			}
			cat("\n")
		}
		return(MATres)
	}

pdf("angle.diff.ridge.plot.pdf",height=6, width=10.8)
	Spine_COL_BMDneverever<-bind_rows(Spine_COL_BMDnever, Spine_COL_BMDever)%>%filter(!is.na(AngleDiff))
	Spine_COL_HEIGHTneverever<-bind_rows(Spine_COL_HEIGHTnever, Spine_COL_HEIGHTever)%>%filter(!is.na(AngleDiff))

	ggplot(Spine_COL_BMDneverever, 
		aes(x = AngleDiff , y = BPtreated, fill=BPtreated)) + 
		geom_density_ridges(color="white")+ scale_x_continuous(breaks =seq(-60,60,by=30), limits=c(-60,60))+
		facet_wrap(~SEX)
	ggplot(Spine_COL_HEIGHTneverever%>%mutate(newGrp=paste0(FiveYearAgeGroup,"-",isMis)), 
		aes(x = SEX, y = AngleDiff)) + 
		geom_violin()+geom_boxplot()

############

	ggplot(Spine_COL_BMDneverever%>%mutate(newGrp=paste0(FiveYearAgeGroup,"-",BPtreated)), 
		aes(x = AngleDiff , y = newGrp, fill=BPtreated)) + 
		geom_density_ridges(color="white")+ scale_x_continuous(breaks =seq(-60,60,by=30), limits=c(-60,60))+
		facet_wrap(~SEX)+xlab("AngleDiff Spine BMD")+theme(legend.position = "none") 

	WILCOX(Spine_COL_BMDneverever, "BPtreated", "SEX")
	WILCOX(Spine_COL_HEIGHTneverever, "BPtreated", "SEX")

	ggplot(Spine_COL_HEIGHTneverever%>%mutate(newGrp=paste0(FiveYearAgeGroup,"-",BPtreated)), 
		aes(x = AngleDiff , y = newGrp, fill=BPtreated)) + 
		geom_density_ridges(color="white")+ scale_x_continuous(breaks =seq(-60,60,by=30), limits=c(-60,60))+
		facet_wrap(~SEX)+xlab("AngleDiff Height")+theme(legend.position = "none") 
############
	ggplot(Spine_COL_BMDneverever%>%mutate(newGrp=paste0(FiveYearAgeGroup,"-",isMis)), 
		aes(x = AngleDiff , y = newGrp, fill=isMis )) + 
		geom_density_ridges(color="white")+ scale_x_continuous(breaks =seq(-60,60,by=30), limits=c(-60,60))+
		facet_wrap(~BPtreated)+xlab("AngleDiff Spine BMD")+theme(legend.position = "none") 

	ggplot(Spine_COL_HEIGHTneverever%>%mutate(newGrp=paste0(FiveYearAgeGroup,"-",isMis)), 
		aes(x = AngleDiff , y = newGrp, fill=isMis )) + 
		geom_density_ridges(color="white")+ scale_x_continuous(breaks =seq(-60,60,by=30), limits=c(-60,60))+
		facet_wrap(~BPtreated)+xlab("AngleDiff Height")+theme(legend.position = "none") 

	WILCOX(Spine_COL_BMDneverever, "BPtreated", "isMis")
	WILCOX(Spine_COL_HEIGHTneverever, "BPtreated", "isMis")

############
	ggplot(Spine_COL_BMDneverever%>%mutate(newGrp=paste0(FiveYearAgeGroup,"-",BPtreated)), 
		aes(x = AngleDiff , y = newGrp, fill=BPtreated)) + 
		geom_density_ridges(color="white")+ scale_x_continuous(breaks =seq(-60,60,by=30), limits=c(-60,60))+
		facet_wrap(~isMis)+xlab("AngleDiff Spine BMD")
	ggplot(Spine_COL_HEIGHTneverever%>%mutate(newGrp=paste0(FiveYearAgeGroup,"-",BPtreated)), 
		aes(x = AngleDiff , y = newGrp, fill=BPtreated)) + 
		geom_density_ridges(color="white")+ scale_x_continuous(breaks =seq(-60,60,by=30), limits=c(-60,60))+
		facet_wrap(~isMis)+xlab("AngleDiff Height")


############
	ggplot(Spine_COL_HEIGHTneverever%>%mutate(newGrp=paste0(FiveYearAgeGroup,"-",isMis)), 
		aes(x = isMis  , y = AngleDiff , group=isMis )) + 
		geom_violin()+geom_boxplot()


dev.off()

##############
##############

ggplot(Spine_COL_HEIGHT%>%filter(BPtreatedSixMoAgo=="yes"), 
		aes(FiveYearAgeGroupSex, AngleDiff)) + 
	geom_violin(width=3)+geom_quasirandom(size=4)
	
	geom_boxplot(aes(color=FiveYearAgeGroupSex), width=1)+
	geom_point(aes(color=SEX_BP, shape=BPtreatedTwoYrAgo), size=4)




boxplot(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_BMD, las=2, ylim=c(-50,50))

boxplot(AngleDiff~FiveYearAgeGroupSex, data=Spine_COL_BMD%>%filter(BPtreatedOneYrAgo=="no"), las=2, ylim=c(-50,50))


boxplot(SlopeDiff~FiveYearAgeGroupSex, data=Spine_COL>%filter(BPtreatedTwoYrAgo=="yes"), las=2)

ggplot(Spine_COL_BMD, aes(ScanAGE, AngleDiff, label=ChiNameOnReport))+
	geom_point() +
	geom_text()

ggplot(Spine_COL_BMD%>%mutate(SEX_BP=paste0(SEX,"__", BPtreatedTwoYrAgo,"__",OIcat2)), 
		aes(FiveYearAgeGroupSex, AngleDiff, group=SEX_BP)) + 
	geom_violin(width=1)+
	geom_boxplot(aes(color=FiveYearAgeGroupSex), width=1)+
	geom_point(aes(color=SEX_BP, shape=BPtreatedTwoYrAgo), size=4)



##########################################
indCOL3<-match(BPTREAT$spineID, Spine_COL$spineINDEX)

allDat<-bind_cols(MATid=seq(nrow(BPTREAT)),BPTREAT, Spine_COL[indCOL3,])%>%
	filter(L1_BMD>0&L2_BMD>0&L3_BMD>0&L4_BMD>0)%>%group_by(ScanDay) %>% slice(which.max(TOT_BMD))


IND<-data.frame("",matrix(rep(0,11),ncol=11))
for(thisID in unique(allDat$IDENTIFIER1)){
	indIDS<-which(allDat$IDENTIFIER1==thisID)
	thisDat<-allDat[indIDS, ]%>%arrange(ScanAGE)
	#table(thisDat$ScanDay, thisDat$TOT_BMD)
	if(length(thisDat$BPage)>0)
		for(j in 1:length(thisDat$BPage)){
			indMin<-max(which(thisDat$ScanAGE<=thisDat$BPage[j]))
			indMax<-min(which(thisDat$ScanAGE>thisDat$BPage[j]))
			indMin2<-indMin; indMax2<-indMax
			if(indMin== -Inf)indMin2<-1
			if(indMax== Inf)indMax2<-length(thisDat$BPage)
			NUMi<-length(unique(thisDat$ScanAGE[seq(indMin2, indMax2)]))
			adjustedAge<-thisDat$BPage[j]
			AdjTYPE<-"notAdj"
			if(is.na(thisDat$BPage[j])){
				1+1
				#next
			}else if(thisDat$BPage[j]<min(thisDat$ScanAGE)){
				thisBMD<-thisDat$TOT_BMD[indMin2]
				adjustedAge<-thisDat$ScanAGE[indMin2]
				AdjTYPE<-"Adj"
			}else if(thisDat$BPage[j]> max(thisDat$ScanAGE)){
					thisBMD<-NA
			}else if(indMin2==indMax2){
				thisBMD<-thisDat$TOT_BMD[indMin2]
			}else if(!is.na(thisDat$BPage[j])){
				y<-thisDat$TOT_BMD[c(indMin2, indMax2)]
				x<-thisDat$ScanAGE[c(indMin2, indMax2)]
				fit1<-lm(y~x)
				thisBMD<-predict(fit1,data.frame(x=thisDat$BPage[j]))
			}else{
				thisBMD<-NA
			}
			IND[nrow(IND)+1,1:2]<- c(thisID, AdjTYPE)
			IND[nrow(IND), -c(1:2)]<- c(indMin, indMax, indMin2, indMax2, NUMi, thisDat$MATid[j], 
				thisDat$BPage[j], adjustedAge, thisBMD, thisDat$TOT_BMD[j])
		}
}
IND<-IND[-1,]
colnames(IND)<-c("IDS", "AdjTYPE", "indMin", "indMax", "indMin2", "indMax2", "NUMi", "MATid", "BPage", "adjustedAge",  "thisBMD", "TOT_BMD")

allDat2<-bind_cols(as.tbl(data.frame(IND[match(allDat$MATid, IND$MATid), c("thisBMD", "AdjTYPE", "adjustedAge")])), allDat)

########################


pdf("BMD.and.Height.pdf",width=16)
	ggplot(Spine_COL, aes(ScanAGE, TOT_BMD, 
			 label=LABEL)) +
		ylim(c(0, 1.0))+xlim(c(0, 30))+
		geom_point(aes(color=BPtreated, shape=SEX),size=4) +
		geom_text(data=Spine_COL%>%group_by(ChiNameOnReport)%>%arrange(desc(TOT_BMD))%>%slice_head(n=1),aes(ScanAGE,TOT_BMD,label=LABEL)) +
		geom_line(aes(color=BPtreated, group=ChiNameOnReport),size=1) + 
		geom_smooth(data=Spine_COL, aes(ScanAGE,TOT_BMD)) +
		geom_smooth(data=Spine_ctrl, aes(ScanAGE,TOT_BMD)) +
		geom_point(data=Spine_ctrl, aes(ScanAGE,TOT_BMD),color="grey", size=4) +
		facet_wrap(~SEX)

	ggplot(Spine_COL, aes(ScanAGE, TOT_BMD, 
			group=ChiNameOnReport, label=ChiNameOnReport)) +
		geom_line(aes(color=BPtreatedTwoYrAgo),size=1) + 
		geom_point(aes(color=BPtreatedTwoYrAgo, shape=SEX),size=4) +
		geom_text(data=Spine_COL%>%group_by(ChiNameOnReport)%>%arrange(desc(TOT_BMD))%>%slice_head(n=1),aes(ScanAGE,TOT_BMD,label=LABEL)) + 
		ylim(c(0, 1.0))+
		facet_wrap(~SEX)

	ggplot(Spine_COL, aes(ScanAGE, TOT_BMD, 
			group=IDENTIFIER1)) +
		geom_line(aes(color=isMis),size=1) + 
		geom_point(aes(color=isMis, shape=SEX),size=4) +
		geom_text(data=Spine_COL%>%group_by(IDENTIFIER1)%>%arrange(desc(TOT_BMD))%>%slice_head(n=1),aes(ScanAGE,TOT_BMD,label=LABEL)) + 
		geom_smooth(data=Spine_ctrl, aes(ScanAGE,TOT_BMD)) +
		geom_point(data=Spine_ctrl, aes(ScanAGE,TOT_BMD),color="grey", size=4) +
		ylim(c(0, 1.0))+xlim(c(0, 30))+
		facet_wrap(~SEX)

###########
	ggplot(Spine_COL, aes(ScanAGE, HEIGHT_76, 
			 label=LABEL)) +
		ylim(c(60, 170))+	xlim(c(0, 30))+
		geom_point(aes(color=BPtreated, shape=SEX),size=4) +
		geom_text(data=Spine_COL%>%group_by(ChiNameOnReport)%>%arrange(desc(HEIGHT_76))%>%slice_head(n=1),aes(ScanAGE,HEIGHT_76,label=LABEL)) +
		geom_line(aes(color=BPtreated, group=ChiNameOnReport),size=1) + 
		geom_smooth(data=Spine_COL, aes(ScanAGE,HEIGHT_76)) +
		geom_smooth(data=Spine_ctrl, aes(ScanAGE,HEIGHT_76)) +
		facet_wrap(~SEX)


	ggplot(Spine_COL, aes(ScanAGE, HEIGHT_76, 
			group=ChiNameOnReport, label=ChiNameOnReport)) +
		geom_line(aes(color=BPtreatedTwoYrAgo),size=1) + 
		geom_point(aes(color=BPtreatedTwoYrAgo, shape=SEX),size=4) +
		geom_text(data=Spine_COL%>%group_by(ChiNameOnReport)%>%arrange(desc(HEIGHT_76))%>%slice_head(n=1),aes(ScanAGE,HEIGHT_76,label=LABEL)) + 
		facet_wrap(~SEX)

	ggplot(Spine_COL, aes(ScanAGE, HEIGHT_76, 
			group=ChiNameOnReport, label=ChiNameOnReport)) +
		geom_line(aes(color=isMis),size=1) + 
		geom_point(aes(color=isMis, shape=SEX),size=4) +
		geom_text(data=Spine_COL%>%group_by(ChiNameOnReport)%>%arrange(desc(HEIGHT_76))%>%slice_head(n=1),aes(ScanAGE,HEIGHT_76,label=LABEL)) + 
		facet_wrap(~SEX)

dev.off()


ggplot(allDat2, aes(ScanAGE, TOT_BMD, 
		group=ChiNameOnReport, label=ChiNameOnReport)) +
	geom_line(aes(color=BPtreatedTwoYrAgo),size=2) + 
	geom_point(aes(color=BPtreatedTwoYrAgo),size=6) +
	geom_text(data=allDat2%>%group_by(ChiNameOnReport)%>%arrange(desc(TOT_BMD))%>%slice_head(n=1),aes(ScanAGE,TOT_BMD,label=LABEL)) + 
	ylim(c(0, 1.0))+
	geom_point(data=allDat2%>%filter(!is.na(thisBMD)&!is.na(BPage)), 
		aes(adjustedAge, thisBMD, shape=AdjTYPE), size=4, color="black") + 
	facet_wrap(~SEX)

pdf("BP.treatment.response.gender.pdf", width=16)
	ggplot(allDat2, aes(ScanAGE, TOT_BMD, 
			group=IDENTIFIER1, label=IDENTIFIER1)) +
		geom_line(aes(color=BPtreated),size=1) + 
		geom_point(aes(color=BPtreated),size=4) +
		geom_text(data=allDat2%>%group_by(IDENTIFIER1)%>%arrange(desc(TOT_BMD))%>%slice_head(n=1),aes(ScanAGE,TOT_BMD,label=LABEL)) + 
		ylim(c(0, 1.0))+
		geom_point(data=allDat2%>%filter(!is.na(thisBMD)&!is.na(BPage)), 
			aes(adjustedAge, thisBMD), size=2, color="black", shape=4) + 
		geom_smooth(data=Spine_ctrl, aes(ScanAGE,TOT_BMD)) +
		geom_point(data=Spine_ctrl, aes(ScanAGE,TOT_BMD),color="grey", size=4) +
		facet_wrap(~SEX)+
		ylim(c(0, 1.0))+xlim(c(0, 30))+theme(legend.position = "none") 

	ggplot(Spine_ctrl, aes(ScanAGE, TOT_BMD)) +
		geom_smooth() +
		facet_wrap(~SEX)+
		ylim(c(0, 1.0))+xlim(c(0, 30))+theme(legend.position = "none") 
###########

	ggplot(allDat2, aes(ScanAGE, HEIGHT_76, 
			group=IDENTIFIER1, label=IDENTIFIER1)) +
		geom_line(aes(color=BPtreated),size=1) + 
		geom_point(aes(color=BPtreated),size=4) +
		geom_text(data=allDat2%>%group_by(IDENTIFIER1)%>%arrange(desc(HEIGHT_76))%>%slice_head(n=1),aes(ScanAGE,HEIGHT_76,label=LABEL)) + 
		ylim(c(0, 1.0))+
		geom_point(data=allDat2%>%filter(!is.na(thisBMD)&!is.na(BPage)), 
			aes(adjustedAge, thisBMD), size=2, color="black", shape=4) + 
		geom_smooth(data=Spine_ctrl, aes(ScanAGE,HEIGHT_76)) +
		geom_point(data=Spine_ctrl, aes(ScanAGE,HEIGHT_76),color="grey", size=4) +
		facet_wrap(~SEX)+
		ylim(c(60, 170))+xlim(c(0, 30))+theme(legend.position = "none") 

	ggplot(Spine_ctrl, aes(ScanAGE, HEIGHT_76)) +
		geom_smooth() +
		facet_wrap(~SEX)+
		ylim(c(60, 170))+xlim(c(0, 30))+theme(legend.position = "none") 
dev.off()


ggplot(allDat2, aes(ScanAGE, HEIGHT_76, 
		group=ChiNameOnReport, label=ChiNameOnReport)) +
	geom_line(aes(color=AffectGenes),size=2) + 
	geom_point(aes(shape=SillenceType),size=4, color="darkgrey") +
	geom_text(data=allDat2%>%group_by(ChiNameOnReport)%>%arrange(desc(HEIGHT_76))%>%slice_head(n=1),aes(ScanAGE,HEIGHT_76,label=ChiNameOnReport)) +
	geom_smooth(data=Spine_COL,  aes(ScanAGE, HEIGHT_76), color="#F8766D") 

	geom_smooth(method="loess",se=T,data=allDat2%>%filter(AffectGenes=="COL1A1"&!is.na(HEIGHT_94)&!is.na(ScanAGE)), aes(ScanAGE, HEIGHT_94),color="#F8766D") +
	geom_smooth(method="loess",se=T,data=allDat2%>%filter(AffectGenes=="COL1A2"&!is.na(HEIGHT_94)&!is.na(ScanAGE)), aes(ScanAGE, HEIGHT_94),color="#00BA38") 








