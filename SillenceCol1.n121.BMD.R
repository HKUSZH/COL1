library(dplyr)
library(ggplot2)

load("BMD.for.COL1.n121.RData")
load("PATIENTS.RData")
########################
BMDName[match(SillenceCOL$HospitalIDs, PATIENTS$IDENTIFIER1)]==toupper(SillenceCOL$EngName)

PATIENTS_COL[["DOB"]]<-SillenceCOL$DOB[match(PATIENTS_COL$IDENTIFIER1, SillenceCOL$HospitalIDs)]

plot((Spine_COL%>%select(L1_BMD,L2_BMD,L3_BMD,L4_BMD,L5_BMD,TOT_BMD)))
table(Spine_COL$SCANID%in%ScanAnalysis_COL$SCANID)

spineBMD<-as.data.frame(Spine_COL2%>%select(L1_BMD,L2_BMD,L3_BMD,L4_BMD,TOT_BMD))
Spine_COL[which(rowSums(spineBMD==0)>0),]
########################
PATIENTS_Spine<-PATIENTS_COL[match(Spine_COL$PATIENT_KEY, PATIENTS_COL$PATIENT_KEY),]
Sillence_Spine<-SillenceCOL[match(PATIENTS_Spine$IDENTIFIER1, SillenceCOL$HospitalIDs),]

Scan_Spine<-ScanAnalysis_COL[match(Spine_COL$SCANID, ScanAnalysis_COL$SCANID),]

DOB<-as.Date(PATIENTS_COL$DOB[match(Spine_COL$PATIENT_KEY, PATIENTS_COL$PATIENT_KEY)],"%d/%m/%Y")

THISDay<-as.Date(Scan_Spine$SCAN_DATE,"%d/%m/%Y")
spineAGE<-as.numeric(difftime(THISDay, DOB, unit="weeks"))/52.25



Spine_COL2<-bind_cols(Spine_COL,Sillence_Spine,PATIENTS_Spine, Scan_Spine)
Spine_COL2[["DateBirth"]]<-DOB
Spine_COL2[["spineAGE"]]<-spineAGE
Spine_COL2[["THISDay"]]<-THISDay
colnames(Spine_COL2)<-gsub("\\.\\.\\.","_",colnames(Spine_COL2))

########################
load("../../REDCAP/REDCAP.2021.sep.24.RData")
BPTREAT<-sapply(strsplit(Sillence_Spine$RedCap_ID,", "),function(ID){
	his1<-nurseBase$BPtreat[match(ID, nurseBase$record_id)]
	x1<-setdiff(unlist(strsplit(his1, " ") ),"")
	x2<-strsplit(gsub("=$","",x1), "=")
	mat1<-sapply(x2,function(x){
		if(length(grep("[a-zA-Z]", x[1]))>0)
			x[1]<-as.character(as.Date(x[1], "%d-%b-%Y") )
		c(ID[1], ID[2], x[1:2])
		})
	t(mat1)
})

MAT<-data.frame(1,Sillence_Spine$HospitalIDs[1], BPTREAT[[1]])
colnames(MAT)<-1:6
for(i in 2:length(BPTREAT)){
	mati<-data.frame(i,Sillence_Spine$HospitalIDs[i], BPTREAT[[i]])
	colnames(mati)<-1:6
	MAT<-rbind(MAT,mati)
}
table(sapply(strsplit(MAT[,5],"-"),length))

########################

names(MAT)<-c("spineID","patID", "red1","red2","BPdate","BPdrug")
Spine_COL3<-bind_cols(spineINDEX=seq(nrow(Spine_COL2)), Spine_COL2)%>%filter(L1_BMD>0&L2_BMD>0&L3_BMD>0&L4_BMD>0)%>%group_by(THISDay) %>% slice(which.max(TOT_BMD))
MAT3<-MAT[MAT$spineID%in% Spine_COL3$spineINDEX, ]
MAT3[["BPage"]]<-as.numeric(difftime(as.Date(MAT3$BPdate, "%Y-%m-%d"), as.Date(PATIENTS_Spine$DOB[MAT3$spineID],"%d/%m/%Y"), unit="weeks"))/52.25
indCOL3<-match(MAT3$spineID, Spine_COL3$spineINDEX)

allDat<-bind_cols(MATid=seq(nrow(MAT3)),as.tbl(MAT3), Spine_COL3[indCOL3,])

IND<-data.frame("",matrix(rep(0,11),ncol=11))
for(thisID in unique(allDat$IDENTIFIER1)){
	indIDS<-which(allDat$IDENTIFIER1==thisID)
	thisDat<-allDat[indIDS, ]%>%arrange(spineAGE)
	#table(thisDat$THISDay, thisDat$TOT_BMD)
	if(length(thisDat$BPage)>0)
		for(j in 1:length(thisDat$BPage)){
			indMin<-max(which(thisDat$spineAGE<=thisDat$BPage[j]))
			indMax<-min(which(thisDat$spineAGE>thisDat$BPage[j]))
			indMin2<-indMin; indMax2<-indMax
			if(indMin== -Inf)indMin2<-1
			if(indMax== Inf)indMax2<-length(thisDat$BPage)
			NUMi<-length(unique(thisDat$spineAGE[seq(indMin2, indMax2)]))
			adjustedAge<-thisDat$BPage[j]
			AdjTYPE<-"notAdj"
			if(is.na(thisDat$BPage[j])){
				1+1
				#next
			}else if(thisDat$BPage[j]<min(thisDat$spineAGE)){
				thisBMD<-thisDat$TOT_BMD[indMin2]
				adjustedAge<-thisDat$spineAGE[indMin2]
				AdjTYPE<-"Adj"
			}else if(thisDat$BPage[j]> max(thisDat$spineAGE)){
					thisBMD<-NA
			}else if(indMin2==indMax2){
				thisBMD<-thisDat$TOT_BMD[indMin2]
			}else if(!is.na(thisDat$BPage[j])){
				y<-thisDat$TOT_BMD[c(indMin2, indMax2)]
				x<-thisDat$spineAGE[c(indMin2, indMax2)]
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
Spine_COL3<-Spine_COL2%>%filter(L1_BMD>0&L2_BMD>0&L3_BMD>0&L4_BMD>0)%>%group_by(THISDay) %>% slice(which.max(TOT_BMD))
	summarise(n = n())

ggplot(Spine_COL2, aes(spineAGE, TOT_BMD, label=ChiNameOnReport, fill=MutationType)) +
	geom_point(aes(color=MutationType),size=4) +
	geom_text() + ylim(c(0, 1.0))

ggplot(Spine_COL2%>%group_by(THISDay) %>% slice(which.max(TOT_BMD)), 
		aes(spineAGE, TOT_BMD, label=ChiNameOnReport, fill=MutationType)) +
	geom_point(aes(color=MutationType),size=4) +
	geom_text() + ylim(c(0, 1.0))


ggplot(Spine_COL2, aes(spineAGE, TOT_BMD, 
		group=ChiNameOnReport, label=ChiNameOnReport)) +
	geom_line(aes(color=MutationType),size=2) + 
	geom_point(aes(shape=SillenceType),size=4, color="darkgrey") +
	geom_text(data=Spine_COL2%>%group_by(ChiNameOnReport)%>%arrange(desc(TOT_BMD))%>%slice_head(n=1),aes(spineAGE,TOT_BMD,label=ChiNameOnReport)) + 
	ylim(c(0, 1.0))


ggplot(allDat2, aes(spineAGE, TOT_BMD, 
		group=ChiNameOnReport, label=ChiNameOnReport)) +
	geom_line(aes(color=MutationType),size=2) + 
	geom_point(aes(shape=SillenceType),size=4, color="darkgrey") +
	geom_text(data=Spine_COL2%>%group_by(ChiNameOnReport)%>%arrange(desc(TOT_BMD))%>%slice_head(n=1),aes(spineAGE,TOT_BMD,label=ChiNameOnReport)) + 
	ylim(c(0, 1.0))+
	geom_point(data=allDat2%>%filter(!is.na(thisBMD)&!is.na(BPage)), 
		aes(adjustedAge, thisBMD, shape=AdjTYPE), size=4, color="black") 




