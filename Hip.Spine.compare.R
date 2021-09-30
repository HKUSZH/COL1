library(dplyr)
library(ggplot2)

load("Hip_BMD2.RData")
load("Spine_BMD2.RData")


table(Spine_BMD2$SCANID_2==Spine_BMD2$SCANID_67)
table(Spine_BMD2$SCANID_2%in%Hip_BMD2$SCANID_67)

table(Spine_BMD2$PATIENT_KEY_34%in%Hip_BMD2$PATIENT_KEY_25)
table(Hip_BMD2$PATIENT_KEY_25%in%Spine_BMD2$PATIENT_KEY_34)


SpineKeyDate<-paste0(Spine_BMD2$PATIENT_KEY_34, "__", as.character(as.Date(Spine_BMD2$SCAN_DATE,"%d/%m/%Y")))
HipKeyDate<-paste0(Hip_BMD2$PATIENT_KEY_25, "__", as.character(as.Date(Hip_BMD2$SCAN_DATE,"%d/%m/%Y")))

table(SpineKeyDate %in%HipKeyDate)
table(HipKeyDate %in% SpineKeyDate)


Spine_BMD3<-Spine_BMD2%>%mutate(ScanDate=as.character(as.Date(SCAN_DATE,"%d/%m/%Y"))) %>%
	mutate(IDScanDate=paste0(PATIENT_KEY_34,"__",ScanDate))%>%
	group_by(IDScanDate) %>% summarize(num = n(), muBMDspine = mean(TOT_BMD, na.rm=TRUE), Age=mean(AGE))

ggplot(Spine_BMD3, aes(Age, muBMD ) )+
		geom_point(size=2) 

Hip_BMD3<-Hip_BMD2%>%mutate(ScanDate=as.character(as.Date(SCAN_DATE,"%d/%m/%Y"))) %>%
	mutate(IDScanDate=paste0(PATIENT_KEY_25,"__",ScanDate))%>%
	group_by(IDScanDate) %>% summarize(num = n(), muBMDhip = mean(HTOT_BMD, na.rm=TRUE), Age=mean(AGE))


table(Spine_BMD3$IDScanDate %in%Hip_BMD3$IDScanDate)
commonIDScanDate<-intersect(Spine_BMD3$IDScanDate, Hip_BMD3$IDScanDate)

str(unique(Spine_BMD3$IDScanDate))
str((Spine_BMD3$IDScanDate))
str(unique(Hip_BMD3$IDScanDate))
str((Hip_BMD3$IDScanDate))

Spine_BMD4<-Spine_BMD3[match(commonIDScanDate, Spine_BMD3$IDScanDate),]
Hip_BMD4<-Hip_BMD3[match(commonIDScanDate, Hip_BMD3$IDScanDate),]

Spine_BMD5<-Spine_BMD2[match(commonIDScanDate, SpineKeyDate),]
Spine_BMD6<-bind_cols(Spine_BMD4,muBMDhip =Hip_BMD4$muBMDhip ,Spine_BMD5)

plot(Spine_BMD4$muBMDspine , Hip_BMD4$muBMDhip)

cor(Spine_BMD4$muBMDspine , Hip_BMD4$muBMDhip, method="pearson")
cor(Spine_BMD4$muBMDspine , Hip_BMD4$muBMDhip, method="spearman")

summary(lm(Hip_BMD4$muBMDhip~ Spine_BMD4$muBMDspine ))
###########################################
pdf("Spine.vs.Hip.corr.pdf",width=9)
	Spine_BMD7<-Spine_BMD6%>%filter(OIcat2!="otherOI")
	summary(lm(muBMDhip~muBMDspine, data=Spine_BMD7))
	cor(Spine_BMD7$muBMDhip, Spine_BMD7$muBMDspine, , method="pearson")
	cor(Spine_BMD7$muBMDspine , Spine_BMD7$muBMDhip, method="spearman")
	summary(Spine_BMD7$muBMDspine )
	sd(Spine_BMD7$muBMDspine )
	summary(Spine_BMD7$muBMDhip)
	sd(Spine_BMD7$muBMDhip)

	ggplot(Spine_BMD7%>%filter(OIcat2!="otherOI"), aes(muBMDspine , muBMDhip) )+
			geom_point(aes(color=OIcat2),size=4) +
		geom_smooth(method="lm", se=T) 


	ggplot(Spine_BMD7%>%group_by(PATIENT_KEY_34)%>%filter(OIcat2!="otherOI"), aes(muBMDspine , muBMDhip, group=PATIENT_KEY_34 ) )+
			geom_point(aes(shape=SEX, color=OIcat2),size=5) +
		geom_line(aes(shape=SEX),size=1) +
		geom_smooth(data=Spine_BMD6%>%group_by(PATIENT_KEY_34), aes(muBMDspine , muBMDhip), se=T) 

#################
	Spine_BMD7<-Spine_BMD6%>%filter(OIcat2!="otherOI"&OIcat2!="ctrl")
	summary(lm(muBMDhip~muBMDspine, data=Spine_BMD7))
	cor(Spine_BMD7$muBMDhip, Spine_BMD7$muBMDspine, , method="pearson")
	cor(Spine_BMD7$muBMDspine , Spine_BMD7$muBMDhip, method="spearman")
	summary(Spine_BMD7$muBMDspine )
	sd(Spine_BMD7$muBMDspine )
	summary(Spine_BMD7$muBMDhip)
	sd(Spine_BMD7$muBMDhip)

	ggplot(Spine_BMD7, aes(muBMDspine , muBMDhip) )+
			geom_point(aes(color=OIcat2),size=4) +
		geom_smooth(method="lm", se=T) 


	ggplot(Spine_BMD7, 
		aes(muBMDspine , muBMDhip, group=PATIENT_KEY_34 ) )+
			geom_point(aes(shape=SEX, color=OIcat2),size=5) +
		geom_line(aes(shape=SEX),size=1) +
		geom_smooth(data=Spine_BMD6%>%group_by(PATIENT_KEY_34), aes(muBMDspine , muBMDhip), se=T) 
#################
	Spine_BMD7<-Spine_BMD6%>%filter(OIcat2=="ctrl")
	summary(lm(muBMDhip~muBMDspine, data=Spine_BMD7))
	cor(Spine_BMD7$muBMDhip, Spine_BMD7$muBMDspine, , method="pearson")
	cor(Spine_BMD7$muBMDspine , Spine_BMD7$muBMDhip, method="spearman")

	ggplot(Spine_BMD7, aes(muBMDspine , muBMDhip) )+
			geom_point(aes(color=OIcat2),size=4) +
		geom_smooth(method="lm", se=T) 


	ggplot(Spine_BMD7, 
		aes(muBMDspine , muBMDhip, group=PATIENT_KEY_34 ) )+
			geom_point(aes(shape=SEX, color=OIcat2),size=5) +
		geom_line(aes(shape=SEX),size=1) +
		geom_smooth(data=Spine_BMD6%>%group_by(PATIENT_KEY_34), aes(muBMDspine , muBMDhip), se=T) 


dev.off()



