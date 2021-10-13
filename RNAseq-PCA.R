	indBat2Col<-which(SAMPLES43$Batch==2 & (grepl("COL|pending",SAMPLES43$Mutation)|SAMPLES43$Index%in%c(14,15, 29, 30))&!SAMPLES43$Index%in%c(20) &maplibsize>1e7)
	PCABat2Col<-prcomp(t(logDatExpr[indNZ3,indBat2Col]))
	varianceBat2Col<-paste0(colnames(PCABat2Col$x)," (",percent(PCABat2Col$sdev^2/sum(PCABat2Col$sdev^2),0.1),")")

	ggplot(data.frame(PCABat2Col$x,SAMPLES43[indBat2Col,], Sillence43[indBat2Col,])%>%
			mutate(INFO=paste0(ChiName,"\n",sampleID,", ",Age, "\n",Mutation,"\n",MutantAllele)),
		aes(PC1,PC2,label=INFO))+
		geom_point(aes(color=Bone),size=4) +
		geom_density_2d() +
		geom_text() + xlab(varianceBat2Col[1]) + ylab(varianceBat2Col[2]) +
		ggtitle("PCA of second batch only (n=9)")
