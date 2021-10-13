load("Spine_BMD2.RData")
##########################################
pdf("weight.vs.height.pdf", width=9)
	ggplot(Spine_BMD2, aes(HEIGHT_76, WEIGHT_75, group=OIcat2)) +
		geom_point(aes(color=OIcat2), size=4) +
		geom_smooth(aes(color=OIcat2)) 

	ggplot(Spine_BMD2, aes(HEIGHT_76, sqrt(WEIGHT_75), group=OIcat2)) +
		geom_point(aes(color=OIcat2), size=4)  +
		geom_smooth(method="lm", aes(color=OIcat2))

	ggplot(Spine_BMD2, aes(ScanAGE, BMI_126, group=OIcat2)) +
		geom_smooth(method="lm", aes(color=OIcat2)) +
		geom_point(aes(color=OIcat2), size=4) 
dev.off()
