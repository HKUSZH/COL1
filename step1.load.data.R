library(dplyr)
library("xlsx")

Sillence<-readr::read_tsv("../../Sillence_subtyping/Summary.of.Sillence.and.Geno.Sep27c.2021.tsv")
#SillenceCOL<-Sillence%>%filter(grepl("COL", AffectGenes))

SillenceCOL<-Sillence%>%filter(AffectGenes%in%c("COL1A1", "COL1A2", "COL1A1, COL1A2"))

#write.xlsx(SillenceCOL[is.na(SillenceCOL$DOB),], "Sillence.OI.COL1.noDOB.xlsx", sheetName = "noDOB", 
#  col.names = TRUE, row.names = TRUE, append = T)

gender<-nurseBase$gender[match(Sillence$ChiNameOnReport, nurseBase$pa_name1)]
ind.nogender<-which(is.na(gender))
gender[ind.nogender]<-c(1, 1,2,2,2,2,2,1,1,1,1)
############################################

PATIENTS<-readr::read_csv("PATIENT.txt")
Hip<-readr::read_csv("Hip.txt")
Forearm<-readr::read_csv("Forearm.txt")
HipHSA<-readr::read_csv("HipHSA.txt")
ObesityIndices<-readr::read_csv("ObesityIndices.txt")
ScanAnalysis<-readr::read_csv("ScanAnalysis.txt")
Spine<-readr::read_csv("Spine.txt")
SubRegionBone<-readr::read_csv("SubRegionBone.txt")
SubRegionComposition<-readr::read_csv("SubRegionComposition.txt")
TenYearFxRisk<-readr::read_csv("TenYearFxRisk.txt")
Wbody<-readr::read_csv("Wbody.txt")
WbodyComposition<-readr::read_csv("WbodyComposition.txt")


#####################################
DOB<-as.Date(PATIENTS$BIRTHDATE,"%d/%m/%Y")
THISDay<-as.Date(PATIENTS$LAST_UPDATE,"%d/%m/%Y")

AGE<-as.numeric(difftime(THISDay, DOB, unit="weeks"))/52.25
PATIENTS[["AGE"]]<-AGE

PATIENTS2<-PATIENTS%>%filter(AGE>0&AGE<100)
#save(PATIENTS,PATIENTS2,file="PATIENTS.RData")

############################################
library(gplots)

BMDName<-paste0(PATIENTS$LAST_NAME,PATIENTS$FIRST_NAME)
BMDName[!grepl("[a-zA-Z]+",BMDName)]
indEng<-which(grepl("[a-zA-Z]+",BMDName)&!is.na(PATIENTS$LAST_NAME)&!is.na(PATIENTS$FIRST_NAME))
BMDName[indEng]<-toupper(BMDName[indEng])

noDupName<-c("caiyuhang", "caozhejia", "cengchengyu", "chenguoguo", "chenpeifen", "chentangxuan", "chenweicheng", "fangzhenxuan", "hejinzeng", "hexiaojin", "hexinyan", "huangsiyu", "jiangmingyue", "jikaide", "liaoxihao", "liaozeguo", "liuhuan", "liuyuping", "pengleying", "suyuxin", "sunruiyu", "tanjingwen", "tangziyan", "tangyibo", "tianzupeng", "wujiebb", "yanzixi", "zhanghanwen", "zhouziqian", "zhuyinqi", "zoumengci", "zhangzhenyi", "xuzhangzhi", "zhuyongchen", "wangchaojia", "shiboxin", "linshiqi", "Ko Pak Sam", "dengzixi", "wangpeiyu", "lijiaqin", "wangwenjun", "dailimin", "yanhongli", "zuziqi", "wusonghao", "yingxianglei", "lixinyan", "guoziran", "zhuxun", "guoyao", "liulinhao", "chenhanhong", "zhangjihaoen", "tongyuhan", "huangjinyu", "yuheng", "weijiahui", "zhangshibo", "liyun", "chenziyang", "chenyouliang", "pengjinyao", "zhangwenyu", "zhangjinshuo", "heyan", "yinhongjuan", "fengyucheng", "tianxiaoxia", "lujia", "wangtao", "huangjixin", "huangyukun", "fujunhao", "xuningze", "yangjiangyuan", "linpeiye", "linzeqin", "lilinwen", "linchangyue", "liumanting", "liumanyan", "liujun", "l???§?1boan", "hongchengming", "zhouyimao", "luoyicheng", "qiuchongxuan", "xudechen", "l???§?1sitong", "huangsijie", "caoyanan", "xieguan", "liujiajian", "tangzihan", "chenqiyue", "liaoyongxin", "wubin", "zhangyangyi", "liyuewen", "yedongchen", "hefeifan", "wushurong", "liangzhixia", "renhaolin", "yangxueying", "wangruojun", "liuxin", "gujunhao", "zhoutianci", "zhangshengtao", "liangyankuan", "wangjian", "wangzihan", "wangxiaoen", "qinyilin", "luoyuchun", "caijinhong", "zoujinxiu", "weijitong", "longyawen", "longyashi", "dengyinyin", "penglei", "liuhao", "wangfanzhuo", "zhangzhijun", "qingzhenyu", "zhangwenhao", "yangxingyue", "mamaileyan", "lizhicheng", "duchengzhi", "gaojiawei", "lintianxing", "wangxiangyu", "liqinghui", "lichenghao", "huangzheng", "chenyitao", "yuanbochao", "zhengziyi", "luositong", "zhanghaosen", "licaiwei", "zouming", "renyuxi", "wangwanqi", "zhentingting", "yangyuxi", "qiuzixiang", "gaozhengqin", "yangjiarui", "zhupengyu", "linjunming", "youxinyu", "liling", "gantingting", "suncong", "qianrongxuan", "zhouxinyi", "zhouxinyue", "zhangpei", "lixianhai", "duzhangyi", "moyuyan", "lipengfei", "tianyuxuan", "zhuyangmei", "weifugang", "yudacheng", "caomingxuan", "maoyunsheng", "fanziyan", "fanzixuan", "chenxipeng", "chenhuangwenyan", "mahaidong", "zhangjingen", "zhouhuifang", "gechenyu", "yecanyang", "meimeirenyuhan", "wuwenyan", "zhangyue", "yuyingying", "dingjunhui", "zhaoyi")
noDupNameUPPER<-toupper(noDupName)

venn(list(noDupNameUPPER,BMDName))
venn(list(toupper(SillenceCOL$EngName),BMDName))
setdiff(toupper(SillenceCOL$EngName), BMDName)

venn(list(SillenceCOL$HospitalIDs,PATIENTS$IDENTIFIER1))
nonInBMD<-setdiff(SillenceCOL$HospitalIDs, PATIENTS$IDENTIFIER1)


sapply(noDupNameUPPER,function(x){
	paste0(unique(PATIENTS$PATIENT_KEY[BMDName%in%x]),collapse=", ")
})
############################################
SillenceCOL[SillenceCOL$HospitalIDs%in%nonInBMD, ]
indnotin<-which(SillenceCOL$HospitalIDs%in%nonInBMD)

 
PATIENTS_COLnotin<-PATIENTS[BMDName%in%toupper(SillenceCOL$EngName[indnotin]),]


PATIENTS_COL<-PATIENTS[which(PATIENTS$IDENTIFIER1%in%setdiff(SillenceCOL$HospitalIDs,NA)),]
BMDName[match(SillenceCOL$HospitalIDs, PATIENTS$IDENTIFIER1)]

table(Spine$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY)
table(Hip$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY)
table(Forearm$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY)
table(Wbody$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY)
table(WbodyComposition$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY)
table(ObesityIndices$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY)

Spine_COL<-Spine[Spine$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY, ]
Hip_COL<-Hip[Hip$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY, ]
Forearm_COL<-Forearm[Forearm$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY, ]
Wbody_COL<-Wbody[Wbody$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY, ]
WbodyComposition_COL<-WbodyComposition[WbodyComposition$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY, ]
ObesityIndices_COL<-ObesityIndices[ObesityIndices$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY, ]
ScanAnalysis_COL<-ScanAnalysis[ScanAnalysis$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY, ] 
TenYearFxRisk_COL<-TenYearFxRisk[TenYearFxRisk$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY, ]
HipHSA_COL<-HipHSA[HipHSA$PATIENT_KEY%in%PATIENTS_COL$PATIENT_KEY, ]


#save(SillenceCOL,PATIENTS_COL,
#	BMDName, ScanAnalysis_COL, TenYearFxRisk_COL, HipHSA_COL,
#	Spine_COL, Hip_COL, Forearm_COL, Wbody_COL, WbodyComposition_COL, ObesityIndices_COL, 
#	file="BMD.for.COL1.n121.RData")


########################################
table(ScanAnalysis$SCAN_TYPE)
table(ScanAnalysis$SCAN_TYPE[ScanAnalysis$SCANID%in%Spine$SCANID])
table(ScanAnalysis$SCAN_TYPE[ScanAnalysis$SCANID%in%Hip$SCANID])
table(ScanAnalysis$SCAN_TYPE[ScanAnalysis$SCANID%in%HipHSA$SCANID])
table(ScanAnalysis$SCAN_TYPE[ScanAnalysis$SCANID%in%Wbody$SCANID])
table(ScanAnalysis$SCAN_TYPE[ScanAnalysis$SCANID%in%WbodyComposition$SCANID])
table(ScanAnalysis$SCAN_TYPE[ScanAnalysis$SCANID%in%SubRegionBone$SCANID])
table(ScanAnalysis$SCAN_TYPE[ScanAnalysis$SCANID%in%Forearm$SCANID])
table(ScanAnalysis$SCAN_TYPE[ScanAnalysis$SCANID%in%TenYearFxRisk$SCANID])
table(ScanAnalysis$SCAN_TYPE[ScanAnalysis$SCANID%in%ObesityIndices$SCANID])

table(Hip$SCANID%in%Spine$SCANID)










