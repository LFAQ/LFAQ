setwd("G:\\InHouseSourceCode\\LFAQ\\LFAQGUI\\LFAQ\\bin\\x64\\Release")
library(ggplot2)
d<-read.table("StandardProteinResults.txt",sep="\t",header = TRUE)
x<-d[d$ExperiemntName=="UPS2_Yeast_R1","SpikedInMols"];
LFAQ<-d[d$ExperiemntName=="UPS2_Yeast_R1","LFAQ"];
LLFAQ<-rep("LFAQ",length(LFAQ))
plotFrame<-data.frame(SpikedIn=x,Label=LLFAQ,LFAQ=LFAQ)
tiff(file="G:\\InHouseSourceCode\\LFAQ\\LFAQGUI\\LFAQ\\bin\\x64\\Release/LFAQOfStandPro.tiff",res=500,width=2000,height=1500,compression = "lzw")
ggplot(plotFrame,aes(SpikedIn,LFAQ,colour=Label))+geom_point()+geom_smooth(method = "lm")+ylab("LFAQ protein intensity(log10)")+xlab("Spiked-in standards amount(log10)")
dev.off()

