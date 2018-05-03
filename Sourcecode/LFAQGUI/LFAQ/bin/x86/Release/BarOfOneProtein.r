setwd("G:\\InHouseSourceCode\\LFAQ\\LFAQGUI\\LFAQ\\bin\\x86\\Release")
library(ggplot2)
d<-read.table("PeptidesOfOneProtein.txt",sep="\t",header = TRUE)
Peptidesindex<-d[,1];
Intensities<-d[,2];
plotFrame<-data.frame(Peptidesindex,Intensities);
tiff(file="G:\\InHouseSourceCode\\LFAQ\\LFAQGUI\\LFAQ\\bin\\x86\\Release/BarOfOneProtein.tiff",res=500,width=2000,height=1500,compression = "lzw")
if(length(Peptidesindex)<10){
ggplot()+geom_bar(aes(Peptidesindex,Intensities),stat = "identity",fill=Peptidesindex)+ ylab("LFAQ peptide intensity")+ xlab("Peptide index")+scale_x_continuous(breaks=Peptidesindex)
}else{
ggplot()+geom_bar(aes(Peptidesindex,Intensities),stat = "identity",fill=Peptidesindex)+ ylab("LFAQ peptide intensity")+ xlab("Peptide index")
}
dev.off()

