celltype<-read.delim("~/Documents/sc_2018/celltype_2018all.txt",header=T,sep="\t",check.names=FALSE)
clonetype<-unique(celltype$clonetype)[3:62]
mat_clonetype<-as.data.frame(matrix(NA,nrow=nrow(celltype),ncol=length(clonetype)))
colnames(mat_clonetype)<-clonetype
celltype<-cbind(celltype,mat_clonetype)
for (i in clonetype){
  celltype[[i]]<-"others"
  celltype[[i]][celltype$clonetype==get("i")]<-get("i")
  celltype[[i]][celltype$clonetype=="N"]<-NA
}
rm(mat_clonetype,clonetype)
