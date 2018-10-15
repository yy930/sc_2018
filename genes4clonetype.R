#use randomForest to select most important genes for clonetype
library(randomForest)
library(caret)

celltype<-read.delim("~/Documents/sc_2018/celltype_clean_clonetype.txt",header=T,sep="\t",
                     check.names=FALSE)

#only include 156 cells with clone_exp
X=as.data.frame(t(out_norm[,QC$clone_exp=="T"]))
Y=as.factor(QC$clonetype[QC$clone_exp=="T"])
Y=droplevels(Y)

#10 fold cross validation to optimize ntree in RF.model
flds <- createFolds(Y, k = 8, list = TRUE, returnTrain = FALSE)
ntree <- seq(100, 1000, length=10)

cvm <- data.frame(matrix(nrow=8*10, ncol=3))
colnames(cvm) <- c("foldid", "ntree", "cvm")
cvm[, 1] <- rep(1:8, each=10)
cvm[, 2] <- factor(rep(ntree, 8))

resamp <- function(x,...){if(length(x)==1) x else sample(x,...)} 

for (i in 1:8)
  {for(j in 1:10){
    train.y <- Y[-flds[[i]]]
    if(any(table(train.y)==0)){
      missing_clonetype=names(which(table(train.y)==0))
      add_flds<-unlist(lapply(missing_clonetype, function(x) resamp(which(Y==x),1)))
      flds[[i]]<-flds[[i]][!flds[[i]]%in%add_flds]
    }
    train.y <- Y[-flds[[i]]]
    test.y <- Y[flds[[i]]]
    trainingset.x <- X[-flds[[i]],]
    testset.x <- X[flds[[i]],]
    mymodel <- randomForest(x=trainingset.x, y=train.y, ntree=ntree[j])
    temp <- predict(mymodel, testset.x)
    cvm[(i-1)*8+j,3] <- sum(temp != test.y)/length(temp)
  }}

min(with(cvm, tapply(cvm, ntree, mean)))

#FIT final model using optimized ntree
RF.mod <- randomForest(x=X, y=Y,ntree = 1000)
imp_genes_clonetype<-RF.mod$importance
imp_genes<-rownames(imp_genes_clonetype)[order(imp_genes_clonetype,decreasing = T)][1:50]
write.table(imp_genes,"/Users/yingy_adm/Documents/sc_2018/50imp_genes_clonetype.txt",col.names = F,row.names = F,quote = F)

library(scater)
pdf("/Users/yingy_adm/Desktop/heatmap_by_clonetype.pdf",height=25,width=25)
M3DropExpressionHeatmap(genes=imp_genes, out_norm[,QC$clone_exp=="T"], cell_labels=Y,key_cells=Y)
dev.off()

pdf("/Users/yingy_adm/Desktop/heatmap_by_batch.pdf",height=25,width=25)
M3DropExpressionHeatmap(genes=imp_genes, out_norm[,QC$clone_exp=="T"], cell_labels=QC$batch[QC$clone_exp=="T"])
dev.off()

heatmap_out <- M3DropExpressionHeatmap(genes=imp_genes, out_norm[,QC$clone_exp=="T"], cell_labels=Y,key_cells=Y)
clusters <- M3DropGetHeatmapCellClusters(heatmap_out, k=53)
for (i in unique(clusters)){
  print (celltype[names(clusters[which(clusters==i)]),"clonetype"])
}

library(SC3)
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(out_norm[imp_genes,QC$clone_exp=="T"]),
    logcounts = log2(as.matrix(out_norm[imp_genes,QC$clone_exp=="T"]) + 1)
  ), 
  colData = celltype[!scdata$drop,][QC$clone_exp=="T",])

sce <- sc3(sce, ks = 4:7, biology = TRUE)
sce<- sc3(sce, ks = 7, biology = TRUE)
sc3_plot_consensus(sce, k = 7, show_pdata = "cell_type")
for (i in as.character(1:7)){
  print (sce$clonetype[which(sce$sc3_7_clusters==i)])
}
