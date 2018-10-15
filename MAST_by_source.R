library(MAST)
library(scater)
library(data.table)
### Hurdle model fitting 
# ZLM (ridge regression for continuous, get standardized deviance residuals)
t_blood_ind<-which(colData(QC)$cell_type=="T" & colData(QC)$source=="Blood")
t_gut_ind <- which(colData(QC)$cell_type=="T" & colData(QC)$source=="Gut")

make_obj<-function(cell_ind_QC){
fData = data.frame(primerid=rownames(out_norm),biotype=rowData(QC)$biotype)
rownames(fData) = rownames(out_norm)

specificity<-factor(colData(QC)$specificity[cell_ind_QC])
specificity<-relevel(specificity,"Tet-T")
patient<-factor(colData(QC)$patient[cell_ind_QC])

#cell_specificity?
cData = data.frame(wellKey=QC$well[cell_ind_QC],specificity=specificity,patient=patient)
rownames(cData) = colnames(out_norm)[cell_ind_QC]

obj <- FromMatrix(as.matrix(log2(out_norm[,cell_ind_QC]+1)), cData, fData)
colData(obj)$cngeneson <- scale(colSums(assay(obj)>0))
return(obj)}

fit_zlm <- function(obj){
  zlm<-zlm(~specificity + patient + cngeneson, obj, method = "bayesglm", 
            ebayes = TRUE, ebayesControl = list(method = "MLE", model = "H1"))
  return (zlm)
}

#function for likelihood ratio test(zlm,factor_tested='specificitySpecific',p=0.05)
lrt_sig<-function(zlm,factor_tested)#,p=0.05)
{
  summarymat <- summary(zlm, doLRT=factor_tested)
  summarymat <- summarymat$datatable
  fcHurdle <- merge(summarymat[contrast==factor_tested & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summarymat[contrast==factor_tested & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  return (fcHurdle)#[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  #fcHurdleSig <- merge(fcHurdle[fdr<p], as.data.table(mcols(obj2)), by='primerid')
  #return (setorder(fcHurdleSig, fdr)) 
}

get_fdr<-function(fcHurdle,FCTHRESHOLD,obj){
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(obj)), by='primerid')
  setorder(fcHurdleSig, fdr)#fdr<.05:304,fdr<.01:215,cutoff:135
}

FCTHRESHOLD <- log2(1.5)

#run testing function and save the DEG results:blood
obj_blood<-make_obj(t_blood_ind)
zlm_blood<-fit_zlm(obj_blood)
fcHurdle.spe_blood<-lrt_sig(zlm_blood,factor_tested='specificityTet+T')#,p=0.05)
fcHurdleSig_blood<-get_fdr(fcHurdle.spe_blood,FCTHRESHOLD,obj_blood)

obj_gut<-make_obj(t_gut_ind)
zlm_gut<-fit_zlm(obj_gut)
fcHurdle.spe_gut<-lrt_sig(zlm_gut,factor_tested='specificityTet+T')#,p=0.05)
fcHurdleSig_gut<-get_fdr(fcHurdle.spe_gut,FCTHRESHOLD,obj_gut)

#alluvial plot compare DEG_blood with DEG_allCells###################################
library(plyr)
library(ggalluvial)
bulk_blood<-read.table(file = "/Users/yingy_adm/Documents/sc_2018/bulk_result/de_genes.csv",header = T, sep = ",")
allDEG<-merge(fcHurdleSig.spe[,c(1,3)],fcHurdleSig_blood[,c(1,3)],by="primerid",all=TRUE)
colnames(allDEG)<-c('gene_name','Tet+/Tet-ALL','Tet+/Tet-Blood')
allDEG[is.na(allDEG)]<-"Not DE"
allDEG$`Tet+/Tet-ALL`[which(as.numeric(allDEG$`Tet+/Tet-ALL`)<0)]<-"Down"
allDEG$`Tet+/Tet-ALL`[which(as.numeric(allDEG$`Tet+/Tet-ALL`)>0)]<-"Up"
allDEG$`Tet+/Tet-Blood`[which(as.numeric(allDEG$`Tet+/Tet-Blood`)<0)]<-"Down"
allDEG$`Tet+/Tet-Blood`[which(as.numeric(allDEG$`Tet+/Tet-Blood`)>0)]<-"Up"

freqDEG<-as.data.frame(table(allDEG$`Tet+/Tet-ALL`,allDEG$`Tet+/Tet-Blood`))
freqDEG<-cbind(freqDEG,freqDEG[,1])
colnames(freqDEG)<-c('ALL','Blood','Freq',"color")

pdf("/Users/yingy_adm/Documents/sc_2018/plots/alluvial_MAST_DEG_ALLvsBlood.pdf",height=5,width=9)
ggplot(data = freqDEG,aes(axis1 = ALL, axis2 = Blood, weight = Freq)) +
  scale_x_discrete(limits = c("Tet+/Tet-ALL", "Tet+/Tet-Blood"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = color)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal() +
  ggtitle("DEG tested using all cells vs using only cells in blood")
dev.off()

D3E_spe<-read.table("/Users/yingy_adm/Documents/sc_2018/D3E_IO/D3E_out_spe_jul.txt",header = T,sep="\t")
D3E_spe<-D3E_spe[order(D3E_spe$p.value),]
D3E_sig<-D3E_spe[which(D3E_spe$p.value<0.01),c(1,16,17)]
colnames(D3E_sig)[c(2,3)]<-c("log2FC_size","log2FC_freq")
D3E_sig$log2FC_size[which(as.numeric(D3E_sig$log2FC_size)<0)]<-"Down"
D3E_sig$log2FC_size[which(as.numeric(D3E_sig$log2FC_size)>0)]<-"Up"
D3E_sig$log2FC_freq[which(as.numeric(D3E_sig$log2FC_freq)<0)]<-"Down"
D3E_sig$log2FC_freq[which(as.numeric(D3E_sig$log2FC_freq)>0)]<-"Up"
#D3E_sig$GeneID[D3E_sig$log2FC_freq!=D3E_sig$log2FC_size]#show gene list
freqDEG_D3E<-as.data.frame(table(D3E_sig$log2FC_size,D3E_sig$log2FC_freq))
freqDEG_D3E<-cbind(freqDEG_D3E,freqDEG_D3E[,1])
colnames(freqDEG_D3E)<-c('log2FC_size','log2FC_freq','Freq',"color")

pdf("/Users/yingy_adm/Documents/sc_2018/plots/alluvial_D3E_DEG_SIZEvsFREQ.pdf",height=5,width=9)
ggplot(data = freqDEG_D3E,aes(axis1 = log2FC_size, axis2 = log2FC_freq, weight = Freq)) +
  scale_x_discrete(limits = c("size", "freq"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = color)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal() +
  ggtitle("DEG by D3E log2FC_size vs log2FC_freq")
dev.off()

#volcano plot comparing DEG_blood labeled with DEG_bulk_blood
colnames(fcHurdle.spe_blood)<-c("primerid","pvalue","logFC","CI_high","CI_low","FDR")
df <- data.frame(primerid=fcHurdle.spe_blood$primerid,logFC=fcHurdle.spe_blood$logFC,adjusted_P=fcHurdle.spe_blood$FDR)
library(ggplot2)
library(ggrepel)
pdf("/Users/yingy_adm/Documents/sc_2018/plots/volcano_BloodCells_bulkBloodLabeled.pdf",height=5,width=9)
ggplot(data = df, aes(x = logFC, y = -log2(adjusted_P), label = primerid)) + theme_bw()+
  geom_point(colour = "black", size = 0.5) +
  geom_point(data=subset(df,primerid%in%bulk_blood$gene_name[bulk_blood$logFC>0]),colour="red",size = 0.5)+
  geom_point(data=subset(df,primerid%in%bulk_blood$gene_name[bulk_blood$logFC<0]),colour="green",size = 0.5)+
  #geom_text(data = subset(df, adjusted_P<.05 & abs(logFC)>log2(1.5)),size=2)+
  geom_hline(yintercept = 4.321928,colour="grey")+#-log2(0.05)
  geom_vline(xintercept=-0.5849625,colour="grey")+#log2(1.5)
  geom_vline(xintercept=0.5849625,colour="grey")+
  ggtitle("DEA using only cells in blood labeled with DEG_bulk_blood")
dev.off()

#alluvial plot compare with DEG using all data and DEG_bulk_blood###################################