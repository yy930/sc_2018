library(rnaseqR)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(stringr)

#data loading
gene_biotype <- read.delim("/Users/yingy_adm/Documents/sc_count/gene_biotype.txt",row.names = 1)
colnames(gene_biotype)[1]<-"Gene"
celltype<-read.delim("~/Documents/sc_2018/celltype_clean_clonetype.txt",header=T,sep="\t",
                     check.names=FALSE)

count2017sep <- read.delim("/Users/yingy_adm/Documents/sc_2017/genecounts_2017sep.txt",
                           header=T,sep=" ",row.names=1,check.names=FALSE)
tpm2017sep <- read.delim("/Users/yingy_adm/Documents/sc_2017/genetpm_2017sep.txt",
                         header=T,sep=" ",row.names=1,check.names=FALSE)
mapping2017sep<- read.table(file = "/Users/yingy_adm/Documents/sc_2017/meta_mapping_2017sep.txt",
                            row.names=1,header = F, sep = "")
count2018feb <- read.delim("/Users/yingy_adm/Documents/sc_2018/genecount_2018feb.txt",
                    header=T,sep=" ",row.names=1,check.names=FALSE)
tpm2018feb <- read.delim("/Users/yingy_adm/Documents/sc_2018/genetpm_2018feb.txt",
                    header=T,sep=" ",row.names=1,check.names=FALSE)
mapping2018feb<- read.table(file = "/Users/yingy_adm/Documents/sc_2018/meta_mapping_2018feb.txt",
                            row.names=1,header = F, sep = "")
tpm2018sep <- read.delim("/Users/yingy_adm/Documents/sc_2018/batch3/genetpm.txt",header=T,sep=" ",
                         row.names=1,check.names=FALSE)


#aggregate transtripts to genes
count1 <- aggregate_tpm(count2017sep)
tpm1 <- aggregate_tpm(tpm2017sep)
count2 <- aggregate_tpm(count2018feb)
tpm2 <- aggregate_tpm(tpm2018feb)
count3 <- aggregate_tpm(count2018sep)
tpm3 <- aggregate_tpm(tpm2018sep)
# import sample information for batch3
samples3 <- import_sample_names('/Users/yingy_adm/Documents/sc_2018/batch3/trans',"samplenames.txt")
# import readcounts using tximport
txi3 <- tx_import('/Users/yingy_adm/Documents/sc_2018/batch3/trans',samples=samples3,file_type='quant.sf',tx2gene=tx2gene_noalt)
txi3$counts->count3
count3<-count3[,substring(celltype$well[1:384],3)]

#add ercc matrix
add_ercc <- function(exprmat,agg_mat,batch){
  rownames(exprmat)=str_split_fixed(rownames(exprmat),"[|]",2)[,1]
  erccmat<-exprmat[grep("ERCC-",rownames(exprmat)),]
  exprmat_sum <- rbind(erccmat,agg_mat) 
  return(exprmat_sum)
}
count1<-add_ercc(count2017sep,count1,batch = 1)
count2<-add_ercc(count2018feb,count2,batch = 2)
erccmat <- count2018feb[grep("ERCC-",rownames(count2018feb)),]#erccmat with all 0
erccmat[erccmat>0] <- 0
colnames(erccmat)<-colnames(count3)
count3 <- rbind(erccmat,count3) 
tpm1 <- add_ercc(tpm2017sep,tpm1,batch = 1)
tpm2 <- add_ercc(tpm2018feb,tpm2,batch = 2)
tpm3 <- add_ercc(tpm2018sep,tpm3,batch = 3)

#combine all 3 batches
allcount <- cbind(count1[1:288],count2,count3)
alltpm <- cbind(tpm1[1:288],tpm2,tpm3)
colnames(allcount) = colnames(alltpm) = celltype$well

#aggregate gene_biotype 
idx_ercc<-grep("ERCC-",rownames(gene_biotype))
gene_biotype_ercc <- data.frame(gene_id=gene_biotype[idx_ercc,1],tx_id=gene_biotype[idx_ercc,1],
                           Gene=gene_biotype[idx_ercc,1],Biotype=gene_biotype[idx_ercc,2,])
gene_biotype_agg <- rbind(gene_biotype_ercc,tx2gene_noalt_hgnc)
allcount[match(gene_biotype_agg$gene_id,rownames(allcount)),]
sum(gene_biotype_agg$gene_id==rownames(allcount))
sum(gene_biotype_agg$gene_id==rownames(alltpm))

#save data
write.table(allcount, file = "~/Documents/sc_2018/allcount_geneid.txt",sep="\t", row.names = T)
write.table(alltpm, file = "~/Documents/sc_2018/alltpm_geneid.txt",sep="\t", row.names = T)
write.table(gene_biotype_agg, file = "~/Documents/sc_2018/agg_geneid_biotype.txt",sep="\t", row.names = T)