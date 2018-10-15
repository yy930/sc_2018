library(SetRank) 
library(GeneSets.Homo.sapiens)
library(GO.db)
library(AnnotationDbi)
library(Homo.sapiens)
#ls("package:GeneSets.Homo.sapiens")

fcHurdleSig.spe$geneID <- gene_biotype$gene_id[match(fcHurdleSig.spe$primerid,rownames(gene_biotype))]
inputData = fcHurdleSig.spe[order(fcHurdleSig.spe$fdr),c(8,3,2,6)]
colnames(inputData) = c("geneID","log2FoldChange","pval","padj") 

ensembl2EntrezID = createIDConverter("Homo.sapiens", "ENSEMBL","ENTREZID")
referenceSet = ensembl2EntrezID(fcHurdle.spe$geneID)

options(mc.cores=1)
#collection = buildSetCollection(GOBP, KEGG, referenceSet = referenceSet, maxSetSize = 500) 
#save(collection, file="/Users/yingy_adm/Documents/sc_2018/collection.Rda")
load("~/Documents/sc_2018/collection.Rda")

inputData$entrezID = ensembl2EntrezID(inputData$geneID, na.rm=FALSE,drop.ambiguous=TRUE)

network = setRankAnalysis(inputData$entrezID, collection, use.ranks = TRUE, setPCutoff = 0.01, fdrCutoff = 0.05)
save(network, file="/Users/yingy_adm/Documents/sc_2018/network.Rda")

load("~/Documents/sc_2018/network.Rda")
load("~/Downloads/Homo.sapiens_18.09.14.Rda")
exportSingleResult(network, inputData$entrezID, collection,
                   "myResults", IDConverter=ensembl2EntrezID, "setrank_results")
