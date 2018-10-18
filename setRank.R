library(SetRank) 
library(GeneSets.Homo.sapiens)
library(GO.db)
library(AnnotationDbi)
library(Homo.sapiens)
#ls("package:GeneSets.Homo.sapiens")

#fcHurdleSig.spe$geneID <- gene_biotype$gene_id[match(fcHurdleSig.spe$primerid,rownames(gene_biotype))]
#inputData = fcHurdleSig.spe[order(fcHurdleSig.spe$fdr),c(8,3,2,6)]
inputData = fcHurdleSig.spe[order(fcHurdleSig.spe$fdr),c(1,3,2,6)]
colnames(inputData) = c("geneID","log2FoldChange","pval","padj") 

ensembl2EntrezID = createIDConverter("Homo.sapiens", "ENSEMBL","ENTREZID")
#referenceSet = ensembl2EntrezID(fcHurdle.spe$geneID)
referenceSet = ensembl2EntrezID(fcHurdle.spe$primerid)

options(mc.cores=2)
collection = buildSetCollection(GOBP, KEGG, referenceSet = referenceSet, maxSetSize = 500) 
save(collection, file="/Users/yingy_adm/Documents/sc_2018/collection1.Rda")
load("~/Documents/sc_2018/collection.Rda")

inputData$entrezID = ensembl2EntrezID(inputData$geneID, na.rm=FALSE,drop.ambiguous=TRUE)

network = setRankAnalysis(inputData$entrezID, collection, use.ranks = TRUE, setPCutoff = 0.05, fdrCutoff = 0.05)
#save(network, file="/Users/yingy_adm/Documents/sc_2018/network.Rda")

load("~/Documents/sc_2018/network.Rda")

exportSingleResult(network, inputData$entrezID, collection,
                   "oct15", IDConverter=ensembl2EntrezID, "setrank_results")

writeMembership <- function(outputFile, selectedGenes, collection, network, 
                            IDConverter=NULL) {
  membershipBool = membershipTable(selectedGenes, collection, network)
  if (is.null(membershipBool)) {
    return(NULL)
  }
  membership = matrix(".", nrow=nrow(membershipBool), 
                      ncol=ncol(membershipBool), dimnames=dimnames(membershipBool))
  membership[membershipBool == TRUE] = "X"
  if (!is.null(IDConverter)) {
    rownames(membership) = IDConverter(rownames(membership))
  }
  membership = orderTable(membership)
  prettyTable = cbind(rownames(membership),membership)
  colnames(prettyTable)[1] = "gene"
  write.table(prettyTable, outputFile, sep="\t", row.names=FALSE, 
              quote=FALSE)
}

orderTable <- function(membershipTable){
  membershipTableMod <- ifelse(membershipTable=="X",1,0)
  if (nrow(membershipTableMod) > 1)  {
    clustering = hclust(dist(membershipTableMod))
    geneOrder = clustering$labels[clustering$order]
  } else {
    geneOrder = 1
  }
  if (ncol(membershipTableMod) > 1) {
    clustering = hclust(dist(t(membershipTableMod)))
    catOrder = clustering$labels[clustering$order]
  } else {
    catOrder = 1
  }
  return(membershipTable[geneOrder,catOrder])
}
writeMembership("/Users/yingy_adm/Documents/sc_2018/setrank_results/myResults_membership.txt",inputData$entrezID, collection, network,IDConverter=createIDConverter("Homo.sapiens", "ENTREZID","SYMBOL"))

load("/Users/yingy_adm/Documents/sc_2018/setrank_results/Homo.sapiens_18.09.14.Rda")
topTables = list(setrank=inputData) 
networks = list(setrank=network)

fieldNames = c("geneID"="entrezID", "symbol"="SYMBOL", "logFC"="log2FoldChange", "p"="padj")
exportGeneNets(topTables, networks, collection, string, "gene_nets", fields = fieldNames)
