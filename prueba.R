
load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3",
  "RISmed",
  "ggplot2",
  "dplyr")

lapply(load.libs, require, character.only = TRUE)

freqs <- read.csv("./frecuencias_alelicas_lab.csv")
exoma <- read.csv("/repositorio/exomas/pipeline/DX010-23/output_dirs/post_process_results/DX010-23_GATK_CODIGO_HPO/DX010-23_GATK_CODIGO_HPO.csv")

genes_univ <- unique(freqs$Gene.refGene)
exoma_goi <- unique(exoma$Gene.refGene)

goi <- clusterProfiler::bitr(exoma_goi,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
universe <- clusterProfiler::bitr(genes_univ,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)



ewp.up <- clusterProfiler::enrichWP(
  goi[,2],
  universe = universe[,2],
  organism = "Homo sapiens",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff; relaxed for demo purposes
)

ewp.up <- DOSE::setReadable(ewp.up, org.Hs.eg.db, keyType = "ENTREZID")

resultados.1 <- ewp.up@result
resultados.1 <- resultados.1[order(resultados.1$pvalue),]

ggplot(resultados.1[1:20,], aes(x=Description, y=Count, fill=pvalue)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "", fill = "p.value") +
  theme(axis.text=element_text(size=11))  

ewp.up <- DOSE::setReadable(ewp.up, org.Hs.eg.db, keyType = "ENTREZID")
resultados.1 <- ewp.up@result
resultados.1 <- resultados.1[order(resultados.1$pvalue),]

goi2 <- unlist(strsplit(resultados.1$geneID,"/"))
goi2 <- clusterProfiler::bitr(goi2,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

dose <- DOSE::enrichDO(goi2[,2])
dose <- DOSE::setReadable(dose, org.Hs.eg.db, keyType = "ENTREZID")

resultado3 <- dose@result[order(dose@result$p.adjust),]
ggplot(resultado3[1:20,], aes(x=Description, y=Count, fill=-log(p.adjust))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_continuous(low="blue", high="red") +
  labs(x = "", y = "", fill = "p.adjust") +
  theme(axis.text=element_text(size=11))  

pathos.list <- list()
pathos <- resultado3
for(l in 1:dim(pathos)[1]){
  repeticion <- length(unlist(strsplit(pathos$geneID[l],"/")))
  pathos.list[[l]] <- data.frame(patologia =rep(pathos$Description[l],repeticion),
                                 Gene.refGene = unlist(strsplit(pathos$geneID[l],"/")),
                                 p_value = rep(pathos$pvalue[l],repeticion),
                                 p_adjust = rep(pathos$p.adjust[l],repeticion))
  
  
}

pathos.df <- (Reduce(rbind,pathos.list))

final <- left_join(exoma,pathos.df,by="Gene.refGene")

# final <- final[which(!is.na(final$CLNDN)),]
# final <- final[which(grepl("benign",ignore.case = T,final$CLNSIG)==F),]
filtrado_1 <- function(X, threshold) {
  X <-
    X[which(X$Func.refGene == "exonic;splicing" |
              X$Func.refGene == "exonic"),]
  X <-
    X[which(
      X$ExonicFunc.refGene == "nonsynonymous SNV" |
        X$ExonicFunc.refGene == "stopgain"|
        X$ExonicFunc.refGene == "exonic" |
        X$ExonicFunc.refGene == "nonsynonymous SNV" |
        X$ExonicFunc.refGene == "frameshift deletion" |
        X$ExonicFunc.refGene == "frameshift insertion" |
        X$ExonicFunc.refGene == "startloss" |
        X$ExonicFunc.refGene == "startgain" |
        X$ExonicFunc.refGene == "stoploss"
    ),]
  
  X <- X[which(X[, grep("freq", colnames(X))] < threshold),]
  
  return(X)
  
}
final <- filtrado_1(final,1)
final <- final[which(!is.na(final$CLNDN)),]
final <- final[which(grepl("benign",ignore.case = T,final$CLNSIG)==F),]

# write.csv(final,"~/Desktop/final_DX004-23.csv")
# 
# final2 <- final


  
  
  