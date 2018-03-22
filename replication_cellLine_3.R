library(broom)

cellLineData <- read.table("rna_celline.tsv", header = T, sep = "\t", stringsAsFactors = F)
head(cellLineData)
table(cellLineData$Unit)
load('Lupus_RF_v2.Rdata')
CDinterestedNumbers <- c(20, 22, 30)
cdNames <- paste0("CD", CDinterestedNumbers)
CDsEntrez <- as.character(CDsannotated[
  CDsannotated$CD_NAME %in% cdNames, "ENTREZ_GENE_ID"])
nCDs <- length(CDsEntrez)
rownames(combinedIDs) <- combinedIDs$ENTREZ_GENE_ID
CDsGeneName <- as.character(CDsannotated[
  CDsannotated$CD_NAME %in% cdNames, "NCBI_NAME"])


for (i in 1:nCDs){
  CDi <- cdNames[i]
  CDiEntrez <- CDsEntrez[i]
  CDiGeneName <- CDsGeneName[i]
  # genesToVerifyiDf <- read.csv(paste0("result61/", CDi, "indMultinomialRanked.csv"), row.names = 1)
  # genesToVerifyi <- rownames(genesToVerifyiDf[genesToVerifyiDf$adjustedPval < 0.05,])
  # genesToVerifyi <- unique(substr(genesToVerifyi, 1, nchar(genesToVerifyi)-2))
  # genesToVerifyiEntrez <-  as.character(combinedIDs[combinedIDs$Gene %in% genesToVerifyi, "ENTREZ_GENE_ID"])
  cdigenes <- read.csv(paste0("resultsMultiSURF/", CDi, "indBinomialRankedLambdamin.csv"), stringsAsFactors = F)[,1]
  cdigenes <- gsub(pattern = "\\**", replacement = "", cdigenes)
  genesOfInterestCDi <- as.character(combinedIDs[
    combinedIDs$Gene %in% cdigenes, "ENTREZ_GENE_ID"])
  # genesOfInterestCD20 <- c(CDs20RelatedEntrez)
  # mapInterest <- genesLymp1[c(genesOfInterestCDi, CDiEntrez),]
  # mapInterest <- mapInterest[!is.na(mapInterest$PROBEID),]
  # # mapInterest <- mapInterest[!duplicated(mapInterest$ENTREZID),]
  # sampleInteresti <- t(exprData1[mapInterest$PROBEID,])
  # colnames(sampleInteresti) <- mapInterest$ENTREZID
  # # sampleInteresti <- t(sampleInteresti)
  # 
  # cellLinei <- cellLineData[cellLineData$Gene.name %in% cdigenes,]
  # dim(cellLinei)
  
  CDiDat <- cellLineData[cellLineData$Gene.name == CDiGeneName, ]
  # CDvali <- CDiDat$Value
  # markCDi <- quantile(CDvali, c(.2, .8), na.rm = T) 
  # # cdtrinary <- factor((CDvali < markCDi[1]) + 2*(CDvali > markCDi[2])) # low 1, high 2
  # cdtrinary <- factor((CDvali < markCDi[1]| CDvali > markCDi[2])) # low 1, high 2
  # levels(cdtrinary) <- c("Non-aberrant", "Aberrant")
  # checkTrinary <- cbind(cdtrinary, CDvali)
  # cdidf <- data.frame(CDiDat, cd = cdtrinary)
  # colnames(cdidf)[-ncol(cdidf)] <- colnames(CDiDat)
  # # hist(cor(importantSamples, myCDs$CD20))
  # xReg <- as.matrix(cdidf[,colnames(cdidf)!="cd"])
  
  
  
  multiSums <- data.frame()
  # impGenesCDi <- setdiff(colnames(sampleInteresti), CDiEntrez)
  impGenesCDi <- cdigenes
  for (j in 1:length(impGenesCDi)){
    genej <- impGenesCDi[j]
    datj <- cellLineData[cellLineData$Gene.name == genej, ]
    if (nrow(datj) > 0){
      datjMerge <- merge(datj, CDiDat, by = "Sample")
      multiSum <- tidy(cor.test(datjMerge$Value.x, datjMerge$Value.y))
      # colnames(datj)[1] <- impGenesCDi[j]
      # uniMultinomial <- glm(cdtrinary ~ gene, data = datj, family = "binomial")
      # multiSum <- tidy(uniMultinomial)
      # multiSum$oddRatio <- exp(multiSum$estimate)
      rownames(multiSum) <- genej
      # # colnames(multiSum)[3] <- "oddRatio" # relative risk, equal exp(coef(uniMultinomial))
      # multiSum$oddRatioL <- multiSum$oddRatio*exp(-1.96*multiSum$std.error)
      # multiSum$oddRatioR <- multiSum$oddRatio*exp(1.96*multiSum$std.error)
      multiSums <- rbind(multiSums, multiSum)
    }

  }
  
  # multiSums$adjustedPval <- multiSums$p.value*length(impGenesCDi)
  multiSums$adjustedPval <- p.adjust(multiSums$p.value, method = "BH")
  # # multiSums[multiSums$y.level == "High", "adjustedPval"] <- 
  # #   multiSums[multiSums$y.level == "High", "p.value"]*length(impGenesCDi)
  # # multiSums[multiSums$y.level == "Low", "adjustedPval"] <- 
  # #   multiSums[multiSums$y.level == "Low", "p.value"]*length(impGenesCDi)
  multiSumsRanked <- multiSums[order(multiSums$p.value),]
  # 
  # multiSumsRanked$needCorrection <- multiSumsRanked$oddRatio < 1
  # multiSumsRanked$correctedOddsRatio <- multiSumsRanked$oddRatio
  # multiSumsRanked[multiSumsRanked$oddRatio < 1, "correctedOddsRatio"] <-
  #   1/multiSumsRanked[multiSumsRanked$oddRatio < 1, "oddRatio"]
  # 
  # multiSumsRanked$correctedOddsRatioL <- multiSumsRanked$oddRatioL
  # multiSumsRanked[multiSumsRanked$oddRatio < 1, "correctedOddsRatioL"] <-
  #   1/multiSumsRanked[multiSumsRanked$oddRatio < 1, "oddRatioR"]
  # 
  # multiSumsRanked$correctedOddsRatioR <- multiSumsRanked$oddRatioR
  # multiSumsRanked[multiSumsRanked$oddRatio < 1, "correctedOddsRatioR"] <-
  #   1/multiSumsRanked[multiSumsRanked$oddRatio < 1, "oddRatioL"]  
  # multiSumsRanked$term <- rownames(multiSumsRanked)
  # rownames(multiSumsRanked)[multiSumsRanked$oddRatio < 1] <- 
  #   paste0(rownames(multiSumsRanked)[multiSumsRanked$oddRatio < 1], "*")
  # 
  # colnames(multiSumsRanked) <- c("term", "estimate", "std.error",  "statistic", "p.value",  "oddRatioUncorrected", "oddRatioLUncorrected", 
  #                                "oddRatioRUncorrected", "adjustedPval", "needCorrection","oddRatio",
  #                                "oddRatioL", "oddRatioR")
  # 
  # 
  write.csv(multiSumsRanked, file = paste0(CDi, "CellLine.csv"))
}
















