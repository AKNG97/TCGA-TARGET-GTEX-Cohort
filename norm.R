## Binded and Individual Norm
#Filtering genes
dataFiltBALLpb <- TCGAanalyze_Filtering(tabDF = rawBALLpb,
                                        method = "quantile",
                                        qnt.cut = 0.25)
dataFiltAML <- TCGAanalyze_Filtering(tabDF = rawAML,
                                     method = "quantile",
                                     qnt.cut = 0.25)
dataFiltNB <- TCGAanalyze_Filtering(tabDF = rawNB,
                                    method = "quantile",
                                    qnt.cut = 0.25)
dataFilBinded <- TCGAanalyze_Filtering(tabDF = readyBinded,
                                    method = "quantile",
                                    qnt.cut = 0.25)

threshold <- round(dim(dataFilBinded)[2]/2)
ridx <- rowSums(dataFilBinded == 0) <= threshold
dataFilBinded <- dataFilBinded[ridx, ]
print(dim(dataFilBinded))
ridx <- rowMeans(dataFilBinded) >= 3.45
dataFilBinded <- dataFilBinded[ridx, ]

print(dim(dataFiltBALLpb))
print(dim(dataFiltAML))
print(dim(dataFiltNB))
print(dim(dataFilBinded))

inter <- intersect(rownames(dataFilBinded), annot$ensembl_gene_id)
f_Binded <- dataFilBinded[rownames(dataFilBinded) %in% inter,]
annotBinded <- annot[annot$ensembl_gene_id %in% inter,]

#Annot & BALL
interBALL <- intersect(rownames(dataFiltBALLpb), annot$ensembl_gene_id)
f_BALL <- dataFiltBALLpb[rownames(dataFiltBALLpb) %in% interBALL,]
annotBALL <- annot[annot$ensembl_gene_id %in% interBALL,]

#Annot & AML
interAML <- intersect(rownames(dataFiltAML), annot$ensembl_gene_id)
f_AML <- dataFiltAML[rownames(dataFiltAML) %in% interAML,]
annotAML <- annot[annot$ensembl_gene_id %in% interAML,]

#Annot & NB
interNB <- intersect(rownames(dataFiltNB), annot$ensembl_gene_id)
f_NB <- dataFiltNB[rownames(dataFiltNB) %in% interNB,]
annotNB <- annot[annot$ensembl_gene_id %in% interNB,]

#Data frame to Matrix
BALL_matrix <- data.matrix(f_BALL)
AML_matrix <- data.matrix(f_AML)
NB_matrix <- data.matrix(f_NB)
Binded_matrix <- data.matrix(f_Binded)

saveRDS(f_Binded, file="f_Binded.RDS")
write.table(annotBinded,"annotBinded.txt",sep="\t")

#Normalizing Binded
ln.data <- withinLaneNormalization(Binded_matrix, annotBinded$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data , annotBinded$GC, which = "full")
norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(factorsB$group))
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n", logtransf = FALSE)
Binded_norm <- exprs(mydata2corr1)

#Normalizing AML
ln.data <- withinLaneNormalization(AML_matrix, annotAML$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data , annotAML$GC, which = "full")
norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
#noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(factorAML$group))
#mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n", logtransf = FALSE, batch = TRUE)
AML_norm_test <- norm.counts

#Normalizing BALL
ln.data <- withinLaneNormalization(BALL_matrix, annotBALL$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data , annotBALL$GC, which = "full")
norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
#noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(factorAML$group))
#mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n", logtransf = FALSE, batch = TRUE)
BALL_norm_test <- norm.counts

#Normalizing NB
ln.data <- withinLaneNormalization(NB_matrix, annotNB$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data , annotNB$GC, which = "full")
norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
#noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(factorAML$group))
#mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n", logtransf = FALSE, batch = TRUE)
NB_norm_test <- norm.counts

dim(NB_norm_test)
dim(AML_norm_test)
dim(BALL_norm_test)
dim(Binded_norm)

#PCAs
mydata_NB = readData(NB_norm_test, factors = as.data.frame(factorNB$group))
myPCA_NB = dat(mydata_NB, type = "PCA")
explo.plot(myPCA_NB)

mydata_NB_before = readData(f_NB, factors = as.data.frame(factorNB$group))
myPCA_NB_before = dat(mydata_NB_before, type = "PCA")
explo.plot(myPCA_NB_before)

## Paired up Norm
# Filtering genes
FiltAMLvBALL <- TCGAanalyze_Filtering(tabDF = AMLvBALL,
                                        method = "quantile",
                                        qnt.cut = 0.25)

threshold <- round(dim(FiltAMLvBALL)[2]/2)
ridx <- rowSums(FiltAMLvBALL == 0) <= threshold
FiltAMLvBALL <- FiltAMLvBALL[ridx, ]
print(dim(FiltAMLvBALL))
ridx <- rowMeans(FiltAMLvBALL) >= 3.45
FiltAMLvBALL <- FiltAMLvBALL[ridx, ]

print(dim(FiltAMLvNB))
print(dim(FiltBALLvNB))
print(dim(FiltAMLvBALL))
