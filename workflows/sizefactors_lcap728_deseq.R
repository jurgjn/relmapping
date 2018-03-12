library("DESeq2")

fp <- snakemake@input[[1]]
countData <- read.table(fp, header=TRUE, sep="\t", row.names=1)

cond <- relevel(factor(c( # Hard-coded(!) unlist(lapply(colnames(countData), function(s) substr(s, 1, nchar(s) - 5)))
    "lcap728_wt_emb", "lcap728_wt_emb",
    "lcap728_wt_l1", "lcap728_wt_l1",
    "lcap728_wt_l2", "lcap728_wt_l2",
    "lcap728_wt_l3", "lcap728_wt_l3",
    "lcap728_wt_l4", "lcap728_wt_l4",
    "lcap728_wt_ya", "lcap728_wt_ya",
    "lcap728_glp1_ya", "lcap728_glp1_ya",
    "lcap728_glp1_d3", "lcap728_glp1_d3",
    "lcap728_glp1_d7", "lcap728_glp1_d7",
    "lcap728_glp1_d10", "lcap728_glp1_d10",
    "lcap728_glp1_d14", "lcap728_glp1_d14")), "lcap728_wt_emb")

colData <- data.frame(row.names=colnames(countData), condition=cond)

dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
dds <- estimateSizeFactors(dds)

fp <- snakemake@output[[1]]
write.table(sizeFactors(dds), fp, sep="\t", quote=FALSE, col.names=F)
