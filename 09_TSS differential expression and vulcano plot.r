library(CAGEfightR)
library(rtracklayer)
library(tidyverse)
library(DESeq2)
library(edgeR)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28


scripts <- c("batchReadTrackData", "normalizeGR", "saveGRangesAsBedGraph", "expandGRtoUnitWidth")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Kindgren_et_al_2019/blob/master/", script, ".R?raw=TRUE") %>% devtools::source_url()
}
devtools::source_url("https://github.com/Maxim-Ivanov/TranscriptomeReconstructoR/blob/main/R/Utility_functions.R?raw=TRUE") 
devtools::source_url("https://github.com/Maxim-Ivanov/Nielsen_et_al_2018/blob/master/Custom_assignTxType_function_for_CAGEfightR.R?raw=TRUE")
devtools::source_url("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/find_matched_control.R?raw=TRUE")


path <- "."  # change to the directory containing merged TSS-seq Bedgraph files obtained from 02_TSS-seq processing.R

triple_fw<- list.files(".", pattern = "fw.bg") # change to the directory containing merged TSS-seq Bedgraph files obtained from 02_TSS-seq processing.sh
triple_rev<- list.files(".", pattern = "rev.bg") # change to the directory containing merged TSS-seq Bedgraph files obtained from 02_TSS-seq processing.sh

batched_fwd <- batchReadTrackData(triple_fw, dir = ".", format = 'bedGraph', strand = "+", seqinfo = seqinfo(txdb)) # change to the directory containing merged TSS-seq Bedgraph files obtained from 02_TSS-seq processing.sh
batched_rev <- batchReadTrackData(triple_rev, dir = ".", format = 'bedGraph', strand = "-", seqinfo = seqinfo(txdb))# change to the directory containing merged TSS-seq Bedgraph files obtained from 02_TSS-seq processing.sh

# Shorten the names:
names(batched_fwd) <- names(batched_fwd) %>% str_replace("_R1_001_trimmed_UMI_Aligned.sortedByCoord.out_mapq10_dedup_fw.bg.gz", "")
names(batched_rev) <- names(batched_rev) %>% str_replace("_R1_001_trimmed_UMI_Aligned.sortedByCoord.out_mapq10_dedup_rev.bg.gz", "")

bg_data <- mapply(function(x, y) { sort(c(x, y)) }, batched_fwd, batched_rev, SIMPLIFY = FALSE)
bg_data_merged_seedlings <- merge_GRanges(bg_data[1:3])
bg_data_merged_siliques <- merge_GRanges(bg_data[4:6])

bg_data_merged_seedlings_norm <- normalizeGR(bg_data_merged_seedlings)
bg_data_merged_siliques_norm <- normalizeGR(bg_data_merged_siliques)

saveGRangesAsBedGraph(bg_data_merged_seedlings_norm, "seedlings-IGV.bedgraph.gz")
saveGRangesAsBedGraph(bg_data_merged_siliques_norm, "siliques-IGV.bedgraph.gz")

# Skip singletons (ranges with a single tag):
bg_data <- lapply(bg_data, function(gr) { return(gr[score(gr) > 1]) })

# Expand ranges with width > 1 to unit width:
expand_bg_data <- lapply(bg_data, expandGRtoUnitWidth)

# Export the expanded GRanges as Plus and Minus BigWig files (to be used as input for CAGEfightR):
for (i in seq_along(expand_bg_data)) {
  data <- expand_bg_data[[i]]
  name <- names(expand_bg_data)[[i]]
  export(data[strand(data) == "+"], paste0(name, "_unitWidth_Plus.bw"), format = "BigWig")
  export(data[strand(data) == "-"], paste0(name, "_unitWidth_Minus.bw"), format = "BigWig")
}

# Call TSS clusters by CAGEfightR:
bw_plus_filenames <- list.files(".", pattern = "unitWidth_Plus.bw$")
bw_minus_filenames <- list.files(".", pattern = "unitWidth_Minus.bw$")
bw_plus <- BigWigFileList(bw_plus_filenames)
bw_minus <- BigWigFileList(bw_minus_filenames)
sample_names <- bw_plus_filenames %>% str_replace("_unitWidth_Plus.bw", "") %>% str_replace("STRIPE-Seq_", "") # better to use %>% piping instead of multiple nested calls
names(bw_plus) <- sample_names
names(bw_minus) <- sample_names

# For DE calling, the "Genotype" column in design_matrix must contain both conditions to compare ("wt" vs "mut"):
design_matrix <- data.frame("Name" = sample_names, 
                            "BigWigPlus" = bw_plus_filenames, 
                            "BigWigMinus" = bw_minus_filenames, 
                            "Genotype" = c("seedlings", "siliques") %>% rep(each = 3) %>% as.factor() %>% relevel(ref = "seedlings"), # better if it is a factor with "wt" being the reference level
                            row.names = sample_names)

# Quantify all tag clusters (TCs):
ctss <- quantifyCTSSs(plusStrand = bw_plus, minusStrand = bw_minus, design = design_matrix, genome = seqinfo(txdb))

# Call candidate TSSs:
tss <- quickTSSs(ctss)

# Annotate TSS and enhancers by genomic features (observe that a custom assignTxType() function is used):
rowRanges(tss)$txType <- suppressWarnings(assignTxType_custom(rowRanges(tss), txdb=txdb, asFactor=TRUE))

# Annotate TCs by gene IDs:
tss <- suppressWarnings(assignGeneID(tss, geneModels=txdb))

# Write all TSSs identified into a bed file
tss_gr <- rowRanges(tss)
export(tss_gr, "All_called_TSS.bed", format = "BED")

# Use DESeq2:
dds <- DESeqDataSet(tss, design = ~ Genotype)
dds <- DESeq(dds)

# Extract DE results:
res <- results(dds, contrast = c("Genotype", "siliques", "seedlings"))
summary(res)

# Make the decisions on differential expression:
log2FoldChange <- res$log2FoldChange
padj <- res$padj

# Add decisions and padj to tss_gr:
mcols(tss_gr)$Log2FC <- log2FoldChange
mcols(tss_gr)$padj <- padj

# Add decisions on the sample specificity of called TSS:
mcols(tss_gr)$nz_seedlings <- rowSums(assay(tss)[,c(1:3)]>1) 
mcols(tss_gr)$nz_siliques <- rowSums(assay(tss)[,c(4:6)]>1)
both <- mcols(tss_gr)$nz_seedlings>1 & mcols(tss_gr)$nz_siliques>1 
seedlings_only <- mcols(tss_gr)$nz_siliques<=1 & mcols(tss_gr)$nz_seedlings>1
siliques_only <- mcols(tss_gr)$nz_seedlings<=1 & mcols(tss_gr)$nz_siliques>1

mcols(tss_gr)$seedlings_siliques <- "ns"
mcols(tss_gr)$seedlings_siliques[both] <- "both"
mcols(tss_gr)$seedlings_siliques[seedlings_only] <- "seedlings_only"
mcols(tss_gr)$seedlings_siliques[siliques_only] <- "siliques_only"

# Calculate cpm/tpm values
cm <- as.data.frame(cpm(assay(tss)))
cm$key <- rownames(cm)
orig <- data.frame("key"=names(rowRanges(tss)))
x <- left_join(orig, cm, by=c("key"))

# Add tpm results to the RSE object:
mcols(tss_gr) <- cbind(mcols(tss_gr), x[-1])

# Calculate mean tpm of samples
mcols(tss_gr) <- cbind(mcols(tss_gr), seedlings_mean = rowMeans(cm[,c(1:3)]))
mcols(tss_gr) <- cbind(mcols(tss_gr), siliques_mean = rowMeans(cm[,c(4:6)]))

write.csv(tss_gr, "all_fold2_padj0.1.csv")

###deduplication - keeping the genes with lowest padj (highest difference)####
dedup <- tss_gr
dedup_a <- dedup[order(dedup$geneID, abs(dedup$padj) ), ] #sort by id and of abs(value)
dedup_si_a <- dedup_a[!duplicated(dedup_a$geneID), ]
dedup_si_b <- dedup_si_a[mcols(dedup_si_a)$seedlings_siliques %in% c("both", "seedlings_only", "siliques_only")]
write.csv(dedup_si_b, "For_Deseq2.csv")

### Generating Volcano Plot
de_res <- read.table("For_Deseq2.csv", sep = ";", dec = ".", header=TRUE)
par(mar=c(5,5,10,4), cex=1.0, cex.main=1.4, cex.axis=1.0, cex.lab=1.4, lwd=2.0)
with(de_res, plot(Log2FC, -log10(padj), pch=16, cex = 0.5, col="#999999", main="Siliques vs Seedling", xlim=c(-17,17)))
with(subset(de_res, padj<.05 & Log2FC<=-2), points(Log2FC, -log10(padj), pch=16, cex = 0.5, col="#003366"))
with(subset(de_res, padj<.05 & Log2FC>=2), points(Log2FC, -log10(padj), pch=16, cex = 0.5, col="#990000"))
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(de_res$padj[de_res$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)

dev.copy2pdf(file = "SiliquesVsSeedling_Volcano.pdf")

### Generating distribution plot
sum(dedup_si_c$txType == "threeUTR" & dedup_si_c$DE == "Up")
sum(dedup_si_c$txType == "intergenic" & dedup_si_c$DE == "Up")
sum(dedup_si_c$txType == "antisense" & dedup_si_c$DE == "Up")
sum(dedup_si_c$txType == "promoter" & dedup_si_c$DE == "Up")
sum(dedup_si_c$txType == "proximal" & dedup_si_c$DE == "Up")
sum(dedup_si_c$txType == "fiveUTR" & dedup_si_c$DE == "Up")
sum(dedup_si_c$txType == "intron" & dedup_si_c$DE == "Up")
sum(dedup_si_c$txType == "exon" & dedup_si_c$DE == "Up")

sum(dedup_si_b$txType == "threeUTR" & dedup_si_b$DE == "Down")
sum(dedup_si_b$txType == "intergenib" & dedup_si_b$DE == "Down")
sum(dedup_si_b$txType == "antisense" & dedup_si_b$DE == "Down")
sum(dedup_si_b$txType == "promoter" & dedup_si_b$DE == "Down")
sum(dedup_si_b$txType == "proximal" & dedup_si_b$DE == "Down")
sum(dedup_si_b$txType == "fiveUTR" & dedup_si_b$DE == "Down")
sum(dedup_si_b$txType == "intron" & dedup_si_b$DE == "Down")
sum(dedup_si_b$txType == "exon" & dedup_si_b$DE == "Down")

#plot
mcols(dedup_si_b) %>% as_tibble() %>% 
  ggplot(aes(y = DE, fill = txType ) ) +
  geom_bar(position = "fill") + # this is a very elegant solution :-)
  labs(y = 'Differential expression', x = 'Fraction of TSS') + 
  theme(aspect.ratio = 0.5)
for (ext in c(".png", ".pdf")) {
  paste0("plot", ext) %>% ggsave(width = 7, height = 7, units = "in")
}

