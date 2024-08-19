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


path <- "."  # change to the directory containing merged TSS-seq Bedgraph files obtained from 02_TSS-seq processing.sh and 03_STRIPE-seq processing.sh


scripts <- c("batchReadTrackData", "normalizeGR", "saveGRangesAsBedGraph", "expandGRtoUnitWidth")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Kindgren_et_al_2019/blob/master/", script, ".R?raw=TRUE") %>% devtools::source_url()
}
devtools::source_url("https://github.com/Maxim-Ivanov/TranscriptomeReconstructoR/blob/main/R/Utility_functions.R?raw=TRUE") # "master" -> "main"
devtools::source_url("https://github.com/Maxim-Ivanov/Nielsen_et_al_2018/blob/master/Custom_assignTxType_function_for_CAGEfightR.R?raw=TRUE"
devtools::source_url("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/find_matched_control.R?raw=TRUE")

triple_fw<- list.files(".", pattern = "fw.bg.gz") # change to the directory containing merged TSS-seq Bedgraph files obtained from 02_TSS-seq processing.sh and 03_STRIPE-seq processing.sh
triple_rev<- list.files(".", pattern = "rev.bg.gz") # change to the directory containing merged TSS-seq Bedgraph files obtained from 02_TSS-seq processing.sh and 03_STRIPE-seq processing.sh

batched_fwd <- batchReadTrackData(triple_fw, dir = ".", format = 'bedGraph', strand = "+", seqinfo = seqinfo(txdb)) # change to the directory containing merged TSS-seq Bedgraph files obtained from 02_TSS-seq processing.sh and 03_STRIPE-seq processing.sh
batched_rev <- batchReadTrackData(triple_rev, dir = ".", format = 'bedGraph', strand = "-", seqinfo = seqinfo(txdb)) # change to the directory containing merged TSS-seq Bedgraph files obtained from 02_TSS-seq processing.sh and 03_STRIPE-seq processing.sh

names(batched_fwd) <- names(batched_fwd) %>% str_replace("_R1_001_trimmed_UMI_Aligned.sortedByCoord.out_mapq10_dedup_fw.bg.gz", "")
names(batched_rev) <- names(batched_rev) %>% str_replace("_R1_001_trimmed_UMI_Aligned.sortedByCoord.out_mapq10_dedup_rev.bg.gz", "")

bg_data <- mapply(function(x, y) { sort(c(x, y)) }, batched_fwd, batched_rev, SIMPLIFY = FALSE)
bg_data_merged <- merge_GRanges(bg_data[1:3])
bg_data_merged_norm <- normalizeGR(bg_data_merged)

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
sample_names <- bw_plus_filenames %>% str_replace("_unitWidth_Plus.bw", "") %>% str_replace("STRIPE-Seq_", "") 
names(bw_plus) <- sample_names
names(bw_minus) <- sample_names

# For DE calling, the "Genotype" column in design_matrix must contain both conditions to compare ("wt" vs "mut"):
design_matrix <- data.frame("Name" = sample_names, 
                            "BigWigPlus" = bw_plus_filenames, 
                            "BigWigMinus" = bw_minus_filenames, 
                            "Genotype" = c("wt") %>% rep(each = 3) %>% as.factor() %>% relevel(ref = "wt"), 
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

# Add decisions on the sample specificity of called TSS:
mcols(tss_gr)$nz_wt <- rowSums(assay(tss)[,c(1:3)]>1) 

# Calculate cpm/tpm values
cm <- as.data.frame(cpm(assay(tss)))
cm$key <- rownames(cm)
orig <- data.frame("key"=names(rowRanges(tss)))
x <- left_join(orig, cm, by=c("key"))
# Add tpm results to the RSE object:
mcols(tss_gr) <- cbind(mcols(tss_gr), x[-1])

# Calculate mean tpm of samples
mcols(tss_gr) <- cbind(mcols(tss_gr), wt_mean = rowMeans(cm[,c(1:3)]))

tss_gr <- tss_gr[mcols(tss_gr)$txType == "promoter", "fiveUTR", "intron", "exon"]
tss_gr <- tss_gr[mcols(tss_gr)$nz_wt > 1]


### Read TI candidates and upstream transcripts 
goi <- read.csv(".csv") # change to the directory containing .csv file obtained from 05_upstream transcripts.r
list_gr <- makeGRangesFromDataFrame(goi, keep.extra.columns = TRUE)

## select TI candidates
gene_list_ti <- list_gr$ti_gene
filtered_granges_ti <- tss_gr[mcols(tss_gr)$geneID %in% gene_list_ti]

# deduplicate
dedup_TI <- filtered_granges_ti[mcols(filtered_granges_ti)$txType == "promoter", "fiveUTR", "intron", "exon"]
dedup_TI_a <- dedup_TI[order(dedup_TI$geneID, -abs(dedup_TI$wt_mean) ), ] #sort by id and reverse of abs(value)
dedup_TI_a <- dedup_TI_a[!duplicated(dedup_TI_a$geneID), ]

write.csv(dedup_TI_a, "dedup_TI_txdb.csv")

goi_ti_tss <- read.csv("./TI_TSS_txdb.csv") # change to the directory containing .csv file obtained from dedup_TI_a
ti_tss_list_gr <- makeGRangesFromDataFrame(goi_ti_tss, keep.extra.columns = TRUE)

## select upstream transcripts
gene_list_up <- up_list_gr$gene_id
filtered_granges_up <- tss_gr[mcols(tss_gr)$geneID %in% gene_list_up]

# deduplicate
dedup_up <- filtered_granges_up[mcols(filtered_granges_up)$txType == "promoter"]
dedup_up_a <- dedup_up[order(dedup_up$geneID, -abs(dedup_up$wt_mean) ), ] #sort by id and reverse of abs(value)
dedup_up_a <- dedup_up_a[!duplicated(dedup_up_a$geneID), ]

write.csv(dedup_up_a, "UP_TSS_txdb.csv")

goi_up_tss <- read.csv("./UP_TSS_txdb.csv") # change to the directory containing .csv file obtained from dedup_up_a
up_tss_list_gr <- makeGRangesFromDataFrame(goi_up_tss, keep.extra.columns = TRUE)

# Violin
devtools::source_url("https://github.com/sa-lee/plyranges/blob/master/R/ranges-bind.R?raw=TRUE")

two <-  bind_ranges(ti = ti_tss_list_gr, up = up_tss_list_gr, .id = "Mutant")
ttl <- "Violin"

two1 <- two[mcols(two)$wt_mean < 50]
two1 %>% mcols() %>% as_tibble() %>%
  #two1 %>% mcols() %>% as_tibble() %>%
  ggplot(aes(x = Mutant, y = wt_mean, fill = Mutant)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.2, outlier.colour = NA, fill = "white") +
  ggtitle(ttl)
ggsave(paste0(ttl, "300.png"), width = 7, height = 7, units = "in")

