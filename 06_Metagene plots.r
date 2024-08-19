library(rtracklayer)
library(tidyverse)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28

scripts <- c("batch_read_track_data", "merge_and_normalize_GRanges", "metagene_matrix", "draw_metagene_plot", "get_overlapping_scores")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}


path <- "." # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh

#Load CUT&Tag data (merged biological replicates):

###seedlings###
###H3K4me1###
H3K4me1.ab <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K4me1.1F4 <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K4me1.epi <- import('I:.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh

H3K4me1.ab <- dropSeqlevels(H3K4me1.ab,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K4me1.1F4 <- dropSeqlevels(H3K4me1.1F4,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K4me1.epi <- dropSeqlevels(H3K4me1.epi,c("ChrM","ChrC"),pruning.mode ="coarse")

H3K4me1.ab <- renameSeqlevels(H3K4me1.ab, gsub("Chr","", seqlevels(H3K4me1.ab)))
H3K4me1.1F4 <- renameSeqlevels(H3K4me1.1F4, gsub("Chr","", seqlevels(H3K4me1.1F4)))
H3K4me1.epi <- renameSeqlevels(H3K4me1.epi, gsub("Chr","", seqlevels(H3K4me1.epi)))

seqinfo(H3K4me1.ab) <- seqinfo(txdb)
seqinfo(H3K4me1.1F4) <- seqinfo(txdb)
seqinfo(H3K4me1.epi) <- seqinfo(txdb)

H3K4me1.ab <- trim(H3K4me1.ab)
H3K4me1.1F4 <- trim(H3K4me1.1F4)
H3K4me1.epi <- trim(H3K4me1.epi)

H3K4me1 <- list(H3K4me1.ab, H3K4me1.1F4, H3K4me1.epi)

###H3K4me2###
H3K4me2.ab <- 
import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K4me2.ab <- dropSeqlevels(H3K4me2.ab,c("ChrM","ChrC"),pruning.mode ="coarse")

H3K4me2.ab <- renameSeqlevels(H3K4me2.ab, gsub("Chr","", seqlevels(H3K4me2.ab)))

seqinfo(H3K4me2.ab) <- seqinfo(txdb)

H3K4me2.ab <- trim(H3K4me2.ab)
H3K4me2 <- list(H3K4me2.ab)

###H3K4me3###
H3K4me3.ab <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K4me3.epi <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh

H3K4me3.ab <- dropSeqlevels(H3K4me3.ab,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K4me3.epi <- dropSeqlevels(H3K4me3.epi,c("ChrM","ChrC"),pruning.mode ="coarse")

H3K4me3.ab <- renameSeqlevels(H3K4me3.ab, gsub("Chr","", seqlevels(H3K4me3.ab)))
H3K4me3.epi <- renameSeqlevels(H3K4me3.epi, gsub("Chr","", seqlevels(H3K4me3.epi)))

seqinfo(H3K4me3.ab) <- seqinfo(txdb)
seqinfo(H3K4me3.epi) <- seqinfo(txdb)

H3K4me3.ab <- trim(H3K4me3.ab)
H3K4me3.epi <- trim(H3K4me3.epi)

H3K4me3 <- list(H3K4me3.ab, H3K4me3.epi)

###H3K36me1###
H3K36me1.ab <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me1.bo <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh

H3K36me1.ab <- dropSeqlevels(H3K36me1.ab,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me1.bo <- dropSeqlevels(H3K36me1.bo,c("ChrM","ChrC"),pruning.mode ="coarse")

H3K36me1.ab <- renameSeqlevels(H3K36me1.ab, gsub("Chr","", seqlevels(H3K36me1.ab)))
H3K36me1.bo <- renameSeqlevels(H3K36me1.bo, gsub("Chr","", seqlevels(H3K36me1.bo)))

seqinfo(H3K36me1.ab) <- seqinfo(txdb)
seqinfo(H3K36me1.bo) <- seqinfo(txdb)

H3K36me1.ab <- trim(H3K36me1.ab)
H3K36me1.bo <- trim(H3K36me1.bo)

H3K36me1 <- list(H3K36me1.ab, H3K36me1.bo)

###H3K36me2###
H3K36me2.ab <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me2.epi <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me2.am <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh

H3K36me2.ab <- dropSeqlevels(H3K36me2.ab,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me2.epi <- dropSeqlevels(H3K36me2.epi,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me2.am <- dropSeqlevels(H3K36me2.am,c("ChrM","ChrC"),pruning.mode ="coarse")

H3K36me2.ab <- renameSeqlevels(H3K36me2.ab, gsub("Chr","", seqlevels(H3K36me2.ab)))
H3K36me2.epi <- renameSeqlevels(H3K36me2.epi, gsub("Chr","", seqlevels(H3K36me2.epi)))
H3K36me2.am <- renameSeqlevels(H3K36me2.am, gsub("Chr","", seqlevels(H3K36me2.am)))

seqinfo(H3K36me2.ab) <- seqinfo(txdb)
seqinfo(H3K36me2.epi) <- seqinfo(txdb)
seqinfo(H3K36me2.am) <- seqinfo(txdb)

H3K36me2.ab <- trim(H3K36me2.ab)
H3K36me2.epi <- trim(H3K36me2.epi)
H3K36me2.am <- trim(H3K36me2.am)

H3K36me2 <- list(H3K36me2.ab, H3K36me2.epi, H3K36me2.am)

###H3K36me3###
H3K36me3.ab <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me3.epi1E2 <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me3.epi1B9 <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me3.ac <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh

H3K36me3.ab <- dropSeqlevels(H3K36me3.ab,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me3.epi1E2 <- dropSeqlevels(H3K36me3.epi1E2,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me3.epi1B9 <- dropSeqlevels(H3K36me3.epi1B9,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me3.ac <- dropSeqlevels(H3K36me3.ac,c("ChrM","ChrC"),pruning.mode ="coarse")

H3K36me3.ab <- renameSeqlevels(H3K36me3.ab, gsub("Chr","", seqlevels(H3K36me3.ab)))
H3K36me3.epi1E2 <- renameSeqlevels(H3K36me3.epi1E2, gsub("Chr","", seqlevels(H3K36me3.epi1E2)))
H3K36me3.epi1B9 <- renameSeqlevels(H3K36me3.epi1B9, gsub("Chr","", seqlevels(H3K36me3.epi1B9)))
H3K36me3.ac <- renameSeqlevels(H3K36me3.ac, gsub("Chr","", seqlevels(H3K36me3.ac)))

seqinfo(H3K36me3.ab) <- seqinfo(txdb)
seqinfo(H3K36me3.epi1E2) <- seqinfo(txdb)
seqinfo(H3K36me3.epi1B9) <- seqinfo(txdb)
seqinfo(H3K36me3.ac) <- seqinfo(txdb)

H3K36me3.ab <- trim(H3K36me3.ab)
H3K36me3.epi1E2 <- trim(H3K36me3.epi1E2)
H3K36me3.epi1B9 <- trim(H3K36me3.epi1B9)
H3K36me3.ac <- trim(H3K36me3.ac)

H3K36me3 <- list(H3K36me3.ab, H3K36me3.epi1E2, H3K36me3.epi1B9, H3K36me3.ac)

###H3K27me3###
H3K27me3.milli <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K27me3.tf <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K27me3.epi <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh

H3K27me3.milli <- dropSeqlevels(H3K27me3.milli,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K27me3.tf <- dropSeqlevels(H3K27me3.tf,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K27me3.epi <- dropSeqlevels(H3K27me3.epi,c("ChrM","ChrC"),pruning.mode ="coarse")

H3K27me3.milli <- renameSeqlevels(H3K27me3.milli, gsub("Chr","", seqlevels(H3K27me3.milli)))
H3K27me3.tf <- renameSeqlevels(H3K27me3.tf, gsub("Chr","", seqlevels(H3K27me3.tf)))
H3K27me3.epi <- renameSeqlevels(H3K27me3.epi, gsub("Chr","", seqlevels(H3K27me3.epi)))

seqinfo(H3K27me3.milli) <- seqinfo(txdb)
seqinfo(H3K27me3.tf) <- seqinfo(txdb)
seqinfo(H3K27me3.epi) <- seqinfo(txdb)

H3K27me3.milli <- trim(H3K27me3.milli)
H3K27me3.tf <- trim(H3K27me3.tf)
H3K27me3.epi <- trim(H3K27me3.epi)

H3K27me3 <- list(H3K27me3.milli, H3K27me3.tf, H3K27me3.epi)

###siliques###
H3K4me1_si <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K4me3_si <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me1_si <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me2_si <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me3_si <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh


H3K4me1_si <- dropSeqlevels(H3K4me1_si,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K4me3_si <- dropSeqlevels(H3K4me3_si,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me1_si <- dropSeqlevels(H3K36me1_si,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me2_si <- dropSeqlevels(H3K36me2_si,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me3_si <- dropSeqlevels(H3K36me3_si,c("ChrM","ChrC"),pruning.mode ="coarse")


H3K4me1_si <- renameSeqlevels(H3K4me1_si, gsub("Chr","", seqlevels(H3K4me1_si)))
H3K4me3_si <- renameSeqlevels(H3K4me3_si, gsub("Chr","", seqlevels(H3K4me3_si)))
H3K36me1_si <- renameSeqlevels(H3K36me1_si, gsub("Chr", "", seqlevels(H3K36me1_si)))
H3K36me2_si <- renameSeqlevels(H3K36me2_si, gsub("Chr", "", seqlevels(H3K36me2_si)))
H3K36me3_si <- renameSeqlevels(H3K36me3_si, gsub("Chr", "", seqlevels(H3K36me3_si)))

seqinfo(H3K4me1_si) <- seqinfo(txdb)
seqinfo(H3K4me3_si) <- seqinfo(txdb)
seqinfo(H3K36me1_si) <- seqinfo(txdb)
seqinfo(H3K36me2_si) <- seqinfo(txdb)
seqinfo(H3K36me3_si) <- seqinfo(txdb)

H3K4me1_si <- trim(H3K4me1_si)
H3K4me3_si <- trim(H3K4me3_si)
H3K36me1_si <- trim(H3K36me1_si)
H3K36me2_si <- trim(H3K36me2_si)
H3K36me3_si <- trim(H3K36me3_si)


### Make lists to be plotted
H3K4me1 <- list(H3K4me1.ab, H3K4me1.1F4, H3K4me1.epi)
H3K4me2 <- list(H3K4me2.ab)
H3K4me3 <- list(H3K4me3.ab, H3K4me3.epi)
H3K36me1 <- list(H3K36me1.ab, H3K36me1.bo)
H3K36me2 <- list(H3K36me2.ab, H3K36me2.epi, H3K36me2.am)
H3K36me3 <- list(H3K36me3.ab, H3K36me3.epi1E2, H3K36me3.epi1B9, H3K36me3.ac)
H3K27me3 <- list(H3K27me3.milli, H3K27me3.tf, H3K27me3.epi)
siliques <- list(H3K36me1_si,H3K36me2_si,H3K36me3_si,H3K4me1_si,H3K4me2_si,H3K4me3_si)

### make genomics intervals ###
win_genes <- genes(txdb, columns = c("gene_id", "tx_type"))
mcols(win_genes)$tx_type <- mcols(win_genes)$tx_type %>% unlist(use.names = FALSE)

# Filter for nuclear protein-coding genes:
win_genes <- win_genes[mcols(win_genes)$tx_type == "protein_coding" & seqnames(win_genes) %in% 1:5]

# Filter by width:
win_genes <- win_genes[width(win_genes) >= 500 & width(win_genes) <= 5000]
win_minus1kb <- flank(win_genes, 1000) %>% trim()
win_plus1kb <- flank(win_genes, 1000, start = FALSE) %>% trim()


### Load and normalize plaNET-seq data:
path_plaNET <- "." # change to the directory containing plaNET-seq Bedgraph file (GSE131733)
wt_file <- list.files(path_plaNET, pattern = "WT_merged_norm1M.bedgraph.gz")
wt_plaNET <- batch_read_track_data(wt_file, dir = path_plaNET, format = 'bedGraph', seqinfo = seqinfo(txdb))
names(wt_plaNET) <- names(wt_plaNET) %>% str_replace('_merged_norm1M.bedgraph.gz', '')

# Count plaNET-seq signal on genes which were chosen for plotting metagenes:
rpm <- win_genes %>% get_overlapping_scores(wt_plaNET, value = "count_matrix") %>% as.numeric() 

# Normalize RPM by the gene width (in Kb):
fpkm <- rpm / width(win_genes) * 1000

# Stratify genes by FPKM values:
grp <- ifelse(fpkm < 1, "Low", ifelse(fpkm < 13, "Medium", "High")) %>% factor(levels = c("Low", "Medium", "High"))

# Add grp as mcols to win_genes:
mcols(win_genes)$grp <- grp
mcols(win_minus1kb)$grp <- grp
mcols(win_plus1kb)$grp <- grp


### Draw metagene plot for medium expressed genes
profile_medium1 <- lapply(#list, metagene_matrix, intervals = win_minus1kb[mcols(win_minus1kb)$grp == "Medium"], scaling = TRUE, skip.zeros = FALSE, skip.outliers = FALSE, matrix.length = 100)
profile_medium2 <- lapply(#list, metagene_matrix, intervals = win_genes[mcols(win_genes)$grp == "Medium"], scaling = TRUE, matrix.length = 500, skip.zeros = FALSE, skip.outliers = FALSE)
profile_medium3 <- lapply(#list, metagene_matrix, intervals = win_plus1kb[mcols(win_plus1kb)$grp == "Medium"], scaling = TRUE, skip.zeros = FALSE, skip.outliers = FALSE, matrix.length = 100)
profile_medium <- mapply(function(x, y, z) { cbind(x, y, z) }, profile_medium1, profile_medium2, profile_medium3, SIMPLIFY = FALSE)
names(profile_medium) <- c(".", ".", ".")
y_limits <- c(0, 4)
draw_metagene_plot(profile_medium, title = ".", x.axis = seq(-99, 600), vline = c(0, 500), width = 14, height = 6, units = "in", ylim = y_limits) 

### Determine highest peak position
avg_h3k4me1_abcam <- apply(profile_h3k4me1$H3K4me1.ab, 2, mean, na.rm = TRUE)
max(avg_h3k4me1_abcam)

avg_h3k4me1_epi <- apply(profile_h3k4me1$H3K4me1.1F4, 2, mean, na.rm = TRUE)
max(avg_h3k4me1_epi) 

avg_h3k4me1_new1F4 <- apply(profile_h3k4me1$H3K4me1.epi, 2, mean, na.rm = TRUE)
max(avg_h3k4me1_new1F4)

H3K4me1_all <-  c(H3K4me1_abcam = avg_h3k4me1_abcam, H3K4me1_epi = avg_h3k4me1_epi, H3K4me1_1F4 = avg_h3k4me1_new1F4)

write.csv(H3K4me1_all,"Metagene_avg_h3k4me1_all.csv")

### Make genomic intervals for selected genes
goi <- read.csv(".") # change to the directory containing .csv with genes of interest
goi_gr <- makeGRangesFromDataFrame(goi, keep.extra.columns=FALSE, ignore.strand=FALSE, seqinfo=NULL, seqnames.field="seqnames", start.field="start", end.field="end", strand.field="strand")
goi_txdb <- transcriptsByOverlaps(txdb, goi_gr, type="any")

win_genes_goi <- goi_txdb
mcols(win_genes_goi)$tx_name <- mcols(win_genes_goi)$tx_name %>% unlist(use.names = FALSE)

### Draw metagene plot for selected genes
win_genes_goi <- win_genes_goi[width(win_genes_goi) >= 500 & width(win_genes_goi) <= 5000]
win_minus1kb_goi <- flank(win_genes_goi, 1000) %>% trim()
win_plus1kb_goi <- flank(win_genes_goi, 1000, start = FALSE) %>% trim()
profile1_goi <- lapply(#list, metagene_matrix, intervals = win_minus1kb_goi, scaling = TRUE, skip.zeros = FALSE, skip.outliers = FALSE, matrix.length = 100)
profile2_goi <- lapply(#list, metagene_matrix, intervals = win_genes_goi, scaling = TRUE, matrix.length = 500, skip.zeros = FALSE, skip.outliers = FALSE)
profile3_goi <- lapply(#list, metagene_matrix, intervals = win_plus1kb_goi, scaling = TRUE, skip.zeros = FALSE, skip.outliers = FALSE, matrix.length = 100)
profile_goi <- mapply(function(x, y, z) { cbind(x, y, z) }, profile1_goi, profile2_goi, profile3_goi, SIMPLIFY = FALSE)
names(profile_goi) <- c(".", ".", ".", ".", ".")
y_limits <- c(0, 4)
draw_metagene_plot(profile_goi, title = ".", x.axis = seq(-99, 600), vline = c(0, 500), width = 14, height = 6, units = "in", ylim = y_limits) #ther

### Draw metagene plot for PAS and 1kb around for selected genes
win_genes_pas <- win_genes_goi %>% resize(0, "end") %>% resize(1000, "center") %>% trim() 
pas_profile_goi <- lapply("list, metagene_matrix, intervals = win_genes_pas, scaling = TRUE, skip.zeros = FALSE, skip.outliers = FALSE, matrix.length = 600)
names(pas_profile_goi) <- c(".", ".", ".", .", ".", ".", ".")
draw_metagene_plot(pas_profile_goi, title = ".", x.axis = seq(-99, 500), vline = 0, width = 8, height = 6, units = "in")

### Load and normalize TSS-seq data:
path_TSS <- "." # change to the directory containing merged TSS-seq Bedgraph files obtained from 02_TSS-seq processing.sh
wt_file_tss <- list.files(path_TSS, pattern = "siliques-IGV.bedgraph")
wt_TSS <- batch_read_track_data(wt_file_tss, dir = path_TSS, format = 'bedGraph', seqinfo = seqinfo(txdb))
names(wt_TSS) <- names(wt_TSS) %>% str_replace('-IGV.bedgraph', '')


# Count TSS-seq signal on genes which were chosen for plotting metagenes:
rpm_tss <- win_genes %>% get_overlapping_scores(wt_TSS, value = "count_matrix") %>% as.numeric() 

# Normalize RPM by the gene width (in Kb):
fpkm_tss <- rpm_tss / width(win_genes) * 1000

# Stratify genes by FPKM values:
grp_tss <- ifelse(fpkm_tss < 0.08, "Low", ifelse(fpkm_tss < 7.2, "Medium", "High")) %>% factor(levels = c("Low", "Medium", "High"))

# Add grp as mcols to win_genes:
mcols(win_genes_tss)$grp_tss <- grp_tss
mcols(win_minus1kb_tss)$grp_tss <- grp_tss
mcols(win_plus1kb_tss)$grp_tss <- grp_tss


### Draw metagene plot for medium expressed genes based on TSS-seq in siliques
# For genes with medium level of expression
profile_siliques1 <- lapply(#list, metagene_matrix, intervals = win_minus1kb[mcols(win_minus1kb)$grp == "Medium"], scaling = TRUE, skip.zeros = FALSE, skip.outliers = FALSE, matrix.length = 100)
profile_siliques2 <- lapply(#list, metagene_matrix, intervals = win_genes[mcols(win_genes)$grp == "Medium"], scaling = TRUE, matrix.length = 500, skip.zeros = FALSE, skip.outliers = FALSE)
profile_siliques3 <- lapply(#list, metagene_matrix, intervals = win_plus1kb[mcols(win_plus1kb)$grp == "Medium"], scaling = TRUE, skip.zeros = FALSE, skip.outliers = FALSE, matrix.length = 100)
profile_siliques <- mapply(function(x, y, z) { cbind(x, y, z) }, profile_siliques1, profile_siliques2, profile_siliques3, SIMPLIFY = FALSE)
names(profile_siliques) <- c(".", ".", ".", ".", ".", ".")
y_limits <- c(0, 4)
draw_metagene_plot(profile_siliques, title = ".", x.axis = seq(-99, 600), vline = c(0, 500), width = 14, height = 6, units = "in", ylim = y_limits)