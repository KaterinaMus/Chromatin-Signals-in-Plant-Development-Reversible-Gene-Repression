library(rtracklayer)
library(tidyverse)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28


win_genes <- genes(txdb, columns = c("gene_id", "tx_type"))
mcols(win_genes)$tx_type <- mcols(win_genes)$tx_type %>% unlist(use.names = FALSE)

### read list of TI candidates obtained from 04_TI ranking
goi_TI <- read.csv(".csv") # change to the directory containing .csv with TI candidates obtained from 04_TI candidate ranking.r
TI_list_gr <- makeGRangesFromDataFrame(goi_TI, keep.extra.columns = TRUE)
gene_list_TI <- TI_list_gr$gene_id
filtered_granges_TI <- win_genes[mcols(win_genes)$gene_id %in% gene_list_TI]

### read annotation with RT tails from Ivanov, M., Sandelin, A., and Marquardt, S. (2021). TrancriptomeReconstructoR: data-driven annotation of complex transcriptomes. BMC Bioinformatics 22, 290.
all_withRTTails <- read.csv(".csv") # change to the directory containing annotation with RT tails (Ivanov, M., et al. 2021)
all_withRTTails_gr <- makeGRangesFromDataFrame(all_withRTTails, keep.extra.columns = TRUE)

### Upstream genes of candidates #TRUE - ignoring strand FALSE - on the same strand
up <- follow(filtered_granges_TI, all_withRTTails_gr, ignore.strand=FALSE)
candidates_upstream <- all_withRTTails_gr[up]

### Distance to the upstream gene
pos_strand <- start(filtered_granges_TI) - end(candidates_upstream) 
neg_strand <- start(candidates_upstream) - end(filtered_granges_TI) 
distance <- pmax(pos_strand, neg_strand)

candidates_upstream$ti_gene <- filtered_granges_TI$gene_id
candidates_upstream$distance <- distance

write.csv(candidates_upstream_all,"upstream_genes_all.csv")

##selection based on median distance <1200bp

candidates_upstream_median <- candidates_upstream[candidates_upstream$distance<1200]

write.csv(candidates_upstream_median,"upstream_genes.csv")

