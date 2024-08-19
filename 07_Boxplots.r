# Load the required libraries:
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(SummarizedExperiment)
library(rtracklayer)
library(ggplot2)
library(ggpubr)
library(tidyverse)


scripts <- c("batch_read_track_data", "merge_and_normalize_GRanges", "metagene_matrix", "draw_metagene_plot", "get_overlapping_scores")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}


###CUT&Tag seedlings
path <- "." # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K4me1 <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K4me3 <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K27me3 <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me1 <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me2 <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me3 <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh

#process them
H3K4me1 <- dropSeqlevels(H3K4me1,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K4me3 <- dropSeqlevels(H3K4me3,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K27me3 <- dropSeqlevels(H3K27me3,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me1 <- dropSeqlevels(H3K36me1,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me2 <- dropSeqlevels(H3K36me2,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me3 <- dropSeqlevels(H3K36me3,c("ChrM","ChrC"),pruning.mode ="coarse")


###to remove Chr1 in seqnames, otherwise it is not compatible with txdb
H3K4me1 <- renameSeqlevels(H3K4me1, gsub("Chr","", seqlevels(H3K4me1)))
H3K4me3 <- renameSeqlevels(H3K4me3, gsub("Chr", "", seqlevels(H3K4me3)))
H3K27me3 <- renameSeqlevels(H3K27me3, gsub("Chr", "", seqlevels(H3K27me3)))
H3K36me1 <- renameSeqlevels(H3K36me1, gsub("Chr", "", seqlevels(H3K36me1)))
H3K36me2 <- renameSeqlevels(H3K36me2, gsub("Chr", "", seqlevels(H3K36me2)))
H3K36me3 <- renameSeqlevels(H3K36me3, gsub("Chr", "", seqlevels(H3K36me3)))

###CUT&Tag siliques
path2 <- "." # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh

H3K4me1_si <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K4me3_si <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me1_si <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me2_si <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me3_si <- import('.',format="bedgraph") # change to the directory containing merged CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh

###process them
H3K4me1_si <- dropSeqlevels(H3K4me1_si,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K4me3_si <- dropSeqlevels(H3K4me3_si,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me1_si <- dropSeqlevels(H3K36me1_si,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me2_si <- dropSeqlevels(H3K36me2_si,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me3_si <- dropSeqlevels(H3K36me3_si,c("ChrM","ChrC"),pruning.mode ="coarse")

###to remove Chr1 in seqnames, otherwise it is not compatible with txdb
H3K4me1_si <- renameSeqlevels(H3K4me1_si, gsub("Chr","", seqlevels(H3K4me1_si)))
H3K4me3_si <- renameSeqlevels(H3K4me3_si, gsub("Chr","", seqlevels(H3K4me3_si)))
H3K36me1_si <- renameSeqlevels(H3K36me1_si, gsub("Chr", "", seqlevels(H3K36me1_si)))
H3K36me2_si <- renameSeqlevels(H3K36me2_si, gsub("Chr", "", seqlevels(H3K36me2_si)))
H3K36me3_si <- renameSeqlevels(H3K36me3_si, gsub("Chr", "", seqlevels(H3K36me3_si)))


## extend function
extend <- function(x, upstream=0, downstream=0)     
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

### Reading selected genes
tss_box <- read.csv(".", sep = ";", dec = ".", header=TRUE) # change to the directory containing .csv with TSS start-end from 08_TSS data validation.r
tss_box_gr <- makeGRangesFromDataFrame(tss_box, keep.extra.columns = TRUE)

# make a window for signal selection 300 bp upstream of TSS
tss_300 <- extend(tss_box_gr,upstream=275,downstream=0)
tss_3 <- extend(tss_300,upstream=50,downstream=0)

#50bp window for overlap
tss_only_3 <- psetdiff(tss_3,tss_300)

# make a window for signal selection 100 bp downstream of TSS
tss_100 <- extend(tss_box_gr,upstream=0,downstream=50)
tss_1 <- extend(tss_100,upstream=0,downstream=50)

#50bp window for overlap
tss_only_1 <- psetdiff(tss_1,tss_100)


### make a window for signal selection 300 bp downstream of PAS for upstream transcripts
goi_up <- read.csv(".") # change to the directory containing .csv with upstream transcripts obtained from 05_upstream transcripts.r
goi_up_gr <- makeGRangesFromDataFrame(goi_up, keep.extra.columns=FALSE, ignore.strand=FALSE, seqinfo=NULL, seqnames.field="seqnames", start.field="start", end.field="end", strand.field="strand")

pas_up_300 <- extend(goi_up_gr,upstream=0,downstream=275)
pas_up <- extend(pas_up_300,upstream=0,downstream=50)

###50bp window for overlap
pas_up_only <- psetdiff(pas_up,pas_up_300)

### make a window for signal selection 300 bp downstream of PAS for 1000 medium expressed genes from 04_Metagene plots.r
goi_medium <- read.csv(".") # change to the directory containing .csv with genes of interest
goi_medium_gr <- makeGRangesFromDataFrame(goi_medium, keep.extra.columns=FALSE, ignore.strand=FALSE, seqinfo=NULL, seqnames.field="seqnames", start.field="start", end.field="end", strand.field="strand")

pas_medium_300 <- extend(goi_medium_gr,upstream=0,downstream=275)
pas_medium <- extend(pas_medium_300,upstream=0,downstream=50)

###50bp window for overlap
pas_medium_only <- psetdiff(pas_medium,pas_medium_300)

### Count the overlaps for a window 300 bp upstream of TSS 
###H3K4me1
se_H3K4me1_up <- countOverlaps(tss_only_3,H3K4me1,type="any")
se_K4me1_up <- tss_3
se_K4me1_up$value <- se_H3K4me1_up

si_H3K4me1_up <- countOverlaps(tss_only_3,H3K4me1_si,type="any")
si_K4me1_up <- tss_3
si_K4me1_up$value <- si_H3K4me1_up

###H3K4me3
se_H3K4me3_up <- countOverlaps(tss_only_3,H3K4me3,type="any")
se_K4me3_up <- tss_3
se_K4me3_up$value <- se_H3K4me3_up

si_H3K4me3_up <- countOverlaps(tss_only_3,H3K4me3_si,type="any")
si_K4me3_up <- tss_3
si_K4me3_up$value <- si_H3K4me3_up

###H3K36me1
se_H3K36me1_up <- countOverlaps(tss_only_3,H3K36me1,type="any")
se_K36me1_up <- tss_3
se_K36me1_up$value <- se_H3K36me1_up

si_H3K36me1_up <- countOverlaps(tss_only_3,H3K36me1_si,type="any")
si_K36me1_up <- tss_3
si_K36me1_up$value <- si_H3K36me1_up


###H3K36me2
se_H3K36me2_up <- countOverlaps(tss_only_3,H3K36me2,type="any")
se_K36me2_up <- tss_3
se_K36me2_up$value <- se_H3K36me2_up

si_H3K36me2_up <- countOverlaps(tss_only_3,H3K36me2_si,type="any")
si_K36me2_up <- tss_3
si_K36me2_up$value <- si_H3K36me2_up

###H3K36me3
se_H3K36me3_up <- countOverlaps(tss_only_3,H3K36me3,type="any")
se_K36me3_up <- tss_3
se_K36me3_up$value <- se_H3K36me3_up

si_H3K36me3_up <- countOverlaps(tss_only_3,H3K36me3_si,type="any")
si_K36me3_up <- tss_3
si_K36me3_up$value <- si_H3K36me3_up

### draw a boxplot
devtools::source_url("https://github.com/sa-lee/plyranges/blob/master/R/ranges-bind.R?raw=TRUE")
two <-  bind_ranges(se_K4me1 = se_K4me1_up, si_K4me1 = si_K4me1_up, se_K4me3 = se_K4me3_up, si_K4me3 = si_K4me3_up, 
                    se_K36me1 = se_K36me1_up, si_K36me1 = si_K36me1_up, 
                    se_K36me2 = se_K36me2_up, si_K36me2 = si_K36me2_up, 
                    se_K36me3 = se_K36me3_up, si_K36me3 = si_K36me3_up, .id = "Mutant")
ttl <- "boxplot"

two <- two[mcols(two)$value < 30]
two%>% mcols() %>% as_tibble() %>%
  ggplot(aes(x = factor (Mutant, level=c('se_K4me3', 'si_K4me3', 'se_K36me3', 'si_K36me3', 'se_K4me1', 'si_K4me1', 'se_K36me2', 'si_K36me2', 'se_K36me1', 'si_K36me1')), y = value, fill = Mutant)) +
  geom_boxplot(width = 0.5, outlier.colour = NA, fill = "white") +
  #ylim (0, 7)
  ggtitle(ttl)
ggsave(paste0(ttl, ".pdf"), width = 7, height = 7, units = "in")

### Count the overlaps for a window 100 bp downstream 
###H3K4me1
se_H3K4me1_down <- countOverlaps(tss_only_3,H3K4me1,type="any")
se_K4me1_down <- tss_3
se_K4me1_down$value <- se_H3K4me1_down

si_H3K4me1_down <- countOverlaps(tss_only_3,H3K4me1_si,type="any")
si_K4me1_down <- tss_3
si_K4me1_down$value <- si_H3K4me1_down

###H3K4me3
se_H3K4me3_down <- countOverlaps(tss_only_3,H3K4me3,type="any")
se_K4me3_down <- tss_3
se_K4me3_down$value <- se_H3K4me3_down

si_H3K4me3_down <- countOverlaps(tss_only_3,H3K4me3_si,type="any")
si_K4me3_down <- tss_3
si_K4me3_down$value <- si_H3K4me3_down

###H3K36me1
se_H3K36me1_down <- countOverlaps(tss_only_3,H3K36me1,type="any")
se_K36me1_down <- tss_3
se_K36me1_down$value <- se_H3K36me1_down

si_H3K36me1_down <- countOverlaps(tss_only_3,H3K36me1_si,type="any")
si_K36me1_down <- tss_3
si_K36me1_down$value <- si_H3K36me1_down


###H3K36me2
se_H3K36me2_down <- countOverlaps(tss_only_3,H3K36me2,type="any")
se_K36me2_down <- tss_3
se_K36me2_down$value <- se_H3K36me2_down

si_H3K36me2_down <- countOverlaps(tss_only_3,H3K36me2_si,type="any")
si_K36me2_down <- tss_3
si_K36me2_down$value <- si_H3K36me2_down

###H3K36me3
se_H3K36me3_down <- countOverlaps(tss_only_3,H3K36me3,type="any")
se_K36me3_down <- tss_3
se_K36me3_down$value <- se_H3K36me3_down

si_H3K36me3_down <- countOverlaps(tss_only_3,H3K36me3_si,type="any")
si_K36me3_down <- tss_3
si_K36me3_down$value <- si_H3K36me3_down

### draw a boxplot
devtools::source_url("https://github.com/sa-lee/plyranges/blob/master/R/ranges-bind.R?raw=TRUE")
two <-  bind_ranges(se_K4me1 = se_K4me1_down, si_K4me1 = si_K4me1_down, se_K4me3 = se_K4me3_down, si_K4me3 = si_K4me3_down, 
                    se_K36me1 = se_K36me1_down, si_K36me1 = si_K36me1_down, 
                    se_K36me2 = se_K36me2_down, si_K36me2 = si_K36me2_down, 
                    se_K36me3 = se_K36me3_down, si_K36me3 = si_K36me3_down, .id = "Mutant")
ttl <- "boxplot"

two <- two[mcols(two)$value < 30]
two%>% mcols() %>% as_tibble() %>%
  ggplot(aes(x = factor (Mutant, level=c('se_K4me3', 'si_K4me3', 'se_K36me3', 'si_K36me3', 'se_K4me1', 'si_K4me1', 'se_K36me2', 'si_K36me2', 'se_K36me1', 'si_K36me1')), y = value, fill = Mutant)) +
  geom_boxplot(width = 0.5, outlier.colour = NA, fill = "white") +
  #ylim (0, 7)
  ggtitle(ttl)
ggsave(paste0(ttl, ".pdf"), width = 7, height = 7, units = "in")


###cound the overlaps for UPstream transcripts
pas_up_H3K4me1 <- countOverlaps(pas_up_only,H3K4me1.ab,type="any")
pas_up_H3K4me3 <- countOverlaps(pas_up_only,H3K4me3.epi,type="any")
pas_up_H3K36me1 <- countOverlaps(pas_up_only,H3K36me1.bo,type="any")
pas_up_H3K36me2 <- countOverlaps(pas_up_only,H3K36me2.ab,type="any")
pas_up_H3K36me3 <- countOverlaps(pas_up_only,H3K36me3.ac,type="any")
pas_up_H3K27me3 <- countOverlaps(pas_up_only,H3K27me3.milli,type="any")

upstream_signal_k4me1 <- pas_up_only
upstream_signal_k4me1$signal <- pas_up_H3K4me1

upstream_signal_k4me3 <- pas_up_only
upstream_signal_k4me3$signal <- pas_up_H3K4me3

upstream_signal_k36me1 <- pas_up_only
upstream_signal_k36me1$signal <- pas_up_H3K36me1

upstream_signal_k36me2 <- pas_up_only
upstream_signal_k36me2$signal <- pas_up_H3K36me2

upstream_signal_k36me3 <- pas_up_only
upstream_signal_k36me3$signal <- pas_up_H3K36me3

upstream_signal_k27me3 <- pas_up_only
upstream_signal_k27me3$signal <- pas_up_H3K27me3


###cound the overlaps for medium expressed transcripts_1000 genes
pas_medium_H3K4me1 <- countOverlaps(pas_medium_only,H3K4me1.ab,type="any")
pas_medium_H3K4me3 <- countOverlaps(pas_medium_only,H3K4me3.epi,type="any")
pas_medium_H3K36me1 <- countOverlaps(pas_medium_only,H3K36me1.bo,type="any")
pas_medium_H3K36me2 <- countOverlaps(pas_medium_only,H3K36me2.ab,type="any")
pas_medium_H3K36me3 <- countOverlaps(pas_medium_only,H3K36me3.ac,type="any")
pas_medium_H3K27me3 <- countOverlaps(pas_medium_only,H3K27me3.milli,type="any")

medium_signal_k4me1 <- pas_medium_only
medium_signal_k4me1$signal <- pas_medium_H3K4me1

medium_signal_k4me3 <- pas_medium_only
medium_signal_k4me3$signal <- pas_medium_H3K4me3

medium_signal_k36me1 <- pas_medium_only
medium_signal_k36me1$signal <- pas_medium_H3K36me1

medium_signal_k36me2 <- pas_medium_only
medium_signal_k36me2$signal <- pas_medium_H3K36me2

medium_signal_k36me3 <- pas_medium_only
medium_signal_k36me3$signal <- pas_medium_H3K36me3

medium_signal_k27me3 <- pas_medium_only
medium_signal_k27me3$signal <- pas_medium_H3K27me3



#boxplot
devtools::source_url("https://github.com/sa-lee/plyranges/blob/master/R/ranges-bind.R?raw=TRUE")
two <-  bind_ranges(K4me1_up = upstream_signal_k4me1, K4me1_med = medium_signal_k4me1, K4me3_up = upstream_signal_k4me3, K4me3_med = medium_signal_k4me3, 
                    K36me1_up = upstream_signal_k36me1, K36me1_med = medium_signal_k36me1, 
                    K36me2_up = upstream_signal_k36me2, K36me2_med = medium_signal_k36me2, 
                    K36me3_up = upstream_signal_k36me3, K36me3_med = medium_signal_k36me3,
                    K27me3_up = upstream_signal_k27me3, K27me3_med = medium_signal_k27me3, .id = "Mutant")
ttl <- "boxplot"

two <- two[mcols(two)$signal < 30]
two%>% mcols() %>% as_tibble() %>%
  #two1 %>% mcols() %>% as_tibble() %>%
  ggplot(aes(x = factor (Mutant, level=c('K4me3_up', 'K4me3_med', 'K36me3_up', 'K36me3_med', 'K4me1_up', 'K4me1_med','K36me2_up', 'K36me2_med', 'K36me1_up', 'K36me1_med', 'K27me3_up', 'K27me3_med')), y = signal, fill = Mutant)) +
  #geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.5, outlier.colour = NA, fill = "white") +
  #ylim (0, 7)
  ggtitle(ttl)
ggsave(paste0(ttl, ".pdf"), width = 7, height = 7, units = "in")

