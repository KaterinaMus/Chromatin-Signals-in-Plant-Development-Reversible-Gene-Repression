library(rtracklayer)
library(tidyverse)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
library(seqinr)
library(GenomicRanges)

win_genes <- genes(txdb, columns = c("gene_id", "tx_type"))
mcols(win_genes)$tx_type <- mcols(win_genes)$tx_type %>% unlist(use.names = FALSE)

scripts <- c("batch_read_track_data", "merge_and_normalize_GRanges", "metagene_matrix", "draw_metagene_plot", "get_overlapping_scores")
for (script in scripts) {
  paste0("https://github.com/Maxim-Ivanov/Utility_functions/blob/main/", script, ".R?raw=TRUE") %>% devtools::source_url()
}


url <- "https://github.com/Maxim-Ivanov/Kindgren_et_al_2019/blob/master"
scripts <- c("batchReadTrackData", "getOverlappingScores", "metageneMatrix", "drawMetagenePlot", "normalizeGR")
for (script in scripts) {
  paste0(url, "/", script, ".R?raw=TRUE") %>% devtools::source_url()
}


path <- "." # change to the directory containing merged ChIP-seq or CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
#CUT&Tag
H3K4me1 <- import('.',format="bedgraph") # change to the directory containing merged ChIP-seq (PMID 28100676) or CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K4me3 <- import('.',format="bedgraph") # change to the directory containing merged ChIP-seq (PMID 28100676) or CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K27me3 <- import('.',format="bedgraph") # change to the directory containing merged ChIP-seq (PMID 22962860) or CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me2 <- import('.',format="bedgraph") # change to the directory containing merged ChIP-seq (PMID 22962860) or CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh
H3K36me3 <- import('.',format="bedgraph") # change to the directory containing merged ChIP-seq (PMID 22962860) or CUT&Tag Bedgraph files obtained from 01-CUT&Tag processing.sh


#process them
H3K4me1 <- dropSeqlevels(H3K4me1,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K4me3 <- dropSeqlevels(H3K4me3,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K27me3 <- dropSeqlevels(H3K27me3,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me2 <- dropSeqlevels(H3K36me2,c("ChrM","ChrC"),pruning.mode ="coarse")
H3K36me3 <- dropSeqlevels(H3K36me3,c("ChrM","ChrC"),pruning.mode ="coarse")

#seqnames(PTM)
H3K4me1.a1 <- resize(H3K4me1, width = 1, fix = "center")
H3K4me3.a3 <- resize(H3K4me3, width = 1, fix = "center")
H3K27me3.1 <- resize(H3K27me3, width = 1, fix = "center")
H3K36me2.a2 <- resize(H3K36me2, width = 1, fix = "center")
H3K36me3.a3 <- resize(H3K36me3, width = 1, fix = "center")


#down select genes up to 1000bp
coding_genes <- genes(txdb, columns = c("gene_id", "tx_type"))
coding_genes_n0big <- coding_genes[width(coding_genes) >= 800, ]
mcols(coding_genes_n0big)$tx_type <- mcols(coding_genes_n0big)$tx_type %>% unlist(use.names = FALSE)

#extend function
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


#shrinks genes to gene bodies
shrinked_genes <- extend(coding_genes_n0big,upstream=-150,downstream=-150)
shrinked_genes <- dropSeqlevels(shrinked_genes,c("C","M","Pt","Mt"),pruning.mode ="coarse")

#define TSSs for the calculation
TSS.2 <- resize(coding_genes_n0big, width = 1 , fix = "start")
extended_TSS.2 <- extend(TSS.2,upstream=150,downstream=150)
extended_TSS.2 <- dropSeqlevels(extended_TSS.2,c("Mt","Pt"),pruning.mode ="coarse")

#create list of histone marks for calculation
data <- list("h3k4me1" = H3K4me1.a1,"h3k4me3" = H3K4me3.a3,"h3k27me3" = H3K27me3.1,"h3k36me2" = H3K36me2.a2,"h3k36me3" = H3K36me3.a3)

# get overlapping scores
getOverlappingScores_v2 <- function(intervals, signal_grl, trim.names=c(0, 0)) {
  require(GenomicRanges)
  before <- trim.names[[1]]; after <- trim.names[[2]]
  names(signal_grl) <- substr(names(signal_grl), before+1, nchar(names(signal_grl))-after)
  results <- vector("list", length(signal_grl))
  names(intervals) <- 1:length(intervals)
  int_plus <- intervals[strand(intervals)=="+"]
  int_minus <- intervals[strand(intervals)=="-"]
  int_star <- intervals[strand(intervals)=="*"]
  intervals_split <- list(int_plus, int_minus, int_star)
  strands <- list("+", "-", c("+", "-"))
  for (i in seq_along(signal_grl)) {
    signal <- signal_grl[[i]]
    message(names(signal_grl)[[i]]); flush.console()
    res <- vector("list", 3)
    for (j in 1:3) {
      curr_strands <- c(strands[[j]], "*")
      curr_signal <- signal[strand(signal) %in% curr_strands]
      curr_intervals <- intervals_split[[j]]
      curr_cov <- coverage(curr_signal)
      cov_by_intervals <- curr_cov[curr_intervals]
      numlist <- as(cov_by_intervals, "NumericList")
      integrated <- lapply(numlist, FUN=sum)
      names(integrated) <- names(curr_intervals)
      res[[j]] <- integrated
    }
    res <- c(res[[1]], res[[2]], res[[3]])
    res <- res[order(as.numeric(names(res)))]
    res <- as.numeric(res)
    results[[i]] <- res
  }
  results <- as.data.frame(t(do.call(rbind, results)))
  names(results) <- names(signal_grl)
  mcols(intervals) <- cbind(mcols(intervals), results)
  return(intervals)
}

#to remove Chr in seqnames

data$h3k4me1 <- renameSeqlevels(data$h3k4me1, gsub("Chr","", seqlevels(data$h3k4me1)))
data$h3k4me3 <- renameSeqlevels(data$h3k4me3, gsub("Chr", "", seqlevels(data$h3k4me3)))
data$h3k27me3 <- renameSeqlevels(data$h3k27me3, gsub("Chr", "", seqlevels(data$h3k27me3)))
data$h3k36me2 <- renameSeqlevels(data$h3k36me2, gsub("Chr", "", seqlevels(data$h3k36me2)))
data$h3k36me3 <- renameSeqlevels(data$h3k36me3, gsub("Chr", "", seqlevels(data$h3k36me3)))

seqinfo(data$h3k4me1) <- seqinfo(txdb)
seqinfo(data$h3k4me3) <- seqinfo(txdb)
seqinfo(data$h3k27me3) <- seqinfo(txdb)
seqinfo(data$h3k36me2) <- seqinfo(txdb)
seqinfo(data$h3k36me3) <- seqinfo(txdb)

trim(data$h3k4me1)
trim(data$h3k4me3)
trim(data$h3k27me3)
trim(data$h3k36me2)
trim(data$h3k36me3)

data$h3k4me1 <- trim(data$h3k4me1)
data$h3k4me3 <- trim(data$h3k4me3)
data$h3k27me3 <- trim(data$h3k27me3)
data$h3k36me2 <- trim(data$h3k36me2)
data$h3k36me3 <- trim(data$h3k36me3)


#calculate the TPMs for each regions of genes
K4ME.tss <- getOverlappingScores_v2(extended_TSS.2,data)
K4ME.genebody <- getOverlappingScores_v2(shrinked_genes,data)

#combine the genebody values to the TSS genomic range

K4ME.tss$h3k4me1.gb <- K4ME.genebody$h3k4me1
K4ME.tss$h3k4me3.gb <- K4ME.genebody$h3k4me3
K4ME.tss$h3k27me3.gb <- K4ME.genebody$h3k27me3
K4ME.tss$h3k36me2.gb <- K4ME.genebody$h3k36me2
K4ME.tss$h3k36me3.gb <- K4ME.genebody$h3k36me3

#calculates FPKM from TPM

K4ME.tss$h3k4me1.gb.norm <- K4ME.tss$h3k4me1.gb/width(K4ME.tss)
K4ME.tss$h3k4me3.gb.norm <- K4ME.tss$h3k4me3.gb/width(K4ME.tss)
K4ME.tss$h3k27me3.gb.norm <- K4ME.tss$h3k27me3.gb/width(K4ME.tss)
K4ME.tss$h3k36me2.gb.norm <- K4ME.tss$h3k36me2.gb/width(K4ME.tss)
K4ME.tss$h3k36me3.gb.norm <- K4ME.tss$h3k36me3.gb/width(K4ME.tss)


### Load TIF-seq data
gr_wt <- import (".") # change to the directory containing TIF-seq Bedgraph file (GSE129523)
gr <- gr_wt


#Measure the diversity of clusters that overlap or not known genes
genes_nO <- coding_genes_n0big

gr$score <- NULL
gr$name  <- as.numeric(gr$name)
gr$cluster_number <- as.character(seq_len(nrow(mcols(gr))))

#finds clusters overlapping coding genes all genes and ncRNA genes

gr$overlap_nO <- countOverlaps(gr,genes_nO,type=c("any"),ignore.strand=FALSE)

intragenic_nO <- gr[gr$overlap_nO != 0]

#identify where overlapping clusters of one orf locate
intragenic_nO$overlap_5 <- countOverlaps(intragenic_nO,extended_TSS.2,type=c("any"),ignore.strand=FALSE)
candidates <- intragenic_nO[intragenic_nO$overlap_5 >= 1 & intragenic_nO$overlap_nO >= 2]

#get overlapping candidates
K4ME.tss$overlap_candidatesTIF <- countOverlaps(K4ME.tss,candidates,type=c("any"),ignore.strand=FALSE)

dat <- as.data.frame(K4ME.tss)

#rank the genes by H3k4me marks ups and down (high H3K4me1 and low H3K4me3 and other repressive marks)
dat$r_h3k4me1 <- NA
dat$r_h3k4me1[order(dat$h3k4me1)] <- 1:nrow(dat)/nrow(dat)
dat$r_h3k4me3 <- NA
dat$r_h3k4me3[order(-dat$h3k4me3)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k27me3 <- NA
dat$r_h3k27me3[order(-dat$h3k27me3)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k36me2 <- NA
dat$r_h3k36me2[order(dat$h3k36me2)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k36me3 <- NA
dat$r_h3k36me3[order(dat$h3k36me3)] <- (1:nrow(dat))/nrow(dat)

dat$r_h3k4me1.gb.norm <- NA
dat$r_h3k4me1.gb.norm[order(dat$h3k4me1.gb.norm)] <- 1:nrow(dat)/nrow(dat)
dat$r_h3k4me3.gb.norm <- NA
dat$r_h3k4me3.gb.norm[order(-dat$h3k4me3.gb.norm)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k27me3.gb.norm <- NA
dat$r_h3k27me3.gb.norm[order(-dat$h3k27me3.gb.norm)] <- (1:nrow(dat))/nrow(dat)
dat$r_h3k36me2.gb.norm <- NA
dat$r_h3k36me2.gb.norm[order(dat$h3k36me2.gb.norm)] <- 1:nrow(dat)/nrow(dat)
dat$r_h3k36me3.gb.norm <- NA
dat$r_h3k36me3.gb.norm[order(-dat$h3k36me3.gb.norm)] <- (1:nrow(dat))/nrow(dat)

#calculate h3k4me1/h3k4me3 ratio and rank it

dat$ra13 <- dat$h3k4me1/dat$h3k4me3
dat$ra36_23 <- dat$h3k36me2/dat$h3k36me3

dat$r_ra13 <- NA
dat$r_ra13[order(dat$ra13)] <- (1:nrow(dat))/nrow(dat)
dat$r_ra36_23 <- NA
dat$r_ra36_23[order(dat$ra36_23)] <- (1:nrow(dat))/nrow(dat)

#Import plaNET-Seq data for gene expression ranking
pNET <- import('.', format ='bedgraph')  # change to the directory containing plaNET-seq Bedgraph file (GSE131733)
pnet <- dropSeqlevels(pNET,c("Pt","Mt"),pruning.mode ="coarse")

seqlevels(pnet) <- c("1", "2", "3", "4", "5")

data <- list("pNET" = pnet)
Norm <- getOverlappingScores_v2(shrinked_genes,data)
txdb2 <- dropSeqlevels(txdb,c("6","7"),pruning.mode ="coarse")
Norm$normpNET <- Norm$pNET/width(Norm)

dat$pNET <- Norm$normpNET
dat$r_pNET <- NA
dat$r_pNET[order(dat$pNET)] <- (1:nrow(dat))/nrow(dat)

#subset genes with high ration and high h3k4me1

dat.2 <- subset(dat,dat$ra13>1)
dat.2 <- subset(dat.2,dat.2$h3k4me1>24) # 24 for CUT&Tag data, 214 for ChIP-seq data

#calculate a score to find best candidates

s1 = 15 #H3K4me1 (X1)
s2 = 10 #h3k4me3 (X2)
s3 = 10 #h3k36me2 (X3)
s4 = 5 #h3k36me3 (X4)
s5 = 20 #ratio H3K4me1/3 (Y1)
s6 = 10 #ratio H3K36me2/3 (Y2)
s7 = 5 #h3k27me3 (Z1)

dat.2$score2 <- (s1*(dat.2$r_h3k4me1) + s2*(dat.2$r_h3k4me3) + s3*(dat.2$r_h3k36me2) + s4*(dat.2$r_h3k36me3) + s5*(dat.2$r_ra13) + s6*(dat.2$r_ra36_23) + s7*(dat.2$r_h3k27me3))
dat.2$tx_type <- "protein_coding"

write.csv(dat.2,".csv")



