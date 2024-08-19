#This script was used to analyse raw TSS-seq fastq files obtained by SE 75bp HighOutput 550 NextSeq sequencing

# Trim Illumina adapter sequences from the reads:
for file in *fastq.gz; do echo $file && trim_galore --illumina $file; done

# Process UMIs (first 8 nt of each read), append them to Fastq, step required for UMI-Tools dedup:
for file in *trimmed.fq.gz; do echo $file && umi_tools extract --stdin=${file} --bc-pattern=NNNNNNNN --stdout=${file/.fq.gz/_UMI.fq.gz}; done

# Align to TAIR10 using STAR:
for file in *trimmed_UMI.fq.gz; do echo $file &&  STAR --runThreadN 16 --runMode alignReads --genomeDir ~/index --readFilesIn $file  --alignEndsType EndToEnd --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 3078189601 --outFileNamePrefix ${file/.fq.gz/_}; done

# Filter out multimapper reads:
for file in *sortedByCoord.out.bam; do echo $file && samtools view -h -q 10 $file -o ${file/.bam/_mapq10.bam}; done

#index uniq bam files
for file in *_mapq10.bam; do echo $file && samtools index $file; done

#remove duplicated reads
for file in *mapq10.bam; do echo $file && ~/umi_tools dedup -I $file -S ${file/.bam/_dedup.bam}; done

#create bedgraph files
for file in *mapq10_dedup.bam; do bedtools genomecov -ibam $file -bg -10 > ${file/.bam/.bedgraph};done

# Make stranded Bedgraph files (only the first base of each read is considered):
for str in "+" "-"; do
  echo $str
  [ "$str" = "+" ] && n="fw" || n="rev"
  for file in *mapq10_dedup.bam; do
    echo $file && bedtools genomecov -ibam $file -bg -5 -strand $str | sort -k 1,1 -k 2,2n > ${file/.bam/}_${n}.bg
  done
done

# Compress the bedgraph files
for file in *.bg; do gzip $file; done
