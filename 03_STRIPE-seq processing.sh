#This script was used to analyse raw STRIPE-seq fastq files obtained by SE 75bp HighOutput 550 NextSeq sequencing

# R1 reads are expected to begin with 8nt UMI followed by TATAGGG 

# Count raw reads:
for file in *fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | 
    awk '{print $1}') / 4 )); 
done

# Run FastQC:
for file in *fastq.gz; do fastqc $file; done

# Process UMI:
for file in *fastq.gz; do 
  echo $file &&
  umi_tools extract --stdin=$file --bc-pattern=NNNNNNNN --stdout=../${file/fastq.gz/UMI.fq.gz};
done

# Trim TATAGGG from 5' end and Illumina adapters from 3' end:
basedir=~/
for file in *UMI.fq.gz; do 
  echo $file &&
  ${basedir}/trim_galore --illumina --clip_R1 7 --length 15 --gzip --no_report_file $file; 
done

# Count trimmed reads:
for file in *trimmed.fq.gz; do 
  echo $file $(( $(zcat $file | wc -l | 
    awk '{print $1}') / 4 )); 
done

# Align to genome using STAR:
for file in *UMI_trimmed.fq.gz; do echo $file &&  STAR --runThreadN 16 --outSAMmultNmax 1 --genomeDir ~/index --readFilesIn $file  --alignEndsType Extend5pOfRead1 --readFilesCommand zcat --outSAMtype BAM Unsorted --sjdbGTFtagExonParentTranscript Parent --outSAMstrandField intronMotif --outFileNamePrefix ${file/trimmed.fq.gz/_}; done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARgenome *STARtmp *out *tab

# Count aligned reads:
for file in *bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' \
    | awk '{print $1}'); 
done

# Sort BAM files and remove low MAPQ reads:with low mapping quality
for file in *_.bam; do 
  echo $file && 
  samtools view -hu -q 10 $file | 
  samtools sort - -o ${file/.bam/mapq10.bam}; 
done

# Count filtered reads:
for file in *mapq10.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Generate indexes:
for file in *mapq10.bam; do samtools index $file; done

# Deduplicate on UMI:
for file in *mapq10.bam; do
  echo $file &&
  umi_tools dedup --stdin=${file} --stdout=${file/.bam/_dedup.bam};
done

# Count deduplicated reads:
for file in *dedup.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Make stranded Bedgraph files (only the first base of each read is considered):
for str in "+" "-"; do
  echo $str
  [ "$str" = "+" ] && n="fw" || n="rev"
  for file in *_dedup.bam; do
    echo $file && bedtools genomecov -ibam $file -bg -5 -strand $str | sort -k 1,1 -k 2,2n > ${file/.bam/}_${n}.bg
  done
done

# Compress the bedgraph files
for file in *.bg; do gzip $file; done
