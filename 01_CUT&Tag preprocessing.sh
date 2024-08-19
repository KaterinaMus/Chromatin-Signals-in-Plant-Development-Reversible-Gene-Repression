log=log.txt

# Count raw reads:
echo "++++ Counting raw reads ++++" >> $log
for file in *.fastq; do 
  echo $file $(( $(cat $file | wc -l | 
    awk '{print $1}') / 4 )) >> $log; 
done

# Trim adapters
for file in *.fastq; do 
  echo $file && 
  trim_galore --nextera --length 15 --no_report_file -o ./ $file
done

# Count adapter trimmed reads:
echo "++++ Counting adapter trimmed reads ++++" >> $log
for file in *_trimmed.fq; do 
  echo $file $(( $(cat $file | wc -l | 
    awk '{print $1}') / 4 )) >> $log; 
done

# Alignment with STAR
genome=/path/to/STARgenome/folder

for file in *_trimmed.fq; do 
  echo $file && 
  STAR --genomeDir $genome --readFilesIn $file \
    --runThreadN 20 --outFileNamePrefix ${file/.fq/_} \
    --outSAMmultNmax 1 \
	--alignEndsType Local \
    --clip3pAdapterSeq "AGATCGGAAGAGC" \
    --outSAMtype BAM Unsorted;
done

rm *out *tab; rmdir *STARtmp
for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done

# Count aligned reads:
echo "++++ Counting aligned reads ++++" >> $log
for file in *bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' \
    | awk '{print $1}') >> $log; 
done

# Sort BAM files and remove low MAPQ reads:
for file in *bam; do 
  echo $file && 
  samtools view -hu -q 10 $file | 
  samtools sort - -o ${file/.bam/_sorted_mapq.bam}; 
done

# Count filtered reads:
echo "++++ Counting mapq filtered reads ++++" >> $log
for file in *mapq.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}') >> $log; 
done

# Deduplicate in PE mode:
for file in *_mapq.bam; do echo $file && samtools markdup -r $file ${file/.bam/_dedup.bam}; done

# Generate indexes:
for file in *_mapq_dedup.bam; do samtools index $file; done

# Count deduplicated reads:
echo "++++ Counting deduplicated reads ++++" >> $log
for file in *sorted_mapq_dedup.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}') >> $log; 
done

#Merge replicates
samtools merge Sample1_Rep11.bam Sample1_Rep2.bam Sample1_Rep.bam Sample1_Merged.bam

# Count reads in merged files to check:
for file in *dedup.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Coverage and Normalization
for file in *sorted_names_fixmate_sorted_dedup.bam; do echo $file && bamCoverage -b $file --binSize 10 --normalizeUsing BPM --effectiveGenomeSize 135000000 --extendReads 76 --outFileFormat bedgraph -o ${file/_dedup.bam/_cov.bg}; done
 
gzip *cov.bg

log=MACS2_log.txt

for file in *_trimmed_sorted_mapq_dedup.bam; do 
	echo "Processing file $file"
	echo "############# $file #############" >> $log
	# Narrow peak calling: -g is genome size; -m is the lower and upper limits of fold change to include; --scale-to large means scale the smaller file to larger file to normalise; -n output file name; -B output bg; -p pValue; -q qValue cutoff
	macs2 callpeak -t $file -f BAM -g 1.35e+08 -p 1e-5 --keep-dup all --scale-to large --SPMR -n ${file/.bam/_macs2} -B --outdir MACS2_Individual 2>> $log
done

cd MACS2_Individual

for file in *.narrowPeak; do
	echo "Processing file $file"
	awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,$7}' $file >> ${file/.narrowPeak/_narrowPeak.bdg}
done
