# Make sure you're outside of raw_gzipped_fastqs_symlinks

for file in raw_gzipped_fastqs_symlinks/*R1.fastq.gz; 
do 
READS=$(expr $(zcat $file | wc -l) / 4)
echo "$file,$READS"
done > read_amounts.csv

