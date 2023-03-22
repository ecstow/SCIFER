module load samtools/1.5
module load bedtools/2.27.1 

prefix=${1%.fastq}  

#Perform unique alignment using bowtie-0.12.8
bowtie -p 10 -m 1 -S -y -v 3 --chunkmbs 8184 <hg38 index> ${prefix}.fastq | samtools view -hbuS - | samtools sort -o ${prefix}_bowtie_hg38_sorted.bam
samtools rmdup -s ${prefix}_bowtie_hg38_sorted.bam ${prefix}_bowtie_hg38_sorted_rmdup.bam

#Extract aligned reads and extract cell-specific alignments based on barcode. Add a colon before each barcode in list_of_barcodes.txt file to ensure the extraction is based on barcode matches alone (instead of reads that potentially contain barcode in sequence).
samtools view -h ${prefix}_bowtie_hg38_sorted_rmdup.bam | awk 'substr($0,1,1) == "@" || $2 == 0 {print}' | samtools view -bS - > ${prefix}_bowtie_hg38_sorted_rmdup_top.bam
samtools view -h ${prefix}_bowtie_hg38_sorted_rmdup.bam | awk 'substr($0,1,1) == "@" || $2 == 16 {print}' | samtools view -bS - > ${prefix}_bowtie_hg38_sorted_rmdup_bottom.bam
samtools merge ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom.bam ${prefix}_bowtie_hg38_sorted_rmdup_top.bam ${prefix}_bowtie_hg38_sorted_rmdup_bottom.bam 
samtools view -h ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom.bam > ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom.sam
while read line in file; do grep $line ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom.sam > ${line}_barcode.sam; done < <list_of_barcodes>

#These files will be used to determine number of mapped reads using
for i in ls *.bam; do echo $(cat ${i} | samtools view -F 0x904 -c)_$i; done

#Extract reads that align to expressed L1 loci and count duplicate barcode-UMIs.
samtools view -b -h -L <expressed_L1_bedfile_plus.bed> ${prefix}_bowtie_hg38_sorted_rmdup_top.bam > ${prefix}_bowtie_hg38_sorted_rmdup_top_expressedL1s.bam
samtools view -b -h -L <expressed_L1_bedfile_minus.bed> ${prefix}_bowtie_hg38_sorted_rmdup_bottom.bam > ${prefix}_bowtie_hg38_sorted_rmdup_bottom_expressedL1s.bam
samtools merge ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s.bam ${prefix}_bowtie_hg38_sorted_rmdup_top_expressedL1s.bam ${prefix}_bowtie_hg38_sorted_rmdup_bottom_expressedL1s.bam
samtools view ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s.bam > ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s.sam

#Cut beginning of header that contains barcode and barcode-UMI. Change as needed for sequence header.
cut -c-51 ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s.sam > ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51.txt

#Count duplicated barcode-UMIs.
sort ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51.txt |uniq -c| sort -n -r > ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51_count.txt

$ awk '{
for ( i=1; i<=NF; i++ )
dict[$i]++;                                                                              
}
END{
for (key in dict)
 if(dict[key]>1)
print dict[key] " : " key
 }' ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51_count.txt

