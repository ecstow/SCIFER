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
 
 #Create a list of unique barcode UMIs and duplicated barcode UMIs. Extract reads from each list and only keep one of each duplicated barcode-UMI.
 grep -w "1" ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51_count.txt > ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51_count_unique.txt
 grep -v -w "1" ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51_count.txt > ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51_count_nonunique.txt 
 cut -b 9-66 ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51_count_unique.txt > ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51_count_unique_cut.txt
 cut -b 9-66 ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51_count_nonunique.txt > ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51_count_nonunique_cut.txt 
 mkdir sam_files
 while read line in file; do grep $line ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s.sam > ${line}_unique.sam; done < ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51_count_unique_cut.txt
 while read line in file; do grep $line ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s.sam > ${line}_unique.sam; done < ${prefix}_bowtie_hg38_sorted_rmdup_top_bottom_expressedL1s_51_count_nonunique_cut.txt 
 mv *_unique.sam sam_files/
 
 #Add header to all output .sam files, convert to .bam, and merge all files.
 for file in sam_files/*.sam; do (head -30 <top30_sam_hg38_MCF7.txt> ; cat $file) > ${file}_header.sam; done
 for file in sam_files/*_header.sam; do samtools view -bS $file> ${file}_merge.bam; done
 samtools merge ${prefix}_L1exp_reads_unique.bam sam_files/*_merge.bam
 samtools view ${prefix}_L1exp_reads_unique.bam > ${prefix}_L1exp_reads_unique.sam
 
 #Extract cell barcodes into individual .sam files, convert to bam, strand separate reads and count reads that overlap with repeat element annotations.
 mkdir sam_files_2
 while read line in file; do grep $line ${prefix}_L1exp_reads_unique.sam
 | head -1  > sam_files_2/${line}_barcode.sam; done < <list_of_cell_barcodes.txt>
 for file in sam_files_2/*.sam; do (head -30 <top30_sam_hg38_MCF7.txt> ; cat $file) > ${file}_header.sam; done
 for file in sam_files_2/*_header.sam; do samtools view -bS $file> sam_files_2/${file}.bam; done
 
 #Strand separate reads and count reads that overlap with repeat element annotations.
 samtools view -h sam_files_2/${file}.bam | awk 'substr($0,1,1) == "@" || $2 == 0 {print}' | samtools view -bS - > sam_files_2/${file}_top.bam
 samtools view -h sam_files_2/${file}.bam | awk 'substr($0,1,1) == "@" || $2 == 16 {print}' | samtools view -bS - > sam_files_2/${file}_bottom.bam
 mkdir read_counts
 bedtools coverage -abam <expressed_L1_bedfile_plus.bed> -b sam_files_2/${file}_top.bam > read_counts/${file}_plus.txt
 bedtools coverage -abam <expressed_L1_bedfile_minus.bed> -b sam_files_2/${file}_bottom.bam > read_counts/${file}_minus.txt
 

