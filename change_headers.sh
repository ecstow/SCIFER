prefix=${1%_[1-2].fastq}  

#Reassign barcode-UMIs in R1 fastq file to header in R2 fastq file.
#Print every other line to obtain only barcode-UMIs/indeces and scores. Delete even lines to obtain only barcode-UMIs.

awk '!(NR%2)' ${prefix}"_1.fastq" > ${prefix}_barcode_score.fastq
sed 'n; d' ${prefix}_barcode_score.fastq > ${prefix}_barcode.fastq

#Print first 26 letters of reads in R1.fastq to obtain barcodes and UMIs. Then print first 16 letters of reads in R1.fastq to obtain barcodes.

cut -c-28 ${prefix}_barcode.fastq > ${prefix}_barcode_28char.fastq
cut -c-16 ${prefix}_barcode.fastq > ${prefix}_barcode_16char.fastq

#Add delimiter to end of barcodes and combine with list of barcode-UMIs.

sed 's/$/\./' ${prefix}_barcode_16char.fastq > ${prefix}_barcode_16char_period.fastq
paste -d'\0' ${prefix}_barcode_16char_period.fastq ${prefix}_barcode_28char.fastq > ${prefix}_barcode_28char_16char_period.fastq

#Print headers from read2 file. Change as needed for different read names.

grep -w “@SRR6129051” ${prefix}"_2.fastq" > ${prefix}_header.fastq
sed 's/233:.*//g' ${prefix}_header.fastq > ${prefix}_header1.fastq

#Paste together file containing the beginning half of header and the file containing barcode-UMIs together.

paste -d'\0' ${prefix}_header1.fastq ${prefix}_barcode_28char_16char_period.fastq > ${prefix}_header1_barcode.fastq
sed '/@SRR6129051/d' ${prefix}"_2.fastq" > ${prefix}_no_header.fastq 
sed 'G;G;G' ${prefix}_barcode_28char_16char_period.fastq > ${prefix}_barcode_28char_16char_period_spaced.fastq
sed '0~3 a\\' ${prefix}_no_header.fastq > ${prefix}_no_header_spaced.fastq 
sed '1i\\' ${prefix}_no_header_spaced.fastq > ${prefix}_no_header_spaced_topline.fastq 
paste -d'\0' ${prefix}_barcode_28char_16char_period_spaced.fastq ${prefix}_no_header_spaced_topline.fastq > ${prefix}_reads_barcodes.fastq
sed 's/ //g' ${prefix}_reads_barcodes.fastq > ${prefix}_reads_barcodes_nospace.fastq
