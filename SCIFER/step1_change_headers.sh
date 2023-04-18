#!/bin/bash

#Usage: for i in *1.fastq; do bash Step1_SCIFER_change_headers.sh $i; done

#This script reassigns barcode-UMIs in R1 fastq file to header in R2 fastq file so that this information is not lost during bowtie alignment of R2.

prefix=${1%_[1-2].fastq}  

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
grep '@NB551227' ${prefix}"_2.fastq" > ${prefix}_header.fastq
sed 's/144:.*//g' ${prefix}_header.fastq > ${prefix}_header1.fastq

#Paste together file containing the beginning half of header and the file containing barcode-UMIs together.
paste -d'\0' ${prefix}_header1.fastq ${prefix}_barcode_28char_16char_period.fastq > ${prefix}_header1_barcode.fastq

#Delete original headers from fastq 2 file
sed '/@NB551227/d' ${prefix}"_2.fastq" > ${prefix}_no_header.fastq 

#Add spaces to fastq file with new headers
sed 'G;G;G' ${prefix}_barcode_28char_16char_period.fastq > ${prefix}_barcode_28char_16char_period_spaced.fastq

#Add spaces to fastq file without header
sed '0~3 a\\' ${prefix}_no_header.fastq > ${prefix}_no_header_spaced.fastq 

#Add top line space
sed '1i\\' ${prefix}_no_header_spaced.fastq > ${prefix}_no_header_spaced_topline.fastq 

#Combine fastq 2 read file with file containing new headers
paste -d'\0' ${prefix}_barcode_28char_16char_period_spaced.fastq ${prefix}_no_header_spaced_topline.fastq > ${prefix}_reads_barcodes.fastq

#Delete spaces in the header to retain barcode information
sed 's/ //g' ${prefix}_reads_barcodes.fastq > ${prefix}_reads_barcodes_nospace.fastq
