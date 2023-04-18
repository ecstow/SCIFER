#!/bin/bash

#Usage: bash Step3_readcount.sh

#Combine read counts from all barcode-specific coverage.txt files for each cluster.
cut -f 10 *plus.txt >> all_read_counts_plus.txt
cut -f 10 *minus.txt >> all_read_counts_minus.txt

#Sort read counts
awk '{a[NR%2483] = a[NR%2483] (NR<=2483 ? "" : ",") $0}
      END{for (i = 1; i <= 2483; i++) print a[i%2483]}' < all_read_counts_plus.txt > all_plus.csv
awk '{a[NR%2480] = a[NR%2480] (NR<=2480 ? "" : ",") $0}
       END{for (i = 1; i <= 2480; i++) print a[i%2480]}' < all_read_counts_minus.txt > all_minus.csv
