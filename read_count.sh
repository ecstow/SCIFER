cd 
#Combine read counts from all barcode-specific coverage.txt files for each cluster.
cut -f 10 *plus.txt >> all_read_counts_plus.txt
cut -f 10 *minus.txt >> all_read_counts_minus.txt

#
awk '{a[NR%305] = a[NR%305] (NR<=305 ? "" : ",") $0}
      END{for (i = 1; i <= 305; i++) print a[i%305]}' < all_read_counts_plus.txt > all_plus.csv
awk '{a[NR%305] = a[NR%305] (NR<=305 ? "" : ",") $0}
       END{for (i = 1; i <= 305; i++) print a[i%305]}' < all_read_counts_minus.txt > all_minus.csv

