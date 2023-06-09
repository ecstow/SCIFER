# SCIFER
Single Cell Implementation to Find Expressed Retrotransposons (SCIFER) aligns scRNA-Seq reads uniquely to the genome and extracts alignments from single cells by cell-specific barcodes. Retrotransposon expression should be confirmed with a matched bulk RNA-Seq sample.

![SCIFER_workflow](https://user-images.githubusercontent.com/108097317/232626714-7f667ba6-46e0-426b-94d9-78ffb633c2eb.png#gh-light-mode-only)
![SCIFER_workflow_dark](https://user-images.githubusercontent.com/108097317/232627059-21786442-3ef5-41e3-b86e-4b5de7c9821d.png#gh-dark-mode-only)

Download from GitHub:
```
git clone https://github.com/ecstow/SCIFER
```

Before running SCIFER, you will need the hg38 reference genome in fasta format, with Bowtie1 index. Add path to the index file in Step2:
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
zcat hg38.fa.gz > hg38.fa
bowtie-build hg38.fa hg38
```

# Additional Details
Our Mobile DNA paper introducing SCIFER: https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-022-00276-0
