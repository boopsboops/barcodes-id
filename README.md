# barcodes-id
taxonomic identification of DNA barcodes

### Identification of cox1 sequences

```bash
# download bold distilled database for blast and sintax and upzip
wget https://us-sea-1.linodeobjects.com/boldistilled/blast.zip -O temp/bold-distilled-blast.zip
wget https://us-sea-1.linodeobjects.com/boldistilled/sintax.zip -O temp/bold-distilled-sintax.zip
wget https://us-sea-1.linodeobjects.com/boldistilled/source.zip -O temp/bold-distilled-source.zip
unzip temp/bold-distilled-blast.zip -d temp
unzip temp/bold-distilled-sintax.zip -d temp
unzip temp/bold-distilled-source.zip -d temp

# run blastn v2.14 local on bold
blastn -task blastn -num_threads 8 -evalue 0.1 -word_size 11 -max_target_seqs 50 -db temp/blast/BOLDistilled_COI_Jan2026_SEQUENCES.fasta -outfmt "6 qseqid sseqid evalue length pident nident bitscore" -out temp/cox1-blast-local.out -query temp/cox1.fasta

# run sintax
vsearch --threads 8 --sintax temp/cox1.fasta --db temp/sintax/BOLDistilled_COI_Jan2026_SEQUENCES_sintax.fasta --sintax_cutoff 0.7 --tabbedout temp/cox1-sintax-out.tsv

# process results
scripts/cox1-process-results.R
```


### Identification of 16S sequences

```bash
wget https://www.reference-midori.info/download/Databases/GenBank269_2025-12-09/SINTAX/uniq/MIDORI2_UNIQ_NUC_GB269_lrRNA_SINTAX.fasta.gz -O temp/MIDORI2_UNIQ_NUC_GB269_lrRNA_SINTAX.fasta.gz
wget https://www.reference-midori.info/download/Databases/GenBank269_2025-12-09/BLAST/uniq/MIDORI2_UNIQ_NUC_GB269_lrRNA_BLAST.zip -O temp/MIDORI2_UNIQ_NUC_GB269_lrRNA_BLAST.zip

unzip temp/MIDORI2_UNIQ_NUC_GB269_lrRNA_BLAST.zip -d temp
gzip -d -k temp/MIDORI2_UNIQ_NUC_GB269_lrRNA_SINTAX.fasta.gz 


# run blastn v2.14 local on midori
blastn -task blastn -num_threads 4 -evalue 0.1 -word_size 11 -max_target_seqs 50 -db temp/MIDORI2_UNIQ_NUC_GB269_lrRNA_BLAST -outfmt "6 qseqid sseqid evalue length pident nident bitscore" -out temp/16s-blast-local.out -query temp/16s.fasta

# run sintax
vsearch --threads 4 --sintax temp/16s.fasta --db temp/MIDORI2_UNIQ_NUC_GB269_lrRNA_SINTAX.fasta --sintax_cutoff 0.7 --tabbedout temp/16s-sintax-out.tsv

# process results
scripts/16s-process-results.R


```