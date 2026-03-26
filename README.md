# barcodes-id
taxonomic identification of DNA barcodes

```sh
# download bold distilled database for blast and sintax and upzip
wget https://us-sea-1.linodeobjects.com/boldistilled/blast.zip -O temp/bold-distilled-blast.zip
wget https://us-sea-1.linodeobjects.com/boldistilled/sintax.zip -O temp/bold-distilled-sintax.zip

wget https://us-sea-1.linodeobjects.com/boldistilled/source.zip -O temp/bold-distilled-source.zip

unzip temp/bold-distilled-blast.zip -d temp
unzip temp/bold-distilled-sintax.zip -d temp
unzip temp/bold-distilled-source.zip -d temp

# run blastn v2.14 local on bold
blastn -task blastn -num_threads 8 -evalue 0.1 -word_size 11 -max_target_seqs 50 -db temp/blast/BOLDistilled_COI_Jan2026_SEQUENCES.fasta -outfmt "6 qseqid sseqid evalue length pident nident bitscore" -out temp/cox1-blast-local.out -query temp/cox1.fasta

# run blast remote
blastn -task blastn -remote -db core_nt -evalue 0.1 -word_size 11 -max_target_seqs 50 -outfmt "6 qseqid sgi evalue length pident nident bitscore" -out temp/cox1-blast-remote.out -query temp/cox1.fasta


# run sintax
vsearch --threads 8 --sintax temp/cox1.fasta --db temp/sintax/BOLDistilled_COI_Jan2026_SEQUENCES_sintax.fasta --sintax_cutoff 0.7 --tabbedout temp/cox1-sintax-out.tsv


# process results
scripts/process-results.R



```