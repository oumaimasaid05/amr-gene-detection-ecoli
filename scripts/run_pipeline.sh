#!/usr/bin/env bash
# Full pipeline, in order. Run from the project root inside the "amr" conda environment.
set -e
mkdir -p data/genomes results/amrfinder results/blast results/figures

# 1. Download genomes (accessions in data/accessions.txt)
datasets download genome accession --inputfile data/accessions.txt --include genome --filename data/genomes.zip
(cd data && unzip -o genomes.zip)
for f in data/ncbi_dataset/data/G*_*/*.fna; do
  cp "$f" data/genomes/$(basename $(dirname "$f")).fna
done

# 2. Genome QC and contig table
seqkit stats -a data/genomes/*.fna | tee results/genome_qc.txt
for f in data/genomes/*.fna; do
  seqkit fx2tab -n -l "$f" | awk -v a="$(basename $f .fna)" '{print a"\t"$0}'
done > results/contigs.tsv

# 3. AMRFinderPlus (primary)
amrfinder -u
for f in data/genomes/*.fna; do
  id=$(basename "$f" .fna)
  amrfinder -n "$f" -O Escherichia --plus --threads 4 -o results/amrfinder/${id}.tsv
done

# 4. CARD + BLAST (validation). Download CARD data first:
#    wget -O data/card/card-data.tar.bz2 https://card.mcmaster.ca/latest/data
#    tar -xvf data/card/card-data.tar.bz2 -C data/card
for f in data/genomes/*.fna; do
  id=$(basename "$f" .fna)
  blastn -query data/card/nucleotide_fasta_protein_homolog_model.fasta -subject "$f" -evalue 1e-10 \
         -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
         -out results/blast/${id}.tsv
done

# 5. Analysis and figures
python scripts/analyze.py
python scripts/locations.py
python scripts/compare_tools.py
