#!/bin/bash

# ----------------------------------------------------------------------------------------
# --- Generate dataset for Snakemake test (subset of Ag1000G data)
# ----------------------------------------------------------------------------------------

# Move into working directory
# cd $SCRATCH/uvri_sm

mkdir -p data

# --- Download Ugandan sample list:

URL_PRE=https://storage.googleapis.com/vo_agam_release/v3/metadata/general

wget $URL_PRE/AG1000G-UG/samples.meta.csv \
    -O data/UG.samples.meta.csv

# --- Download the catalog of run accessions and alignments

wget https://storage.googleapis.com/vo_agam_release/v3/metadata/ena_runs.csv \
    -O data/ena_runs.csv

wget https://storage.googleapis.com/vo_agam_release/v3/metadata/ena_alignments.csv \
    -O data/ena_alignments.csv

# --- Download FASTQ reads

#    FTP=ftp://ftp.sra.ebi.ac.uk/vol1
#
#    for IND in `cut -f1 -d','  data/UG.samples.meta.csv | grep -v "sample"`; do
#
#        SEQ_IDS=($(grep $IND data/ena_runs.csv | cut -f 2 -d','))
#
#        # For this, only download first file
#
#        THIS_SEQ_ID=${SEQ_IDS[1]}
#
#        PREFIX=$FTP/fastq/${THIS_SEQ_ID:0:6}/$THIS_SEQ_ID/${THIS_SEQ_ID}
#        wget ${PREFIX}_1.fastq.gz -O data/${IND}_1.fastq
#        wget ${PREFIX}_1.fastq.gz -O data/${IND}_2.fastq
#
#    done

# --- Download BAM files

# For now, only download first 20

for IND in `cut -f1 -d','  data/UG.samples.meta.csv | grep -v "sample" | sed -n '1,20p'`; do

    SEQ_ID=`grep $IND data/ena_alignments.csv | cut -f 2 -d',' | sed 's/\r$//'`

    # Download BAM files

    CMD="cd `pwd`; \
        wget ftp://ftp.sra.ebi.ac.uk/vol1/${SEQ_ID:0:6}/$SEQ_ID/$IND.bam \
            -O data/$IND.bam"

    echo $CMD | qsub -V -M cxb585@psu.edu -A open -m abe \
        -l nodes=1:ppn=1,walltime=4:00:00 -l feature=rhel7 -N $IND

done

# --- Pull out gene of interest (Vgsc / kdr)

# Gene info from VEuPathDB:
#    AgamP4_2L   VEuPathDB   protein_coding_gene 2358158 2431617 .   +   .
#    ID=AGAP004707;description=voltage-gated sodium channel

COORDS="2L:2358158-2431617"

module load samtools/1.2

for IND in `cut -f1 -d','  data/UG.samples.meta.csv | grep -v "sample" | sed -n '1,20p'`; do

    samtools index data/$IND.bam
    samtools view -h data/$IND.bam "$COORDS" > data/$IND.subset.bam
done

# --- Convert to FASTQ format

module load samtools/1.2
module load bedtools

for IND in `cut -f1 -d','  data/UG.samples.meta.csv | grep -v "sample" | sed -n '1,20p'`; do

    samtools sort -n data/$IND.subset.bam data/$IND.subset.sort
    bedtools bamtofastq -i data/$IND.subset.sort.bam \
        -fq data/$IND.R1.fastq -fq2 data/$IND.R2.fastq
done

# --- Also download reference genome, AgamP4, but just chromosome 2L

# Script to download Anopheles gambiae genome, PEST, AgamP4

mkdir -p genomes/AgamP4
cd genomes/AgamP4

GENOME_FA=AgamP4_2L.fa

AGAM_URL=ftp://ftp.ensemblgenomes.org/pub/metazoa/release-51/fasta
AGAM_URL=${AGAM_URL}/anopheles_gambiae/dna/
AGAM_URL=${AGAM_URL}/Anopheles_gambiae.AgamP4.dna_sm.chromosome.2L.fa.gz

wget $AGAM_URL \
    -O ${GENOME_FA}.gz
gunzip ${GENOME_FA}.gz

cd ../..

# --- And index the reference "genome" (actually just chromosome arm 2L)

module load bwa

bwa index genomes/AgamP4/AgamP4_2L.fa

exit
