#!/usr/bin/env bash
# RASHER 0.1.0 benchmark runs on the public data already on disk (no installs, no downloads):
#   1. human mt: 32 1KGP chrM Illumina BAMs -> paired fastq -> RASHER templated vs rCRS (NC_012920.1)
#   2. Arabidopsis accessions: 24 HiSeq PE runs vs plastome NC_000932.1, templated and SPAdes
# Usage: bash run_rasher.sh [THREADS]
set -euo pipefail
B=/media/8tb/analyses/bacon_benchmark
R=$B/rasher
T=${1:-24}

echo "[$(date +%T)] converting chrM BAMs to paired fastq"
mkdir -p $R/mt_fastq
for bam in $B/human_mt/illumina/*.chrM.bam; do
    s=$(basename "$bam" .chrM.bam)
    [ -f "$R/mt_fastq/${s}_R1.fastq.gz.ok" ] && continue
    samtools collate -u -O -@ 4 "$bam" | samtools fastq -@ 4 -1 "$R/mt_fastq/${s}_R1.fastq.gz" \
        -2 "$R/mt_fastq/${s}_R2.fastq.gz" -0 /dev/null -s /dev/null -n -
    touch "$R/mt_fastq/${s}_R1.fastq.gz.ok"
done

echo "[$(date +%T)] RASHER on human mt (templated)"
/usr/bin/time -v -o $R/mt_samtools.time \
    rasher -r $B/refs/NC_012920.1.fasta -i $R/mt_fastq -o $R/mt_samtools -t "$T" -p 4 \
    > $R/mt_samtools.log 2>&1

echo "[$(date +%T)] RASHER on Arabidopsis accessions, plastome, templated"
/usr/bin/time -v -o $R/acc_cp_samtools.time \
    rasher -r $B/refs/NC_000932.1.fasta -i $B/accessions/illumina -o $R/acc_cp_samtools -t "$T" -p 3 \
    > $R/acc_cp_samtools.log 2>&1

echo "[$(date +%T)] RASHER on Arabidopsis accessions, plastome, SPAdes"
/usr/bin/time -v -o $R/acc_cp_spades.time \
    rasher -r $B/refs/NC_000932.1.fasta -i $B/accessions/illumina -o $R/acc_cp_spades -a spades -t "$T" -p 3 \
    > $R/acc_cp_spades.log 2>&1

echo "[$(date +%T)] all done"
