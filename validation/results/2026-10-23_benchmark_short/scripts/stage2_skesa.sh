#!/usr/bin/env bash
# SKESA's real-data shot (decision rule from the user: drop -a skesa if never better than SPAdes).
set -euo pipefail
B=/media/8tb/analyses/bacon_benchmark
R=$B/rasher
T=${1:-24}
echo "[$(date +%T)] RASHER on human mt (SKESA)"
/usr/bin/time -v -o $R/mt_skesa.time \
    rasher -r $B/refs/NC_012920.1.fasta -i $R/mt_fastq -o $R/mt_skesa -a skesa -t "$T" -p 4 \
    > $R/mt_skesa.log 2>&1
echo "[$(date +%T)] RASHER on Arabidopsis accessions, plastome, SKESA"
/usr/bin/time -v -o $R/acc_cp_skesa.time \
    rasher -r $B/refs/NC_000932.1.fasta -i $B/accessions/illumina -o $R/acc_cp_skesa -a skesa -t "$T" -p 3 \
    > $R/acc_cp_skesa.log 2>&1
echo "[$(date +%T)] stage 2 done"
