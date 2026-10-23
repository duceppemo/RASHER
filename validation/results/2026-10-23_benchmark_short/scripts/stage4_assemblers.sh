#!/usr/bin/env bash
# More SKESA-vs-SPAdes evidence before fixing the de novo recommendation (user request 2026-10-23):
# 1. human mt SPAdes (the missing direct comparison on dataset 2)
# 2. Arabidopsis MITOGENOME (NC_037304.1, ~367 kb, repeat-rich -- the hard case): templated, SPAdes, SKESA
# 3. depth sweep on the already-baited plastome reads (input = acc_cp_samtools/1_extracted, so baiting
#    is near-instant): SPAdes and SKESA at --target-depth 30 and 15 (100 is the existing runs)
# Usage: bash stage4_assemblers.sh [THREADS]
set -euo pipefail
B=/media/8tb/analyses/bacon_benchmark
R=$B/rasher
T=${1:-24}

echo "[$(date +%T)] RASHER on human mt (SPAdes)"
/usr/bin/time -v -o $R/mt_spades.time \
    rasher -r $B/refs/NC_012920.1.fasta -i $R/mt_fastq -o $R/mt_spades -a spades -t "$T" -p 4 \
    > $R/mt_spades.log 2>&1

for a in samtools spades skesa; do
    echo "[$(date +%T)] RASHER on Arabidopsis accessions, mitogenome, $a"
    /usr/bin/time -v -o $R/acc_mt_${a}.time \
        rasher -r $B/refs/NC_037304.1.fasta -i $B/accessions/illumina -o $R/acc_mt_${a} -a "$a" -t "$T" -p 3 \
        > $R/acc_mt_${a}.log 2>&1
done

for d in 30 15; do
    for a in spades skesa; do
        echo "[$(date +%T)] depth sweep: plastome, $a, target-depth $d"
        /usr/bin/time -v -o $R/acc_cp_${a}_d${d}.time \
            rasher -r $B/refs/NC_000932.1.fasta -i $R/acc_cp_samtools/1_extracted -o $R/acc_cp_${a}_d${d} \
            -a "$a" --target-depth "$d" -t "$T" -p 6 \
            > $R/acc_cp_${a}_d${d}.log 2>&1
    done
done
echo "[$(date +%T)] stage 4 done"
