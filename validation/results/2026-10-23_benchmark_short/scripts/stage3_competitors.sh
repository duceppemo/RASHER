#!/usr/bin/env bash
# Competitors on the same inputs: GetOrganelle, NOVOPlasty (assemblies), Snippy (SNPs).
# Envs: bench_getorganelle, bench_novoplasty, bench_snp. Sequential per sample, threads within tools.
set -uo pipefail   # no -e: one sample's failure must not stop the benchmark; failures are visible in the logs
B=/media/8tb/analyses/bacon_benchmark
R=$B/rasher
T=${1:-16}
source ~/miniconda3/etc/profile.d/conda.sh

acc_samples() {  # name<TAB>R1<TAB>R2
    for r1 in $B/accessions/illumina/*_1.fastq.gz; do
        printf '%s\t%s\t%s\n' "$(basename "$r1" _1.fastq.gz)" "$r1" "${r1%_1.fastq.gz}_2.fastq.gz"
    done
}
mt_samples() {
    for r1 in $R/mt_fastq/*_R1.fastq.gz; do
        printf '%s\t%s\t%s\n' "$(basename "$r1" _R1.fastq.gz)" "$r1" "${r1%_R1.fastq.gz}_R2.fastq.gz"
    done
}

echo "[$(date +%T)] GetOrganelle"
conda activate bench_getorganelle
mkdir -p $R/getorganelle/acc $R/getorganelle/mt
while IFS=$'\t' read -r s r1 r2; do
    [ -f "$R/getorganelle/acc/$s/.ok" ] && continue
    /usr/bin/time -v -o "$R/getorganelle/acc/$s.time" \
        get_organelle_from_reads.py -1 "$r1" -2 "$r2" -o "$R/getorganelle/acc/$s" -F embplant_pt -t "$T" \
        > "$R/getorganelle/acc/$s.log" 2>&1 && touch "$R/getorganelle/acc/$s/.ok"
done < <(acc_samples)
while IFS=$'\t' read -r s r1 r2; do
    [ -f "$R/getorganelle/mt/$s/.ok" ] && continue
    /usr/bin/time -v -o "$R/getorganelle/mt/$s.time" \
        get_organelle_from_reads.py -1 "$r1" -2 "$r2" -o "$R/getorganelle/mt/$s" -F animal_mt -t "$T" \
        > "$R/getorganelle/mt/$s.log" 2>&1 && touch "$R/getorganelle/mt/$s/.ok"
done < <(mt_samples)
conda deactivate

echo "[$(date +%T)] NOVOPlasty"
conda activate bench_novoplasty
mkdir -p $R/novoplasty/acc $R/novoplasty/mt
np_config() {  # name seed range_min range_max out r1 r2
    cat << CFG
Project:
-----------------------
Project name          = $1
Type                  = $([ "$3" -lt 100000 ] && echo mito || echo chloro)
Genome Range          = $3-$4
K-mer                 = 33
Max memory            = 60
Extended log          = 0
Save assembled reads  = no
Seed Input            = $2
Extend seed directly  = no
Reference sequence    =
Variance detection    =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = 150
Insert size           = 350
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = $6
Reverse reads         = $7

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = no
Output path           = $5/
CFG
}
while IFS=$'\t' read -r s r1 r2; do
    out=$R/novoplasty/acc/$s
    [ -f "$out/.ok" ] && continue
    mkdir -p "$out"
    np_config "$s" "$B/refs/NC_000932.1.fasta" 120000 180000 "$out" "$r1" "$r2" > "$out/config.txt"
    /usr/bin/time -v -o "$out.time" NOVOPlasty.pl -c "$out/config.txt" > "$out.log" 2>&1 && touch "$out/.ok"
done < <(acc_samples)
while IFS=$'\t' read -r s r1 r2; do
    out=$R/novoplasty/mt/$s
    [ -f "$out/.ok" ] && continue
    mkdir -p "$out"
    np_config "$s" "$B/refs/NC_012920.1.fasta" 15000 18000 "$out" "$r1" "$r2" > "$out/config.txt"
    /usr/bin/time -v -o "$out.time" NOVOPlasty.pl -c "$out/config.txt" > "$out.log" 2>&1 && touch "$out/.ok"
done < <(mt_samples)
conda deactivate

echo "[$(date +%T)] Snippy"
conda activate bench_snp
mkdir -p $R/snippy/acc $R/snippy/mt
while IFS=$'\t' read -r s r1 r2; do
    [ -f "$R/snippy/acc/$s/.ok" ] && continue
    /usr/bin/time -v -o "$R/snippy/acc/$s.time" \
        snippy --outdir "$R/snippy/acc/$s" --ref $B/refs/NC_000932.1.fasta --R1 "$r1" --R2 "$r2" \
        --cpus "$T" --force > "$R/snippy/acc/$s.log" 2>&1 && touch "$R/snippy/acc/$s/.ok"
done < <(acc_samples)
(cd $R/snippy/acc && snippy-core --ref $B/refs/NC_000932.1.fasta */ > core.log 2>&1)
while IFS=$'\t' read -r s r1 r2; do
    [ -f "$R/snippy/mt/$s/.ok" ] && continue
    /usr/bin/time -v -o "$R/snippy/mt/$s.time" \
        snippy --outdir "$R/snippy/mt/$s" --ref $B/refs/NC_012920.1.fasta --R1 "$r1" --R2 "$r2" \
        --cpus "$T" --force > "$R/snippy/mt/$s.log" 2>&1 && touch "$R/snippy/mt/$s/.ok"
done < <(mt_samples)
(cd $R/snippy/mt && snippy-core --ref $B/refs/NC_012920.1.fasta */ > core.log 2>&1)
conda deactivate
echo "[$(date +%T)] stage 3 done"
