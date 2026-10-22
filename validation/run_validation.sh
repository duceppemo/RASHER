#!/usr/bin/env bash
# Run RASHER on the simulated scenarios with every assembler and SNP method, then score against the truth.
#
# Usage: bash validation/run_validation.sh OUTPUT_FOLDER [THREADS] [all|new]
#
# Run it in the RASHER conda environment (environment.yml). The simulated reads are generated first (fixed seed:
# the same files every time); each run writes into OUTPUT_FOLDER/<scenario>/<run> and the scores go to
# OUTPUT_FOLDER/SCORES_<scenario>.md and .json. See validation/README.md.
# all (the default): plastid, linear and rdna, with every assembler and SNP method; new: heteroplasmy and numt,
# with the templated assembly and Flye (SKA2).
set -euo pipefail

if [ ! -x /usr/bin/time ] || ! realpath -m / >/dev/null 2>&1; then
    echo "error: needs GNU time (/usr/bin/time; Debian/Ubuntu package 'time') and GNU realpath" >&2
    exit 1
fi
out=$(realpath -m "${1:?output folder}")
threads=${2:-16}
group=${3:-all}
here=$(cd "$(dirname "$0")" && pwd)
case $group in
    all) scenarios=(plastid linear rdna) ;;
    new) scenarios=(heteroplasmy numt) ;;
    *) echo "error: the third argument is all or new, not '$group'" >&2; exit 1 ;;
esac

if [ -d "$out/${scenarios[0]}" ]; then  # Resumed runs would skip steps and make the timings meaningless
    echo "error: $out already has results; use a new folder (simulated reads can be reused: copy $out/sim)" >&2
    exit 1
fi
mkdir -p "$out"
if [ ! -f "$out/sim/${scenarios[-1]}/truth.json" ]; then
    python "$here/simulate.py" "$out/sim" --scenario "$group"
fi

rasher_run() {  # rasher_run SCENARIO RUN [rasher options...]
    local scenario=$1 run=$2
    shift 2
    local dir="$out/$scenario/$run"
    mkdir -p "$out/$scenario"
    echo "[$(date +%T)] $scenario / $run: rasher $*"
    # One line per rasher call, appended: the first is the full run, later ones only the steps that changed.
    /usr/bin/time -a -f "%e s, %M kB: rasher $*" -o "$dir.time" \
        rasher -r "$out/sim/$scenario/reference.fasta" -i "$out/sim/$scenario/reads" -o "$dir" \
        -t "$threads" -p 4 "$@" >> "$dir.log" 2>&1 || echo "    FAILED (see $dir.log)"
}

parsnp_default() {  # Parsnp without -c (its default genome filtering), on existing assemblies
    local scenario=$1 dir="$out/$1/flye/4_compared/parsnp_default"
    echo "[$(date +%T)] $scenario / flye: parsnp without -c"
    rm -rf "$dir"
    mkdir -p "$(dirname "$dir")"
    local asm=("$out/$scenario/flye/3_assembled/all_assemblies/"*.fasta)
    parsnp -r "$out/$scenario/flye/reference.fasta" -d "${asm[@]}" -o "$dir" -p "$threads" > "$dir.log" 2>&1 || true
    if [ -f "$dir/parsnp.ggr" ]; then
        harvesttools -i "$dir/parsnp.ggr" -S "$dir/snps.raw.fasta" >> "$dir.log" 2>&1
        python - "$dir" "$out/$scenario/flye" <<'PY'
import sys
from pathlib import Path
from rasher.compare import clean_alignment, snp_distances, write_distances
d, run = Path(sys.argv[1]), Path(sys.argv[2])
rename = {p.name: p.stem for p in (run / "3_assembled" / "all_assemblies").glob("*.fasta")}
rename.update({"reference.fasta.ref": "Reference", "reference.fasta": "Reference"})
recs = clean_alignment(d / "snps.raw.fasta", d / "parsnp.snps.fasta", rename)
write_distances(d / "snp_distances.tsv", *snp_distances(recs))
PY
    fi
}

{
    echo "# Validation run $(date -Iseconds)"
    echo "rasher: $(rasher --version)"
    for tool in minimap2 filtlong flye myloasm samtools ska harvesttools; do
        echo "$tool: $("$tool" --version 2>&1 | head -1)"
    done
    echo "parsnp: $(parsnp --version 2>&1 | grep -o 'Parsnp [0-9.]*' | head -1)"
    echo "FastTree: $(FastTree -expert 2>&1 | grep -o 'FastTree [0-9.]*' | head -1)"
} > "$out/versions.txt"

for scenario in "${scenarios[@]}"; do
    # The defaults (--snp-method auto): templated assembly compared with minimap2, Flye with SKA2
    rasher_run "$scenario" samtools
    rasher_run "$scenario" flye -a flye
    # The other method on the same assemblies (the assembly step is reused from the checkpoint)
    rasher_run "$scenario" samtools --snp-method ska
    rasher_run "$scenario" flye -a flye --snp-method minimap2
    if [ "$group" = all ]; then
        rasher_run "$scenario" myloasm -a myloasm
        # Other comparisons of the Flye assemblies (the assembly step is reused from the checkpoint).
        rasher_run "$scenario" flye -a flye --snp-method parsnp
        rasher_run "$scenario" flye -a flye --ska-min-freq 0.5
        parsnp_default "$scenario"
    fi
    python "$here/score.py" "$out/sim/$scenario" "$out/$scenario"/*/ --markdown "$out/SCORES_$scenario.md" \
        --json "$out/SCORES_$scenario.json" > /dev/null
    echo "Scores: $out/SCORES_$scenario.md"
done
