# Validation

RASHER's methods are chosen and tested against simulated datasets with a known truth and against public real
data. The validation suite lives in
[`validation/`](https://github.com/duceppemo/RASHER/blob/main/validation/README.md): the simulator
(plastid, linear, rdna, heteroplasmy and NUMT/NUPT scenarios, now as paired 2x150 bp reads), the runner and
the scoring scripts, with a dated record of every validation under `validation/results/`.

The first short-read record — the simulated scenarios with the templated assembly, SPAdes and SKESA, and a
benchmark against GetOrganelle, NOVOPlasty and Snippy on public Illumina data — is in progress and will be
summarised here. In particular it will establish:

- which assembler should be the default, and whether SKESA (conservative; it breaks at repeats such as the
  plastome inverted repeat) earns its place;
- how the templated consensus behaves at low depth and inside repeats with 150 bp pairs (preliminary: reads
  in a long inverted repeat map with ambiguous placement, so low-depth samples call `N` there);
- how NUPT/NUMT reads affect baiting and the consensus (preliminary: nearly all plastome-like nuclear reads
  are baited; the consensus calls `N` rather than taking the nuclear allele);
- timings on real datasets.

The parent tool's long-read validation is documented in
[BACoN's wiki](https://github.com/duceppemo/BACoN/wiki/Validation); its records remain under
`validation/results/` with their original BACoN names. Those results are about Nanopore reads and do not
carry over to RASHER.
