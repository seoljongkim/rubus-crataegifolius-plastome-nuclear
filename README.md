# Analysis code — *Rubus crataegifolius* plastome and nuclear resequencing

Code accompanying **"Discordant chloroplast and nuclear divergence among eight *Rubus crataegifolius*
(Rosaceae) cultivars."** *Genes* **2026**, *17*(10), 1287. https://doi.org/10.3390/genes17101287

Eight cultivars were sequenced on Illumina NovaSeq 6000 (151 bp paired-end, 19.9–23.8 M read pairs each).
Plastomes: GenBank **PZ143672–PZ143679**. Raw sequencing reads are available from the corresponding
author on reasonable request.

**This repository contains code only.** Input and intermediate data files are not included; every
script expects the directory layout of the original analysis (see *Paths*).

## Pipeline

| Step | Directory | What it does |
|---|---|---|
| 1 | `01_assembly_annotation/` | Plastome standardisation (SSC orientation, rotation to *trnH*-GUG) and GenBank submission files |
| 2 | `02_cp_analysis/` | Variant detection from the MAFFT alignment (`structural_variant_analysis.py`), SSR (MISA, conservative thresholds), tandem and dispersed repeats, RSCU, IR boundaries |
| 3 | `03_selection/` | Length-normalised substitution rates and Nei–Gojobori Ka/Ks, within cultivars and against three congeners |
| 4 | `04_phylogeny/` | ML tree figure from the IQ-TREE output (21 taxa, TVM+F+R3) |
| 5 | `05_nuclear/` | Reference comparison, mapping and joint SNP calling, marker panel selection, distance robustness checks; heterozygosity QC / strict SNP set and panel re-check added at revision (`08_`, `09_`) |
| 6 | `06_figures/` | Figure generation |
| 7 | `07_its_leafy/` | nrDNA ITS assembly, read-level check and ML tree figure. The LEAFY intron 2 check reported in the Discussion was run interactively and has no script here; its inputs and alignments are described in Supplementary Note 2 |

## Cultivar codes

Sample codes map to cultivars as follows: `DAHK_2`/`DAHE` = 'Dongakheuk', `DAH_B1`/`DAHO` =
'Dongakhong', `DB_A1`/`DB` = 'Dabok', `DS_A1`/`DS` = 'Daesung', `GDCH_1`/`GDC` = 'Geumdongcheon',
`GDM_A1`/`GDM` = 'Goldmoon', `GDW_5`/`GDW` = 'Geumdongwang', `HGB_A1`/`HGB` = 'Hwangeumboll'.

Two further accessions, `GB048` and `GB066`, appear in the assembly, annotation and repeat scripts
because ten plastomes were assembled initially. **They were excluded from every analysis reported in
the paper** and are filtered out downstream (see `03_selection/05_kaks.py` and `06_figures/`).

## Paths

Every script resolves its input and output locations from a single base directory:

```sh
export RUBUS_BASE=/path/to/your/analysis/root
```

If `RUBUS_BASE` is unset, the base defaults to the directory three levels above the script, which is
correct when this repository sits inside the analysis root (`<root>/code_release/<step>/<script>`).
The expected layout under the base is `standardized/`, `misa_output/`, `trf_output/`,
`reputer_output/`, `rscu_output/`, `irscope_output/`, `variant_output/`, `kaks_output/`,
`bankit/`, `bankit_gb/`, `phylogeny21/`, `its_analysis/`, `trimmed/`, `figures/` and
`home_wsl_package/`.

One exception: `05_nuclear/06_marker_selection.py` uses paths relative to the working directory and
must be run from inside `home_wsl_package/nuclear_snp/results/`.

## External software

Trimmomatic (`ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`),
GetOrganelle v1.7.7.1, GeSeq, MISA, Tandem Repeats Finder v4.10, BLAST+ v2.17.0, MAFFT v7,
IQ-TREE v3.0.1, bwa-mem2 v2.3, SAMtools v1.24, BCFtools v1.24, PLINK v1.90 and v2.00,
Biopython v1.86, SciPy.

## Notes on reproducing the nuclear analysis

- Reference: *Rubus idaeus* `GCA_965153425.1` (drRubIdae1.hap1.1). No chromosome-level assembly exists
  for *R. crataegifolius*; this reference gave the higher mapping rate of the two candidates tested
  (`05_nuclear/03a_ref_compare.sh`).
- Joint calling is split into 5 Mb windows and run in parallel (`03b_nuclear_snp_parallel.sh`); the
  windows are concatenated afterwards. **Do not pass `--threads` to `bcftools concat`** — in our
  environment it produced a CRC32-corrupted output.
- SNPs at which all eight cultivars carry the same non-reference genotype are fixed interspecific
  differences and are removed (`plink2 --mac 1`) before distance estimation.
- 1 − IBS distances depend on heterozygosity; `07_distance_robustness.py` repeats the analysis on
  homozygous-only sites and checks for confounding by depth and mapping rate.
- Revision (Genes R1): `08_het_qc.py` tabulates heterozygous calls by mapping quality, depth and allele
  balance and writes the strict SNP set (MQ ≥ 50, per-sample depth 8–40×, allele balance 0.2–0.8, no
  missing genotype); `09_panel_strict_check.py` re-checks the marker panel against that set. RcNu-04 of
  the original submission failed the check and was replaced by OZ232403.1:12,029,618 (same genotype
  pattern); `06_marker_selection.py` therefore reproduces the original panel, not the revised one.

## Citation and archive

Archived release: **v1.0.0**, Zenodo DOI [10.5281/zenodo.23364695](https://doi.org/10.5281/zenodo.23364695)
(concept DOI for all versions: 10.5281/zenodo.23364694). Code is released under the MIT License.

Article:

Kim, S.-J.; Kim, Y.-Y.; Kwon, S.-H. Discordant Chloroplast and Nuclear Divergence Among Eight *Rubus crataegifolius*
(Rosaceae) Cultivars. *Genes* **2026**, *17*, 1287. https://doi.org/10.3390/genes17101287

Code:

Kim, S.-J.; Kim, Y.-Y.; Kwon, S.-H. Analysis code: Discordant chloroplast and nuclear divergence among eight
*Rubus crataegifolius* (Rosaceae) cultivars. Zenodo, 2026. https://doi.org/10.5281/zenodo.23364695
