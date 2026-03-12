# rrn copy number assignment with rrnDB

Pipeline for assigning 16S rRNA operon copy numbers (rrn) to ASVs using [rrnDB](https://rrndb.umms.med.umich.edu/) as a reference database, via phylogenetic placement with EPA-ng and maximum parsimony ancestral state reconstruction with PICRUSt2's `hsp.py`.

Copy number is used as a proxy for trophic strategy: organisms with few rrn copies tend to be oligotrophs (slow-growing, efficient in low-nutrient environments), while those with many copies tend to be copiotrophs (rapid response to nutrient pulses) ([Klappenbach et al., 2000](https://doi.org/10.1128/AEM.66.4.1328-1333.2000); [Gao & Wu, 2019](https://doi.org/10.1101/350348)).

---

## Overview

```
ASVs (FASTA)
    │
    ├── hmmalign (align to reference HMM)
    │
    ├── place_seqs.py (PICRUSt2, --min_align 0.8)
    │       └── EPA-ng phylogenetic placement onto rrnDB reference tree
    │
    └── hsp.py (Maximum Parsimony)
            └── rrn_copies_estimated.tsv
```

---

## Requirements

- [PICRUSt2](https://github.com/picrust/picrust2) (≥2.5) — provides `place_seqs.py` and `hsp.py`
- [EPA-ng](https://github.com/Pbdas/epa-ng)
- [HMMER](http://hmmer.org/) — `hmmbuild`, `hmmalign`
- [seqkit](https://bioinf.shenlab.ac.cn/seqkit/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) *(only if building your own tree)*
- [IQ-TREE2](http://www.iqtree.org/) *(only if building your own tree)*

---

## Reference tree

A pre-built reference tree for the **V3-V4 region** is provided in `tree/`. It was constructed from rrnDB v5.10 sequences (see [Building your own reference tree](#building-your-own-reference-tree) below if your amplicons target a different region).

The `tree/` directory contains:

| File | Description |
|------|-------------|
| `rrndb_V3V4.treefile` | IQ-TREE2 phylogenetic tree (best model: TVM+F+I+R10) |
| `rrndb_unique_V3V4_filtered_alignment.fasta` | Reference alignment (43,140 seqs, 3,152 columns) |
| `ref.hmm` | HMM profile built from the reference alignment |
| `rrndb_copies_filtered.tsv` | rrn copy numbers per genome ID from rrnDB v5.10 |

### How the reference tree was built

```bash
# 1. Download rrnDB
#    rrnDB-5.10_16S_rRNA.fasta  (265,316 sequences)
#    rrnDB-5.10.tsv             (copy numbers per genome ID)

# 2. Extract V3-V4 region with HVRLocator, keep one sequence per genome ID
#    → 48,451 sequences

# 3. Filter by minimum length 430 bp
seqkit seq -m 430 rrndb_unique_V3V4.fasta > rrndb_unique_V3V4_filtered.fasta
# → 43,140 sequences, avg length 494.3 bp

# 4. Align with MAFFT
mafft --auto --thread 6 rrndb_unique_V3V4_filtered.fasta \
    > rrndb_unique_V3V4_filtered_alignment.fasta

# 5. Build tree with IQ-TREE2
iqtree2 -s rrndb_unique_V3V4_filtered_alignment.fasta -m MFP -T AUTO \
    --prefix rrndb_V3V4

# 6. Build HMM profile
hmmbuild ref.hmm rrndb_unique_V3V4_filtered_alignment.fasta
```

---

## Running the pipeline

### 1. Align your ASVs to the reference HMM

```bash
hmmalign --trim --dna \
    --mapali rrndb_unique_V3V4_filtered_alignment.fasta \
    --informat FASTA \
    -o query_align.stockholm \
    ref.hmm your_asvs.fasta
```

This produces a joint alignment (reference + query sequences) of 3,152 columns.

### 2. Phylogenetic placement + copy number estimation (PICRUSt2)

```bash
place_seqs.py \
    -s your_asvs.fasta \
    --ref_dir rrndb_ref/ \
    -o placed_seqs.tre \
    --min_align 0.8 \
    --chunk_size 500 \
    --processes 4 \
    --intermediate placement_working/ \
    --verbose

hsp.py \
    -t placed_seqs.tre \
    --observed_trait_table rrndb_copies_filtered.tsv \
    -m mp \
    -o rrn_copies_estimated.tsv \
    -n
```

> **Note:** `--min_align 0.8` excludes ASVs whose aligned fraction with reference sequences is below 80%. This is recommended — ASVs that align poorly introduce noise into the tree topology and degrade parsimony reconstruction for neighboring nodes.

### Output

`rrn_copies_estimated.tsv` contains one row per ASV with the estimated copy number (`metadata_mean_16S_rRNA_Count`) and the NSTI value (`metadata_NSTI`, phylogenetic distance to the nearest reference sequence).

---

## Interpreting the results

Copy numbers can be used to classify trophic strategy following [Klappenbach et al. (2000)](https://doi.org/10.1128/AEM.66.4.1328-1333.2000) and [Gao & Wu (2019)](https://doi.org/10.1101/350348):

| Category | rrn copies |
|----------|-----------|
| Oligotroph | 1–2 |
| Intermediate | 3–4 |
| Copiotroph | ≥5 |

The intermediate category is retained to minimize misclassification in the transition zone between strategies.

**NSTI thresholds** (from PICRUSt2 recommendations):

| NSTI | Reliability |
|------|-------------|
| < 0.15 | High |
| 0.15 – 0.50 | Moderate |
| > 0.50 | Low (interpret with caution) |

High NSTI values are expected for samples from poorly characterized environments (e.g., freshwater systems, soils from underrepresented regions) where rrnDB coverage is limited.

---

## Building your own reference tree

If your sequences target a **different hypervariable region** (e.g., V4, V1-V3), you should build a region-specific reference tree instead of using the V3-V4 tree provided here. The steps are the same — just use HVRLocator to extract the appropriate region from rrnDB sequences before aligning and building the tree (steps 2–6 above).

The minimum length filter in step 3 may need adjustment depending on the expected length of your target region.

---

## Example analysis

See `example_analysis/` for an R Markdown workflow that demonstrates:
- Loading the `hsp.py` output
- Classifying ASVs into trophic strategies
- Evaluating placement quality (NSTI distribution)
- Visualizing rrn copy number distributions by phylum

---

## References

- Klappenbach JA, Dunbar JM, Schmidt TM. 2000. rRNA operon copy number reflects ecological strategies of bacteria. *Appl Environ Microbiol* 66:1328–1333.
- Gao Y, Wu M. 2019. Free-living bacterial communities are mostly dominated by oligotrophs. *bioRxiv* doi:10.1101/350348.
- Douglas GM et al. 2020. PICRUSt2 for prediction of metagenome functions. *Nat Biotechnol* 38:685–688.
- Neufeld JD et al. rrnDB: improved tools for interpreting rRNA gene abundance in bacteria and archaea and a new foundation for future development. *Nucleic Acids Res* 2016.
