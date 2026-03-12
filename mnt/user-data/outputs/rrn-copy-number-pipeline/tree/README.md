# Reference tree files

This directory should contain the pre-built V3-V4 reference tree and associated files:

| File | Description |
|------|-------------|
| `rrndb_V3V4.treefile` | IQ-TREE2 phylogenetic tree (best model: TVM+F+I+R10) |
| `rrndb_unique_V3V4_filtered_alignment.fasta` | Reference alignment (43,140 seqs, 3,152 columns) |
| `ref.hmm` | HMM profile built from the reference alignment |
| `rrndb_copies_filtered.tsv` | rrn copy numbers per genome ID (rrnDB v5.10) |

> **Note:** These files are not tracked by Git due to their size. See the main README for instructions on how to rebuild them, or download from [Releases](../../releases).

## For use with place_seqs.py

Place all four files into a directory (e.g., `rrndb_ref/`) and pass it to
`place_seqs.py` via `--ref_dir`. PICRUSt2 expects the following naming convention
inside the reference directory:

```
rrndb_ref/
├── ref_seqs.fna          ← reference sequences (FASTA)
├── ref_tree.tre           ← reference tree
├── ref_msa.fasta          ← reference alignment
└── ref_metadata.tsv       ← trait table (rrn copies)
```

Rename or symlink the files accordingly before running `place_seqs.py`.
