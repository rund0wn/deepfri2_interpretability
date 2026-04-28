# GPSite_dataset

Raw ligand splits and processed exports for the GPSite binding-site dataset.

## Provenance

The archive was downloaded from Zenodo record **10845362** and unzipped:

https://zenodo.org/records/10845362/files/GPSite_dataset.zip?download=1

Ligand-specific train/test text files live in [`GPSite_rawdata/`](GPSite_rawdata/). Each file is named `{LIGAND}_{train|test}.txt`. Every entry is three lines: a `>{PDB}_{chain}` header, the amino-acid sequence, and a same-length `0`/`1` mask (`1` = binding site for that ligand).

## Processed outputs (`GPSite_processed/`)

Regenerate with:

```bash
cd GPSite_dataset
uv run python scripts/build_gpsite_processed.py
```

Optional flags: `--raw-dir` (default: `GPSite_rawdata`), `--out-dir` (default: `GPSite_processed`), `--biolip-nr` (default: [`BioLip_rawdata/BioLiP_nr.txt`](BioLip_rawdata/BioLiP_nr.txt); columns per [`BioLip_rawdata/readme.txt`](BioLip_rawdata/readme.txt)).

| File | Description |
|------|-------------|
| `gpsite_master.fasta` | Unique `>{protein_id}` sequences merged across all raw files (sequence must agree if the same ID appears twice). |
| `gpsite_binding_records.pkl` | `list[dict]` with `protein_id`, `ligand`, `split`, `sequence`, `sequence_length`, `binding_positions` (1-based indices where the mask is `1`). One dict per raw row (ligand- and split-specific). |
| `gpsite_master_lt1020.fasta` | Same as master FASTA, only proteins with `len(sequence) < 1020`. |
| `gpsite_binding_records_lt1020.pkl` | Same as full pickle, plus `uniprot_id` (or `null`): PDB+chain matched to BioLiP column 18 UniProt on [`BioLiP_nr.txt`](BioLip_rawdata/BioLiP_nr.txt). Only lt1020 rows. |
| `gpsite_ligand_counts.tsv` | Per-ligand: `n_records_before`, `n_records_after`, and `n_proteins_uniprot` (distinct lt1020 `protein_id` with a BioLiP UniProt hit for that ligand). |
| `gpsite_ligand_summary.md` | Markdown table of the same columns; total row omits summed UniProt (not additive across ligands); footer gives distinct mapped lt1020 proteins globally. |

## Scripts

| Script | Role |
|--------|------|
| [`scripts/build_gpsite_processed.py`](scripts/build_gpsite_processed.py) | Builds everything in `GPSite_processed/` from `GPSite_rawdata/`. |


