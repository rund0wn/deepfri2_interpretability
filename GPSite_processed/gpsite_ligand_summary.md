# GPSite ligand record counts

Rows = train + test entries per ligand. After = `len(sequence) < 1020`. `n_proteins_uniprot` = distinct `protein_id` in the lt1020 split with a UniProt from BioLiP `(pdb, chain)` lookup.

| ligand | n_records_before | n_records_after | n_proteins_uniprot |
|--------|------------------|-----------------|--------------------|
| ATP | 426 | 416 | 78 |
| CA | 1737 | 1737 | 1045 |
| DNA | 807 | 786 | 458 |
| HEM | 224 | 223 | 48 |
| MG | 1964 | 1964 | 1005 |
| MN | 604 | 604 | 300 |
| PEP | 1486 | 1469 | 768 |
| PRO | 710 | 710 | 80 |
| RNA | 1035 | 1003 | 405 |
| ZN | 1857 | 1857 | 1262 |
| **total** | **10850** | **10769** | **—** |

Distinct lt1020 `protein_id` with UniProt (BioLiP): **5039** / **10150** unique proteins under length cap.
