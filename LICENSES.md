# Data Source Licenses

This file catalogs every external data source used by RMS-ISP. Pipeline code is MIT-licensed (`LICENSE`); documentation is CC-BY-4.0. External data sources retain their own licenses; we redistribute only what each license permits, and otherwise ship retrieval scripts that pull from the source.

Default posture: when a license is unclear, treat as non-redistributable. Users fetch raw files themselves through documented retrieval scripts.

## Reference genome and core annotation

| Source | License | Redistributable? | Citation |
|---|---|---|---|
| GRCh38 / GENCODE | LGPL-2.1 | Yes (with attribution) | Frankish et al., Nucleic Acids Res 2021 |
| Ensembl VEP | Apache 2.0 | Yes | McLaren et al., Genome Biol 2016, PMID 27268795 |
| ANNOVAR | Free for academic use | No (academic license required) | Wang et al., Nucleic Acids Res 2010, PMID 20601685 |
| gnomAD v4 | CC-BY 4.0 | Yes (with attribution) | Karczewski et al., Nature 2020, PMID 32461654 |
| ClinVar | Public domain (NCBI) | Yes | Landrum et al., Nucleic Acids Res 2014, PMID 24234437 |
| COSMIC Cancer Gene Census | Free for academic use | No (academic license required) | Sondka et al., Nat Rev Cancer 2018 |
| OncoKB | Free for academic use | No (signup required) | Chakravarty et al., JCO PO 2017 |
| CADD | Free for academic use | No | Kircher et al., Nat Genet 2014, PMID 24487276 |
| REVEL | Free for academic use | No | Ioannidis et al., AJHG 2016, PMID 27666373 |
| AlphaMissense | CC-BY-NC-SA 4.0 | Conditional | Cheng et al., Science 2023, PMID 37733863 |

## Pediatric cancer molecular data

| Source | License / Access | Redistributable? | Notes |
|---|---|---|---|
| TARGET-RT processed (GDC) | Open access | Yes (with attribution) | dbGaP phs000720; raw tier requires DAR |
| TARGET-RT raw | Controlled (dbGaP) | No | Requires DAR; not used in pilot |
| OpenPedCan v15+ | CC-BY 4.0 | Yes | CHOP D3b Center |
| Shern 2014 (GSE66533) | NCBI GEO terms | Yes | Cancer Discov 2014 |
| Shern 2021 supplementary tables | Publisher terms | Per publisher | J Clin Oncol 2021 |
| Wei 2022 scRNA-seq | NCBI GEO terms | Yes | Dev Cell 2022 |
| DeMartino 2023 | NCBI GEO terms | Yes | Cancer Discov 2023 |
| Danielli 2024 | NCBI GEO terms | Yes | Cancer Cell 2024 |
| Gryder 2017-2020 | NCBI GEO/SRA terms | Yes | Cancer Discov 2017, Mol Cell 2019, Nat Genet 2020 |
| MCI (dbGaP phs002790) | Controlled (DAR) | No | Pilot defers MCI access; parallel DAR workstream |
| INSTRuCT | Per consortium | No | Summary stats only for pilot |

## Cell line and PDX resources

| Source | License | Redistributable? | Notes |
|---|---|---|---|
| DepMap Public 24Q2+ | CC-BY 4.0 | Yes | Broad Institute |
| CCLE | CC-BY 4.0 | Yes | Broad Institute |
| PPTC public tier | Per study | Per study | Houghton/Kurmasheva |
| ITCC-P4 public subsets | Per study | Per study | Full panel needs MTA |
| St. Jude CSTN | Per study | Per study | Cell line MTAs available |

## Drug and target databases

| Source | License | Redistributable? | Notes |
|---|---|---|---|
| DGIdb v5 | MIT | Yes | API access |
| OpenTargets | CC0 1.0 (most data) | Yes | Verify per dataset |
| Molecular Targets Platform | NCI/CCDI terms | Per terms | API access |
| ChEMBL v34+ | CC-BY-SA 3.0 | Yes (with attribution) | EBI |
| DrugBank Open Data | CC-BY-NC 4.0 | Conditional (non-commercial) | |
| ClinicalTrials.gov | Public domain | Yes | Official API |
| CMap / LINCS L1000 | Per platform | Per dataset | clue.io / lincsportal |

## Protein structure resources

| Source | License | Redistributable? | Notes |
|---|---|---|---|
| AlphaFold DB | CC-BY 4.0 | Yes | Jumper et al., Nature 2021 |
| RCSB PDB | Public domain | Yes | |
| Boltz-1 | MIT | Yes (model + weights) | Wohlwend et al., MIT |
| Chai-1 | Apache 2.0 (code); CC-BY-NC 4.0 (weights) | Conditional | Chai Discovery |
| ESMFold | MIT | Yes | Meta AI |

## Pipeline and computational tools

| Tool | License | Notes |
|---|---|---|
| Nextflow | Apache 2.0 | Di Tommaso et al., Nat Biotechnol 2017 |
| Docker | Apache 2.0 | |
| Singularity / Apptainer | BSD-3 | |
| FoldX | Free for academic use | Schymkowitz et al., Nucleic Acids Res 2005 |
| AutoDock Vina | Apache 2.0 | |
| STRING v12 | CC-BY 4.0 | Szklarczyk et al., NAR 2023 |
| clusterProfiler | Artistic 2.0 | Yu et al., OMICS 2012 |
| Cytoscape | LGPL 2.1 | Shannon et al., Genome Res 2003 |
| MSigDB | CC-BY 4.0 | Liberzon et al., Cell Syst 2015 |
| Reactome | CC-BY 4.0 | Gillespie et al., NAR 2022 |

## Adding a source

When integrating a new source:

1. Confirm license at the source URL on the date of integration.
2. Add a row to the appropriate table above.
3. Note the access date in `corpus_manifest.json`.
4. If unsure, treat as non-redistributable: include a retrieval script, not the data.
