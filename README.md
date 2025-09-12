# Nitrogen Responsive Genes in Sugarcane (NRGSC)

[![DOI:10.6084/m9.figshare.30100633](https://img.shields.io/badge/dataset-figshare-orange)](https://figshare.com/projects/NitrogeGene_co-expression_networks_highlight_key_nodes_associated_with_ammonium_nitrate_response_in_sugarcane/242225)  
[![NCBI BioProject: PRJNA1176579](https://img.shields.io/badge/data-NCBI%20BioProject%20PRJNA1176579-blue)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1176579)

## Overview

This repository accompanies the study:

**Gene co-expression networks highlight key nodes associated with ammonium nitrate response in sugarcane**  
Jorge Mario Muñoz-Pérez, Denis Bassi, Lucia Mattiello, Marcelo Menossi, Diego Mauricio Riaño-Pachón  
([Contact: diego.riano@cena.usp.br](mailto:diego.riano@cena.usp.br))

This resource investigates the molecular mechanisms underlying nitrogen use efficiency (NUE) in sugarcane by integrating RNA-seq, transcriptome assembly, differential expression, and gene co-expression network analysis. The repository provides all scripts and pipelines for annotation, differential expression analysis, and gene network construction.

---

## Repository Structure

```
NRGSC/
│
├── annotation/              # Transcriptome annotation
│   ├── blast/               # BLAST-based annotation
│   │   ├── code/
│   │   │   ├── fasta-splitter.pl
│   │   │   ├── run_blastParts.sh
│   │   │   └── split_reference.sh
│   │   ├── data/
│   │   └── results/
│   │ 
│   └── transdecoder/
│       └── run_transcoder.sh
│
├── diff_expression/         # Differential expression analysis
│   └── code/
│       ├── DESEQ_DEA_GENE.r
│       ├── DESEQ_DEA_TRANSCRIPTS.r
│       ├── make_full_matrixes_for_trinotate_web.r
│       ├── pheatmap_by_groups.r
│       └── up_set_plot.r
│
├── networks/                # Gene co-expression network analysis
│   └── code/
│       ├── corr_modules.r
│       ├── format_modules.sh
│       ├── get_adjTOlist.r
│       ├── get_hub.r
│       ├── get_tf.sh
│       ├── GO_DEG.r
│       ├── inflation_plot.r
│       ├── make_graph_diff_corrected.r
│       ├── make_heatmap_for_tf.r
│       ├── make_heatmaps_for_all_modules.r
│       ├── make_heatmaps_for_modules.r
│       ├── make_histogram.r
│       ├── make_hp_module_20.r
│       ├── run_graph.sh
│       ├── run_himcl2.sh
│       ├── run_mcl.sh
│       ├── TF_hipergeometric.r
│       └── topGO.r
│
└── README                   # (This file)
```

---

## Workflow Description

### 1. **Transcriptome Annotation (`annotation/`)**

#### a. **ORF Prediction**
- `transdecoder/run_transcoder.sh`: Detects candidate open reading frames (ORFs) in assembled transcripts using TransDecoder.

#### b. **BLAST Annotation**
- `blast/code/fasta-splitter.pl`: Splits large FASTA files into manageable chunks for parallel BLAST.
- `blast/code/split_reference.sh`: Wrapper for splitting the transcriptome or predicted peptides into parts for annotation.
- `blast/code/run_blastParts.sh`: Runs BLASTX and BLASTP in parallel on the split files for high-throughput annotation.

---

### 2. **Differential Expression Analysis (`diff_expression/`)**

#### a. **Gene- and Transcript-level Analysis**
- `code/DESEQ_DEA_GENE.r`: Performs differential expression at the gene level using DESeq2. Generates result tables, MA plots, volcano plots, and filtered DEG lists.
- `code/DESEQ_DEA_TRANSCRIPTS.r`: Performs differential expression at the transcript (isoform) level.

#### b. **Downstream Processing and Visualization**
- `code/make_full_matrixes_for_trinotate_web.r`: Builds complete expression matrices for integration with Trinotate web reports.
- `code/pheatmap_by_groups.r`: Creates heatmaps of expression patterns by group.
- `code/up_set_plot.r`: Generates UpSet plots to visualize intersections of DEG sets across conditions.

---

### 3. **Gene Co-Expression Network Analysis (`networks/`)**

#### a. **Network Construction and Module Detection**
- `code/make_graph_diff_corrected.r`: Constructs the gene co-expression network based on correlation matrices.
- `code/run_graph.sh`: Runs the main network analysis pipeline.
- `code/run_himcl2.sh`, `run_mcl.sh`: Run the Markov Clustering Algorithm (MCL) and/or HiMCL for module identification.

#### b. **Module Characterization and Visualization**
- `code/format_modules.sh`: Formats MCL output into usable modules.
- `code/corr_modules.r`: Correlates modules (eigengenes) with external traits (genotype, nitrogen, metabolites).
- `code/get_adjTOlist.r`: Converts adjacency matrices to edge lists.
- `code/get_hub.r`: Identifies hub genes within modules.
- `code/GO_DEG.r`, `topGO.r`: Performs GO enrichment analysis of DEGs and network modules.
- `code/make_heatmaps_for_modules.r`, `make_heatmap_for_tf.r`, `make_heatmaps_for_all_modules.r`: Generates heatmaps for candidate modules, transcription factor modules, and all modules.
- `code/make_histogram.r`, `inflation_plot.r`: Visualizes module sizes and MCL inflation parameter effects.
- `code/make_hp_module_20.r`: Dedicated script for visualizing and analyzing the MYB-enriched module 20.
- `code/get_tf.sh`, `TF_hipergeometric.r`: Identifies transcription factors and tests for TF enrichment in modules.

---

## How to Use

1. **Clone the repository**
   ```sh
   git clone https://github.com/labbces/NRGSC.git
   cd NRGSC
   ```

2. **Annotation**
   - Run `annotation/transdecoder/run_transcoder.sh` to predict ORFs.
   - Use `annotation/blast/code/split_reference.sh` and `run_blastParts.sh` for protein/transcript annotation.
   - Organize results in `annotation/blast/results/`.

3. **Differential Expression**
   - Prepare Salmon quantification files and sample metadata.
   - Use `diff_expression/code/DESEQ_DEA_GENE.r` and/or `DESEQ_DEA_TRANSCRIPTS.r` for DE analysis.
   - Use `make_full_matrixes_for_trinotate_web.r`, `pheatmap_by_groups.r`, and `up_set_plot.r` for downstream analysis and visualization.

4. **Co-Expression Network Analysis**
   - Build correlation matrices and networks using `networks/code/make_graph_diff_corrected.r`, then run clustering (`run_mcl.sh`, `run_himcl2.sh`).
   - Use `corr_modules.r`, `get_hub.r`, `GO_DEG.r`, `topGO.r`, and visualization scripts for module-trait correlation, hub gene detection, functional enrichment, and heatmap generation.

---

## Requirements

- R (with Bioconductor packages: `DESeq2`, `tximport`, `topGO`, `WGCNA`, `tidyverse`, etc.)
- Perl, Bash, BLAST+, TransDecoder, Trinity (for transcriptome assembly, not included)
- Linux environment, access to compute cluster for parallel BLAST and network analysis

---

## Data Availability

- **Raw RNA-seq**: [NCBI BioProject PRJNA1176579](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1176579)
- **Processed data, assemblies, networks, functional enrichment**: [Figshare project](https://figshare.com/projects/NitrogeGene_co-expression_networks_highlight_key_nodes_associated_with_ammonium_nitrate_response_in_sugarcane/242225)
- **Analysis code and workflows**: This repository

---

## FAIR Compliance

- All data and code are publicly available.
- Processed matrices, co-expression networks, and functional enrichment tables are provided for reuse.

---

## Citation

If you use this resource, please cite:

Muñoz-Pérez JM, Bassi D, Mattiello L, Menossi M, Riaño-Pachón DM.  
Gene co-expression networks highlight key nodes associated with ammonium nitrate response in sugarcane.  
**[BioRxiv, 2025-05]**  
NCBI BioProject: PRJNA1176579  
Figshare: 10.6084/m9.figshare.30100633

---

## Acknowledgements

Funded by CNPq, FAPESP, and CAPES. See manuscript for details.

---

**Contact:** diego.riano@cena.usp.br  
**Lab:** https://github.com/labbces
