# SRS Short-Read Tumor-Normal Pipeline

**Somatic Alignment and Variant Calling Workflow for the SRS Initiative**

**Research Use Only — Not a Production or Clinical Pipeline**

This repository contains the short-read tumor-normal workflow used in
the **Somatic Reference Samples (SRS) Initiative** to generate
alignments and somatic variant callsets for engineered HG002 clones and
FFPE mixtures described in:

Daniels et al., *Somatic Reference Sample Development and Evaluation
Using Unedited and CRISPR-Cas9 Edited Human Cell Lines* (in preparation)

This workflow was used to generate orthogonal short-read somatic
callsets reported in the manuscript (Methods; Data & Code Availability).

This repository reflects the research workflow used during the SRS
initiative and is provided for transparency and reproducibility.
It is **not intended for clinical use, regulated environments, or production deployment.**

------------------------------------------------------------------------

## Overview

This AWS HealthOmics WDL workflow performs:

- Tumor–normal alignment (BWA-MEM, GRCh38)
- Somatic SNV and INDEL calling using:
  - **Sentieon TNscope (202112.07)**
  - **DeepSomatic (v1.8, Google)**
- Multisample batch execution
- Support for clone tumor-normal runs and FFPE tumor-only analysis (modified parameters)

The workflow was applied to:

- Illumina WGS (30×, 120×)
- Illumina WGS FFPE (600× aggregated)
- Illumina WES (600×)
- RNA-seq (handled separately via DRAGEN RNA pipeline)

------------------------------------------------------------------------

## Repository Contents

- `SRS-ShortRead-Map-VC-Multisample.wdl`  
  Primary workflow definition (WDL)

- `SRS_Shortread_Map_VC_Multisample_parameters_definition.json`  
  Input schema and parameter definitions

- `SRS_Shortread_Map_VC_Multisample_parameters.json`  
  Example input file

- `shortread_parameter_creation.sh`  
  Helper script to generate workflow parameter JSON from sample sheet

- `USAGE.md`  
  Execution guidance for AWS HealthOmics deployment

------------------------------------------------------------------------

## Workflow Architecture

### Alignment

- **Aligner:** BWA-MEM  
- **Reference:** GRCh38 (DRAGEN hg38 multigenome v4)  
- **Output:** Coordinate-sorted BAM/CRAM  

### Somatic Variant Calling

Two independent callers are executed:

1. **Sentieon TNscope (v202112.07)** — Sensitive somatic SNV/INDEL detection  
2. **DeepSomatic (v1.8)** — Deep learning–based somatic variant caller  

Variants are retained per caller and downstream merged/unified as
described in the manuscript (Repun-based unification).

------------------------------------------------------------------------

## Computational Environment

- Platform: AWS HealthOmics  
- Execution backend: WDL  
- Compute: dynamically allocated (f2 instances typical for SRS runs)  
- Storage: S3 streaming  

------------------------------------------------------------------------

## Reproducibility Notes

- Reference: GRCh38 (DRAGEN hg38 multigenome v4)  
- Chromosomes analyzed: chr1–22, X, Y  
- PASS variants extracted for downstream analysis  
- VAFs derived directly from FORMAT fields (allele support / total depth)  
- Structural variants handled separately in DRAGEN pipeline  

Variant representation unification and benchmark construction are
described in the Repun pipeline:  
https://github.com/nate-d-olson/mdic_repun

------------------------------------------------------------------------

## Relationship to the SRS Manuscript

This workflow generated:

- Orthogonal short-read somatic callsets (30×, 120×)
- FFPE validation datasets (600×)
- Data used for:
  - Allele frequency stability analysis
  - Clone-specific variant detection
  - Benchmark tier construction
  - Cross-caller concordance analysis

------------------------------------------------------------------------

## Citation

If you use this workflow, please cite:

Daniels CA, Abdulkadir A, Cleveland MH, et al.  
*Somatic Reference Sample Development and Evaluation Using Unedited and CRISPR-Cas9 Edited Human Cell Lines*. (in preparation)

------------------------------------------------------------------------

## License

(Add license here)
