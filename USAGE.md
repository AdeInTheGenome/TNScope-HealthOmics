# Usage Guide

**Research Workflow — Not Production Hardened**

This workflow was developed for the Somatic Reference Samples (SRS)
research initiative and is not intended for clinical diagnostics,
regulated pipelines, or production environments.

Users are responsible for validating performance and compliance in
their own environment.

------------------------------------------------------------------------

## Platform Requirements

- AWS HealthOmics
- IAM roles with:
  - S3 read/write permissions
  - HealthOmics workflow execution permissions
- Sentieon license (canonical_user_id required)

------------------------------------------------------------------------

## Required Inputs

### Mandatory Parameters

- `cohort` — Tumor-normal FASTQ group definitions
- `reference_name` — Genome reference (`hg38`, `t2t_mat`, `t2t_pat`)
- `canonical_user_id` — Required for Sentieon license
- `sentieon_docker` — Sentieon container URI
- `deepsomatic_docker` — DeepSomatic container URI

### Optional Parameters

- `n_threads` (default: 32)
- `memory` (default: 64 GiB)

------------------------------------------------------------------------

## Example Input JSON

```json
{
  "cohort": {
    "samples": [
      {
        "tumor_sample_name": "BRAF-V600E_clone",
        "normal_sample_name": "HG002_parent",
        "r1_tumor_fastqs": ["s3://bucket/tumor_R1.fastq.gz"],
        "r2_tumor_fastqs": ["s3://bucket/tumor_R2.fastq.gz"],
        "r1_normal_fastqs": ["s3://bucket/normal_R1.fastq.gz"],
        "r2_normal_fastqs": ["s3://bucket/normal_R2.fastq.gz"]
      }
    ]
  },
  "reference_name": "hg38",
  "canonical_user_id": "YOUR_ID",
  "sentieon_docker": "URI",
  "deepsomatic_docker": "URI"
}

### Running the Workflow (AWS HealthOmics)

aws omics start-run \
  --workflow-id <workflow_id> \
  --name "SRS_shortread_$(date +%Y%m%d)" \
  --output-uri s3://<output-bucket>/ \
  --parameters file://parameters.json \
  --workflow-version-name "v6"
