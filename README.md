# SRS ShortRead Map & Variant Calling (Multisample) - AWS HealthOmics Workflow

This AWS HealthOmics private workflow performs read mapping and variant calling on paired tumor-normal short read FASTQ files. It supports multiple samples in a batch run and uses industry-standard variant callers (Sentieon's TNScope and Google's DeepSomatic).

---

## Overview

- **Workflow ID:** `1651565`
- **Latest Successful Run ID:** `2429213`
- **Output Location:** [Run 2429213 Outputs in S3](https://s3.console.aws.amazon.com/s3/buckets/srs-poc-test/omics-test-ade/shortread_pipeline_multisample_test/2429213/?region=us-east-1&tab=objects)

This pipeline:
- Aligns tumor and normal FASTQ files using Sentieon
- Calls variants using multiple callers including DeepSomatic
- Supports multisample input with shared or unique normal samples

---

## Repository Contents

- `SRS-ShortRead-Map-VC-Multisample.wdl` – Workflow definition
- `SRS_Shortread_Map_VC_Multisample_parameters_definition.json` – Input schema and descriptions
- `SRS_Shortread_Map_VC_Multisample_parameters.json` – Sample input JSON file

---

## Required Inputs

### Parameters

| Name                  | Type    | Description                                                              | Required |
|-----------------------|---------|--------------------------------------------------------------------------|----------|
| `cohort`              | Array   | Tumor-normal FASTQ file groups with sample metadata                      | Yes      |
| `reference_name`      | String  | Genome reference (`hg38`, `t2t_mat`, `t2t_pat`)                           | Yes      |
| `canonical_user_id`   | String  | AWS canonical user ID (for Sentieon license)                             | Yes      |
| `sentieon_docker`     | String  | Sentieon Docker image URI                                                | Yes      |
| `deepsomatic_docker`  | String  | DeepSomatic Docker image URI                                             | Yes      |
| `n_threads`           | Int     | Number of vCPUs (default: 32)                                            | No       |
| `memory`              | String  | Memory per task (default: 64 GiB)                                        | No       |

---

## Sample Input Structure

```json
{
  "cohort": {
    "samples": [
      {
        "tumor_sample_name": "test_sample1",
        "normal_sample_name": "chr4_chr7_normal",
        "r1_tumor_fastqs": [
          "s3://srs-poc-test/fastqs/compressed/chr4_chr7_tumor_R1.fastq.gz"
        ],
        "r2_tumor_fastqs": [
          "s3://srs-poc-test/fastqs/compressed/chr4_chr7_tumor_R2.fastq.gz"
        ],
        "r1_normal_fastqs": [
          "s3://srs-poc-test/fastqs/compressed/chr4_chr7_normal_R1.fastq.gz"
        ],
        "r2_normal_fastqs": [
          "s3://srs-poc-test/fastqs/compressed/chr4_chr7_normal_R2.fastq.gz"
        ],
        "tumor_read_groups": [
          "@RG\\tID:chr4_chr7_tumor_test1\\tSM:chr4_chr7_tumor_test1\\tPL:ILLUMINA"
        ],
        "normal_read_groups": [
          "@RG\\tID:normal_sample\\tSM:chr4_chr7_normal\\tPL:ILLUMINA"
        ]
      }
    ]
  },
  "reference_name": "hg38",
  "canonical_user_id": "YOUR_CANONICAL_ID",
  "sentieon_docker": "YOUR_SENTIEON_IMAGE_URI",
  "deepsomatic_docker": "YOUR_DEEPSOMATIC_IMAGE_URI"
}
```
