# SRS ShortRead Map & Variant Calling (Multisample) - AWS HealthOmics Workflow

This AWS HealthOmics private workflow (written in WDL) performs read mapping and variant calling on paired tumor-normal short read FASTQ files. It supports multiple samples in a batch run and uses industry-standard variant callers (Sentieon's TNScope and Google's DeepSomatic).

---

## Overview

- **Workflow ID:** `1651565`
- **Workflow version** v6
- Aligns tumor and normal FASTQ files using Sentieon DNASeq
- Calls variants using TNScope and DeepSomatic
- Supports multisample input with shared or unique normal samples

---

## Repository Contents

- `SRS-ShortRead-Map-VC-Multisample.wdl` – Workflow definition
- `SRS_Shortread_Map_VC_Multisample_parameters_definition.json` – Input schema and descriptions
- `SRS_Shortread_Map_VC_Multisample_parameters.json` – Sample input JSON file

---

## Required Inputs

### Parameters

cohort: An array that defines tumor-normal FASTQ file groups along with sample metadata. Required.
reference_name: The genome reference to use (options include hg38, t2t_mat, and t2t_pat). Required.
canonical_user_id: Your AWS canonical user ID, needed for obtaining a Sentieon license. Required.
sentieon_docker: URI for the Sentieon Docker image to use in the workflow. Required.
deepsomatic_docker: URI for the DeepSomatic Docker image. Required.
n_threads: Number of vCPUs to allocate (default is 32). Optional.
memory: Amount of memory to allocate for tasks (default is 64 GiB). Optional.

---
## Parameter Json Creation
The input parameter file can be created manually or using the shortread_parameter_creation.sh script. 

Create a sample sheet csv copied from the WGS_Metadata_Table. Keep all headers from the WGS_Metadata_Table and copy only the rows for the corresponding samples you want to run. The csv MUST contain at least the following headers:

  - Sample ID (B1 or B2 refers to batch)
  - Sample Type
  - Read 1 Fastq S3 Path
  - Read 2 Fastq S3 Path 

Run the command as shown below

```sh
./shortread_parameter_creation.sh -i WGS_Metadata_Table.csv -o test.json
```

## Example Input Parameter Structure

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


## Example Command to run the HealthOmics Pipeline

```sh
aws omics start-run --role-arn "arn:aws:iam::860660336427:role/service-role/OmicsWorkflow-20240124114364" --workflow-id 1651565 --name "<Input Descriptive Name Here> $(date +%Y%m%d-%H%M%S)" --output-uri s3://srs-hg002-seq/main-broad-082025/analysis/short-read/ --storage-capacity 5000 --parameters file://SRS_Shortread_Map_VC_Multisample_parameters.json --workflow-version-name 'v6'
```
