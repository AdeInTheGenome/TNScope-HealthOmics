version 1.0

workflow run_deepsomatic {
    input {
        File aligned_reads
        File normal_reads
        File aligned_index
        File normal_index
        File ref_fasta
        File ref_fasta_index
        Array[File] contigs
        Int threads = 16
        Int preemptible_tries
        String tumor_name
        String normal_name
        String deepsomatic_docker
    }

    scatter (ctg in contigs){
        call call_DeepSomatic {
            input:
                aligned_reads = aligned_reads,
                normal_reads = normal_reads,
                aligned_index = aligned_index,
                normal_index = normal_index,
                ref_fasta = ref_fasta,
                ref_fai = ref_fasta_index,
                contig = ctg,
                tumor_name = tumor_name,
                normal_name = normal_name,
                threads = threads,
                deepsomatic_docker = deepsomatic_docker
        }
    }


    call gather_deepsomatic {
        input:
            deepsomatic_vcfs = call_DeepSomatic.DeepSomatic_calls_vcf,
            deepsomatic_gvcfs = call_DeepSomatic.DeepSomatic_calls_gvcf,
            threads = threads,
            tumor_name = tumor_name,
            normal_name = normal_name
    }


    output {
        File DeepSomatic_calls_vcf = gather_deepsomatic.DeepSomatic_calls_vcf
        File DeepSomatic_calls_vcf_tbi = gather_deepsomatic.DeepSomatic_calls_vcf_tbi
        File DeepSomatic_calls_gvcf = gather_deepsomatic.DeepSomatic_calls_gvcf
        File DeepSomatic_calls_gvcf_tbi = gather_deepsomatic.DeepSomatic_calls_gvcf_tbi
        }
}

task call_DeepSomatic {
    input {
        # Input tumor bam files
        File aligned_reads
        File aligned_index
        String tumor_name

        # Input normmal bam files
        File normal_reads
        File normal_index
        String normal_name

        # Reference genome files
        File ref_fasta
        File ref_fai

        File? contig 

        # Execution
        String threads = "16"
        Int preemptible_tries = 3
        String deepsomatic_docker
    }

    Float file_size = ceil(size(aligned_reads, "GB") * 2 + size(normal_reads, "GB") * 2 + size(ref_fasta, "GB") + size(contig, "GB") + 20)

    command <<<
    set -xv
    set -exvuo pipefail

    run_deepsomatic \
    --model_type=WGS \
    --ref="~{ref_fasta}" \
    --reads_normal="~{normal_reads}" \
    --reads_tumor="~{aligned_reads}" \
    --output_vcf=sample_deepsomatic.vcf.gz \
    --output_gvcf=sample_deepsomatic.g.vcf.gz \
    --sample_name_tumor="~{tumor_name}" \
    --sample_name_normal="~{normal_name}" \
    --num_shards="~{threads}" \
    --logging_dir=logs \
    ~{"--regions=" + contig}


    # Filter for PASS
    bcftools --version
    bcftools view \
        -f PASS -Oz \
        -o "sample_deepsomatic.PASS.vcf.gz"\
        "sample_deepsomatic.vcf.gz"
    tabix -p vcf "sample_deepsomatic.PASS.vcf.gz"
    >>>

    output {
    File DeepSomatic_calls_vcf = "sample_deepsomatic.vcf.gz"
    File DeepSomatic_calls_vcf_tbi = "sample_deepsomatic.vcf.gz.tbi"
    File DeepSomatic_calls_gvcf = "sample_deepsomatic.g.vcf.gz"
    File DeepSomatic_calls_gvcf_tbi = "sample_deepsomatic.g.vcf.gz.tbi"
  }

    runtime {
        docker: deepsomatic_docker
        cpu: threads
        memory: "~{threads * 8} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
        acceleratorCount: 1
        acceleratorType: "nvidia-tesla-t4"
    }
}

task gather_deepsomatic {
    input {
        Array[File] deepsomatic_vcfs
        Array[File] deepsomatic_gvcfs
        Int threads
        String tumor_name
        String normal_name
    }

    Float file_size = ceil((size(deepsomatic_vcfs, "GB") * 2) + (size(deepsomatic_gvcfs, "GB") * 2) + 20)

    command <<<
    set -euxo pipefail

    echo "Merging deepsomatic vcfs for ~{tumor_name}"

    bcftools --version

    echo '~{sep="\n" deepsomatic_vcfs}' > vcf.list
    echo '~{sep="\n" deepsomatic_gvcfs}' > gvcf.list

    # Sometimes the indices failed to be localized properly. 
    # Reindex here to make it more robust...
    count=0
    for file in $(cat vcf.list)
    do
        cp ${file} ${count}.vcf.gz
        tabix ${count}.vcf.gz
        count=$((count+1))
    done
    countmax_vcf=${count}

    count=0
    for file in $(cat gvcf.list)
    do
        cp ${file} ${count}.gvcf.gz
        tabix ${count}.gvcf.gz
        count=$((count+1))
    done
    countmax_gvcf=${count}

    # Create new VCF list from the copied files
    ls *.vcf.gz > vcf.list
    ls *.gvcf.gz > gvcf.list

    bcftools concat -a -f vcf.list |\
        # Set GT to 0/1 for somatic variants as DeepSomatic set 1/1 which will not be phased by HiPhase
        bcftools +setGT -- -t q -n c:"0/1" -i 'FMT/GT="1/1"' |\
        bcftools sort -Oz -o sample_deepsomatic.vcf.gz
    bcftools concat -a -f gvcf.list |\
        bcftools sort -Oz -o sample_deepsomatic.g.vcf.gz

    tabix -p vcf sample_deepsomatic.vcf.gz
    tabix -p vcf sample_deepsomatic.g.vcf.gz

    # Remove intermediate files using countmax
    for i in $(seq 0 $((countmax_vcf-1)))
    do
        rm ${i}.vcf.gz ${i}.vcf.gz.tbi
    done

    for i in $(seq 0 $((countmax_gvcf-1)))
    do
        rm ${i}.gvcf.gz ${i}.gvcf.gz.tbi
    done
    >>>

    output {
    File DeepSomatic_calls_vcf = "sample_deepsomatic.vcf.gz"
    File DeepSomatic_calls_vcf_tbi = "sample_deepsomatic.vcf.gz.tbi"
    File DeepSomatic_calls_gvcf = "sample_deepsomatic.g.vcf.gz"
    File DeepSomatic_calls_gvcf_tbi = "sample_deepsomatic.g.vcf.gz.tbi"
    }

    runtime {
        docker: "860660336427.dkr.ecr.us-east-1.amazonaws.com/hifisomatic:bcftools"
        cpu: threads
        memory: "~{threads * 4} GB"
        disk: file_size + " GB"
        maxRetries: 2
        preemptible: 1
        }
}