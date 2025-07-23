version 1.1

# Somatic variant calling with Sentieon TNScope and Google's Deepsomatic

struct ShortReadSample {
    String tumor_sample_name
    String normal_sample_name
    Array[File] r1_tumor_fastqs
    Array[File] r2_tumor_fastqs
    Array[File] r1_normal_fastqs
    Array[File] r2_normal_fastqs
    Array[String] tumor_read_groups
    Array[String] normal_read_groups
}

struct Cohort {
    Array[ShortReadSample] samples
}

workflow short_read_variant_calling {
  input {
    # Individual tumor sample inputs as defined in the JSON.
    Cohort cohort

    # Reference and other parameters
    String reference_name

    String canonical_user_id
    String sentieon_license = "aws-omics.sentieon.com:9011"

    # Execution
    String n_threads_deepsomatic = "16"
    String n_threads_sentieon = "32"
    String memory_deepsomatic = "~{n_threads_deepsomatic * 4} GB"
    String memory_sentieon = "~{n_threads_sentieon * 2} GB"
    Int preemptible_tries = 3
    String deepsomatic_docker
    String sentieon_docker
    # Default number of threads for misc tasks (4GB per thread assigned for almost all steps)
    Int def_threads = 2
    # Scatter small variants calling into equal chunk per chromosome to make use of multiple nodes. Default of 75 Mbp per chromosome (total of 42 chunks for hg38)
    Int chunk_size = 75000000

    # process arguments
    String bwa_xargs = ""
    String bwa_karg = "10000000"
    String sort_xargs = "--bam_compression 1"
    String lc_xargs = ""
    String alnstat_xargs = "--adapter_seq \'\'" # Corrected quoting
    String dedup_xargs = "--cram_write_options version=3.0,compressor=rans"
    String calling_driver_xargs = ""
    String calling_algo_xargs = ""

    # Scatter small variants calling into equal chunk per chromosome to make use of multiple nodes. Default of 75 Mbp per chromosome (total of 42 chunks for hg38)
    Int chunk_size = 75000000
  }

  # Perform a license check
  call SentieonLicense {
    input:
      canonical_user_id = canonical_user_id,
      sentieon_license = sentieon_license,
      sentieon_docker = sentieon_docker
  }

  # Download the reference genome
  call DownloadReference {
    input:
      reference_name = reference_name,
      sentieon_docker = sentieon_docker
  }

  # Process each tumor sample against the single normal sample.
  scatter(individual in cohort.samples) {

    call SentieonMapping {
      input:
        r1_fastq = individual.r1_tumor_fastqs,
        r2_fastq = individual.r2_tumor_fastqs,
        tumor_name = individual.tumor_sample_name,
        read_groups = individual.tumor_read_groups,

        normal_r1_fastq = individual.r1_normal_fastqs,
        normal_r2_fastq = individual.r2_normal_fastqs,
        normal_name = individual.normal_sample_name,
        normal_read_groups = individual.normal_read_groups,

        ref_fasta = DownloadReference.ref_fasta,
        ref_fai = DownloadReference.ref_fai,
        ref_bwt = DownloadReference.ref_bwt,
        ref_sa = DownloadReference.ref_sa,
        ref_amb = DownloadReference.ref_amb,
        ref_ann = DownloadReference.ref_ann,
        ref_pac = DownloadReference.ref_pac,

        bwa_xargs = bwa_xargs,
        bwa_karg = bwa_karg,
        sort_xargs = sort_xargs,
        lc_xargs = lc_xargs,
        alnstat_xargs = alnstat_xargs,
        dedup_xargs = dedup_xargs,
        calling_driver_xargs = calling_algo_xargs,

        license_ok = SentieonLicense.license_ok,
        canonical_user_id = canonical_user_id,
        sentieon_license = sentieon_license,

        n_threads = n_threads_sentieon,
        memory = memory_sentieon,
        preemptible_tries = preemptible_tries,
        sentieon_docker = sentieon_docker
    }


    call SentieonVariantCalling {
      input:
        aligned_reads = SentieonMapping.aligned_reads,
        aligned_index = SentieonMapping.aligned_index,
        tumor_name = individual.tumor_sample_name,

        normal_reads = SentieonMapping.normal_reads,
        normal_index = SentieonMapping.normal_index,
        normal_name = individual.normal_sample_name,

        ref_fasta = DownloadReference.ref_fasta,
        license_ok = SentieonLicense.license_ok,
        canonical_user_id = canonical_user_id,
        sentieon_license = sentieon_license,

        n_threads = n_threads_sentieon,
        memory = memory_sentieon,
        preemptible_tries = preemptible_tries,
        sentieon_docker = sentieon_docker,
        calling_driver_xargs = calling_driver_xargs,
        calling_algo_xargs = calling_algo_xargs
      }

    call split_contigs {
      input:
        ref_fasta_index = DownloadReference.ref_fai,
        chunk_size = chunk_size,
        threads = def_threads
    }

    call DeepSomatic {
      input:
        aligned_reads = SentieonMapping.aligned_reads,
        aligned_index = SentieonMapping.aligned_index,
        tumor_name = individual.tumor_sample_name,
        normal_reads = SentieonMapping.normal_reads,
        normal_index = SentieonMapping.normal_index,
        normal_name = individual.normal_sample_name,
        ref_fasta = DownloadReference.ref_fasta,
        ref_fai = DownloadReference.ref_fai,
        n_threads = n_threads_deepsomatic,
        memory = memory_deepsomatic,
        preemptible_tries = preemptible_tries,
        deepsomatic_docker = deepsomatic_docker,
        contigs = split_contigs.contigs
    }
  }

  output {
    Array[File] aligned_reads = SentieonMapping.aligned_reads
    Array[File] aligned_index = SentieonMapNormal.aligned_index
    Array[File?]? normal_reads = SentieonMapping.normal_reads
    Array[File?]? normal_index = SentieonMapping.normal_index

    Array[File] Sentieon_calls_vcf = SentieonVariantCalling.Sentieon_calls_vcf
    Array[File] Sentieon_calls_vcf_tbi = SentieonVariantCalling.Sentieon_calls_vcf_tbi

    Array[File] Sentieon_dedup_metrics = SentieonMapping.Sentieon_dedup_metrics
    Array[File] Sentieon_mq_metrics = SentieonMapping.Sentieon_mq_metrics
    Array[File] Sentieon_qd_metrics = SentieonMapping.Sentieon_qd_metrics
    Array[File] Sentieon_gc_summary = SentieonMapping.Sentieon_gc_summary
    Array[File] Sentieon_gc_metrics = SentieonMapping.Sentieon_gc_metrics
    Array[File] Sentieon_as_metrics = SentieonMapping.Sentieon_as_metrics
    Array[File] Sentieon_is_metrics = SentieonMapping.Sentieon_is_metrics

    Array[File] Sentieon_mq_plot = SentieonMapping.Sentieon_mq_plot
    Array[File] Sentieon_qd_plot = SentieonMapping.Sentieon_qd_plot
    Array[File] Sentieon_gc_plot = SentieonMapping.Sentieon_gc_plot
    Array[File] Sentieon_is_plot = SentieonMapping.Sentieon_is_plot

    Array[File] DeepSomatic_calls_vcf = DeepSomatic.DeepSomatic_calls_vcf
    Array[File] DeepSomatic_calls_vcf_tbi = DeepSomatic.DeepSomatic_calls_vcf_tbi
    Array[File] DeepSomatic_calls_gvcf = DeepSomatic.DeepSomatic_calls_gvcf
    Array[File] DeepSomatic_calls_gvcf_tbi = DeepSomatic.DeepSomatic_calls_gvcf_tbi

    Array[File] DeepSomatic_visual_report = DeepSomatic.DeepSomatic_visual_report
  }
}

task SentieonLicense {
  input {
    String canonical_user_id
    String sentieon_license
    String sentieon_docker
  }
  command <<<
    set -exvuo pipefail
    source /opt/sentieon/omics_credentials.sh "~{sentieon_license}" "~{canonical_user_id}"
    sentieon licclnt ping && echo "Ping is OK"
    echo "License OK" > license_ok.txt
  >>>
  runtime {
    docker: sentieon_docker
    memory: "1 GB"
    cpu: 1
    preemptible: 3
  }
  output {
    File license_ok = "license_ok.txt"
  }
}

task DownloadReference {
  input {
    String reference_name
    String sentieon_docker
  }
  command <<<
    set -exvuo pipefail

    ref_build=~{reference_name}
    s3_bucket_basename="s3://srs-poc-test"
    ref_idx_files=(
      ""
      ".fai"
      ".amb"
      ".ann"
      ".bwt"
      ".pac"
      ".sa"
    )

    case "$ref_build" in
      t2t_mat)
        ref_base="ref/t2t-q100/v1.0.1-maternal/uncompressed/hg002v1.0.1.mat_Y_EBV_MT.fasta"
        has_alt=false
        has_dbSNP=false
        has_knownsites=false
        ;;
      t2t_pat)
        ref_base="ref/t2t-q100/v1.0.1-paternal/uncompressed/hg002v1.0.1.pat_X_EBV_MT.fasta"
        has_alt=false
        has_dbSNP=false
        has_knownsites=false
        ;;
      hg38_alt|hg38_gatk)
        ref_base="hg38/Homo_sapiens_assembly38.fasta"
        has_alt=true
        has_dbSNP=true
        has_knownsites=true
        VCFS=(
          "hg38/dbsnp_138.hg38.vcf.gz"
          "hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
          "hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        )
        ;;
      hg38|hg38_noalt|hs38|GRCh38)
        ref_base="ref/hg38_giab_v3/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta"
        has_alt=false
        has_dbSNP=true
        has_knownsites=true
        ;;
      b37_gatk)
        ref_base="b37/Homo_sapiens_assembly19.fasta"
        has_alt=false
        has_dbSNP=true
        has_knownsites=true
        VCFS=(
          "b37/dbsnp_138.b37.vcf.gz"
          "b37/1000G_phase1.indels.b37.vcf.gz"
          "b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
        )
        ;;
      b37|hs37d5)
        ref_base="b37/hs37d5.fa"
        has_alt=false
        has_dbSNP=true
        has_knownsites=true
        VCFS=(
          "b37/dbsnp_138.b37.vcf.gz"
          "b37/1000G_phase1.indels.b37.vcf.gz"
          "b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
        )
        ;;
      hg19|ucsc_hg19)
        ref_base="hg19/ucsc.hg19.fasta"
        has_alt=false
        has_dbSNP=true
        has_knownsites=true
        VCFS=(
          "hg19/dbsnp_138.hg19.vcf.gz"
          "hg19/1000G_phase1.indels.hg19.sites.vcf.gz"
          "hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
        )
        ;;
      quickstart)
        ref_base="quickstart/ucsc.hg19_chr22.fasta"
        has_alt=false
        has_dbSNP=true
        has_knownsites=true
        VCFS=(
          "quickstart/dbsnp_135.hg19_chr22.vcf.gz"
          "quickstart/1000G_phase1.snps.high_confidence.hg19_chr22.sites.vcf.gz"
          "quickstart/Mills_and_1000G_gold_standard.indels.hg19_chr22.sites.vcf.gz"
        )
        ;;
      *)
        echo "ERROR: unknown genome build"
        exit 1
        ;;
    esac

    for idx in "${ref_idx_files[@]}"; do
      aws s3 cp "$s3_bucket_basename/$ref_base$idx" "./reference.fa$idx"
    done
    if [[ "$has_alt" == true ]]; then
      aws s3 cp "$s3_bucket_basename/${ref_base}.alt" "./reference.fa.alt"
    fi
    if [[ "$has_dbSNP" == true ]]; then
      aws s3 cp "$s3_bucket_basename/${ref_base}.alt" "./reference.fa.alt"
      for i in $(seq 1 "${#VCFS[@]}"); do
      i=$((i - 1))
      vcf="${VCFS[$i]}"
      aws s3 cp "$s3_bucket_basename/$vcf" "./sites_${i}.vcf.gz"
      aws s3 cp "$s3_bucket_basename/${vcf}.tbi" "./sites_${i}.vcf.gz.tbi"
    done
    fi
    if [[ "$has_knownsites" == true ]]; then
      aws s3 cp "$s3_bucket_basename/${ref_base}.alt" "./reference.fa.alt"
      for i in $(seq 1 "${#VCFS[@]}"); do
      i=$((i - 1))
      vcf="${VCFS[$i]}"
      aws s3 cp "$s3_bucket_basename/$vcf" "./sites_${i}.vcf.gz"
      aws s3 cp "$s3_bucket_basename/${vcf}.tbi" "./sites_${i}.vcf.gz.tbi"
    done
    fi
  >>>
  runtime {
    preemptible: 3
    docker: sentieon_docker
    memory: "4 GB"
    cpu: 1
  }
  output {
    File ref_fasta = "reference.fa"
    File ref_fai   = "reference.fa.fai"
    File ref_bwt   = "reference.fa.bwt"
    File ref_sa    = "reference.fa.sa"
    File ref_amb   = "reference.fa.amb"
    File ref_ann   = "reference.fa.ann"
    File ref_pac   = "reference.fa.pac"
  }
}

# Use bedtools to split contigs
task split_contigs {
  input {
    File ref_fasta_index
    Int chunk_size = 75000000
    Int threads
  }

  Float file_size = ceil(size(ref_fasta_index, "GB") + 10)

  command <<<
  set -euxo pipefail

  bedtools --version

  echo "Splitting contigs for ~{ref_fasta_index}"
  bedtools makewindows -g ~{ref_fasta_index} -w ~{chunk_size} > contigs.bed
  grep -v -E "random|chrUn|chrM|chrEBV" contigs.bed > noalt.bed
  # Split the contig bed files into one file for each line
  split -l 1 noalt.bed contigs_split.
  # Add .bed to all the contigs_split file
  for file in $(ls contigs_split.*); do mv $file $file.bed; done
  >>>

  output {
    Array[File] contigs = glob("contigs_split.*.bed")
  }

  runtime {
    docker: "860660336427.dkr.ecr.us-east-1.amazonaws.com/hifisomatic:bedtools"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

task SentieonMapping {
  input {
    # Input tumor fastq files
    Array[File] r1_fastq = []
    Array[File] r2_fastq = []
    Array[String] read_groups = []
    String tumor_name

    # Input normal fastq files
    Array[File] normal_r1_fastq = []
    Array[File] normal_r2_fastq = []
    Array[String] normal_read_groups = []
    String normal_name
    
    
    # Reference genome files
    File ref_fasta
    File ref_fai
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac

    # Optional process arguments
    String bwa_xargs = ""
    String bwa_karg = "10000000"
    String sort_xargs = "--bam_compression 1"
    String lc_xargs = ""
    String alnstat_xargs = "--adapter_seq ''"
    String dedup_xargs = "--cram_write_options version=3.0,compressor=rans"
    String qcal_xargs = ""
    String calling_driver_xargs = ""
    String calling_algo_xargs = ""

    # Sentieon license configuration
    File license_ok
    String canonical_user_id
    String sentieon_license = "aws-omics.sentieon.com:9011"

    # Execution
    String n_threads = "32"
    String memory = "64 GiB"
    Int preemptible_tries = 3
    String sentieon_docker

  }
  command <<<
    set -xv
    source /opt/sentieon/omics_credentials.sh "~{sentieon_license}" "~{canonical_user_id}"
    set -exvuo pipefail

    numa_nodes=$(lscpu | grep "NUMA node(s):" | sed 's/^NUMA node.* //')
    numa_cpulist=()
    for i in $(seq 1 "$numa_nodes"); do
        i=$((i - 1))
        numa_cpulist+=($(lscpu | grep "NUMA node$i CPU" | sed 's/^NUMA.* //'))
    done

    nt=$(nproc)
    n_threads=$((nt / numa_nodes))

    mem_kb=$(cat /proc/meminfo | grep "MemTotal" | awk '{print $2}')
    bwt_mem=$((mem_kb / 1024 / 1024 / numa_nodes - 10))
    export bwt_max_mem="$bwt_mem"G

    r1_fastq=("~{sep='" "' r1_fastq}")
    r2_fastq=("~{sep='" "' r2_fastq}")
    read_groups=("~{sep='" "' read_groups}")

    normal_r1_fastq=("~{sep='" "' normal_r1_fastq}")
    normal_r2_fastq=("~{sep='" "' normal_r2_fastq}")
    normal_read_groups=("~{sep='" "' normal_read_groups}")

    # Sanity check the input
    if [[ ${#r1_fastq[@]} -ne ${#read_groups[@]} ]]; then
      echo "The number of readgroups does not match the number of input fastq files"
      exit 1
    fi
    if [[ ${#r1_fastq[@]} -ne ${#r2_fastq[@]} ]]; then
      echo "The number of r1 fastq does not equal the number of r2 fastq"
      exit 1
    fi
    if [[ ${#normal_r1_fastq[@]} -ne ${#normal_read_groups[@]} ]]; then
      echo "The number of normal readgroups does not match the number of input fastq files"
      exit 1
    fi
    if [[ ${#normal_r1_fastq[@]} -ne ${#normal_r2_fastq[@]} ]]; then
      echo "The number of normal r1 fastq does not equal the number of r2 fastq"
      exit 1
    fi

    # Alignment with BWA
    alignment_output=()
    for i in $(seq 1 "$numa_nodes"); do
      i=$((i - 1))
      r1="${r1_fastq[$i]}"
      r2="${r2_fastq[$i]}"
      rg="${read_groups[$i]}"

      for j in $(seq 1 "$numa_nodes"); do 
        j=$((j - 1))
        cpulist="${numa_cpulist[$j]}"

        # Alignment command
        perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)';
        taskset -c "$cpulist" sentieon bwa mem -R "$rg" \
          ~{bwa_xargs} -K ~{bwa_karg} -t $n_threads "~{ref_fasta}" \
          <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
            sentieon fqidx extract -F "$j"/"$numa_nodes" -K ~{bwa_karg} \
            <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
              igzip -dc "$r1") \
            <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
              igzip -dc "$r2")) | \
          taskset -c "$cpulist" sentieon util sort -t $n_threads --sam2bam \
          -o "sample_sorted_${i}_${j}.bam" -i - ~{sort_xargs} &
        alignment_output+=("sample_sorted_${i}_${j}.bam")
      done
      wait

      (rm "$r1" "$r2" || (exit 0)) &  # Try removing the input files to save disk space
    done

    # Alignment with BWA - normal
    normal_alignment_output=()
    for i in $(seq 1 ${#normal_r1_fastq[@]}); do
      i=$((i - 1))
      r1="${normal_r1_fastq[$i]}"
      r2="${normal_r2_fastq[$i]}"
      rg="${normal_read_groups[$i]}"

      for j in $(seq 1 "$numa_nodes"); do
        j=$((j - 1))
        cpulist="${numa_cpulist[$j]}"

        # Alignment command
        perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)';
        taskset -c "$cpulist" sentieon bwa mem -R "$rg" \
          ~{bwa_xargs} -K ~{bwa_karg} -t $n_threads "~{ref_fasta}" \
          <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
            sentieon fqidx extract -F "$j"/"$numa_nodes" -K ~{bwa_karg} \
            <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
              igzip -dc "$r1") \
            <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
              igzip -dc "$r2")) | \
          taskset -c "$cpulist" sentieon util sort -t $n_threads --sam2bam \
          -o "normal_sorted_${i}_${j}.bam" -i - ~{sort_xargs} &
        normal_alignment_output+=("normal_sorted_${i}_${j}.bam")
      done
      wait

      (rm "$r1" "$r2" || (exit 0)) &  # Try removing the input files to save disk space
    done
    

    # Add the bam input files to the bam_str
    bam_str=()
    for i in $(seq 1 ${#alignment_output[@]}); do
      i=$((i - 1))
      bam_str+=("-i")
      bam_str+=("${alignment_output[$i]}")
    done
    normal_bam=()
    for i in $(seq 1 ${#normal_alignment_output[@]}); do
      i=$((i - 1))
      normal_bam+=("-i")
      normal_bam+=("${normal_alignment_output[$i]}")
    done

    # Dedup and Metrics Calc - tumor
    tumor_deduped="sample_aligned.cram"
    sentieon driver "${bam_str[@]}" -r ~{ref_fasta} \
      --algo LocusCollector ~{lc_xargs} "Sentieon_sample_score.txt.gz" \
      --algo MeanQualityByCycle "Sentieon_sample_mq_metrics.txt" \
      --algo QualDistribution "Sentieon_sample_qd_metrics.txt" \
      --algo GCBias --summary "Sentieon_sample_gc_summary.txt" "Sentieon_sample_gc_metrics.txt" \
      --algo AlignmentStat ~{alnstat_xargs} "Sentieon_sample_aln_metrics.txt" \
      --algo InsertSizeMetricAlgo "Sentieon_sample_is_metrics.txt"

    sentieon driver "${bam_str[@]}" -r ~{ref_fasta} --algo Dedup \
      ~{dedup_xargs} --score_info "Sentieon_sample_score.txt.gz" \
      --metrics "Sentieon_sample_dedup_metrics.txt" \
      "$tumor_deduped"

    # Plot the metrics output
    sentieon plot GCBias -o "Sentieon_sample_gc-report.pdf" "Sentieon_sample_gc_metrics.txt" &
    sentieon plot QualDistribution -o "Sentieon_sample_qd-report.pdf" "Sentieon_sample_qd_metrics.txt" &
    sentieon plot MeanQualityByCycle -o "Sentieon_sample_mq-report.pdf" "Sentieon_sample_mq_metrics.txt" &
    sentieon plot InsertSizeMetricAlgo -o "Sentieon_sample_is-report.pdf" "Sentieon_sample_is_metrics.txt" &

    (rm "${alignment_output[@]}" || (exit 0)) &   # Remove intermediate files to save space

    # Dedup and Metrics Calc - normal
    normal_deduped="normal_aligned.cram"
    sentieon driver "${normal_bam[@]}" -r ~{ref_fasta} \
      --algo LocusCollector ~{lc_xargs} "Sentieon_normal_score.txt.gz" \
      --algo MeanQualityByCycle "Sentieon_normal_mq_metrics.txt" \
      --algo QualDistribution "Sentieon_normal_qd_metrics.txt" \
      --algo GCBias --summary "Sentieon_normal_gc_summary.txt" "Sentieon_normal_gc_metrics.txt" \
      --algo AlignmentStat ~{alnstat_xargs} "Sentieon_normal_aln_metrics.txt" \
      --algo InsertSizeMetricAlgo "Sentieon_normal_is_metrics.txt"

    sentieon driver "${normal_bam[@]}" -r ~{ref_fasta} --algo Dedup \
      ~{dedup_xargs} --score_info "Sentieon_normal_score.txt.gz" \
      --metrics "Sentieon_normal_dedup_metrics.txt" \
      "$normal_deduped"

    # Plot the metrics output
    sentieon plot GCBias -o "Sentieon_normal_gc-report.pdf" "Sentieon_normal_gc_metrics.txt" &
    sentieon plot QualDistribution -o "Sentieon_normal_qd-report.pdf" "Sentieon_normal_qd_metrics.txt" &
    sentieon plot MeanQualityByCycle -o "Sentieon_normal_mq-report.pdf" "Sentieon_normal_mq_metrics.txt" &
    sentieon plot InsertSizeMetricAlgo -o "Sentieon_normal_is-report.pdf" "Sentieon_normal_is_metrics.txt" &

    (rm "${normal_alignment_output[@]}" || (exit 0)) &   # Remove intermediate files to save space

    # Extract the samples
    tumor_sm=$(samtools samples "$tumor_deduped" | cut -f 1 | head -n 1)
    normal_sm=$(samtools samples "$normal_deduped" | cut -f 1 | head -n 1)

    wait
    unset http_proxy
    exit 0
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: sentieon_docker
    memory: memory
    cpu: n_threads
  }
  output {
    # Core alignment files
    File aligned_reads = "sample_aligned.cram"
    File aligned_index = "sample_aligned.cram.crai"
    File? normal_reads = "normal_aligned.cram"
    File? normal_index = "normal_aligned.cram.crai"

    # QC output metrics
    File Sentieon_dedup_metrics = "Sentieon_sample_dedup_metrics.txt"
    File Sentieon_mq_metrics = "Sentieon_sample_mq_metrics.txt"
    File Sentieon_qd_metrics = "Sentieon_sample_qd_metrics.txt"
    File Sentieon_gc_summary = "Sentieon_sample_gc_summary.txt"
    File Sentieon_gc_metrics = "Sentieon_sample_gc_metrics.txt"
    File Sentieon_as_metrics = "Sentieon_sample_aln_metrics.txt"
    File Sentieon_is_metrics = "Sentieon_sample_is_metrics.txt"

    # QC output plots
    File Sentieon_mq_plot = "Sentieon_sample_mq-report.pdf"
    File Sentieon_qd_plot = "Sentieon_sample_qd-report.pdf"
    File Sentieon_gc_plot = "Sentieon_sample_gc-report.pdf"
    File Sentieon_is_plot = "Sentieon_sample_is-report.pdf"

    # Normal QC output metrics
    File? Sentieon_normal_dedup_metrics = "Sentieon_normal_dedup_metrics.txt"
    File? Sentieon_normal_mq_metrics = "Sentieon_normal_mq_metrics.txt"
    File? Sentieon_normal_qd_metrics = "Sentieon_normal_qd_metrics.txt"
    File? Sentieon_normal_gc_summary = "Sentieon_normal_gc_summary.txt"
    File? Sentieon_normal_gc_metrics = "Sentieon_normal_gc_metrics.txt"
    File? Sentieon_normal_as_metrics = "Sentieon_normal_aln_metrics.txt"
    File? Sentieon_normal_is_metrics = "Sentieon_normal_is_metrics.txt"

    # Normal QC output plots
    File? Sentieon_normal_mq_plot = "Sentieon_normal_mq-report.pdf"
    File? Sentieon_normal_qd_plot = "Sentieon_normal_qd-report.pdf"
    File? Sentieon_normal_gc_plot = "Sentieon_normal_gc-report.pdf"
    File? Sentieon_normal_is_plot = "Sentieon_normal_is-report.pdf"
  }
}

task SentieonVariantCalling {
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
  }
  command <<<
    # Set the BQSR and calling intervals
    first_chrom=$(head -n 1 ~{ref_fai} | cut -f 1)
    case "$first_chrom" in
      chr10_MATERNAL)
        BQSR_INTERVALS="chr1_MATERNAL,chr2_MATERNAL,chr3_MATERNAL,chr4_MATERNAL,chr5_MATERNAL,chr6_MATERNAL,chr7_MATERNAL,chr8_MATERNAL,chr9_MATERNAL,chr10_MATERNAL,chr11_MATERNAL,chr12_MATERNAL,chr13_MATERNAL,chr14_MATERNAL,chr15_MATERNAL,chr16_MATERNAL,chr17_MATERNAL,chr18_MATERNAL,chr19_MATERNAL,chr20_MATERNAL,chr21_MATERNAL,chr22_MATERNAL"
        CALLING_INTERVALS="chrM,chr1_MATERNAL,chr2_MATERNAL,chr3_MATERNAL,chr4_MATERNAL,chr5_MATERNAL,chr6_MATERNAL,chr7_MATERNAL,chr8_MATERNAL,chr9_MATERNAL,chr10_MATERNAL,chr11_MATERNAL,chr12_MATERNAL,chr13_MATERNAL,chr14_MATERNAL,chr15_MATERNAL,chr16_MATERNAL,chr17_MATERNAL,chr18_MATERNAL,chr19_MATERNAL,chr20_MATERNAL,chr21_MATERNAL,chr22_MATERNAL,chrEBV,chrX_MATERNAL"
        ;;
      chr10_PATERNAL)
        BQSR_INTERVALS="chr1_PATERNAL,chr2_PATERNAL,chr3_PATERNAL,chr4_PATERNAL,chr5_PATERNAL,chr6_PATERNAL,chr7_PATERNAL,chr8_PATERNAL,chr9_PATERNAL,chr10_PATERNAL,chr11_PATERNAL,chr12_PATERNAL,chr13_PATERNAL,chr14_PATERNAL,chr15_PATERNAL,chr16_PATERNAL,chr17_PATERNAL,chr18_PATERNAL,chr19_PATERNAL,chr20_PATERNAL,chr21_PATERNAL,chr22_PATERNAL"
        CALLING_INTERVALS="chrM,chr1_PATERNAL,chr2_PATERNAL,chr3_PATERNAL,chr4_PATERNAL,chr5_PATERNAL,chr6_PATERNAL,chr7_PATERNAL,chr8_PATERNAL,chr9_PATERNAL,chr10_PATERNAL,chr11_PATERNAL,chr12_PATERNAL,chr13_PATERNAL,chr14_PATERNAL,chr15_PATERNAL,chr16_PATERNAL,chr17_PATERNAL,chr18_PATERNAL,chr19_PATERNAL,chr20_PATERNAL,chr21_PATERNAL,chr22_PATERNAL,chrEBV,chrY_PATERNAL"
        ;;
      chrM)
        BQSR_INTERVALS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
        CALLING_INTERVALS="chrM,chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
        ;;
      1)
        BQSR_INTERVALS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
        CALLING_INTERVALS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT"
        ;;
      chr1)
        BQSR_INTERVALS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
        CALLING_INTERVALS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
        ;;
      chr22)
        BQSR_INTERVALS="chr22"
        CALLING_INTERVALS="chr22"
        ;;
      *)
        echo "ERROR: unknown first reference chromosome"
        exit 1
        ;;
    esac

    # Extract the samples
    tumor_sm=$(samtools samples "$tumor_deduped" | cut -f 1 | head -n 1)
    normal_sm=$(samtools samples "$normal_deduped" | cut -f 1 | head -n 1)

    # Variant calling with TNScope
    sentieon driver -r "~{ref_fasta}" \
      -i "~{aligned_reads}" \
      ${normal_deduped:+-i "~{normal_reads}"} \
      --interval "$CALLING_INTERVALS" \
      --algo TNscope\
        --tumor_sample "$tumor_sm" \
        ${normal_sm:+--normal_sample "$normal_sm"} \
        "sample_tnscope.vcf.gz" \

    wait
    unset http_proxy
    exit 0
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: sentieon_docker
    memory: memory
    cpu: n_threads
  }
  output {
    # VCF outputs
    File Sentieon_calls_vcf = "sample_tnscope.vcf.gz"
    File Sentieon_calls_vcf_tbi = "sample_tnscope.vcf.gz.tbi"
  }
task split_contigs {
  input {
    File ref_fasta_index
    Int chunk_size = 75000000
    Int threads
  }

  Float file_size = ceil(size(ref_fasta_index, "GB") + 10)

  command <<<
  set -euxo pipefail

  bedtools --version

  echo "Splitting contigs for ~{ref_fasta_index}"
  bedtools makewindows -g ~{ref_fasta_index} -w ~{chunk_size} > contigs.bed
  grep -v -E "random|chrUn|chrM|chrEBV" contigs.bed > noalt.bed
  # Split the contig bed files into one file for each line
  split -l 1 noalt.bed contigs_split.
  # Add .bed to all the contigs_split file
  for file in $(ls contigs_split.*); do mv $file $file.bed; done
  >>>

  output {
    Array[File] contigs = glob("contigs_split.*.bed")
  }

  runtime {
    docker: "860660336427.dkr.ecr.us-east-1.amazonaws.com/hifisomatic:bedtools"
    cpu: threads
    memory: "~{threads * 4} GB"
    disk: file_size + " GB"
    maxRetries: 2
    preemptible: 1
  }
}

task DeepSomatic {
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

    Array[File] contigs

    # Execution
    String n_threads = "16"
    Int preemptible_tries = 3
    String deepsomatic_docker
  }
  command <<<
    set -xv
    set -exvuo pipefail

    scatter (ctg in contigs){

    run_deepsomatic \
    --model_type=WGS \
    --ref="~{ref_fasta}" \
    --reads_normal="~{normal_reads}" \
    --reads_tumor="~{aligned_reads}" \
    --output_vcf=sample_deepsomatic.vcf.gz \
    --output_gvcf=sample_deepsomatic.g.vcf.gz \
    --sample_name_tumor="~{tumor_name}" \
    --sample_name_normal="~{normal_name}" \
    --num_shards="~{n_threads}" \
    --logging_dir=logs

    ls -l

    wait
    unset http_proxy
    exit 0 
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: deepsomatic_docker
    cpu: n_threads
    memory: "~{n_threads * 4} GB"
  }
  output {
    File DeepSomatic_calls_vcf = "sample_deepsomatic.vcf.gz"
    File DeepSomatic_calls_vcf_tbi = "sample_deepsomatic.vcf.gz.tbi"
    File DeepSomatic_calls_gvcf = "sample_deepsomatic.g.vcf.gz"
    File DeepSomatic_calls_gvcf_tbi = "sample_deepsomatic.g.vcf.gz.tbi"
    File DeepSomatic_visual_report = "sample_deepsomatic.visual_report.html"
  }
}
