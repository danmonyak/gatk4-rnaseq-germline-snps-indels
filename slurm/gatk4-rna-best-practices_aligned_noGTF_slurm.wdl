## Copyright Broad Institute, 2019
##
## Workflows for processing RNA data for germline short variant discovery with GATK (v4) and related tools
##
## Requirements/expectations :
## - BAM 
##
## Output :
## - A BAM file and its index.
## - A VCF file and its index. 
## - A Filtered VCF file and its index. 
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs. 
 
 workflow RNAseq_array {

	File inputAlignedBam
	String sampleName = basename(inputAlignedBam,".bam")

	File refFasta
	File refFastaIndex
	File refDict

	String? gatk_path_override
	String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])
	
	Array[File] knownVcfs
	Array[File] knownVcfsIndices
	
	File dbSnpVcf
	File dbSnpVcfIndex

	Int? minConfidenceForVariantCalling

 	## Directory of scatter interval list files from previous run
	String ScatterIntervalListDir
	#Array[File] ScatterIntervalList_manual = glob(ScatterIntervalListDir + "/*.interval_list")
	call glob_task{input: silDir=ScatterIntervalListDir}

	## Optional user optimizations
	Int? haplotypeScatterCount
	Int scatterCount = select_first([haplotypeScatterCount, 6])

	Int? preemptible_tries
	Int preemptible_count = select_first([preemptible_tries, 3])

	call MarkDuplicates {
	    input:
	        input_bam = inputAlignedBam,
	        base_name = sampleName + ".dedupped",
	        preemptible_count = preemptible_count,
	        #docker = gatk4_docker,
	        gatk_path = gatk_path
	}


	call SplitNCigarReads {
	input:
	    input_bam = MarkDuplicates.output_bam,
	    input_bam_index = MarkDuplicates.output_bam_index,
	    base_name = sampleName + ".split",
	    ref_fasta = refFasta,
	    ref_fasta_index = refFastaIndex,
	    ref_dict = refDict,
	    preemptible_count = preemptible_count,
	    #docker = gatk4_docker,
	    gatk_path = gatk_path
	}


	call BaseRecalibrator {
		input:
			input_bam = SplitNCigarReads.output_bam,
			input_bam_index = SplitNCigarReads.output_bam_index,
			recal_output_file = sampleName + ".recal_data.csv",
  			dbSNP_vcf = dbSnpVcf,
  			dbSNP_vcf_index = dbSnpVcfIndex,
  			known_indels_sites_VCFs = knownVcfs,
  			known_indels_sites_indices = knownVcfsIndices,
  			ref_dict = refDict,
  			ref_fasta = refFasta,
  			ref_fasta_index = refFastaIndex,
  			preemptible_count = preemptible_count,
			#docker = gatk4_docker,
			gatk_path = gatk_path
	}

	call ApplyBQSR {
		input:
			input_bam =  SplitNCigarReads.output_bam,
			input_bam_index = SplitNCigarReads.output_bam_index,
			base_name = sampleName + ".aligned.duplicates_marked.recalibrated",
			ref_fasta = refFasta,
			ref_fasta_index = refFastaIndex,
			ref_dict = refDict,
			recalibration_report = BaseRecalibrator.recalibration_report,
			preemptible_count = preemptible_count,
			#docker = gatk4_docker,
			gatk_path = gatk_path
	}

	
	scatter (interval in glob_task.out) {
        call HaplotypeCaller {
		input:
		input_bam = ApplyBQSR.output_bam,
		input_bam_index = ApplyBQSR.output_bam_index,
		base_name = sampleName + ".hc",
		interval_list = interval,
		ref_fasta = refFasta,
		ref_fasta_index = refFastaIndex,
		ref_dict = refDict,
		dbSNP_vcf = dbSnpVcf,
		dbSNP_vcf_index = dbSnpVcfIndex,
		stand_call_conf = minConfidenceForVariantCalling,
		preemptible_count = preemptible_count,
		#docker = gatk4_docker,
		gatk_path = gatk_path
        }

		File HaplotypeCallerOutputVcf = HaplotypeCaller.output_vcf
		File HaplotypeCallerOutputVcfIndex = HaplotypeCaller.output_vcf_index
	}

	call MergeVCFs {
	input:
	    input_vcfs = HaplotypeCallerOutputVcf,
	    input_vcfs_indexes =  HaplotypeCallerOutputVcfIndex,
	    output_vcf_name = sampleName + ".g.vcf.gz",
	    preemptible_count = preemptible_count,
	    #docker = gatk4_docker,
	    gatk_path = gatk_path
	}

	call VariantFiltration {
	input:
		input_vcf = MergeVCFs.output_vcf,
		input_vcf_index = MergeVCFs.output_vcf_index,
		base_name = sampleName + ".variant_filtered.vcf.gz",
		ref_fasta = refFasta,
		ref_fasta_index = refFastaIndex,
		ref_dict = refDict,
		preemptible_count = preemptible_count,
		#docker = gatk4_docker,
		gatk_path = gatk_path
	}

	output {
		File recalibrated_bam = ApplyBQSR.output_bam
		File recalibrated_bam_index = ApplyBQSR.output_bam_index
		File merged_vcf = MergeVCFs.output_vcf
		File merged_vcf_index = MergeVCFs.output_vcf_index
		File variant_filtered_vcf = VariantFiltration.output_vcf
		File variant_filtered_vcf_index = VariantFiltration.output_vcf_index
	}
}

task MarkDuplicates {

 	File input_bam
 	String base_name

	String gatk_path

	#String docker
 	Int preemptible_count

 	command <<<
 	    ${gatk_path} \
 	        MarkDuplicates \
 	        --INPUT ${input_bam} \
 	        --OUTPUT ${base_name}.bam  \
 	        --CREATE_INDEX true \
 	        --VALIDATION_STRINGENCY SILENT \
 	        --METRICS_FILE ${base_name}.metrics
 	>>>

 	output {
 		File output_bam = "${base_name}.bam"
 		File output_bam_index = "${base_name}.bai"
 		File metrics_file = "${base_name}.metrics"
 	}

	runtime {
		requested_memory_mb_per_core: 200
	}
}

task SplitNCigarReads {

  File input_bam
  File input_bam_index
  String base_name

  File ref_fasta
  File ref_fasta_index
  File ref_dict

	String gatk_path
	#String docker
        Int preemptible_count

    command <<<
        ${gatk_path} \
                SplitNCigarReads \
                -R ${ref_fasta} \
                -I ${input_bam} \
                -O ${base_name}.bam \
		-sequence-dictionary ${ref_dict}
    >>>

        output {
                File output_bam = "${base_name}.bam"
                File output_bam_index = "${base_name}.bai"
        }

    runtime {
        #disks: "local-disk " + sub(((size(input_bam,"GB")+1)*5 + size(ref_fasta,"GB")),"\\..*","") + " HDD"
        #docker: docker
        #memory: "5 GB"
        #preemptible: preemptible_count
    }
}

task BaseRecalibrator {

    File input_bam
    File input_bam_index
    String recal_output_file

    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String gatk_path

    #String docker
    Int preemptible_count

    command <<<
        ${gatk_path} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --use-original-qualities \
            -O ${recal_output_file} \
            -known-sites ${dbSNP_vcf} \
            -known-sites ${sep=" --known-sites " known_indels_sites_VCFs} \
	    -sequence-dictionary ${ref_dict}
    >>>

    output {
        File recalibration_report = recal_output_file
    }

    runtime {
        #memory: "5 GB"
        #disks: "local-disk " + sub((size(input_bam,"GB")*3)+30, "\\..*", "") + " HDD"
        #docker: docker
        #preemptible: preemptible_count
    }
}


task ApplyBQSR {

    File input_bam
    File input_bam_index
    String base_name
    File recalibration_report

    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String gatk_path

    #String docker
    Int preemptible_count

    command <<<
        ${gatk_path} \
            --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R ${ref_fasta} \
            -I ${input_bam} \
            --use-original-qualities \
            -O ${base_name}.bam \
            --bqsr-recal-file ${recalibration_report} \
	    -sequence-dictionary ${ref_dict}
    >>>

    output {
        File output_bam = "${base_name}.bam"
        File output_bam_index = "${base_name}.bai"
    }

    runtime {
        #memory: "5 GB"
        #disks: "local-disk " + sub((size(input_bam,"GB")*4)+30, "\\..*", "") + " HDD"
        #preemptible: preemptible_count
        #docker: docker
    }
}

task HaplotypeCaller {

	File input_bam
	File input_bam_index
	String base_name

	File interval_list

	File ref_dict
	File ref_fasta
	File ref_fasta_index

	File dbSNP_vcf
	File dbSNP_vcf_index

	String gatk_path
	#String docker
	Int preemptible_count

	Int? stand_call_conf

	command <<<
		${gatk_path} --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
		HaplotypeCaller \
		-R ${ref_fasta} \
		-I ${input_bam} \
		-L ${interval_list} \
		-O ${base_name}.vcf.gz \
		-dont-use-soft-clipped-bases \
		--standard-min-confidence-threshold-for-calling ${default=20 stand_call_conf} \
		--dbsnp ${dbSNP_vcf} \
	    	-sequence-dictionary ${ref_dict}
	>>>

	output {
		File output_vcf = "${base_name}.vcf.gz"
		File output_vcf_index = "${base_name}.vcf.gz.tbi"
	}

	runtime {
		#docker: docker
		#memory: "7 GB"
		#disks: "local-disk " + sub((size(input_bam,"GB")*2)+30, "\\..*", "") + " HDD"
		#preemptible: preemptible_count
	}
}

task VariantFiltration {

	File input_vcf
	File input_vcf_index
	String base_name

 	File ref_dict
 	File ref_fasta
 	File ref_fasta_index

	String gatk_path
	#String docker
 	Int preemptible_count

	command <<<
		 ${gatk_path} \
		    VariantFiltration \
			--R ${ref_fasta} \
			--V ${input_vcf} \
			--window 35 \
			--cluster 3 \
			--filter-name "FS" \
			--filter "FS > 30.0" \
			--filter-name "QD" \
			--filter "QD < 2.0" \
			-O ${base_name} \
	    		-sequence-dictionary ${ref_dict}
	>>>

	output {
    	File output_vcf = "${base_name}"
    	File output_vcf_index = "${base_name}.tbi"
	}

	runtime {
		#docker: docker
		#memory: "5 GB"
		#disks: "local-disk " + sub((size(input_vcf,"GB")*2)+30, "\\..*", "") + " HDD"
		#preemptible: preemptible_count
	}
}

task MergeVCFs {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name

    Int? disk_size = 5

    String gatk_path

    #String docker
    Int preemptible_count

    # Using MergeVcfs instead of GatherVcfs so we can create indices
    # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
    command <<<
        ${gatk_path} --java-options "-Xms2000m"  \
            MergeVcfs \
            --INPUT ${sep=' --INPUT ' input_vcfs} \
            --OUTPUT ${output_vcf_name}
    >>>

    output {
        File output_vcf = output_vcf_name
        File output_vcf_index = "${output_vcf_name}.tbi"
    }

    runtime {
        #memory: "5 GB"
        #disks: "local-disk " + disk_size + " HDD"
        #docker: docker
        #preemptible: preemptible_count
    }
}

task glob_task {
    String silDir
    command {}
    output { Array[File] out = glob(silDir + "/out/*/*.interval_list") }
}
