#!/bin/bash
#SBATCH --job-name=haplotypeCaller
#SBATCH -o haplotypeCaller.out
#SBATCH -e haplotypeCaller.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --mail-user=danmonyak@gmail.com
#SBATCH --mail-type=END

sampleName="SRR17843648_unmapped"
#BSQR_base_name = sampleName + ".aligned.duplicates_marked.recalibrated"


gatk_path="/hpc/group/ryserlab/gatk-4.3.0.0/gatk"
ref_fasta="/hpc/group/ryserlab/working/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
input_bam="/hpc/group/ryserlab/working/gatk4/gatk4-rnaseq-germline-snps-indels/cromwell-executions/RNAseq/2617d9f2-b516-4dc6-b4e7-990f2ccdc12e/call-ApplyBQSR/execution/SRR17843648_CR.aligned.duplicates_marked.recalibrated.bam"
interval_list="/hpc/group/ryserlab/working/gatk4/gatk4-rnaseq-germline-snps-indels/cromwell-executions/RNAseq/2617d9f2-b516-4dc6-b4e7-990f2ccdc12e/call-ScatterIntervalList/execution/out/temp_0004_of_6/4scattered.interval_list"
base_name=sampleName".hc"
dbSNP_vcf="/hpc/group/ryserlab/working/gatk4/gatk4-rnaseq-germline-snps-indels/realdata/Homo_sapiens_assembly38.dbsnp138.vcf"
ref_dict="/hpc/group/ryserlab/working/refdata-gex-GRCh38-2020-A/fasta/genome.dict"

${gatk_path} --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
		HaplotypeCaller \
		-R ${ref_fasta} \
		-I ${input_bam} \
		-L ${interval_list} \
		-O ${base_name}.vcf.gz \
		-dont-use-soft-clipped-bases \
		--standard-min-confidence-threshold-for-calling 20 \
		--dbsnp ${dbSNP_vcf} \
	    	-sequence-dictionary ${ref_dict}

echo "Done!"
