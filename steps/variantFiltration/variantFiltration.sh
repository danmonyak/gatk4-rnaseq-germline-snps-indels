#!/bin/bash
#SBATCH --job-name=variantFiltration
#SBATCH -o variantFiltration.out
#SBATCH -e variantFiltration.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=20G
#SBATCH --mail-user=danmonyak@gmail.com
#SBATCH --mail-type=END

sampleName="SRR17843648_CR"

gatk_path="/hpc/group/ryserlab/gatk-4.3.0.0/gatk"
ref_fasta="/hpc/group/ryserlab/working/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
input_vcf="/hpc/group/ryserlab/working/gatk4/gatk4-rnaseq-germline-snps-indels/steps/mergeVCFs/SRR17843648_CR.g.vcf.gz"
base_name=$sampleName".variant_filtered.vcf.gz"
ref_dict="/hpc/group/ryserlab/working/refdata-gex-GRCh38-2020-A/fasta/genome.dict"

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

echo "Done!"
