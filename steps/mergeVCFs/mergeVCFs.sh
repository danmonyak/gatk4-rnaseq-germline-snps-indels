#!/bin/bash
#SBATCH --job-name=mergeVCFs
#SBATCH -o mergeVCFs.out
#SBATCH -e mergeVCFs.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH --mail-user=danmonyak@gmail.com
#SBATCH --mail-type=END

sampleName="SRR17843648_unmapped"

input_vcfs=$HaplotypeCallerOutputVcf
output_vcf_name=$sampleName".g.vcf.gz"

/hpc/group/ryserlab/working/gatk4/gatk4-rnaseq-germline-snps-indels

gatk_path="/hpc/group/ryserlab/gatk-4.3.0.0/gatk"

${gatk_path} --java-options "-Xms2000m"  \
            MergeVcfs \
            --INPUT ${sep=' --INPUT ' input_vcfs} \
            --OUTPUT ${output_vcf_name}

echo "Done!"
