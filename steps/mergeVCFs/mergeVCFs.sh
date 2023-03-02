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

gatkDir="/hpc/group/ryserlab/working/gatk4/gatk4-rnaseq-germline-snps-indels/cromwell-executions/RNAseq/2617d9f2-b516-4dc6-b4e7-990f2ccdc12e/call-HaplotypeCaller"
manualDir="/hpc/group/ryserlab/working/gatk4/gatk4-rnaseq-germline-snps-indels/steps/haplotypeCaller"

vcfFilename="SRR17843648_CR.hc.vcf.gz"

gatk_path="/hpc/group/ryserlab/gatk-4.3.0.0/gatk"

${gatk_path} --java-options "-Xms2000m"  \
            MergeVcfs \
            --INPUT $gatkDir/shard-0/execution/$vcfFilename \
            --INPUT $gatkDir/shard-1/execution/$vcfFilename \
            --INPUT $gatkDir/shard-2/execution/$vcfFilename \
            --INPUT $manualDir/$vcfFilename \
            --INPUT $gatkDir/shard-4/execution/$vcfFilename \
            --INPUT $gatkDir/shard-5/execution/$vcfFilename \
            --OUTPUT ${output_vcf_name}

echo "Done!"
