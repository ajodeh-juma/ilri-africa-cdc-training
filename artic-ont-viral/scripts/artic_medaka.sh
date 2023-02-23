#!/bin/bash

#SBATCH --partition batch
#SBATCH -w compute06
#SBATCH -n 8
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --job-name=align_minimap2
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=J.Juma@cgiar.org

# clear environment
# module purge

# load modules
#module load minimap2/2.22
#module load samtools/1.9

SRC_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/scripts"

# define inputs/outputs
DATA_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/output/dataset-002/minimap2"
OUTPUT_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/output/dataset-002/medaka"
primer_bed="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/primer-schemes/DENV2/DENV2.primer.bed"
reference="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/genomes/dengue/DENV2/DENV2.fasta"


medaka_model="r941_min_high_g360"

if [ ! -d "${OUTPUT_DIR}" ]; then
    mkdir -p "${OUTPUT_DIR}"
fi

cd "${OUTPUT_DIR}"


# samtools sort -o ${prefix}.sorted.bam -
# Trim alignments from an amplicon scheme
for bam in "${DATA_DIR}"/*.bam;
do
    echo -e "$bam"
    sample=$(basename "${bam}" | cut -d . -f 1)
    # log_file="${OUTPUT_DIR}/${sample}_trim_align.log"

    artic-tools align_trim \
      --normalise 200 \
      --start \
      --remove-incorrect-pairs \
      --report "$sample.alignreport.txt" \
      <${bam} \
      ${primer_bed} | samtools sort - -o ${sample}.trimmed.rg.sorted.bam

      samtools index "${sample}.trimmed.rg.sorted.bam"


    artic-tools align_trim \
      --normalise 200 \
      --remove-incorrect-pairs \
      --report "${sample}.alignreport.txt" \
      <${sample}.trimmed.rg.sorted.bam \
      ${primer_bed} | samtools sort - -o ${sample}.primertrimmed.rg.sorted.bam

    samtools index "${sample}.primertrimmed.rg.sorted.bam"
 done




# collect the primer pools, values to be used as read groups in medaka consensus step
pools=(1 2)


### generate consensus for each read group
for p in ${pools[*]};
do
 for bamfile in ${OUTPUT_DIR}/*.primertrimmed.rg.sorted.bam;
 do
   sample=$(basename "${bamfile}" | cut -d . -f 1)
   output_prefix=${sample}.${p}
   hdf=${output_prefix}.hdf
   vcf=${output_prefix}.vcf

   if [ ! -f "${hdf}" ];then
    echo -e "Medaka consensus alignment for read group ${p} on bam file ${sample}"
    medaka consensus --model ${medaka_model} --threads 8 --RG ${p} $bamfile ${hdf}
   fi

   if [ ! -f "${vcf}" ];then
    echo -e "Medaka variant calling for read group ${p} on bam file ${sample}"
    medaka variant ${reference} ${hdf} ${vcf}
   fi
 done
done


# merge, compress and index vcf files
for i in ./*.1.vcf*;
do
  sample=$(basename "${i}" | cut -d . -f 1)
  merged=${sample}.merged.vcf
  j="${i%*.1.vcf}.2.vcf"
    if [[ -e "$j" ]]; then
      echo "merging vcfs in $i and $j files"
      artic_vcf_merge ${sample} ${primer_bed} 2>${sample}.primersitereport.txt 2:${j} 1:${i}
      echo "compressing merged vcf file $merged"
      bgzip -f ${sample}.merged.vcf
      echo "indexing compressed vcf file $merged"
      tabix -f -p vcf ${sample}.merged.vcf.gz
    fi
done

# call SNVs, filter, compress and filter vcfs
for bamfile in ./*.primertrimmed.rg.sorted.bam;
do
  sample=$(basename "${bamfile}" | cut -d . -f 1)
  echo "calling SNVs in $bamfile"
  longshot -P 0 \
    -F -A --no_haps \
    --bam ${bamfile} \
    --ref ${reference} \
    --out ${sample}.merged.vcf \
    --potential_variants ${sample}.merged.vcf.gz


  echo -e "filtering variants in sample : $sample"
  artic_vcf_filter --medaka ${sample}.merged.vcf ${sample}.pass.vcf ${sample}.fail.vcf

  echo "compressing passed variants in vcf sample : $sample"
  bgzip -f ${sample}.pass.vcf

  echo "indexing passed variants in vcf sample : $sample"
  tabix -p vcf ${sample}.pass.vcf.gz
done


# compute coverage, mask low depth regions and generate consensus
for bamfile in ./*.primertrimmed.rg.sorted.bam;
do
  sample=$(basename "${bamfile}" | cut -d . -f 1)

  echo "computing coverage per position in bam file $bamfile"
  artic_make_depth_mask \
    --store-rg-depths \
    ${reference} \
    ${bamfile} \
    ${sample}.coverage_mask.txt


  echo "masking regions with low coverage in $sample"
  artic_mask \
    ${reference} \
    ${sample}.coverage_mask.txt \
    ${sample}.fail.vcf ${sample}.preconsensus.fasta


  echo "generating consensus sequence in $sample"
  bcftools consensus \
    -f ${sample}.preconsensus.fasta \
    ${sample}.pass.vcf.gz \
    -m ${sample}.coverage_mask.txt \
    -o ${sample}.consensus.fasta


  echo "renaming header for consensus sequence in $sample"
  artic_fasta_header ${sample}.consensus.fasta "${sample}/ARTIC/medaka"


  # align to reference
  cat ${sample}.consensus.fasta ${reference} > ${sample}.muscle.in.fasta
  muscle -in ${sample}.muscle.in.fasta -out ${sample}.muscle.out.fasta




done


