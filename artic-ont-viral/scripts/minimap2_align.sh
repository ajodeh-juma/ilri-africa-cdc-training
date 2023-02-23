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
module purge

# load modules
module load minimap2/2.22
module load samtools/1.9

# define inputs/outputs
#WORK_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023"
DATA_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/output/dataset-002/artic-guppyplex"
OUTPUT_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/output/dataset-002/minimap2"
reference="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/genomes/dengue/DENV2/DENV2.fasta"


if [ ! -d "${OUTPUT_DIR}" ]; then
    mkdir -p "${OUTPUT_DIR}"
fi

cd "${OUTPUT_DIR}"

for fastq in ${DATA_DIR}/*.fastq.gz;
do
    sample=$(basename $fastq | cut -d . -f 1)
    echo -e "aligning reads with minimap2 on: ${sample}"
    minimap2 -a -x map-ont -t 8 ${reference} ${fastq} | samtools view -bS -F 4 - | samtools sort -o ${sample}.sorted.bam -
    samtools index ${sample}.sorted.bam
done

       
