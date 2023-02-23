#!/bin/bash

#SBATCH --partition batch
#SBATCH -w compute06
#SBATCH -n 8
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --job-name=nanoqc
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=J.Juma@cgiar.org


#DATA_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/data/dataset-2"
DATA_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/data/dataset-002"
OUTPUT_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/output/dataset-002/fastqc"

if [ ! -d "${OUTPUT_DIR}" ]; then
    mkdir -p "${OUTPUT_DIR}"
fi


for fastq in "${DATA_DIR}"/*/*.fastq;
do
    sample=$(basename $fastq | cut -d . -f 1)
    echo -e "quality check with FastQC for sample: ${sample}"
    fastqc -o ${OUTPUT_DIR} -f fastq --extract $fastq
done


MULTIQC_OUTPUT_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/output/dataset-002/multiqc"

echo -e "summarizing quality reports with MultiQC from FastQC output: ${OUTPUT_DIR}"
multiqc -o ${MULTIQC_OUTPUT_DIR} ${OUTPUT_DIR}