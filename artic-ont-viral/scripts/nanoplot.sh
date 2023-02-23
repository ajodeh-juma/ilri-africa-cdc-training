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
OUTPUT_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/output/dataset-002/nanoplot"

if [ ! -d "${OUTPUT_DIR}" ]; then
    mkdir -p "${OUTPUT_DIR}"
fi


for fastq in "${DATA_DIR}"/*/*.fastq;
do
    sample=$(basename $fastq | cut -d . -f 1)
    echo -e "quality check with nanoplot for sample: ${sample}"
    NanoPlot -t 8 --fastq $fastq --outdir "${OUTPUT_DIR}" --prefix $sample --tsv_stats --plots kde hex
done