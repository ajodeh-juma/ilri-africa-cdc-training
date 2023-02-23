#!/bin/bash

#SBATCH --partition batch
#SBATCH -w compute06
#SBATCH -n 8
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --job-name=nanoqc
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=J.Juma@cgiar.org


DATA_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/data/rvfv"
OUTPUT_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/output/rvfv/pycoqc"

if [ ! -d "${OUTPUT_DIR}" ]; then
    mkdir -p "${OUTPUT_DIR}"
fi

echo -e "Quality control of data using PycoQC"
pycoQC -f ${DATA_DIR}/sequencing_summary_FAR90564_0e97e5b0.txt -o ${OUTPUT_DIR}/rvfv.html