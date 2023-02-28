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
module load artic/1.2.3

# define inputs/outputs
DATA_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/data/dataset-002"
OUTPUT_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/output/dataset-002/artic-guppyplex"


if [ ! -d "${OUTPUT_DIR}" ]; then
    mkdir -p "${OUTPUT_DIR}"
fi

cd "${OUTPUT_DIR}"


for dir in "${DATA_DIR}"/*;
do
  for fastq in "${dir}"/*;
  do
    sample=$(basename $fastq | cut -d . -f 1)
    output="${OUTPUT_DIR}/${sample}.fastq"
    log_file="${OUTPUT_DIR}/${sample}-guppyplex.log"
    echo -e "Running artic guppyplex to filter reads in ${sample}"
    artic guppyplex \
      --skip-quality-check \
      --min-length 356 \
      --max-length 588 \
      --directory "${dir}" \
      --output $output >>$log_file 2>&1
    pigz --no-name --processes 8 $output
  done
done

       
