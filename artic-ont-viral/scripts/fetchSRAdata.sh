#!/bin/bash

#SBATCH --partition batch
#SBATCH -w compute06
#SBATCH -n 8
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --job-name=fetchings
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=J.Juma@cgiar.org

# clear environment
module purge

# load required module(s)
module load sra-tools/3.0.0


# WORK_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023"
# DATA_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/data"
# OUTPUT_DIR="${DATA_DIR}/dataset-001"

# accessions_dataset="${HOME}/trainings/africacdc-ilri-aslm-2023/metadata/accessions-001.txt"

if [ $# -eq 0 ]; then 
	  echo -e "$0 : You must provide two arguments (1 the path to the accessions file and 2, the path to the output directory" 
	  exit 1 
fi  

if [ ! -f $1 ]; then
  echo -e "Please provide a file with list of accessions as one accession per line"
fi

if [ ! -d $2 ];then
    mkdir -p $2
fi




# fetch raw reads from SRA/ENA
while IFS=$'\n' read -r line;
do
    if [[ "${line}" =~ \#.* ]]; then
        echo -e "comment line:$line"
    else
      output_dir=$2/${line}
      

      outfile="${output_dir}/${line}.fastq"
      if [ ! -f "${outfile}" ]; then
          echo -e "fetching read data for accession: ${line}"
          fasterq-dump --split-3 --threads 4 --outdir ${output_dir} ${line}
      fi
    fi
done < $1
# done < ${accessions_dataset}