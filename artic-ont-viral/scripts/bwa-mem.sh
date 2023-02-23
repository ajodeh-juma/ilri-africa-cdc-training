#!/bin/bash

#SBATCH --partition batch
#SBATCH -w compute06
#SBATCH -n 8
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --job-name=align_minimap2
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=J.Juma@cgiar.org

# load modules
#module load bwa/0.7.17
#module load samtools/1.9

#--------------------------------------------------------------------------#
#          build genome index and align/map reads to the genome            #
#--------------------------------------------------------------------------#


# specify the data directory
#DATA_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/output-dir/dataset-2/artic-guppyplex"
DATA_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/output/dataset-002/artic-guppyplex"


# create output directory if not exists
OUTPUT_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/output/dataset-002/bwa"


# index the reference genome
reference="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/genomes/dengue/DENV2/DENV2.fasta"
# reference="${HOME}/trainings/africacdc-ilri-aslm-2023/genomes/ebola/EBOV.fasta"


ext="${reference##*.}"
basename=$(basename ${reference} \.${ext})



# create output directory if not exists
# INDEX="${HOME}/trainings/africacdc-ilri-aslm-2023/genomes/ebola"
INDEX="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/genomes/dengue/DENV2"


if [ ! -d "${INDEX}" ]; then
 mkdir -p $INDEX
fi


bwa index -p ${INDEX}/${basename} $reference


# align reads to the genome




if [ ! -d "${OUTPUT_DIR}" ]; then
    mkdir -p "${OUTPUT_DIR}"
fi

# change to output directory
cd "${OUTPUT_DIR}"


# loop through the files in the data directory and run bowtie2
for fastq in ${DATA_DIR}/*.fastq.gz;
do
    sample=$(basename ${fastq} | cut -d . -f 1)
    echo -e "aligning reads with bwa-mem on: ${sample}"
    bwa mem -x ont2d -t 8 ${INDEX}/${basename} $fastq | samtools view -bS -F 4 - | samtools sort -o ${sample}.sorted.bam -
done