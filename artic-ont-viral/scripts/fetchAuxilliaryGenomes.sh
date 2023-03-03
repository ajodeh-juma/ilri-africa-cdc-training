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
module load edirect/7.80

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

infile=$1
outdir=$2

# fetch genomes from NCBI nucleotide database

# while IFS=$'\n' read -r line;
# do
#     echo -e "Retrieving reference genome for accession: ${line}"
#     esearch -db nucleotide -query $line | efetch -format fasta > ${outdir}/${line}.fasta
# done < $1


IFS=$'\n'; 
for accession in $(cat ${infile}); 
do
    outfile=${outdir}/${accession}.fasta
    
    if [ ! -f "${outfile}" ]; then
        echo -e "Retrieving reference genome for accession: ${accession}"
        esearch -db nucleotide -query $accession | efetch -db nucleotide -format fasta > ${outdir}/${accession}.fasta;
    fi
done

# cat $1 | while read line 
# do
#     echo -e "Retrieving reference genome for accession: ${line}"
#     esearch -db nucleotide -query "$line" | efetch -format fasta > $outdir/$line.fasta
# done
