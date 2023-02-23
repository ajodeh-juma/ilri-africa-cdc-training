#!/bin/bash

#SBATCH --partition batch
#SBATCH -w compute06
#SBATCH -n 8
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --job-name=fetch_reference
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=J.Juma@cgiar.org

# clear environment
module purge

# load required module(s)
module load edirect/7.80


WORK_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023"
REFERENCE_DIRS="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/genomes/dengue"
PRIMER_SCHEMES_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/artic-ont-viral/primer-schemes/"

mkdir -p "${REFERENCE_DIRS}"/{DENV1,DENV2,DENV3,DENV4}


accessions="NC_001477.1 NC_001474.2 NC_001475.2 NC_002640.1"

for accession in $accessions;
do
  if [ $accession == "NC_001477.1" ]; then
    prefix="DENV1"
    echo -e "Retrieving reference genome for accession: $accession"
    esearch -db nucleotide -query "$accession" | efetch -format fasta > "$REFERENCE_DIRS"/$prefix/$prefix.fasta
    cp "$REFERENCE_DIRS"/$prefix/$prefix.fasta "$PRIMER_SCHEMES_DIR"/$prefix/$prefix.reference.fasta
  fi

  if [ "$accession" == "NC_001474.2" ]; then
    prefix="DENV2"
    echo -e "Retrieving reference genome for accession: $accession"
    # esearch -db nucleotide -query "$accession" | efetch -format fasta > "$REFERENCE_DIRS"/$prefix/$prefix.fasta
    cp "$REFERENCE_DIRS"/$prefix/$prefix.fasta "$PRIMER_SCHEMES_DIR"/$prefix/$prefix.reference.fasta
  fi

  if [ "$accession" == "NC_001475.2" ]; then
    prefix="DENV3"
    echo -e "Retrieving reference genome for accession: $accession"
    esearch -db nucleotide -query "$accession" | efetch -format fasta > "$REFERENCE_DIRS"/$prefix/$prefix.fasta
    cp "$REFERENCE_DIRS"/$prefix/$prefix.fasta "$PRIMER_SCHEMES_DIR"/$prefix/$prefix.reference.fasta
  fi

  if [ "$accession" == "NC_002640.1" ]; then
    prefix="DENV4"
    echo -e "Retrieving reference genome for accession: $accession"
    esearch -db nucleotide -query "$accession" | efetch -format fasta > "$REFERENCE_DIRS"/$prefix/$prefix.fasta
    cp "$REFERENCE_DIRS"/$prefix/$prefix.fasta "$PRIMER_SCHEMES_DIR"/$prefix/$prefix.reference.fasta
  fi
done

#if [ ! -d ${REFERENCE_DIR} ]; then
#    mkdir ${REFERENCE_DIR}
#fi
#
#cd "${REFERENCE_DIR}"

#wget -c https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Dengue_virus/latest_assembly_versions/GCF_000871845.1_ViralProj20183/GCF_000871845.1_ViralProj20183_genomic.fna.gz
#wget -c https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Dengue_virus/latest_assembly_versions/GCF_000871845.1_ViralProj20183/GCF_000871845.1_ViralProj20183_genomic.gff.gz
#
#
#gunzip *.gz
#
#mv GCF_000871845.1_ViralProj20183_genomic.fna DENV2.fasta
#mv GCF_000871845.1_ViralProj20183_genomic.gff DENV2.gff

#esearch -db nucleotide -query "NC_001474.2" | efetch -format fasta > DENV2.fasta
#esearch -db nucleotide -query "NC_001474.2" | efetch -format gff > DENV2.gff


# ebola virus
#
#REFERENCE_DIR="${HOME}/trainings/africacdc-ilri-aslm-2023/genomes/ebola"
#
#if [ ! -d ${REFERENCE_DIR} ]; then
#    mkdir ${REFERENCE_DIR}
#fi
#
#cd "${REFERENCE_DIR}"

# esearch -db nucleotide -query "KR063671.1" | efetch -format fasta > EBOV.fasta

# accessory files required for amplicon sequence data analysis
# modules: artic-tools
#primer_scheme_dir="${HOME}/trainings/africacdc-ilri-aslm-2023/primer-schemes/ZaireEbola"
#
#if [ ! -d ${primer_scheme_dir} ]; then
#    mkdir ${primer_scheme_dir}
#fi
#artic-tools get_scheme --schemeVersion 1 --outDir ${primer_scheme_dir} ebola

