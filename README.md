---
title: README.md
tags: ["SARS-CoV-2", "Genomics", "Bioinformatics", "Metadata", "Linux", "Analysis"]
---
# **Building capacity in SARS-CoV-2 genomics in Africa**
---
###### ***Trainers***: [John Juma](https://github.com/ajodeh-juma), [Ouso Daniel](https://github.com/ousodaniel) & [Gilbert Kibet](https://github.com/kibet-gilbert)
---

- [Introduction](#introduction)
- [Scope](#scope)
- [Background](#background)
- [Prerequisite](#prerequisite)
- [Set-Up](#setup)
- [Preparations](#preparations)
    - [Log into the HPC](#log-into-the-HPC)
    - [Project Organisation](#project-organisation)
    - [Data retrieval and integrity checks](#data-retrieval-and-integrity-checks)
- [Analysis](#analysis)
    - [Loading modules](#loading-modules)
    - [Prepare the reference genome](#prepare-the-reference-genome)
    - [Quality assessment](#quality-assessment)
    - [Quality and Adapter filtering](#quality-and-adapter-filtering)
    - [Decontamination](#decontamination)
    - [Alignment](#mappingalignment)
    - [Sort and Index alignment map](#sort-and-index-alignment-map)
    - [Primer trimming](#primer-trimming)
    - [Compute coverage](#compute-coverage)
    - [Variant calling](#variant-calling)
    - [Variant annotation](#variant-annotation)
    - [Consensus genome assembly](#consensus-genome-assembly)
    - [Pangolin: Lineage assignment](#pangolin-lineage-assignment)
    - [Nextclade: Clade assignment](#nexclade-clade-assignment)
- [Summarize results](#summarize-results)
- [Download reports](#download-reports)
- [Data Retrieval and Review](#data-retrieval-and-review)
    - [Transfer of data: MiSeq](#transfer-of-data-miseq)
    - [Transfer of data: MinION](#transfer-of-data-minion)
    - [Reviewing data: Illumina](#reviewing-data-illumina)
    - [Reviewing data: ONT](#reviewing-data-ont)
- [Working with metadata](#working-with-metadata)


## Introduction
In early January 2020, the novel coronavirus (SARS-CoV-2) responsible for a
pneumonia outbreak in Wuhan, China, was identified using next-generation
sequencing (NGS) and readily available bioinformatics pipelines. In addition to
virus discovery, these NGS technologies and bioinformatics resources are
currently being employed for ongoing genomic surveillance of SARS-CoV-2
worldwide, tracking its spread, evolution and patterns of variation on a global
scale.

## Scope
In this short workshop we will tackle, hands-on, the basic principles employed by the numerous bioinformatic pipelines:
to generate consensus genome sequences of SARS-CoV-2 and identify variants using
an actual dataset generated in our facility.

> **Note**

> This is part of the initiative fronted by the [Africa
CDC](https://africacdc.org/) with generous support from the [Rockeffeler
foundation](https://www.rockefellerfoundation.org/) to build capacity in pathogen genomics
in Africa.


## Background

We will use a dataset comprising of raw sequence reads of SARS-CoV-2 samples
obtained from a sequencing run on NextSeq 550 platorm at [ILRI](www.ilri.org).
NextSeq 550 flowcell uses 4 lanes; and so, 4 reads of data per sequenced sample corresponding to the 4
lanes are generated with suffixes L001, L002, L003 and L004. The dataset we are
using in this tutorial comprises of already concatenated sequences. These reads
can be combined/concatenated into a single file bearing in mind the type of library
sequencing either ```single``` or ```paired-end```.

## Prerequisite

This module will come after the introductory Linux module and therefore assumes familiarity with basic Linux command-line use. It also assumes you have an account and are operating in the ILRI computing cluster from a local Linux environment.

>**Note**

>Once inside the `hpc`, all instances of ```$USER``` will be equivalent to the hpc username that you were assigned, for example `Bio4Info$$`. Your username, by default, is stored in a variable called `USER`. By using it, you will not have to type-in your username, rather, your shell will automatically pick your username which is the value stored in the `USER` variable. The `$` (dollar) character-prefix to a variable name is used to call the value of that variable.

### Set-Up
We will use the computer lab at ILRI, which is already equipped with Linux-operating desktop computers. Since we will be working from the remote servers, we will not need special setup for personal laptops. However, toward the end of the program, we can look into access to a Linux server from a Windows PC; or how to install a Linux (sub)system for any interested persons.

### Preparation

#### ***Log into the HPC***
From the terminal (or equvalent tool) of your local computer, you can log into the HPC using the folowing command line, followed by pressing <ENTER>. You will be promted to type-in your password (the password will not be visible as you type it; just have faith). On a Linux system, you can use the `Ctrl-Alt-T` 7keyboard shortcut to open a terminal.
```
ssh <user_name>@hpc.ilri.cgiar.org
```
The HPC head node has 4 CPUs and we need to access more CPUs/resources in other compute nodes.
You will have to move from the cluster's master node into the node where we will be working from (it is called `compute05`). Use the following command; `-w` requests (a) specific list of host(s).
```
interactive -w compute05
```

`ssh` allows you to securely connect to the remote computer over internet, while `interactive` allows you to reserve resources to work interactively in a specified node within the computing cluster using the `-w` flag.
>**Note**
>When running a job interactively, the time limit is 8 hours and Default number of CPU is 1.

#### ***Project organisation***

1. We then change into the `compute05` `scratch` directory to create our project directory. Using the option`-p` (parent) `mkdir` will create any missing intermediate directories.
    ```
    cd /var/scratch/
    mkdir -p $USER/AfricaCDC_training
    cd $USER/AfricaCDC_training
    ```
2. The `assets, databases, primer-schemes, scripts` directories will be linked to the project directory, to limit redundancy. `-s` (soft) means that we are creating a soft link.
    ```
    ln -s /var/scratch/global/AfricaCDC_training/[adps]* .
    ```
3. We will create the directories `data`, `results` and `genome` to store raw data in ```fastq``` format, output and reference genomes respectively. Intermediate output files per `tool/software` will be created within the `results` directory. We will exploit the bash array data structure to create all the directories at once.
    ```
    mkdir data genome results
    mkdir -p results/{fastqc,fastp,kraken,samtools,ivar,snpeff,pangolin,nextclade,multiqc,bowtie2,bedtools}
    ```
4. Change into the `data` directory, from where we will retrieve our ```fastq``` files.
    ```
    cd data
    ls
    ```
#### ***Data retrieval and integrity checks***
1. While there are specialised tools for data retieval from nucleotide sequence databases, universal `Unix` command (`wget`) can be used to download data over internet.
    ```
    wget --no-check-certificate https://hpc.ilri.cgiar.org/~douso/AfricaCDC_training/fastq-metadata.tar.gz
    ```
2. After downloading your data, say from a sequencing facility site, it is often good practice to verify that your data was not intentionally/accidentally tampered with. To do this, your data service provider will likely accompany your data with a file containing a verification code: `checksum_file` (***will be provided***). The `md5sum` command, using the `-c` (check) tag, allows for checking the integrity of a file downloaded or acquired from a different source.
    ```
    wget --no-check-certificate https://hpc.ilri.cgiar.org/~douso/AfricaCDC_training/fastq-metadata.md5
    ls
    md5sum -c fastq-metadata.md5
    ```
3. Next, we will unzip the file using `tar` with the `-xf` (extract, file; respectively) tags, which tells `tar` extract the given file.
    ```
    tar -xf fastq-metadata.tar.gz
    ls
    ```
4.  Download SARS-CoV-2 reference genome and the genome annotation file.

    We will retrieve SARS-CoV-2 reference genome and the annotation from [NCBI](https://www.ncbi.nlm.nih.gov/).
    1. On a web browser, open the link [NCBI](https://www.ncbi.nlm.nih.gov/).
    2. Type 'SARS-CoV-2' on the search box and select 'Genome' database.
    3. Select the [Genbank](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/) hyperlink.
    4. Select the genome version [GCA_009858895.3_ASM985889v3](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/).
    5. Right click on the genome [FASTA](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.fna.gz) and select 'copy link'.
    6. Change into the ```genome``` directory using the command
    ```cd ../genome```.
    7. Use ```wget``` to fetch the file.
    8. Retrieve the feature annotation file [GFF](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gff.gz) using ```wget``` command.
    9. Dowload the [md5checksum](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/md5checksums.txt)  using `wget` command and check for integrity of your reference genome (FASTA) and annotation (GFF) files.

        ```
        echo "$(grep *GCA_009858895.3_ASM985889v3_genomic.fna.gz* md5checksums.txt | cut -f1 -d' ') GCA_009858895.3_ASM985889v3_genomic.fna.gz" | md5sum -c -
        echo "$(grep *GCA_009858895.3_ASM985889v3_genomic.gff.gz* md5checksums.txt | cut -f1 -d' ') GCA_009858895.3_ASM985889v3_genomic.gff.gz" | md5sum -c -
        ```

    11. If integrity check of the files has passed (`OK`), Uncompress the ```.gz``` files
       ```
       gunzip *.gz
       ```
    12. Rename the `FASTA` and `GFF` files
        ```
        mv GCA_009858895.3_ASM985889v3_genomic.fna nCoV-2019.fasta
        mv GCA_009858895.3_ASM985889v3_genomic.gff nCoV-2019.gff
        ```

## Analysis

#### ***Loading modules***
1. Clear the environment.
    ```
    module purge
    ```
2. Load modules using the `module load <tool-name>`command.
    ```
    module load fastqc/0.11.7
    module load fastp/0.22.0
    module load kraken/2.0.8-beta
    module load bowtie2/2.3.4.1
    module load samtools/1.11
    module load ivar/1.3.1
    module load bedtools/2.29.0
    module load R/3.6
    module load bcftools/1.11
    module load snpeff/4.1g
    module load multiqc/1.12
    module load nextclade/1.11.0
    module load python/3.9
    ```

    **Optional**
    The above modules can also be loaded using a single command
    ```
    module load fastqc/0.11.7 fastp/0.22.0 \
    kraken/2.0.8-beta bowtie2/2.3.4.1 samtools/1.11 ivar/1.3.1 \
    bedtools/2.29.0 R/3.6 bcftools/1.11 snpeff/4.1g multiqc/1.12 \
    nextclade/1.11.0 python/3.9
    ```
3. To list the loaded modules, type the below command.
    ```
    module list
    ```

#### ***Prepare the reference genome***


1. While still in the `genome` directory, we will index the reference sequence using samtools' `faidx`. Indexing produces a `.fai` file consisting of five tab-separated columns: `chrname, seqlength, first-base offset, seqlinewidth` without `\n` (newline character) and `seqlinewidth` with`\n`. This is essential for samtools' operations.

    ```
    samtools faidx nCoV-2019.fasta
    ```
    The above command generates the index for reference genome with the name `nCoV-2019.fasta.fai`.
2. We can take a sneak-view of the generated file and manipulate it for fun, say, to extract the genome size of reference fasta. This can be extracted from the `faidx`-indexed genome file using the ```cut``` command. The ```-f``` specifies the field(s) of interest.
    ```
    cut -f 1,2 nCoV-2019.fasta.fai > nCoV-2019.fasta.sizes
    ```
3. In order to allow easy access of genome regions during read mapping we will index the reference genome using ```bowtie2-build``` command.

    ```
    mkdir /var/scratch/$USER/AfricaCDC_training/genome/bowtie2
    ```

    ```
    bowtie2-build \
      --threads 1 \
      /var/scratch/$USER/AfricaCDC_training/genome/nCoV-2019.fasta \
      /var/scratch/$USER/AfricaCDC_training/genome/bowtie2/nCoV-2019
    ```
    The above command generates index files with the suffix `.bt2` for the reference genome with the prefix `nCoV-2019.`
4. Build SnpEff database for the reference genome

    [SnpEff](http://pcingola.github.io/SnpEff/se_introduction/), a variant annotation and predictor needs a database to perform genomic annotations. There are pre-built databases for thousands of genomes, so chances are that your organism of choice already has a SnpEff database available.

    >**Note** We will use pre-built SARS-CoV-2 SnpEff database


    ***Optional***
    In the (unlikely?) event that you need to build one yourself, you can build one using the commands found [here](http://pcingola.github.io/SnpEff/se_buildingdb/)

#### ***Quality assessment***
[`FastQC`](https://www.youtube.com/watch?v=bz93ReOv87Y)  is a common tool for Illumina read quality checks. The basic statistics from this report include `total sequences`, `sequence length` and `%GC`. Another 10 measures of quality are also graphically represented. Your experimental design will be crirical in interpreting `FastQC` reports. This step is very important for the subsequent data processes, especially at initial optimisation steps.


1. Change into the results ```fastqc``` directory
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/fastqc/
    ```
2. Run ```fastqc```
    ```
    fastqc \
        -t 1 \
        -o . \
        /var/scratch/$USER/AfricaCDC_training/data/COVM02379_R1.fastq.gz \
        /var/scratch/$USER/AfricaCDC_training/data/COVM02379_R2.fastq.gz
    ```
    ***Optional***
        Run step 3. above for the other 2 samples.

#### ***Quality and adapter filtering***
The preceeding step will guide us on the possible filtering and trimming operations to subject our data to. Depending on your study design, it is important to minimise noise as much as to zero, if possible. However, the latter case may be practically impossible.


1. Change into the output ```fastp``` directory.
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/fastp/
    ```

2. Run ```fastp```. `i,I` (input(s)) are for read1, read2; respectively. `o,O` (output(s)) are for the respective read1, read2; respectively. The `2>` construct redirects the standard error channel for saving as a log file.

    ```
    fastp \
        -w 1 \
        -i /var/scratch/$USER/AfricaCDC_training/data/COVM02379_R1.fastq.gz \
        -I /var/scratch/$USER/AfricaCDC_training/data/COVM02379_R2.fastq.gz \
        -o COVM02379_R1.trim.fastq.gz \
        -O COVM02379_R2.trim.fastq.gz \
        -h COVM02379.fastp.html \
        -j COVM02379.fastp.json \
        2> COVM02379.fastp.log
    ```

    ***Optional***
        Run steps 3 and 4 above for the other 2 samples.

#### ***Decontamination***

At times, sequencing experients will pick up non-target nucleic acids: for instance, host genetic material in SARS-CoV-2  sequencing. Such may obscure our signal of interest in the data; therefore, it is important to minimise or remove such sources of noise (unwanted background).
There are several bioinformatics tools and databases which can be used in querying the reads data in order to remove such noise. A commonly used tool is [Kraken2](https://github.com/DerrickWood/kraken2/wiki).
Kraken2 is a fast and memory efficient tool for taxonomic assignment of metagenomics sequencing reads. ```Kraken2``` can allow us to query the composition of our samples by searching for sequence reads against a pre-formatted database ("contaminant").

>**Note**
>Preformatted Kraken 2 and Bracken indexes can be found here: https://benlangmead.github.io/aws-indexes/k2 and downloaded without need of building new ones from scractch.


In this tutorial, we will use pre-formatted kraken2 ```human``` database to identify human-derived reads in our samples. This may give us an indication of contamination from host reads.

**Quiz:** *What type of contaminants would you think of in a SARS-CoV-2 sequencing experiment?*

---
<details close>
  <summary>Answer</summary>
  Host DNA,
  Host RNA and
  Internal control (PhiX)
</details>

---

1. Human database search
    - Change into the `kraken` directory results
        ```
        cd /var/scratch/$USER/AfricaCDC_training/results/kraken
        ```
    - Run `kraken2`
        ```
        kraken2 \
              -db /var/scratch/$USER/AfricaCDC_training/databases/kraken2-human-db \
              --threads 1 \
              --unclassified-out COVM02379.unclassified#.fastq \
              --classified-out COVM02379.classified#.fastq \
              --report COVM02379.kraken2.report.txt \
              --output COVM02379.kraken2.out \
              --gzip-compressed \
              --report-zero-counts \
              --paired /var/scratch/$USER/AfricaCDC_training/results/fastp/COVM02379_R1.trim.fastq.gz \
              /var/scratch/$USER/AfricaCDC_training/results/fastp/COVM02379_R2.trim.fastq.gz
        ```
    **Quiz:** *How many reads have hit the human genome as targets in the sample(s)?*

    ---
    <details close>
      <summary>Answer</summary>
      COVM02379  -  60529 reads    (13.00%)<br>
    </details>

    ---

    **Quiz:** *What percent of the sequencing reads are classified as SARS-CoV-2?*

    More information on output formats can be found [here](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats).

    ***Optional***
        Run steps 1 above for the other 2 samples.


    **Quiz:** *Based on the results, which sample do you think will give us good results in terms of genome coverage?*


#### ***Alignment***
Aligning sequence reads to a reference genome is the first step in many
comparative genomics pipelines, including pipelines for variant calling,
isoform quantitation and differential gene expression. In many cases, the alignment step is the slowest. This is because for each read the aligner must solve a difficult computational problem: determining the read's likely point of origin with respect to a reference genome. This is always non trivial for several reasons:
- The reference genome is often very big. Searching big things is harder than
  searching small things.
- You aren’t always looking for exact matches in the reference genome–or, at
least, probably not.

Here we use [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.


1. Change to the ```bowtie2``` directory.
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/bowtie2
    ```

2. Run the ```bowtie2``` command to align reads to the reference genome.

    ```
    bowtie2 \
          -x /var/scratch/$USER/AfricaCDC_training/genome/bowtie2/nCoV-2019 \
          -1 /var/scratch/$USER/AfricaCDC_training/results/kraken/COVM02379.unclassified_1.fastq \
          -2 /var/scratch/$USER/AfricaCDC_training/results/kraken/COVM02379.unclassified_2.fastq \
          --threads 1 \
          --un-conc-gz COVM02379.unmapped.fastq.gz \
          --local \
          --very-sensitive-local \
          2> COVM02379.bowtie2.log \
          | samtools view -@ 1 -F4 -bhS -o COVM02379.trim.dec.bam -
    ```
    ***Optional***
        Run steps 1 and 2. above for the other 2 samples.

#### ***Sort and Index alignment map***

Alignments can often be manipulated using [```samtools```](http://www.htslib.org/) using several sub commands


1. Sort the converted binary alignment (```.bam```)
    ```
    samtools sort -@ 1 -o COVM02379.sorted.bam -T COVM02379 COVM02379.trim.dec.bam
    ```

2. Index the sorted alignment
    ```
    samtools index -@ 1 COVM02379.sorted.bam
    ```
    ***Optional***
        Run steps 1 and 2. above for the other 2 samples.

#### ***Primer trimming***
Trim amplicon primers using [```iVar```](https://andersen-lab.github.io/ivar/html/manualpage.html)
iVar uses primer positions supplied in a BED file to soft clip primer sequences from an aligned and sorted BAM file. Following this, the reads are trimmed based on a quality threshold (Default: 20). To do the quality trimming, iVar uses a sliding window approach (Default: 4). The windows slides from the 5' end to the 3' end and if at any point the average base quality in the window falls below the threshold, the remaining read is soft clipped. If after trimming, the length of the read is greater than the minimum length specified (Default: 30), the read is written to the new trimmed BAM file.

1. Change to the output directory ```ivar```
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/ivar/
    ```

2. Run the command to trim primers

    ```
    ivar trim \
        -i /var/scratch/$USER/AfricaCDC_training/results/bowtie2/COVM02379.sorted.bam \
        -b /var/scratch/$USER/AfricaCDC_training/primer-schemes/V3/nCoV-2019.primer.bed \
        -p COVM02379.primertrimmed \
        -m 30 \
        -q 20 > COVM02379.ivar.log
    ```
3. Sort the primer trimmed alignment
    ```
    samtools sort \
          -@ 1 \
          -o COVM02379.primertrimmed.sorted.bam \
          -T COVM02379 COVM02379.primertrimmed.bam
    ```
4. Index the sorted primer trimmed alignment
    ```
    samtools index -@ 1 COVM02379.primertrimmed.sorted.bam
    ```


#### ***Compute coverage***
Here we will use [bedtools](https://github.com/arq5x/bedtools2), your swiss-army knife for genomic arithmetic and interval manipulation.

1. Change to the output directory ```bedtools```
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/bedtools/
    ```

2. Compute coverage
    ```
    bedtools \
        genomecov \
        -d \
        -ibam \
        /var/scratch/$USER/AfricaCDC_training/results/ivar/COVM02379.primertrimmed.sorted.bam \
        > COVM02379.coverage
    ```
3. Plot to visualize

    ```
    Rscript /var/scratch/$USER/AfricaCDC_training/scripts/plotGenomecov.R COVM02379.coverage
    ```



#### ***Variant calling***
`iVar` uses the output of the ```samtools mpileup``` command to call variants - single nucleotide variants(SNVs) and indels.


Pileup format consists of TAB-separated lines, with each line representing the pileup of reads at a single genomic position.

Several columns contain numeric quality values encoded as individual ASCII characters. Each character can range from "!" to "~" and is decoded by taking its ASCII value and subtracting 33; e.g., "A" encodes the numeric value 32.

The first three columns give the position and reference:

1. Chromosome name.
2. 1-based position on the chromosome.
3. Reference base at this position (this will be "N" on all lines if ```-f``` or ```--fasta-ref``` has not been used)

In generating the mpileup, we will use the flags:
```--count-orphans```: Do not skip anomalous read pairs in variant calling. Anomalous read pairs are those marked in the FLAG field as paired in sequencing but without the properly-paired flag set.
```--ignore-overlaps```: Disable read-pair overlap detection
```--no-BAQ```: Disable base alignment quality (BAQ) computation

The ```tee``` command, used with a pipe, reads standard input from ```samtools mpileup```, then writes the output of the program to standard output and simultaneously copies it into the specified file ```.mpileup```


In order to call variants correctly, the reference file used for alignment must be passed to `iVar` using the ```-r``` flag. The output of samtools pileup is piped into ivar variants to generate a ```.tsv``` file with the variants. There are two parameters that can be set for variant calling using `iVar` - minimum quality (Default: 20) and minimum frequency (Default: 0.03). Minimum quality is the minimum quality for a base to be counted towards the ungapped depth to calculate iSNV frequency at a given position. For insertions, the quality metric is discarded and the mpileup depth is used directly. Minimum frequency is the minimum frequency required for a SNV or indel to be reported.

`iVar` can identify codons and translate variants into amino acids using a [GFF](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) format containing the required coding regions (CDS). In absence of a GFF file, iVar will not perform the translation and "NA" will be added to the output file in place of the reference and alternate codons and amino acids.

1. Change to the output directory ```ivar```
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/ivar/
    ```

2. Call variants

    ```
    samtools mpileup \
            --ignore-overlaps \
            --count-orphans \
            --no-BAQ \
            --max-depth 0 \
            --min-BQ 0 \
            --reference /var/scratch/$USER/AfricaCDC_training/genome/nCoV-2019.fasta \
            /var/scratch/$USER/AfricaCDC_training/results/ivar/COVM02379.primertrimmed.sorted.bam \
            | tee COVM02379.mpileup \
            | ivar \
                variants \
                -t 0.25 \
                -q 20 \
                -m 10 \
                -g /var/scratch/$USER/AfricaCDC_training/genome/nCoV-2019.gff \
                -r /var/scratch/$USER/AfricaCDC_training/genome/nCoV-2019.fasta \
                -p COVM02379.variants
    ```
3. Convert the variants from ```.tsv``` to ```.vcf``` (Variant Call Format)

    ```
    python /var/scratch/$USER/AfricaCDC_training/scripts/ivar_variants_to_vcf.py \
      COVM02379.variants.tsv \
      COVM02379.vcf \
      --pass_only \
      --allele_freq_thresh 0.75 > COVM02379.variant.counts.log
    ```
    ##### **VCF file format**

    The header begins the file and provides metadata describing the body     of the file.
    Header lines are denoted as starting with `#`.
    Special keywords in the header are denoted with `##`.
    Recommended keywords   include fileformat, fileDate and reference.

    The header contains keywords that optionally semantically and      syntactically describe the fields used in the body of the file, notably INFO, FILTER, and FORMAT.


    |   |      Name    |  Brief description (see the specification for details)[VCF](https://samtools.github.io/hts-specs/VCFv4.1.pdf).  |
    |---|:-------------|:---------------------------------------------------------|
    | 1 |  CHROM       |The name of the sequence (typically a chromosome) on which the variation is being called.                                                           |
    | 2 |  POS         |The 1-based position of the variation on the given sequence.                                                          |
    | 3 |  ID          |The identifier of the variation, e.g. a dbSNP rs identifier, or if unknown a ".". Multiple identifiers should be separated by semi-colons without white-space.                                                           |
    | 4 |  REF         |The reference base (or bases in the case of an indel) at the given position on the given reference sequence.                                                          |
    | 5 |  ALT         |The list of alternative alleles at this position.                                                           |
    | 6 |  QUAL        |A quality score associated with the inference of the given alleles.                                                          |
    | 7 |  FILTER      |A flag indicating which of a given set of filters the variation has failed or PASS if all the filters were passed successfully.                                                          |
    | 8 |  INFO        |An extensible list of key-value pairs (fields) describing the variation.                                                          |
    | 9 |  FORMAT      |An (optional) extensible list of fields for describing the samples                                                          |
    | + |  SAMPLES     |For each (optional) sample described in the file, values are given for the fields listed in FORMAT                                                           |

4. Compress vcf file
    ```
    bgzip -c COVM02379.vcf > COVM02379.vcf.gz
    ```

5. Create tabix index from a sorted bgzip tab-delimited genome file
    ```
    tabix -p vcf -f COVM02379.vcf.gz
    ```
6. Generate stats from VCF file
    ```
    bcftools stats COVM02379.vcf.gz > COVM02379.stats.txt
    ```


#### ***Variant annotation***

We will use [SnpEff](http://pcingola.github.io/SnpEff/se_introduction/). It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes). It requires a configured SnpEff database with the annotation or features of the genome.

1. Change to the output directory ```snpeff```
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/snpeff/
    ```

2. Annotate and predict variants

    ```
    java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar \
        nCoV-2019 \
        -c /var/scratch/$USER/AfricaCDC_training/databases/snpeff_db/snpeff.config \
        -dataDir /var/scratch/$USER/AfricaCDC_training/databases/snpeff_db/data \
        /var/scratch/$USER/AfricaCDC_training/results/ivar/COVM02379.vcf.gz \
        > COVM02379.ann.vcf
    ```
3. Compress vcf file
    ```
    bgzip -c COVM02379.ann.vcf > COVM02379.ann.vcf.gz
    ```
4. Rename the ```summary.html``` and ```genes.txt``` file
    ```
    mv snpEff_summary.html COVM02379.summary.html
    mv snpEff_genes.txt COVM02379.genes.txt
    ```

5. Create tabix index from a sorted bgzip tab-delimited genome file
    ```
    tabix -p vcf -f COVM02379.ann.vcf.gz
    ```

6. Generate stats from VCF file
    ```
    bcftools stats COVM02379.ann.vcf.gz > COVM02379.stats.txt
    ```

7. Filter variants
    [SnpSift](http://pcingola.github.io/SnpEff/ss_introduction/) annotates genomic variants using databases, filters, and manipulates genomic annotated variants. Once you annotated your files using SnpEff, you can use SnpSift to help you filter large genomic datasets in order to find the most significant variants for your experiment.
    ```
      java -Xmx4g -jar /export/apps/snpeff/4.1g/SnpSift.jar \
            extractFields \
            -s "," \
            -e "." \
            COVM02379.ann.vcf.gz \
            CHROM POS REF ALT \
            "ANN[*].GENE" "ANN[*].GENEID" \
            "ANN[*].IMPACT" "ANN[*].EFFECT" \
            "ANN[*].FEATURE" "ANN[*].FEATUREID" \
            "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \
            "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \
            "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \
            "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \
            "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \
            > COVM02379.snpsift.txt
    ```


#### ***Consensus genome assembly***
To generate a consensus sequence iVar uses the output of samtools mpileup command. The mpileup output must be piped into ivar consensus. There are five parameters that can be set:
- minimum quality ```-q``` (Default: 20).
- minimum frequency threshold ```-t``` (Default: 0).
- minimum depth to call a consensus ```-m``` (Default: 10).
- a flag ```-n``` to exclude nucleotides from regions with depth less than the minimum depth and a character to call in regions with coverage lower than the speicifed minimum depth (Default: 'N').

Minimum quality is the minimum quality of a base to be considered in calculations of variant frequencies at a given position. Minimum frequency threshold is the minimum frequency that a base must match to be called as the consensus base at a position. If one base is not enough to match a given frequency, then an ambigious nucleotide is called at that position. Minimum depth is the minimum required depth to call a consensus. If ```-k``` flag is set then these regions are not included in the consensus sequence. If ```-k``` is not set then by default, a 'N' is called in these regions. You can also specfy which character you want to add to the consensus to cover regions with depth less than the minimum depth. This can be done using ```-n``` option. It takes one of two values: ```-``` or ```N```.

1. Change to the output directory ```ivar```
    ```
    cd /var/scratch/$USER/AfricaCDC_training/results/ivar/
    ```

2. Generate pileup and consensus genome sequences

    ```
    samtools \
            mpileup \
            --reference /var/scratch/$USER/AfricaCDC_training/genome/nCoV-2019.fasta \
            --count-orphans \
            --no-BAQ \
            --max-depth 0 \
            --min-BQ 0 \
            -aa \
            /var/scratch/$USER/AfricaCDC_training/results/ivar/COVM02379.primertrimmed.sorted.bam \
            | tee COVM02379.mpileup \
            | ivar \
                consensus \
                -t 0.75 \
                -q 20 \
                -m 10 \
                -n N \
                -p COVM02379.cons
    ```
The ```tee``` command reads from the standard input and writes to both standard output and one or more files at the same time. ```tee``` is mostly used in combination with other commands through piping.

3. Rename the consensus genome header

    ```
    sed -i '/^>/s/Consensus_\(.*\)_threshold.*/\1/' COVM02379.cons.fa
    ```



#### ***Nextclade: Clade assignment***
[**Nextclade**](https://docs.nextstrain.org/projects/nextclade/en/stable/) is a tool within the [**Nextrain**](https://nextstrain.org/) collection that uses sequence differences for their assignment to [clades](https://clades.nextstrain.org/). It also reports suspect quality issues with such sequences. There are both [web-](https://clades.nextstrain.org/) and [command-line-interfaces](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli.html) for *nextclade*. To run it in the command-line, we need some reference files: genome, feature map, origin tree, primers and quality configurations. Luckily, for SARS-CoV-2, these can be easily retrieved using the same tool, otherwise, you will have to create/retrieve accordingly.

1. Get the reference dataset
    ```
    nextclade dataset get --name 'sars-cov-2' --reference 'MN908947' --output-dir /var/scratch/$USER/AfricaCDC_training/nextclade_db
    ```
2. Perform clade assignment
    ```
    nextclade \
       --input-fasta /var/scratch/$USER/AfricaCDC_training/results/ivar/COVM02379.cons.fa \
       --input-dataset /var/scratch/$USER/AfricaCDC_training/nextclade_db \
       --input-root-seq /var/scratch/$USER/AfricaCDC_training/nextclade_db/reference.fasta \
       --genes E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
       --input-gene-map /var/scratch/$USER/AfricaCDC_training/nextclade_db/genemap.gff \
       --input-tree /var/scratch/$USER/AfricaCDC_training/nextclade_db/tree.json \
       --input-qc-config /var/scratch/$USER/AfricaCDC_training/nextclade_db/qc.json \
       --input-pcr-primers /var/scratch/$USER/AfricaCDC_training/nextclade_db/primers.csv \
       --output-csv /var/scratch/$USER/AfricaCDC_training/results/nextclade/COVM02379.csv \
       --output-tree /var/scratch/$USER/AfricaCDC_training/results/nextclade/COVM02379.auspice.json \
       --output-dir /var/scratch/$USER/AfricaCDC_training/results/nextclade/ \
       --output-basename COVM02379.cons \
       2> /var/scratch/$USER/AfricaCDC_training/results/nextclade/COVM02379.nextclade.log
    ```

3. Visualization: The output of Nextclade includes a phylogenetic tree in `.json` format. This tree can be visualized in [Auspice](https://auspice.us/). First let us download the the `.json` file:
```
cp /var/scratch/$USER/AfricaCDC_training/results/nextclade/COVM02379.auspice.json ~/
```
In your local computer use `scp` to copy the file to any desired destination:
```
scp <user_name>@hpc.ilri.cgiar.org:/home/<user_name>/*.json <destination_folder>
```
Open [Auspice](https://auspice.us/) and drag and drop the `.json` file in the [Auspice](https://auspice.us/). Now edit the tree.
  - In `Dataset` click the drop down arrow and select `ncov`, below it select `open` and below it select `global`.
  - In `Color By` click the drop down arrow and select `clade`.
  - Do any other adjustments as you wish.


#### ***Pangolin Lineage Assignment***

[Phylogenetic Assignment of Named Global Outbreak Lineages (Pangolin)](https://cov-lineages.org/resources/pangolin.html) implements a dynamic nomenclature of SARS-CoV-2 lineages, known as the Pango nomenclature. To assign [Pangolin Lineages](https://cov-lineages.org/lineage_list.html), we will use [the web version of pangolin](https://pangolin.cog-uk.io/). It also has a robust [command-line version](https://github.com/cov-lineages/pangolin) that we will look into later. With the web version we must retrieve out consensus genome from the analysis server (HPC) to our local computer. Follow the following steps to assign your query sequences pangolin lineages. Also here is a [***tutorial***](https://cov-lineages.org/resources/pangolin/tutorial.html).

1. In your bowser open [Pangolin Web Application](https://pangolin.cog-uk.io/). This is the online version of [Pangolin](https://github.com/cov-lineages/pangolin).
2. Now copy-paste/drag-drop your consensus file to the site and click `Start Analysis`
3. Once done, download the results and if you need to, copy the names of sequences that failed the analysis to a file. The download is called `results.csv`.

Here is how pangolin performs the analysis:
![alt text](https://cov-lineages.org/assets/images/pangolin_pipelines.svg "Pangolin Analysis Workflow")

Alternatively to perform the commandline analysis for Pangolin, let us proceed as follows. We will need to use a singularity image (think of a singularity image as a ready-to-use container within which we have packaged all the software needed to do a certain task) in this case packaging pangolin softwares*.
1. Let us create a directory to store our image:
```
mkdir /var/scratch/$USER/AfricaCDC_training/singularity
```
2. Download the image:
```
singularity pull --dir /var/scratch/$USER/AfricaCDC_training/singularity/ \
                    --force docker://staphb/pangolin:latest
```
3. Download Pangolin's referemce data: Downloads/updates to latest release of pangoLEARN and constellations
```
singularity run /var/scratch/$USER/AfricaCDC_training/singularity/pangolin_latest.sif \
          pangolin --update-data \
          --datadir /var/scratch/$USER/AfricaCDC_training/pangolin_db
```
4. Conduct Pangolin Lineage assignment:
```
singularity run /var/scratch/$USER/AfricaCDC_training/singularity/pangolin_latest.sif
          pangolin /var/scratch/$USER/AfricaCDC_training/results/ivar/COVM02379.cons.fa \
          --alignment \
          --usher \
          --max-ambig 0.3 \
          --min-length 25000 \
          --outdir /var/scratch/$USER/AfricaCDC_training/results/pangolin/ \
          --outfile COVM02379.pangolin.usher.csv \
          --datadir /var/scratch/$USER/AfricaCDC_training/pangolin_db/ \
          2> /var/scratch/$USER/AfricaCDC_training/results/pangolin/COVM02379.pangolin.usher.log
```


##Running [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline.
The [nf-core](https://nf-co.re/) -"A community effort to collect a curated set of analysis pipelines built using Nextflow" - has [many pipelines that can be easily setup and used to analyse genomics data. We will be using [nf-core/viralrecon](https://nf-co.re/viralrecon) for our analysis. The most recent verion being [2.4.1](https://nf-co.re/viralrecon/2.4.1).

To set it up we will follow the workflow in [Launch pipeline](https://nf-co.re/launch?id=1649311271_b34e829e5e3f).

Our data is stored in `/var/scratch/global/ilri_AuCDC/miseq`. We need to:
- SSH into HPC:

```
ssh <username>@hpc.ilri.cgiar.org
```
- Go into ineractive mode in compute05
```
interactive -w compute05
```
- Symbolicly link our data to the to a directory in `scratch`
```
mkdir /var/scratch/user10/miseq_analysis/
cd /var/scratch/user10/miseq_analysis/
ln -s /var/scratch/global/ilri_AuCDC/miseq ./
```
- Now let us go back to [Launch pipeline](https://nf-co.re/launch?id=1649311271_b34e829e5e3f) and do step by step set up.

- Finally let us transfer `the parameters JSON to a file` to the HPC.

- We can now launch the analysis as follows:
  - Load modules and set some Java options
```
module load nextflow/21.10
NXF_OPTS='-Xms1g -Xmx4g'
```
  - Launch the analysis: Use your user name in <user##> and the you will see your run name in <nextflow run nf-core/viralrecon -r 2.4.1 -name AfricaCDC_sarscov2 -profile singularity -resume -params-file nf-params.json>:
```
srun --partition=batch -w compute05 -J <user##> -n 3 nextflow run nf-core/viralrecon -r 2.4.1 -name <AfricaCDC_sarscov2> -profile singularity -resume -params-file nf-params.json
```
- Wait for the run to be completed.
- Let us download the multiQC file and visualize as follows:
  - Check the `results` directory and `results/multiqc` directory:
  ```
  ls ./results
  ls ./results/multiqc
  ```
  - Download the `./results/multiqc/multiqc_report.html` and visualize
  - Clue: use the command `scp`


## Summarize results
Aggregate results from bioinformatics analyses across many samples into a single report with [MultiQC](https://multiqc.info/)

1. Change to the output directory `multiqc_out`
    ```
    $ cd /var/scratch/$USER/AfricaCDC_training/results/multiqc/
    ```
2. Aggregate tools' outputs
    ```
    $ multiqc \
        --force \
        --title SARS-CoV-2 \
        --export \
        --outdir . \
        --config /var/scratch/$USER/AfricaCDC_training/assets/multiqc_config.yaml \
        /var/scratch/$USER/AfricaCDC_training/results/
    ```

## Download reports

1. To help with file transfers, we will create a user-specific directory inside the `global` temporary directory of the head-node.
    ```
    $ mkdir /var/scratch/global/AfricaCDC_training_outputs/$USER
    ```
2. Copy output files of interest (multiqc, primer-trimmed bam).
    ```
    cp /var/scratch/$USER/AfricaCDC_training/results/ivar/*.sorted.bam* \
        /var/scratch/global/AfricaCDC_training_outputs/$USER/
    cp -r SARS-CoV-2_multiqc_report_plots *.html \
        /var/scratch/global/AfricaCDC_training_outputs/$USER/
    ```
3. On your local computer, open a terminal and create a directory in the `Downloads` directory

    ```
    mkdir -p /mnt/c/Downloads/AfricaCDC_training/results
    cd /mnt/c/Downloads/AfricaCDC_training/results/
    ```

4. Copy all the contents of the `/var/scratch/global/AfricaCDC_training_outputs/<username>/` directory in the HPC to the local outputs directory you created in the previous step.
    >**Note**

    >Replace ```username``` with the actual provided HPC account username
    ```
    rsync \
        -avP \
        <username>@hpc.ilri.cgiar.org:/var/scratch/global/AfricaCDC_training_outputs/<username>/* \
        .
    ```
    OR
    ```
        scp <username>@hpc.ilri.cgiar.org:/var/scratch/global/AfricaCDC_training_outputs/<username>/* \
        .
    ```

## Data Retrieval and Review
Having sequenced our samples in both MiSeq (Illumina) and MinION (ONT) we can now transfer the sequence output to the HPC where we will conduct the bioinformatics analysis.

#### ***Transfer of data: Illumina***
MiSeq is based on Windows and data transfer will be done by copy-pasting the data to HPC through the network. Ensure the transfer is completed successfully without errors.

#### ***Transfer of data: MinION***
The MinION sequencer stores its sequencing output in a Linux based computer. To transfer the data we logged into computer and transferred the data on the command line as follows.
```
rsync -avP <path-to-the-directory-with_sequencing-ouput>/ <username>:<path-to-the-directory-to-store-sequencing-ouput>/
```
Replaced `<path-to-the-directory-with_sequencing-ouput>` with the path to the directory storing the sequencing output. Replaced `<HPC-login-username>` with your hpc login username (i.e user##@hpc.ilri.cgiar.org) and `<path-to-the-directory-to-store-sequencing-output>` with path to the directory you want to store the data in the HPC.
Example:`
rsync -avP /media/SeqData_LTS/20220405_1121_MN2816_FAH91436_2b8e9827/ <username>@hpc.ilri.cgiar.org:/var/scratch/global/$USER/20220405_1121_MN2816_FAH91436_2b8e9827`

#### ***Reviewing data: Illumina***
Change working directory into the directory that stores the Illumina dataset
```
cd /var/scratch/global/miseq/
```
Now let us view what is the sequencing output. The output is a FASTA format. Note: in the second command replace `<one-of-the-fastq.gz>`
 with the name of the fastq.gz file
 ```
ls fastq/
less -S fastq/<one-of-the-fastq.gz>
```

#### ***Reviewing data: ONT***
Change working directory into directory that stores the ONT data. Note: replace `<name-of-run-folder>` with the name of the run folder
```
cd /var/scratch/global/ONT/
ls <name-of-run-folder>
```

## Working with metadata
Metadata is the data associated with your main data; it describes your data in a manner that can allow drawing information upon analysing the data. Often, the importance of `metadata` is ignored; things like how to capture, store, encode and organise metadata, especially, having the downstream analyses and interpretation in mind.

> **IMPORTANT**
> Data without metadata is, mostly, garbage.

We will take a look at an example metadata to highlight some concerns with metadata, and the reason(s) why they are important in data analyses workflows.

test_lab|case_id|lab_id|loc|age|gender|occup|samp_type|symp|vacc_state|coll_dt|confir_dt|recep_dt
|---|---|---|---|---|---|---|---|---|---|---|---|---|
<strong>TESTING_LAB</strong>|<strong>CASE-ID</strong>|<strong>SampleNumber</strong>|<strong>LOCATION</strong>|<strong>AGE</strong>|<strong>GENDER (M/F)</strong>|<strong>OCCUPATION</strong>|<strong>SAMPLE TYPE</strong>|<strong>SYMPTOMS SHOWN (COUGH;FEVER;ETC)</strong>|<strong>VACCINATION STATUS</strong>|<strong>DATE OF SAMPLE COLLECTION</strong>|<strong>LAB CONFIRMATION DATE</strong>|<strong>DATE SAMPLE RECEIVED IN THE LAB</strong>
LABA|place/id/date|COVD0308|Some Place|80|F|Food handler|NP & OP Swab|Asymptomatic|Yes|5th June 2020|08-Jun-20|08/06/20
LABB|COM/SARS001/2022|COVD360|Comoros|38|F|None|NP-OP Swab|FC;CO;H|Yes|5th June 2020|08-Jun-20|08/06/20
LABC|DJI/SARS001/2023|COVD 273|Djibouti|30|M|Business|NP-OP Swab|Asymptomatic|No|5th June 2020|08-Jun-20|08/06/20
LABD|SWZ/SARS001/2024|COVD154|Eswatini|23|Female|Food seller|OP-NP Swap|No symptoms|Yes|5th June 2020|08-Jun-20|08/06/20
LABE|ETH/SARS001/2025|COVD0875|Ethiopia|34|M|Targeted testing|NP-OP Swab|Asymptomatic|Yes|5th June 2020|08-Jun-20|08/06/20
LABF|LBY/SARS001/2027|COVD00672|Libya|18|Female|Food handler|NP-OP Swab|Fever/Chills, Cough, Headache|Yes|5th June 2020|08-Jun-20|08/06/20
LABG|MDG/SARS001/2028|COVD499|Madagascar|25|F|Targeted testing|NP Swab|Asymptomatic|No|5th June 2020|08-Jun-20|08/06/20
LABH|MUS/SARS001/2029|COVD078|Mauritius|25|F|Not indicated|NP-OP Swab|Asymptomatic|No|5th June 2020|08-Jun-20|08/06/20
LABI|SYC/SARS001/2030|COVD579|Seychelles|2 years|Male|Targeted testing|NP Swab|Asymptomatic|No|5th June 2020|08-Jun-20|08/06/20
LABJ|SOM/SARS001/2031|300|Somalia|26|F|Food handler|NP-OP Swab|Asymptomatic|No|5th June 2020|08-Jun-20|08/06/20
LABK|SSD/SARS001/2031|COVD00381|South Sudan|45|M|NA|NP-OP Swab|Asymptomatic|Yes|5th June 2020|08-Jun-20|08/06/20

What are some of the issues you notice with the above metadata?

---
<details close>
  <summary>Answer</summary>

  - Inconsistent header naming
  - Future dates
  - Mixed data types within column
  - Inconsistent date formats
  - Inconsistent sample type capturing
  - Inconsistent symptoms
  - Spaces in headers
  - Use of commas
  - ...
</details>


---

***Computation-wise, whichever way (poor/good) you choose to organise your data, ensure consistency***.
##Galaxy Tutorial
https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/sars-cov-2-variant-discovery/tutorial.html

---
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
