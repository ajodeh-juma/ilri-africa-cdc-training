---
title: README.md
---

# **Building capacity in SARS-CoV-2 genomics in Africa**

- [Introduction](#introduction)
- [Background](#background)
- [Genome preparation](#genome-preparation)
- [Data retrieval and integrity checks](#data-retrieval-and-integrity-checks)
- [Bioinformatic Analysis](#bioinformatic-analysis)
- [Quality assessment](#quality-assessment)
- [Decontamination](#decontamination)
- [Mapping/Alignment](#mappingalignment)
- [Sort and Index alignment](#sort-and-index-alignment)
- [Primer trimming](#primer-trimming)
- [Compute coverage](#compute-coverage)
- [Variant calling](#variant-calling)
- [Variant annotation](#variant-annotation)
- [Consensus genome](#consensus-genome)
- [Summarize results](#summarize-results)
- [Download MultiQC](#download-multiqc)
- [Pangolin Lineage assignment](#pangolin-lineage-assignment)
- [Nexclade analysis](#nexclade-analysis)
- [Data Retrieval and Review](#data-retrieval-and-review)



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

>In this tutorial, replace all instances of ```<$user>``` with the hpc username that you were assigned, for example `student1`. Alternatively, you may store your username in a variable like `user=student1`. This way, you will not have to do the replacements, rather, your shell will automatically pick your username as the value stored in the `$user` variable.
>**Note**
>If you get disconnected or logout from the hpc, you WILL have to repeat the above preparatory step.

### Preparation

#### ***Log into the computing cluster***
From the terminal (or equvalent tool) of your local computer, you can log into the HPC using the folowing command line, followed by pressing <ENTER>. You will be promted to type-in your password:
```
$ ssh $user@hpc.ilri.cgiar.org
```
The HPC head node has 4 CPUs and we need to access more CPUs/resources in other compute nodes.
You will have to move from the cluster's master node into the node where we will be working from (it is called `compute05`). Use the following command:
```
$ interactive -w compute05
```

`ssh` allows you to securely connect to the remote computer over internet, while `interactive` allows you to reserve resources to work interactively in a specified node within the computing cluster using the `-w` flag.
>**Note**
>When running a job interactively, the time limit is 8 hours and Default number of cpus is 1

#### ***Project organisation***
1. Change into the global temporary directory and create a project directory and `cd` into it.
    ```
    $ cd /var/scratch/global/
    $ mkdir -p $user/AfricaCDC_training
    $ cd $user/AfricaCDC_training
    ```
2. The `assets, databases, primer-schemes, scripts` directories will be linked to the project directory.
    ```
    $ ln -s /var/scratch/global/AfricaCDC_training/[adps]* .
    ```
3. We will create the directories `data`, `results` and `genome` to store raw data in ```fastq``` format, output and reference genomes respectively. Intermediate output files per `tool/software` will be created within the `results` directory
    ```
    $ mkdir data results genome
    $ mkdir -p results/{fastqc,fastp,kraken/{human,standard},samtools,bowtie2,ivar/{trim,variants,consensus,coverage},snpeff,multiqc}

    ```  
4. Change into the `data` directory, from where we will deposit our retrieved ```fastq``` files.
    ```
    $ cd data
    ```
#### ***Data retrieval and integrity checks***
1. While there are specialised tools for data retieval from NA databases, universal `Unix` command (`wget`) can be used to download data.
    ```
     $ wget https://hpc.ilri.cgiar.org/~douso/AfricaCDC_training/enc_fastq.tar.gz

    ```
2. After downloading your data, say from a sequencing facility site, it is often good practice to verify that your data was not intentionally/accidentally tampered with. To do this, your data service provider will likely accompany your data with a file containing a verification code: `checksum_file`. The `md5sum` command allows for checking the integrity of a file downloaded or acquired from a different source.
    ```
     $ md5sum -c <checksum_file.md5>
    ``` 
3.  Download SARS-CoV-2 genome and annotation
 
    We will retrieve SARS-CoV-2 reference genome and annotation from [NCBI](https://www.ncbi.nlm.nih.gov/).
    1. On a web browser, open the link [NCBI](https://www.ncbi.nlm.nih.gov/).
    2. Type 'sars-cov-2' on the search box and select 'Genome' database
    3. Select the [Genbank](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/) hyperlink.
    4. Select the genome version [GCA_009858895.3_ASM985889v3](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/)
    5. Right click on the genome [FASTA](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.fna.gz) and select 'copy link'
    6. Change into the ```genome``` directory using the command
    ```cd ../genome```
    7. Use ```wget``` to fetch the file.
    8. Retrieve the feature annotation file [GFF](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.gff.gz) using ```wget``` command
    9. Dowload the [md5checksum](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Severe_acute_respiratory_syndrome-related_coronavirus/latest_assembly_versions/GCA_009858895.3_ASM985889v3/md5checksums.txt)  using `wget` command and check for integrity of your files using the command:
        `md5sum GCA_009858895.3_ASM985889v3_genomic.fna.gz`
        `md5sum GCA_009858895.3_ASM985889v3_genomic.gff.gz`
    11. If integrity check of the files has passed (equal to values in the 'md5checksum.txt', Uncompress the ```.gz``` files
        `gunzip *.gz`
    13. Rename the `FASTA` and `GFF` files
        ```
        mv GCA_009858895.3_ASM985889v3_genomic.fna nCoV-2019.fasta
        mv GCA_009858895.3_ASM985889v3_genomic.gff nCoV-2019.gff
        ```

## Analysis

#### ***Making available the modules to be used***  
1. Clear the environment.
    ```
    $ module purge
    ```
2. Load modules using the `module load <tool-name>`command. 
    ```
    $ module load fastqc/0.11.7
    $ module load fastp/0.22.0
    $ module load kraken/2.0.8-beta
    $ module load bowtie2/2.3.4.1
    $ module load samtools/1.11
    $ module load ivar/1.3.1
    $ module load bedtools/2.29.0
    $ module load R/3.6
    $ module load bcftools/1.11
    $ module load snpeff/4.1g
    $ module load multiqc/1.6
    ```
3. To list the loaded modules, type the command
    ```
    module list
    ```

#### ***Prepare the reference genome***


1. While still in the `genome` directory, index reference sequence in the FASTA format or extract subsequence from indexed reference sequence. This is essential for the variants calling step.

    ```
    $ samtools faidx nCoV-2019.fasta
    ```  
    The above command generates the index for reference genome with the name `nCoV-2019.fasta.fai`
2. We will need the genome size for downstream step . This can be extracted from the `faidx`-indexed genome file using the ```cut``` command. The ```-f``` specifies the field(s) of interest. 
    ```
    $ cut -f 1,2 nCoV-2019.fasta.fai > nCoV-2019.fasta.sizes
    ```
3. In order to allow easy access of genome regions during read mapping we will index the reference genome using ```bowtie2-build``` command.  

    ```
    $ mkdir /var/scratch/global/$user/AfricaCDC_training/genome/bowtie2
    ```

    ```
    $ bowtie2-build \
      --threads 1 \
      /var/scratch/global/$user/AfricaCDC_training/genome/nCoV-2019.fasta \
      /var/scratch/global/$user/AfricaCDC_training/genome/bowtie2/nCoV-2019
    ```
    The above command generates index files with the suffix `.bt2` for the reference genome with the prefix `nCoV-2019.`
4. Build SnpEff database for the reference genome

    [SnpEff](http://pcingola.github.io/SnpEff/se_introduction/), a variant annotation and predictor needs a database to perform genomic annotations. There are pre-built databases for thousands of genomes, so chances are that your organism of choice already has a SnpEff database available. In the (unlikely?) event that you need to build one yourself, you can build one. We will build a database for the SARS-CoV-2 genome.
    
    ```
    mkdir -p snpeff_db/genomes/
    cd snpeff_db/genomes/
    ln -s /var/scratch/global/$user/AfricaCDC_training/genome/nCoV-2019.fasta nCoV-2019.fa
    ```
    
    ```
    cd ../../
    mkdir -p snpeff_db/nCoV-2019/
    cd snpeff_db/nCoV-2019/
    ln -s /var/scratch/global/$user/AfricaCDC_training/genome/nCoV-2019.gff genes.gff
    ```

    ```
    cd ../../
    echo "nCoV-2019.genome : nCoV-2019" > snpeff.config
    java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar \
        build \
        -config snpeff.config \
        -dataDir ./snpeff_db \
        -gff3 \
        -v \
        nCoV-2019
    ```

#### ***Quality assessment***
[`FastQC`](https://www.youtube.com/watch?v=bz93ReOv87Y)  is a common tool for Illumina read quality checks. The basic statistics from this report include `total sequences`, `sequence length` and `%GC`. Another 10 measures of quality are also graphically represented. Your experimental design will be crirical in interpreting `FastQC` reports. This step is very important for the subsequent data processes, especially at initial optimisation steps.


1. Change into the output ```fastqc_out``` directory
    ```
    $ cd /var/scratch/global/$user/AfricaCDC_training/results/fastqc/
    ```
2. Run ```fastqc```  
    ```
    $ fastqc \
        -t 1 \
        -o . \
        /var/scratch/global/$user/AfricaCDC_training/data/COVM02379_R1.fastq.gz \
        /var/scratch/global/$user/AfricaCDC_training/data/COVM02379_R2.fastq.gz
    ```
3. ***Optional***
        Run step 3. above for the other 2 samples.
        
#### ***Quality and adapter filtering***
The preceeding step will guide us on the possible filtering and trimming operations to subject our data to. Depending on your study design, it is important to minimise noise as much as to zero, if possible. However, the latter case may be practically impossible.


1. Change into the output ```fastp``` directory
    ```
    $ cd /var/scratch/global/$user/AfricaCDC_training/results/fastp/
    ```

2. Run ```fastp``` 
    
    ```
    $ fastp \
        -i /var/scratch/global/$user/AfricaCDC_training/data/COVM02379_R1.fastq.gz \
        -I /var/scratch/global/$user/AfricaCDC_training/data/COVM02379_R2.fastq.gz \
        -o COVM02379_R1.trim.fastq.gz \
        -O COVM02379_R2.trim.fastq.gz \
        2> COVM02379.log
    ```
3. Rename `.json` and `.html` files
    ```
    mv fastp.json COVM02379.fastp.json
    mv fastp.html COVM02379.fastp.html
    ```
4. ***Optional***
        Run steps 3 and 4 above for the other 2 samples.

#### ***Metagenomics***

At times, sequencing experients will pick up non-target nucleic acids: for instance, host genetic material in SARS-CoV-2  sequencing. Such may obscure our signal of interest in the data; therefore, it is important to minimise or remove such sources of noise (unwanted background).
There are several bioinformatics tools and databases which can be used in querying the reads data in order to remove such noise. A commonly used tool is [Kraken2](https://github.com/DerrickWood/kraken2/wiki). 
Kraken 2 is a fast and memory efficient tool for taxonomic assignment of metagenomics sequencing reads. ```Kraken2``` can allow us to query the composition of our samples by searching for sequence reads against a pre-formatted database (contaminant).

>**Note**
>Preformatted Kraken 2 and Bracken indexes can be found here: https://benlangmead.github.io/aws-indexes/k2 and downloaded without need of building new ones from scractch.


In this tutorial, we will use pre-formatted kraken2 ```human``` and ```standard``` databases to identify human-derived reads in our samples. This may give us an indication of contamination from host reads.

**Quiz:** *What type of contaminants would you think of in a SARS-CoV-2 sequencing experiment?*

---
<details close>
  <summary>Answer</summary>
  Host DNA
  Host RNA
  Internal control (PhiX)
</details>

---

1. Human database search
    - Change into the ```kraken/human``` directory
        ```
        $ cd /var/scratch/global/$user/AfricaCDC_training/results/kraken/human/
        ```
    - Run ```kraken2```
        ```
        $ kraken2 \
              -db /var/scratch/global/AfricaCDC_training/databases/kraken2-human-db \
              --threads 1 \
              --unclassified-out COVM02379.unclassified#.fastq \
              --classified-out COVM02379.classified#.fastq \
              --report COVM02379.kraken2.report.txt \
              --output COVM02379.kraken2.out \
              --gzip-compressed \
              --report-zero-counts \
              --paired /var/scratch/global/$user/AfricaCDC_training/data/COVM02379_R1.fastq.gz \
              /var/scratch/global/$user/AfricaCDC_training/data/COVM02379_R2.fastq.gz
        ```
    **Quiz:** *How many reads have hit the human genome as targets in the samples?*

    ---
    <details close>
      <summary>Answer</summary>
      COVM02379  -  66438 reads    (13.87%)<br>
      COVM02524  -  60686 reads    (57.23%)<br>
      COVM02720  -  245351 reads    (83.17%)
    </details>

    ---
        
1. Standard database search
    - Change into the ```kraken/standard``` directory
        ```
        $ cd /var/scratch/global/$user/AfricaCDC_training/results/kraken/standard/
        ```
    - Run ```kraken2```
        ```
        $ kraken2 \
                  -db /export/data/bio/kraken2/db/standard/ \
                  --threads 1 \
                  --unclassified-out COVM02379.unclassified#.fastq \
                  --classified-out COVM02379.classified#.fastq \
                  --report COVM02379.kraken2.report.txt \
                  --output COVM02379.kraken2.out \
                  --gzip-compressed \
                  --report-zero-counts \
                  --paired /var/scratch/global/$user/AfricaCDC_training/data/COVM02379_R1.fastq.gz \
                  /var/scratch/global/$user/AfricaCDC_training/data/COVM02379_R2.fastq.gz
        ```

    **Quiz:** *What percent of the sequencing reads are classified as SARS-CoV-2?*
    
    More information on output formats can be found at https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats

3. ***Optional***
        Run steps 1 and 2. above for the other 2 samples.


    **Quiz:** *Based on the results, which sample do you think will give us good results in terms of genome coverage?*


## Alignment
Aligning sequence reads to a reference genome is the first step in many
comparative genomics pipelines, including pipelines for variant calling,
isoform quantitation and differential gene expression. In many cases, the alignment step is the slowest. This is because for each read the aligner must solve a difficult computational problem: determining the read's likely point of origin with respect to a reference genome. This is always non trivial for several reasons:
- The reference genome is really big. Searching big things is harder than
  searching small things.
- You aren’t always looking for exact matches in the reference genome–or, at
least, probably not.

Here we use [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.


1. Change to the ```bowtie2``` directory
    ```
    $ cd /var/scratch/global/$user/AfricaCDC_training/results/bowtie2
    ```  

2. Run the ```bowtie2``` command to align reads to the reference genome

    ```
    $ bowtie2 \
          -x /var/scratch/global/$user/AfricaCDC_training/genome/bowtie2/nCoV-2019 \
          -1 /var/scratch/global/$user/AfricaCDC_training/results/fastp/COVM02379_R1.trim.fastq.gz \
          -2 /var/scratch/global/$user/AfricaCDC_training/results/fastp/COVM02379_R2.trim.fastq.gz \
          --threads 1 \
          --un-conc-gz COVM02379.unmapped.fastq.gz \
          --local \
          --very-sensitive-local \
          2> COVM02379.bowtie2.log \
          | samtools view -@ 1 -F4 -bhS -o COVM02379.bam -
    ```
3. ***Optional***
        Run steps 1 and 2. above for the other 2 samples.
    
## Sort and Index alignment

Alignments can often be manipulated using [```samtools```](http://www.htslib.org/) using several sub commands


1. Sort the converted binary alignment (```.bam```)
    ```
    $ samtools sort -@ 1 -o COVM02379.sorted.bam -T COVM02379 COVM02379.bam
    ```
    
2. Index the sorted alignment
    ```
    $ samtools index -@ 1 COVM02379.sorted.bam
    ```
3. ***Optional***
        Run steps 1 and 2. above for the other 2 samples.

## Primer trimming
Trim amplicon primers using [```iVar```](https://andersen-lab.github.io/ivar/html/manualpage.html)
iVar uses primer positions supplied in a BED file to soft clip primer sequences from an aligned and sorted BAM file. Following this, the reads are trimmed based on a quality threshold (Default: 20). To do the quality trimming, iVar uses a sliding window approach (Default: 4). The windows slides from the 5' end to the 3' end and if at any point the average base quality in the window falls below the threshold, the remaining read is soft clipped. If after trimming, the length of the read is greater than the minimum length specified (Default: 30), the read is written to the new trimmed BAM file.

1. Change to the output directory ```ivar/trim```
    ```
    $ cd /var/scratch/global/$user/AfricaCDC_training/results/ivar/trim/
    ```  

2. Run the command to trim primers

    ```
    $ ivar trim \
        -i /var/scratch/global/$user/AfricaCDC_training/results/bowtie2/COVM02379.sorted.bam \
        -b /var/scratch/global/AfricaCDC_training/primer-schemes/V3/nCoV-2019.primer.bed \
        -p COVM02379.primertrimmed \
        -m 30 \
        -q 20 > COVM02379.ivar.log
    ```
3. Sort the primer trimmed alignment
    ```
    $ samtools sort \
          -@ 1 \
          -o COVM02379.primertrimmed.sorted.bam \
          -T COVM02379 COVM02379.primertrimmed.bam
    ```
4. Index the sorted primer trimmed alignment
    ```
    samtools index -@ 1 COVM02379.primertrimmed.sorted.bam
    ```


## Compute coverage
Here we will use [bedtools](https://github.com/arq5x/bedtools2), your swiss-army knife for genomic arithmetic and interval manipulation.

1. Change to the output directory ```ivar/coverage```
    ```
    $ cd /var/scratch/global/$user/AfricaCDC_training/results/ivar/coverage/
    ```  

2. Compute coverage
    ```
    $ bedtools \
        genomecov \
        -d \
        -ibam \
        /var/scratch/global/$user/AfricaCDC_training/results/ivar/trim/COVM02379.primertrimmed.sorted.bam \
        > COVM02379.coverage
    ```
3. Plot to visualize

    ```
    Rscript /var/scratch/global/AfricaCDC_training/scripts/plotGenomecov.R COVM02379.coverage
    ```



## Variant calling
iVar uses the output of the ```samtools mpileup``` command to call variants - single nucleotide variants(SNVs) and indels. 


Pileup format consists of TAB-separated lines, with each line representing the pileup of reads at a single genomic position.

Several columns contain numeric quality values encoded as individual ASCII characters. Each character can range from "!" to "~" and is decoded by taking its ASCII value and subtracting 33; e.g., "A" encodes the numeric value 32.

The first three columns give the position and reference:

1. Chromosome name.
2. 1-bases position on the chromosome
3. Reference base at this position (this will be "N" on all lines if ```-f``` or ```--fasta-ref``` hs not been used)

In generating the mpileup, we will use the flags:
```--count-orphans```: Do not skip anomalous read pairs in variant calling. Anomalous read pairs are those marked in the FLAG field as paired in sequencing but without the properly-paired flag set.
```--ignore-overlaps```: Disable read-pair overlap detection
```--no-BAQ```: Disable base alignment quality (BAQ) computation

The ```tee``` command, used with a pipe, reads standard input from ```samtools mpileup```, then writes the output of the program to standard output and simultaneously copies it into the specified file ```.mpileup```


In order to call variants correctly, the reference file used for alignment must be passed to iVar using the ```-r``` flag. The output of samtools pileup is piped into ivar variants to generate a ```.tsv``` file with the variants. There are two parameters that can be set for variant calling using iVar - minimum quality (Default: 20) and minimum frequency (Default: 0.03). Minimum quality is the minimum quality for a base to be counted towards the ungapped depth to calculate iSNV frequency at a given position. For insertions, the quality metric is discarded and the mpileup depth is used directly. Minimum frequency is the minimum frequency required for a SNV or indel to be reported.

iVar can identify codons and translate variants into amino acids using a [GFF](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) format containing the required coding regions (CDS). In absence of a GFF file, iVar will not perform the translation and "NA" will be added to the output file in place of the reference and alternate codons and amino acids.

1. Change to the output directory ```ivar/variants```
    ```
    $ cd /var/scratch/global/$user/AfricaCDC_training/results/ivar/variants/
    ```  

2. Call variants

    ```
    $ samtools mpileup \
            --ignore-overlaps \
            --count-orphans \
            --no-BAQ \
            --max-depth 0 \
            --min-BQ 0 \
            --reference /var/scratch/global/$user/AfricaCDC_training/genome/nCoV-2019.fasta \
            /var/scratch/global/$user/AfricaCDC_training/results/ivar/trim/COVM02379.primertrimmed.sorted.bam \
            | tee COVM02379.mpileup \
            | ivar \
                variants \
                -t 0.25 \
                -q 20 \
                -m 10 \
                -g /var/scratch/global/$user/AfricaCDC_training/genome/nCoV-2019.gff \
                -r /var/scratch/global/$user/AfricaCDC_training/genome/nCoV-2019.fasta \
                -p COVM02379
    ```
3. Convert the variants from ```.tsv``` to ```.vcf``` (Variant Call Format)

    ```
    $ python /var/scratch/global/AfricaCDC_training/scripts/ivar_variants_to_vcf.py \
      COVM02379.tsv \
      COVM02379.vcf \
      --pass_only \
      --allele_freq_thresh 0.75 > COVM02379.variant.counts.log
    ```
    VCF file format
    
    The header begins the file and provides metadata describing the body     of the file. 
    Header lines are denoted as starting with #. 
    Special keywords in the header are denoted with ##.
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
    $ bgzip -c COVM02379.vcf > COVM02379.vcf.gz
    ```

5. Create tabix index from a sorted bgzip tab-delimited genome file
    ```
    $ tabix -p vcf -f COVM02379.vcf.gz
    ```
6. Generate stats from VCF file
    ```
    $ bcftools stats COVM02379.vcf.gz > COVM02379.stats.txt
    ```


## Variant annotation

We will use [SnpEff](http://pcingola.github.io/SnpEff/se_introduction/). It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes). It requires a configured SnpEff database with the annotation or features of the genome.

1. Change to the output directory ```snpeff```
    ```
    $ cd /var/scratch/global/$user/AfricaCDC_training/results/snpeff/
    ```  

2. Annotate and predict variants

    ```
    $ java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar \
        nCoV-2019 \
        -c /var/scratch/global/jjuma/AfricaCDC_training/genome/snpeff.config \
        -dataDir /var/scratch/global/$user/AfricaCDC_training/genome/snpeff_db/ \
        /var/scratch/global/$user/AfricaCDC_training/results/ivar/variants/COVM02379.vcf.gz \
        > COVM02379.vcf
      
    ```
3. Compress vcf file 
    ```
    bgzip -c COVM02379.vcf > COVM02379.vcf.gz
    ```
4. Rename the ```summary.html``` and ```genes.txt``` file
    ```
    $ mv snpEff_summary.html COVM02379.summary.html
    $ mv snpEff_genes.txt COVM02379.genes.txt
    ```

5. Create tabix index from a sorted bgzip tab-delimited genome file
    ```
    $ tabix -p vcf -f COVM02379.vcf.gz
    ```

6. Generate stats from VCF file
    ```
    $ bcftools stats COVM02379.vcf.gz > COVM02379.stats.txt
    ```

7. Filter variants
    [SnpSift](http://pcingola.github.io/SnpEff/ss_introduction/) annotates genomic variants using databases, filters, and manipulates genomic annotated variants. Once you annotated your files using SnpEff, you can use SnpSift to help you filter large genomic datasets in order to find the most significant variants for your experiment. 
    ```
      $ java -Xmx4g -jar /export/apps/snpeff/4.1g/SnpSift.jar \
            extractFields \
            -s "," \
            -e "." \
            COVM02379.vcf.gz \
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


## Consensus genome
To generate a consensus sequence iVar uses the output of samtools mpileup command. The mpileup output must be piped into ivar consensus. There are five parameters that can be set: 
- minimum quality ```-q``` (Default: 20)
- minimum frequency threshold ```-t``` (Default: 0)
- minimum depth to call a consensus ```-m``` (Default: 10)
- a flag ```-n``` to exclude nucleotides from regions with depth less than the minimum depth and a character to call in regions with coverage lower than the speicifed minimum depth (Default: 'N'). 

Minimum quality is the minimum quality of a base to be considered in calculations of variant frequencies at a given position. Minimum frequency threshold is the minimum frequency that a base must match to be called as the consensus base at a position. If one base is not enough to match a given frequency, then an ambigious nucleotide is called at that position. Minimum depth is the minimum required depth to call a consensus. If ```-k``` flag is set then these regions are not included in the consensus sequence. If ```-k``` is not set then by default, a 'N' is called in these regions. You can also specfy which character you want to add to the consensus to cover regions with depth less than the minimum depth. This can be done using ```-n``` option. It takes one of two values: ```-``` or ```N```.

1. Change to the output directory ```ivar/consensus```
    ```
    $ cd /var/scratch/global/$user/AfricaCDC_training/results/ivar/consensus/
    ```  

2. Generate pileup and consensus genome sequences

    ```
      $ samtools \
            mpileup \
            --reference /var/scratch/global/$user/AfricaCDC_training/genome/nCoV-2019.fasta \
            --count-orphans \
            --no-BAQ \
            --max-depth 0 \
            --min-BQ 0 \
            -aa \
            /var/scratch/global/$user/AfricaCDC_training/results/ivar/trim/COVM02379.primertrimmed.sorted.bam \
            | tee COVM02379.mpileup \
            | ivar \
                consensus \
                -t 0.75 \
                -q 20 \
                -m 10 \
                -n N \
                -p COVM02379
    ```
The ```tee``` command reads from the standard input and writes to both standard output and one or more files at the same time. ```tee``` is mostly used in combination with other commands through piping.

3. Rename the consensus genome header

    ```
    sed -i '/^>/s/Consensus_\(.*\)_threshold.*/\1/' COVM02379.fa
    ```


## Summarize results
Aggregate results from bioinformatics analyses across many samples into a single report with [MultiQC](https://multiqc.info/)

1. Change to the output directory `multiqc_out`
    ```
    $ cd /var/scratch/global/$user/AfricaCDC_training/results/multiqc/
    ```
2. Aggregate tools' outputs
    ```
    $ multiqc \
        --force \
        --title SARS-CoV-2 \
        --export \
        --outdir . \
        --config /var/scratch/global/AfricaCDC_training/assets/multiqc_config.yaml \
        /var/scratch/global/$user/AfricaCDC_training/results/
    ```
    
## Download MultiQC 

1. On your local computer, open a terminal and create a directory in the `Downloads` directory

    ```
    $ mkdir -p ~/Downloads/AfricaCDC_training/summary
    $ cd ~/Downloads/AfricaCDC_training/summary/
    ```
2. Copy the `MultiQC` reports, summaries and plots
    >**Note**
    >Replace ```user``` with the actual provided hpc account username
    ```
    $ rsync \
        -avP \
        --partial user@hpc.ilri.cgiar.org:/var/scratch/global/user/AfricaCDC_training/results/multiqc/SARS-CoV-2_multiqc_report.html \
        .
    
    rsync \
        -avP \
        --partial \
        user@hpc.ilri.cgiar.org:/var/scratch/global/user/AfricaCDC_training/results/ivar/coverage/*.pdf \
        .
    ```

## Pangolin Lineage Assignment
To assign [**Pangolin Lineages**](https://cov-lineages.org/lineage_list.html) to our consensus genomes we need to download the data from the HPC to your local computer.
1. Create a directory/folder on your local machine to store the data
    ```
    $ mkdir consensus
    $ rsync -r --progress --size-only <HPC-login-username>:/<path-to-the-directory-to-source-consensus-genomes>/ ./consensus/
    ```
On your browser open [Pangolin Web Application](https://pangolin.cog-uk.io/). This is the online version of [Pangolin](https://github.com/cov-lineages/pangolin)
Now copy and paste your consensus file to the site and click `Start Analysis`
Once done, download the results and if you need to, copy the names of sequences that failed the analysis to a file.


## Nextclade Analysis

## Data Retrieval and Review
Having sequenced our samples in both MiSeq (Illumina) and MinION (ONT) we can now transfer the sequence output to the HPC where we will conduct the bioinformatics analysis.
***
**Transfer of data from MiSeq:**  
MiSeq is based on Windows and data transfer will be done by copy-pasting the data to HPC through the network. Ensure the transfer is completed successfully without errors.  
***
**Transfer of data from MinION:**  
The MinION sequencer stores its sequencing output in a Linux based computer. To transfer the data we logged into computer and transferred the data on the command line as follows.
```
rsync -r --progress --size-only /<path-to-the-directory-with_sequencing-ouput>/ <HPC-login-username>:/<path-to-the-directory-to-store-sequencing-ouput>/
```
Replaced `<path-to-the-directory-with_sequencing-ouput>` with the path to the directory storing the sequencing output. Replaced `<HPC-login-username>` with your hpc login username (i.e user##@hpc.ilri.cgiar.org) and `<path-to-the-directory-to-store-sequencing-output>` with path to the directory you want to store the data in the HPC.
Example:`
rsync -r --progress /media/SeqData_LTS/20220405_1121_MN2816_FAH91436_2b8e9827/ user10@hpc.ilri.cgiar.org:/var/scratch/global/user10/20220405_1121_MN2816_FAH91436_2b8e9827`
***
**Reviewing the Illumina and ONT data:**
Log in to your HPC account.
```
ssh user##@hpc.ilri.cgiar.org
```
***
**Reviewing the Illumina data:**
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
***
**Reviewing the ONT data:**
Change working directory into directory that stores the ONT data. Note: replace `<name-of-run-folder>` with the name of the run folder
```
cd /var/scratch/global/ONT/
ls <name-of-run-folder>
```
