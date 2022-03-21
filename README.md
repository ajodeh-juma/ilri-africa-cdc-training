**Building capacity in SARS-CoV-2 genomics in Africa**.

- [Introduction](#introduction)
  - [Background](#background)
  - [Genome preparation](#genome-preparation)
  - [Metagenomics](#metagenomics)
  - [Mapping/Alignment](#mapping/alignment)
  - [Sorting and indexing](#sorting-and-indexing)
  - [Primer trimming](#primer-trimming)
  - [Compute coverage](#compute-and-plot-genome-sequencing-coverage)
  - [Variant calling](#variant-calling)
  - [Variant annotation](#variant-annotation)
  - [Consensus genome](#consensus-genome)



## Introduction
In early January 2020, the novel coronavirus (SARS-CoV-2) responsible for a 
pneumonia outbreak in Wuhan, China, was identified using next-generation
sequencing (NGS) and readily available bioinformatics pipelines. In addition to 
virus discovery, these NGS technologies and bioinformatics resources are
currently being employed for ongoing genomic surveillance of SARS-CoV-2
worldwide, tracking its spread, evolution and patterns of variation on a global
scale.

Here we apply basic principles employed by the numerous bioinformatic pipelines
to generate consensus genome sequences of SARS-CoV-2 and identfy variants using
a dataset obtained by sequencing samples from Kenya.

This is part of the initiative fronted by the [Africa
CDC](https://africacdc.org/) with generous support from the [Rockeffeler 
foundation](https://www.rockefellerfoundation.org/) to build pathogen genomics
in Africa.

## Background

We will use a dataset comprising of raw sequence reads of SARS-CoV-2 samples 
obtained from a sequencing run on NextSeq 550 platorm at [ILRI](www.ilri.org). 
NextSeq 550 flowcell uses 4 lanes and so 4 reads data per sequenced sample corresponding to the 4
lanes are generated with suffixes L001, L002, L003 and L004. The dataset we are 
using in this tutorial comprises of already concatenated sequences. These reads
can be combined/concatenated into a single file bearing in mind the type of library
sequencing either ```single``` or ```paired-end```.

**NOTE**  
In this tutorial, replace all instances of ```<user>``` with the output from ```whoami``` command. You can save the output of ```whoami``` as a variable and user it in subsequent steps
```
user=`whoami`
echo $user
```

For example to create a directory with my user name, I can simply execute these commands:
```
user=`whoami`
mkdir $user
```


### Prep

1. Project structure 

    data directory -
    ```/var/scratch/global/AfricaCDC_training/testdata_nextseq/testdata_nextseq```  
    1. Change to the global temporary directory
        ```
        cd /var/scratch/global/
        ```
    2. Find your user id with the command
        ```
        whoami
        ```
    3. Make a directory with the output from ```whoami``` command by storing the
       output in a variable(here we are storing it as ```user``` but you can use
       any name without spaces) and using it to create a directory
        ```
        user=`whoami`
        mkdir $user
        ```
    4. change to the newly created subdirectory in step 3 and make a project directory named ```sarscov2``` 
        ```
        cd $user
        mkdir sarscov2
        ```  
    5. Make a subdirectory in the ```sarscov2``` dircetory to store data
        ```
        mkdir data
        ```  
    6. Make a subdirectory in the parent dircetory to store reference genome  
        ```
        mkdir genome
        ```  
    7. Make a subdirectory in the parent dircetory to store results  
        ```
        mkdir results
        ```  

        **Optional**  
        The above command can be executed using the command below.  
        The ```mkdir``` command allows one to create intermediate directories with the option ```-p```.  
        ```
          cd /var/scratch/global/
          user=`whoami`
          mkdir -p $user/sarscov2/{data,genome,results}
        ```  
      8. Change to the project directory ```sarscov2``` and copy the raw
         sequence data using ```cp``` command to the ```data``` subdirectory.

          ```
          cd /var/scratch/global/$user/sarscov2
          cp /var/scratch/global/AfricaCDC_training/testdata_nextseq/testdata_nextseq/* data/
          ```

2. Requirements  
    Clear the environment
    ```
    module purge
    ```

    Load required modules
    ```
    module load kraken/2.0.8-beta
    module load bowtie2/2.3.4.1
    module load samtools/1.11
    module load ivar/1.3.1
    module load bcftools/1.11
    module load snpeff/4.1g
    ```


3.  SARS-CoV-2 genome  

    genomes directory - ```/var/scratch/global/AfricaCDC_training/genomes/sars-cov-2```  
    - fasta - ```nCoV-2019.reference.fasta```  
    - gff - ```GCA_009858895.3_ASM985889v3_genomic.200409.gff```  
    - Copy the ```fasta``` and the ```gff3``` files using ```cp``` command to the ```genome``` subdirectory you created.

## Genome preparation
Prepare reference genome for several downstream processes  

- index reference genome with ```samtools```. This is essential for variants calling  

    ```
    samtools faidx /var/scratch/global/$user/sarscov2/genome/nCoV-2019.reference.fasta
    ```  

    ```
    cut -f 1,2 /var/scratch/global/$user/sarscov2/genome/nCoV-2019.reference.fasta.fai > /var/scratch/global/<user>/sarscov2/genome/nCoV-2019.reference.fasta.sizes
    ``` 
- index reference genome with ```bowtie2```.  

    ```
    mkdir /var/scratch/global/$user/sarscov2/genome/bowtie2
    ```

    ```
    bowtie2-build \
      --threads 1 \
      /var/scratch/global/$user/sarscov2/genome/nCoV-2019.reference.fasta \
      /var/scratch/global/$user/sarscov2/genome/bowtie2/nCoV-2019.reference
    ```


## Metagenomics

At times, one may be interested in finding out the composition of reads within a
sequenced sample. This is often the case when the sequencing library was
prepared using a metagenomics approach where the total nucleic acid material is
sequenced without a bias (no enrichment of the target pathogen or material).
There are several bioinformatics tools and databases which can be used in
querying the read data in order to ascertain the taxonomic composition of a
particular sample. One of the commonly used tool is
[Kraken2](https://github.com/DerrickWood/kraken2/wiki). ```Kraken2``` can allow
us to query the composition of our samples by searching for sequence reads
against a pre-formatted database.

#### Kraken2 database
In this tutorial, we will use a pre-formatted human database to see the level of
human-derived reads in our samples. This may give us an indication of
contamination from host reads. 

1. load module  
    ```
    module load kraken/2.0.8-beta
    ```  
2.  make a directory in the ```results``` subdirectory within the parent
  ```sars-cov-2-genomics``` directory  
    ```
    mkdir /var/scratch/global/$user/sarscov2/results/kraken2_output
    ```  

3. run the following commands for individual samples:  

    ```
    cd /var/scratch/global/$user/sarscov2/results/kraken2_output
    ```
    ```
    kraken2 \
      -db /var/scratch/global/AfricaCDC_training/databases/kraken2-human-db \
      --threads 1 \
      --unclassified-out COVM02379_S37.unclassified#.fastq \
      --classified-out COVM02379_S37.classified#.fastq \
      --report COVM02379_S37.kraken2.report.txt \
      --output COVM02379_S37.kraken2.out.txt \
      --gzip-compressed \
      --report-zero-counts \
      --paired /var/scratch/global/$user/sarscov2/data/COVM02379_S37_con_R1_001.fastq.gz \
      /var/scratch/global/$user/sarscov2//data/COVM02379_S37_con_R2_001.fastq.gz
    ```

Run step 3. above for the other 2 samples.

Q1. How many reads have hit the human genome as targets in the samples?
<details close>
  <summary>Answer</summary>
  COVM02379_S37: 66438 reads (13.87%)
  COVM02524_S209: 60686 reads (57.23%)
  COVM02720_S365: 245351 reads (83.17%)
</details>
Q2. Based on the results, which sample do you think will give us good results in terms of genome coverage?


## Mapping/Alignment
Aligning sequencing reads to a reference genome is the first step in many
comparative genomics pipelines, including pipelines for variant calling, 
isoform quantitation and differential gene expression. In many cases, the
alignment step is the slowest. This is because for each read the aligner must 
solve a difficult computational problem: determining the read's likely point of
origin with respect to a reference genome. This is always trivial for several
reasons:
- The reference genome is really big. Searching big things is harder than
  searching small things.
- You aren’t always looking for exact matches in the reference genome–or, at 
least, probably not.
 
Here we use [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml),
an ultrafast and memory-efficient tool for aligning sequencing reads to long
reference sequences.


1. make a directory in the ```results``` subbdirectory to store the output and
  change to the directory
    ```
    mkdir /var/scratch/global/$user/sarscov2/results/bowtie2_output
    cd /var/scratch/global/$user/sarscov2/results/bowtie2_output
    ```  

2. run the command to align reads to the reference genome

    ```
    bowtie2 \
      -x /var/scratch/global/$user/sarscov2/genome/bowtie2/nCoV-2019.reference \
      -1 /var/scratch/global/$user/sarscov2/data/COVM02379_S37_con_R1_001.fastq.gz \
      -2 /var/scratch/global/$user/sarscov2/data/COVM02379_S37_con_R2_001.fastq.gz \
      --threads 1 \
      --un-conc-gz COVM02379_S37.unmapped.fastq.gz \
      --local \
      --very-sensitive-local \
      --seed 1 \
      2> COVM02379_S37.bowtie2.log \
      | samtools view -@ 1 -F4 -bhS -o COVM02379_S37.bam -
    ```

Run step 2. above for the other 2 samples.


## Sorting and indexing

Sort and index the alignment with ```samtools```

1. make a directory in the ```results``` subbdirectory to store the output and
  change to the directory
    ```
    cd /var/scratch/global/$user/sarscov2/results/bowtie2_output
    ```  

2. run the command

    ```
    samtools sort -@ 1 -o COVM02379_S37.sorted.bam -T COVM02379_S37 COVM02379_S37.bam
    samtools index -@ 1 COVM02379_S37.sorted.bam
    ```

Run step 2. above for the other 2 samples.



## Primer trimming
Trim amplicon primers using ```ivar```

1. make a directory in the ```results``` subbdirectory to store the output and
  change to the directory
    ```
    mkdir /var/scratch/global/$user/sarscov2/results/ivar_output/primer_trimmed
    cd /var/scratch/global/$user/sarscov2/results/ivar_output/primer_trimmed
    ```  

2. run the command

    ```
    ivar trim \
      -i /var/scratch/global/$user/sarscov2/results/bowtie2_output/COVM02379_S37.sorted.bam \
      -b /var/scratch/global/AfricaCDC_training/primer-schemes/V3/nCoV-2019.primer.bed \ 
      -p COVM02379_S37.primertrimmed \
      -m 30 \
      -q 20 > COVM02379_S37.ivar.log
        samtools sort \
          -@ 1 \
          -o COVM02379_S37.primertrimmed.sorted.bam \
          -T COVM02379_S37 COVM02379_S37.primertrimmed.bam
    ```

Run step 2. above for the other 2 samples.



## Variant calling
Call variants with ivar 

1. make a directory in the ```results``` subbdirectory to store the output and
  change to the directory
    ```
    mkdir /var/scratch/global/$user/sarscov2/results/ivar_output/variants
    cd /var/scratch/global/$user/sarscov2/results/ivar_output/variants
    ```  

2. run the command

    ```
    samtools mpileup \
            --ignore-overlaps \
            --count-orphans \
            --no-BAQ \
            --max-depth 0 \
            --min-BQ 0 \
            --reference /var/scratch/global/$user/sarscov2/genome/nCoV-2019.reference.fasta \
            /var/scratch/global/$user/sarscov2/results/ivar_output/primer_trimmed/COVM02379_S37.primertrimmed.sorted.bam \
            | tee COVM02379_S37.mpileup \
            | ivar \
                variants \
                -t 0.25 \
                -q 20 \
                -m 10 \
                -g ${gff} \
                -r ${fasta} \
                -p COVM02379_S37
    ```

    ```
    python /var/scratch/global/AfricaCDC_training/scripts/ivar_variants_to_vcf.py \
      COVM02379_S37.tsv \
      COVM02379_S37.vcf \
      --pass_only \
      --allele_freq_thresh 0.75 > COVM02379_S37.variant.counts.log

    # compress vcf file
    bgzip -c COVM02379_S37.vcf > COVM02379_S37.vcf.gz

    # create tabix index from a sorted bgzip tab-delimited genome file
    tabix -p vcf -f COVM02379_S37.vcf.gz

    generates stats from VCF file
    bcftools stats COVM02379_S37.vcf.gz > COVM02379_S37.bcftools.stats.txt
    ```

Run step 2. above for the other 2 samples.



## Compute and plot genome sequencing coverage

1. make a directory in the ```results``` subbdirectory to store the output and
  change to the directory
    ```
    mkdir -p /var/scratch/global/$user/sarscov2/results/mosdepth/genome
    cd /var/scratch/global/$user/sarscov2/results/mosdepth/genome
    ```  

2. run the command
    ```
    mosdepth \
        --by 200 \
        --fast-mode \
        COVM02379_S37 \
        /var/scratch/global/$user/sarscov2/results/ivar_output/primer_trimmed/COVM02379_S37.primertrimmed.sorted.bam
    ```

cd /var/scratch/global/$user/sarscov2/results/ivar_output/primer_trimmed


## Variant annotation

Annotate variants to predict variants effects

1. make a directory in the ```results``` subbdirectory to store the output and
  change to the directory
    ```
    mkdir /var/scratch/global/AfricaCDC_training/results/snpeff_output
    cd /var/scratch/global/AfricaCDC_training/results/snpeff_output
    ```  

2. run the command

    ```
    java -Xmx4g -jar /export/apps/snpeff/4.1g/snpEff.jar \
      nCoV-2019.reference \
      -c /var/scratch/global/AfricaCDC_training/databases/sars-cov-2-snpeff_db/snpeff.config \
      -dataDir /var/scratch/global/AfricaCDC_training/databases/sars-cov-2-snpeff_db/data \
      /var/scratch/global/AfricaCDC_training/results/snpeff_output/COVM02379_S37.vcf.gz \
      | bgzip -c > COVM02379_S37.snpeff.vcf.gz

      mv snpEff_summary.html COVM02379_S37.snpeff.summary.html
      mv snpEff_genes.txt COVM02379_S37.snpeff.genes.txt

      # create tabix index from a sorted bgzip tab-delimited genome file
      tabix -p vcf -f COVM02379_S37.snpeff.vcf.gz

      # generates stats from VCF file
      bcftools stats COVM02379_S37.snpeff.vcf.gz > COVM02379_S37.bcftools.stats.txt

      # Filter variants
      java -Xmx4g -jar /export/apps/snpeff/4.1g/SnpSift.jar \
            extractFields \
            -s "," \
            -e "." \
            COVM02379_S37.snpeff.vcf.gz \
            CHROM POS REF ALT \
            "ANN[*].GENE" "ANN[*].GENEID" \
            "ANN[*].IMPACT" "ANN[*].EFFECT" \
            "ANN[*].FEATURE" "ANN[*].FEATUREID" \
            "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \
            "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \
            "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \
            "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \
            "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \
            > COVM02379_S37.snpsift.txt
    ```

Run step 2. above for the other 2 samples


## Consensus genome

1. make a directory in the ```results``` subbdirectory to store the output and
  change to the directory
    ```
    mkdir /var/scratch/global/$user/sarscov2/results/ivar_output/consensus
    cd /var/scratch/global/$user/sarscov2/results/ivar_output/consensus
    ```  

2. run the command

    ```
      samtools \
            mpileup \
            --reference /var/scratch/global/$user/sarscov2/genome/nCoV-2019.reference.fasta \
            --count-orphans \
            --no-BAQ \
            --max-depth 0 \
            --min-BQ 0 \
            -aa \
            /var/scratch/global/$user/sarscov2/results/ivar_output/primer_trimmed/COVM02379_S37.primertrimmed.sorted.bam \
            | tee COVM02379_S37.mpileup \
            | ivar \
                consensus \
                -t 0.75 \
                -q 20 \
                -m 10 \
                -n N \
                -p COVM02379_S37
    ```

Run step 2. above for the other 2 samples





