# README for NGSgenotyp

## Overview

Loci with extremely high levels of molecular polymorphism such as the self-incompatibility locus (S-locus) of Brassicaceae have remained recalcitrant to genotyping with NGS technologies based on short reads, as they are typically challenging to assemble de novo as well as to align to a given reference.
NGSgenotyp is an efficient pipeline to map raw reads from individual outcrossing Arabidopsis genomes against a dataset of multiple reference sequences of the pistil specificity determining gene of the Brassicaceae S-locus (SRK) and determine individual S-genotypes. In line with the important trans-specific polymorphism observed in this genetic system, we show that this approach can be first used to successfully obtain S-locus genotypes in related Brassicaceae genera, even if reference sequences from the species are not included in the initial database. We further show that this approach can be used to specifically assemble full-length individual S-allele sequences, and even discover new allelic sequences that were not initially present in the database.
This pipeline can in principle be adapted to other highly polymorphic loci, given datasets of relevant reference sequences are available.

## Reference list of dependencies

NGSgenotyp is written with Python 3 (>=3.9.6) and requires following packages installed before running:

	   biopython - >= 1.79
	   ete3 - >= 3.1.2
	   matplotlib >= 3.4.2
	   scipy - >= 1.7.0
	   numpy - >= 1.21.1
	   pysam - >= 0.16.0.1
	   pyyaml - >= 5.4.1
	   psutil - >= 5.9.4

NGSgenotyp also includes the following embedded tools (no need to install them):

	   bowtie - 2-2.2.6
	   CAP3
	   FastQC - 0.11.4
	   muscle - 3.8.31
	   PhyML  - 3.1
	   quast - 4.5
	   samtools - 1.4
	   SPAdes - 3.11.0
	   sratoolkit - 2.8.2-1
	   yass - 1.14

## Installing NGSgenotyp form the git repository

### 1a) Create a conda environment

Install conda : https://conda.io/projects/conda/en/latest/user-guide/install/index.html

```
sh create_conda_environment.sh
```

### 1b) Or install Python packages

```
pip install biopython ete3 matplotlib scipy numpy pysam pyyaml psutil
```

### 2) Create a repository for the pipeline and download it inside:

```
wget https://github.com/mathieu-genete/NGSgenotyp/archive/master.zip
unzip master.zip
```

## Running NGSgenotyp

If the conda environment is available, use the ***NGSgenotyp*** script. If not, use the ***NGSgenotyp2.py*** script.

```
usage: NGSgenotyp [feature]

features list:

	help            	-- show this help message

	version         	-- program version

	genotyp         	-- execute genotyp pipeline from raw NGS reads data

	haploAsm        	-- Assembly pipeline from genotyping results

```

### NGSgenotyp genotyp

```
usage: genotyp [-h] [-V] [-v] [-f] [-k] [-ks KMERSIZE] [-m MISMATCHTHRLD] [-pdf] [-s] [-sh] [-sm] [-kid] [-T MAXPARALLELSJOBS]
               [-e ERRORATETHRLD] [-M READSMAPPEDTHRLD] [-A ALIGNMENTMODE] [-S ALIGNMENTSENSITIVITY] -o OUTFOLDER -i READSINFO -d REFDATABASE
               [-x REFEXCLUDE] [-on OUTNBSHEETPERXLS]

NGSgenotyp v2 -- genotyping pipeline

optional arguments:
  -h, --help            show this help message and exit
  -V, --verbose         full verbose
  -v, --tinyverbose     verbose
  -f, --force
  -k, --kmerfilter      filtering fastq raw data with kmers dictionnary generated from references sequences. Should be use if your input
                        fastq are not yet filtered (significatively reduces compute time)
  -ks KMERSIZE, --kmerSize KMERSIZE
                        kmer size for kmer filtering - default = 20
  -m MISMATCHTHRLD, --mismatchthrld MISMATCHTHRLD
                        mismatch threshold on aligned reads (remove reads whose mitmatch value is greater than threshold)
  -pdf, --pdfreports    generate PDF reports (can take long time)
  -s, --statsonly       do stats only
  -sh, --sharedreads    generates shared reads file
  -sm, --statsmismatchthrld
                        do stats with filtered bam with mismatchthrld value
  -kid, --keepindel     force keep reads with insertion/deletion during bam filtering step
  -T MAXPARALLELSJOBS, --MaxParallelsJobs MAXPARALLELSJOBS
                        max number of parallels jobs to run - default= see config file
  -e ERRORATETHRLD, --ErroRateThrld ERRORATETHRLD
                        Force error rate threshold - default= see config file
  -M READSMAPPEDTHRLD, --readsMappedThrld READSMAPPEDTHRLD
                        Take into account only alignments with 'reads mapped' greater then threshold (default = 10)
  -A ALIGNMENTMODE, --alignmentMode ALIGNMENTMODE
                        bowtie2 reads alignments set to local or end-to-end (default = end-to-end)
  -S ALIGNMENTSENSITIVITY, --alignmentSensitivity ALIGNMENTSENSITIVITY
                        bowtie2 reads alignments sensitivity set to very-fast, fast, sensitive or very-sensitive (default = sensitive)
  -o OUTFOLDER, --outfolder OUTFOLDER
                        destination folder (create it if not exist)
  -i READSINFO, --readsinfo READSINFO
                        Configuration file with reads informations if reads ares paired add [format=paired] parameter
  -d REFDATABASE, --refdatabase REFDATABASE
                        reference database in fasta format (see documentation)
  -x REFEXCLUDE, --refexclude REFEXCLUDE
                        simple text file contains for each line, references names to exclude for current analysis
  -on OUTNBSHEETPERXLS, --outnbsheetperxls OUTNBSHEETPERXLS
                        maximum sheet number by output xls files (default = 15)
```

### NGSgenotyp haploasm

```
usage: haploasm [-h] [-V] [-v] [-f] [-p] [-y] [-x] [-T THREADED] [-st SPADESTHREADS] [-l INDIVLIST] [-m CONTIGMINLEN] [-M CONTIGMAXLEN] -g GENOTYPFOLDER -s ASMSUFFIX [-a ASSEMBLER]

NGSgenotyp Haplotyp Assembly

optional arguments:
  -h, --help            show this help message and exit
  -V, --verbose         full verbose
  -v, --tinyverbose     verbose
  -f, --force
  -p, --paired          orignals reads files with paired-end reads forward and reverse files distincts
  -y, --includeParalogsyass
                        include paralogs for yass analyse
  -x, --includeParalogstree
                        include paralogs for phylogeny
  -T THREADED, --threaded THREADED
                        define threads number used for phylogeny. By default 8 cpus
  -st SPADESTHREADS, --spadesThreads SPADESTHREADS
                        define threads number used for spades. By default 50% of cpu numbers
  -l INDIVLIST, --indivlist INDIVLIST
                        list of haplotype name for assembly
  -m CONTIGMINLEN, --contigminlen CONTIGMINLEN
                        keep contigs up from minimum lenght (default=500)
  -M CONTIGMAXLEN, --contigmaxlen CONTIGMAXLEN
                        keep contigs down to maximum lenght
  -g GENOTYPFOLDER, --genotypfolder GENOTYPFOLDER
                        genotyping folder
  -s ASMSUFFIX, --asmsuffix ASMSUFFIX
                        assembly folder suffix
  -a ASSEMBLER, --assembler ASSEMBLER
                        Program used to assemble alleles (spades,minia) defaut: spades
```

---

# Pipeline quick use

## Demo dataset

Datas format:

      sample DRS032518: concatenated reads - use the split option
      sample DRS032528: paired end datas - if split option is set, add ,format=paired,nosplit=True to the reads list file (see below)
      sample DRS032540: concatenated reads - use the split option

Expected genotype for all three samples

	 sample DRS032518: AhSRK12 and AhSRK01
	 sample DRS032528: AhSRK01
	 sample DRS032540: AhSRK02 and AhSRK19


format for the reads list file:

```
       /path/to/the/fastq/folder
       sample_name_1,fastq_file_1, fastq_file_2,...
       sample_name_2,fastq_Fwd_1,fastq_Rev_1,...,fastq_Fwd_n,fastq_Rev_n,format=paired,nosplit=True
       .
       .
       .
```

You can add a comment, only after the first line. The first line must always contain the absolute path to the fastq files folder.

- Add `format=paired` when you have paired data for a sample.

- To prevent files from being split, add `nosplit=True`.

- When you have paired data, use both parameters simultaneously: `...,format=paired,nosplit=True`.

## Databank format

The databank containing all your references must be a single file in fasta format. Avoid using spaces, commas, or slashes in your sequence IDs; use underscores instead.

Example:

```
>CgrSRK46|specie=Capsella_grandiflora|grpRef=H3-10
ACTCCAAATTTAATCAATCAAATGGATTTTTGTGGCAGAGC... ... ...
>Carubv10016249_Aly9|Paralog=1
ATGAGAGGCGCATTACCAAACTCTTACCATTCTTACACTTT... ... ...
>CgSRK07
GAGTGGAGAGATAGAGAGATGAGAAGTGAAGGACCAAAC... ... ...
...
...
```

The description of your sequences must follow this format:

```
>Sequence_ID|Param1=xxx|Param2=xxx... ...
```

- Sequence ID must be unique and always appear after the '>'.

- (Optional) Parameters must be separated by '|'. There is no precise order for their statement.

To indicate that your reference corresponds to a paralogue, add the parameter : Paralog=1 (example above with "Carubv10016249_Aly9")

You can specify the species with the parameter: `specie=xxxxx` (replace spaces with underscores).

**grpRef** uses two pieces of information in the following format: `H<group_nbr>-<group_id>` (e.g., H4-2: group number 4 and group ID number 2).

**grpId** is the allelic class ranging from 0 to infinity.

**HaploId** is the orthologous identifier within the allelic class (ranging from 1 to infinity), with numbering being arbitrary.

### *Example*:

*Table 1.*

| grpRef     | halleri       | lyrata        | Capsella grandiflora |
| ---------- | ------------- | ------------- | -------------------- |
| H0-1       | AhSRK27       | --            | CgrSRK40             |
| H1-1       | AhSRK01       | AlSRK01       | CgrSRK03             |
| H2-1       | AhSRK03       | AlSRK03       | --                   |
| H2-2       | AhSRK08       | --            | CgrSRK10             |
| **_H2-3_** | **_AhSRK09_** | **_AlSRK14_** | --                   |

In the table 1, hortologous references AhSRK09 (*A. halleri*) and AlSRK14 (*A. lyrata*) have the same grpRef:

They belong to group 2 (allelic class): gprId = 2

Because they ar hortologous, they have the same HaploId: HaploId = 3
 ==> grpRef = H2-3

So in the fasta file they will be noted as follows:

```
>AhSRK09|grpRef=H2-3
séquence... ...
>AlSRK14|grpRef=H2-3
séquence... ...
```

### you can define an allele from multiple sequences:

An allele can be composed of several sequences using the **allelePart** parameter. In this case, the sequence alone will no longer be considered as an allele.

example if you define an allele called **CgrSRK01** represented by 2 sequences (CgrSRK01P1,CgrSRK01P2) :

```
>CgrSRK01P1|allelePart=CgrSRK01|grpRefPart=H0-11
>CgrSRK01P2|allelePart=CgrSRK01|grpRefPart=H0-11
```

***grpRefPart*** and ***allelePart*** will be the same for each sequences

**CgrSRK01P1** and **CgrSRK01P2** will not appear in the results, only **CgrSRK01** (which combines these two sequences) will be displayed.

A sequence can be part of one or more alleles :

```
>ID_Sequence|allelePart=ID_Allele_1,ID_Allele_2,ID_Allele_n|grpRefPart=ID_groupe_1,ID_groupe_2,ID_groupe_n
>ID_Sequence2|allelePart=ID_Allele_2|grpRefPart=ID_groupe_2
>ID_Sequence3|allelePart=ID_Allele_1|grpRefPart=ID_groupe_1
```

**ID_Allele_1** contain **ID_Sequence** + **ID_Sequence3**

**ID_Allele_2** contain **ID_Sequence** + **ID_Sequence2**

---

## Basic command line usage for genotyping

Before using the demo dataset, edit the **reads_list_reduced** file in the demo_datas folder. Update the demo data path on the first line (absolute path required). Inside the NGSgenotyp folder, launch the following command:

```
./NGSgenotyp2.py genotyp -v -k -S -o demo_datas/demo_results -i demo_datas/reads_list_reduced -d demo_datas/Ahalleri_SRK_Database.fa
```

**Parameters:**

- **genotyp**: Use the genotyping feature of the pipeline.
- **-v**: Force NGSgenotyp to be verbose and print progress.
- **-k**: Filter raw fastq files using a kmer dictionary generated from the reference fasta.
- **-S**: Split concatenated reads from SRA into forward and reverse fastq files.
- **-o**: Specify the output folder for the results.
- **-i**: Provide the reads list file containing the absolute paths to fastq files, sample names, and fastq filenames (see format for the reads list file).
- **-d**: Specify the reference databank in fasta format (see Databank format).

### outputs

NGSgenotyp genotyp feature generates the following output folders:

- **BWT_Results**: Contains raw alignment data for each sample.
- **FQ_filtered**: Contains filtered fastq files generated by kmerRefFilter (if the -k option was set).
- **FQ_Splited**: Contains split fastq files (if the -S option was set).
- **logs**: Contains logs for NGSgenotyp and kmerRefFilter.
- **Results**: Contains all genotyping result files.

Results files:

- **Genotyp_Stats_XXXXXXXX.xls**: Excel file containing genotyping results.
- **putative_alleles.txt**: Tabular text file containing putative alleles for each sample.
- **putative_groups.txt**: Tabular text file containing putative groups for each sample.
- **xls_samples_index.txt**: Text file to locate a sample in the 'Genotyp_Stats_XXXXXXXX.xls' results files.
- **samtools_stats.yaml**: File contains raw results from the pipeline

## Genotyp_stats_XXXXXXXX.xls format

### Color Coding:

- **Orange**: Paralogous alleles, used as references for genotype score calculation.
- **Green**: Putative alleles.

Each sample has its own sheet.

### Column Descriptions:

1. **Allele Name**
2. **Genotype Score (gs) [0,1]**: The closer the genotype score is to 1, the better the prediction. It depends on the error rate and allele coverage.
3. **Error Rate**: This corresponds to the ratio of bases mapped to mismatches.
4. **Mean Covered Depth**: The average depth of reads that aligned with the reference. For a correct prediction, we expect an average depth of at least 5 for the paralogous sequences.
5. **Normalized Covered Depth**: This depth is normalized by the average depth of paralogous sequences. Comparing paralogous depth with sample depth can help determine heterozygosity.
   - If the depth of the allele matches the level of the paralogs, it indicates a homozygote (normalized values: ~1 for the paralogs and ~1 for the allele).
   - If there are two alleles and their depths are each about half of the paralogs' depth, it indicates a heterozygote (normalized values: ~1 for the paralogs and ~0.5 for each allele).
   - If there is only one predicted allele and its depth is about half of the paralogs, it indicates a heterozygote, but the second allele might not be present in the database.

6. **Homolog IDs**: This group corresponds to all orthologous sequences for an allele as defined in the database. Same group = same allele.
7. **Reference Length (bp)**
8. **Bases Mapped**
9. **Mismatches**
10. **Region Coverage**: Value between 0 and 1; 0 means the region is not covered, 1 means the aligned reads cover the entire region.
11. **Reads <X Mismatches**: Ratio of reads aligned with X mismatches. By default X=0 but the threshold can be modified with the `--mismatchthrld` option.
12. **Mean Covered Depth Mismatch Ratio Corrected**: This is the mean covered depth multiplied by the ratio of reads with **X** mismatches. It corrects the mean covered depth to reflect the alignment of reads under the mismatch threshold.
13. **Normalized Covered Depth Mismatch Ratio Corrected**: This is the normalized covered depth multiplied by the ratio of reads with **X** mismatches. It adjusts the normalized covered depth to account for the alignment of reads under the mismatch threshold.
14. **Warnings**: For example, if the reference length is shorter than 500 bp.

---

## Basic command line usage for haplotypes Assembly

The assembly pipeline is intended to assemble haplotypes individually. **To run haplotypes assembly, you need completed genotyping results.**

Launch haplotypes assembly for demo dataset:

```
NGSgenotyp2.py haploasm -g demo_datas/demo_results -s asmv1 -T 5
```

**Parameters:**

- **haploasm**: Use the NGSgenotyp haplotypes assembly feature.
- **-g**: Specify the genotyping folder from the 'genotyp' command (e.g., 'demo_datas/demo_results').
- **-s**: Provide the suffix for the assembly output folder (e.g., 'Assembly_asmv1').
- **-T**: Define the number of threads used for phylogeny (e.g., 5 CPUs).

### outputs

NGSgenotyp haploAsm feature generates the following output folders:

- **log**: Contains haploAsm and SPAdes log files.
- **XXXXXX**: Sample name folders containing alignment quality folders, SPAdes folders, and phylogeny folders.

## Configuration file

Some pipeline parameters can be adjusted in the **configs/genotyp_config.yaml** file.

- **MaxParallelsJobs**: Default is set to 20, which is the number of parallel jobs to run. Check your configuration to modify this value.

## Contact Information

Mathieu Genete

Email: mathieu.genete@univ-lille.fr

## Licence Agreement

This software is covered by GNU General Public License, version 3 (GPL-3.0).