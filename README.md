# README for NGSgenotyp

## Overview
Loci with extremely high levels of molecular polymorphism such as the self-incompatibility locus (S-locus) of Brassicaceae have remained recalcitrant to genotyping with NGS technologies based on short reads, as they are typically challenging to assemble de novo as well as to align to a given reference.
NGSgenotyp is an efficient pipeline to map raw reads from individual outcrossing Arabidopsis genomes against a dataset of multiple reference sequences of the pistil specificity determining gene of the Brassicaceae S-locus (SRK) and determine individual S-genotypes. In line with the important trans-specific polymorphism observed in this genetic system, we show that this approach can be first used to successfully obtain S-locus genotypes in related Brassicaceae genera, even if reference sequences from the species are not included in the initial database. We further show that this approach can be used to specifically assemble full-length individual S-allele sequences, and even discover new allelic sequences that were not initially present in the database.
This pipeline can in principle be adapted to other highly polymorphic loci, given datasets of relevant reference sequences are available.

## Reference list of dependencies
NGSgenotyp is written with Python 3 and requires following packages installed before running:

	   biopython - >= 1.79
	   ete3 - >= 3.1.2
	   matplotlib >= 3.4.2
	   scipy - >= 1.7.0
	   numpy - >= 1.21.1
	   pysam - >= 0.16.0.1
	   pyyaml - >= 5.4.1
	   psutil - >= 5.9.4

NGSgenotyp also use embeded following tools (no need to install them):

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

```
conda create -n NGSgenotyp2_env python=3.9.6 matplotlib biopython scipy numpy pysam ete3 pyyaml xlrd xlwt psutil
```
### 1b) Or install Python packages
```
biopython, ete3, matplotlib, scipy, numpy, pysam, pyyaml, psutil
```

### 2) Create a repository for the pipeline and download it inside:
```
wget https://github.com/mathieu-genete/NGSgenotyp/archive/master.zip
unzip master.zip
```

## Running NGSgenotyp
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
### NGSgenotyp haploAsm
```
usage: haploAsm [-h] [-V] [-v] [-f] [-p] [-y] [-x] [-T THREADED] [-st SPADESTHREADS] [-l INDIVLIST] [-m CONTIGMINLEN] [-M CONTIGMAXLEN] -g
                GENOTYPFOLDER -s ASMSUFFIX [-a ASSEMBLER]

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

add "format=paired" when you have paired datas for a sample.
to prevent files from being split add "nosplit=True"
When you have paired datas, you should use both parameters at the same time: ....,format=paired,nosplit=True
## Databank format
The databank containing all your references must be a single file in fasta format.

You should not use spaces, commas, slash in your sequences id, use the character "_" instead.

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

sequence id must be unique and always be after the sign ">".
the (optional) parameters must be separated by the "|". There is no precise order in their statement.

To indicate that your reference corresponds to a paralogue you must add the parameter: Paralog=1 (example above with "Carubv10016249_Aly9")

You can specify the species with the parameter (be careful, replace spaces with "_"): specie=xxxxx

grpRef uses 2 informations in the following format: Hg-h (with g=gprId et h=HaploId)

gprId is the allelic class ranging from 0 to infinite

HaploId is the hortologous identifiant inside allelic class (ranging from 1 to infinite), the numbering choice of your HaploId is arbitrary.

### *Example*:

*Table 1.*

|  grpRef |  halleri 	| lyrata  | Capsella grandiflora  |
|		:-:		|		:-:			|	:-:			|		:-:									|
|  H0-1	|   AhSRK27 |		--		|	CgrSRK40							|
|  H1-1	|   AhSRK01	| AlSRK01	| CgrSRK03 							|
|  H2-1	|   AhSRK03	| AlSRK03	|		--				        	|
|  H2-2	|   AhSRK08	|		--		|	CgrSRK10							|
|  **_H2-3_**	|   **_AhSRK09_**	| **_AlSRK14_**	|   --									|

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
## Basic command line usage for genotyping
Before use demo dataset, edit **reads_list_reduced** file in demo_datas folder. Update demo datas path on the first line (absolute path required).
Inside NGSgenotyp folder launch the following command:
```
./NGSgenotyp.py genotyp -pdf -v -k -S -o demo_datas/demo_results -i demo_datas/reads_list_reduced -d demo_datas/Ahalleri_SRK_Database.fa

genotyp => use genotype feature of the pipeline
-pdf => generate pdf plots coverage informations for each samples
-v => force NGSgenotyp to be verbose and print progression
-k => filter raw fastq with kmers dictionnary generated from the references fasta
-S => reads from SRA are concatenated. Split reads in half for each fastq file and generate 2 fastq file (forward and reverse)
-o => output folder for the demo_results
-i => reads list: text file contains fastq files absolute path, samples name and fastq file names (see format for the reads list file)
-d => reference databank in fasta (see Databank format)
```
### outputs

NGSgenotyp genotyp feature generates following output folders:

- BWT_Results => contains raw alignments data for each samples
- FQ_filtered => contains filtered fastq generated by kmerRefFilter (if the -k option was set)
- FQ_Splited => contains splited fastq (if the -S option was set)
- logs => contains NGSgenotyp and kmerRefFilter logs
- QC_SamplesRef => contains FastQC results for each samples
- Results => contains all genotyping results files (see below)

Results files:

- Coverage_sample_xxxxx.pdf => plot coverage for each sample (if -pdf option is set). By default, only aligments with an error rate lower than 0.08 are displayed (threshold can be modified in configuration file)
- ErrCovDepth_plot.pdf => scatterplots for all aligments results / Allele probability plots
- ErrRegionCov_Plot.pdf => coverage/error rate scatterplots for each sample
- Genotyp_Stats_XXXXXXXX.xls => xls file contains genotyping results
- putative_alleles.txt => tabular text file contains putative alleles for each sample
- shared_reads.xls => contains numbers of shared reads between aligments for each sample

## Basic command line usage for haplotypes Assembly

To run haplotypes assembly, you need genotyping results.

for remote SSH access, **ete3** python module need export display for pdf generation (ssh -X)

Launch haplotypes assembly for demo dataset:

```
./NGSgenotyp.py haploAsm -p -d -x -y -db demo_datas/Ahalleri_SRK_Database.fa -o demo_datas/demo_results_2/Assembly -m 900 -filFQ demo_datas/demo_results_2/FQ_Splited/Ahalleri_SRK_Database_splited_reads_list.yaml -s demo_datas/demo_results_2/Results/samview_stats.p

haploAsm => use NGSgenotyp haplotypes assembly feature
-p => originals reads files with paired-end reads forward and reverse files distincts
-d => use DipSpades for diploid highly polymorphic genomes
-x => exclude paralogous sequences (identified as paralog in database)
-y => generates phylogenetic trees one for each samples and global tree
-db => reference database used for genotyping
-o => output folder for results
-m => keep contigs up from minimum length (900 bp)
-filFQ => text or yaml file contains fastq filtered files absolute path, samples name and fastq filtered file names. NGSgenotyp generates a yaml file for filtered fastq, which can be used by haploAsm assembly feature
-s => samview_stats.p file from genotyping results
```
### outputs

NGSgenotyp haploAsm feature generates following output folders:

- log => folder contains haploAsm and spades logs files
- XXXXXX => sample names folders
	- align_qual => folder
	- DipSpades_OUT/Spades_OUT => spades folders
	- XXXXXX_phylogeny

## Configuration file
Some pipeline parameters can be adjust in the **config.yaml** file.

- MaxParallelsJobs: by default set to 10, is the number of parallels jobs to run. Check your configuration to modifie this value

 If necessary, you can change MaxParallelsJobs value during pipeline execution by using:
 ```
 ./NGSgenotyp chMaxJob -T <number of max parallels jobs>
 ```
## Contact Information
Mathieu Genete

Email: mathieu.genete@univ-lille.fr

## Licence Agreement
This software is covered by GNU General Public License, version 3 (GPL-3.0).
