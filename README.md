# README for NGSgenotyp

## Overview
Loci with extremely high levels of molecular polymorphism such as the self-incompatibility locus (S-locus) of Brassicaceae have remained recalcitrant to genotyping with NGS technologies based on short reads, as they are typically challenging to assemble de novo as well as to align to a given reference.
NGSgenotyp is an efficient pipeline to map raw reads from individual outcrossing Arabidopsis genomes against a dataset of multiple reference sequences of the pistil specificity determining gene of the Brassicaceae S-locus (SRK) and determine individual S-genotypes. In line with the important trans-specific polymorphism observed in this genetic system, we show that this approach can be first used to successfully obtain S-locus genotypes in related Brassicaceae genera, even if reference sequences from the species are not included in the initial database. We further show that this approach can be used to specifically assemble full-length individual S-allele sequences, and even discover new allelic sequences that were not initially present in the database.
This pipeline can in principle be adapted to other highly polymorphic loci, given datasets of relevant reference sequences are available.

## Reference list of dependencies
NGSgenotyp is written with Python 2.7.5 and requires following packages installed before running:

	   biopython - >= 1.68
	   ete3 - >= 3.0.0b35
	   matplotlib >= 2.0.0
	   scipy - >= 0.19.1
	   numpy - >= 1.13.1
	   pysam - >= 0.8.2.1

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
Create a repository for the pipeline and download it inside:
```
wget https://github.com/mathieu-genete/NGSgenotyp/archive/master.zip
unzip master.zip
```

## Running NGSgenotyp
```
usage: NGSgenotyp [feature]

features list:
	RapidFastQC     	-- A qualitative tool for quick quality check on fastq raw files (.fastq or .gz)

	help            	-- show this help message

	kmerRefFilter   	-- fastq raw files filtering by kmers dictionnary generated from references sequences

	extractFromFasta	-- extract sequences from multiple haploAsm contigs files

	splitFastq      	-- Split reads of one or more fastq file(s) in 2 files (_R1_RSplit.fastq and _R2_RSplit.fastq)

	genotyp         	-- execute genotyp pipeline from raw NGS reads data

	graphs_genotyp  	-- draw genotyps score plots in pdf file(s)

	SRAGetDatas     	-- download sequencing data and metadata from SRA number

	chMaxJob        	-- change the number of background jobs launched during genotyp pipeline execution.

	version         	-- program version

	databases       	-- list available databases

	haploAsm        	-- haplotypes assembly from genotyps mapped reads

```
### NGSgenotyp SRAGetDatas
```
usage: NGSgenotyp SRAGetDatas [-h] [-m] [-f] [-s SRAQUERY] [-l LISTOFSRA] -d DESTDIR
                              [-o CSVOUTPUTFILE] [-x ADDFILTER] [-F FILTERFROMREF]

Get SRA raw files and Metadatas

optional arguments:
  -h, --help            show this help message and exit
  -m, --metaDataOnly    get metadata only for given SRA
  -f, --localSratoolkit
  -s SRAQUERY, --sraquery SRAQUERY
                        SRA Query
  -l LISTOFSRA, --listOfSra LISTOFSRA
                        text file with list Of SRA Query (one per line)
  -d DESTDIR, --destdir DESTDIR
                        destination directory
  -o CSVOUTPUTFILE, --CSVoutputfile CSVOUTPUTFILE
                        csv output file
  -x ADDFILTER, --addfilter ADDFILTER
                        download data wich correspond to filer: header=value
  -F FILTERFROMREF, --filterFromRef FILTERFROMREF
                        directly filter reads while download with reference
```
### NGSgenotyp RapidFastQC
```
usage: NGSgenotyp RapidFastQC [-h] [-v] [-d OUTPUTDIR] [-n NUMBEROFREADS] -i
                              [FASTQFILE [FASTQFILE ...]]

A qualitative tool for quick quality check on fastq raw files.

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         show results on output
  -d OUTPUTDIR, --outputdir OUTPUTDIR
                        output directory for report
  -n NUMBEROFREADS, --numberOfReads NUMBEROFREADS
                        reads sampling. By default 2000
  -i [FASTQFILE [FASTQFILE ...]], --fastqfile [FASTQFILE [FASTQFILE ...]]
                        Fastq file(s) to analyse (fastq or gz file). For gz
                        file, reads number is not available
```
### NGSgenotyp kmerRefFilter
```
usage: NGSgenotyp kmerRefFilter [-h] [-v] [-prog] [-a] [-d] [-y] [-knf] [-o OUTPUTDIR]
                               [-k KMERSIZE] [-m MINMATCHKEEPSEQ] [-z MINSHENTROPY]
                               [-q MAXRATIOAMBIGOUS] [-e [EXTREMITY5P3P]]
                               [-r [REFERENCESFASTA [REFERENCESFASTA ...]]]
                               [-i KMERINFILE] [-p KMEROUTPUT]
                               [-x [EXCLUDEFASTA [EXCLUDEFASTA ...]]] [-maxRD MAXREADS]
                               [-f [FASTQFILE [FASTQFILE ...]]] [-l FASTQLIST]
                               [-1 [FWDFASTQ [FWDFASTQ ...]]]
                               [-2 [REVFASTQ [REVFASTQ ...]]]
                               [-u [FASTQURL [FASTQURL ...]]] [-u1 FWDFASTQURL]
                               [-u2 REVFASTQURL] [-ugzip] [-s STREAMINPUT]
                               [-c FWDREVCHOICE] [-kc KMERFWDREV]

Tool for raw reads kmer-based filtering.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -prog, --printprogress
                        print progession if set
  -a, --append          if fastq filtered exists, add new reads at the end of
                        file
  -d, --ambigousDNA     generate all possible kmers from an ambigous kmer
  -y, --yamlout         generate yaml output on stdout
  -knf, --keepNotFiltered
                        keep not filtered reads in a file _KEEPED.fastq
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        output directory for filtered fastq files
  -k KMERSIZE, --kmersize KMERSIZE
                        size for kmers - default=20
  -m MINMATCHKEEPSEQ, --minMatchKeepSeq MINMATCHKEEPSEQ
                        minimum number of ref reads matched to keep fastq read
                        -- default=1
  -z MINSHENTROPY, --minShEntropy MINSHENTROPY
                        minimum Shannon Entropy for kmers (0 to 2.0)--
                        defaut=0.8
  -q MAXRATIOAMBIGOUS, --maxratioAmbigous MAXRATIOAMBIGOUS
                        maximum ambigous bases accepted in kmers in % --
                        defaut=0.2
  -e [EXTREMITY5P3P], --extremity5p3p [EXTREMITY5P3P]
                        compare 5' and 3' read extremity only -- nbr of bases
                        to test from extremities default=1
  -r [REFERENCESFASTA [REFERENCESFASTA ...]], --referencesfasta [REFERENCESFASTA [REFERENCESFASTA ...]]
                        fasta files with references sequences
  -i KMERINFILE, --kmerinfile KMERINFILE
                        load kmer dictionary from a saved dictionary file
  -p KMEROUTPUT, --kmeroutput KMEROUTPUT
                        only save kmer dictionary in a file and exit without
                        filtering
  -x [EXCLUDEFASTA [EXCLUDEFASTA ...]], --excludefasta [EXCLUDEFASTA [EXCLUDEFASTA ...]]
                        fasta files with sequences to exclude
  -maxRD MAXREADS, --maxreads MAXREADS
                        maximum number of read to analyze
  -f [FASTQFILE [FASTQFILE ...]], --fastqfile [FASTQFILE [FASTQFILE ...]]
                        Fastq to filter (fastq, gz, bz2)
  -l FASTQLIST, --fastqlist FASTQLIST
                        text file with all fastq (fastq, gz, bz2) path -- 1
                        per line
  -1 [FWDFASTQ [FWDFASTQ ...]], --fwdfastq [FWDFASTQ [FWDFASTQ ...]]
                        forward Fastq to filter (fastq, gz, bz2) respect file
                        order with reverse files
  -2 [REVFASTQ [REVFASTQ ...]], --revfastq [REVFASTQ [REVFASTQ ...]]
                        reverse Fastq to filter (fastq, gz, bz2) respect file
                        order with forward files
  -u [FASTQURL [FASTQURL ...]], --fastqurl [FASTQURL [FASTQURL ...]]
                        url to Fastq to filter (fastq, gz)
  -u1 FWDFASTQURL, --fwdfastqurl FWDFASTQURL
                        url to forward Fastq to filter (fastq, gz)
  -u2 REVFASTQURL, --revfastqurl REVFASTQURL
                        url to reverse Fastq to filter (fastq, gz)
  -ugzip, --urlgzip     indicate that url's given in -u or in -u1 and -u2 are
                        gz files
  -s STREAMINPUT, --streamInput STREAMINPUT
                        use stdind as input - specify output filename
  -c FWDREVCHOICE, --fwdRevChoice FWDREVCHOICE
                        for paired filtering: filtering on forward reads only
                        (FWD), reverse reads only (REV) or both (BOTH - by
                        default)
  -kc KMERFWDREV, --kmerFwdRev KMERFWDREV
                        kmders dictionary generation: use forward kmers only
                        (FWD), reverse kmers only (REV) or both (BOTH - by
                        default)
```
### NGSgenotyp genotyp
```
usage: NGSgenotyp genotyp [-h] [-V] [-v] [-f] [-k] [-ks KMERSIZE] [-pdf] [-q] [-s]
                          [-r REGION REGION] [-T MAXPARALLELSJOBS] [-e ERRORATETHRLD]
                          [-S [SPLITREADS]] -o OUTFOLDER -i READSINFO -d REFDATABASE
                          [-x REFEXCLUDE] [--config CONFIG]

genotyping pipeline

optional arguments:
  -h, --help            show this help message and exit
  -V, --verbose         full verbose
  -v, --tinyverbose     verbose
  -f, --force
  -k, --kmerfilter      filtering fastq raw data with kmers dictionnary
                        generated from references sequences. Should be use if
                        your input fastq are not yet filtered (significatively
                        reduces compute time)
  -ks KMERSIZE, --kmerSize KMERSIZE
                        kmer size for kmer filtering - default = 20
  -pdf, --pdfreports    generate PDF reports (can take long time)
  -q, --samplesqc       samples quality control only
  -s, --statsonly       do stats only
  -r REGION REGION, --region REGION REGION
                        Region to analyse in sequences <min> <max>
  -T MAXPARALLELSJOBS, --MaxParallelsJobs MAXPARALLELSJOBS
                        max number of parallels jobs to run - default= see
                        config file
  -e ERRORATETHRLD, --ErroRateThrld ERRORATETHRLD
                        Force error rate threshold - default= see config file
  -S [SPLITREADS], --splitreads [SPLITREADS]
                        split fastq files reads in half use this option if
                        reads are concatenated - add [nosplit=True] parameter in
                        reads information file for sample not need to be split
  -o OUTFOLDER, --outfolder OUTFOLDER
                        destination folder (create it if not exist)
  -i READSINFO, --readsinfo READSINFO
                        Configuration file with reads informations if reads
                        ares paired add [format=paired] parameter
  -d REFDATABASE, --refdatabase REFDATABASE
                        reference database in fasta format (see documentation)
  -x REFEXCLUDE, --refexclude REFEXCLUDE
                        simple text file contains for each line, references
                        names to exclude for current analysis
  --config CONFIG       config file
```
### NGSgenotyp graphs_genotyp
```
usage: NGSgenotyp graphs_genotyp [-h] [-v] [-g] [-S] -s SAMTVIEWSTATS
                                  [-d GRAPHSPERDOCUMENT] -o OUTPDF
                                  [-f [REFFILTERNAME [REFFILTERNAME ...]]]

graphic for genotyp

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -g, --sortbygroup     sort by group
  -S, --ShowScore       print genotyp scores on plots
  -s SAMTVIEWSTATS, --samtviewstats SAMTVIEWSTATS
                        samview_stats.p file to analyse
  -d GRAPHSPERDOCUMENT, --graphsPerDocument GRAPHSPERDOCUMENT
                        number of graphs per documents (15 by default)
  -o OUTPDF, --outpdf OUTPDF
                        out pdf filename
  -f [REFFILTERNAME [REFFILTERNAME ...]], --refFilterName [REFFILTERNAME [REFFILTERNAME ...]]
                        filter for reference
```
### NGSgenotyp haploAsm
```
usage: NGSgenotyp haploAsm [-h] [-v] [-f] [-i] [-p] [-d] [-x] [-y] [-t THREADED]
                            [-st SPADESTHREADS] -db REFDATABASE [-k KMERSIZE] -o OUTPUTDIR
                            [-l INDIVLIST] [-m CONTIGMINLEN] [-M CONTIGMAXLEN]
                            [-mdc MAXDISTCONTIGS] -filFQ FILTEREDFQ -s SAMVIEWSTATS
                            [--config CONFIG]

NGSgenotyp Haplotyp Assembly

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -f, --force           force
  -i, --concatenated      orignals reads files with concatenated forward and
                        reverse paired-end reads
  -p, --paired          orignals reads files with paired-end reads forward and
                        reverse files distincts
  -d, --diploid         use with diploid highly polymorphic genomes
  -x, --excludeParalogs
                        exclude paralogs from yass results
  -y, --phylotree       generate phylogenetic tree
  -t THREADED, --threaded THREADED
                        define threads number used for phylogeny. By default
                        10% of cpu numbers
  -st SPADESTHREADS, --spadesThreads SPADESTHREADS
                        define threads number used for spades. By default 50%
                        of cpu numbers
  -db REFDATABASE, --refdatabase REFDATABASE
                        reference database used for genotyping
  -k KMERSIZE, --kmerSize KMERSIZE
                        comma-separated list of k-mer sizes for Spades (must
                        be odd and less then 128)
  -o OUTPUTDIR, --outputdir OUTPUTDIR
                        output directory for assembly
  -l INDIVLIST, --indivlist INDIVLIST
                        list of haplotype name for assembly
  -m CONTIGMINLEN, --contigminlen CONTIGMINLEN
                        keep contigs up from minimum lenght
  -M CONTIGMAXLEN, --contigmaxlen CONTIGMAXLEN
                        keep contigs down to maximum lenght
  -mdc MAXDISTCONTIGS, --maxDistContigs MAXDISTCONTIGS
                        contig maximum distance from nearest node - default =
                        see maxDistContigs value in configuration file
  -filFQ FILTEREDFQ, --filteredFQ FILTEREDFQ
                        Configuration file with filtered reads informations
  -s SAMVIEWSTATS, --samviewstats SAMVIEWSTATS
                        samview_stats file from genotyp results
  --config CONFIG       config file
```
### NGSgenotyp extractFromFasta
```
usage: NGSgenotyp extractFromFasta [-h] [-v] -f INFASTASPATH -o OUTFOLDER -l CONSTIGSLIST

Extract sequences from multiple NGSgenotyp haploAsm contigs files

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -f INFASTASPATH, --infastaspath INFASTASPATH
                        haploAsm results folder
  -o OUTFOLDER, --outfolder OUTFOLDER
                        output folder for fasta files
  -l CONSTIGSLIST, --constigslist CONSTIGSLIST
                        contigs list to extract in text format -- sample for a
                        line: out_fasta_filename_(without
                        exention),contig_Id_1,contig_Id_2,...,contig_Id_n
```
### NGSgenotyp chMaxJob
```
usage: NGSgenotyp chMaxJob [-h] [-u UTILSPARAM] -T MAXPARALLELSJOBS

Change Max parallels number of jobs runing

optional arguments:
  -h, --help            show this help message and exit
  -u UTILSPARAM, --utilsparam UTILSPARAM
                        'XXXX_utils_params' file
  -T MAXPARALLELSJOBS, --MaxParallelsJobs MAXPARALLELSJOBS
                        max number of parallels jobs to run
```
### NGSgenotyp databases
```
usage: NGSgenotyp databases [-h] [--config CONFIG] [-r REMOVEDB]

database list

optional arguments:
  -h, --help            show this help message and exit
  --config CONFIG       config file
  -r REMOVEDB, --removedb REMOVEDB
                        database to remove
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
>CgrSRK46|HaploId=10|gprId=3|specie=Capsella_grandiflora|grpRef=H3010
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

The following 3 parameters HaploId, gprId and grpRef are linked:

gprId is the allelic class ranging from 0 to 9

Chaque groupe contient des homologues HaploId dont la numérotation est arbitraire (entier entre 1 et 100)
HaploId is the homologous identifiant inside allelic class (ranging from 1 to 100).

grpRef uses these 2 informations in the following format: Hghhh (with g=gprId et hhh=HaploId)

### *Example*:

*Table 1.*

|  grpRef |  halleri 	| lyrata  | Capsella grandiflora  |
|		:-:		|		:-:			|	:-:			|		:-:									|
|  H0001	|   AhSRK27 |		--		|	CgrSRK40							|
|  H1001	|   AhSRK01	| AlSRK01	| CgrSRK03 							|
|  H2001	|   AhSRK03	| AlSRK03	|		--				        	|
|  H2002	|   AhSRK08	|		--		|	CgrSRK10							|
|  H2003	|   AhSRK09	| AlSRK14	|   --									|

In the table 1, hortologous references AhSRK09 (A. halleri) and AlSRK14 (A. lyrata) have the same grpRef:

They belong to group 2 (allelic class): gprId = 2
Because they ar hortologous, they have the same HaploId: HaploId = 3
 ==> grpRef = H2003

So in the fasta file they will be noted as follows:

```
>AhSRK09|gprId=2|HaploId=3|grpRef=H2003
séquence... ...
>AlSRK14|gprId=2|HaploId=3|grpRef=H2003
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
