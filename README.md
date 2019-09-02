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

NGSgenotyp also use embeded following tools:

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
                        reads are interlaced - add [nosplit=True] parameter in
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
  -i, --interlaced      orignals reads files with interlaced forward and
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

      sample DRS032518: interlaces data - use the split option
      sample DRS032528: paired end datas - if split option is set, add ,format=paired,nosplit=True to the reads list file (see below)
      sample DRS032540: interlaces data - use the split option

Expected genotype for all three samples

	 sample DRS032518: AhSRK12 and AhSRK01
	 sample DRS032528: AhSRK01
	 sample DRS032540: AhSRK02 and AhSRK19


format for the reads list file:
       /path/to/the/fastq/folder
       sample_name_1,fastq_file_1, fastq_file_2,...
       sample_name_2,fastq_Fwd_1,fastq_Rev_1,...,fastq_Fwd_n,fastq_Rev_n,format=paired,nosplit=True
       .
       .
       .

You can add a comment, only after the first line. The first line must always contain the absolute path to the fastq files folder.

add "format=paried" when you have paired datas for a sample.
to prevent files from being split add "nosplit=True"
When you have paired datas, you should use both parameters at the same time: ....,format=paired,nosplit=True
## Databank format

## Basic command line usage for genotyping


### Contact Information
Mathieu Genete

Email: mathieu.genete@univ-lille.fr

### Licence Agreement
This software is covered by GNU General Public License, version 3 (GPL-3.0).
