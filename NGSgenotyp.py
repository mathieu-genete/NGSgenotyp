#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

AppPath = os.path.realpath(__file__)
App_Folder = AppPath[:AppPath.rfind('/')]
tools_folder = os.path.join(App_Folder,'TOOLS')

#Library path
PyLibPath = "PyLib/"
sys.path.insert(0, os.path.join(App_Folder,PyLibPath))

import argparse
import subprocess
import collections
import yaml
import xlrd
import xlwt
import pysam
from Bio import SeqIO
from Bio import SeqUtils
import utils
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import random
import math
import numpy as np
from scipy.stats import gaussian_kde
from itertools import cycle
import pickle
import hashlib
import shutil
from datetime import datetime
import time
import threading
import ctypes
import urlparse
from urllib2 import urlopen
from operator import itemgetter

plt.switch_backend('agg')
utils.set_paramFileName('NGSgenotyp')

#Globals variables
__version__="v1.4"
AppName = "NGSgenotyp"

config = None
fastqParams = None
args = None
HomoParaFromRef = {}
StderrPaths = []
samtools_path=''
bowtie2_path = ''
fastqc_path = ''

#*************************
#* genotyp main function *
#*      Workflow         *
#*************************

def genotyp(ArgsVal):

	global config
        global fastqParams
	global args
	global HomoParaFromRef
	global StderrPaths
	global samtools_path
	global bowtie2_path
	global fastqc_path

	startTime = time.time()
	
	#ArgsParser
	description = """ genotyping pipeline """
	
	parser = argparse.ArgumentParser(prog="genotyp",description=description)
	
	parser.add_argument("-V", "--verbose", help="full verbose", action="store_const", const="-V", default="")
	parser.add_argument("-v", "--tinyverbose", help="verbose",action="store_const", const="-v", default="")
	parser.add_argument("-t", "--test", help=argparse.SUPPRESS, action="store_const", const="-t", default="")
	parser.add_argument("-f", "--force", action="store_const", const="-f", default="")
	parser.add_argument("-c", "--checks", help=argparse.SUPPRESS, action="store_const", const="-c", default="")
	parser.add_argument("-k", "--kmerfilter", help="filtering fastq raw data with kmers dictionnary generated from references sequences. Should be use if your input fastq are not yet filtered (significatively reduces compute time)", action="store_const", const="-k", default="")
        parser.add_argument("-ks","--kmerSize", help="kmer size for kmer filtering - default = 20", type=int, default=20)
	parser.add_argument("-pdf", "--pdfreports",help="generate PDF reports (can take long time)", action="store_const", const="-pdf", default="")
	
	parser.add_argument("-q","--samplesqc", action="store_const", const="-q", default="", help="samples quality control only")
	parser.add_argument("-s","--statsonly", action="store_const", const="-s", default="", help="do stats only")
	parser.add_argument("-r","--region", help="Region to analyse in sequences <min> <max>",nargs=2)
	parser.add_argument("-T","--MaxParallelsJobs", help="max number of parallels jobs to run - default= see config file", type=int)
        parser.add_argument("-e","--ErroRateThrld", help="Force error rate threshold - default= see config file", type=float)
	parser.add_argument("-S","--splitreads", help="split fastq files reads in half use this option if reads are interlaced - add [nosplit=True] parameter in reads information file for sample not need to be split", nargs='?', const=-1, type=int)
	
	parser.add_argument("-o","--outfolder", help="destination folder (create it if not exist)", required=True)
	parser.add_argument("-i","--readsinfo", help="Configuration file with reads informations if reads ares paired add [format=paired] parameter", required=True)
	parser.add_argument("-d","--refdatabase", help="reference database in fasta format (see documentation)", required=True)
	parser.add_argument("-x","--refexclude", help="simple text file contains for each line, references names to exclude for current analysis")
	
	parser.add_argument("--config", help="config file")
	args = parser.parse_args(ArgsVal)
	
	if args.config:
		config_file = args.config
	else:
		config_file = os.path.join(App_Folder,"config.yaml")
		
	config = yaml.load(open(config_file))

        #Force error rate threshold
        if args.ErroRateThrld:
                config['genotyp_def_ErrorRate']=args.ErroRateThrld

        config['refexclude']=[]

        if args.refexclude:
                if os.path.exists(args.refexclude):
                        tmp=[]
                        exfile=open(args.refexclude,'r')
                        for exref in exfile:
                                tmp.append(exref.strip())
                        exfile.close()
                        config['refexclude']=tmp
	
	currentDateTime = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
	
	config['log_file'] = "{appn}_Log_{dt}.txt".format(appn=AppName,dt=datetime.now().strftime("%d_%m_%Y-%H_%M_%S"))
	
	readsInfo_filename = args.readsinfo
	
	if check_FileExtensions(readsInfo_filename,['yaml','yml']):
		readsConfig = yaml.load(open(readsInfo_filename))
	else:
		readsConfig = load_readsList(readsInfo_filename)
        
	config['reads_dir'] = readsConfig['reads_dir']
	config['reads_files'], fastqParams = extract_paramsFromReadsConfig(readsConfig['reads_files'])

	keysList=['data_dir','samtools','bowtie2']
	
	add_AppFolderToYamlConfig(keysList)
	
	if args.MaxParallelsJobs:
		maxPjobs = args.MaxParallelsJobs
	else:
		maxPjobs = config['MaxParallelsJobs']
		
	stdout_print("Max parallels job set to {}".format(maxPjobs))
	
	utils.setMaxParallelJobsToConfig(maxPjobs)
	
	#Add out folder to config
	config['outfolder'] = args.outfolder
	
	dirFile_exist(config['outfolder'],1)
	dirFile_exist(config['data_dir'],1)
	
	if args.pdfreports:
		config['generate_pdf_reports']=True
	
	cmdline = " ".join(sys.argv)
	
	#log file write
	startText = "{appn} {appv} -- {cdate}".format(appn=AppName,appv=__version__,cdate=currentDateTime)
	starDateline = "*"*len(startText)
	addTextToLogFile("{}\n{}\n{}".format(starDateline,startText,starDateline),False)
	addTextToLogFile("user command: {}".format(cmdline))
	
	if config['samtools']!=None:
		samtools_path = os.path.join(App_Folder,config['samtools'])
		addTextToLogFile("USE samtools : {}".format(samtools_path))
	else:
		samtools_path = ''
	
	if config['bowtie2']!=None:
		bowtie2_path = os.path.join(App_Folder,config['bowtie2'])
		addTextToLogFile("USE bowtie2 : {}".format(bowtie2_path))
	else:
		bowtie2_path = ''

        if config['fastqc']!=None:
		fastqc_path = os.path.join(App_Folder,config['fastqc'])
		addTextToLogFile("USE fastqc : {}".format(fastqc_path))
	else:
		fastqc_path = ''
        
        if len(config['refexclude'])>0:
                addTextToLogFile("References to exluce : {}".format(",".join(config['refexclude'])))

	#Database fasta file
	IsRefInCatalog, config['ref_file'] = refDatabase_in_catalog(args.refdatabase,check_DatabaseFastaInputName(args.refdatabase))
	
	#check fastqFiles
	checkFQ = check_fastqFiles()

	if checkFQ['url'] and not args.kmerfilter:
                errormsg="ERROR: url should be use with -k (--kmerfilter) option"
                exitProg(errormsg)  

	load_HomoParaFromRef()
	
	# ----- MAIN WORKFLOW -----

	if args.samplesqc:
		addTextToLogFile("Sample QC only")
		stdout_print("-- Sample QC only ---")
		#Quality control for samples (reads)
		stdout_print("Samples QC control...")
		samples_quality_cont()

	else:
	
		if not args.statsonly:
					
			if args.kmerfilter:
				stdout_print("Fastq filtering with kmerRefFilter")
				FilteredReads_Config = Run_kmerRefFilter(args.kmerSize)
				config['reads_dir'] = FilteredReads_Config['reads_dir']
				config['reads_files'] = FilteredReads_Config['reads_files']
				
			#should be after kmerfilter
			if args.splitreads:
				stdout_print("Split fastq reads")
				splitedReads_Config = split_reads(args.splitreads,fastqParams)
				config['reads_dir'] = splitedReads_Config['reads_dir']
				config['reads_files'] = splitedReads_Config['reads_files']

                        finalReadsDic={'reads_dir':config['reads_dir'],'reads_files':config['reads_files']}
                        create_YamlFromDict(os.path.join(config['outfolder'],config['finalReads_yaml']), finalReadsDic)

			#Create fasta with ref and ref table correspondances and quality control for references
			stdout_print("Split fasta...")
			split_fasta(IsRefInCatalog)
		
			#Quality control for samples (reads)
			stdout_print("Samples QC control...")
			samples_quality_cont()
			
			#Create indexes files for Bowtie2
			if not IsRefInCatalog or args.force:
				stdout_print("Bowtie2 index...")
				bowtie_index()
	
			#Paired-end alignement with Bowtie
			stdout_print("Bowtie2 align...")
			bowtie_align()
			
		else:
			msg = "-- Stats and plots only ---"
			stdout_print(msg)
			addTextToLogFile("\t{}".format(msg))
			
		addTextToLogFile("-- STATS --")
		
		#Sort and indexes bam file with samtools and plots
		stdout_print("Samtools stats...")
		AllBamDepthCoverage = samtools_SortIndexStats()
		
		if (config['generate_pdf_reports']):
			stdout_print("Plot Coverage...")
			addTextToLogFile("\tGenerate PDF reports")
			bamPlotCoverage(AllBamDepthCoverage)
			
	#delete utils param file
	utils.delete_paramFileName()
	
	totalTime = time.time() - startTime
	
        endMessage = "-- ENDED IN {} --".format(time.strftime('%Hh%Mm%Ss', time.gmtime(totalTime)))
        stdout_print(endMessage)
	addTextToLogFile(endMessage)
	
#----- END OF MAIN WORKFLOW -----
#=================================================================================================================================

#*********************
#* genotyp Functions *
#*********************
def split_reads(cutpos,fastqParams):

	addTextToLogFile(" -- split Reads -- ")
	
	dt=datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
	
	splited_reads = {}
	
	jblst = utils.Jobslist("kmerRefFilter for all fastq files")
	
	reads_dir = config['reads_dir']
	
	reads_files = config['reads_files']
	
	out_folder = config['outfolder']
	
	splited_dir = config['out_splited']
	
	log_folder = os.path.join(out_folder,config['log_folder'])
	
	output = os.path.join(out_folder,splited_dir)
	
	dirFile_exist(output,1)
	
	splited_reads['reads_files']={}
	
	splited_reads['reads_dir'] = os.path.join(os.getcwd(),output)
	
	targetList = []
	
	cutpos_str = ""

        extList=['FASTQ','FQ']

	
	if cutpos>0:
		cutpos_str = " -c {}".format(cutpos)
		
	for k,reads in reads_files.items():	
		tmp = []		
		tmp.append(reads[0])
                params=fastqParams[reads[0]]
                NoSplit=False
                if "nosplit" in params.keys():
                        NoSplit=True
		for fread in reads[1:]:
                        readFile = os.path.join(reads_dir,fread)
			if "_RSplit" not in fread and not NoSplit:
				cmd = "{pth}/splitFastq_v1/splitFastq{ct} -o {out} {fq} >>{log}".format(ct=cutpos_str, pth=tools_folder,out=output,fq=readFile,log=os.path.join(log_folder,"splitFastq_{dt}.txt".format(dt=dt)))
				target = os.path.join(output,"{}_R1_RSplit.fastq".format(remove_fileExt(fread, extList)))
				targetList.append(target)
				
				tmp.append("{}_R1_RSplit.fastq".format(remove_fileExt(fread, extList)))
				tmp.append("{}_R2_RSplit.fastq".format(remove_fileExt(fread, extList)))
				
				if not os.path.exists(target):
					jblst.add_a_job(cmd,"splitFastq {}".format(fread),target)
                        elif NoSplit:
                                stdout_print("[{}] nosplit parameter for {}".format(reads[0],fread),True)
                                nosplitFname = "{}_RSplit.fastq".format(remove_fileExt(fread, extList))
                                slinkName = os.path.join(output,nosplitFname)
                                if os.path.exists(slinkName):
                                        cmd = "rm {}".format(slinkName)
                                        os.system(cmd)
                                cmd = "ln -s {src} {name}".format(src=readFile,name=slinkName)
                                os.system(cmd)
                                tmp.append(nosplitFname)
					
		splited_reads['reads_files'][k] = tmp
	
	stdout_print("{} fastq files to split".format(jblst.get_SizeOfJoblist()),True)
	
	dbfile = os.path.basename(config['ref_file'])
	dbname = dbfile[:dbfile.rfind('.')]
	create_YamlFromDict(os.path.join(output,"{}_{}".format(dbname,config['out_splitedReads_yaml'])),splited_reads)
	
	if jblst.get_SizeOfJoblist()>0:
		utils.trun(args, jblst)
	else:
		stdout_print("All files already splited",True)
	
	for f in targetList:
		if not os.path.exists(f):
			exitProg("ERROR \"{}\" NOT FOUND. File was not generated during splitFastq. See log file".format(f))
	return splited_reads

def remove_fileExt(filename,fileExtList):
        if is_url(filename):
                filename = os.path.basename(filename)
        for ext in fileExtList:
                if filename.upper().rfind(".{}".format(ext))>0:
                        return filename[0:filename.upper().rfind(".{}".format(ext))]
        return filename

def Run_kmerRefFilter(kmerSize):
	
	addTextToLogFile(" -- kmerRefFilter --")
	
	dt=datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
	
	filtered_reads = {}
	
	jblst = utils.Jobslist("kmerRefFilter for all fastq files")
	
	reads_dir = config['reads_dir']
	
	reads_files = config['reads_files']
	
	out_folder = config['outfolder']
	
	filtered_dir = config['out_filtered']
	
	log_folder = os.path.join(out_folder,config['log_folder'])
	
	output = os.path.join(out_folder,filtered_dir)
	
	kmerFilter_Yaml_File = os.path.join(output,config['kmerFilter_Yaml'])
	
	dirFile_exist(output,1)
	
	filtered_reads['reads_files']={}
	
	filtered_reads['reads_dir'] = os.path.join(os.getcwd(),output)
	
	targetList = []

	NbrFilesToFilter = 0

        extList=['FASTQ','FQ']

	for k,reads in reads_files.items():	
		tmp = []
		idName = reads[0]
                freadList = reads[1:]
		tmp.append(idName)
                parms = fastqParams[idName]
                
                if 'format' in parms.keys() and parms['format']=='paired' and len(freadList)%2==0:
                        for i in range(0,len(freadList),2):
                                fwdread = freadList[i]
                                revread = freadList[i+1]
                                if "_filtered" not in fwdread and "_filtered" not in revread:
                                        stdout_print("paired filtering for indiv: {} -- fwd:{} rev:{}".format(idName,fwdread,revread),True)
                                        if is_url(fwdread) and is_url(revread):
                                                fwdreadFile = fwdread
                                                revreadFile = revread
                                                gzipParam=""
                                                if check_FileExtensions(fwdreadFile,['gz']) and check_FileExtensions(revreadFile,['gz']):
                                                        gzipParam="-ugzip "

                                                cmd="{pth}kmerRefFilter.py -y -k {ksize} -u1 {FwdFq} -u2 {RevFq} -r {fastaref} {gzipval}-o {outdir} 2>>{log} 1>>{yamlKmer}".format(gzipval=gzipParam,pth=os.path.join(App_Folder,PyLibPath),FwdFq=fwdreadFile,RevFq=revreadFile,outdir=output,fastaref=config['ref_file'],log=os.path.join(log_folder,"kmerRefFilter_log_{dt}.txt".format(dt=dt)),yamlKmer=kmerFilter_Yaml_File,ksize=kmerSize)
                                        else:
                                                fwdreadFile = os.path.join(reads_dir,fwdread)
                                                revreadFile = os.path.join(reads_dir,revread)
                                                cmd="{pth}kmerRefFilter.py -y -k {ksize} -1 {FwdFq} -2 {RevFq} -r {fastaref} -o {outdir} 2>>{log} 1>>{yamlKmer}".format(pth=os.path.join(App_Folder,PyLibPath),FwdFq=fwdreadFile,RevFq=revreadFile,outdir=output,fastaref=config['ref_file'],log=os.path.join(log_folder,"kmerRefFilter_log_{dt}.txt".format(dt=dt)),yamlKmer=kmerFilter_Yaml_File,ksize=kmerSize)

                                        target = os.path.join(output,"{}FWD_filtered.fastq".format(remove_fileExt(fwdread, extList)))
                                        targetList.append(target)
                                        tmp.append("{}FWD_filtered.fastq".format(remove_fileExt(fwdread, extList)))
                                        tmp.append("{}REV_filtered.fastq".format(remove_fileExt(revread, extList)))
                                        NbrFilesToFilter += 2
                                        if not os.path.exists(target):
                                                jblst.add_a_job(cmd,"kmerRefFilter paired {} {}".format(fwdread,revread),target)
                else:

                        for fread in freadList:
                                if "_filtered" not in fread:
                                        if is_url(fread):
                                                readFile = fread
                                                gzipParam=""
                                                if check_FileExtensions(readFile,['gz']):
                                                        gzipParam="-ugzip "
                                                cmd="{pth}kmerRefFilter.py -y -k {ksize} -u {fq} -r {fastaref} {gzipval}-o {outdir} 2>>{log} 1>>{yamlKmer}".format(gzipval=gzipParam,pth=os.path.join(App_Folder,PyLibPath),fq=readFile,outdir=output,fastaref=config['ref_file'],log=os.path.join(log_folder,"kmerRefFilter_log_{dt}.txt".format(dt=dt)),yamlKmer=kmerFilter_Yaml_File,ksize=kmerSize)
                                        else:
                                                readFile = os.path.join(reads_dir,fread)
                                                cmd="{pth}kmerRefFilter.py -y -k {ksize} -f {fq} -r {fastaref} -o {outdir} 2>>{log} 1>>{yamlKmer}".format(pth=os.path.join(App_Folder,PyLibPath),fq=readFile,outdir=output,fastaref=config['ref_file'],log=os.path.join(log_folder,"kmerRefFilter_log_{dt}.txt".format(dt=dt)),yamlKmer=kmerFilter_Yaml_File,ksize=kmerSize)
                                        target = os.path.join(output,"{}_filtered.fastq".format(remove_fileExt(fread, extList)))
                                        targetList.append(target)
                                        tmp.append("{}_filtered.fastq".format(remove_fileExt(fread, extList)))
                                        NbrFilesToFilter += 1
                                        if not os.path.exists(target):
                                                jblst.add_a_job(cmd,"kmerRefFilter {}".format(fread),target)
					
		filtered_reads['reads_files'][k] = tmp

	#stdout_print("{} fastq files to filter".format(jblst.get_SizeOfJoblist()),True)
        stdout_print("{} fastq files to filter".format(NbrFilesToFilter),True)
	
	dbfile = os.path.basename(config['ref_file'])
	dbname = dbfile[:dbfile.rfind('.')]
	create_YamlFromDict(os.path.join(output,"{}_{}".format(dbname,config['out_reads_yaml'])),filtered_reads)
	
	if jblst.get_SizeOfJoblist()>0:
	
		if args.tinyverbose:
			t2 = threading.Thread(target=progress_FastqFilter, args=(jblst.get_SizeOfJoblist(),))
			t2.start()
		
		utils.trun(args, jblst)
	
		if args.tinyverbose:
			Wait_ProgressBarFinished(t2)
			
	else:
		stdout_print("All files already filtered",True)
	
	for f in targetList:
		if not os.path.exists(f):
			exitProg("ERROR \"{}\" NOT FOUND. File was not generated during filtering. See kmerRefFilter logs".format(f))
	
	kmerFilter_Yaml = yaml.load(open(kmerFilter_Yaml_File))
	
	filtered_reads['reads_Nbr'] = {}
	
	for k,v in kmerFilter_Yaml.items():
		filtered_reads['reads_Nbr'][k] = v['reads Nbr']
	
	return filtered_reads
	
def check_fastqFiles():

	reads_files = config['reads_files']
	reads_dir = config['reads_dir']
	
        url=False

	for read in reads_files.values():
		for fq in read[1:]:
                        if is_url(fq):
                                addTextToLogFile("{} is URL".format(fq))
                                url=True
                                #if not check_url_exists(fq):
                                #        errormsg="ERROR: \"{}\" URL not exists".format(fq)
                                #        exitProg(errormsg)  
                        else:
                                if not os.path.exists(os.path.join(reads_dir,fq)):
                                        errormsg="ERROR: \"{}\" not found".format(os.path.join(reads_dir,fq))
                                        exitProg(errormsg)
				
			if args.kmerfilter:
			
				if not check_FileExtensions(fq,['fastq','gz','bz2']):
					errormsg="ERROR: \"{}\" file extension should be .fastq (gz or bz2 with -k option)".format(fq)
					exitProg(errormsg)
			else:
			
				if not check_FileExtensions(fq,['fastq']):
					errormsg="ERROR: \"{}\" file extension should be .fastq".format(fq)
					exitProg(errormsg)
        return {'url':url}	
	
def stdout_print(txt, printInLog=False):
	if args.tinyverbose:
		print txt
	if printInLog:
		addTextToLogFile(txt)

def check_DatabaseFastaInputName(databaseName):

	if check_FileExtensions(databaseName,['fasta','fa']):
		return 'fasta'
	
	if databaseName.count('/')>0:
		errormsg="ERROR: Database entry should be a fasta file (.fasta, .fa) OR a database name only"
		exitProg(errormsg)
		
	return 'DB'

def check_readsInfoPaired():
	readsList = config['reads_files']
	for read in readsList.values():
		if (len(read)!=3): return False
		
	return True

def shared_reads(picklestats_list):

	out_folder = config['outfolder']
	all_results_dir = os.path.join(out_folder,config['all_results_dir'])
	dir_results = os.path.join(out_folder,config['bwt_dirResult'])
	
	sharedReads = {}
	
	for sample, SValues in picklestats_list.items():
		
		sharedReads[sample] = {}
		
		for refName, StatVal in SValues.items():
			if StatVal['error rate']<config['sharedReads_errorRatelim'] and StatVal['reads mapped']>0:
				Sample_path = os.path.join(dir_results,sample)
				aln_name = "Aln_{}_VS_{}".format(sample,refName)
				aln_path = os.path.join(Sample_path,aln_name)
				sorted_bam_file = os.path.join(aln_path, "{}_SORTED.bam".format(aln_name))
				samfile = pysam.AlignmentFile(sorted_bam_file, "rb")
				if args.region:
					fetched_SamFile = samfile.fetch(refName, int(args.region[0]), int(args.region[1]))
				else:
					fetched_SamFile = samfile.fetch()
				
				tmp = []
				for read in fetched_SamFile:
					hashval = read.query_name+read.seq+read.qual
					tmp.append(hashval)
					
				sharedReads[sample][refName] = tmp
				
				samfile.close()
				
		rslt = {}
		
		for sample, val in sharedReads.items():
			rslt[sample] = {}
			i=0
			for k1,v1 in val.items():
				rslt[sample][k1] = {}
				for k2,v2 in val.items()[i:]:
					rslt[sample][k1][k2] = int(len(set(v1) & set(v2)))
			i += 1
			
		generateXlsSharedReads(rslt, os.path.join(all_results_dir,config['xls_shared_reads']), picklestats_list)

def generateXlsSharedReads(ShReadsdict, xls_save_file, picklestats_list):
	
	#output xls file generation
	
	style = xlwt.easyxf('alignment: horizontal center, vertical center;')
	styleYellow = xlwt.easyxf('alignment: horizontal center, vertical center;pattern: pattern solid, fore_colour yellow')
	styleRed = xlwt.easyxf('alignment: horizontal center, vertical center;pattern: pattern solid, fore_colour red')
	styleHeaderTop = xlwt.easyxf('align: rotation 90, horizontal center, vertical center;')
	
	workbook = xlwt.Workbook()
	
	for k1,sample in sorted(ShReadsdict.items(),key=itemgetter(0)):
		if len(sample.keys())>0:
			if len(k1)>31:
				stdout_print("Shared reads: sheet name lenght '{}' to long truncate to '{}'".format(k1,k1[:31]))
			sheet = workbook.add_sheet(k1[:31])
			row=0
			col=1
			for h in sample.keys():
				sheet.write(row,col,"({ln}) {name}".format(name=h,ln=getRealRefLen(picklestats_list[k1][h]['RefLen'])),styleHeaderTop)
				sheet.col(col).width=1000
				col += 1
		
			sheet.row(0).height = max([len(v) for v in sample.keys()])*150
			sheet.col(0).width = max([len(v) for v in sample.keys()])*340
		
			row = 1
			col = 0
			i=0
			for ref, statRef in sample.items():
				sheet.write(row,0,"{name} ({ln})".format(name=ref,ln=getRealRefLen(picklestats_list[k1][ref]['RefLen'])),style)
			
				col =1
			
				for h in statRef.keys():
					if h==ref:
						cstyle = styleYellow
					else:
						cstyle = style
						if statRef[h]==statRef[ref]:
							cstyle = styleRed
					
					if statRef[h]>0:
						ccontent = statRef[h]
					else:
						ccontent = ''
					sheet.write(row,col,ccontent,cstyle)
					col += 1
				row += 1
				i+= 1
		
			region = ''
		
			if args.region:
				region = "from {} to {}\n".format(args.region[0],args.region[1])
			
			sheet.write(0,0,"{}alignments with error rate < {}".format(region, config['sharedReads_errorRatelim']),style)
		
	workbook.save(xls_save_file)

def getRealRefLen(RefLen):
	if args.region:
		if RefLen > int(args.region[1]):
			return int(args.region[1]) - int(args.region[0])
		else:
			return RefLen - int(args.region[0])
	else:
	
		return RefLen

def ScatterPlot_ERate_MeanCov(out_pdf,picklestats_list):
	
	colormap = cm.plasma
	
	region =''
	
	if args.region:
		region = "(from {rmin} To {rmax} bp)".format(rmin=args.region[0],rmax=args.region[1])
		
	with PdfPages(out_pdf) as pdf:
	
		for sample,ref in picklestats_list.items():
			MeanCov = [v['mean cover'] for i,v in ref.items()]
			ErrRate = [v['error rate'] for i,v in ref.items()]
			NormMed = [v['Norm med'] for i,v in ref.items()]

                        #avoid division by zero
			if max(NormMed)>0 and max(MeanCov)>0:

                                T = np.arctan2([v/max(ErrRate)-0.5 for v in ErrRate],[v/max(NormMed)-0.5 for v in NormMed])
                                
                                norm = Normalize()
                                norm.autoscale(T)
                                
                                fig = plt.figure()
                                
                                plot = plt.scatter(NormMed, NormMed, c = NormMed, cmap = 'plasma')
                                plt.clf()
                                
                                ax = fig.add_subplot(111)
                                
                                ax.set_ylim([-0.05,max(MeanCov)])
                                
                                plt.scatter(ErrRate,MeanCov, color=colormap(norm(T)))
                                plt.title("Sample {title} {region}".format(title=sample,region=region))
                                plt.xlabel('Error rate',fontsize=12)
                                plt.ylabel('Mean coverage',fontsize=12)
                                
                                cbar = plt.colorbar(plot)
                                cbar.set_label('Normalize Median')
                                pdf.savefig(plt.gcf())
                                plt.close('all')

def ScatterPlot_ERate_RegionCov(out_pdf,picklestats_list):
	
	region =''
        
        pointColor={True:{True:'blue',False:'blue'},False:{True:'red',False:'black'}}
	
	if args.region:
		region = "(from {rmin} To {rmax} bp)".format(rmin=args.region[0],rmax=args.region[1])
		
	with PdfPages(out_pdf) as pdf:
	
		for sample,ref in picklestats_list.items():
			RegionCov = [v['Region Cov'] for i,v in ref.items() if v['mean cover']>0]
			ErrRate = [v['error rate'] for i,v in ref.items() if v['mean cover']>0]
			Colors = [pointColor[v['IsParalog']][v['IsPositiv']] for i,v in ref.items() if v['mean cover']>0]

                        plt.clf()                                
                        fig, ax = plt.subplots()
                        ax.set_ylim([0,max(RegionCov)*(1+10.0/100.0)])
                        ax.set_xlim([0,max(ErrRate)*(1+10.0/100.0)])

                        plt.scatter(ErrRate,RegionCov,color=Colors,marker='x')
                        ax.axvline(x=config['genotyp_def_ErrorRate'],color='red', linewidth=0.5)
                        plt.title("Sample {title} {region}".format(title=sample,region=region))
                        plt.xlabel('Error rate',fontsize=12)
                        plt.ylabel('Region coverage',fontsize=12)
                        red_patch = mpatches.Patch(color='red', label='positiv samples')
                        blue_patch = mpatches.Patch(color='blue', label='paralogs')
                        black_patch = mpatches.Patch(color='black', label='negativ samples')
                        plt.legend(handles=[red_patch,black_patch,blue_patch],loc=1,framealpha=0.5,fontsize=6)
                        pdf.savefig(plt.gcf())
                        plt.close('all')


def ScatterPlot_HomoPara(out_pdf,picklestats_list):
	
	colormap = cm.hsv
	
	region =''
	
	if args.region:
		region = "(from {rmin} To {rmax} bp)".format(rmin=args.region[0],rmax=args.region[1])
		
	with PdfPages(out_pdf) as pdf:
	
		for sample,ref in picklestats_list.items():
			MeanCov = [v['mean cover'] for i,v in ref.items()]
			ErrRate = [v['error rate'] for i,v in ref.items()]
			Paralog = [v['IsParalog'] for i,v in ref.items()]
			GrpID = list(set([v['Group ID'] for i,v in ref.items()]))
			
			ColorVal = [get_HomoPara_Parameters(i)['grpColor'] for i,v in ref.items()]
			GrpNbr = [get_HomoPara_Parameters(i)['gprId'] for i,v in ref.items()]
			HaploId = [get_HomoPara_Parameters(i)['HaploId'] for i,v in ref.items()]
			grpRef = [get_HomoPara_Parameters(i)['grpRef'] for i,v in ref.items()]
			
			Labels = []
			GrpColors = {GrpID[i]:colormap(int(float(i)/len(GrpID)*255)) for i in range(0,len(GrpID),1)}

			fig = plt.figure()
			
			ax = fig.add_subplot(111)
			
			ax.set_ylim([-0.05,max(MeanCov)])
			
			for ErrRt, MeanC, GrpN, HapId, IsParalog, grpRf in zip(ErrRate, MeanCov, GrpNbr, HaploId, Paralog, grpRef):
				
				ColorV = GrpColors[grpRf]
				
				if GrpN!='':
					tmp = "Group {}".format(GrpRomanNumber(GrpN))
					marker = "o"
					backcolor = ColorV
				else:
					if IsParalog:
						tmp = "Paralog"
						marker = "s"
						ColorV = 'black'
						backcolor = 'black'
					else:
						tmp = "No group"
						marker = "o"
						ColorV = 'black'
						backcolor = 'None'
				
				if MeanC>0:
					alpha = 1
				else:
					alpha = 0.1
				
				Labels.append(tmp)
				
				plt.scatter(ErrRt,MeanC, color=backcolor, marker=marker, edgecolors=ColorV, label=tmp, alpha=alpha)
				
			LabelsList = ['Paralog','Homolog','No information']
			pltLgd = []
			pltLgd.append(Line2D(range(1), range(1), color="white", marker='s', markersize=8, markerfacecolor = 'black', markeredgecolor = 'black'))
			pltLgd.append(Line2D(range(1), range(1), color="white", marker='o', markersize=8, markerfacecolor = 'red', markeredgecolor = 'red'))
			pltLgd.append(Line2D(range(1), range(1), color="white", marker='o', markersize=8, markerfacecolor = 'white', markeredgecolor = 'black'))
			
			plt.title("Sample {title} {region}".format(title=sample,region=region))
			plt.xlabel('Error rate',fontsize=12)
			plt.ylabel('Mean coverage',fontsize=12)
			plt.legend(pltLgd, LabelsList, scatterpoints=1, loc=1, fontsize=8)
			pdf.savefig(plt.gcf())
			plt.close('all')
			
def get_HomoPara_Parameters(refName, HomoPara_grp = None):

	if HomoPara_grp==None:
		HomoPara_grp = HomoParaFromRef

	for k,v in HomoPara_grp.items():
		v['RefName'] = refName
		keylist = v.keys()
		if refName.upper()==k.upper():
			IsKeyPresentInDic('HaploId',v,'')
			IsKeyPresentInDic('grpColor',v,'#000000')
			IsKeyPresentInDic('gprId',v,'')
			IsKeyPresentInDic('specie',v,'')
			IsKeyPresentInDic('grpRef',v,'')
			return v

def Load_Paralogs():
	list_paralogs = [ k for k,v in HomoParaFromRef.items() if 'Paralog' in v]
	return list_paralogs

def bamDepthCoverage(bamfile, aln_path,aln_name, ref_name):
	
	cmd = "{path}samtools view -h -b {bam} | {path}samtools depth -a - | cut -f 2-".format(path=samtools_path,bam=bamfile)

	depthrstl = utils.run(args,"samtools depth", '', cmd, True)

	xdepth = []
	ydepth = []

	for line in depthrstl[0].split("\n"):
		if len(line.split("\t"))>1:
			xdepth.append(int(line.split("\t")[0]))
			ydepth.append(int(line.split("\t")[1]))
	
	fxdepth = xdepth
	fydepth = ydepth

	if args.region:
		rmin = int(args.region[0])
		rmax = int(args.region[1])
		xdepth = xdepth[rmin:rmax+1]
		xdepth[:] = [int(x) - rmin for x in xdepth]
		ydepth = ydepth[rmin:rmax+1]

	if (len(ydepth)>0 and sum(ydepth)>0):
	
		mean = sum(ydepth)/float(len(ydepth))
		gmean = sum(fydepth)/float(len(fydepth))
		pSeqCov = sum([1 for v in fydepth if v>0])/float(len(fydepth))
		
		revCumSum = []
		for i in range(0,int(max(ydepth))):
			pos = 0
			tmp = 0
			for val in ydepth:
				if val >= i:
					tmp += 1
				pos += 1
			revCumSum.append(tmp/float(xdepth[-1]))
		
		covfreq = getFrequence(ydepth)
		
		maxdepth = int(max(covfreq[0]))

		med = np.median(ydepth)

		q75, q25 = np.percentile(ydepth, [75 ,25])

		iqr = q75 - q25

		riqr=iqr/maxdepth
		rmed=med/maxdepth
	
	else:
		mean = 0
		gmean = 0
		pSeqCov = 0
		xdepth = []
		xdepth = []
		revCumSum = []
		fxdepth = []
		fydepth = []
		covfreq = []
		med = 0
		q25 = 0
		q75 = 0
		maxdepth = 0
		iqr = 0
		riqr = 0
		rmed = 0
	
	return {'mean':mean,'gmean':gmean,'pSeqCov':pSeqCov,'xdepth':xdepth,'ydepth':ydepth,'revCumSum':revCumSum,'aln_name':aln_name,'ref_name':ref_name,'fxdepth':fxdepth,'fydepth':fydepth,'covfreq':covfreq,'median':med,'q25':q25,'q75':q75,'maxdepth':maxdepth,'iqr':iqr, 'riqr':riqr, 'rmed':rmed}

def plotDephCoverage(depthCov, pp, title):
	
	depthCov = [v for v in depthCov if len(v['fydepth'])>0]
	
	linesNbr = config['MaxgraphPerFigure']
	colNbr = 3
	gpos = 1
	
	if linesNbr>0:
		plt.figure(figsize=(20,5*(linesNbr+1)))
	
		plt.subplot(linesNbr,colNbr,2)
		plt.title("{}\n".format(title),fontsize=25)
	
		for value in depthCov:
			mean = value['mean']
			gmean = value['gmean']
			xdepth = value['xdepth']
			ydepth = value['ydepth']
			revCumSum = value['revCumSum']
			aln_name = value['aln_name']
			ref_name = value['ref_name']
			err_rate = value['error rate']
			fxdepth = value['fxdepth']
			fydepth = value['fydepth']
			freq = value['covfreq']
			med = value['median']
			q75 = value['q75']
			q25 = value['q25']
			maxdepth = value['maxdepth']
			iqr = value['iqr']
			riqr = value['riqr']
			rmed = value['rmed']
                        IsParalog = value['IsParalog']
                        IsPositiv = value['IsPositiv']
		
                        if IsParalog:
                                textColor='orange'
                        elif not IsParalog and IsPositiv:
                                textColor='green'
                        else:
                                textColor='black'

			maxfdepth = int(max(fydepth))
		
			axes = plt.subplot(linesNbr,colNbr,gpos)
			axes.spines['right'].set_visible(False)
			axes.spines['top'].set_visible(False)
			axes.set_ylim([-0.5,maxfdepth+maxfdepth*0.2])
			plt.ylabel('depth',fontsize=12)
			plt.xlabel('Ref sequence in bp',fontsize=10)
			plt.text(0,maxfdepth+maxfdepth*0.1, "Ref.: {title}\nerror rate = {er}".format(title=ref_name, er=round(err_rate,4)),fontsize=9,fontweight='bold',color = textColor,backgroundcolor='white')
			plt.plot(fxdepth,fydepth,color = '#ff7e7e')
			if args.region:
				rmin = int(args.region[0])
				rmax = int(args.region[1])
				if rmax>len(fxdepth): rmax =len(fxdepth)
				plt.axvspan(rmin, rmax, alpha=0.2, ymax=0.85, color='#fff300')
				plt.plot([rmin,rmax],[mean,mean],linestyle='--',color='#ff00f7',label="region mean={mean}".format(mean=round(mean,2)))
			plt.plot([0,fxdepth[-1]],[gmean,gmean],'g--',label="mean={gmean}".format(gmean=round(gmean,2)))
			plt.legend(loc=1,fontsize=9)
	
			axes = plt.subplot(linesNbr,colNbr,gpos+1)
			axes = plt.gca()
			axes.set_ylim([-0.05,1.05])
			axes.set_xlim([0,max(ydepth)])
			axes.spines['right'].set_visible(False)
			axes.spines['top'].set_visible(False)
			plt.ylabel('fraction of bases\nsampled >= coverage',fontsize=8)
			plt.xlabel('coverage (#reads per bp)',fontsize=10)
			plt.plot(revCumSum,color = 'blue')
		
			ybar = int(max(freq[1]))/2
			axes = plt.subplot(linesNbr,colNbr,gpos+2)
			axes = plt.gca()
			axes.set_ylim([0,int(max(freq[1]))])
			axes.set_xlim([0,maxdepth])
			axes.spines['right'].set_visible(False)
			axes.spines['top'].set_visible(False)
			plt.xlabel('coverage (#reads per bp)',fontsize=8)
			plt.ylabel('Number of Bases',fontsize=10)
			plt.bar(freq[0],freq[1],color = '#e8e8e8')
			plt.errorbar(med, ybar, xerr=np.array([[med-q25 ,q75-med]]).T, capsize=5, fmt='ko',label="Norm IQR: {riqr}\nNorm med: {rmed}\nmed: {med} IQR: {iqr}".format(med=med,iqr=iqr,riqr=round(riqr,5),rmed=round(rmed,5)))
			plt.legend(loc=1,fontsize=9)
		
			gpos += colNbr
		
		pp.savefig(plt.gcf())
		plt.close('all')
	
def bamPlotCoverage(AllBamDepthCoverage):

	graphPerFigure = config['MaxgraphPerFigure']
	
	out_folder = config['outfolder']
	dir_pdfout = os.path.join(out_folder,config['all_results_dir'])
	dir_results = os.path.join(out_folder,config['bwt_dirResult'])
	all_results_dir = os.path.join(out_folder,config['all_results_dir'])
	picklestats_file = config['picklestats_file']
	
	bwt_yam_reflist = os.path.join(os.path.join(out_folder,config['bwt_IndexdirResult']), config['bwt_yam_reflist'])
	
	picklestats_list = load_DicFromPickleFile(os.path.join(all_results_dir,picklestats_file))
	
	plot_params = config['plot_selection_param']
	
	str_plot_params = ' - '.join([str(val[0])+str(val[1])+str(val[2]) for val in plot_params])
	
	if config['ERate_MeanCov_PDF']:
                ScatterPlot_ERate_MeanCov(os.path.join(dir_pdfout,"ErrMeanCov_plot.pdf"),picklestats_list)
	
        if config['Group_ErrMeanCov_PDF']:
                ScatterPlot_HomoPara(os.path.join(dir_pdfout,"Group_ErrMeanCov_plot.pdf"),picklestats_list)

        if config['ErrRegionCov_PDF']:
	        ScatterPlot_ERate_RegionCov(os.path.join(dir_pdfout,"ErrRegionCov_plot.pdf"),picklestats_list)

	region =''
	
	if args.region:
		region = "(from {rmin} To {rmax} bp)".format(rmin=args.region[0],rmax=args.region[1])
	
	for read,references in picklestats_list.items():
		depthCov = []
		
		with PdfPages(os.path.join(dir_pdfout,"Coverages_sample_{}.pdf".format(read))) as pp:
			page = 1
			
			for ref,stats in references.items():
								
				title = "Coverages for sample {sp} {region} - page {p}\n{params}".format(sp=read,p=page, region=region,params=str_plot_params)
				
				if checkSelection_param(stats,'plot_selection_param') and stats['max depth']>0:
					Sample_path = os.path.join(dir_results,read)
				
					aln_name = "Aln_{}_VS_{}".format(read,ref)
					aln_path = os.path.join(Sample_path,aln_name)
					sorted_bam_file = os.path.join(aln_path, "{}_SORTED.bam".format(aln_name))
					target = os.path.join(aln_path, "{}_coverage.png".format(aln_name))
				
					if os.path.exists(sorted_bam_file):
						depthCov.append(dict(AllBamDepthCoverage[aln_name].items() + stats.items()))
					
					if len(depthCov)%graphPerFigure==0:
						plotDephCoverage(depthCov,pp,title)
						depthCov = []
						page += 1
				
				
			if len(depthCov)>0:
				plotDephCoverage(depthCov,pp,title)

def samtools_SortIndexStats():
	
	DB_dir = config['DB_dir']
	
	out_folder = config['outfolder']
	
	DB_Path = os.path.join(config['data_dir'],DB_dir)
	
	dir_results = os.path.join(out_folder,config['bwt_dirResult'])
	
	dir_allResults = os.path.join(out_folder,config['all_results_dir'])
	
	bwt_yam_reflist = os.path.join(os.path.join(DB_Path,config['bwt_IndexdirResult']), config['bwt_yam_reflist'])
	
	References_list = yaml.load(open(bwt_yam_reflist))
	
	ReferencesQC = yaml.load(open(os.path.join(DB_Path,config['yml_qc_ref'])))
	
	ReadsCounts = yaml.load(open(os.path.join(dir_allResults,config['yamlReads_file'])))
	
	picklestats_file = config['picklestats_file']
	
	picklestats_params = config['picklestats_params']
	
	picklestat_path = os.path.join(dir_allResults,picklestats_file)
	
	xlsstat_filename = str(config['outfolder']).rstrip('/')

        ErrCovDensityPlot_path = os.path.join(dir_allResults,"ErrCovDepth_plot.pdf")
	
	rfind_fn = xlsstat_filename.rfind('/')
	
	if rfind_fn>=0:
		xlsstat_filename = xlsstat_filename[rfind_fn+1:]
	
	xlsstat_path = os.path.join(dir_allResults,"Genotyp_Stats_{}.xls".format(xlsstat_filename))

	outTXT = os.path.join(dir_allResults,config['putative_alleles'])

	remove_File(picklestat_path)
	
	remove_File(xlsstat_path)
	
	dirFile_exist(dir_allResults,1)
		
	stats_rslt = {}
	AllBamDepthCoverage = {}
	nbaln = 0
	
	errorStderrFiles = get_errorStderrFiles()
	
	for reads in config['reads_files']:
	
		read_file = config['reads_files'][reads]
		Sample_path = os.path.join(dir_results,read_file[0])
		
		stats_rslt [read_file[0]] = {}
		headers = []
		
		readsNbr = sum(ReadsCounts[read_file[0]].values())
		
                references_Names = [v for v in References_list.values() if v not in config['refexclude']]

		for reference in references_Names:
		
			aln_name = "Aln_{}_VS_{}".format(read_file[0],reference)
			aln_path = os.path.join(Sample_path,aln_name)
			
			bam_file = os.path.join(aln_path, "{}.bam".format(aln_name))
			stderr_file = os.path.join(aln_path, "stderr_{}".format(aln_name))
			sorted_bam_file = os.path.join(aln_path, "{}_SORTED.bam".format(aln_name))
					
			if not Is_stderrFileInError(stderr_file):
			
				if os.path.exists(bam_file):
					cmd = "{path}samtools sort {inbam} -o {outbam}".format(path=samtools_path,inbam=bam_file,outbam=sorted_bam_file)
					utils.run(args,"samtools sort",sorted_bam_file,cmd)
				
				if os.path.exists(sorted_bam_file):
				
					nbaln += 1
				
					target = os.path.join(aln_path, "{}_SORTED.bam.bai".format(aln_name))
					cmd = "{path}samtools index {inbam}".format(path=samtools_path,inbam=sorted_bam_file)
					utils.run(args,"samtools index",target,cmd)
				
					region = ''
				
					stats_rslt [read_file[0]][reference] = {}
				
					if args.region:
						region = "{ref}:{rmin}-{rmax} ".format(ref=reference,rmin=args.region[0],rmax=args.region[1])
						stats_rslt [read_file[0]][reference]['Region'] = "{rmin}-{rmax}".format(rmin=args.region[0],rmax=args.region[1])
				
					cmd = "{path}samtools stats {inbam} {region}| grep ^SN | cut -f2-3".format(path=samtools_path,inbam=sorted_bam_file,region=region)
					stdout, stderr = utils.run(args,"samtools stats",picklestat_path,cmd, True)
				
								
					stats_rslt [read_file[0]][reference]['RefLen'] = ReferencesQC[reference]['seq_len']
			
					tab_stat = stdout.split('\n')
					for val in tab_stat:
						tmp = val.split(':\t')
						if len(tmp)==2 and (tmp[0] in picklestats_params):
							headers.append(tmp[0])
							stats_rslt [read_file[0]][reference][tmp[0]] = castStrNbr_FloarOrInt(tmp[1])
						
					rslt_BDC = bamDepthCoverage(sorted_bam_file, aln_path, aln_name, reference)
				
					stats_rslt [read_file[0]][reference]['max depth'] = float(rslt_BDC['maxdepth'])
					stats_rslt [read_file[0]][reference]['mean cover'] = float(rslt_BDC['mean'])
					stats_rslt [read_file[0]][reference]['NormDepth'] = 0
					stats_rslt [read_file[0]][reference]['cov median'] = float(rslt_BDC['median'])
					stats_rslt [read_file[0]][reference]['cov IQR'] = float(rslt_BDC['iqr'])
					stats_rslt [read_file[0]][reference]['Norm IQR'] = float(rslt_BDC['riqr'])
					stats_rslt [read_file[0]][reference]['Norm med'] = float(rslt_BDC['rmed'])
					stats_rslt [read_file[0]][reference]['Region Cov'] = round(float(rslt_BDC['pSeqCov']),3)
                                        stats_rslt [read_file[0]][reference]['gscore'] = 0
                                        stats_rslt [read_file[0]][reference]['Warnings'] = []
					AllBamDepthCoverage[aln_name] = (rslt_BDC)
						
	stats_rslt = analyse_StatsResults(stats_rslt,ErrCovDensityPlot_path)

        generatePutativeAllelesTXTFile(outTXT,stats_rslt)

	create_PickleFromDict(picklestat_path,stats_rslt)
	
	headers.append('RefLen')
	headers.append('mean cover')
        headers.append('NormDepth')
	headers.append('cov median')
	headers.append('cov IQR')
	headers.append('Norm IQR')
	headers.append('Norm med')
	headers.append('max depth')
	headers.append('Region Cov')
	headers.append('gscore')
        headers.append('Warnings')
	
	headers = list(set(headers))
	
	#does config already contain headers	
        if ('result_stats_headers' in config):
                headers = config['result_stats_headers']
		
	if args.region:
		headers.append('Region')
	
	generateXlsStatsFromDict(stats_rslt, ReadsCounts, headers, xlsstat_path)
	
	shared_reads(stats_rslt)
	
	msg = "stats on {} alignements".format(str(nbaln))
	stdout_print(msg)
	addTextToLogFile("\t{}".format(msg))
	
	return AllBamDepthCoverage

#                 Score calculation
#------------------------------------------------------
def loi_Exponentielle(x,l):
        if(x>0):
                return l*math.exp(-l*x)
        else:
                return l

def loi_LogNormal(x,mu,sigma):
        if(x>0):
                p1 = 1/(x*sigma*math.sqrt(2*math.pi))
                p2 = math.exp(-((math.log(x)-mu)**2/(2*sigma**2)))
                return (p1*p2)
        else:
                return 0.0

def loi_Normal(x,mu,sigma):
        p1 = 1/(sigma*math.sqrt(2*math.pi))
        p2 = math.exp(-0.5*((x-mu)/sigma)**2)
        return p1*p2
        
def loi_ExpoFixed(x):
        l=config['lbda']
        return loi_Exponentielle(1-x,l)

def loi_LogNormalFixed(x):
        mu=config['mu']
        sigma=config['sigma']
        return loi_LogNormal(x,mu,sigma)

def loi_NormalFixed(x):
        mu=config['mu']
        sigma=config['sigma']
        return loi_Normal(x,mu,sigma)

def calculIntegrale(f,a,b):
        if callable(f):
                s=0
                stepVal=float(config['simpson_step'])
                for x in np.arange(a,b,stepVal):
                        s += (stepVal/6.0)*(f(x)+4*f((x+x+stepVal)/2.0)+f(x+stepVal))
                return s
        else:
                print "f is not a function"

def get_ErrProb(err):
        return calculIntegrale(loi_NormalFixed,err,1)
        #return calculIntegrale(loi_LogNormalFixed,err,1)

def get_CovProb(cov):
        return calculIntegrale(loi_ExpoFixed,0,cov)

def allele_prob(err,cov):
        errProb = get_ErrProb(err)
        covProb = get_CovProb(cov)
        allProb = eval(config['alleleProb_formula'])
        return {'alleleP':allProb,'errP':errProb,'covP':covProb}

def get_maxMeanCover(stats_rslt,sname):
        maxMeanCover = 0.00000000000001
        if sname in stats_rslt.keys():
               maxMeanCover = max([statRef['mean cover'] for statRef in stats_rslt[sname].values()])
        return maxMeanCover

def ErrCov_Density_plots(ErrDatas,CovDatas,IsPositiv,THRLD_range,maxerr=0):
        colorPos = {True:'red',False:'black'}
        colors = [colorPos[v] for v in IsPositiv]
        plt.clf()
        fig, ax1 = plt.subplots()
        ax1.xaxis.set_major_locator(ticker.AutoLocator())
        ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax1.set_axisbelow(True)
        ax1.xaxis.grid(True, linestyle=':',color='lightgray',which='minor')
        ax1.xaxis.grid(True, linestyle='-',color='lightgray',which='major')
        ax1.scatter(ErrDatas,CovDatas,marker=".",color=colors)
        err_THRLD_range = [v[0] for v in THRLD_range]
        cov_THRLD_range = [v[1] for v in THRLD_range]
        ax1.plot(err_THRLD_range,cov_THRLD_range,linestyle=':',color='orange')
        if maxerr>0:
                ax1.set_xlim([0,maxerr])
        ax1.set_ylim([0,1.0])
        ax1.text(max(err_THRLD_range),max(cov_THRLD_range),"genotyp score thresold: {}".format(config['genotyp_alleleProb_THRLD']),fontsize=6,color='orange')
        ax1.set_ylabel(r"Coverage", fontsize=12)
        ax1.set_xlabel("Error rate", fontsize=12)
        return fig

def ErrDepth__plot(ErrDatas,DepthDatas,IsPositiv,maxcov=0):
        colorPos = {True:'red',False:'black'}
        colors = [colorPos[v] for v in IsPositiv]
        plt.clf()
        fig, ax1 = plt.subplots()
        ax1.yaxis.set_major_locator(ticker.AutoLocator())
        ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax1.set_axisbelow(True)
        ax1.yaxis.grid(True, linestyle=':',color='lightgray',which='minor')
        ax1.yaxis.grid(True, linestyle='-',color='lightgray',which='major')
        ax1.scatter(DepthDatas,ErrDatas,marker=".",color=colors)
        if maxcov>0:
                ax1.set_ylim([0,maxcov])
        ax1.set_xlabel(r"Mean Depth", fontsize=12)
        ax1.set_ylabel("Error rate", fontsize=12)
        return fig

def LogalleleProb__plot(sortedLogAllelesProb,ProbThrld,title):
        plt.clf()
        fig, ax1 = plt.subplots()
        ax1.yaxis.set_major_locator(ticker.AutoLocator())
        ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax1.set_axisbelow(True)
        ax1.yaxis.grid(True, linestyle=':',color='lightgray',which='minor')
        ax1.yaxis.grid(True, linestyle='-',color='lightgray',which='major')
        ax1.set_title(title,fontsize=8)
        val = [v[1] for v in sortedLogAllelesProb]
        ax1.plot(val,marker=".")
        #ax1.set_yticks(np.arange(math.floor(min(val)),0,step=1))
        #ax1.tick_params(axis='y',labelsize=6)
        ax1.axhline(y=math.log10(ProbThrld),color='red', linewidth=1)
        ax1.set_xlabel(r"Alignments", fontsize=12)
        ax1.set_ylabel("log(Allele probability)", fontsize=12)
        return fig

def ErrCovDepth_plot3d(ErrDatas,CovDatas,Depth,IsPositiv):
    col={True:'red',False:'black'}
    colors=[col[v] for v in IsPositiv]
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("Error rate VS Coverage VS Depth (without paralogs)",fontsize=8)
    ax.scatter(Depth,CovDatas,ErrDatas,c=colors)
    ax.set_ylabel('Reference coverage')
    ax.set_xlabel('Mean depth')
    ax.set_zlabel('Error rate')
    return fig

def ErrCovScore_plot3d(ErrDatas,CovDatas,NormScore,IsPositiv):
    col={True:'red',False:'black'}
    colors=[col[v] for v in IsPositiv]
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("Error rate VS Coverage VS Genotyp score (without paralogs)",fontsize=8)
    ax.scatter(NormScore,CovDatas,ErrDatas,c=colors)
    ax.set_ylabel('Reference coverage')
    ax.set_xlabel('Genotyp score')
    ax.set_zlabel('Error rate')
    return fig

def ErrorRate_Density(ErrList):
        ErrRatesList=sorted(ErrList)
        density = gaussian_kde(ErrRatesList)
        xs=np.linspace(0,max(ErrRatesList),len(ErrRatesList))
        density._compute_covariance()       
        Adensity = np.array(density(xs))
        return {'xs':xs,'density':Adensity}

def get_errorratesCov(stats,List_Paralogs,onlyParalogs=False):
        ErrorRateCovList={'err':[],'cov':[],'IsPositiv':[],'depth':[],'gscore':[]}
        for indiv,v in stats.items():
                for ref,w in v.items():
                        if onlyParalogs:
                                b=ref in List_Paralogs
                        else:
                                b=ref not in List_Paralogs

                        if w['mean cover']>0 and b:
                                ErrorRateCovList['err'].append(w['error rate'])
                                ErrorRateCovList['cov'].append(w['Region Cov'])
                                ErrorRateCovList['IsPositiv'].append(w['IsPositiv'])
                                ErrorRateCovList['depth'].append(w['mean cover'])
                                ErrorRateCovList['gscore'].append(w['gscore'])
        return ErrorRateCovList

def analyse_StatsResults(stats_rslt,ErrCovDensityPlot_path):
	
	List_Paralogs = Load_Paralogs()
        allele_prob_THRLD = float(config['genotyp_alleleProb_THRLD'])
        LogAllelesProbList=[]
        LogParalogsAllelesProbList=[]

	for sname, values in stats_rslt.items():
		for ref, statRef in values.items():

                        if stats_rslt[sname][ref]['RefLen']<config['seq_min_len']:
                                stats_rslt[sname][ref]['Warnings'].append(config['Warn_refLen'])

                        HomoParaInfo = get_HomoPara_Parameters(ref, HomoParaFromRef)
			stats_rslt[sname][ref]['Group ID'] = str(HomoParaInfo['grpRef'])

			stats_rslt[sname][ref]['IsPositiv'] = False

			if ref in List_Paralogs:
				stats_rslt[sname][ref]['IsParalog'] = True
				stats_rslt[sname][ref]['Group ID'] = 'Paralog'
			else:
				stats_rslt[sname][ref]['IsParalog'] = False

                #NormDepth calculation
                parDepth = np.array([v['mean cover'] for v in values.values() if (v['IsParalog'] and v['error rate']<config['genotyp_def_ErrorRate'] and v['mean cover']>0)])
                if len(parDepth)>0:
                        parMedDepth = np.median(parDepth)

                parCov = np.array([v['Region Cov'] for v in values.values() if (v['IsParalog'] and v['error rate']<config['genotyp_def_ErrorRate'] and v['mean cover']>0)])
                if len(parCov)>0:
                        parMedCov = np.median(parCov)

                #Allele probability calculation
                for ref, statRef in values.items():
                        if config['normalized_cov']:
                                covValue = statRef['Region Cov']*(1.0/float(parMedCov))
                        else:
                                covValue = statRef['Region Cov']

                        AlleleProb = allele_prob(statRef['error rate'],covValue)
                        stats_rslt[sname][ref]['alleleProb'] = AlleleProb['alleleP']
                        if AlleleProb['alleleP']>0:
                                log10AlleleP = math.log10(AlleleProb['alleleP'])
                        else:
                                log10AlleleP = 0                                
                        stats_rslt[sname][ref]['LOGalleleProb'] = log10AlleleP
			stats_rslt[sname][ref]['gscore'] = AlleleProb['alleleP']
			stats_rslt[sname][ref]['errProb'] = AlleleProb['errP']
			stats_rslt[sname][ref]['covProb'] = AlleleProb['covP']

                        if len(parDepth)>0:
                                statRef['NormDepth']= statRef['mean cover']/parMedDepth
                        else:
                                statRef['NormDepth']= statRef['mean cover']

                LogAllelesProbList = LogAllelesProbList + [(v['alleleProb'],v['LOGalleleProb']) for v in values.values() if (v['alleleProb']>0  and not v['IsParalog'])]
                LogParalogsAllelesProbList = LogParalogsAllelesProbList + [(v['alleleProb'],v['LOGalleleProb']) for v in values.values() if (v['alleleProb']>0 and v['IsParalog'])]

                #define putatives alleles
                for ref, statRef in values.items():
                        if (stats_rslt[sname][ref]['alleleProb']>=config['genotyp_alleleProb_THRLD']):
                                stats_rslt[sname][ref]['IsPositiv'] = True
                        else:
                                stats_rslt[sname][ref]['IsPositiv'] = False

        THRLD_range = get_threshold_limits(0.1,0.5,config['genotyp_alleleProb_THRLD'])

        sortedLogAllelesProb = sorted(LogAllelesProbList,key=itemgetter(0),reverse=True)
        sortedLogParalogsAllelesProbList = sorted(LogParalogsAllelesProbList,key=itemgetter(0),reverse=True)

        ErrorRateCovList = get_errorratesCov(stats_rslt,List_Paralogs)
        ErrorRateParalogsCovList = get_errorratesCov(stats_rslt,List_Paralogs,True)
        ErrorRateDensityDatas = ErrorRate_Density(ErrorRateCovList['err'])

        with PdfPages(ErrCovDensityPlot_path) as pdf:
                pdf.savefig(ErrCov_Density_plots(ErrorRateCovList['err'],ErrorRateCovList['cov'],ErrorRateCovList['IsPositiv'],THRLD_range,0.35))
                pdf.savefig(ErrCov_Density_plots(ErrorRateParalogsCovList['err'],ErrorRateParalogsCovList['cov'],ErrorRateParalogsCovList['IsPositiv'],THRLD_range,0.35))
                pdf.savefig(ErrDepth__plot(ErrorRateCovList['err'],ErrorRateCovList['depth'],ErrorRateCovList['IsPositiv'],0.3))
                pdf.savefig(ErrDepth__plot(ErrorRateParalogsCovList['err'],ErrorRateParalogsCovList['depth'],ErrorRateParalogsCovList['IsPositiv'],0.3))
                pdf.savefig(ErrCovDepth_plot3d(ErrorRateCovList['err'],ErrorRateCovList['cov'],ErrorRateCovList['depth'],ErrorRateCovList['IsPositiv']))
                pdf.savefig(ErrCovScore_plot3d(ErrorRateCovList['err'],ErrorRateCovList['cov'],ErrorRateCovList['gscore'],ErrorRateCovList['IsPositiv']))
                pdf.savefig(LogalleleProb__plot(sortedLogAllelesProb,config['genotyp_alleleProb_THRLD'],"Sorted Allele probability for references alleles (exclude paralogs)"))
                pdf.savefig(LogalleleProb__plot(sortedLogParalogsAllelesProbList,config['genotyp_alleleProb_THRLD'],"Sorted Allele probability for paralogs only"))
        return stats_rslt

def get_threshold_limits(maxErr,minCov,alleleProbTHRLD):
        THRLD_range=[]
        tmp=0
        for err in np.arange(0,maxErr,step=0.005):
                for cov in np.arange(1,minCov,step=-0.01):
                        alleleP = allele_prob(err,cov)['alleleP']
                        if alleleP<alleleProbTHRLD:
                                if cov>tmp:
                                        THRLD_range.append((err,cov))
                                        tmp=cov
                                        break
        return THRLD_range


def first_max(values):
        tmp=(0,0)
        for v in values:
                if v[1]<tmp[1] and tmp[1]>10: break
                if v[1]>tmp[1]: tmp=v
        return tmp
             
def generatePutativeAllelesTXTFile(outTXT,stats_rslt):
        out=""
        fout = open(outTXT,'w')
        pAlleles = {}
        for sname in sorted(stats_rslt.keys()):
                values = stats_rslt[sname]
                pAlleles[sname]=[]
		for ref, statRef in values.items():
                        if stats_rslt[sname][ref]['IsPositiv'] and not stats_rslt[sname][ref]['IsParalog']:
                                pAlleles[sname].append(ref)
        refList=[]
        for ref in pAlleles.values():
                refList = refList + ref
        refList = sorted(list(set(refList)))
        
        out+="Sample\t{}\tNb Alleles\n".format("\t".join(refList))
        allelesSum={}
        for sample in sorted(pAlleles.keys()):
                refsample=pAlleles[sample]
                isRef = []
                for ref in refList:
                        if ref not in allelesSum.keys():
                                allelesSum[ref]=0
                        if ref in refsample:
                                isRef.append('1')
                                allelesSum[ref]+=1
                        else:
                                isRef.append('0')
                out+="{}\t{}\t{}\n".format(sample,"\t".join(isRef),sum([int(v) for v in isRef]))

        sumAlleles=[]
        for ref in refList:
                sumAlleles.append(str(allelesSum[ref]))
        out+="sum\t{}".format("\t".join(sumAlleles))

        fout.write(out)
        fout.close()

def generateXlsStatsFromDict(statdict, ReadsCounts, headers, xls_save_file):
	
	#output xls file generation
	
	style = xlwt.easyxf('alignment: horizontal center, vertical center;')
	hstyle = xlwt.easyxf('font: bold on; alignment: horizontal center, vertical center;')
	greenStyle = xlwt.easyxf('alignment: horizontal center, vertical center;pattern: pattern solid, fore_colour light_green;')
	orangeStyle = xlwt.easyxf('alignment: horizontal center, vertical center;pattern: pattern solid, fore_colour light_orange;')

	saveName = xls_save_file[:xls_save_file.rfind('.')]

	workbook = xlwt.Workbook()
	sheetNbr = 0
	fileNbr = 1

	for k1 in sorted(statdict.keys()):
                sample = statdict[k1]
		if len(k1)>31:
			stdout_print("Stats xls file: sheet name lenght '{}' to long truncate to '{}'".format(k1,k1[:31]))
		sheet = workbook.add_sheet(k1[:31])
		sheetNbr += 1

		readsNbr = sum(ReadsCounts[k1].values())
		
		sheet.write(0,0,"Reads count: {}".format(formatMillionNumbers(readsNbr)),style)
		
		row=0
		col=1
                xlsTitles = config['xls_headers_titles']
		for h in headers:
			sheet.write(row,col,xlsTitles[h],hstyle)
			col += 1
		
		row = 1
		col = 0
		
		items = sorted(sample.items(), key = SortKey, reverse=config['sortOrderReverse'])
		
		for ref, statRef in items:
			cstyle = style
			if statRef['mean cover']>0:
			
				if (statRef['IsPositiv']):
					cstyle = greenStyle
					if (statRef['IsParalog']):
						cstyle = orangeStyle
						
				sheet.write(row,0,ref,cstyle)
				col = 1
				for h in headers:
                                        cellContent = statRef[h]
                                        if isinstance(cellContent,list):
                                                cellContent = ",".join(cellContent)
					sheet.write(row,col,cellContent,cstyle)
					col += 1
				row += 1
                if(sheetNbr%config['MaxSheetNbr']==0):
                        remove_File("{}_{}.xls".format(saveName,fileNbr))
                        workbook.save("{}_{}.xls".format(saveName,fileNbr))
                        fileNbr += 1
                        workbook = xlwt.Workbook()
                        sheetNbr = 0
        if(sheetNbr>0):
                remove_File("{}_{}.xls".format(saveName,fileNbr))
                workbook.save("{}_{}.xls".format(saveName,fileNbr))
                        


def SortKey(tup):
	key, d = tup
	return d[config['SortKey']]

def formatMillionNumbers(number, separator="."):
	value = str(number)
	
	rslt = ''
	i = 1
	for v in range(len(value),0,-1):
		if i%3==0:
			rslt = separator + value[v-1] + rslt
		else:
			rslt = value[v-1] + rslt
		i+=1
			
	return rslt.lstrip(separator)

def bowtie_align():

	DB_dir = config['DB_dir']
	
	out_folder = config['outfolder']
	
	DB_Path = os.path.join(config['data_dir'],DB_dir)
	
	dir_IndexRef = os.path.join(DB_Path,config['bwt_IndexdirResult'])
	
	bwt_yam_reflist = os.path.join(dir_IndexRef, config['bwt_yam_reflist'])
	
	References_list = yaml.load(open(bwt_yam_reflist))
	
	dir_results = os.path.join(out_folder,config['bwt_dirResult'])
	
	dirFile_exist(dir_results,1)
	
	params = ' '.join(config['bowtie2_args'])
	
	jblst = utils.Jobslist("bowtie2 for all")
	
	addTextToLogFile("-- ALIGNMENTS WITH BOWTIE2 --")
	
	for reads in config['reads_files']:

		read_file = config['reads_files'][reads]
		rslt_path = os.path.join(dir_results,read_file[0])
		dirFile_exist(rslt_path,1)
		
                fastq_list = []
                for filename in read_file[1:]:
                        fastq_list.append(os.path.join(config['reads_dir'],filename))
		
                #unpaired allignment
                samples = "-U {fastqFiles}".format(fastqFiles=','.join(fastq_list))
			
		references_Names = [v for v in References_list.values() if v not in config['refexclude']]

		for reference in references_Names:
			aln_name = "Aln_{}_VS_{}".format(read_file[0],reference)
			aln_path = os.path.join(rslt_path,aln_name)
			
			target = os.path.join(aln_path, "{}.bam".format(aln_name))
			
			dirFile_exist(aln_path,1)
			
                        mapped_unmapped_out = "--al {mapped}".format(mapped=os.path.join(aln_path, "mapped_{}".format(aln_name)))
			
			cmd = "{bpath}bowtie2 {params} -x {RefIndex} {samples} {map_unmap} 2>{err_bowtie} | {path}samtools view -S -b -F 0x4 - > {bamFile}".format(bpath=bowtie2_path,path=samtools_path,params=params,RefIndex=os.path.join(dir_IndexRef,reference),samples=samples, map_unmap=mapped_unmapped_out, err_bowtie=os.path.join(aln_path, "stderr_{}".format(aln_name)), bamFile=os.path.join(aln_path, "{}.bam".format(aln_name)))
			
			#TO REMOVE START
#			test = random.choice(['','A','A','A','A'])
#			cmd = "bowtie2 {params} -x {RefIndex}{test} {samples} {map_unmap} 2>{err_bowtie} | samtools view -S -b -F 0x4 - > {bamFile}".format(test=test,params=params,RefIndex=os.path.join(dir_IndexRef,reference),samples=samples, map_unmap=mapped_unmapped_out, err_bowtie=os.path.join(aln_path, "stderr_{}".format(aln_name)), bamFile=os.path.join(aln_path, "{}.bam".format(aln_name)))
			#TO REMOVE END
			
			StderrPaths.append(os.path.join(aln_path, "stderr_{}".format(aln_name)))
			
			jblst.add_a_job(cmd,"bowtie2 {}".format(aln_name),target)
					
	addTextToLogFile("\tLaunch {} alignments".format(jblst.get_SizeOfJoblist()))
	
	if args.tinyverbose:
		t1 = threading.Thread(target=progress_StderrFiles)
		t1.start()
		
	utils.trun(args, jblst)
	
	if args.tinyverbose:
		Wait_ProgressBarFinished(t1)
		
	errorFileList = get_errorStderrFiles()
	
	msg = "{}/{} Alignements performed correctly".format(len(StderrPaths)-len(errorFileList),len(StderrPaths))
	stdout_print(msg)
	addTextToLogFile("\t{}".format(msg))
	
	if len(errorFileList)>0:
		msg = "{} Alignements not performed correctly see log file".format(len(errorFileList))
		stdout_print(msg)
		addTextToLogFile("\t{}".format(msg))
		
		for f in errorFileList:
			addTextToLogFile("\t\tSee: {}".format(f))

def Is_stderrFileInError(stderrFile):
	with open(stderrFile, "r") as f:
		data = f.read()
		if "(ERR)" in data:
			return True
	return False
	
def get_errorStderrFiles():
	
	errFileList=[]
	
	for ferr in StderrPaths:
		with open(ferr, "r") as f:
			data = f.read()
			if "(ERR)" in data:
				errFileList.append(ferr)
				
	return errFileList
	
def Wait_ProgressBarFinished(thrd):

	while(thrd.isAlive()):
		time.sleep(2)
	
def count_FinishedStderrFiles():
	c=0
	for f in StderrPaths:
		if os.path.exists(f) and os.path.getsize(f)>0:
			c+=1
	return c
	
def progress_FastqFilter(FilesNbr):
	c=0
	percent = 0

	out_folder = config['outfolder']
	filtered_dir = config['out_filtered']
	output = os.path.join(out_folder,filtered_dir)
	kmerFilter_Yaml_File = os.path.join(output,config['kmerFilter_Yaml'])
	
	while c<FilesNbr:
		time.sleep(2)
		
		if os.path.exists(kmerFilter_Yaml_File) and os.path.getsize(kmerFilter_Yaml_File)>0:
			f=yaml.load(open(kmerFilter_Yaml_File,"r"))
			c = len(f.keys())
		
		percent = int(100*float(c)/float(FilesNbr))
		
		progstr = "="*int(percent/2)+" "*(50-int(percent/2))
		
		printProgress("Kmer Filter progress: [{}] {}%".format(progstr,percent))
	print ""

def progress_StderrFiles():
	c=0
	tm = time.time()
	percent = 0
	estime = 0
	
	while c<len(StderrPaths):
		
		time.sleep(2)
		
		c = count_FinishedStderrFiles()
		
		tmpPercent = percent
		percent = int(round(100*float(c)/float(len(StderrPaths)),0))
		
		progstr = "="*int(percent/2)+" "*(50-int(percent/2))
		
		eltime = time.time() - tm
		
		if percent>0:
			if percent != tmpPercent: estime = (100-percent)*eltime/(percent)
			str_estime = " -- {}".format(time.strftime('%Hh%Mm%Ss left', time.gmtime(estime)))
		else:
			str_estime = ""
			
		printProgress("Alignments progress: [{}] {}%{}".format(progstr,percent,str_estime))
	print ""

def bowtie_index():

	DB_dir = config['DB_dir']

	out_folder = os.path.join(config['data_dir'],DB_dir)
	
	refs_dir = os.path.join(out_folder,config['ref_dir'])
	
	yml_qc_ref = os.path.join(out_folder,config['yml_qc_ref'])
	
	bwt_yam_reflist = config['bwt_yam_reflist']
	
	dir_results = os.path.join(out_folder,config['bwt_IndexdirResult'])
	
	bwt_out = config['bwt_build_out']
	
	dirFile_exist(dir_results,1)
	
	References_info = yaml.load(open(yml_qc_ref))
	
	jblst = utils.Jobslist("bowtie2-build for all")
	
	addTextToLogFile("-- INDEX A NEW REFERENCE --")
	
	ref_list = {}
	c=0
	for key, RefFiles in References_info.items():
		c += 1
		filename = os.path.join(refs_dir,key+'.fasta')
		if os.path.exists(filename):
			
			target = os.path.join(dir_results,"{}{}".format(key,'.1.bt2'))
			cmd = '{bpath}bowtie2-build {reffile} {refname} >{bwtout}'.format(bpath=bowtie2_path,reffile = filename,refname = os.path.join(dir_results,key),bwtout = os.path.join(dir_results,bwt_out))
			addTextToLogFile("\tBowtie2-build: {}".format(cmd))
			jblst.add_a_job(cmd,"bowtie2-build {}".format(key),target)
			ref_list['Ref_'+str(c)]= key
		else:
			errormsg = "ERROR: BOWTIE2_Index {} NOT FOUND\n".format(filename)
			exitProg(errormsg)
		
	create_YamlFromDict(os.path.join(dir_results,bwt_yam_reflist), ref_list)
	
	utils.trun(args, jblst)
		

#Quality control for samples
def samples_quality_cont():
	
	addTextToLogFile("-- DO SAMPLE QUALITY CONTROL --")
	
	samples_info = {}
	
	out_folder = config['outfolder']
	
	reads_files = config['reads_files']
	
	dir_results = os.path.join(out_folder,config['sampleQC_dir'])
	
	dir_allResults = os.path.join(out_folder,config['all_results_dir'])
	
	dir_reads = config['reads_dir']
	
	fastQC_stdout = os.path.join(dir_results,config['fastqc_stdout'])
	fastQC_stderr = os.path.join(dir_results,config['fastqc_stderr'])
			
	remove_File(fastQC_stdout)
	remove_File(fastQC_stderr)
	
	yamlReads_file = os.path.join(dir_allResults,config['yamlReads_file'])
	
	yamlReadsSamples = {}
	
	if os.path.exists(yamlReads_file):
		yamlReadsSamples = yaml.load(open(yamlReads_file))
	
	dirFile_exist(dir_results,1)
	dirFile_exist(dir_allResults,1)
	
	num_reads = 0
	
	NbrReads = {}
	
	for read in reads_files.values():

		sname = read[0]

		samples_info[sname]=read[1:]
	
	jblst = utils.Jobslist("FASTQC for all")
	
	for k,sample in samples_info.items():
	
		dirFile_exist(os.path.join(dir_results,k),1)
		
		fastq_files = sample
		
		NbrReads[k] = {}
		
		for filename in fastq_files:
			if is_url(filename):
                                stdout_print("\tNO FASTQC FOR URL: {}".format(filename),True)
                        else:
                                outputdir = os.path.join(dir_results,k)
                                target = os.path.join(outputdir,"{}_fastqc.html".format(filename[:filename.find('.')]))
                                cmd = "{pth}fastqc -o {output} {fastqfile} 2>> {stderr} 1>> {stdout}".format(output=outputdir,fastqfile=os.path.join(dir_reads, filename),stderr=fastQC_stderr,stdout=fastQC_stdout,pth=fastqc_path)
                                addTextToLogFile("\tFASTQC: {}".format(cmd))
			
                                if (k not in yamlReadsSamples.keys()) or args.force:
                                        num_reads = bufcount(os.path.join(dir_reads, filename))/4
                                        stdout_print("count Ok for {} : {}".format(filename,num_reads))
                                else:
                                        num_reads = yamlReadsSamples[k][filename]
				
			
                                NbrReads[k][filename] = num_reads
			
                                jblst.add_a_job(cmd,filename,target)
			
	for k,v in NbrReads.items():
		addTextToLogFile("\tSample: {Sample}\tfiles: [{files}]\tTotal reads: {totReads}".format(Sample=k,files=",".join(v.keys()),totReads=sum(v.values())))
	create_YamlFromDict(yamlReads_file, NbrReads)
	
	utils.trun(args, jblst, False)

def get_paramsFromSeqID(seqID):
	tabID = seqID.split("|")
	params = {}
	
	if len(tabID)>1:
		for p in tabID[1:]:
			kv = p.split("=")
			if len(kv)==2:
				params[kv[0]]=kv[1]
		
	return tabID[0],params

def HomoParaFromRef_AddGrpColor():
	GrpColorDic={}
	r = lambda: random.randint(0,255)
	
	for k,v in HomoParaFromRef.items():
		if 'grpRef' in v:
			GrpColorDic[v['grpRef']]=('#%02X%02X%02X' % (r(),r(),r()))
	for k,v in HomoParaFromRef.items():
		if 'grpRef' in v:
			v.update({'grpColor':GrpColorDic[v['grpRef']]})

def exitProg(msg=''):
	if msg!='':
		addTextToLogFile("PROGRAM EXIT\t{}".format(msg))
	sys.exit(msg)

def refDatabase_in_catalog(dbname,DBFormat=''):

	if DBFormat=='DB':
		orig_dbname = dbname
		tmp = os.path.join(config['data_dir'],dbname)
		dbname = os.path.join(tmp,"{}.fa".format(dbname))
		if not os.path.exists(dbname):
			exitProg("{} DATABASE NOT FOUND...".format(orig_dbname))
		
	filemd5 = getFileMD5(dbname)
	
	bname = os.path.basename(dbname)
	DB_dir = bname[:bname.rfind('.')]
	refCatalogFile = os.path.join(config['data_dir'],config['refs_catalog'])
	
	if os.path.exists(refCatalogFile):

		refCatalogDic = load_DicFromPickleFile(refCatalogFile)
		
		for k,dbVal in refCatalogDic.items():
			if not os.path.exists(os.path.join(config['data_dir'],dbVal)):
				del refCatalogDic[k]
		
		if filemd5 in refCatalogDic:
				config['DB_dir'] = refCatalogDic[filemd5]
				addTextToLogFile("Found Database : {}".format(config['DB_dir']))
				return os.path.exists(os.path.join(config['data_dir'],refCatalogDic[filemd5])), dbname
				
		else:
			
			if config['remove_OldDB']:
				if DB_dir in refCatalogDic.values():
					stdout_print("DELETE \"{}\" DATABASE FOLDER".format(DB_dir))
					shutil.rmtree(os.path.join(config['data_dir'],DB_dir), ignore_errors=True)
					refCatalogDic = dict((k,v) for k,v in refCatalogDic.items() if v!=DB_dir)
					addTextToLogFile("Database \"{}\" folder removed".format(DB_dir))

			refCatalogDic[filemd5] = DB_dir
			create_PickleFromDict(refCatalogFile,refCatalogDic)
			config['DB_dir'] = DB_dir
			addTextToLogFile("Create Database : {}".format(DB_dir))
			return False, dbname
	else:
	
		create_PickleFromDict(refCatalogFile,{filemd5 : DB_dir})
		
		config['DB_dir'] = DB_dir
		addTextToLogFile("Create Database : {}".format(DB_dir))
		return False, dbname

def split_fasta(IsRefInCatalog):
	
	addTextToLogFile("-- REFERENCE QUALITY CONTROL --")
	
	n_seq = 0
	n_bp = 0
	
	seq_info = {}
	
	change_list = []
	
	stderr_bases = ''
	stderr_gc = ''
	stderr_minlen = ''
	
	refs_file = config['ref_file']
	
	data_dir = config['data_dir']

	ref_name = config['ref_id']
	
	DB_dir = config['DB_dir']
	
	out_folder= os.path.join(data_dir,DB_dir)
	
	dir_results = os.path.join(out_folder,config['ref_dir'])
	
	species_dict_file = config['species_dict']
	
	yamlqc_file = os.path.join(out_folder,config['yml_qc_ref'])
	
	seqLen_distrib_graph = os.path.join(out_folder,config['seqLen_distrib_graph'])
	
	xls_save_file = os.path.join(out_folder,config['xls_results_qc'])
	
	addTextToLogFile("\tReference file: {}".format(refs_file))
	
	if not IsRefInCatalog or args.force:
		#test directories and files
		dirFile_exist(data_dir)
		dirFile_exist(refs_file)
		dirFile_exist(os.path.join(data_dir, species_dict_file))
		check_FastaFile(refs_file)
		dirFile_exist(dir_results,1)
		
		BnameRefs_file = os.path.basename(refs_file)
		
		copyFiles(refs_file,"{}/{}.fa".format(os.path.join(out_folder),BnameRefs_file[:BnameRefs_file.rfind('.')]))
		
		species_dict = yaml.load(open(os.path.join(data_dir, species_dict_file)))
	
		count_seqid = []
	
		seqLen_distrib = []
	
		for seq in SeqIO.parse(refs_file, "fasta"):

			n_seq += 1
			n_bp += len(seq)
		
			seqLen_distrib.append(len(seq))
		
			seqID, params = get_paramsFromSeqID(seq.id)
		
			count_seqid.append(seqID)
		
			check_duplicates(count_seqid,"in {}".format(refs_file),1)

			NewRef_id = str_normalize(seqID,change_list)
		
			fasta_Filename = dir_results + '/'+ NewRef_id+'.fasta'
		
			seq_info[NewRef_id] = {}
		
			seq_info[NewRef_id]['specie'] = getSpecieFromSeqName(seqID, species_dict)
	
			#inFile':val,'seqName':col_seqname[r],'species':col_species[r],'complete':col_complete[r],'accessNbr':col_accessNbr[r]
			seq_info[NewRef_id]['seq_id'] = seqID
			seq_info[NewRef_id]['seq_len'] = len(seq)
			seq_info[NewRef_id]['seq_description'] = seq.description.replace(seqID,'').strip(' ')
	
			seq_info[NewRef_id]['fragments'] = {}
	
			split_fragment = seq.seq.split(config['sep_motif'])
	
			n_fragment = 0
	
			for cont in split_fragment:
		
				if len(cont)>0:
					n_fragment += 1
					pos_cont = seq.seq.find(cont)
			
					seq_info[NewRef_id]['fragments']['F'+str(n_fragment)] = {}
			
					seq_info[NewRef_id]['fragments']['F'+str(n_fragment)]['pos'] = pos_cont
					seq_info[NewRef_id]['fragments']['F'+str(n_fragment)]['len'] = len(cont)
					seq_info[NewRef_id]['fragments']['F'+str(n_fragment)]['GC'] = round(SeqUtils.GC(cont),1)
                                        seq_info[NewRef_id]['fragments']['F'+str(n_fragment)]['ShEntropy'] = get_ShannonEntropy(cont)
			
					seq_info[NewRef_id]['fragments']['F'+str(n_fragment)]['bases_count'] = dict()
			
					base_count = []
			
					for i in ['N','U','W','S','M','K','R','Y','B','D','H','V','Z']:
			
						if (cont.count(i)>0):
				
							base_count.append(i+':'+str(cont.count(i)))
					
					seq_info[NewRef_id]['fragments']['F'+str(n_fragment)]['bases_count'] = ','.join(base_count)						
			
					if (len(seq_info[NewRef_id]['fragments']['F'+str(n_fragment)]['bases_count'])>0):
			
						stderr_bases += "Id:{}\tfragment:C{}\t{}\n".format(NewRef_id, str(n_fragment), ', '.join(base_count))
	
					if (SeqUtils.GC(cont)<config['GC_interval'][0] or SeqUtils.GC(cont)>config['GC_interval'][1]):
		
						stderr_gc += "WARNING:\t{}\tGC:{}%\n".format(NewRef_id,round(SeqUtils.GC(cont),1))
			
					if len(cont)<config['seq_min_len']:
							stderr_minlen += "WARNING: length of fragment:C{} ({}bp) in sequence >{} ({}bp)\n".format(str(n_fragment),len(cont),seqID,len(seq))
				
		
			seq.id = NewRef_id
			seq.description = ''
			with open(fasta_Filename, 'w') as output_handle:
				SeqIO.write(seq, output_handle, "fasta")
	
	
		create_YamlFromDict(yamlqc_file, seq_info)
	
		creatHist(seqLen_distrib,"References length distribution\n({} references)".format(len(seqLen_distrib)),'Frequency','bp',seqLen_distrib_graph,False)
	
		#stderr outputs
	
		resume = stderr_write("WARNING: DNA Alphabet on {} sequences".format(len(stderr_bases.split('\n'))-1), stderr_bases)
	
		resume += stderr_write("WARNING: GC% on {} sequences out of interval: [{}% < GC% < {}%]".format(len(stderr_gc.split('\n'))-1,config['GC_interval'][0],config['GC_interval'][1]), stderr_gc)
	
		resume += stderr_write("WARNING: Sequence length < {}bp on {} item(s)".format(config['seq_min_len'], len(stderr_minlen.split('\n'))-1), stderr_minlen)
	
		resume += stderr_write("WARNING: String normalization on {} item(s)".format(len(change_list)),'\n'.join(change_list) + '\n')
	
		resume += stderr_write("INFORMATION: Sequences per species",seq_per_species(seq_info))
	
		resume += stderr_write("END","{} seq - {} bp\n".format(n_seq, n_bp))
	
		xlsgen_refQC(yamlqc_file, xls_save_file, resume)

	else:
		stdout_print("{refname} Found...".format(refname=DB_dir))
	
	copyFiles(xls_save_file,os.path.join(config['outfolder'],"{}_{}".format(DB_dir,config['xls_results_qc'])))
	copyFiles(seqLen_distrib_graph,os.path.join(config['outfolder'],"{}_{}".format(DB_dir,config['seqLen_distrib_graph'])))

def addTextToLogFile(logtxt,printTime=True):

	if printTime:
		ttime = "{}\t".format(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"))
	else:
		ttime=""
	
	logFolder = os.path.join(config['outfolder'],config['log_folder'])

	dirFile_exist(logFolder,1)

	logFile = os.path.join(logFolder,config['log_file'])
	lf = open(logFile,'a+')
	if logtxt.count('\n')>0:
		for l in logtxt.split('\n'):
			lf.write("{tm}{log}\n".format(tm=" "*len(ttime),log=l))
	else:
		lf.write("{tm}{log}\n".format(tm=ttime,log=logtxt))
	lf.close()
	
def load_HomoParaFromRef():
	
	refs_file = config['ref_file']
	
	for seq in SeqIO.parse(refs_file, "fasta"):
		
		seqID, params = get_paramsFromSeqID(seq.id)
		
		HomoParaFromRef[str_normalize(seqID)]= params
	
	HomoParaFromRef_AddGrpColor()

def getSpecieFromSeqName(seqName, species_dict):
        if type(species_dict) is dict:
                for k in species_dict.keys():
                        if seqName.find(k)==0:
                                return species_dict[k]
	
	return 'Other'
	
def check_FastaFile(fasta_Filename):

	fasta = open(fasta_Filename)
	
	rslt = True
	
	error = ""
	
	for line in fasta:
		if line[0]=='>':
			if line.find(' ')>-1:
				rslt = False
				error += line
	if rslt:
		return True
	else:
		errormsg = "ERROR: spaces in fasta sequences name -- {fname}\n{errors}".format(fname=fasta_Filename,errors=error)
		exitProg(errormsg)
			
			
def check_duplicates(datas, message, exitprog=0):
		
	c = 0
	
	for k,nbr in collections.Counter(datas).items():
		if nbr>1:
			stderr_infile_dup = "{} - present {} times\n".format(k,str(nbr));
			c += 1
			
	if c>0:
		stderr_write("ERROR duplicates: {}".format(message), stderr_infile_dup,exitprog)
	
	return c
	
# References quality controle XLS generation from yaml
def xlsgen_refQC(yamlfile, xls_save_file, resume=''):
	
	refQC_Yaml = yaml.load(open(yamlfile))
	
	#output xls file generation
	
	style = xlwt.easyxf('alignment: horizontal center, vertical center;')
	
	workbook = xlwt.Workbook()
	
	sheet = workbook.add_sheet('S-LocGen_RefQC')
	sheet_resume = workbook.add_sheet('Resume')
	
	#resume sheet
	rownbr = 0
	for text in resume.split("\n"):
		sheet_resume.write(rownbr,0,text)
		rownbr += 1
	
	rownbr = 0
	
	colnbr = 0
	
	keyslist = ['specie', 'seq_id', 'seq_len']
	
	fragmentkeyslist = []
	
	#col headers
	maxcont = 0
	fragmentIdnbr = 0
	for k,v in refQC_Yaml.items():
	
		nbrcont = len(v['fragments'])
		
		if (fragmentIdnbr == 0 and 'F1' in v['fragments']):
			fragmentIdnbr = len(v['fragments']['F1'])
			fragmentkeyslist = v['fragments']['F1'].keys()
			
		if (nbrcont>maxcont):
			maxcont = nbrcont
	
	for headers in keyslist:
		sheet.write_merge(rownbr, rownbr+1, colnbr, colnbr, headers,style)
		colnbr += 1
	
	for c in range(1,maxcont+1):
		sheet.write_merge(rownbr, rownbr, colnbr, colnbr-1+fragmentIdnbr,'fragment F'+str(c),style)
		
		coltmp = colnbr
		for items in fragmentkeyslist:
			sheet.write(rownbr+1,coltmp,items,style)
			coltmp += 1
		
		colnbr += fragmentIdnbr
	
	#fill row with datas
	rownbr = 2
	
	for k,v in refQC_Yaml.items():
		
		colnbr = 0
		
		v['ref_id'] = k
		
		for keys in keyslist:
			
			if keys in v:
				cell_cont = v[keys]
			else:
				cell_cont = ''
				
			sheet.write(rownbr,colnbr,cell_cont,style)
			colnbr += 1
	
		for kc,kv in v['fragments'].items():
			
			for fragmentkey in fragmentkeyslist:
				sheet.write(rownbr,colnbr,kv[fragmentkey],style)
				colnbr += 1
	
		rownbr += 1
	
	workbook.save(xls_save_file)


def seq_per_species(seq_info):
	
	rslt = {}
	
	for seq_id in seq_info:
		if 'specie' in seq_info[seq_id]:
			specie = seq_info[seq_id]['specie']
		else:
			specie = 'not specified'
			
		if not specie in rslt:
			rslt[specie] = 1
		else:
			rslt[specie] += 1

	str_rslt = "Number of species {}  -- total number of sequences {}\n".format(len(rslt),sum(rslt.values()))
	
	str_rslt += "specie\tsequence(s)\n"
	
	for spe_val in rslt:
		str_rslt += "{}\t{}\n".format(spe_val,rslt[spe_val])
		
	return str_rslt

#stderr output
def stderr_write(title,stderr_str,exitprog=0):
	
	stderr_text = ''
	
	if len(stderr_str.split('\n'))>1 and args.tinyverbose:
		stderr_text = "----- " + title + " -----\n"
		stderr_text += stderr_str + "\n"
		
		sys.stderr.write(stderr_text)
		addTextToLogFile(stderr_text)
		
	if exitprog:
		sys.exit("Exit...")
		
	return stderr_text

def creatHist(data,Title,xlabel,ylabel,pngfilename,cumulative):
	fig, ax = plt.subplots()
	x = np.array(data)
	plt.clf()
	histval = plt.hist(x,histtype='stepfilled', cumulative=cumulative)
	plt.grid(True)
	plt.title(Title)
	plt.ylabel(xlabel);
	plt.xlabel(ylabel)
	plt.xticks(range(0, max(data),1000))
	plt.savefig(pngfilename)

def getFrequence(data,excludeZero=False):

	if excludeZero:
		data = remove_values_from_list(data,0)
	
	x = np.array(data)
	
	unique, counts = np.unique(x, return_counts=True)
	
	return [unique, counts]

def checkSelection_param(stats,YamlSelectParam):
	plot_selection_param = config[YamlSelectParam]
	
	for param in plot_selection_param:
		eval_str = "bool({})".format(str(stats[param[0]])+str(param[1])+str(param[2]))
		if not eval(eval_str): return False
	
	return True

def GrpRomanNumber(number):
	romain = {'1':'I','2':'II','3':'III','4':'IV','5':'V','6':'VI','7':'VII','8':'VIII','9':'IX','0':'?'}
	if str(number) in romain:
		return romain[str(number)]
	else:
		return str(number)

def add_AppFolderToYamlConfig(KeysList):
	for key in KeysList:
		if config[key]!=None:
			config[key] = os.path.join(App_Folder, config[key])

def load_readsList(readsList):
	rslt = {}
	
	rl = open(readsList,'r')
	
	rslt['reads_dir'] = rl.readline().strip()
	
	rslt['reads_files'] = {}
	
	i=0
	for read in rl:
                commentID = read.find('#')
                if commentID<0:
                        tmp = read.strip()
                elif commentID>0:
                        tmp = read[:commentID].strip()
                if commentID!=0:
                        read_info = tmp.split(',')
                        if len(read_info)>=2:
                                rslt['reads_files']["reads_{}".format(str(i))] = read_info
                                i+=1
	
	rl.close()
	
	return rslt

#=================================================================================================================================

#********************
#* Varous Functions *
#********************

def check_url_exists(url):
        try:
                code = urlopen(url).code
                return True
        except:
                return False

def is_url(url):
        return urlparse.urlparse(url).scheme in ('http','https','ftp', 'ftps')

def get_ShannonEntropy(dna,alphabet="ATCG"):
        ShEntropy = 0
        Udna = dna.upper()
        for base in alphabet:
                freq = float(Udna.count(base))/float(len(Udna))
                if freq==0: freq = 0.00000001
                ShEntropy = ShEntropy + freq*math.log(freq,2)
        return -ShEntropy

def extract_paramsFromReadsConfig(readsConfig):
        params = {}
        for line in readsConfig.values():
                idName = line[0]
                if idName not in params.keys():
                        params[idName]={}
                toRemove=[]
                for v in line:
                        if '=' in v:
                                tmp = v.split('=')
                                params[idName][tmp[0]]=tmp[1]
                                toRemove.append(v)
                for r in toRemove:
                        line.remove(r)

        return readsConfig, params

def copyFiles(src,dest):
	srcBname = os.path.basename(src)
	if not os.path.exists(os.path.join(dest,srcBname)):
		try:
			shutil.copy2(src,dest)
		except:
			stdout_print("Except CopyFiles src={} -- dst={}".format(src,dest),True) 

def remove_File(filename):
	if os.path.exists(filename):
		try:
			os.remove(filename)
			return True
		except:
			return False
	else:
		return False

def remove_values_from_list(the_list, val):
	return [value for value in the_list if value != val]

def str_normalize(cstr,change_list=None):
	
	orig = cstr
		
	patterns = {')':'_', '(':'_', '|':"\t", ' ':'_','>':'','.':'_'}

	for key,val in patterns.items():
		cstr = cstr.replace(key,val)
	
	if (cstr != orig) and change_list != None:
		change_list.append("orig:{}\tnorm:{}".format(orig,cstr))

	return str(cstr)

#Test if a file or directory exist if creatdir set to 0: exit
# if creatdir set to 1: create new directory
def dirFile_exist(path,creatdir = 0):

	if not os.path.exists(path):
		if (creatdir == 0):
			exitProg("ERROR: {} not found".format(path))
		else:
			dirmsg = "-- create {} directory".format(path)
			stdout_print(dirmsg)
			os.makedirs(path)

#use for yaml add_representer function
def quoted_presenter(dumper, data):
	return dumper.represent_scalar('tag:yaml.org,2002:str', data, style="'")	

def getFileMD5(filename):
	return hashlib.md5(open(filename, 'rb').read()).hexdigest()

def IsKeyPresentInDic(key,dic,ValueIfNorPresent):
	if key not in dic:
		dic[key]=ValueIfNorPresent

def check_FileExtensions(filename, extList):
	fileExt = filename[filename.rfind('.')+1:]
	
	if fileExt.upper() in [ext.upper() for ext in extList]:
		return True
	
	return False

def create_PickleFromDict(pickleFilePath, dic):
	pf = open(pickleFilePath,'wb')
	pickle.dump(dic,pf,protocol=pickle.HIGHEST_PROTOCOL)
	pf.close()

def load_DicFromPickleFile(pickleFilePath):
	pf = open(pickleFilePath,'rb')
	retVal = pickle.load(pf)
	pf.close()
	return retVal

def castStrNbr_FloarOrInt(strNbr):
	
	try:
	
		if strNbr.find('.')>=0 or strNbr.upper().find('E')>0:
			#float
			tmp = float(strNbr)
		else:
			tmp = int(strNbr)
			
	except:
	
		tmp = strNbr
		
	return tmp

def create_YamlFromDict(yamlOutput_file, ymldict, overwrite=True, add_representer = True):

		if not os.path.exists(yamlOutput_file) or overwrite:
			#add single quotes for str variables
			if add_representer: yaml.add_representer(str, quoted_presenter)
		
			with open(yamlOutput_file, 'w') as yaml_file:
				yaml.dump(ymldict, yaml_file, default_flow_style=False)

def printProgress(txt):
	sys.stdout.write("\r"+str(txt))
	sys.stdout.flush()

def bufcount(filename):

	number_lines = sum(1 for line in open(filename))
	
	return number_lines

#=================================================================================================================================

#***************************
#* Others program features *
#***************************

def databases_list(ArgsVal):

	description = """ database list """
	
	parser = argparse.ArgumentParser(prog="databases",description=description)
		
	parser.add_argument("--config", help="config file")
	parser.add_argument("-r","--removedb", help="database to remove")
	args = parser.parse_args(ArgsVal)
	
	if args.config:
		config_file = args.config
	else:
		config_file = os.path.join(App_Folder,"config.yaml")
		
	config = yaml.load(open(config_file))

	tmp = os.path.join(config['data_dir'],config['refs_catalog'])
	refCatalogFile = os.path.join(App_Folder,tmp)
	
	if os.path.exists(refCatalogFile):
		catalog = load_DicFromPickleFile(refCatalogFile)
		
		if args.removedb:
			dbfound = False
			for k,v in catalog.items():
				if v==args.removedb:
					dbfound = True
					confirm = raw_input("Delete database {} (y/n): ".format(args.removedb))
					if confirm.upper()=='Y':
						del catalog[k]
						cf = open(refCatalogFile,"w")
						pickle.dump(catalog,cf)
						cf.close()
						shutil.rmtree(os.path.join(config['data_dir'],args.removedb), ignore_errors=True)
						print "{} deleted\n\n".format(args.removedb)
			if not dbfound:
				print "database \'{}\' not found... \n".format(args.removedb)
		
		print "\n-- List of available databases -- \n"
		print "type database name after -d option in genotyp feature\n"
		print "{} database(s) found:".format(len(catalog))
		
		for v in catalog.values():
			print "-\t{}\n".format(v)
	else:
		print "Catalog file not found. It while be created after first use"

def change_MaxJobConfig(ArgsVal):
	
	description = """ Change Max parallels number of jobs runing """
	
	parser = argparse.ArgumentParser(prog="chMaxJob",description=description)
		
	parser.add_argument("-u", "--utilsparam", help="'XXXX_utils_params' file")
	parser.add_argument("-T","--MaxParallelsJobs", help="max number of parallels jobs to run", type=int, required=True)
	args = parser.parse_args(ArgsVal)
	
	if args.utilsparam:
		yamlOutput_file = args.utilsparam
	else:
		yamlOutput_file = utils.params_file

	if os.path.exists(yamlOutput_file):
	
		config = yaml.load(open(yamlOutput_file))

		oldvalue = config['MaxParallelJobs']
		newvalue = args.MaxParallelsJobs
		
		if oldvalue!=newvalue:

			config['MaxParallelJobs'] = newvalue

			with open(yamlOutput_file, 'w') as yaml_file:
				yaml.dump(config, yaml_file, default_flow_style=False)
				
			print "MaxParallelJobs changed from {} to {}".format(oldvalue,config['MaxParallelJobs'])
		else:
			print "MaxParallelJobs already to {}".format(oldvalue)
	else:
		print "'{}' NOT FOUND. You can specifie a location by using -u option (see help)".format(yamlOutput_file)
	
def show_program_version():
	print "{} {}".format(AppName,__version__)

#Use splitFastq library, see C source in TOOLS directory
def run_splitFastq(ArgsVal):
	description = """ Split one or more fastq file(s) in 2 files (_R1_RSplit.fastq and _R2_RSplit.fastq) """
	
	parser = argparse.ArgumentParser(prog="splitFastq",description=description)
	parser.add_argument("-v", "--verbose", help=argparse.SUPPRESS,action="store_const", const="-v", default="")
	parser.add_argument("-t", "--test", help=argparse.SUPPRESS, action="store_const", const="-t", default="")
	parser.add_argument("-F", "--force", help=argparse.SUPPRESS, action="store_const", const="-F", default="")
	parser.add_argument("-C", "--checks", help=argparse.SUPPRESS, action="store_const", const="-c", default="")
	
	parser.add_argument("-o","--outfolder", help="output directory for splited fastq")
	parser.add_argument("-f","--fastqfiles", help="reference database in fasta format (see documentation)",nargs='*' , required=True)
	parser.add_argument("-c","--cutpos", help="cut position: if this option is not set the read is half cut", type=int)
	args = parser.parse_args(ArgsVal)
	
	outf=""
	cutpos = 0
	
	if args.outfolder:
		outf = args.outfolder
	
	if args.cutpos:
		cutpos = args.cutpos
	
	fasqL = args.fastqfiles
	
	splitFastq(cutpos, outf, fasqL, args = args)

def splitFastq(cutpos, outf, fasqL, args = args):
	
	cutpos_str = ""
	
	cur_dir = os.getcwd()
	
	if outf=="":
		outf = os.getcwd()
		
	elif outf[0]!="/":
		outf = os.path.join(cur_dir,outf)
	
	if not os.path.exists(outf):
		os.mkdir(outf)
	
	outf_str=" -o {}".format(outf)
	
	if cutpos>0:
		cutpos_str = " -c {}".format(cutpos)
	
	fq = []
	for f in fasqL:
		if f[0]!="/":
			fq.append(os.path.join(cur_dir,f))
		else:
			fq.append(f)
	cmd = "{pth}/splitFastq_v1/splitFastq{ct}{out} {fq}".format(pth=tools_folder,ct=cutpos_str,out=outf_str,fq=" ".join(fq))
	os.system(cmd)
	
def RapidFastQC(ArgsVal):
	import RapidFastQC
	RapidFastQC.run(ArgsVal)

def SraGetDatas(ArgsVal):
	import sra_getmeta
        db_argPath('-F',ArgsVal)
	sra_getmeta.run(ArgsVal)

def RunHaploAsm(ArgsVal):
	import haploAsm
        db_argPath('-db',ArgsVal)
	haploAsm.run(ArgsVal)

def RunGraphs_genotyp(ArgsVal):
        import graphs_genotyp
        graphs_genotyp.run(ArgsVal)

def RunExtractFromFasta(ArgsVal):
        import extractFromFasta
        extractFromFasta.run(ArgsVal)

def db_argPath(argOpt,ArgsVal):

	if argOpt in ArgsVal:
	
		config_file = os.path.join(App_Folder,"config.yaml")
		config = yaml.load(open(config_file))

		db_dir = os.path.join(App_Folder,config['data_dir'])
		
		dbIdx = ArgsVal.index(argOpt)+1
		dbname = ArgsVal[dbIdx]
		
		refFasta_dir = os.path.join(db_dir,dbname)
		refFasta = os.path.join(refFasta_dir,"{}.fa".format(dbname))
		
		if dbname.find('/')<0 and not check_FileExtensions(dbname, ['fa','fasta']):
			if os.path.exists(refFasta):
				ArgsVal[dbIdx] = refFasta
	
#Show help message for NGSgenotyp
def show_help(ProgFeatures):
	print "usage: {} [feature]\n".format(AppName)
	print "features list:"
	maxLength = max([len(k) for k in ProgFeatures.keys()])
	
	for feature,value in ProgFeatures.items():
		featureTitle = "{}{}".format(feature," "*(maxLength - len(feature)))
		print "\t{}\t-- {}\n".format(featureTitle,value['description'])


#=================================================================================================================================

#***********************
#* Program ENTRY POINT *
#***********************

if	__name__ == '__main__':
	
	#Select functionality
	ProgFeatures = {'help':{'Nbr':0,'description':"show this help message"}\
	,'genotyp':{'Nbr':1, 'description':"execute genotyp pipeline from raw NGS reads data"}\
	,'databases':{'Nbr':2, 'description':"list available databases"}\
	,'chMaxJob':{'Nbr':3, 'description':"change the number of background jobs launched during genotyp pipeline execution."}\
	,'kmerRefFilter':{'Nbr':4, 'description':"fastq raw files filtering by kmers dictionnary generated from references sequences"}\
	,'version':{'Nbr':5, 'description':"program version"}\
	,'splitFastq':{'Nbr':6, 'description':"Split reads of one or more fastq file(s) in 2 files (_R1_RSplit.fastq and _R2_RSplit.fastq)"}\
	,'RapidFastQC':{'Nbr':7, 'description':"A qualitative tool for quick quality check on fastq raw files (.fastq or .gz)"}\
	,'SRAGetDatas':{'Nbr':8, 'description':"download sequencing data and metadata from SRA number"}\
        ,'haploAsm':{'Nbr':9, 'description':"haplotypes assembly from genotyps mapped reads"}\
        ,'graphs_genotyp':{'Nbr':10, 'description':"draw genotyps score plots in pdf file(s)"}\
        ,'extractFromFasta':{'Nbr':11, 'description':"extract sequences from multiple haploAsm contigs files"}}

	func = None
	
	if len(sys.argv)>1:
	
		func = sys.argv[1]
		
		if (func in ProgFeatures):
			choice = ProgFeatures[func]['Nbr']
		
			if choice == 0:
				show_help(ProgFeatures)
			
			elif choice == 1:
				genotyp(sys.argv[2:])
			
			elif choice == 2:
				databases_list(sys.argv[2:])
			
			elif choice == 3:
				change_MaxJobConfig(sys.argv[2:])
		
			elif choice == 4:
				import kmerRefFilter
				kmerRefFilter.kmerRefFilter(sys.argv[2:])
				
			elif choice == 5:
				show_program_version()
				
			elif choice == 6:
				run_splitFastq(sys.argv[2:])
				
			elif choice == 7:
				RapidFastQC(sys.argv[2:])
				
			elif choice == 8:
				SraGetDatas(sys.argv[2:])
				
			elif choice == 9:
				RunHaploAsm(sys.argv[2:])

			elif choice == 10:
				RunGraphs_genotyp(sys.argv[2:])

                        elif choice == 11:
				RunExtractFromFasta(sys.argv[2:])

		else:
			print "{} : unknown feature\n".format(func)
			show_help(ProgFeatures)
	else:
		show_help(ProgFeatures)
	
