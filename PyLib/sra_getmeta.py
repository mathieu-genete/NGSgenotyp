#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# ref: https://www.ncbi.nlm.nih.gov/books/NBK242621/
# ref: https://www.ncbi.nlm.nih.gov/books/NBK56913/
# sratoolkit: https://github.com/ncbi/sra-tools/wiki/Downloads

import sys
import os

AppPath = os.path.realpath(__file__)
App_Folder = AppPath[:AppPath.rfind('/')]
sratoolkit_folder = os.path.join(App_Folder,'../TOOLS/sratoolkit/bin')

import csv
import argparse
import subprocess
import urllib2
import time
import signal
from datetime import datetime

__version__ = "v1.4"

args = None
UseLocalSraToolKit = False

nohupMode = False

# if need url can be changed, just replace query by: {qry}
metadataURL = "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term={qry}"

def run(ArgsVal):

	global args
        global nohupMode
	
	description = """ Get SRA raw files and Metadatas """

	parser = argparse.ArgumentParser(description=description)

	parser.add_argument("-m", "--metaDataOnly", help="get metadata only for given SRA", action="store_const", const="-m", default="")
	parser.add_argument("-f", "--localSratoolkit", help="", action="store_const", const="-f", default="")
	parser.add_argument("-s","--sraquery", help="SRA Query")
	parser.add_argument("-l","--listOfSra", help="text file with list Of SRA Query (one per line)")
	parser.add_argument("-d","--destdir", help="destination directory", required=True)
	parser.add_argument("-o","--CSVoutputfile", help="csv output file")
        parser.add_argument("-x","--addfilter", help="download data wich correspond to filer: header=value")
	parser.add_argument("-F","--filterFromRef", help="directly filter reads while download with reference")
	parser.add_argument("-ks","--kmerSize", help="kmer size for kmer filtering - default = 20", type=int, default=20)
	parser.add_argument("-mk","--minMatchKeepSeq", help="minimum number of ref reads matched to keep fastq read -- default=1",default=1, type=int)
	parser.add_argument("-z","--minShEntropy", help="minimum Shannon Entropy for kmers (0 to 2.0)-- defaut=0.8",default=0.8, type=float)
	parser.add_argument("-q","--maxratioAmbigous", help="maximum ambigous bases accepted in kmers in %% -- defaut=0.2",default=0.2, type=float)
	args = parser.parse_args(ArgsVal)
		
	kmersize=args.kmerSize
        minMatchKeepSeq=args.minMatchKeepSeq
        minShEntropy=args.minShEntropy
        maxratioAmbigous=args.maxratioAmbigous

	query = args.sraquery
	destdir = args.destdir
        absdestDir = os.path.abspath(destdir)

        launch_datetime=datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
        outRawFileName = os.path.join(destdir,"raw_data_list_{}".format(launch_datetime))
	proclist = []
	
        noCacheConfiguration()

        if not signal.getsignal(signal.SIGHUP) == signal.SIG_DFL:
                print "nohup mode"
                nohupMode=True

	if not is_tool("fastq-dump") or args.localSratoolkit:
		if os.path.exists(os.path.join(sratoolkit_folder,"fastq-dump")):
			UseLocalSraToolKit = True
		else:
			sys.exit("ERROR: fastq-dump not found or not installed...")
	
	if args.CSVoutputfile:
		outcsv = args.CSVoutputfile
                if ".csv" not in outcsv:
                        outcsv="{}.csv".format(outcsv)
	else:
		outcsv = "SRAs_Metadatas_{}.csv".format(launch_datetime)
	
	if query: #and query[:3].upper() in ["SRA","SRP","SRX","DRA","DRP","DRX","PRJ","ERA","ERP","ERX"]:
		csvrslt = GetFrom_URL(metadataURL.format(qry=query))
		csvtab = [l.split(',') for l in csvrslt.split("\n") if len(l)>0]
		ListFile = [l[csvtab[0].index("Run")] for l in csvtab[1:]]
	
	if args.listOfSra:
		ListFile = [l.strip() for l in open(args.listOfSra)]

	MetadataList = {}
	ErrorSRA = []
	LineNbr = len(ListFile)
	
	print "{} samples to download".format(LineNbr)
	
	i = 0
	n = 0
	
	printProgress("GET Metadatas 0%")
	
	for SRANbr in ListFile:
		SRANbr = SRANbr.rstrip()
		tmpMD = GetSRAMetadatas(SRANbr, destdir, None)

		if len(tmpMD)>1:
                        if args.filterFromRef:
                                tmpMD['Fastq File'] = "{}_filtered.fastq".format(tmpMD['Run'])
                        else:
                                tmpMD['Fastq File'] = "{}.fastq.gz".format(tmpMD['Run'])
			MetadataList[SRANbr] = tmpMD
			n += 1
		else:
			ErrorSRA.append(SRANbr)
		i += 1
		printProgress("GET Metadatas {}%".format(int(100*i/LineNbr)))
		
	print "\n\nget {} metadatas\n".format(n)
	
	if (i-n>0):
		print "*** ERROR for {} SRA number. See errors log file ***\n".format(i-n)
		
	KeysList = []
	for values in MetadataList.values():
		for k in values.keys():
			KeysList.append(k)
	
	totalSize_MB = sum(int(v['size_MB']) for v in MetadataList.values())
	
        if args.filterFromRef:
                print "Filtering fastq with '{}' -- kmers size = {}\n".format(os.path.basename(args.filterFromRef),kmersize)
        else:
                print "about {} MB while be downloaded\n".format(totalSize_MB)
	
	Headers = list(set(KeysList))
	
	generate_csv_File(os.path.join(destdir,outcsv), Headers, MetadataList)

        if args.addfilter:
                filterTxt = str(args.addfilter).split('=')
                headFilter = filterTxt[0]
                valFilter = filterTxt[1]
                print "filter: [{}]={}".format(headFilter,valFilter)
	outRawFile={}
	for val in MetadataList.values():
                if args.addfilter:
                        if headFilter in val.keys() and val[headFilter]==valFilter:
                                launchDwn = True
                        else:
                                launchDwn = False
                else:
                        launchDwn = True

                if launchDwn:
                        if val['SampleName'] not in outRawFile.keys():
                                outRawFile[val['SampleName']]=[]
                        outRawFile[val['SampleName']].append(val['Fastq File'])

                        proc = downloadFastqFromSRA(destdir, val['Run'],kmersize,minMatchKeepSeq,minShEntropy,maxratioAmbigous)
                        proclist.append([proc, val['Run']])
	
	if len(ErrorSRA)>0:
		errfile = open(os.path.join(destdir,"sra_getmeta_errors.log"),'w')
		for sranbr in ErrorSRA:
			errfile.write("No metadatas for SRA: {}".format(sranbr))
		errfile.close()

        outRfile = open(outRawFileName,"w")
        outRfile.write("{}\n".format(absdestDir))
        for k,v in outRawFile.items():
                outRfile.write("{},{}\n".format(k,",".join(v)))
        outRfile.close()
	waitProcTermination(proclist,destdir)

def noCacheConfiguration(ncbiConfigfolder=".ncbi",userSettings="user-settings.mkfg"):
        homefolder = os.path.expanduser('~')
        ncbiConfigPath = os.path.join(homefolder,ncbiConfigfolder)
        userSettingPath = os.path.join(ncbiConfigPath,userSettings)
        cacheDisabledLine="/repository/user/cache-disabled"
        cdLineExists = False
        if not os.path.exists(ncbiConfigPath):
                os.mkdir(ncbiConfigPath)

        allLines=[]
        if os.path.exists(userSettingPath):
                setFile = open(userSettingPath,'r')
                for line in setFile:
                        if cacheDisabledLine in line:
                                allLines.append('{} = "true"'.format(cacheDisabledLine))
                        else:
                                allLines.append(line.strip())
                setFile.close()
                os.remove(userSettingPath)
        else:
                allLines.append('{} = "true"'.format(cacheDisabledLine))
        
        outSetFile = open(userSettingPath,'w')
        for line in allLines:
                outSetFile.write("{}\n".format(line))
        outSetFile.close()

def waitProcTermination(proclist,destdir):
	if not args.metaDataOnly:
                sraTerminatedLog = os.path.join(destdir,"sra_getmeta_terminated")
		procTermn = 0
                if os.path.exists(sraTerminatedLog):
                        os.remove(sraTerminatedLog)

		os.mknod(sraTerminatedLog)
		while procTermn < len(proclist):
		
			time.sleep(1)
			
			oldProcTermn = procTermn
			
			procTermn = sum(1 for p in proclist if p[0].poll()!=None)
			
			if oldProcTermn != procTermn:
				trfile = open(sraTerminatedLog,"w")
				trfile.write("\n".join(["{}.fastq.gz".format(p[1]) for p in proclist if p[0].poll()!=None]))
				trfile.close()
				
			printProgress("Download progression {} / {} file(s)".format(procTermn,len(proclist)))
		
		print "Terminated..."
			
def printProgress(txt):
        if not nohupMode:
                sys.stdout.write("\r"+str(txt))
                sys.stdout.flush()
	
def generate_csv_File(csv_filename, headers, datas):
	with open(csv_filename, 'wb') as csvfile:
		csvwt = csv.writer(csvfile)
		csvwt.writerow(headers)
		for line in datas.values():
			tmp = []
			for k in headers:
				tmp.append(line[k])
			csvwt.writerow(tmp)

def replace_comma_betweenQote(data):
        rslt = ''.join(x if i%2==0 else x.replace(',','_') for i,x in enumerate(data.split('"')))
        return rslt

def GetSRAMetadatas(query, destdir, csvfilename=None):

	url = metadataURL.format(qry = query)
	
	SRA_MetaData = replace_comma_betweenQote(GetFrom_URL(url))

	if SRA_MetaData:
	
		if csvfilename:
			f = open(os.path.join(destdir, csvfilename), 'w')
			f.write(SRA_MetaData)
			f.close()
			
		headers = SRA_MetaData.split("\n")[0].split(",")
		values = SRA_MetaData.split("\n")[1].split(",")
	
		rslt = dict_metadataLine(headers,values)
	
		return rslt
		
	else:
		return False

def dict_metadataLine(headers,values):
	rslt = {}
	if len(headers) == len(values):
		i=0
		for heads in headers:
			rslt[heads]=values[i]
			i+=1

	return rslt

def downloadFastqFromSRA(destinationFolder, Run ,kmersize,minMatchKeepSeq,minShEntropy,maxratioAmbigous):

	if not args.metaDataOnly:
			
                if args.filterFromRef:
                        print "Launch Download and filter for {}".format(Run)
                else:
                        print "Launch Download for {}".format(Run)
		
                cmdFld = ""

		if UseLocalSraToolKit:
		
			if destinationFolder[0]=="/":
				destFold = destinationFolder
			else:
				destFold = os.path.join(os.getcwd(),destinationFolder)
			cmdFld="{}/".format(sratoolkit_folder)				

                else:
                        destFold = destinationFolder

                if args.filterFromRef:
                        cmd = "{stkpath}fastq-dump -Z {qry} | {appfld}/kmerRefFilter.py -k {kmsize} -m {minMatch} -z {minSh} -q {maxratio} -r {ref} -s {qry} -o {outf} 2>> {outf}/kmerRefilter_stdout".format(stkpath=cmdFld,qry=Run, outf=destFold, appfld=App_Folder,ref=args.filterFromRef,kmsize=kmersize,minMatch=minMatchKeepSeq,minSh=minShEntropy,maxratio=maxratioAmbigous)
                else:
                        cmd = "{stkpath}fastq-dump --gzip -O {outf} {qry}".format(stkpath=cmdFld,qry=Run, outf=destFold)
		proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		return proc
		
	return None

def GetFrom_URL(url):

	attempts = 0
	
	while attempts < 3:
		try:
			response = urllib2.urlopen(url, timeout = 15)
			content = response.read()
			
			return content
			
		except urllib2.URLError as e:
			attempts += 1
			sys.stderr.write("ERREUR {}".format(e))
			return None

def dlfile(url,destFld):
	# Open the url
	try:
		f = urllib2.urlopen(url)
		print "downloading " + url
		CHUNK = 16*1024
		# Open our local file for writing
		with open(os.path.join(destFld, os.path.basename(url)), "wb") as local_file:
			while True:
				chunk = f.read(CHUNK)
				if not chunk:
					break
				local_file.write(chunk)

	#handle errors
	except urllib2.URLError as e:
		sys.stderr.write("{}".format(e))

def is_tool(name):
	try:
		devnull = open(os.devnull)
		subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
	except OSError as e:
		if e.errno == os.errno.ENOENT:
			return False
	return True

if	__name__ == '__main__':

	run(sys.argv[1:])

