#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import random
import subprocess
import numpy
import math
import sys
import time
import os
import argparse
import io
from Bio import SeqIO

#Globals variables
__version__="v1.3"
AppName = "RapidFastQC"
args = None
graphYsize = 10

def run(ArgsVal):

	global args
	
	description = """ A qualitative tool for quick quality check on fastq raw files. """

	parser = argparse.ArgumentParser(prog=AppName,description=description)

	parser.add_argument("-v","--verbose", help="show results on output", action="store_const", const="-v", default="")
	parser.add_argument("-d","--outputdir", help="output directory for report")
	parser.add_argument("-n","--numberOfReads", help="reads sampling. By default 2000",type=int)
	parser.add_argument("-i","--fastqfile", help="Fastq file(s) to analyse (fastq or gz file). For gz file, reads number is not available", nargs='*', required=True)
	args = parser.parse_args(ArgsVal)
	
	NbrTestReads = 2000
	
	if args.numberOfReads:
		NbrTestReads = args.numberOfReads
	
	if args.outputdir:
		outdir = args.outputdir
		if not os.path.exists(outdir):
			os.mkdir(outdir)
	else:
		outdir = ''
	
	NbrFastqFile = len(args.fastqfile)
	
	reportFileName = os.path.join(outdir,"RapidFastQC_REPORT_{}.txt".format(time.strftime('%d_%m_%Y_%Hh%Mm%S')))
	repfile = open(reportFileName,'w')
	
	
	head = "{} fastq file(s) input\n{} reads sample for each file\n\n".format(NbrFastqFile,NbrTestReads)
	
	if args.verbose:
		print head
		
	repfile.write(head)
	
	n = 1
	telapsed = 0
	str_estime = '...'
	
	for fastqFile in args.fastqfile:
		tstart = time.time()
		
		FastqBasename=os.path.basename(fastqFile)
	
		if args.verbose:
			print "\n-- Start analyse {} -- {}/{}  ete: {}".format(FastqBasename,n,NbrFastqFile,str_estime)
		
		if check_FileExtensions(fastqFile, ['fastq','fq','gz']):
		
			Rslt = GetFastqEstim(fastqFile,NbrTestReads,True)
			IsGzipFile = Rslt['gzipFile']
		
			if Rslt['meanReadNbr']==0:
				strEstReads = "Estimated reads numbers N/A"
			else:
				strEstReads = "Estimated reads numbers {}M ±{}".format(Rslt['meanReadNbr'],Rslt['ErrStdReadNbr'])
		
			if IsGzipFile: strEstReads += " (gzip)"
		
			telapsed += time.time() - tstart
		
			str_estime = time.strftime('%Hh%Mm%Ss', time.gmtime((telapsed/n)*(NbrFastqFile-n)))
		
			report = TXTFrame("File:{}\n{}\nEstimated Reads size {}pb ±{}\nMean quality = {}".format(os.path.basename(FastqBasename),strEstReads, Rslt['meanReadsSize'],Rslt['ErrStdReadsSize'], Rslt['meanqal']))
			report += "\n{}\n\n".format(Rslt['qalplot'])
		else:
		
			report = "ERROR: \'{}\' IS NOT A FASTQ FILE".format(fastqFile)
		
		if args.verbose:
				print report
		
		repfile.write(report)
		
		n += 1
		
	repfile.close()

def ExtractMeanQual(fastqVar):

	fastq = io.BytesIO(fastqVar)
	
	qal = []
	
	maxQalLen = 0
	
	for seq in SeqIO.parse(fastq,"fastq"):
		quality = seq.letter_annotations['phred_quality']
		if len(quality)>maxQalLen:
			maxQalLen = len(quality)
		qal.append(quality)
	
	for v in qal:
		if len(v)<maxQalLen:
			for i in range(0,maxQalLen-len(v)):
				v.append(0)
	
	qualA = numpy.array(qal)
	
	meanqal = numpy.mean(qualA, axis=0)

	meanqal = [((float(i)+float(j))/2) for i,j in zip(meanqal[0::2], meanqal[1::2])]

	gMean = numpy.mean(meanqal)
	
	plot = cmd_plot(meanqal,graphYsize)
	
	return [round(gMean,1), plot]

def cmd_plot(datas,ysize):
	minval = min(datas)
	maxval = max(datas)-minval
	tmp = [(d-minval)/maxval for d in datas]
	rslt = "^ Qual (min={} - max={})\n".format(round(min(datas),1),round(max(datas),1))
	valOnx=[]
	for i in range(ysize,-1,-1):
			line="|"
			for j,v in enumerate(tmp):
					if v>=(float(i)/float(ysize)) and j not in valOnx:
							line+="+"
							valOnx.append(j)
					else:
							line+=" "
			rslt += line+"\n"
	xline = ""

	x=0
	while x<(len(datas)+1):
		if x%10==0:
			xline += "|"
		else:
			xline += "-"
		x+=1
	rslt += xline+"-->\n"
	
	xLabelline = ""
	x=0
	while x<(len(datas)+1):
		if x%10==0:
			xLabelline += "{}".format(str(x*2))
			x+=len(str(x*2))-1
		else:
			xLabelline += " "
		x+=1
	
	rslt += xLabelline+" (bp)\n"
	return rslt
		
def GetFastqEstim(fastqFile, NbrTestReads,ShowprintProgress=False):

	random.seed(time.time())
	IsGzipFile = False
	FileSize = 0
	
	if is_gzip(fastqFile):
	
		IsGzipFile = True
		
		f = gzip.open(fastqFile, "r")
		
	else:
	
		f = open(fastqFile,"r")
		FileSize = GetFileSize(f)

	rslt = ''

	LineSize = []
	ReadsLen = []

	for i in range(0,NbrTestReads):
	
		sumLine = 0
		
		if IsGzipFile:
			for _ in range(0,randomInt_mult(0,9000,4)):
				if f.readline()=='':
					f.seek(0)
		else:
			f.seek(RandomFindHeadLine(f))
	
		Hline = f.readline()
		rslt += Hline
		Seqline = f.readline()
		rslt += Seqline
	
		ReadsLen.append(len(Seqline.rstrip()))
	
		sumLine += len(Hline)
		sumLine += len(Seqline)
	
		for l in range(0,2):
			line = f.readline()
			rslt += line
			sumLine += len(line)
	
		LineSize.append(FileSize/sumLine)
		
		if ShowprintProgress and args.verbose:
			printProgress("{} %".format(int(100*i/(NbrTestReads-1))))
	
	if args.verbose:
		print ""
	
	a = numpy.array(LineSize)
	b = numpy.array(ReadsLen)
	meanReadNbr = round((a.mean())/1000000,2)
	ErrStdReadNbr = round((numpy.std(a, ddof = 1)/math.sqrt(len(LineSize))/1000000)*2,2)
	meanReadsSize = round(b.mean(),2)
	ErrStdReadsSize = round((numpy.std(b, ddof = 1)/math.sqrt(len(ReadsLen)))*2,2)
	
	meanqal = ExtractMeanQual(rslt)
	
	return {'meanReadNbr': meanReadNbr, 'ErrStdReadNbr' : ErrStdReadNbr, 'meanReadsSize' : meanReadsSize, 'ErrStdReadsSize' : ErrStdReadsSize, 'FastqData' : rslt, 'meanqal' : meanqal[0], 'qalplot' : meanqal[1], 'gzipFile' : IsGzipFile}

def OffsetFirstCharLine(f,offset):

	if offset<0:
			return 0
	f.seek(offset)
	f.readline()
	
	return f.tell()
	
def EndOffset(f):
	tmp = f.tell()
	f.seek(0,2)
	endf = f.tell()
	
	while f.read(1)!="\n":
		endf -= 1
		f.seek(endf)
	
	f.seek(tmp)
	
	return endf-1

def RandomFindHeadLine(f):
	Retoffset = 0
	
	endf = EndOffset(f)

	headLine = False
	
	offset = random.randint(0,endf)

	loffset = OffsetFirstCharLine(f,offset)
	f.seek(loffset)

	while headLine == False:
		Retoffset = f.tell()
		if f.read(1) != '@':
			f.readline()
		else:
			f.seek(Retoffset)
			line = f.readline()
			if line.find(' ')>0:
				headLine = True
		
			
	return Retoffset
	

def GetFileSize(f):
	tmp = f.tell()
	
	f.seek(0,2)
	
	size = f.tell()
	
	f.seek(tmp)
	
	return size

def TXTFrame(text):
	updn = ''
	core = ''
	
	maxlength = max(len(i.replace('±','A')) for i in text.split('\n'))
	
	for line in text.split('\n'):
		core += "* {} ".format(line)
		core += " "*(maxlength - len(line.replace('±','A')))
		core += "*\n"
	
	updn = "*"*(maxlength+4)
	
	return "{}\n{}{}\n".format(updn,core,updn)

def printProgress(txt):
	sys.stdout.write("\r"+str(txt))
	sys.stdout.flush()

def find_Magic_Bytes(filename,magic_bytes):
	
	with open(filename) as infile:
		file_start = infile.read(len(magic_bytes))
		
	if file_start.startswith(magic_bytes):
		return True
		
	return False

def is_gzip(filename):

	magic_bytes = "\x1f\x8b\x08"
	
	return find_Magic_Bytes(filename,magic_bytes)

def is_b2z(filename):

	magic_bytes = "\x42\x5a\x68"
	
	return find_Magic_Bytes(filename,magic_bytes)

def randomInt_mult(imin,imax,mult):
	nbr = random.randint(imin,imax)
	
	while nbr%mult!=0:
		nbr = random.randint(imin,imax)
	
	return nbr

def check_FileExtensions(filename, extList):
	fileExt = filename[filename.rfind('.')+1:]
	
	if fileExt.upper() in [ext.upper() for ext in extList]:
		return True
	
	return False

if	__name__ == '__main__':

	run(sys.argv[1:])
