#!/usr/bin/env python
# -*- coding: utf-8 -*-

### Author: Mathieu Genete
### 2018.12.17

import gzip
import bz2
import argparse
import random
import math
import time
from datetime import datetime
from Bio.SeqUtils import GC
from Bio import SeqIO
from Bio import Seq
from itertools import product
import sys
import cPickle
import stat
import os
import yaml

__version__= "v1.20"

args = None
deltaTprogress = None

random.seed(time.time())

def stderr_print(txt):
	sys.stderr.write(str(txt)+"\n")

def CompSeqAlphabet(seq,alphabet):
	bases=list(set(seq))
	
	for b in bases:
		if b not in alphabet:
			return False
	
	return True

def get_ratioAmbigous(seq):
	ud = Seq.IUPAC.IUPACData.unambiguous_dna_letters
	nbAmb = sum(seq.count(b) for b in ud)
	ratioAmb = 1.0-(float(nbAmb)/float(len(seq)))
	return ratioAmb
	
def get_ambiguous_dna(seq):
	d = Seq.IUPAC.IUPACData.ambiguous_dna_values
	return list(map("".join, product(*map(d.get, seq))))

def get_ShannonEntropy(dna,alphabet="ATCG"):
        ShEntropy = 0
        Udna = dna.upper()
        for base in alphabet:
                freq = float(Udna.count(base))/float(len(Udna))
                if freq==0: freq = 0.00000001
                ShEntropy = ShEntropy + freq*math.log(freq,2)
        return -ShEntropy

def RatioMaxBaseOccurence(seq):
	bases=list(set(seq))
	bocc=[]

	for b in bases:
		bocc.append(seq.count(b))

	ratio=float(max(bocc))/float(len(seq))
	
	return ratio

def NbrSuccessifBase(base,seq):
	l = len(seq)
	while l>0:
		q=str(base)*l
		nbr = seq.count(q)
		if nbr>1:
			return l
		l-= 1		
	return 0

def filterkmers(kmerList,kmerExcludeList,minShEntropy=0.8,MaxratioAmb=0.2,alphabet='ATCGN'):
	return_kmers=set()
	nbrAmbigous = 0
	for k in kmerList:
		if get_ratioAmbigous(k)<MaxratioAmb and CompSeqAlphabet(k,alphabet) and get_ShannonEntropy(k)>=minShEntropy and k not in kmerExcludeList:
			if args.ambigousDNA:
				ambDNA = get_ambiguous_dna(k)
				if len(ambDNA)>1: nbrAmbigous += len(ambDNA)
				for ambSeq in ambDNA:
					return_kmers.add(ambSeq)
			else:
				return_kmers.add(k)
				
	return return_kmers, nbrAmbigous

def generateKmerDic(refFiles,excludeFiles,kmerFwdRev,kmerSize,alphabet,minShEntropy,MaxratioAmb):

	kmerList = set()
        kmerExcludeList = set()
	totalLen = 0
	nbrSeq = 0
	
	stderr_print(TitleFrame("construct Kmers dictionnary"))
        #construct kmer Exclude List
        if len(excludeFiles)>0:
                stderr_print(TitleFrame("\t=> construct Kmers exclude list"))
                for exfile in excludeFiles:
                        if os.path.exists(exfile):
                                exFileType="fasta"
                                if ".fastq" in exfile or ".fq" in exfile:
                                        exFileType="fastq"
                                for exRec in SeqIO.parse(exfile,exFileType):
                                        exSeq = str(exRec.seq).upper()
                                        exRcSeq = str(exRec.seq.reverse_complement())
                                        for i in range(0,len(exSeq)-kmerSize,1):
                                                if kmerFwdRev=="BOTH":
                                                        kmerExcludeList.add(exSeq[i:i+kmerSize].upper())
                                                        kmerExcludeList.add(exRcSeq[i:i+kmerSize].upper())
                                                elif kmerFwdRev=="FWD":
                                                        kmerExcludeList.add(exSeq[i:i+kmerSize].upper())
                                                elif kmerFwdRev=="REV":
                                                        kmerExcludeList.add(exRcSeq[i:i+kmerSize].upper())

                stderr_print(TitleFrame("\t=> {} Kmers to exclude".format(len(kmerExcludeList))))
        #construct kmer list
	for rfile in refFiles:
		if os.path.exists(rfile):
                        refFileType="fasta"
                        if ".fastq" in rfile or ".fq" in rfile:
                                refFileType="fastq"
			for record in SeqIO.parse(rfile,refFileType):

				seq = str(record.seq).upper()
				rcSeq = str(record.seq.reverse_complement())
				nbrSeq +=1
				totalLen += len(seq)
	
				for i in range(0,len(seq)-kmerSize,1):
                                        if kmerFwdRev=="BOTH":
                                                kmerList.add(seq[i:i+kmerSize].upper())
                                                kmerList.add(rcSeq[i:i+kmerSize].upper())
                                        elif kmerFwdRev=="FWD":
                                                kmerList.add(seq[i:i+kmerSize].upper())
                                        elif kmerFwdRev=="REV":
                                                kmerList.add(rcSeq[i:i+kmerSize].upper())
		else:
			sys.exit("ERROR: file '{}' not found".format(rfile))

	s_kmerList, nbrAmbigous = filterkmers(kmerList,kmerExcludeList,minShEntropy,MaxratioAmb,alphabet)
	ambPrint=""
	if nbrAmbigous>0:
		ambPrint=" - {} ambigous sequences generated".format(nbrAmbigous)
	stderr_print("{nbfiles} ref file(s) analysed\n{nbref} ref seq traited\n{nbkmer} kmers ({rmkmer} filtered{amb})\n".format(nbfiles=len(refFiles),nbref=nbrSeq,nbkmer=len(s_kmerList),rmkmer=(len(kmerList)-len(s_kmerList)+nbrAmbigous),amb=ambPrint))
	
	return s_kmerList

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

def creat_FastqReadsDictionary(fastqFile):
	rslt=set()
	fastq = open(fastqFile,'r')
	
	while True:
		l1 = fastq.readline()
		l2 = fastq.readline()
		l3 = fastq.readline()
		l4 = fastq.readline()
		seq = l2.strip()
		rslt.add(seq)
		
		if not l4: break
		
	fastq.close()
	
	return rslt

def printProgress(txt):
	sys.stderr.write("\r"+str(txt))
	sys.stderr.flush()

def kmerRefFilterProc(dicReads,kmerSize,fastqFile,minMatchKeepSeq,keepNotFiltered,namedPipe,kmerstats):

	start = time.time()
	
	stderr_print(TitleFrame("Launch reads filtering"))
	if namedPipe:
                stderr_print("Filtering on download stream...")
	
	bname = os.path.basename(fastqFile)

	outFastq = "{}_filtered.fastq".format(get_headFname(bname,['FASTQ','FQ']))
	keepoutFastq = "{}_KEEPED.fastq".format(get_headFname(bname,['FASTQ','FQ']))

	if args.outputdir:
		try:
			if not os.path.exists(args.outputdir):
				stderr_print("-- create {}".format(args.outputdir))
				os.mkdir(args.outputdir)
		except OSError as err:
			stderr_print(err)
			
		outFastq = os.path.join(args.outputdir,outFastq)
                keepoutFastq = os.path.join(args.outputdir,keepoutFastq)

	if args.append:
		writeMode="a"
	else:
		writeMode="w"

	if args.streamInput:
                fastq = sys.stdin
        else:
                fastq = open_fastqFile(fastqFile,namedPipe)

	outReadsDic=set()
	
	if os.path.exists(outFastq):
	
		outReadsDic = creat_FastqReadsDictionary(outFastq)
	
	outf = open(outFastq,writeMode)
	
        if keepNotFiltered:
                keepoutf = open(keepoutFastq,writeMode)
                
	tot = 0
	t=0
	matched=0
	appendNbr=0
        keepedNotFilt=0
        
        timeprint=time.time()

	while True:

		l1 = fastq.readline()
		l2 = fastq.readline()
		l3 = fastq.readline()
		l4 = fastq.readline()
		
		seq = l2.strip()
		
		if not l4: break
		
		m = 0
		
		
		if args.extremity5p3p:
			ext=args.extremity5p3p
			readRange = range(0,ext)+range(len(seq)-kmerSize-ext,len(seq)-kmerSize)
		else:
			readRange = range(0,len(seq)-kmerSize,1)

                readFound=False
		for i in readRange:
			km = seq[i:i+kmerSize].upper()
			if km in dicReads:
					m+=1
                                        add_kmerInStats(km,kmerstats)
			if m>=minMatchKeepSeq:
				matched+=1
				if (seq not in outReadsDic) or (not args.append):
                                        readFound=True
					appendNbr+=1
					outf.write(l1)
					outf.write(l2)
					outf.write(l3)
					outf.write(l4)
				break
		t += 1
		tot +=1

                if not readFound and keepNotFiltered:
                        keepedNotFilt += 1
                        keepoutf.write(l1)
                        keepoutf.write(l2)
                        keepoutf.write(l3)
                        keepoutf.write(l4)

                if args.maxreads>0 and t>=args.maxreads:
                        break

                if args.printprogress and time.time()-timeprint>deltaTprogress:
                        printProgress("{tr} reads analysed ({rps} reads per sec) -- {mr} reads kept ({pm}%) -- elapsed time {etime} min".format(etime=round(float(time.time()-start)/60.0,1),rps=round((float(tot)/float(time.time()-start)),1),tr=tot,mr=matched,pm=round(100*(float(matched)/float(tot)),3)))
                        timeprint = time.time()
		
	outf.close()

        if keepNotFiltered:
                keepoutf.close()

	if not args.streamInput:
                fastq.close()

        elapTime = time.strftime("%Hh%Mm%Ss", time.gmtime(time.time()-start))
	
	appendInfo=""
	notFilteredInfo=""

        if keepNotFiltered:
                notFilteredInfo="-- {} Not Filtered reads Kept ".format(keepedNotFilt)
	
	if args.append:
		appendInfo="-- {} reads append ".format(appendNbr)
		
	stderr_print("\n-- complete for {Fqfile}: {tot} total reads analysed {apnd}-- {pm}% reads kept ({match} reads) {notFilt}- in {tm}".format(Fqfile=os.path.basename(fastqFile),pm=round(100*(float(matched)/float(tot)),3),match=matched,tot=tot,tm=elapTime,apnd=appendInfo,notFilt=notFilteredInfo))
	
	return {'reads Nbr': tot,'reads kept': matched,'out file':[os.path.abspath(outFastq)]}

def kmerRefFilterPairedProc(dicReads,kmerSize,fastqFwdFile,fastqRevFile,minMatchKeepSeq,fwdRevChoice,keepNotFiltered,namedPipe,kmerstats):

	start = time.time()
	
	stderr_print(TitleFrame("Launch paired reads filtering"))
	if namedPipe:
                stderr_print("Filtering on download stream...")
	fwdbname = os.path.basename(fastqFwdFile)
        revbname = os.path.basename(fastqRevFile)

	outFwdFastq = "{}FWD_filtered.fastq".format(get_headFname(fwdbname,['FASTQ','FQ']))
	outRevFastq = "{}REV_filtered.fastq".format(get_headFname(revbname,['FASTQ','FQ']))

	outKeepFwdFastq = "{}FWD_KEEPED.fastq".format(get_headFname(fwdbname,['FASTQ','FQ']))
	outKeepRevFastq = "{}REV_KEEPED.fastq".format(get_headFname(revbname,['FASTQ','FQ']))

	
	if args.outputdir:
		try:
			if not os.path.exists(args.outputdir):
				stderr_print("-- create {}".format(args.outputdir))
				os.mkdir(args.outputdir)
		except OSError as err:
			stderr_print(err)
			
		outFwdFastq = os.path.join(args.outputdir,outFwdFastq)
                outRevFastq = os.path.join(args.outputdir,outRevFastq)

		outKeepFwdFastq = os.path.join(args.outputdir,outKeepFwdFastq)
                outKeepRevFastq = os.path.join(args.outputdir,outKeepRevFastq)

	if args.append:
		writeMode="a"
	else:
		writeMode="w"
	
        fwdfastq = open_fastqFile(fastqFwdFile,namedPipe)
        revfastq = open_fastqFile(fastqRevFile,namedPipe)

	outFwdReadsDic=set()
        outRevReadsDic=set()

	if os.path.exists(outFwdFastq) and os.path.exists(outRevFastq):
	
		outFwdReadsDic = creat_FastqReadsDictionary(outFwdFastq)
		outRevReadsDic = creat_FastqReadsDictionary(outRevFastq)
	
	Fwdoutf = open(outFwdFastq,writeMode)
        Revoutf = open(outRevFastq,writeMode)

	if keepNotFiltered:
                KeepFwdoutf = open(outKeepFwdFastq,writeMode)
                KeepRevoutf = open(outKeepRevFastq,writeMode)

	tot = 0
	t=0
	matched=0
	appendNbr=0
        keepedNotFilt=0
        FWDFound=0
        REVFound=0

        timeprint=time.time()

	while True:

		Fwdl1 = fwdfastq.readline()
		Fwdl2 = fwdfastq.readline()
		Fwdl3 = fwdfastq.readline()
		Fwdl4 = fwdfastq.readline()
		
		Fwdseq = Fwdl2.strip()
                
                Revl1 = revfastq.readline()
		Revl2 = revfastq.readline()
		Revl3 = revfastq.readline()
		Revl4 = revfastq.readline()
		
		Revseq = Revl2.strip()
		
                FwdL = [Fwdl1,Fwdl2,Fwdl3,Fwdl4]
                RevL = [Revl1,Revl2,Revl3,Revl4]

		if not Fwdl4 or not Revl4: break
		
		
		if args.extremity5p3p:
			ext=args.extremity5p3p
			FwdreadRange = range(0,ext)+range(len(Fwdseq)-kmerSize-ext,len(Fwdseq)-kmerSize)
                        RevreadRange = range(0,ext)+range(len(Revseq)-kmerSize-ext,len(Revseq)-kmerSize)
		else:
			FwdreadRange = range(0,len(Fwdseq)-kmerSize,1)
			RevreadRange = range(0,len(Revseq)-kmerSize,1)

                Rslt = {'RF':False}

                if fwdRevChoice in ['FWD','BOTH']:
                        Rslt = search_kmerInRead(FwdreadRange,kmerSize,dicReads,minMatchKeepSeq,outFwdReadsDic,Fwdseq,Fwdoutf,Revoutf,FwdL,RevL,kmerstats)
                        if Rslt['RF']:
                                FWDFound +=1

                if fwdRevChoice=='REV' or ( not Rslt['RF'] and fwdRevChoice=='BOTH'):
                        Rslt = search_kmerInRead(RevreadRange,kmerSize,dicReads,minMatchKeepSeq,outRevReadsDic,Revseq,Fwdoutf,Revoutf,FwdL,RevL,kmerstats)
                        if Rslt['RF']:
                                REVFound +=1

                matched += Rslt['matched']
                appendNbr += Rslt['appendNbr']

		t += 1
		tot +=1
		
                if keepNotFiltered:
                        keepVal = keep_Not_Filtered_PairedReads(Rslt['RF'],KeepFwdoutf,KeepRevoutf,FwdL,RevL)
                        if keepVal:
                                keepedNotFilt += 1

                if args.maxreads>0 and t>=args.maxreads:
                        break

                if args.printprogress and time.time()-timeprint>deltaTprogress:
                        printProgress("{tr} reads pair analysed ({rps} reads pair per sec) -- {mr} reads pair kept ({pm}%) -- elapsed time {etime} min".format(etime=round(float(time.time()-start)/60.0,1),rps=round((float(tot)/float(time.time()-start)),1),tr=tot,mr=matched,pm=round(100*(float(matched)/float(tot)),3)))
                        timeprint = time.time()

	Fwdoutf.close()
	Revoutf.close()
        fwdfastq.close()
        revfastq.close()

	if keepNotFiltered:
                KeepFwdoutf.close()
                KeepRevoutf.close()


        elapTime = time.strftime("%Hh%Mm%Ss", time.gmtime(time.time()-start))
	
	appendInfo=""
	notFilteredInfo=""

        if keepNotFiltered:
                notFilteredInfo="-- {} Not Filtered paired reads Kept ".format(keepedNotFilt)

	if args.append:
		appendInfo="-- {} reads append ".format(appendNbr)

        strFqFiles = "[{} - {}]".format(fwdbname,revbname)
	stderr_print("-- PAIRED complete for {Fqfile}: {tot} total reads per pair analysed {apnd}-- {pm}% reads per pair kept ({match} reads F:{fwdFnd} + R:{revFnd}) {notFilt}- in {tm}".format(Fqfile=strFqFiles,pm=round(100*(float(matched)/float(tot)),3),match=matched,tot=tot,tm=elapTime,apnd=appendInfo,notFilt=notFilteredInfo,fwdFnd=FWDFound,revFnd=REVFound))

	return {'reads Nbr': tot,'reads kept': matched,'out file':[os.path.abspath(outFwdFastq),os.path.abspath(outRevFastq)]}

def get_headFname(filename,fileExtList):
        for ext in fileExtList:
                if filename.upper().rfind(".{}".format(ext))>0:
                        return filename[0:filename.upper().rfind(".{}".format(ext))]
        return filename

def keep_Not_Filtered_PairedReads(ReadFoundValue,KeepFwdoutf,KeepRevoutf,FwdL,RevL):

        readKept=False

        if not ReadFoundValue:

                readKept=True

                for Fline in FwdL:
                        KeepFwdoutf.write(Fline)

                for Rline in RevL:
                        KeepRevoutf.write(Rline)

        return readKept

def search_kmerInRead(readRange,kmerSize,dicReads,minMatchKeepSeq,outReadsDic,seq,Fwdoutf,Revoutf,FwdL,RevL,kmerstats):
        m=0
        matched = 0
        appendNbr = 0
        Readfiltered = False

        for i in readRange:
                km = seq[i:i+kmerSize].upper()
                if km in dicReads:
                        m+=1
                        add_kmerInStats(km,kmerstats)

                if m>=minMatchKeepSeq:
                        matched+=1
                        if (seq not in outReadsDic) or (not args.append):
                                appendNbr+=1
                                Readfiltered = True
                                for Fline in FwdL:
                                        Fwdoutf.write(Fline)

                                for Rline in RevL:
                                        Revoutf.write(Rline)
                                
                        break

        outVal = {'RF':Readfiltered,'matched':matched,'appendNbr':appendNbr}

        return outVal

def add_kmerInStats(km,kmerstats):
        if km not in kmerstats.keys():
                kmerstats[km]={'count':0}
        kmerstats[km]['count']+=1

def create_PickleFromDict(pickleFilePath, dic):
	pf = gzip.open(pickleFilePath,'wb')
	cPickle.dump(dic,pf)
	pf.close()

def load_DicFromPickleFile(pickleFilePath):
	pf = gzip.open(pickleFilePath,'rb')
	retVal = cPickle.load(pf)
	pf.close()
	return retVal

def open_fastqFile(fastqFile,namedPipe):

        if  namedPipe:
                fastq = open(fastqFile)
        else:
                if is_gzip(fastqFile):
	
                        fastq = gzip.open(fastqFile, "r")
	
                elif is_b2z(fastqFile):
	
                        fastq = bz2.BZ2File(fastqFile, "r")
	
                else:
	
                        fastq = open(fastqFile,"r")

        return fastq

def TitleFrame(title):
	cnt = "* - {} - *".format(title)
	updn = "*"*len(cnt)
	
	return "\n{}\n{}\n{}".format(updn,cnt,updn)

def truncateStr(txt,maxlen):
	if len(txt)>maxlen:
		return "{} ... {}".format(txt[0:maxlen], txt[-1])
	else:
		return txt

def create_fifo(fifoname):
        if os.path.exists(fifoname):
                #isFifo = stat.S_ISFIFO(os.stat(fifoname).st_mode)
                os.remove(fifoname)
        os.mkfifo(fifoname)

def kmerRefFilter(ArgsVal):

	global args
	global deltaTprogress

        deltaTprogress = 3

	description = """ Tool for raw reads kmer-based filtering. """

	parser = argparse.ArgumentParser(prog="kmerRefFilter",description=description)
	
	parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

        parser.add_argument("-prog","--printprogress", help="print progession if set", action="store_const", const=True, default=False)
	parser.add_argument("-a","--append", help="if fastq filtered exists, add new reads at the end of file", action="store_const", const="-a", default="")
	parser.add_argument("-d","--ambigousDNA", help="generate all possible kmers from an ambigous kmer", action="store_const", const="-d", default="")
	parser.add_argument("-y","--yamlout", help="generate yaml output on stdout", action="store_const", const="-y", default="")
	parser.add_argument("-knf","--keepNotFiltered", help="keep not filtered reads in a file _KEEPED.fastq", action="store_const", const=True, default=False)
	
	parser.add_argument("-o","--outputdir", help="output directory for filtered fastq files")

	parser.add_argument("-k","--kmersize", help="size for kmers - default=20",default=20, type=int)
        parser.add_argument("-ks","--kmerstats", help="generate kmer statistics output table USED or ALL (USED by default)", nargs='?',const='USED')
	parser.add_argument("-m","--minMatchKeepSeq", help="minimum number of ref reads matched to keep fastq read -- default=1",default=1, type=int)
	parser.add_argument("-z","--minShEntropy", help="minimum Shannon Entropy for kmers (0 to 2.0)-- defaut=0.8",default=0.8, type=float)
	parser.add_argument("-q","--maxratioAmbigous", help="maximum ambigous bases accepted in kmers in %% -- defaut=0.2",default=0.2, type=float)
	
	parser.add_argument("-e","--extremity5p3p", help="compare 5' and 3' read extremity only -- nbr of bases to test from extremities default=1", nargs='?', const=1, type=int)
	
	parser.add_argument("-r","--referencesfasta", help="fasta files with references sequences", nargs='*')
	parser.add_argument("-i","--kmerinfile", help="load kmer dictionary from a saved dictionary file")
	parser.add_argument("-p","--kmeroutput", help="only save kmer dictionary in a file and exit without filtering")
	parser.add_argument("-x","--excludefasta", help="fasta files with sequences to exclude", nargs='*',default=[])
	parser.add_argument("-maxRD","--maxreads", help="maximum number of read to analyze",default=0, type=int)

	parser.add_argument("-f","--fastqfile", help="Fastq to filter (fastq, gz, bz2)", nargs='*')
	parser.add_argument("-l","--fastqlist", help="text file with all fastq (fastq, gz, bz2) path -- 1 per line")
        parser.add_argument("-1","--fwdfastq", help="forward Fastq to filter (fastq, gz, bz2) respect file order with reverse files", nargs='*')
        parser.add_argument("-2","--revfastq", help="reverse Fastq to filter (fastq, gz, bz2) respect file order with forward files", nargs='*')
        parser.add_argument("-u","--fastqurl", help="url to Fastq to filter (fastq, gz)", nargs='*')
        parser.add_argument("-u1","--fwdfastqurl", help="url to forward Fastq to filter (fastq, gz)")
        parser.add_argument("-u2","--revfastqurl", help="url to reverse Fastq to filter (fastq, gz)")
        parser.add_argument("-ugzip","--urlgzip", help="indicate that url's given in -u or in -u1 and -u2 are gz files", action="store_const", const=True, default=False)
	parser.add_argument("-s","--streamInput", help="use stdind as input - specify output filename")
	parser.add_argument("-c","--fwdRevChoice", help="for paired filtering: filtering on forward reads only (FWD), reverse reads only (REV) or both (BOTH - by default)", default='BOTH')
	parser.add_argument("-kc","--kmerFwdRev", help="kmders dictionary generation: use forward kmers only (FWD), reverse kmers only (REV) or both (BOTH - by default)", default='BOTH')
	
	args = parser.parse_args(ArgsVal)
	
        if not args.outputdir and not args.kmeroutput:
                stderr_print("No output directory defined (-o)\nexit...")
                sys.exit(0)

        if not args.referencesfasta and not args.kmerinfile:
                stderr_print("No reference fasta file or dictionary file set as input (-r or -i)\nexit...")
                sys.exit(0)

        if args.fwdRevChoice not in ['FWD','REV','BOTH']:
                stderr_print("-c (--fwdRevChoice) option not correctly set. Use \"FWD\", \"REV\" or \"BOTH\"\nexit...")
                sys.exit(0)

        if args.kmerFwdRev not in ['FWD','REV','BOTH']:
                stderr_print("-kc (--kmerFwdRev) option not correctly set. Use \"FWD\", \"REV\" or \"BOTH\"\nexit...")
                sys.exit(0)

        if args.kmerstats and args.kmerstats not in ['USED','ALL']:
                stderr_print("-ks (--kmerstats) option not correctly set. Use \"USED\" or \"ALL\"\nexit...")
                sys.exit(0)

	progVersion = "{pname} {ver}".format(pname=os.path.basename(sys.argv[0]),ver=__version__)
	
	stderr_print(TitleFrame(progVersion))
	
	refFiles = args.referencesfasta
        excludeFiles = args.excludefasta
	
        paired = False
        namedPipe = False

	if args.fastqfile:
		fastqFile = args.fastqfile
		
	elif args.fastqlist:
		flist=open(args.fastqlist,"r")
		fastqFile = [f.strip() for f in flist]
		flist.close()
        elif args.fwdfastq and args.revfastq:

                fastqFwdFile = args.fwdfastq
                fastqRevFile = args.revfastq

                if len(fastqFwdFile) != len(fastqRevFile):
                        stderr_print("forward and reverse file number are not identical")
                        sys.exit()
                paired = True

        elif args.streamInput:
                fastqFile = [args.streamInput]

        elif args.fwdfastqurl and not args.fastqurl:
                fastqFwdurl = args.fwdfastqurl
                fastqRevurl = args.revfastqurl

                fastqFwdFile = ["/tmp/{}".format(os.path.basename(fastqFwdurl).strip(".gz"))]
                fastqRevFile = ["/tmp/{}".format(os.path.basename(fastqRevurl).strip(".gz"))]
                create_fifo(fastqFwdFile[0])
                create_fifo(fastqRevFile[0])

                if args.urlgzip:
                        wget_cmd="wget -q -O - {url} | gzip -d -c > {fifo} &"
                else:
                        wget_cmd="wget -q -O - {url} > {fifo} &"

                cmd1=wget_cmd.format(fifo=fastqFwdFile[0],url=fastqFwdurl)
                cmd2=wget_cmd.format(fifo=fastqRevFile[0],url=fastqRevurl)

                os.system(cmd1)
                os.system(cmd2)

                paired = True
                namedPipe = True

        elif args.fastqurl and not args.fwdfastqurl:
                fastqFile=[]
                fastqListUrl = args.fastqurl
                for fastqFileurl in fastqListUrl:
                        fifo = "/tmp/{}".format(os.path.basename(fastqFileurl).strip(".gz"))
                        fastqFile.append(fifo)
                        create_fifo(fifo)

                        if args.urlgzip:
                                wget_cmd="wget -q -O - {url} | gzip -d -c > {fifo} &"
                        else:
                                wget_cmd="wget -q -O - {url} > {fifo} &"

                        cmd=wget_cmd.format(fifo=fifo,url=fastqFileurl)
                        os.system(cmd)

                namedPipe = True

        elif args.kmeroutput:
                stderr_print("no input fastq files required to save kmer dictionnary")

	else:
		stderr_print("no input fastq files")
		sys.exit()
		
	if args.ambigousDNA:
		alphabet=Seq.IUPAC.IUPACData.ambiguous_dna_letters
	else:
		alphabet = 'ATCGN'
	
	kmerSize = args.kmersize
	
	if (args.maxratioAmbigous*kmerSize)>8:
		args.maxratioAmbigous=8.0/float(kmerSize)

        kmerFwdRev = args.kmerFwdRev
	MaxratioAmb=args.maxratioAmbigous
	minMatchKeepSeq = args.minMatchKeepSeq
	minShEntropy = args.minShEntropy
	
	maxStrParamLen= 50
	parameters="\n".join(["- "+str(k)+": "+truncateStr(str(v),maxStrParamLen) for k,v in vars(args).items() if v])
	stderr_print("Parameters:\n{}\n".format(parameters))
	
        if not args.kmerinfile:
                dicReads = generateKmerDic(refFiles,excludeFiles,kmerFwdRev,kmerSize,alphabet,minShEntropy,MaxratioAmb)
        else:
                dicReads = load_DicFromPickleFile(args.kmerinfile)
                kmerSize = len(list(dicReads)[0])
                stderr_print(TitleFrame("Load Kmers dictionnary from file"))
                stderr_print("file: {}".format(args.kmerinfile))
                stderr_print("{} kmers -- kmers size = {}\n".format(len(dicReads),kmerSize))

        if args.kmeroutput:
                stderr_print("Create kmer dictionnary file : {}\n".format(args.kmeroutput))
                create_PickleFromDict(args.kmeroutput,dicReads)
                stderr_print("\nkmer dictionnary saved\nexit...")
                sys.exit(0)

	start = time.time()
	
	yamlResume={}
        kmerstats={}

	fasqFileNbr = 0

        if paired:
                fasqFileNbr = len(fastqFwdFile)
                for i in range(0,len(fastqFwdFile)):
                        yamlKey = "{} - {}".format(fastqFwdFile[i],fastqRevFile[i])
                        yamlResume[yamlKey] = kmerRefFilterPairedProc(dicReads,kmerSize,fastqFwdFile[i],fastqRevFile[i],minMatchKeepSeq,args.fwdRevChoice,args.keepNotFiltered,namedPipe,kmerstats)
                        if namedPipe:
                                os.remove(fastqFwdFile[i])
                                os.remove(fastqRevFile[i])
        else:
                fasqFileNbr = len(fastqFile)
                for fastq in fastqFile:
                        yamlResume[os.path.basename(fastq)] = kmerRefFilterProc(dicReads,kmerSize,fastq,minMatchKeepSeq,args.keepNotFiltered,namedPipe,kmerstats)
                        if namedPipe:
                                os.remove(fastq)
	
	if args.yamlout:
		stderr_print(TitleFrame("Yaml output on stdout"))
		yaml.dump(yamlResume, sys.stdout, default_flow_style=False)

        
        stderr_print(TitleFrame("kmer statistics generation"))
        #kmuse=sum([1 for km in dicReads if km in kmerstats.keys()])
        kmuse=len(kmerstats.keys())

        if args.kmerstats:
                dt=datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
                kmerstatsFileName=os.path.join(args.outputdir,"kmerstats_{}".format(dt))
                with open(kmerstatsFileName,'w') as kmfile:
                        for km in dicReads:
                                if km in kmerstats.keys():
                                        kmfile.write("{kmer}\t{cnt}\t{entro}\n".format(kmer=km,cnt=kmerstats[km]['count'],entro=get_ShannonEntropy(km)))
                                elif args.kmerstats=='ALL':
                                        kmfile.write("{kmer}\t{cnt}\t{entro}\n".format(kmer=km,cnt=0,entro=get_ShannonEntropy(km)))

        stderr_print("Total kmer: {} -- kmers used nbr {} -- kmers nevers used {}\n".format(len(dicReads),kmuse,len(dicReads)-kmuse))
        stderr_print("30 most used kmers:\n")
        tenSortedkmuse=sorted(kmerstats.iteritems(), key=lambda x : x[1], reverse=True)[:30]
        stderr_print("\tkmer\tcount\tentropy")
        for km,v in tenSortedkmuse:
                stderr_print("\t{kmer}\t{cnt}\t{entro}".format(kmer=km,cnt=v['count'],entro=get_ShannonEntropy(km)))
                                
                        
	elapTime = time.strftime("%Hh%Mm%Ss", time.gmtime(time.time()-start))
	
	stderr_print("\nFiltering terminated in {} - {} fastq file(s) analysed".format(elapTime,fasqFileNbr))

if	__name__ == '__main__':
	kmerRefFilter(sys.argv[1:])
