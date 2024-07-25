#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Author: Mathieu Genete
### 2018.12.17

import gzip
import bz2
import argparse
import random
import math
import numpy as np
import time
from datetime import datetime
from Bio import SeqIO
from Bio import Seq
from itertools import product
import sys
import pickle as cPickle
import stat
import os
import psutil
import yaml
from multiprocessing import Pool
from functools import partial

__version__= "v1.5"

args = None
logfile = None
deltaTprogress = None
readscount = 0

random.seed(time.time())

def find_read_ORF(seq,minPercORF,maxPercORF):
    outCoding={'IsCoding':False,'pctORF':0,'strand':{'value':'','phase':0,'prot':'','protShE':0,'protFreq':{'letters':'','freq':-1}},'dnaShE':0}
    if len(seq)>=6:
        theoricProtSize=(len(seq)-(len(seq)%3))/3
        dnaSeq=Seq.Seq(seq)
        RCdnaSeq=dnaSeq.reverse_complement()
        for phase in range(0,3):
            pstrand=dnaSeq[phase:]
            nstrand=RCdnaSeq[phase:]
            phased_pstrand=pstrand[:(len(pstrand)-(len(pstrand)%3))]
            phased_nstrand=nstrand[:(len(nstrand)-(len(nstrand)%3))]
            dnaPahsed={'+':phased_pstrand,'-':phased_nstrand}
            for strand,pdna in dnaPahsed.items():
                protSeq=str(pdna.translate())
                maxprotSize=float(max([len(subpseq) for subpseq in protSeq.replace("X","*").split("*")]))
                pctORF=maxprotSize/theoricProtSize
                s=str(pdna)
                if pctORF>=minPercORF and pctORF<=maxPercORF:
                    outCoding['IsCoding']=True
                    outCoding['pctORF']=pctORF
                    outCoding['strand']['value']=strand
                    outCoding['strand']['phase']=phase
                    outCoding['strand']['prot']=protSeq
                    outCoding['strand']['protFreq']=repeats_frequency(protSeq)
    return outCoding

def is_read_ORF_PAIRED(Fwdseq,Revseq,minPercORF,maxPercORF,maxprotrepeatfreq,sbcoding):
    FWDcoding=find_read_ORF(Fwdseq,minPercORF,maxPercORF)
    REVcoding=find_read_ORF(Revseq,minPercORF,maxPercORF)
    
    FwdFreqUnderThrld=(FWDcoding['strand']['protFreq']['freq']<=maxprotrepeatfreq and FWDcoding['strand']['protFreq']['freq']>0)
    RevFreqUnderThrld=(REVcoding['strand']['protFreq']['freq']<=maxprotrepeatfreq and REVcoding['strand']['protFreq']['freq']>0)

    if sbcoding=="BOTH":
        TotalpctORF=((FWDcoding['pctORF']+REVcoding['pctORF'])/2)
        TotalpctORFThrld=(TotalpctORF>=minPercORF)
        FwdRevFreqUnderThrld=(FwdFreqUnderThrld and RevFreqUnderThrld)
        FwdRevSameStrand=(FWDcoding['strand']['value']==REVcoding['strand']['value'] and test_strand(FWDcoding['strand']['value']))
        boolTest=(FWDcoding['IsCoding'] and REVcoding['IsCoding'] and FwdRevSameStrand and FwdRevFreqUnderThrld and TotalpctORFThrld)

    elif sbcoding=="SINGLE":
        boolTest=((FWDcoding['IsCoding'] and FwdFreqUnderThrld and test_strand(FWDcoding['strand']['value'])) or (REVcoding['IsCoding'] and RevFreqUnderThrld and test_strand(REVcoding['strand']['value'])))
    elif sbcoding=="SINGEXBOTH":
        boolTest=((FWDcoding['IsCoding']!=REVcoding['IsCoding']) and ((FWDcoding['IsCoding'] and FwdFreqUnderThrld and test_strand(FWDcoding['strand']['value'])) or (REVcoding['IsCoding'] and RevFreqUnderThrld and test_strand(REVcoding['strand']['value']))))
    else:
        boolTest=False
    return boolTest

def is_read_ORF_SINGLE(readseq,minPercORF,maxPercORF,maxprotrepeatfreq,sbcoding):
    Readcoding=find_read_ORF(readseq,minPercORF,maxPercORF)
    
    if Readcoding['IsCoding']:
        FreqUnderThrld=(Readcoding['strand']['protFreq']['freq']<=maxprotrepeatfreq and Readcoding['strand']['protFreq']['freq']>0)
        TotalpctORFThrld=(Readcoding['pctORF']>=minPercORF)
        boolTest=(FreqUnderThrld and TotalpctORFThrld and test_strand(Readcoding['strand']['value']))
    else:
        boolTest=False
    return boolTest

def test_strand(strandValue):
    strandintdic={'+':1,'-':2}
    strandint=strandintdic[strandValue]
    if args.strandorf!=0:
        return (args.strandorf==strandint)
    else:
        return True

def get_procmem():
    process = psutil.Process(os.getpid())
    return process.memory_info()

def check_max_memoryGo_used():
    memused=0
    if args.maxmemory>0:
            maxmem=args.maxmemory*1000000000
            memused=get_procmem().rss
            if  memused > maxmem:
                    stderr_print("ERROR: maximum memory allowed ({} Go) used.({})".format(args.maxmemory,memused))
                    sys.exit(-5)                
    return memused

def stderr_print(txt):
    if not args.noverbose:
        sys.stderr.write(str(txt)+"\n")
    log_write(txt)

def log_write(txt):
    if logfile and not args.nologfile :
        with open(logfile,'a') as hlogf:
            hlogf.write("{}\n".format(str(txt)))

def CompSeqAlphabet(seq,alphabet):
    s1=set(seq)
    s2=set(alphabet)
    if len(s1.difference(s2))>0:
            return False
    return True

def get_ratioAmbigous(seq):
    ud = ['G','A','T','C']
    nbAmb = sum(seq.count(b) for b in ud)
    ratioAmb = 1.0-(float(nbAmb)/float(len(seq)))
    return ratioAmb

def get_ShannonEntropy(dna,alphabet="ATCG"):
    ShEntropy = 0
    Udna = dna.upper()
    for base in alphabet:
            freq = float(Udna.count(base))/float(len(Udna))
            if freq>0:
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

def build_kmers(sequence, ksize,muted=False):
    kmers = set()
    n_kmers = len(sequence) - ksize + 1
    for i in range(n_kmers):
            kmer = sequence[i:i + ksize]
            if muted:
                    kmers.update(mutedkmers(kmer))
            else:
                    kmers.add(kmer)
    return kmers

def repeats_frequency(sequence,winmax=4):
    maxfreq={'letters':'','freq':0}
    for i in range(2,winmax+1):
            alphabet=set([sequence[s:s+i] for s in range(0,len(sequence)-(i-1))])
            for letters in alphabet:
                    freq=float(sequence.count(letters))/(float(len(sequence))/float(i))
                    if freq>maxfreq['freq']:
                            maxfreq['freq']=freq
                            maxfreq['letters']=letters
    return maxfreq

def filterkmers(kmerList,kmerExcludeList,minShEntropy=0.8,MaxratioAmb=0.2,maxRepeatsFreq=0.5,alphabet='ATCGN'):
    return_kmers=set()
    nbrAmbigous = 0
    rawsize=len(kmerList)
    dt=time.time()-6
    startclean=time.time()
    if args.multiproc>0:
            p = Pool(args.multiproc)
            test_kmerSeq_k=partial(test_kmerSeq,kmerExcludeList=kmerExcludeList,minShEntropy=minShEntropy,MaxratioAmb=MaxratioAmb,maxRepeatsFreq=maxRepeatsFreq,alphabet=alphabet)
            outpoutmp=[]
            psync=p.map_async(test_kmerSeq_k,kmerList,callback=outpoutmp.extend)
            p.close()
            while not psync.ready():
                    time.sleep(1)
                    check_max_memoryGo_used()
                    if args.printprogress:
                        remaining=100*float(abs(rawsize-(psync._number_left * psync._chunksize)))/float(rawsize)
                        printProgress("Cleaning progress {}%".format(round(remaining,1)))
            if None in outpoutmp: outpoutmp.remove(None)
            return_kmers=set(outpoutmp)
    else:
            
            for c,k in enumerate(kmerList):
                    return_kmers.add(test_kmerSeq(k,kmerExcludeList,minShEntropy,MaxratioAmb,maxRepeatsFreq,alphabet))

                    if args.printprogress and time.time()-dt>5:
                            check_max_memoryGo_used()
                            pcomplete=100*(float(c+1)/rawsize)
                            printProgress("Cleaning progress {}%".format(round(pcomplete,1)))
                            dt=time.time()
    printProgress("Cleaning progress 100.0% in {} sec".format(round(time.time()-startclean,2)))
            
    stderr_print("\n")
    return return_kmers, nbrAmbigous

def test_kmerSeq(k,kmerExcludeList,minShEntropy,MaxratioAmb,maxRepeatsFreq,alphabet):
    if k not in kmerExcludeList and CompSeqAlphabet(k,alphabet) and get_ratioAmbigous(k)<MaxratioAmb and get_ShannonEntropy(k)>=minShEntropy and repeats_frequency(k)['freq']<=maxRepeatsFreq:
        return k

def generateKmerDic(refFiles,excludeFiles,kmerFwdRev,kmerSize,alphabet,minShEntropy,MaxratioAmb,mutedkmers,maxRepeatsFreq):

    kmerList = set()
    kmerExcludeList = set()
    totalLen = 0
    nbrSeq = 0
    dt=time.time()-6
    stderr_print(TitleFrame("construct Kmers dictionnary"))
    startbuild=time.time()
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
                    #Check memomry used every 5 sec
                    if time.time()-dt>=5:
                        check_max_memoryGo_used()
                        dt=time.time()
                    if kmerFwdRev=="BOTH":
                        kmerExcludeList.update(build_kmers(exSeq,kmerSize))
                        kmerExcludeList.update(build_kmers(exRcSeq,kmerSize))
                    elif kmerFwdRev=="FWD":
                        kmerExcludeList.update(build_kmers(exSeq,kmerSize))
                    elif kmerFwdRev=="REV":
                        kmerExcludeList.update(build_kmers(exRcSeq,kmerSize))

        stderr_print(TitleFrame("\t=> {} Kmers to exclude".format(len(kmerExcludeList))))
    #construct kmer list
    for rfile in refFiles:
        if os.path.exists(rfile):
            refFileType="fasta"
            if ".fastq" in rfile or ".fq" in rfile:
                refFileType="fastq"
            for record in SeqIO.parse(rfile,refFileType):

                seq = str(record.seq).upper()
                rcSeq = str(record.seq.reverse_complement()).upper()
                nbrSeq +=1
                totalLen += len(seq)
                #Check memomry used every 5 sec and if set print progess
                if time.time()-dt>=5:
                    check_max_memoryGo_used()
                    dt=time.time()
                    if args.printprogress:
                            printProgress("constructing raw kmer dictionary ({} kmers)".format(len(kmerList)))
                if kmerFwdRev=="BOTH":
                    kmerList.update(build_kmers(seq,kmerSize,mutedkmers))
                    kmerList.update(build_kmers(rcSeq,kmerSize,mutedkmers))
                elif kmerFwdRev=="FWD":
                    kmerList.update(build_kmers(seq,kmerSize,mutedkmers))
                elif kmerFwdRev=="REV":
                    kmerList.update(build_kmers(rcSeq,kmerSize,mutedkmers))

        else:
            sys.exit("ERROR: file '{}' not found".format(rfile))

    memused=check_max_memoryGo_used()
    stderr_print("\n{} raw kmers generated in {} sec -- (memory usage {} Mo)".format(len(kmerList),round(time.time()-startbuild,2),round(float(memused)/1000000,2)))
    if mutedkmers:
        stderr_print("add all one base divergent kmers -- {}% of divergence for kmers size of {}".format(round(100.0/float(kmerSize),2),kmerSize))
    if args.nokmerclean:
        stderr_print("Disable kmers cleaning")
        s_kmerList=kmerList
    else:
        stderr_print("Launch kmers cleaning...")
        s_kmerList, nbrAmbigous = filterkmers(kmerList,kmerExcludeList,minShEntropy,MaxratioAmb,maxRepeatsFreq,alphabet)
        ambPrint=""
        if nbrAmbigous>0:
            ambPrint=" - {} ambigous sequences generated".format(nbrAmbigous)
        stderr_print("{nbfiles} ref file(s) analysed\n{nbref} ref seq traited\n{nbkmer} kmers ({rmkmer} filtered{amb})\n".format(nbfiles=len(refFiles),nbref=nbrSeq,nbkmer=len(s_kmerList),rmkmer=(len(kmerList)-len(s_kmerList)+nbrAmbigous),amb=ambPrint))
    
    return s_kmerList

def mutedkmers(kmseq):
    bases = 'ATCG'
    splits = [(kmseq[:i], kmseq[i:]) for i in range(len(kmseq) + 1)]
    replaces = [L + c + R[1:] for L, R in splits if R for c in bases]
    return set(replaces)

def find_Magic_Bytes(filename,magic_bytes):
    with open(filename,'r',encoding="ISO-8859-1") as infile:
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
    if not args.noverbose:
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
    SumReadsSize=0
    timeprint=time.time()

    while True:

        l1 = fastq.readline()
        l2 = fastq.readline()
        l3 = fastq.readline()
        l4 = fastq.readline()
        
        seq = l2.strip()

        kmseq=build_kmers(seq, kmerSize)

        if not l4: break
        
        m = 0
        
        
        if args.extremity5p3p:
            ext=args.extremity5p3p
            readRange = range(0,ext)+range(len(seq)-kmerSize-ext,len(seq)-kmerSize)
        else:
            readRange = range(0,len(seq)-kmerSize,1)

        readFound=False
        intersections = kmseq.intersection(dicReads)
        m=len(intersections)

        if args.minporf>0 or args.maxporf<1:
            keepReadBool=(m>0 and m>=minMatchKeepSeq and is_read_ORF_SINGLE(seq,args.minporf,args.maxprotfreq,args.maxporf,args.singbothcoding))
        else:
            keepReadBool=(m>0 and m>=minMatchKeepSeq)

        if keepReadBool and ((seq not in outReadsDic) or (not args.append)):
            matched+=1
            add_kmerInStats(intersections,kmerstats)
            readFound=True
            appendNbr+=1
            SumReadsSize+=len(seq)
            outf.write(l1)
            outf.write(l2)
            outf.write(l3)
            outf.write(l4)

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
        
    if args.minporf>0 or args.maxporf<1:
        singbothTXT={'SINGLE':'(Fwd OR Rev)','BOTH':'(Fwd AND Rev)','SINGEXBOTH':'(Fwd OR Rev exclude Fwd AND Rev)'}
        appendInfo+="-- {} contains at least {}% to {}% ORF ".format(singbothTXT[args.singbothcoding],args.minporf*100,args.maxporf*100)

    if matched>0:
        meanRdSize=round(float(SumReadsSize)/float(matched),1)
    else:
        meanRdSize=0

    stderr_print("\r-- complete for {Fqfile}: {tot} total reads analysed {apnd}-- {pm}% reads kept ({match} reads -- mean size {rdsize} bp) {notFilt}- in {tm}".format(Fqfile=os.path.basename(fastqFile),pm=round(100*(float(matched)/float(tot)),3),match=matched,tot=tot,tm=elapTime,apnd=appendInfo,notFilt=notFilteredInfo,rdsize=meanRdSize))
    
    return {'reads Nbr': tot,'reads kept': matched,'out file':[os.path.abspath(outFastq)]}

def kmerRefFilterPairedProc(dicReads,kmerSize,fastqConcatPairedFile,fastqFwdFile,fastqRevFile,minMatchKeepSeq,fwdRevChoice,keepNotFiltered,namedPipe,kmerstats):

    start = time.time()
    concatPairedFQ=(fastqConcatPairedFile!="")

    if concatPairedFQ:
        txttitle="Launch paired concatenated reads filtering"
    else:
        txttitle="Launch paired reads filtering"

    stderr_print(TitleFrame(txttitle))


    if namedPipe:
        stderr_print("Filtering on download stream...")

    if concatPairedFQ:
        fwdbname = os.path.basename(fastqConcatPairedFile)
        revbname = os.path.basename(fastqConcatPairedFile)
    else:
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
    
    if concatPairedFQ:
        fastqCP = open_fastqFile(fastqConcatPairedFile,namedPipe)
    else:
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
    SumReadsSize=0

    timeprint=time.time()

    while True:
        if concatPairedFQ:
            FwdL,RevL=read_ConcatenatedPairedFastqLine(fastqCP)
        else:
            FwdL,RevL=read_pairedFastqLine(fwdfastq,revfastq)
        Fwdseq = FwdL[1].strip()
        Revseq = RevL[1].strip()
        

        if not FwdL[3] or not RevL[3]: break
        
        
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
        SumReadsSize += Rslt['meanSize']

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

    if concatPairedFQ:
        fastqCP.close()
    else:
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

    if args.minporf>0 or args.maxporf<1:
        singbothTXT={'SINGLE':'(Fwd OR Rev)','BOTH':'(Fwd AND Rev)','SINGEXBOTH':'(Fwd OR Rev exclude Fwd AND Rev)'}
        appendInfo+="-- {} contains at least {}% to {}% ORF ".format(singbothTXT[args.singbothcoding],args.minporf*100,args.maxporf*100)

    if concatPairedFQ:
        strFqFiles = "[{} - concatenated ]".format(fwdbname)
    else:
        strFqFiles = "[{} - {}]".format(fwdbname,revbname)

    if matched>0:
        meanRdSize=round(float(SumReadsSize)/float(matched),1)
    else:
        meanRdSize=0

    stderr_print("\r-- PAIRED complete for {Fqfile}: {tot} total reads per pair analysed {apnd}-- {pm}% reads per pair kept ({match} reads F:{fwdFnd} + R:{revFnd} -- mean size {rdsize} bp) {notFilt}- in {tm}".format(Fqfile=strFqFiles,pm=round(100*(float(matched)/float(tot)),3),match=matched,tot=tot,tm=elapTime,apnd=appendInfo,notFilt=notFilteredInfo,fwdFnd=FWDFound,revFnd=REVFound,rdsize=meanRdSize))

    return {'reads Nbr': tot,'reads kept': matched,'out file':[os.path.abspath(outFwdFastq),os.path.abspath(outRevFastq)]}

def read_pairedFastqLine(hfwd,hrev):
    FwdLines=[hfwd.readline() for i in range(0,4)]
    RevLines=[hrev.readline() for i in range(0,4)]
    return FwdLines,RevLines

def read_ConcatenatedPairedFastqLine(hcp):
    FwdLines=[]
    RevLines=[]
    HcpLines=[hcp.readline() for i in range(0,4)]
    seq=HcpLines[1].strip()
    qual=HcpLines[3].strip()
    FwdLines.append(HcpLines[0])
    RevLines.append(HcpLines[0])
    FwdLines.append("{}\n".format(seq[:len(seq)/2]))
    RevLines.append("{}\n".format(seq[len(seq)/2:]))
    FwdLines.append(HcpLines[2])
    RevLines.append(HcpLines[2])
    if not qual:
        FwdLines.append("")
        RevLines.append("")
    else:
        FwdLines.append("{}\n".format(qual[:len(qual)/2]))
        RevLines.append("{}\n".format(qual[len(qual)/2:]))
    return FwdLines,RevLines

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
    SumReadsSize=0
    Readfiltered = False

    kmseq=build_kmers(seq, kmerSize)
    intersections = kmseq.intersection(dicReads)
    m=len(intersections)

    if args.minporf>0 or args.maxporf<1:
        keepReadBool=(m>0 and m>=minMatchKeepSeq and is_read_ORF_PAIRED(FwdL[1].strip(),RevL[1].strip(),args.minporf,args.maxporf,args.maxprotfreq,args.singbothcoding))
    else:
        keepReadBool=(m>0 and m>=minMatchKeepSeq)

    if keepReadBool and ((seq not in outReadsDic) or (not args.append)):
        matched+=1
        add_kmerInStats(intersections,kmerstats)
        appendNbr+=1
        Readfiltered = True
        SumReadsSize+=(float(len(FwdL[1].strip()))+float(len(RevL[1].strip())))/2
        for Fline in FwdL:
            Fwdoutf.write(Fline)
            
        for Rline in RevL:
            Revoutf.write(Rline)
    if matched>0:
        rdmeanSize=float(SumReadsSize)/float(matched)
    else:
        rdmeanSize=0

    outVal = {'RF':Readfiltered,'matched':matched,'appendNbr':appendNbr,'meanSize':rdmeanSize}

    return outVal

def add_kmerInStats(kmers,kmerstats):
    global readscount
    for km in kmers:
        if km not in kmerstats.keys():
            kmerstats[km]={'count':0,'idxreads':[]}
        kmerstats[km]['count']+=1
        kmerstats[km]['idxreads'].append(readscount)
    readscount+=1

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
            fastq = gzip.open(fastqFile, 'rt', encoding='utf-8')
        elif is_b2z(fastqFile):    
            fastq = bz2.open(fastqFile, 'rt', encoding='utf-8')    
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

def test_arg_range(val,allowedRange,argname):
    if  (allowedRange['min']!=None and val<allowedRange['min']) or (allowedRange['max']!=None and val>allowedRange['max']):
        stderr_print("ERROR: argument value out of range for {} : allowed range: [{},{}]".format(argname,allowedRange['min'],allowedRange['max']))
        return False
    return True

def clean_filteredFQ(yamlResume,kmerstats,thresholdMinMax,paired):
    yamlvalues=list(yamlResume.values())[0]
    outfiles=yamlvalues['out file']
    retainedkmerstats={}
    kmidxList_toremove=set()
    kmidxList_outside=set()
    for km,kmvalues in kmerstats.items():
        if (kmvalues['count']<thresholdMinMax[0] and thresholdMinMax[0]>=0) or (kmvalues['count']>thresholdMinMax[1] and thresholdMinMax[1]>0):
            kmidxList_toremove.update(set(kmvalues['idxreads']))
        else:
            kmidxList_outside.update(set(kmvalues['idxreads']))
            retainedkmerstats[km]=kmvalues
    kmidxList=kmidxList_toremove.difference(kmidxList_outside)
    rdcount=0
    keepreads=0
    if paired:
        outclean=["{}_cleaned.fastq".format(v[:v.rfind('.')]) for v in outfiles]
        with open(outclean[0],'w') as outClFwd:
            with open(outclean[1],'w') as outClRev:
                with open(outfiles[0],'r') as fwdFile:
                    with open(outfiles[1],'r') as revFile:
                        while True:
                            fwdLines,RevLines = read_pairedFastqLine(fwdFile,fwdFile)
                            if fwdLines[1]=="": break
                            if rdcount not in kmidxList:
                                keepreads+=1
                                write_lines(outClFwd,fwdLines)
                                write_lines(outClRev,RevLines)
                            rdcount+=1
    else:
        outclean="{}_cleaned.fastq".format(outfiles[0][:outfiles[0].rfind('.')])
        with open(outclean,'w') as outCl:
            with open(outfiles[0],'r') as fqfile:
                while True:
                    fqlines=[fqfile.readline() for i in xrange(4)]
                    if fqlines[1]=="" : break
                    if rdcount not in kmidxList:
                        keepreads+=1
                        write_lines(outCl,fqlines)
                    rdcount+=1
    total_reads=yamlvalues['reads Nbr']
    percent_retained=100*float(keepreads)/float(total_reads)
    return keepreads,total_reads,percent_retained,retainedkmerstats

def write_lines(fh,lines):
    for l in lines:
        fh.write(l)

def kmerRefFilter(ArgsVal):

    global args
    global deltaTprogress
    global logfile
    global readscount

    deltaTprogress = 3

    description = """ Tool for raw reads kmer-based filtering. """

    parser = argparse.ArgumentParser(prog="kmerRefFilter",description=description)
    
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

    parser.add_argument("-prog","--printprogress", help="print progession if set", action="store_const", const=True, default=False)
    parser.add_argument("-nc","--nokmerclean", help="Disable kmers cleaning step", action="store_const", const=True, default=False)
    parser.add_argument("-a","--append", help="if fastq filtered exists, add new reads at the end of file", action="store_const", const="-a", default="")
    parser.add_argument("-y","--yamlout", help="generate yaml output on stdout", action="store_const", const="-y", default="")
    parser.add_argument("-knf","--keepNotFiltered", help="keep not filtered reads in a file _KEEPED.fastq", action="store_const", const=True, default=False)

    parser.add_argument("-mem","--maxmemory", help="maximum memory used in Go (default = 5 Go) - set to 0 for unlimited memory usage",default=5, type=int)
    parser.add_argument("-P","--multiproc", help="use multiprocessing for kmer dictionary generation -- number of CPU",default=0, type=int)
    parser.add_argument("-o","--outputdir", help="output directory for filtered fastq files")
    parser.add_argument("-nolog","--nologfile", help="Disable log file output", action="store_const", const=True, default=False)
    parser.add_argument("-nov","--noverbose", help="Disable output log on stderr", action="store_const", const=True, default=False)

    parser.add_argument("-k","--kmersize", help="size for kmers - default=20",default=20, type=int)
    parser.add_argument("-ks","--kmerstats", help="generate kmer statistics output table USED or ALL (USED by default)", nargs='?',const='USED')
    parser.add_argument("-m","--minMatchKeepSeq", help="minimum number of ref reads matched to keep fastq read -- default=1",default=1, type=int)
    parser.add_argument("-z","--minShEntropy", help="minimum Shannon Entropy for kmers (0 to 2.0)-- defaut=0.8",default=0.8, type=float)
    parser.add_argument("-q","--maxratioAmbigous", help="maximum frequency of ambigous bases accepted in kmers -- defaut=0.2",default=0.2, type=float)
    parser.add_argument("-j","--maxRepeatsFreq", help="maximum frequency of 2-4 mers repetitions accepted in kmers -- defaut=1.0",default=1.0, type=float)
    parser.add_argument("-b","--mutedkmers", help="generates all single bases replacement possibilities for each dictionary kmers (huge memory consumption)", action="store_const", const=True, default=False)    
    parser.add_argument("-e","--extremity5p3p", help="compare 5' and 3' read extremity only -- nbr of bases to test from extremities default=1", nargs='?', const=1, type=int)
    
    parser.add_argument("-r","--referencesfasta", help="fasta files with references sequences", nargs='*')
    parser.add_argument("-i","--kmerinfile", help="load kmer dictionary from a saved dictionary file")
    parser.add_argument("-p","--kmeroutput", help="only save kmer dictionary in a file and exit without filtering")
    parser.add_argument("-x","--excludefasta", help="fasta files with sequences to exclude", nargs='*',default=[])
    parser.add_argument("-maxRD","--maxreads", help="maximum number of read to analyze",default=0, type=int)

    parser.add_argument("-f","--fastqfile", help="Fastq to filter (fastq, gz, bz2)")
    parser.add_argument("-ct","--concatenatedfastq", help="input fastq files are concatenated (FWD+REV in one read). Split half reads and quality values. Should be use with -f or -l or -u or -s options only", action="store_const", const=True, default=False)
    parser.add_argument("-1","--fwdfastq", help="forward Fastq to filter (fastq, gz, bz2) respect file order with reverse files")
    parser.add_argument("-2","--revfastq", help="reverse Fastq to filter (fastq, gz, bz2) respect file order with forward files")
    parser.add_argument("-u","--fastqurl", help="url to Fastq to filter (fastq, gz)")
    parser.add_argument("-u1","--fwdfastqurl", help="url to forward Fastq to filter (fastq, gz)")
    parser.add_argument("-u2","--revfastqurl", help="url to reverse Fastq to filter (fastq, gz)")
    parser.add_argument("-ugzip","--urlgzip", help="indicate that url's given in -u or in -u1 and -u2 are gz files", action="store_const", const=True, default=False)
    parser.add_argument("-s","--streamInput", help="use stdind as input - specify output filename")
    parser.add_argument("-c","--fwdRevChoice", help="for paired filtering: filtering on forward reads only (FWD), reverse reads only (REV) or both (BOTH - by default)", default='BOTH')
    parser.add_argument("-kc","--kmerFwdRev", help="kmders dictionary generation: use forward kmers only (FWD), reverse kmers only (REV) or both (BOTH - by default)", default='BOTH')
    
    parser.add_argument("-morf","--minporf", help="minimum percent of ORF in read",default=0, type=float)
    parser.add_argument("-Morf","--maxporf", help="maximum percent of ORF in read",default=1, type=float)

    parser.add_argument("-mkc","--minkmercount", help="minimum kmer count for recruited reads",default=0, type=int)
    parser.add_argument("-Mkc","--maxkmercount", help="maximum kmer count for recruited reads",default=0, type=int)

    parser.add_argument("-sorf","--strandorf", help="force strand for ORF filtering (both=0, +=1, -=2) -- default=0", default=0, type=int)
    parser.add_argument("-pfreq","--maxprotfreq", help="maximum frequency of AA repetitions in proteics sequences from reads (default = 0.6)",default=0.8, type=float)
    parser.add_argument("-sbc","--singbothcoding", help="apply ORF threshold for single (SINGLE); both (BOTH) forward and reverse reads to keep the pair or single exclude both (SINGEXBOTH) (SINGLE - by default)", default='SINGLE')

    args = parser.parse_args(ArgsVal)
    
    if not args.outputdir and not args.kmeroutput:
        stderr_print("No output directory defined (-o)\nexit...")
        sys.exit(1)

    if not args.referencesfasta and not args.kmerinfile:
        stderr_print("No reference fasta file or dictionary file set as input (-r or -i)\nexit...")
        sys.exit(2)

    if args.fwdRevChoice not in ['FWD','REV','BOTH']:
        stderr_print("-c (--fwdRevChoice) option not correctly set. Use \"FWD\", \"REV\" or \"BOTH\"\nexit...")
        sys.exit(3)

    if args.kmerFwdRev not in ['FWD','REV','BOTH']:
        stderr_print("-kc (--kmerFwdRev) option not correctly set. Use \"FWD\", \"REV\" or \"BOTH\"\nexit...")
        sys.exit(4)

    if args.kmerstats and args.kmerstats not in ['USED','ALL']:
        stderr_print("-ks (--kmerstats) option not correctly set. Use \"USED\" or \"ALL\"\nexit...")
        sys.exit(5)

    if args.singbothcoding not in ['SINGLE','BOTH','SINGEXBOTH']:
        stderr_print("-sbc (--singbothcoding) option not correctly set. Use \"SINGLE\" or \"BOTH\"\nexit...")
        sys.exit(6)

    if (args.fwdfastq and not args.revfastq) or (not args.fwdfastq and args.revfastq):
        stderr_print("-1 or -2 fastq files are missing\nexit...")
        sys.exit(7)

    if (args.fwdfastqurl and not args.revfastqurl) or (not args.fwdfastqurl and args.revfastqurl):
        stderr_print("-u1 or -u2 fastq urls are missing\nexit...")
        sys.exit(8)
    
    if args.concatenatedfastq and not args.fastqfile and not args.fastqurl and not args.streamInput:
        stderr_print("-ct option must be associated to one of these options: -f or -l or -u or -s\nexit...")
        sys.exit(9)

    #test arguments values:
    argToTest=[(args.maxmemory,{'min':0,'max':None},"--maxmemory"),\
               (args.kmersize,{'min':5,'max':100},"--kmersize"),\
               (args.minMatchKeepSeq,{'min':1,'max':None},"--minMatchKeepSeq"),\
               (args.minShEntropy,{'min':0.0,'max':2.0},"--minShEntropy"),\
               (args.maxratioAmbigous,{'min':0.0,'max':1.0},"--maxratioAmbigous"),\
               (args.maxRepeatsFreq,{'min':0.0,'max':1.0},"--maxRepeatsFreq"),\
               (args.multiproc,{'min':0,'max':None},"--multiproc"),\
               (args.minporf,{'min':0.0,'max':1.0},"--minporf"),\
               (args.maxporf,{'min':0.0,'max':1.0},"--maxporf"),\
               (args.maxprotfreq,{'min':0.0,'max':1.0},"--maxprotfreq"),\
               (args.strandorf,{'min':0,'max':2},"--strandorf"),\
               (args.minkmercount,{'min':0,'max':None},"--minkmercount"),\
               (args.maxkmercount,{'min':0,'max':None},"--maxkmercount")]

    kmerCountThrld=(args.minkmercount,args.maxkmercount)        
    testFailed=False
    for argt in argToTest:
        testrslt=test_arg_range(argt[0],argt[1],argt[2])
        if not testrslt:
            testFailed=True
    if testFailed:
        sys.exit(10)

    if args.maxporf<args.minporf:
        stderr_print("--minporf should be lesser than --maxporf")
        sys.exit(11)

    progVersion = "{pname} {ver} {dat}".format(pname=os.path.basename(sys.argv[0]),ver=__version__,dat=datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"))
    
    logfile=os.path.join(args.outputdir,"log_kmerRefFilter_{}".format(datetime.now().strftime("%d_%m_%Y-%H_%M_%S")))

    stderr_print(TitleFrame(progVersion))
    
    refFiles = args.referencesfasta
    excludeFiles = args.excludefasta

    paired = False
    concatenated=False
    namedPipe = False

    if args.fastqfile:
        fastqFile = args.fastqfile
        if args.concatenatedfastq:
            paired = True
            concatenated=True
    elif args.fwdfastq and args.revfastq:

        fastqFwdFile = args.fwdfastq
        fastqRevFile = args.revfastq

        if len(fastqFwdFile) != len(fastqRevFile):
            stderr_print("forward and reverse file number are not identical")
            sys.exit()
        paired = True

    elif args.streamInput:
        fastqFile = [args.streamInput]
        if args.concatenatedfastq:
            paired = True
            concatenated=True

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
        if args.concatenatedfastq:
            paired = True
            concatenated=True

    elif args.kmeroutput:
        stderr_print("no input fastq files required to save kmer dictionnary")

    else:
        stderr_print("no input fastq files")
        sys.exit()
        
    alphabet = 'ATCGN'
    kmerSize = args.kmersize
    
    if (args.maxratioAmbigous*kmerSize)>8:
        args.maxratioAmbigous=8.0/float(kmerSize)

    kmerFwdRev = args.kmerFwdRev
    MaxratioAmb=args.maxratioAmbigous
    minMatchKeepSeq = args.minMatchKeepSeq
    minShEntropy = args.minShEntropy
    maxRepeatsFreq = args.maxRepeatsFreq
    mutedkmers = args.mutedkmers

    maxStrParamLen= 50
    parameters="\n".join(["- "+str(k)+": "+truncateStr(str(v),maxStrParamLen) for k,v in vars(args).items() if v])
    stderr_print("Parameters:\n{}\n".format(parameters))
    
    if not args.kmerinfile:
        dicReads = generateKmerDic(refFiles,excludeFiles,kmerFwdRev,kmerSize,alphabet,minShEntropy,MaxratioAmb,mutedkmers,maxRepeatsFreq)
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

    if paired:
        if concatenated:
            yamlResume[os.path.basename(fastqFile)] = kmerRefFilterPairedProc(dicReads,kmerSize,fastqFile,"","",minMatchKeepSeq,args.fwdRevChoice,args.keepNotFiltered,namedPipe,kmerstats)
        else:
            yamlKey = "{} - {}".format(fastqFwdFile,fastqRevFile)
            yamlResume[yamlKey] = kmerRefFilterPairedProc(dicReads,kmerSize,"",fastqFwdFile,fastqRevFile,minMatchKeepSeq,args.fwdRevChoice,args.keepNotFiltered,namedPipe,kmerstats)
            if namedPipe:
                os.remove(fastqFwdFile)
                os.remove(fastqRevFile)
    else:
        yamlResume[os.path.basename(fastqFile)] = kmerRefFilterProc(dicReads,kmerSize,fastqFile,minMatchKeepSeq,args.keepNotFiltered,namedPipe,kmerstats)
        if namedPipe:
            os.remove(fastqFile)
    
    if args.yamlout:
        stderr_print(TitleFrame("Yaml output on stdout"))
        yaml.dump(yamlResume, sys.stdout, default_flow_style=False)

        
    stderr_print(TitleFrame("kmer statistics generation"))
    kmuse=len(kmerstats.keys())
    kmcountList=[kmvalues['count'] for km,kmvalues in kmerstats.items()]
    #print([v for k,v in kmerstats.items()])
    if args.kmerstats:
        dt=datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
        kmerstatsFileName=os.path.join(args.outputdir,"kmerstats_{}".format(dt))
        with open(kmerstatsFileName,'w') as kmfile:
            for km in dicReads:
                if km in kmerstats.keys():
                    kmfile.write("{kmer}\t{cnt}\t{entro}\t{repeats}\n".format(kmer=km,cnt=kmerstats[km]['count'],entro=get_ShannonEntropy(km),repeats=repeats_frequency(km)))
                elif args.kmerstats=='ALL':
                    kmfile.write("{kmer}\t{cnt}\t{entro}\t{repeats}\n".format(kmer=km,cnt=0,entro=get_ShannonEntropy(km),repeats=repeats_frequency(km)))
        
    stderr_print("Total kmer: {} -- kmers used nbr {} -- kmers nevers used {}\n".format(len(dicReads),kmuse,len(dicReads)-kmuse))
    stderr_print("stats: mean count {} -- standard deviation {} -- median count {}\n".format(np.mean(kmcountList),np.std(kmcountList),np.median(kmcountList)))
    hist=np.histogram(kmcountList)
    stderr_print("Distribution:\nFreq\t{}\nCounts\t{}\n".format("\t".join([str(v) for v in hist[0]]),"\t".join([str(round(v,1)) for v in hist[1]])))
    stderr_print("30 most used kmers:\n")
    tenSortedkmuse=sorted(kmerstats.items(), key=lambda x : x[1]['count'], reverse=True)[:30]
    stderr_print("\tkmer\tcount\tentropy\trepetitions")
    for km,v in tenSortedkmuse:
        rprslt=repeats_frequency(km)
        stderr_print("\t{kmer}\t{cnt}\t{entro}\t{repeats}".format(kmer=km,cnt=v['count'],entro=round(get_ShannonEntropy(km),2),repeats=(rprslt['letters'],round(rprslt['freq'],2))))

    if args.minkmercount>=0 and args.maxkmercount>0:
        stderr_print(TitleFrame("Clean filtered fastq files on kmer counts [{},{}]".format(args.minkmercount,args.maxkmercount)))
        nbrClean,total_reads,percent_retained,retainedkmerstats=clean_filteredFQ(yamlResume,kmerstats,kmerCountThrld,paired)
        stderr_print("\n{nbc} reads keeped afer cleaning ({pm}%)".format(nbc=nbrClean,pm=round(percent_retained,3)))
        stderr_print("\n30 most used kmers after kmer count cleaning:\n")
        tenSortedkmuse=sorted(retainedkmerstats.items(), key=lambda x : x[1]['count'], reverse=True)[:30]
        stderr_print("\tkmer\tcount\tentropy\trepetitions")
        for km,v in tenSortedkmuse:
            rprslt=repeats_frequency(km)
            stderr_print("\t{kmer}\t{cnt}\t{entro}\t{repeats}".format(kmer=km,cnt=v['count'],entro=round(get_ShannonEntropy(km),2),repeats=(rprslt['letters'],round(rprslt['freq'],2))))
                        
    elapTime = time.strftime("%Hh%Mm%Ss", time.gmtime(time.time()-start))
    
    stderr_print("\nFiltering terminated in {}".format(elapTime))
    stderr_print("\n{}".format(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]")))

if    __name__ == '__main__':
    kmerRefFilter(sys.argv[1:])
