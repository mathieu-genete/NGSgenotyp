#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
from Bio import SeqIO
import sys
import os
import time
import datetime
import argparse

def split(sample_infile,outfld,outcompress=False):
    sample=sample_infile[0]
    infile=sample_infile[1]
    
    logout=""
    bnameInfile=os.path.basename(infile)

    logout+="== Split fastq files for sample{} ==\n".format(sample)

    if outcompress:
        filesuffix=".gz"
        logout+="Compress out fastq activated\n"
    else:
        filesuffix=".fastq"

    fwdsplited=os.path.join(outfld,"{}_R1_splited{}".format(bnameInfile[:bnameInfile.rfind('.')],filesuffix))
    revsplited=os.path.join(outfld,"{}_R2_splited{}".format(bnameInfile[:bnameInfile.rfind('.')],filesuffix))
    if not os.path.exists(fwdsplited) or not os.path.exists(revsplited):
        logout+="File to split: {}\n".format(infile)
        logout+="Output splited files:\n\t{}\n\t{}\n".format(fwdsplited,revsplited)
        dt=0
        nbr=0
        start=time.time()
        if outcompress:
            outR1=gzip.open(fwdsplited,'wb')
            outR2=gzip.open(revsplited,'wb')
        else:
            outR1=open(fwdsplited,'wb')
            outR2=open(revsplited,'wb')
        ingz=open_fastqFile(infile)
        for i,rec in enumerate(ingz):
            line=rec.strip()
            if i%4==0:
                nbr+=1
            if i%4==1 or i%4==3:
                recR2="{}\n".format(line[int(len(line)/2):])
                recR1="{}\n".format(line[:int(len(line)-len(line)/2)])
                outR1.write(recR1.encode())
                outR2.write(recR2.encode())
            else:
                recLine="{}\n".format(line)
                outR1.write(recLine.encode())
                outR2.write(recLine.encode())
        ingz.close()                
        outR1.close()
        outR2.close()
        logout+="Split Ended in {} sec\n".format(round(time.time()-start,1))
    else:
        logout+="File already splited: {}\n".format(infile)
        
    return {'sample':sample,'fqFWD':fwdsplited,'fqREV':revsplited,'log':logout}

def open_fastqFile(fastqFile):
    if is_gzip(fastqFile):
        fastq = gzip.open(fastqFile, 'rt', encoding='utf-8')
    else:
        fastq = open(fastqFile,"r")
    return fastq

def find_Magic_Bytes(filename,magic_bytes):
    with open(filename,'r',encoding="ISO-8859-1") as infile:
        file_start = infile.read(len(magic_bytes))
    if file_start.startswith(magic_bytes):
        return True		
    return False

def is_gzip(filename):
    magic_bytes = "\x1f\x8b\x08"
    return find_Magic_Bytes(filename,magic_bytes)
