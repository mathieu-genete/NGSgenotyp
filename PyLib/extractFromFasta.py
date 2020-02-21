#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os
import argparse
import sys
from Bio import SeqIO

__version__ = "1.0"

def run(ArgsVal):

    description = """ Extract sequences from multiple NGSgenotyp haploAsm contigs files """	

    parser = argparse.ArgumentParser(prog="extractFromFasta",description=description)

    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument("-f","--infastaspath", help="haploAsm results folder", required=True)
    parser.add_argument("-o","--outfolder", help="output folder for fasta files", required=True)
    parser.add_argument("-l","--constigslist", help="contigs list to extract in text format -- sample for a line: out_fasta_filename_(without exention),contig_Id_1,contig_Id_2,...,contig_Id_n", required=True)

    args = parser.parse_args(ArgsVal)

    yassContigsFiles=[]
    contigslist={}
    output_folder = args.outfolder

    stderr_print("{} - version {}".format(os.path.basename(__file__),__version__))
    stderr_print(" ".join(sys.argv))

    if not os.path.exists(args.constigslist):
        stderr_print("ERROR: {} NOT FOUND...".format(args.constigslist))
        sys.exit(2)

    fcontigslist=open(args.constigslist,'r')
    for line in fcontigslist:
        tmp=line.strip().split(',')
        fasta_fname = tmp[0]
        if fasta_fname not in contigslist.keys():
            contigslist[fasta_fname]=tmp[1:]
    fcontigslist.close()

    for root,dirs,files in os.walk(args.infastaspath):
        for file in files:
            if file.endswith('yass_oriented_contigs.fasta'):
                yassContigsFiles.append(os.path.join(root,file))

    if len(yassContigsFiles)==0:
        stderr_print("No contigs files found. Verify your folder (included subfolders) contains files ended with 'yass_oriented_contigs.fasta'")
        sys.exit(0)

    stderr_print("{} 'yass_oriented_contigs.fasta' files found in folder and subfolders\n".format(len(yassContigsFiles)))
        
    outFasta={}
    for fastaName in yassContigsFiles:
        fasta = SeqIO.parse(fastaName,'fasta')
        for rec in fasta:
            for OutFastaFname,contigs in contigslist.items():
                if OutFastaFname not in outFasta.keys():
                    outFasta[OutFastaFname]=""
                if str(rec.id) in contigs or len(contigs)==0:
                    stderr_print("- Extract [{}] FROM {} - output: {}\n".format(rec.id,fastaName,OutFastaFname))
                    outFasta[OutFastaFname] += ">{}\n{}\n".format(rec.id,rec.seq)

    dirFile_exist(output_folder,True)

    for outfastaFile,seqs in outFasta.items():
        outname = os.path.join(output_folder,"{}.fasta".format(outfastaFile))
        outFaFile = open(outname,'w')
        outFaFile.write(seqs)
        outFaFile.close()

def stderr_print(txt):
	sys.stderr.write(str(txt)+"\n")

#Test if a file or directory exist if creatdir set to 0: exit
# if creatdir set to 1: create new directory
def dirFile_exist(path,creatdir = False):
    if not os.path.exists(path):
        if not creatdir:
            stderr_print("ERROR: {} not found".format(path))
            sys.exit(1)
        else:
            dirmsg = "-- create {} directory".format(path)
            stderr_print(dirmsg)
            os.makedirs(path)
