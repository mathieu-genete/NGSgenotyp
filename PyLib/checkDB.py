#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import encodings
import subprocess
from Bio import SeqIO
from Bio import Seq
import codecs
import uuid

args=None

def run(ArgsVal):
    global args

    description = """ NGSgenotyp Check Fasta Database """
    parser = argparse.ArgumentParser(prog="checkDB",description=description)
    parser.add_argument("-f","--infasta", help="input fasta database file", required = True)
    parser.add_argument("-o","--outputfasta", help="output fasta normalized file")
    args = parser.parse_args(ArgsVal)
    file_type(args.infasta)
    outencoded=convertToUnix(args.infasta)
    encodedfasta="".join(outencoded)
    fastadict=parse_fasta(encodedfasta)
    analyse_fasta(fastadict)
    if args.outputfasta:
        save_encoded_file(outencoded,args.outputfasta,encod='utf-8')

def analyse_fasta(fastadict):
    NotAllowedChar=[',','.','/',' ',';','#','@']
    DNAalphabet=Seq.IUPAC.IUPACData.ambiguous_dna_values.keys()
    AllowedParameters=['Paralog','specie','grpRef','HaploId','gprId']
    #?:{'title':'?','err':[]}
    errorsout={1:{'title':'unsupported characters in sequence id:','err':[]},2:{'title':'unsupported characters in sequence string:','err':[]},3:{'title':'sequence id format','err':[]},4:{'title':'Gene ID duplicates','err':[]}}
    IdList=[]
    for idnbr,value in fastadict.items():
        #===== search Id errors=====
        idCharERROR=False
        for c in NotAllowedChar:
            if c in value['description']:
                cpos=[str(i) for i,ltr in enumerate(value['description']) if ltr==c]
                positionStr=list("."*len(value['description']))
                for p in cpos:
                    positionStr[int(p)]="^"
                lineheader="\n\tid number: {idnb}\t".format(idnb=idnbr)
                errorsout[1]['err'].append("{head}{seqid} -- character: '{ch}' -- position(s): {pos}\n\t{lnh} {posstr}".format(head=lineheader,seqid=value['description'],ch=c,pos=",".join(cpos),posstr="".join(positionStr),lnh=" "*len(lineheader)))
                idCharERROR=True
        
        parseID=value['description'].split("|")
        geneId=parseID[0]
        if geneId not in IdList:
            IdList.append(parseID[0])
        else:
            errorsout[4]['err'].append("\n\tDuplicate found for id: '{}' -- id number: {}".format(geneId,idnbr))

        #===== search parametters errors=====
        paramsSTR=parseID[1:]
        for p in paramsSTR:
            tmp=p.split('=')
            if len(tmp)==1 and tmp[0]=="":
                errorsout[3]['err'].append("\n\tid number: '{}' -- id: {} -- EMPTY parameter: '||'".format(idnbr,geneId))
            elif len(tmp)==1 or len(tmp)>2:
                errorsout[3]['err'].append("\n\tid number: '{}' -- id: {} -- param error: {}".format(idnbr,geneId,p))
            else:
                if tmp[0] not in AllowedParameters:
                    errorsout[3]['err'].append("\n\tid number: '{}' -- id: {} -- paramameter not allowed: {}".format(idnbr,geneId,tmp[0]))

        #===== search sequences errors=====
        sequence=value['seq']
        SequenceNotAllowedAlphabet=[b for b in set(list(sequence.upper())) if b not in DNAalphabet]

        for base in SequenceNotAllowedAlphabet:
            bpos=[str(i) for i,ltr in enumerate(sequence) if ltr.upper()==base]
            errorsout[2]['err'].append("\n\tid number: '{}' -- id: {} -- base '{}' not in allowed alphabet -- position: {} bp".format(idnbr,geneId,base,",".join(bpos)))

    for err in errorsout.values():
        if len(err['err'])>0:
            print("\n{}".format(err['title']))
            for chaine in err['err']:
                print("{}".format(chaine))

def parse_fasta(outencoded):
    fastadict={}
    tmpname=os.path.join("/tmp",str(uuid.uuid4()))
    with open(tmpname,"w") as tmpfasta:
        for line in outencoded:
            tmpfasta.write(line)
    fasta=SeqIO.parse(tmpname,'fasta')
    seqNbr=0
    for rec in fasta:
        seqNbr+=1
        fastadict[seqNbr]={'id':str(rec.id),'description':str(rec.description),'seq':str(rec.seq)}
    os.remove(tmpname)
    return fastadict

def file_type(infile):
    AllowedTypes=set(encodings.aliases.aliases.values())
    cmd="file -b {}".format(infile)
    job = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,sdterr = job.communicate()
    ftype=None
    for t in AllowedTypes:
        if t.upper() in stdout.strip().upper():
            ftype=t
    crlf=("CRLF" in stdout.strip().upper())
    print(stdout.strip())
    return ftype,crlf

def convertToUnix(infile):
    MsDOS_endline='\r\n'
    convertedFasta=[]
    with codecs.open(infile,'r') as file:
        for line in file.readlines():
            IS_MsDOS_EL=line[-2:]==MsDOS_endline
            if IS_MsDOS_EL:
                convertedFasta.append("{}\n".format(line[:-2]))
            else:
                convertedFasta.append("{}\n".format(line.strip()))
    return convertedFasta

def save_encoded_file(outencoded,outfile,encod='utf-8'):
    print("output file {} with encoding {}".format(outfile,encod))
    with codecs.open(outfile,'w',encoding=encod) as wfile:
        for line in outencoded:
            wfile.writelines(line.decode(encod))

if __name__ == '__main__':
    ArgsVal = sys.argv[1:]
    run(ArgsVal)