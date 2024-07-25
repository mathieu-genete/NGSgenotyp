#!/usr/bin/env python3
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

__version__ = "0.5"
AppName = "checkDB"

def run(infasta):
    outcheck=""
    outcheck+="=== {} v{} ===\n".format(AppName,__version__)
    outcheck+="input fasta: {}\n".format(infasta)
    fileOut,ftype,crlf=file_type(infasta)
    outencoded=convertToUnix(infasta)
    encodedfasta="".join(outencoded)
    fastadict=parse_fasta(encodedfasta)
    outcheck+="=> START CHECK\n"
    fastaOK,logout=analyse_fasta(fastadict)
    outcheck+=logout
    return fastaOK,outcheck
    
def analyse_fasta(fastadict):
    logout=""
    fastaOK=True
    NotAllowedChar=[',','.','/',' ',';','#','@',':']
    #NotAllowedChar=['.','/',' ',';','#','@']
    DNAalphabet=['A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'T', 'W', 'V', 'Y', 'X']
    AllowedParameters=['Paralog','specie','grpRef','HaploId','gprId','allelePart','grpRefPart']
    comma_exception=['allelePart','grpRefPart']
    deprecatedParameters=['HaploId','gprId']
    warnerr={'E':'ERROR - ','W':'WARNING - '}

    #?:{'title':'?','err':[]}
    errorsout={1:{'title':'unsupported characters in sequence id:','err':[]},2:{'title':'unsupported characters in sequence string:','err':[]},3:{'title':'sequence id format','err':[]},4:{'title':'Gene ID duplicates','err':[]}}
    IdList=[]

    logout+="Sequences number: {}\n".format(len(fastadict))

    for idnbr,value in fastadict.items():
        
        parseID=value['description'].split("|")
        geneId=parseID[0]
        
        #===== search Id errors=====
        idCharERROR=False
        for c in NotAllowedChar:
            if c in value['description']:
                param_with_comma=[v.split('=')[0] for v in parseID if ',' in v and v.split('=')[0] not in comma_exception]
                if len(param_with_comma)==0:
                    continue
                cpos=[str(i) for i,ltr in enumerate(value['description']) if ltr==c]
                positionStr=list("."*len(value['description']))
                for p in cpos:
                    positionStr[int(p)]="^"
                lineheader="\n\t{we}id number: {idnb}\t".format(we=warnerr['E'],idnb=idnbr)
                errorsout[1]['err'].append("{head}{seqid} -- character: '{ch}' -- position(s): {pos}\n\t{lnh} {posstr}".format(head=lineheader,seqid=value['description'],ch=c,pos=",".join(cpos),posstr="".join(positionStr),lnh=" "*len(lineheader)))
                idCharERROR=True
                fastaOK=False
        
        if geneId not in IdList:
            IdList.append(parseID[0])
        else:
            errorsout[4]['err'].append("\t{}Duplicate found for id: '{}' -- id number: {}".format(warnerr['E'],geneId,idnbr))
            fastaOK=False

        #===== search parametters errors=====
        paramsSTR=parseID[1:]
        for p in paramsSTR:
            tmp=p.split('=')
            if len(tmp)==1 and tmp[0]=="":
                errorsout[3]['err'].append("\t{}id number: '{}' -- id: {} -- EMPTY parameter: '||'".format(warnerr['E'],idnbr,geneId))
                fastaOK=False
            elif len(tmp)==1 or len(tmp)>2:
                errorsout[3]['err'].append("\t{}id number: '{}' -- id: {} -- param error: {}".format(warnerr['E'],idnbr,geneId,p))
                fastaOK=False
            else:
                if tmp[0] not in AllowedParameters:
                    errorsout[3]['err'].append("\t{}id number: '{}' -- id: {} -- paramameter not allowed: {}".format(warnerr['E'],idnbr,geneId,tmp[0]))
                    fastaOK=False
                if tmp[0] in deprecatedParameters:
                    errorsout[3]['err'].append("\t{}id number: '{}' -- id: {} -- deprecated parameter: {}".format(warnerr['W'],idnbr,geneId,tmp[0]))
                if tmp[0]=='grpRef':
                    grpRefVal=tmp[1]
                    grpRefSplit=grpRefVal.split("-")
                    if len(grpRefSplit)==0 or len(grpRefSplit)>2:
                        errorsout[3]['err'].append("\t{}id number: '{}' -- id: {} -- grpRef value: {}".format(warnerr['E'],idnbr,geneId,grpRefVal))
                        fastaOK=False
                    elif len(grpRefSplit)==1 and len(grpRefVal)==5:
                        errorsout[3]['err'].append("\t{}id number: '{}' -- id: {} -- grpRef value deprecated use 'Hg-h': {}".format(warnerr['W'],idnbr,geneId,grpRefVal))

        #===== search sequences errors=====
        sequence=value['seq']
        SequenceNotAllowedAlphabet=[b for b in set(list(sequence.upper())) if b not in DNAalphabet]

        for base in SequenceNotAllowedAlphabet:
            bpos=[str(i) for i,ltr in enumerate(sequence) if ltr.upper()==base]
            errorsout[2]['err'].append("\t{}id number: '{}' -- id: {} -- base '{}' not in allowed alphabet -- position: {} bp".format(warnerr['E'],idnbr,geneId,base,",".join(bpos)))
            fastaOK=False

    for err in errorsout.values():
        if len(err['err'])>0:
            logout+="{}\n".format(err['title'])
            for chaine in err['err']:
                logout+="{}\n".format(chaine)

    return fastaOK,logout

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
    stdout=stdout.decode("utf-8")
    ftype=None
    for t in AllowedTypes:
        if t.upper() in stdout.strip().upper():
            ftype=t
    crlf=("CRLF" in stdout.strip().upper())
    return stdout.strip(),ftype,crlf

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
