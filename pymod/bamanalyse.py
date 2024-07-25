#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
import numpy as np
from collections import Counter

def read_bam(bamfilename):
    """
    This function extract informations from input bam file
    Input: bam file path
    Return a dictionary:
    'readsNb': reads number in bam
    'ratio_readsNoInDel': ratio of reads contains no InDel
    'readsNoInDelNbr': number of reads contains no InDel
    'mismatchs_list': list with all numbers of mismatches in the alignment
    'mismatch_mean': mean of mismatches in the alignment
    """
    inbam=pysam.AlignmentFile(bamfilename,'rb')
    #MD_list=[]

    #XM tag from bowtie2 manual:
    #XM:i:<N> The number of mismatches in the alignment. Only present if SAM record is for an aligned read.

    #MAPQ: MAPping Quality. It equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer. A value 255 indicates that the mapping quality is not available.
    XM_list=[]
    XM_Pct=[]
    MQ_list=[]
    readsNoInDelNbr=0
    ratio_readsNoInDel=0
    readsNb=0
    read_list=[r for r in inbam]
    if len(read_list)>0:
        for read in read_list:
            tags=tagsToDict(read.tags)
            #print(read.reference_start,tags['XM'],tags['MD'], read.query_name,read.cigarstring,containsInDel(read))
            readsNb+=1
            if not containsInDel(read):
                #MD_list.append(tags['MD'])
                XM_list.append(tags['XM'])
                XM_Pct.append(float(tags['XM'])/float(read.reference_length))
                MQ_list.append(read.mapping_quality)
        readsNoInDelNbr=len(XM_list)
        ratio_readsNoInDel=float(readsNoInDelNbr)/float(readsNb)
        #print(readsNb,readsNoInDelNbr)
        mutatedfraction={m:float(Counter(XM_list)[m])/float(readsNoInDelNbr) for m in sorted(Counter(XM_list).keys())}
    inbam.close()
    return {'readsNb':readsNb,'ratio_readsNoInDel':ratio_readsNoInDel,'readsNoInDelNbr':readsNoInDelNbr,'mismatchs_list':XM_list,'MAPQ_list':MQ_list,'mismatch_mean':float(np.mean(XM_list))}

def filter_bam(bamfilename,outbamname,mutThrld,filterInDel=True):
    """
    This function filter a bam file on the XM value with a mutation threshold
    Input:
    bamfilename : path to the input bam file
    outbamname : path to the filtered bam file
    mutThrld : mutation threshold, keep read have XM value <= mutThrld
    filterInDel : filter read with insertion/deletions

    Output:
    create the filtered bam file and the "..._OTHER.bam" file
    contains reads not kept.
    """
    inbam=pysam.AlignmentFile(bamfilename,'rb')
    outbamkept=pysam.AlignmentFile(outbamname,'wb', template=inbam)

    bamothername="{}_OTHER.bam".format(outbamname[:outbamname.rfind('.')])
    outbamother=pysam.AlignmentFile(bamothername,'wb', template=inbam)
    
    for read in inbam:
        tags=tagsToDict(read.tags)
        if tags['XM']<=mutThrld and (not containsInDel(read) and filterInDel):
            outbamkept.write(read)
        else:
            outbamother.write(read)
            
    inbam.close()
    outbamkept.close()
    outbamother.close()
    
def tagsToDict(tagslist):
    outdic={}
    for t in tagslist:
        outdic[t[0]]=t[1]
    return outdic

def containsInDel(read):
    return("I" in read.cigarstring or "D" in read.cigarstring)

def analyse_MD(MD_list):
    mutspaces_list=[]
    for md in MD_list:
        md_sub=re.sub(r'[A-Z\^+]','X',md)
        size=len(md_sub.split('X'))
        mutspaces=[int(v) for v in md_sub.split('X') if v]
        if not bool(re.match('^[A-Z\^+]',md)):
            mutspaces=mutspaces[1:]
        if not bool(re.match('[A-Z\^+]$',md)):
            mutspaces=mutspaces[:-1]
            
        if len(mutspaces)>1:
            a=float(sum(mutspaces))/float((size-2))
            b=(a**2)
            dist=np.mean([v**2 for v in mutspaces])/b
            print(md,mutspaces,sum(mutspaces),(size-1),dist)
            mutspaces_list.append(mutspaces)
