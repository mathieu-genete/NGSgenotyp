#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Mathieu Genete

@NGS-feature
@NGS-description: execute genotyp pipeline from raw NGS reads data
"""

#=== imports libraries ===
import os
import sys
import argparse
import subprocess
import collections
import xlrd
import xlwt
import pysam
from Bio import SeqIO
from Bio import SeqUtils
from mpl_toolkits.mplot3d import Axes3D
import random
import math
import numpy as np
from scipy.stats import gaussian_kde
from itertools import cycle
import hashlib
import shutil
from datetime import datetime
import time
import threading
import ctypes
from urllib.parse import urlparse
from urllib.request import urlopen
from operator import itemgetter

#Multimodal test
from unidip import UniDip
import dip

#Multiprocessing
from multiprocessing import Pool
from functools import partial

#home libraries
import utils
import checkDB
import splitfastq
import bamanalyse
import genotyp_graphs
import calculate_genotypes as cg
import genotyp_depth as gd
import load_dump as ld
import genotyp_output as go

utils.set_paramFileName('NGSgenotyp2_genotyp')

__version__="2.0.1"
AppName = "NGSgenotyp genotyp"

config = None
fastqParams = None
args = None
samtools_path=''
bowtie2_path = ''
kmerrefilter_path=''

def main(ArgsVal,main_configs):

    global config
    global args
    global samtools_path
    global bowtie2_path
    global kmerrefilter_path
    
    #TO REMOVE
    #print("genotyp feature: ",__name__)
    #print(main_configs)

    startTime = time.time()
    
    #ArgsParser
    description = """ NGSgenotyp v2 -- genotyping pipeline """
    
    parser = argparse.ArgumentParser(prog="genotyp",description=description)
    
    parser.add_argument("-V", "--verbose", help="full verbose", action="store_const", const="-V", default="")
    parser.add_argument("-v", "--tinyverbose", help="verbose",action="store_const", const="-v", default="")
    parser.add_argument("-t", "--test", help=argparse.SUPPRESS, action="store_const", const="-t", default="")
    parser.add_argument("-f", "--force", action="store_const", const="-f", default="")
    parser.add_argument("-c", "--checks", help=argparse.SUPPRESS, action="store_const", const="-c", default="")
    parser.add_argument("-k", "--kmerfilter", help="filtering fastq raw data with kmers dictionnary generated from references sequences. Should be use if your input fastq are not yet filtered (significatively reduces compute time)", action="store_const", const="-k", default="")
    parser.add_argument("-ks","--kmerSize", help="kmer size for kmer filtering - default = 20", type=int, default=20)
    parser.add_argument("-m","--mismatchthrld", help="mismatch threshold on aligned reads (remove reads whose mitmatch value is greater than threshold)", type=int,default=0)
    parser.add_argument("-pdf", "--pdfreports",help="generate PDF reports (can take long time)", action="store_const", const="-pdf", default="")
	
    parser.add_argument("-s","--statsonly", action="store_const", const=True, default=False, help="do stats only")
    parser.add_argument("-sh","--sharedreads", action="store_const", const=True, default=False, help="generates shared reads file")
    parser.add_argument("-sm","--statsmismatchthrld", action="store_const", const=True, default=False, help="do stats with filtered bam with mismatchthrld value")
    parser.add_argument("-kid","--keepindel", action="store_const", const=True, default=False, help="force keep reads with insertion/deletion during bam filtering step")
    parser.add_argument("-T","--MaxParallelsJobs", help="max number of parallels jobs to run - default= see config file", type=int)
    parser.add_argument("-e","--ErroRateThrld", help="Force error rate threshold - default= see config file", type=float)
    parser.add_argument("-M","--readsMappedThrld", help="Take into account only alignments with 'reads mapped' greater then threshold (default = 10)", type=int,default=10)
    parser.add_argument("-A","--alignmentMode", help="bowtie2 reads alignments set to local or end-to-end (default = end-to-end)", type=str,default="end-to-end")
    parser.add_argument("-S","--alignmentSensitivity", help="bowtie2 reads alignments sensitivity set to very-fast, fast, sensitive or very-sensitive (default = sensitive)", type=str,default="sensitive")
    parser.add_argument("-o","--outfolder", help="destination folder (create it if not exist)", required=True)
    parser.add_argument("-i","--readsinfo", help="Configuration file with reads informations if reads ares paired add [format=paired] parameter", required=True)
    parser.add_argument("-d","--refdatabase", help="reference database in fasta format (see documentation)", required=True)
    parser.add_argument("-x","--refexclude", help="simple text file contains for each line, references names to exclude for current analysis")

    parser.add_argument("-on","--outnbsheetperxls", help="maximum sheet number by output xls files (default = 15)", type=int,default=15)

    args = parser.parse_args(ArgsVal)
    
    #Check arguments
    alignmentMode=['local','end-to-end']
    alignmentSensitivity=['very-fast','fast','sensitive','very-sensitive']

    if args.alignmentMode not in alignmentMode:
        exitProg(msg='Bad Alignment mode...')

    if args.alignmentSensitivity not in alignmentSensitivity:
        exitProg(msg='Bad Alignment sensitivity...')

    #Check if files exists
    files_exits([args.readsinfo,args.refdatabase,args.refexclude])
    
    #genotyp and tools configurations
    config = ld.yaml_loadDic(os.path.join(main_configs['configs_folder'],'genotyp_config.yaml'))
    samtools_cfg=main_configs['tools_avail']['samtools']
    bowtie2_cfg = main_configs['tools_avail']['bowtie2']
    kmerrefilter_cfg = main_configs['tools_avail']['kmerRefFilter']

    if args.ErroRateThrld:
        config['genotyp_def_ErrorRate']=args.ErroRateThrld

    if args.outnbsheetperxls>=1:
        config['MaxSheetNbr']=args.outnbsheetperxls
    else:
        config['MaxSheetNbr']=1

    config['reads_mapped_THRLD']=args.readsMappedThrld
    
    #update config variable
    config['outfolder']=args.outfolder
    config['log_folder']=os.path.join(config['outfolder'],'logs')
    
    #Create output and logs folder
    for folder in [config['outfolder'],config['log_folder']]:
        if not os.path.exists(folder):
            os.mkdir(folder)
    
    #genotyp log informations
    cmdline = " ".join(sys.argv)
    currentDateTime = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    config['file_time']=datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
    config['log_file'] = "{appn}_Log_{dt}.txt".format(appn=AppName,dt=config['file_time'])
    startText = "{appn} {appv} -- {cdate}".format(appn=AppName,appv=__version__,cdate=currentDateTime)
    starDateline = "*"*len(startText)
    addTextToLogFile("{}\n{}\n{}".format(starDateline,startText,starDateline),False)
    addTextToLogFile("user command: {}".format(cmdline))

    #Configure tools variable (local use or system)
    if samtools_cfg['uselocal']:
        samtools_path=samtools_cfg['folder']
        addTextToLogFile("USE local samtools : {}".format(samtools_path))
    else:
        samtools_path=''
        
    if bowtie2_cfg['uselocal']:
        bowtie2_path = bowtie2_cfg['folder']
        addTextToLogFile("USE local bowtie2 : {}".format(bowtie2_path))
    else:
        bowtie2_path=''

    if kmerrefilter_cfg['uselocal']:
        kmerrefilter_path = kmerrefilter_cfg['folder']
        addTextToLogFile("USE local KmerRefFilter : {}".format(kmerrefilter_path))
    else:
        kmerrefilter_path = ''

    #Set MaxParallelsJobs value
    if args.MaxParallelsJobs:
        config['MaxParallelsJobs'] = args.MaxParallelsJobs
    stdout_print("Max parallels job set to {}".format(config['MaxParallelsJobs']))
    utils.setMaxParallelJobsToConfig(config['MaxParallelsJobs'])

    #create output folders
    for f,v in config['OutFolders'].items():
        v=os.path.join(config['outfolder'],v)
        dir_create(v)
        config['OutFolders'][f]=v

    #path configuration for output files
    config['splited_reads']=os.path.join(config['OutFolders']['out_splited'],config['splited_reads'])
    refbname=os.path.basename(args.refdatabase)
    config['kmer_Ref_PickleFile']=os.path.join(config['OutFolders']['reference_DB_fld'],"KRFDB_{}_k{}.pk".format(refbname[:refbname.rfind('.')],args.kmerSize))
    config['readsconfig']=os.path.join(config['OutFolders']['configs_fld'],config['readsconfig'])

    #***********************************************************
    #*=== MAIN WORKFLOW ===*
    #***********************

    #==============================
    #check files (reads, fasta,...)
    
    #Check DB Fasta
    fastaOK,outcheck=checkDB.run(args.refdatabase)
    for line in outcheck.strip().split('\n'):
        addTextToLogFile(line)

    if not fastaOK:
        exitProg("ERROR in database file '{}' - see log file for more informations".format(args.refdatabase))
    else:
        stdout_print("Check database '{}' OK".format(os.path.basename(args.refdatabase)),printInLog=True)

    #Load and check reads file
    stdout_print("Parse reads file '{}'".format(args.readsinfo),printInLog=True)
    if check_FileExtensions(args.readsinfo,['yaml','yml']):
        readsConfig = ld.yaml_loadDic(args.readsinfo)
    else:
        readsConfig = load_readsList(args.readsinfo)
        
    #Load excluded references names
    config['refexclude']=[]
    if args.refexclude:
        if os.path.exists(args.refexclude):
            tmp=[]
            with open(args.refexclude,'r') as exfile:
                for exref in exfile:
                    tmp.append(exref.strip())
            config['refexclude']=tmp
        if len(config['refexclude'])>0:
            addTextToLogFile("References to exluce : {}".format(",".join(config['refexclude'])))

    #Database file parsing
    pk_refdatabase=parse_referenceDB(args.refdatabase)
    #create index for references
    yaml_index=index_references(pk_refdatabase)

    #Split reads for concatenated reads
    # if kmerfilter not used, split raw reads (if taged in read txt file)
    
    #filter using kmerRefFilter if set
    if args.kmerfilter:
        readsConfig = run_kmerRefFilter(args.kmerSize,readsConfig)
        #Split filtered fastq if split is set
        readsConfig = split_reads(readsConfig,compress=False)
    else:
        readsConfig = split_reads(readsConfig)

    #save readsConfig
    ld.yaml_dumpDic(config['readsconfig'],readsConfig)

    #Alignment using Bowtie2 (unpaired alignments)
    bwt_results,reads_nbr=run_bowtie2(readsConfig,yaml_index,config['refexclude'],args.alignmentMode,args.alignmentSensitivity,pairedAln=False)
    
    #Metrics evaluation and statistics
    outstats,depthfiles = alignments_stats(bwt_results,pk_refdatabase)

    stdout_print("Combine Alleles...",printInLog=True)
    combined_stats=cg.combine_alleles_stats(outstats,pk_refdatabase)

    stdout_print("Determine genotypes...",printInLog=True)
    cg.determine_genotypes(combined_stats,pk_refdatabase,config['lbda'],config['mu'],config['sigma'],config['genotyp_alleleProb_THRLD'])

    stdout_print("Normalize depth...",printInLog=True)
    gd.normalize_depth(combined_stats,pk_refdatabase)

    go.write_xls_output(config['outfolder'],config['OutFolders']['main_results_dir'],combined_stats,reads_nbr,pk_refdatabase,args.mismatchthrld,config)

    go.write_genotypTXT(config['outfolder'],config['OutFolders']['main_results_dir'],combined_stats,pk_refdatabase)
    print(config['OutFolders']['main_results_dir'])
    stdout_print("Save stats...",printInLog=True)
    samtools_stats=save_stats(combined_stats)
    assembly_config_file(pk_refdatabase)
    
    #***************************
    #*=== END MAIN WORKFLOW ===*
    #********************************************************************

    #delete utils param file
    utils.delete_paramFileName()

    totalTime = time.time() - startTime
    endMessage = "-- ENDED IN {} --".format(time.strftime('%Hh%Mm%Ss', time.gmtime(totalTime)))
    stdout_print(endMessage)
    addTextToLogFile(endMessage)
    #END =================================================================


def assembly_config_file(pk_refdatabase):
    samtools_stats_yaml = os.path.join(config['OutFolders']['main_results_dir'],config['samtools_stats'])
    readsconfig="/".join(config['readsconfig'].split("/")[1:])
    asmconfig={'readsconfig':readsconfig,'samtools_stats':"/".join(samtools_stats_yaml.split("/")[1:]),'ref':"/".join(pk_refdatabase.split("/")[1:])}
    asm_configfile=os.path.join(config['OutFolders']['configs_fld'],config['assembly_config'])
    ld.yaml_dumpDic(asm_configfile,asmconfig)
    
def alignments_stats(bwt_results,pk_refdatabase):
    stdout_print("Do Alignments stats...",printInLog=True)
    tps_stat=time.time()
    outstats={}
    refinfo=ld.pickle_loadDic(pk_refdatabase)
    depthfiles={}
    sample_nbr=len(bwt_results.keys())
    c=0
    tps1=0
    ellapsed_times=[]
    for sample, refdatas in bwt_results.items():
        outstats[sample]={}
        depthfiles[sample]={}
        c+=1
        if tps1>0:
            ellapsed_times.append(time.time()-tps1)
            ETE=float(np.mean(ellapsed_times))*(sample_nbr-(c-1))
        else:
            ETE=0
        stdout_print("\t{}/{} => do stats for {} -- ETE {} s".format(c,sample_nbr,sample,round(ETE,2)),printInLog=True)
        tps1 = time.time()
        for refID, idval in refdatas.items():
            stats_file=idval['stats']
            reduced_stats_file=idval['reduced_stats']
            reflength=refinfo[refID]['len_seq']
            stats_results=parse_stats_file(stats_file)
            reduced_stats_results=parse_stats_file(reduced_stats_file)
            depthyaml=os.path.join(idval['aln_path'],'depthinfos_{}'.format(idval['aln_name']))
            reduced_depthyaml=os.path.join(idval['aln_path'],'reduced_depthinfos_{}'.format(idval['aln_name']))
            if stats_results['reads mapped']>config['reads_mapped_THRLD']:
                #print(idval['aln_name'])
                depthfiles[sample][refID]={'sorted':depthyaml,'reduced':reduced_depthyaml}
                outstats[sample][refID]={'stats':stats_results}
                outstats[sample][refID].update({'reduced_stats':reduced_stats_results})
                outstats[sample][refID].update(analyse_bam(idval['bam_sorted']))
                outstats[sample][refID].update({'ref_length':reflength})
                outstats[sample][refID].update({'coverage':extract_coverage(idval['bam_sorted'],depthyaml,reflength)})
                outstats[sample][refID].update({'reduced_coverage':extract_coverage(idval['reduced_sorted'],reduced_depthyaml,reflength)})
    stdout_print("END Alignments stats in {} s".format(round(time.time()-tps_stat,2)),printInLog=True)
    return outstats,depthfiles

def save_stats(outstats):
    samtools_stats=os.path.join(config['OutFolders']['main_results_dir'],config['samtools_stats'])
    ld.yaml_dumpDic(samtools_stats,outstats)
    return samtools_stats
    
def analyse_bam(bamfile):
    #{'ratio_readsNoInDel':ratio_readsNoInDel,'readsNoInDelNbr':readsNoInDelNbr,'mismatchs_list':XM_list,'mismatch_mean':float(np.mean(XM_list))}

    mismatchthrld=args.mismatchthrld
    
    bamdatas=bamanalyse.read_bam(bamfile)
    XM_list=bamdatas['mismatchs_list']
    
    if len(XM_list)>0:
        n_bins=max(XM_list)-min(XM_list)
        if n_bins==0: n_bins=1
        #if 0 in XM_list: print(bamfile)
        distrib_XM_values=[[int(i) for i in v] for v in np.histogram(XM_list,bins=n_bins)]
    else:
        distrib_XM_values=[]
    Nb_below_Thrld=sum([1 for v in XM_list if v<=mismatchthrld])
    ratio_below_Thrld = 0
    if len(XM_list)>0:
        ratio_below_Thrld = float(Nb_below_Thrld)/float(len(XM_list))
        
    #MAPQ: MAPping Quality. It equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer. A value 255 indicates that the mapping quality is not available.
    MQ_list=[v for v in bamdatas['MAPQ_list'] if v!=255]
    if len(MQ_list)>0:
        mean_MQ = float(np.mean(MQ_list))
    else:
        mean_MQ = 0
    #prob all mapping positions are wrong = (1/10)^(PHRED_SCORE/10)
    MQ_prob=0.1**(mean_MQ/10)
    
    return {'readsNb':bamdatas['readsNb'],'ratio_readsNoInDel':bamdatas['ratio_readsNoInDel'],'readsNoInDelNbr':bamdatas['readsNoInDelNbr'],'mismatch_mean':bamdatas['mismatch_mean'],'distrib_XM_values':distrib_XM_values,'MQ_prob':MQ_prob}

#TO DO -- REMOVE gname
#DO Only for positiv alleles...
def test_mutltimodal(data):
    data = np.msort(data)
    test=dip.diptst(data)
    pvalue=test[1]
    return pvalue

def extract_coverage(bamfile,depthyaml,reflength):
    cmd = "{path}samtools depth {bam} | cut -f 2-".format(path=samtools_path,bam=bamfile)
    depthrstl = utils.run(args,"samtools depth", '', cmd, retOUT=True)
    depthstdout=depthrstl[0].decode('utf-8')
    
    depth_datas={}
    for depthline in depthstdout.split('\n'):
        depth_array=depthline.split('\t')
        if len(depth_array)==2:
            depth_datas[int(depth_array[0])]=int(depth_array[1])
            
    depth_mean=float(sum(depth_datas.values()))/float(reflength)
    depth_covered_mean=float(np.mean(list(depth_datas.values())))
    covered_pos=len(depth_datas)
    covered_pct=float(covered_pos)/float(reflength)
    min_depth=0
    max_depth=0
    if reflength==len(depth_datas):
        min_depth=min(depth_datas.values())
    if len(depth_datas.values())>0:
        max_depth=max(depth_datas.values())
    ld.yaml_dumpDic(depthyaml,{'ref_length':reflength,'depth_datas':depth_datas})
    
    return {'depth_mean':depth_mean,'depth_covered_mean':depth_covered_mean,'covered_pos':covered_pos,'covered_pct':covered_pct,'min_depth':min_depth,'max_depth':max_depth}
        
    
def parse_stats_file(stats_file):
    if os.path.exists(stats_file):
        filestats={}
        with open(stats_file) as stfile:
            for line in stfile:
                tmp=[v.strip() for v in line.split(':')]
                keyid=tmp[0]
                value=float(tmp[1])
                filestats[keyid]=value
            if filestats['bases mapped']>0:
                error_rate=filestats['mismatches']/filestats['bases mapped']
                return {'raw sequences':filestats['raw total sequences'],'reads mapped':filestats['reads mapped'],'bases mapped':filestats['bases mapped'],'mismatches':filestats['mismatches'],'average length':filestats['average length'],'error rate':error_rate}
    return {'reads mapped':0,'bases mapped':0,'mismatches':0,'average length':0,'error rate':0}
                
def run_bowtie2(readsConfig,yaml_index,refexclude,alignmentMode,alignmentSensitivity,pairedAln):
    addTextToLogFile("-- ALIGNMENTS WITH BOWTIE2 --")
    index_dic=ld.yaml_loadDic(yaml_index)
    aln_name_template="Aln_{refID}_VS_{sample}"
    align_outdic={}
    errBWTfile_list=[]
    jblst = utils.Jobslist("bowtie2 for all")

    if len(refexclude)>0:
        stdout_print("Exclude reference file set",printInLog=True)
    for xref in refexclude:
        stdout_print("\t exclude reference: {}".format(xref),printInLog=True)
        
    for sample,values in readsConfig.items():
        sample_bwt_path=os.path.join(config['OutFolders']['bwt_dirResult'],sample)
        dir_create(sample_bwt_path)
        
        #Paired or single alignments
        if pairedAln and len(values['reads'])>1 and len(values['reads'])%2==0:
            stdout_print("Do paired alignments",printInLog=True)
            forward_files=[s for i,s in enumerate(values['reads']) if i%2==0]
            reverse_files=[s for i,s in enumerate(values['reads']) if i%2!=0]
            samples="-1 {fwd} -2 {rev}".format(fwd=','.join(forward_files),rev=','.join(reverse_files))
            flags="-f 0x3"
        else:
            stdout_print("Do single alignments",printInLog=True)
            samples = "-U {}".format(','.join(values['reads']))
            flags="-F 0x4"
            
        align_outdic[sample]={}
        for refID,indexval in index_dic.items():
            if refID not in refexclude:
                aln_name=aln_name_template.format(sample=sample,refID=refID)
                aln_path=os.path.join(sample_bwt_path,aln_name)
                dir_create(aln_path)

                errBWTfile=os.path.join(aln_path, "stderr_{}".format(aln_name))
                errBWTfile_list.append(errBWTfile)
                outbam=os.path.join(aln_path,"{}.bam".format(aln_name))
                reducedbam=os.path.join(aln_path,"{}_reduced.bam".format(aln_name))
                reducedbam_index=os.path.join(aln_path,"{}_reduced.bam.bai".format(aln_name))
                outbam_sorted=os.path.join(aln_path,"{}_SORTED.bam".format(aln_name))
                outbam_index=os.path.join(aln_path,"{}_SORTED.bam.bai".format(aln_name))
                #mapped_unmapped_out = os.path.join(aln_path, "mapped_{}".format(aln_name))
                outstats=os.path.join(aln_path,"{}_samtools_stats".format(aln_name))
                reduced_outstats=os.path.join(aln_path,"{}_samtools_stats_reduced".format(aln_name))
                    
                bwt_cmd="{bpath}bowtie2 --phred33 --{alnmode} --{alnsens} -x {RefIndex} {samples} 2>{err_bowtie} | {path}samtools view -S -b {flags} - > {bamFile}".format(bpath=bowtie2_path,path=samtools_path,RefIndex=indexval['index'],samples=samples, err_bowtie=errBWTfile, bamFile=outbam,flags=flags,alnmode=alignmentMode,alnsens=alignmentSensitivity)
                align_outdic[sample][refID]={'aln_name':aln_name,'aln_path':aln_path,'aln_path':aln_path,'errBWTfile':errBWTfile,'bam':outbam,'orig_sorted':outbam_sorted,'reduced_sorted':reducedbam,'reduced_index':reducedbam_index,'bam_sorted':outbam_sorted,'bam_index':outbam_index,'orig_index':outbam_index,'bwt_cmd':bwt_cmd,'stats':outstats,'orig_stats':outstats,'reduced_stats':reduced_outstats,'align_ok':False}
                addTextToLogFile("bowtie2 command for {} : {}".format(aln_name,bwt_cmd))
                jblst.add_a_job(bwt_cmd,"bowtie2 unpaired {}".format(aln_name),outbam)

    align_nbr=sum(len(v) for v in align_outdic.values())

    stdout_print("Launch {} alignments".format(align_nbr),printInLog=True)

    #Launch all bowtie2 commands

    thread_list=[]
    #show progress bar if verbose option is set
    if args.tinyverbose:
        t1 = threading.Thread(target=progress_bar,args=['Alignments progress',len(errBWTfile_list),thread_list])
        t1.start()

    utils.trun(args, jblst,thrd=thread_list)

    if args.tinyverbose:
        wait_progressbar_finished(t1)

    #check alignments completion
    errnbr,align_outdic=check_alignments(align_outdic)
    
    stdout_print("Alignments results: {} completed / {} misaligned".format(align_nbr-errnbr,errnbr),printInLog=True)
    if errnbr>0:
        stdout_print("\tSee log file for alignments errors".format(align_nbr-errnbr,errnbr))

    stdout_print("\nSort, index and stats for bam files")

    list_stats=[]
    list_stats_reduced=[]
    list_index=[]
    list_index_reduced=[]
    list_sort=[]
    reads_nbr={}
    for sample,refAlnDatas in align_outdic.items():
        for refID, aligndatas in refAlnDatas.items():
            if align_outdic[sample][refID]['align_ok']:
                list_stats.append((align_outdic[sample][refID]['bam_sorted'],align_outdic[sample][refID]['stats']))
                list_stats_reduced.append((align_outdic[sample][refID]['reduced_sorted'],align_outdic[sample][refID]['reduced_stats']))
                list_index.append((align_outdic[sample][refID]['bam_sorted'],align_outdic[sample][refID]['bam_index']))
                list_index_reduced.append((align_outdic[sample][refID]['reduced_sorted'],align_outdic[sample][refID]['reduced_index']))
                list_sort.append((align_outdic[sample][refID]['bam'],align_outdic[sample][refID]['bam_sorted']))
                if sample not in reads_nbr.keys():
                    reads_nbr[sample]=get_reads_nbr(align_outdic[sample][refID]['errBWTfile'])
                
    launch_samtools_sort(list_sort)
    launch_samtools_index(list_index)

    launch_samtools_stats(list_stats)
    
    #Filtrer avant samtools stats => à revoir
    align_outdic=read_filter_bam(align_outdic,args.mismatchthrld)
    
    launch_samtools_index(list_index_reduced)
                
    launch_samtools_stats(list_stats_reduced)
    
    return align_outdic,reads_nbr

def get_reads_nbr(stderr_file):
    with open(stderr_file,'r') as sfile:
        line = sfile.readline()
        try:
            return int(line.split()[0])
        except:
            return 0
    
def read_filter_bam(align_outdic,mutThrld):
    stdout_print("Start reduce",printInLog=True)
    for sample, refdatas in align_outdic.items():
        for refID, idval in refdatas.items():
            bam_file=idval['bam_sorted']
            reduced_bam=idval['reduced_sorted']
            
            if not os.path.exists(reduced_bam):
                filtInDel=(not args.keepindel)
                bamanalyse.filter_bam(bam_file,reduced_bam,mutThrld,filterInDel=filtInDel)
            
            if args.statsmismatchthrld:
                align_outdic[sample][refID]['bam_sorted']=reduced_bam
                align_outdic[sample][refID]['bam_index']=idval['reduced_index']
                align_outdic[sample][refID]['stats']=idval['reduced_stats']

    stdout_print("END reduce",printInLog=True)
                
    return align_outdic

def launch_samtools_stats(list_stats,title=""):
    stats_jblst = utils.Jobslist("bam stats for all")
    for bamfile,statsfile in list_stats:
        #stats command
        stats_cmd = "{path}samtools stats {inbam} | grep ^SN | cut -f2-3 > {outstats}".format(path=samtools_path,inbam=bamfile,outstats=statsfile)
        stats_jblst.add_a_job(stats_cmd,"samtools stats {}".format(bamfile),statsfile)

    thread_list=[]
    if args.tinyverbose:
        t3 = threading.Thread(target=progress_bar,args=['Samtools stats progress {}'.format(title),len(stats_jblst.get_joblist()),thread_list])
        t3.start()
        
    stdout_print("Launch samtools stats {}".format(title),printInLog=True)
    utils.trun(args, stats_jblst, thrd=thread_list)

    if args.tinyverbose:
        wait_progressbar_finished(t3)

def launch_samtools_index(list_index,title=""):
    index_jblst = utils.Jobslist("bam index for all")
    
    for bamfile,indexfile in list_index:
        #index command
        index_cmd= "{path}samtools index {inbam}".format(path=samtools_path,inbam=bamfile)
        index_jblst.add_a_job(index_cmd,"samtools index {}".format(bamfile),indexfile)

    thread_list=[]
    if args.tinyverbose:
        t2 = threading.Thread(target=progress_bar,args=['Samtools index progress {}'.format(title),len(index_jblst.get_joblist()),thread_list])
        t2.start()
        
    stdout_print("Launch samtools index {}".format(title),printInLog=True)
    utils.trun(args, index_jblst, thrd=thread_list)

    if args.tinyverbose:
        wait_progressbar_finished(t2)
        
def launch_samtools_sort(list_sort,title=""):
    sort_jblst = utils.Jobslist("bam sort for all")

    for bamfile,sortfile in list_sort:
        #sort command
        sort_cmd = "{path}samtools sort {inbam} -o {outbam}".format(path=samtools_path,inbam=bamfile,outbam=sortfile)
        sort_jblst.add_a_job(sort_cmd,"samtools sort {}".format(bamfile),sortfile)

    thread_list=[]
    if args.tinyverbose:
        t1 = threading.Thread(target=progress_bar,args=['Samtools sort progress {}'.format(title),len(sort_jblst.get_joblist()),thread_list])
        t1.start()
        
    stdout_print("Launch samtools sort {}".format(title),printInLog=True)
    utils.trun(args, sort_jblst, thrd=thread_list)

    if args.tinyverbose:
        wait_progressbar_finished(t1)
    
def check_alignments(align_outdic):
    errlist=[]
    for sample,refAlnDatas in align_outdic.items():
        for refID, aligndatas in refAlnDatas.items():
            errVal=is_stderrFileInError(aligndatas['errBWTfile'])
            align_outdic[sample][refID]['align_ok']=not(errVal)
            errlist.append(errVal)
            if errVal:
                addTextToLogFile("\tERROR in alignment {} VS {}".format(sample,refID))
    return sum(errlist),align_outdic
            
def progress_bar(title,nbr,thread_list):
    finished=0
    all_thrd=[]
    while finished<nbr:
        time.sleep(2)
        if len(thread_list)>0:
            for obj in thread_list:
                if obj not in all_thrd:
                    all_thrd.append(obj)
            #finished=sum(1 for f in file_list if os.path.exists(f) and os.path.getsize(f)>0 )
            finished=utils.getNumberJobTerminated(all_thrd)
            percent=round(100*float(finished)/float(nbr),1)

            progstr="="*(int(percent/2)-1)+">"+" "*(50-int(percent/2))
            if percent>=100:
                progstr="="*int(percent/2)

            printProgress("{}: [{}] {}%".format(title,progstr,percent))
    print("")

def is_stderrFileInError(stderrFile):
    with open(stderrFile,'r') as sf:
        fcontent=sf.read()
        return "(ERR)" in fcontent
    
def wait_progressbar_finished(pbThrd):
    while(pbThrd.is_alive()):
        time.sleep(2)

def printProgress(txt):
    sys.stdout.write("\r"+str(txt))
    sys.stdout.flush()
    
def run_kmerRefFilter(kmerSize,readsConfig):
    stdout_print(" -- kmerRefFilter --",printInLog=True)

    outKRF_dic={}
    
    outFQfiltered=config['OutFolders']['out_filtered']
    yaml_resume_file=os.path.join(outFQfiltered,config['KRF_yamlFile'])
    
    out_template_FWD="{}FWD_filtered.fastq"
    out_template_REV="{}REV_filtered.fastq"
    out_template_Single="{}_filtered.fastq"

    out_log_PK_KRF=os.path.join(config['log_folder'],"kmerRefFilter_pickle_{dt}.txt".format(dt=config['file_time']))
    out_log_KRF=os.path.join(config['log_folder'],"kmerRefFilter_log_{dt}.txt".format(dt=config['file_time']))

    #Create kmerRefFilter pickle database
    stdout_print("\tconstruct kmer DB for filtering : {}".format(config['kmer_Ref_PickleFile']),printInLog=True)
    cmd="{pth}kmerRefFilter.py -nolog -P {cpu} -p {outpickle} -k {ksize} -r {fastaref} -o . 2>{log} 1>&2".format(pth=kmerrefilter_path,outpickle=config['kmer_Ref_PickleFile'],cpu=config['MaxParallelsJobs'],fastaref=args.refdatabase,log=out_log_PK_KRF,ksize=kmerSize)
    if not os.path.exists(config['kmer_Ref_PickleFile']) or args.force:
        os.system(cmd)
    
    jblst = utils.Jobslist("kmerRefFilter for all")
    for sample,values in readsConfig.items():
        outKRF_dic[sample]={}
        #Check if read already filtered is set
        extlist=[".fastq",".fq",".gz",".bz2"]
        if not values['params']['filtered']:
            readsConfig[sample]['params']['filtered']=True
            sample_reads=readsConfig[sample]['reads']
            readsConfig[sample]['reads']=[]
            if values['params']['format']=='paired':
                PairedList=get_paired_reads_from_list(sample_reads)
                for pair in PairedList:
                    fq_fwd_bname=get_fileNoExt(pair[0],extlist)
                    fq_rev_bname=get_fileNoExt(pair[1],extlist)

                    out_filt_fwd=os.path.join(outFQfiltered,out_template_FWD.format(fq_fwd_bname))
                    out_filt_rev=os.path.join(outFQfiltered,out_template_REV.format(fq_rev_bname))

                    readsConfig[sample]['reads'].append(out_filt_fwd)
                    readsConfig[sample]['reads'].append(out_filt_rev)
                    cmd="{pth}kmerRefFilter.py -nolog -y -k {ksize} -1 {FwdFq} -2 {RevFq} -i {pkref} -o {outdir} 2>>{log} 1>>{yamlKmer}".format(pth=kmerrefilter_path,FwdFq=pair[0],RevFq=pair[1],outdir=outFQfiltered,pkref=config['kmer_Ref_PickleFile'],log=out_log_KRF,yamlKmer=yaml_resume_file,ksize=kmerSize)
                    addTextToLogFile("\tkmerRefFilter paired: {}".format(cmd))
                    jblst.add_a_job(cmd,"kmerRefFilter paired {}".format(sample),out_filt_fwd)
                    
            elif values['params']['format']=='single':
                for fqfile in sample_reads:
                    fq_bname=get_fileNoExt(fqfile,extlist)
                    out_filt_single=os.path.join(outFQfiltered,out_template_Single.format(fq_bname))

                    readsConfig[sample]['reads'].append(out_filt_single)
                    cmd="{pth}kmerRefFilter.py -nolog -y -k {ksize} -f {Fq} -i {pkref} -o {outdir} 2>>{log} 1>>{yamlKmer}".format(pth=kmerrefilter_path,Fq=fqfile,outdir=outFQfiltered,pkref=config['kmer_Ref_PickleFile'],log=out_log_KRF,yamlKmer=yaml_resume_file,ksize=kmerSize)
                    addTextToLogFile("\tkmerRefFilter single: {}".format(cmd))
                    jblst.add_a_job(cmd,"kmerRefFilter paired {}".format(sample),out_filt_single)

    thread_list=[]
    if args.tinyverbose:
        t1 = threading.Thread(target=progress_bar,args=['KmerRefFilter progress',len(jblst.get_joblist()),thread_list])
        t1.start()
        
    #Launch kmerRefFilter process
    utils.trun(args, jblst,thrd=thread_list)

    if args.tinyverbose:
        wait_progressbar_finished(t1)

    stdout_print(" -- END kmerRefFilter --",printInLog=True)

    return readsConfig

def get_fileNoExt(filename,extlist):
    filebname=os.path.basename(filename)
    for ext in extlist:
        ext_index=filebname.upper().rfind(ext.upper())
        if ext_index>0:
            filebname=filebname[:ext_index]
    return filebname

def get_paired_reads_from_list(reads_list):
    PairedList=[]
    pair=[None,None]
    for c,reads in enumerate(reads_list):
        pair[c%2]=reads
        if (c%2)==1:
            PairedList.append(pair)
            pair=[None,None]
    return PairedList
    
def split_reads(readsConfig,compress=True):
    sample_infile=[]
    outfld=config['OutFolders']['out_splited']
    for sample,values in readsConfig.items():
        if values['params']['split']:
            for fqfile in values['reads']:
                sample_infile.append((sample,fqfile))
            readsConfig[sample]['reads']=[]
            readsConfig[sample]['params']['split']=False

    outlog_list=[]
    if len(sample_infile)>0:
        #Parallelize split function
        #splitfastq.split(sample_infile,outfld,outcompress=compress)
        p=Pool(config['MaxParallelsJobs'])
        partial_split_reads=partial(splitfastq.split,outfld=outfld,outcompress=compress)
        split_rslt=[]
        psync=p.map_async(partial_split_reads,sample_infile,callback=split_rslt.extend)
        p.close()
        #Wait for split termination
        while not psync.ready():
            time.sleep(1)

        #get and format results from Pool
        for rslt in split_rslt:
            sample=rslt['sample']
            readsConfig[sample]['reads'].append(rslt['fqFWD'])
            readsConfig[sample]['reads'].append(rslt['fqREV'])
            outlog_list.append(rslt['log'])
            
    ld.yaml_dumpDic(config['splited_reads'],readsConfig)

    for logline in outlog_list:
        for line in logline.strip().split('\n'):
            addTextToLogFile(line)
            
    return readsConfig

def index_references(pk_refdatabase):
    stdout_print(" -- Index reference file --",printInLog=True)
    index_dic={}
    bwt_out = os.path.join(config['log_folder'],config['bwt_build_out'])
    yaml_out=os.path.join(config['OutFolders']['configs_fld'],config['index_cfg'])

    if os.path.exists(yaml_out) and not args.force:
        return yaml_out
    
    #Create single fasta files to index
    refDB_dic=ld.pickle_loadDic(pk_refdatabase)
    for seqID,seqVal in refDB_dic.items():
        singlefasta=os.path.join(config['OutFolders']['single_fasta_ref_fld'],"{}.fa".format(seqID))
        indexfile=os.path.join(config['OutFolders']['bwt_IndexdirResult'],seqID)
        if not seqVal['combined']:
            index_dic[seqID]={'fasta':singlefasta,'index':indexfile}
            
        if (not os.path.exists(singlefasta) or args.force) and not seqVal['combined']:
            with open(singlefasta,'w') as outfasta:
                outfasta.write(">{}\n{}".format(seqID,seqVal['seq']))
                
    #Index all fasta files
    jblst = utils.Jobslist("bowtie2-build for all")
    for seqID,indexval in index_dic.items():
        target = "{}{}".format(indexval['index'],'.1.bt2')
        cmd = '{bpath}bowtie2-build {reffile} {refname} >>{bwtout}'.format(bpath=bowtie2_path,reffile = indexval['fasta'],refname = indexval['index'],bwtout = bwt_out)
        addTextToLogFile("\tBowtie2-build: {}".format(cmd))
        jblst.add_a_job(cmd,"bowtie2-build {}".format(seqID),target)
    ld.yaml_dumpDic(yaml_out,index_dic)
    utils.trun(args, jblst)
    return yaml_out

def parse_referenceDB(infastafile):
    stdout_print(" -- parse reference file {} --".format(infastafile),printInLog=True)
    dictfasta={}
    allelePart={}
    subAlleleList=[]
    fastafile_basename=os.path.basename(infastafile)
    
    #TO REMOVE YAML FILE
    out_yamlfile=os.path.join(config['OutFolders']['reference_DB_fld'],"{}.yaml".format(fastafile_basename[:fastafile_basename.rfind('.')]))
    
    out_pkfile=os.path.join(config['OutFolders']['reference_DB_fld'],"{}.pk".format(fastafile_basename[:fastafile_basename.rfind('.')]))
    #Parameters: Paralog, specie, grpRef
    Allowed_parameters={'Paralog':'bool', 'specie':'str', 'grpRef':'str','allelePart':'list','grpRefPart':'list'}

    if os.path.exists(out_pkfile) and not args.force:
        return out_pkfile
    
    reffasta=SeqIO.parse(infastafile,'fasta')
    for rec in reffasta:
        seqID_list=str(rec.id).split('|')
        seqID=seqID_list[0]
        dictfasta[seqID]={}
        dictfasta[seqID]['seq']=str(rec.seq)
        dictfasta[seqID]['len_seq']=len(str(rec.seq))
        dictfasta[seqID]['combined']=False
        for p,t in Allowed_parameters.items():
            dictfasta[seqID][p]=convert_value(None,t)
        
        if len(seqID_list)>1:
            for values in seqID_list[1:]:
                if values.count("=")==1:
                    tmp=values.split('=')
                    param_name=tmp[0]
                    param_val=tmp[1]

                    if param_name in Allowed_parameters.keys():
                        dictfasta[seqID][param_name]=convert_value(param_val,Allowed_parameters[param_name])
    for refname,refparam in dictfasta.items():
        grpID,haploID=parse_GrpRef(refparam['grpRef'],refname)
        if grpID is not None and haploID is not None:
            refparam['grpRef']="H{}-{}".format(grpID,haploID)
        dictfasta[refname]['grpIDPart']=None
        dictfasta[refname]['haploIDPart']=None
        
        if dictfasta[refname]['grpRefPart'] is not None and len(dictfasta[refname]['grpRefPart'])>0:
            subAlleleList.append(refname)
            dictfasta[refname]['grpID']=None
            dictfasta[refname]['haploID']=None
            dictfasta[refname]['grpIDPart']=[]
            dictfasta[refname]['haploIDPart']=[]
            for grpRefPart in dictfasta[refname]['grpRefPart']:
                grpIDPart,haploIDPart=parse_GrpRef(grpRefPart,refname)
                dictfasta[refname]['grpIDPart'].append(grpIDPart)
                dictfasta[refname]['haploIDPart'].append(haploIDPart)
        else:
            dictfasta[refname]['grpID']=grpID
            dictfasta[refname]['haploID']=haploID

    for subAllele in subAlleleList:
        for allele,group in zip(dictfasta[subAllele]['allelePart'],dictfasta[subAllele]['grpRefPart']):
            if allele not in allelePart.keys():
                allelePart[allele]={'subAllele':[],'len_seq':0,'subGrp':None,'subGrpID':None,'subHaploID':None}
            allelePart[allele]['subAllele'].append(subAllele)
            allelePart[allele]['subGrp']=group
            grpID,haploID=parse_GrpRef(group,allele)
            allelePart[allele]['subGrpID']=grpID
            allelePart[allele]['subHaploID']=haploID
            allelePart[allele]['len_seq']+=dictfasta[subAllele]['len_seq']
            
    for allele, alleles_val in allelePart.items():
        dictfasta[allele]={'combined':True,'Paralog': False, 'allelePart': None, 'grpID': alleles_val['subGrpID'], 'grpIDPart': None, 'grpRef': alleles_val['subGrp'], 'grpRefPart': None, 'haploID': alleles_val['subHaploID'], 'haploIDPart': None, 'len_seq': alleles_val['len_seq'], 'seq': None, 'specie': None}
        
        
    #TO REMOVE YAML FILE
    ld.yaml_dumpDic(out_yamlfile,dictfasta)
    
    ld.pickle_dumpDic(out_pkfile,dictfasta)
        
    return out_pkfile
            
def parse_GrpRef(grpRef,refname):
    grpID=None
    haploID=None
    if grpRef and grpRef[0].upper()=="H":
        if grpRef.count('-')==1:
            tmp=grpRef[1:].split('-')
            grpID=int(tmp[0])
            haploID=int(tmp[1])
        elif grpRef.count('-')==0:
            grpID=int(grpRef[1])
            haploID=int(grpRef[2:])
        else:
            stdout_print("WARNING grpRef parameter '{}' for reference '{}' not correctly formated".format(grpRef,refname),printInLog=True)
    return grpID,haploID
            
        
def convert_value(value,vtype):
    if vtype=='bool':
        if value=="True": return True
        if value=="False": return False
        return bool(value)
    if value is not None:
        if vtype=='str':
            return str(value)
        elif vtype=='int':
            return int(value)
        elif vtype=='float':
            return float(value)
        elif vtype=='list':
            return list(value.split(','))

    return value
        
def load_readsList(readsList):
    readlist_error=False
    rslt = {}
    reads_dir=""
    Allowed_params={'format':'str','split':'bool','filtered':'bool'}
    default_value={'format':'single','split':False,'filtered':False}
    with open(readsList,'r') as rl:
        for read in rl:
            #nosplit value compatibility
            read=read.replace("nosplit=False","split=True")
            read=read.replace("nosplit=True","split=False")
            commentID = read.find('#')
            if commentID<0:
                tmp = read.strip()
            elif commentID>0:
                tmp = read[:commentID].strip()
            if commentID!=0:
                read_info = tmp.split(',')
                if len(read_info)>=2:
                    readlist_error=(reads_dir=="")
                    infos_dict={}
                    sample_name=read_info[0]
                    infos_dict['reads']=[]
                    infos_dict['params']={}
                    for p,t in Allowed_params.items():
                        infos_dict['params'][p]=convert_value(default_value[p],t)
                    for val in read_info[1:]:
                        if val.count("=")==1:
                            tmp_val=val.split('=')
                            if tmp_val[0] in Allowed_params.keys():
                                infos_dict['params'][tmp_val[0]]=convert_value(tmp_val[1],Allowed_params[tmp_val[0]])
                        else:
                            readfile=os.path.join(reads_dir,val)
                            if os.path.exists(readfile):
                                infos_dict['reads'].append(readfile)
                            else:
                                stdout_print("ERROR '{}' not exists".format(readfile),printInLog=True)
                                readlist_error=True
                                
                    rslt[sample_name] = infos_dict
                elif os.path.isdir(tmp):
                    reads_dir=tmp
    if readlist_error:
        exitProg("ERROR in read file '{}'".format(readsList))
        
    return rslt
    
def check_FileExtensions(filename, extList):
    fileExt = filename[filename.rfind('.')+1:]
    if fileExt.upper() in [ext.upper() for ext in extList]:
        return True
    return False
    
def files_exits(file_list):
    errfiles=[]
    for f in [v for v in file_list if v]:
        if not os.path.exists(f):
            errfiles.append(f)
    if len(errfiles)>0:
        print("ERROR - files listed not exists:")
        print("\n".join(errfiles))
        print("EXIT...")
        sys.exit(1)
        
def addTextToLogFile(logtxt,printTime=True):
    if printTime:
        ttime = "{}\t".format(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"))
    else:
        ttime=""
	
    logFolder = config['log_folder']
    
    dir_create(logFolder)
    
    logFile = os.path.join(logFolder,config['log_file'])
    lf = open(logFile,'a+')
    if logtxt.count('\n')>0:
        for l in logtxt.split('\n'):
            lf.write("{tm}{log}\n".format(tm=" "*len(ttime),log=l))
    else:
        lf.write("{tm}{log}\n".format(tm=ttime,log=logtxt))
    lf.close()

def dir_create(path,printOnScreen=True):
    if not os.path.exists(path):
        dirmsg = "-- create {} directory".format(path)
        if printOnScreen:
            stdout_print(dirmsg)
        try:
            os.makedirs(path)
        except:
            exitProg("ERROR can't create destination directory {}".format(path))

def stdout_print(txt, printInLog=False):
    if args.tinyverbose:
        print(txt)
    if printInLog:
        addTextToLogFile(txt)

def exitProg(msg=''):
    if msg!='':
        addTextToLogFile("PROGRAM EXIT\t{}".format(msg))
    sys.exit(msg)
