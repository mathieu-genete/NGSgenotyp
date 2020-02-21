#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import os

AppPath = os.path.realpath(__file__)
App_Folder = AppPath[:AppPath.rfind('/')]

import yaml
import pickle
import subprocess
import numpy as np
import argparse
import threading
import gzip
import time
from datetime import datetime
import multiprocessing
from Bio import SeqIO
from Bio import Seq
from operator import itemgetter
from ete3 import Tree, TreeStyle, TextFace
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

__version__ = "1.0.1"
AppName = "haploAsm"

args = None
config = None
spades_folder = None
quast_folder = None
yass_folder = None
cap3_folder = None
muscle_folder = None
phyml_folder = None
bowtie_folder = None
samtools_folder=None

cpu_number =  multiprocessing.cpu_count()

def run(ArgsVal):

	global args
	global spades_folder
	global quast_folder
	global yass_folder
	global cap3_folder
	global muscle_folder
	global phyml_folder
        global bowtie_folder
        global samtools_folder
	global config
	global cpu_number

	interlaced = False
	paired = False
	
	description = """ NGSgenotyp Haplotyp Assembly """

	parser = argparse.ArgumentParser(prog="haploAsm",description=description)
	
	parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
	parser.add_argument("-f","--force",help="force",action="store_const", const="-f", default="")
	parser.add_argument("-i","--concatenated",help="orignals reads files with concatenated forward and reverse paired-end reads",action="store_const", const="-i", default="")
	parser.add_argument("-p","--paired",help="orignals reads files with paired-end reads forward and reverse files distincts",action="store_const", const="-p", default="")
	parser.add_argument("-d","--diploid",help="use with diploid highly polymorphic genomes",action="store_const", const="-d", default="")
	parser.add_argument("-x","--excludeParalogs",help="exclude paralogs from yass results",action="store_const", const="-x", default="")
	parser.add_argument("-y","--phylotree",help="generate phylogenetic tree",action="store_const", const="-y", default="")
        parser.add_argument("-yx","--excludeParaPhylo",help="exclude paralogs from phylogeny trees",action="store_const", const=True, default=False)

	parser.add_argument("-t","--threaded",help="define threads number used for phylogeny. By default 10%% of cpu numbers",type=int)	
	parser.add_argument("-st","--spadesThreads",help="define threads number used for spades. By default 50%% of cpu numbers",type=int)	
	parser.add_argument("-db","--refdatabase", help="reference database used for genotyping", required = True)
	parser.add_argument("-k","--kmerSize", help="comma-separated list of k-mer sizes for Spades (must be odd and less then 128)")
	parser.add_argument("-o","--outputdir", help="output directory for assembly", required = True)
	parser.add_argument("-l","--indivlist", help="list of haplotype name for assembly")
	parser.add_argument("-m","--contigminlen", help="keep contigs up from minimum lenght",type=int)
	parser.add_argument("-M","--contigmaxlen", help="keep contigs down to maximum lenght",type=int)
        parser.add_argument("-mdc","--maxDistContigs", help="contig maximum distance from nearest node - default = see maxDistContigs value in configuration file",type=float)
	parser.add_argument("-filFQ","--filteredFQ",help="Configuration file with filtered reads informations", required = True)
	parser.add_argument("-s","--samviewstats", help="samview_stats file from genotyp results", required = True)
	parser.add_argument("--config", help="config file", default=os.path.join(App_Folder,'haploAsm_config.yaml'))
	
	args = parser.parse_args(ArgsVal)
	
	config = yaml.load(open(args.config,"r"))

	outfolder = args.outputdir

	config['outfolder'] = outfolder

        config['logFolder'] = os.path.join(config['outfolder'],'log')

        config['currentDate'] = datetime.now().strftime("%d_%m_%Y-%H_%M_%S")

        config['log_file'] = "{appn}_Log_{dt}.txt".format(appn=AppName,dt=config['currentDate'])

	bwt_dirResult = config['bwt_dirResult']
	
	spades_folder = path_InConfig('spades',App_Folder)
	quast_folder = path_InConfig('quast',App_Folder)
	yass_folder = path_InConfig('yass',App_Folder)
	cap3_folder = path_InConfig('cap3',App_Folder)
        muscle_folder = path_InConfig('muscle',App_Folder)
	phyml_folder = path_InConfig('phyml',App_Folder)
        bowtie_folder = path_InConfig('bowtie',App_Folder)
        samtools_folder = path_InConfig('samtools',App_Folder)

	print "\n\n\t{} -- v{}\n\n".format(description,__version__)

	#log file write
        currentDateTime = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        cmdline = " ".join(sys.argv)
	startText = "{appn} {appv} -- {cdate}".format(appn=AppName,appv=__version__,cdate=currentDateTime)
	starDateline = "*"*len(startText)
	addTextToLogFile("{}\n{}\n{}".format(starDateline,startText,starDateline),False,False)
	addTextToLogFile("user command: {}".format(cmdline),False)
	
        if args.maxDistContigs:
                config['maxDistContigs']=args.maxDistContigs
                addTextToLogFile("maxDistContigs user define: {}".format(config['maxDistContigs']))

        if checkX11():
                addTextToLogFile("X11 display detected...")
        else:
                addTextToLogFile("WARNING: Can't print PDF files, should use X11 display...")

        if args.diploid:
                dipStr="diploid_"
        else:
                dipStr=""
                
	if args.paired:
		paired = True
		
	if args.concatenated:
		interlaced = True
		paired = False
		
        readsInfoFilFQ_filename = args.filteredFQ
	
	ref_database = args.refdatabase
	print "database: {}\n".format(ref_database)

        paralogs_List = get_paralogs_list(ref_database)
        if args.excludeParaPhylo:
                addTextToLogFile("Remove Paralogs on phylogenetics trees : {}".format(" ".join(paralogs_List)))

	if not os.path.exists(outfolder):
		os.mkdir(outfolder)
			
        ReadsFilFQPath,ReadsFilFQInfos = get_readsConfigs(readsInfoFilFQ_filename)
	
	samviewstatsFName = args.samviewstats
	
	AnalyseFolder = samviewstatsFName[:samviewstatsFName.rfind('/',0,samviewstatsFName.rfind('/'))]
	
	BWTResult_Folder = os.path.join(AnalyseFolder,bwt_dirResult)
	
	fileExt = samviewstatsFName[samviewstatsFName.rfind('.')+1:]
	
	if (fileExt!='p'):
		samviewstatsFile = yaml.load(open(samviewstatsFName,"r"))
	else:
		samviewstatsFile = pickle.load(open(samviewstatsFName,"r"))
	
	AlignSelect = {}
	
	if args.indivlist and os.path.exists(args.indivlist):
		lf = open(args.indivlist,"r")
		IndList = lf.read()
		HaploNames = IndList.split('\n')
		HaploNames.remove('')
		lf.close()
	else:
                HaploNames = ReadsFilFQInfos.keys()
                
	
	for hname in HaploNames:
		AlignSelect[hname]={}
		for ref, val in samviewstatsFile[hname].items():
			if val['mean cover']>0 and val['IsParalog']==False:
				AlignSelect[hname][ref]=val
				
	addTextToLogFile("{} sample(s) to analyse:".format(len(AlignSelect.keys())))

        addTextToLogFile("Samples list:",False)
	for sample in AlignSelect.keys():
		print "\t - {}".format(sample)
                addTextToLogFile("\t{}".format(sample),False)

        thrdPhyList = []
	for IndivName,Refs in AlignSelect.items():
                t1 = launch_analyse(IndivName,Refs,outfolder,BWTResult_Folder,ReadsFilFQPath,ReadsFilFQInfos,interlaced,paired,ref_database,paralogs_List)
                if t1:
                        thrdPhyList.append(t1)
        
        if len(thrdPhyList)>0:
                if args.threaded:
                        ThrdNbr = args.threaded
                else:
                        ThrdNbr = (cpu_number * 10)/100
                if ThrdNbr<=0:
                        ThrdNbr = 1

                addTextToLogFile("Phylogeny Threaded: {} -- {} samples to launch".format(ThrdNbr,len(thrdPhyList)))

                c=0
                for thrd in thrdPhyList:
                        while get_NbrAliveThreads(thrdPhyList)>=ThrdNbr:
                                time.sleep(1)
			thrd.start()
                        c+=1
                        addTextToLogFile("launch thread {}/{}".format(c,len(thrdPhyList)))

                while get_NbrAliveThreads(thrdPhyList)>0:
                        time.sleep(1)

                AllTruncContigsFasta = ""
                for IndivName in AlignSelect.keys():
                        indivFolder = os.path.join(outfolder,IndivName)
                        phyloFolder = os.path.join(indivFolder,"{}_phylogeny".format(IndivName))
                        treefileName = os.path.join(phyloFolder,"MUSCLE_{}_ContigsRef.phy_phyml_tree.txt".format(IndivName))
                        yass_contigsFileName = os.path.join(indivFolder,"{idn}_{dip}yass_oriented_contigs.fasta".format(idn=IndivName,dip=dipStr))
                        yass_contigsTruncFileName = os.path.join(indivFolder,"{idn}_{dip}yass_truncated_contigs.fasta".format(idn=IndivName,dip=dipStr))
                        AllTruncContigsFasta += DrawPhyloTree(ref_database,IndivName,treefileName,yass_contigsFileName,yass_contigsTruncFileName,phyloFolder)

                refFastaFile = SeqIO.parse(ref_database,"fasta")
                for rec in refFastaFile:
                        id = str(rec.id)
                        seq = str(rec.seq)

                        limitPos = id.find('|')

                        if limitPos>0:
                                seqName = id[:limitPos]
                        else:
                                seqName = id

                        if args.excludeParaPhylo:
                                if seqName not in paralogs_List:
                                        AllTruncContigsFasta += ">{}_length_{}\n{}\n".format(seqName,len(seq),seq)
                        else:
                                AllTruncContigsFasta += ">{}_length_{}\n{}\n".format(seqName,len(seq),seq)

                AllTruncContigsFileName = os.path.join(outfolder,"AllTruncContigsRefs.fasta")
                AllTruncContigsFile = open(AllTruncContigsFileName,"w")
                AllTruncContigsFile.write(AllTruncContigsFasta)
                AllTruncContigsFile.close()

                if len(AlignSelect.keys())>1:
                        addTextToLogFile("Launch phylogeny for ALL contigs")
                        Pylo_From_Fasta(AllTruncContigsFileName,outfolder)

        addTextToLogFile("-- ENDED --")

def get_paralogs_list(ref_database):
        paralogsList=[]
        fastaF=SeqIO.parse(ref_database,'fasta')
        for rec in fastaF:
                idseq=str(rec.id)
                refname=idseq.split("|")[0]
                if "Paralog=1" in idseq:
                        paralogsList.append(refname)
        return paralogsList

def path_InConfig(configKey,App_Folder):
        if config[configKey]!=None:
                return os.path.join(App_Folder,config[configKey])
        else:
                return ''

def Pylo_From_Fasta(inFasta,outfolder):

        MuscleFilename = os.path.join(outfolder,"MUSCLE_AllContigsRef.fasta")
        PhylipFilename = os.path.join(outfolder,"MUSCLE_AllContigsRef.phy")

        if not os.path.exists(MuscleFilename) or args.force:
                cmd = "{path}muscle -in {infasta} -out {outMsl}".format(path=muscle_folder,infasta=inFasta,outMsl=MuscleFilename)
                os.system(cmd)
        else:
                addTextToLogFile("Muscle alignent already exists -- PASS THIS STEP")

        if fastaToPhylip(MuscleFilename,PhylipFilename):
                if not os.path.exists("{}_phyml_tree.txt".format(PhylipFilename)) or os.path.getsize("{}_phyml_tree.txt".format(PhylipFilename))==0 or args.force:
                        cmd = "{path}phyml --quiet -i {phy} -d nt".format(path=phyml_folder,phy=PhylipFilename)
                        os.system(cmd)
                else:
                        addTextToLogFile("Phyml tree already exists -- PASS THIS STEP")

                treefileName = "{}_phyml_tree.txt".format(PhylipFilename)
                outname = "AllContigs_Phylo"

                if checkX11():
                        DrawSimplePhyloTree(treefileName,outname,outfolder)

def DrawSimplePhyloTree(treefileName,outName,outFolder):

        t = Tree(treefileName)
        ts = TreeStyle()
        ts.scale = 500
        ts.branch_vertical_margin = 10
        ts.show_leaf_name = False
        ts.show_branch_support = True
        ts.title.add_face(TextFace("Phylogeny for {}".format(treefileName), fsize=12), column=0)
        for leaf in t.iter_leaves():
                if "Contig" in leaf.name:
                        name_face = TextFace(leaf.name, fgcolor="red")
                else:
                        name_face = TextFace(leaf.name, fgcolor="black")

                leaf.add_face(name_face, column=0, position='branch-right')

        render_Phylo(os.path.join(outFolder,outName),t,ts)

def DrawPhyloTree(ref_database,indivName,treefileName,yass_contigsFileName,yass_contigsTruncFileName,outFolder):

    if not os.path.exists(outFolder):
        return ""

    refID = get_refID(ref_database)
    groupIds = get_GroupIds(refID)
    colorId=config['IDgroupColors']
    paralogColor="#D8D8D8"
    outContigs = open(os.path.join(outFolder,"Contigs_candidates.txt"),"w")
    fastaContigs = open(os.path.join(outFolder,"Contigs_candidates.fasta"),"w")
    outFastaContigsIds = []
    addTextToLogFile(yass_contigsFileName)
    outName = "MUSCLE_{}_ContigsRef_phy_phyml".format(indivName)

    t = Tree(treefileName)
    ts = TreeStyle()
    ts.scale = 500
    ts.branch_vertical_margin = 10
    ts.show_leaf_name = False
    ts.show_branch_support = True
    ts.title.add_face(TextFace("Phylogeny for {}".format(treefileName), fsize=12), column=0)
    keep_nodes = list()   
    leaves_Annotation(t,refID)
    contigId = ['C']

    for gId in groupIds:
        vals=[str(gId)]+contigId
        keep_nodes = list(set(Group_Monophyl(t,vals,colorId[gId])+keep_nodes))

    
    paralog_nodes = list(set(Group_Monophyl(t,['P']+contigId,paralogColor)))
    AllNodes = list(set(t))

    #ContigList = [{'name':v.name,'dist':v.dist} for v in keep_nodes if "Contig" in v.name]
    #ContigList = [{'name':v.name,'dist':v.dist} for v in AllNodes if ("Contig" in v.name and v.name not in [pnd.name for pnd in paralog_nodes])]

    MonophyleticContigList = [{'name':v.name,'dist':v.dist,'allele':'MONOPHYL. GPR'} for v in keep_nodes if "Contig" in v.name]
    ContigList = get_ContigsCandidates(t,ref_database,config['maxDistContigs'])
    ContigListName = [v['name'] for v in ContigList]
    ContigList += [v for v in MonophyleticContigList if v['name'] not in ContigListName]

    SortedContigList = sorted(ContigList, key=lambda k: k['dist'])
    if os.path.exists(yass_contigsFileName):
        contFasta = SeqIO.parse(yass_contigsFileName,"fasta")
        contigs_names = ["{}_{}".format(indivName,v['name']) for v in SortedContigList]
        for rec in contFasta:
                if rec.id in contigs_names:
                        fastaLine = ">{}\n{}\n".format(rec.id,rec.seq)
                        outFastaContigsIds.append(rec.id)
                        fastaContigs.write(fastaLine)
        outTxt = ""
        outTxt += "File: {}\n".format(treefileName)
        outTxt += "max dist value: {}\n".format(config['maxDistContigs'])
        for v in SortedContigList:
                outTxt += "{}\tnearest node:{}\tdist={}\n".format(v['name'],v['allele'],v['dist'])
        outTxt += "\n"
        addTextToLogFile(outTxt)
        outContigs.write(outTxt)

    keep_nodes = list(set(Group_Monophyl(t,['P']+contigId,paralogColor)+keep_nodes))
    
    if checkX11():
        if not os.path.exists(os.path.join(outFolder,"{}.pdf".format(outName))) or args.force:
                render_Phylo(os.path.join(outFolder,"{}_ALL".format(outName)),t,ts)
                #t.prune(keep_nodes)
                #render_Phylo(os.path.join(outFolder,outName),t,ts)        
        else:
                addTextToLogFile("{}.pdf already exists -- PASS THIS STEP".format(outName))
                    
    
    outContigs.close()
    fastaContigs.close()
    
    outTruncContigs = ""
    TruncFname = SeqIO.parse(yass_contigsTruncFileName,"fasta")
    for rec in TruncFname:
            if rec.id in outFastaContigsIds:
                    outTruncContigs += ">{}\n{}\n".format(rec.id,rec.seq)
    return outTruncContigs

def get_ContigsCandidates(tree,ref_database,maxDistContigs):
        allelsList=[]
        contigList=[]

        ref = SeqIO.parse(ref_database,"fasta")
        for rec in ref:
                id=str(rec.id)
                if "Paralog=1" not in id:
                        allelsList.append(id.split('|')[0])

        allelTree = [j for j in tree if j.name in allelsList]

        for i in tree:
                if "Contig" in i.name:
                        tmpDist=10
                        alleleName=""
                        for j in allelTree:
                                if i.name!=j.name:
                                        d = tree.get_distance(i.name,j.name)
                                        if d<tmpDist:
                                                tmpDist=d
                                                alleleName=j.name
                        if tmpDist<=maxDistContigs:
                                contigList.append({'name':i.name,'allele':alleleName,'dist':tmpDist})

        return contigList

def checkX11():
        if 'DISPLAY' in os.environ.keys():
                return True

        return False

def render_Phylo(outPdf, inTree, TreeStyle):
    inTree.sort_descendants()
    OutName = "{}.pdf".format(outPdf)
    inTree.render(OutName,tree_style=TreeStyle)

def get_GroupIds(refID):
    groupIds = list()
    for k,v in refID.items():
        if 'gprId' in v.keys():
            groupIds = list(set([v['gprId']]+groupIds))
    return groupIds

def Group_Monophyl(inTree,monoValues,color=""):

    NodesOut = list()

    tm = inTree.get_monophyletic(values=monoValues,target_attr="gprId")
    
    for node in tm:
        if color!="":
            node.img_style["bgcolor"]=color
        [NodesOut.append(v) for v in node]

    return NodesOut

def get_refID(ref_database):
    refDB = SeqIO.parse(args.refdatabase,"fasta")

    refID = {}
    for rec in refDB:
        tmp = str(rec.id).split('|')
        refID[tmp[0]]={v.split('=')[0]:v.split('=')[1] for v in tmp[1:]}

    return refID

def leaves_Annotation(inTree,refID):

    for leaf in inTree.iter_leaves():
        leaf.add_feature('gprId',"X")

        if "Contig" in leaf.name:
                leaf.add_features(isContig=True)
                leaf.add_feature('gprId',"C")
                name_face = TextFace(leaf.name, fgcolor="red")
        else:
                name_face = TextFace(leaf.name, fgcolor="black")
                leaf.add_features(isContig=False)
                if leaf.name in refID.keys():
                    for k,v in refID[leaf.name].items():
                        leaf.add_feature(k,v)
                        if k=='Paralog':
                            leaf.add_feature('gprId','P')

        leaf.add_face(name_face, column=0, position='branch-right')
        if leaf.gprId in get_GroupIds(refID):
            leaf.add_face(TextFace(" - Grp {}".format(leaf.gprId), fgcolor="black"), column=1, position='branch-right')


def get_NbrAliveThreads(thrdList):
        aliveNbr = 0
        for thrd in thrdList:
                if thrd.is_alive():
                        aliveNbr += 1
        return aliveNbr

def get_readsConfigs(readsInfo_filename):

	if check_FileExtensions(readsInfo_filename,['yaml','yml']):
		readsConfig = yaml.load(open(readsInfo_filename))
	else:
		readsConfig = load_readsList(readsInfo_filename)
	
	ReadsPath = readsConfig['reads_dir']
	ReadsInfos = {d[0]:d[1:] for d in readsConfig['reads_files'].values()}
	
	return ReadsPath,ReadsInfos
	
def launch_analyse(IndivName,Refs,outfolder,BWTResult_Folder,ReadsFilFQPath,ReadsFilFQInfos,interlaced,paired,ref_database,paralogs_List):

	outIndivFolder = os.path.join(outfolder,IndivName)
	
	if not os.path.exists(outIndivFolder):
		os.mkdir(outIndivFolder)
	
	Aligns_Folder = os.path.join(BWTResult_Folder,IndivName)
	addTextToLogFile(TitleFrame("Assembly for {}".format(IndivName)))
	addTextToLogFile("\tAlignFolder: {}\n\treads path: {}\n\treads files: {}\n\tSelected mapped reads to Ref : {}\n".format(Aligns_Folder, ReadsFilFQPath, ",".join(ReadsFilFQInfos[IndivName]), ",".join(Refs.keys())))
	
	ALL_mapped = ""
		
        All_filtered = ""
                
        ALL_filtered_fileName = os.path.join(outIndivFolder,"ALL_filtered_{}.fastq".format(IndivName))

        FWD_filtered_fileName = os.path.join(outIndivFolder,"ALL_filtered_{}_R1.fastq".format(IndivName))
        REV_filtered_fileName = os.path.join(outIndivFolder,"ALL_filtered_{}_R2.fastq".format(IndivName))
                
        if os.path.exists(ALL_filtered_fileName) and not args.force:
                addTextToLogFile("\tfile \'{}\' allready exists\n".format(ALL_filtered_fileName))
                mfile = open(ALL_filtered_fileName,"r")
                All_filtered = mfile.read()
                mfile.close()
                ALL_mapped = All_filtered
        else:
                parity = 0
                for read in ReadsFilFQInfos[IndivName]:
                        rdFile = open(os.path.join(ReadsFilFQPath,read),"r")
                        fileFastqContent = rdFile.read()
                        All_filtered += fileFastqContent
                        rdFile.close()

                        if paired:
                                if parity%2==0:
                                        outPairedFilename = FWD_filtered_fileName
                                else:
                                        outPairedFilename = REV_filtered_fileName
                                
                                fltFile = open(outPairedFilename,"a")
                                fltFile.write(fileFastqContent)
                                fltFile.close()
                                parity+=1
                                
                if interlaced:
                        ALL_mapped, fwd_mapped, rev_mapped = split_fastq(All_filtered)

                        fltFile = open(FWD_filtered_fileName,"a")
                        fltFile.write(fwd_mapped)
                        fltFile.close()

                        fltFile = open(REV_filtered_fileName,"a")
                        fltFile.write(rev_mapped)
                        fltFile.close()

                else:
                        ALL_mapped = All_filtered
                                
                fltFile = open(ALL_filtered_fileName,"a")
                fltFile.write(ALL_mapped)
                fltFile.close()
                        
                InSpadesFastq = ALL_filtered_fileName
		
	if args.kmerSize:
		usedKmerSize = args.kmerSize
	else:
		usedKmerSize = "21,41,81"

	if args.spadesThreads:
                spThrd = args.spadesThreads
        else:
                spThrd = (cpu_number * 50/100)

        if interlaced or paired:
                inputSpadesFQ = "-1 {fwdFQ} -2 {revFQ}".format(fwdFQ=FWD_filtered_fileName,revFQ=REV_filtered_fileName)
        else:
                inputSpadesFQ = "-s {fastq}".format(fastq=InSpadesFastq)

        spadesLogName = "spades_log_file_{ind}_{dt}.txt".format(ind=IndivName,dt=config['currentDate'])
        spadesLogFile = os.path.join(config['logFolder'],spadesLogName)

        addTextToLogFile("LAUNCH Assembly for {} -- spades log file: {}".format(IndivName,spadesLogName))

	if args.diploid:
                addTextToLogFile("Diploid mode for assembly")
		spadesOutFolder = os.path.join(outIndivFolder,"DipSpades_OUT")
		
		if not os.path.exists(spadesOutFolder) and not args.force:
			cmd = "{path}dipspades.py -t {thrds} --careful -k {kmerS} -o {outf} {fastq} > {logf}".format(path=spades_folder,outf=spadesOutFolder,kmerS=usedKmerSize,fastq=inputSpadesFQ,thrds=spThrd,logf=spadesLogFile)
                        addTextToLogFile("dipspades command: {}".format(cmd))
			os.system(cmd)
                else:
                        addTextToLogFile("Assembly already exists -- PASS ASSEMBLY STEP")
		
		contigsFastaFile = os.path.join(spadesOutFolder,"dipspades/consensus_contigs.fasta")
		
	else:
	
		spadesOutFolder = os.path.join(outIndivFolder,"Spades_OUT")
	
		if not os.path.exists(spadesOutFolder) and not args.force:
			cmd = "{path}spades.py -t {thrds} --careful -k {kmerS} -o {outf} {fastq} > {logf}".format(path=spades_folder,outf=spadesOutFolder,kmerS=usedKmerSize,fastq=inputSpadesFQ,thrds=spThrd,logf=spadesLogFile)
                        addTextToLogFile("spades command: {}".format(cmd))
			os.system(cmd)
                else:
                        addTextToLogFile("Assembly already exists -- PASS ASSEMBLY STEP")
	
		contigsFastaFile = os.path.join(spadesOutFolder,"contigs.fasta")
	
		if not os.path.exists(contigsFastaFile):
			kmers = usedKmerSize.split(',')
			kmers.sort(reverse=True)
		
			for k in kmers:
				contigsFastaFileKmer = os.path.join(spadesOutFolder,"K{}/final_contigs.fasta".format(k))
			
				if os.path.exists(contigsFastaFileKmer):
					kfile = open(contigsFastaFileKmer,"r")
					cfile = open(contigsFastaFile,"w")
					cfile.write(kfile.read())
					cfile.close()
					kfile.close()
					break
                else:
                        addTextToLogFile("{} already exists -- PASS THIS STEP".format(contigsFastaFile))

	if os.path.exists(contigsFastaFile):
	
                fastaAllContigs,fastaUnfilteredAllContigs = apply_cap3(contigsFastaFile,IndivName)

		tmpFastaFile = os.path.join(outIndivFolder,"tmp.fasta")
                outFastaContigs = open(tmpFastaFile,"w")
		outFastaContigs.write(fastaAllContigs)
		outFastaContigs.close()
                fastaAllContigs,NoReturn = apply_cap3(tmpFastaFile,IndivName,True)
                os.remove(tmpFastaFile)
		
                if fastaAllContigs=="":
                        addTextToLogFile("WARNING: NO CONTIGS ASSEMBLED OR SELECTED FOR {}".format(IndivName))
                        return False

		dipStr = ''
		
		if args.diploid:
			dipStr = "diploid_"
			
			QuastFolder = os.path.join(outIndivFolder,"quast_dip_results")
		else:
		
			QuastFolder = os.path.join(outIndivFolder,"quast_results")

		outFastaContigsFilename = os.path.join(outIndivFolder,"{iname}_{dip}contigs.fasta".format(iname=IndivName,dip=dipStr))
		outFastaContigs = open(outFastaContigsFilename,"w")
		outFastaContigs.write(fastaAllContigs)
		outFastaContigs.close()
		
		outFastaUnfilteredContigsFilename = os.path.join(outIndivFolder,"{iname}_{dip}contigs_unfiltered.fasta".format(iname=IndivName,dip=dipStr))
		outFastaUnfilteredContigs = open(outFastaUnfilteredContigsFilename,"w")
		outFastaUnfilteredContigs.write(fastaUnfilteredAllContigs)
		outFastaUnfilteredContigs.close()
		
                if not os.path.exists(os.path.join(QuastFolder,"report.txt")) or args.force:
                        cmd = "{path}quast.py -o {outd} {contigf} 1> /dev/null 2> /dev/null".format(path=quast_folder,outd=QuastFolder,contigf=outFastaUnfilteredContigsFilename)
                        addTextToLogFile("quast command: {}".format(cmd))
                        os.system(cmd)
		else:
                        addTextToLogFile("Quast files already exists -- PASS THIS STEP")

                yass_stdout = os.path.join(outIndivFolder,"yass_stdout")
                yass_stderr = os.path.join(outIndivFolder,"yass_stderr")
                if not os.path.exists(yass_stdout) or args.force:
                        cmd = "{path}yass -d 2 {params}{contigFile} {refFile} > {ystdout} 2> {ystderr}".format(ystderr=yass_stderr,ystdout=yass_stdout,path=yass_folder,params=config['yass_params'],contigFile=outFastaContigsFilename,refFile=ref_database)
                        addTextToLogFile("yass command: {}".format(cmd))
                        os.system(cmd)
                else:
                        addTextToLogFile("Yass files already exists -- PASS THIS STEP")

		ystdoutFile = open(yass_stdout,"r")
                stdout = ystdoutFile.read()
                ystdoutFile.close()

                YassOut = stdout.split('\n')
		YassOut.remove('')
		
		refDict = SeqIO.to_dict(SeqIO.parse(ref_database, "fasta"))
		
		YassHead = [y.strip() for y in YassOut[0].split(',')]
		YassHead[0]='Query id'
		YassHead.append('RefLen')
		YassHead.append('Error rate')
		YassHead.append('% ref cov')
		#YassHead.append('Distance')
		YassHead.append('Contig orientation')
		YassHead.append('% contig cov')

                SubjectID_head=['Ref name','HaploId','gprId','specie','grpRef','Paralog']
                for hd in SubjectID_head:
                        YassHead.append(hd)

		YassRslt = []
		FwRvInfo = {}
                QstartEnd = {}

		# Yass Fields: [0] Query id,[1] Subject id, [2] % identity, [3] alignment length, [4] mismatches, [5] gap openings, [6] q. start, [7] q. end, [8] s. start, [9] s. end,[10] e-value, [11] bit score

		for YassVal in YassOut[1:]:
			YassLine = [check_num(yl.strip()) for yl in YassVal.split('\t')]
			YassLine.append(len(refDict[YassLine[1]].seq))
			
			errorRate = float(YassLine[4])/float(YassLine[3])
			YassLine.append(errorRate)
			
			YassLine.append(round((float(YassLine[3])/float(YassLine[12]))*100,2))
			
			# distance = 1 - ((%cov/100) * (1 - Error rate))
			#distance = 1-((float(YassLine[3])/float(YassLine[12]))*(1-errorRate))
			#YassLine.append(distance)
			
			orient = YassLine[9]-YassLine[8]
                        if YassLine[0] not in QstartEnd.keys():
                                QstartEnd[YassLine[0]]={'start':[],'end':[]}
                        QstartEnd[YassLine[0]]['start'].append(YassLine[6])
                        QstartEnd[YassLine[0]]['end'].append(YassLine[7])

			if YassLine[0] not in FwRvInfo.keys():
				FwRvInfo[YassLine[0]] = 0
				
			if orient>0:
				YassLine.append('Fwd')
				FwRvInfo[YassLine[0]] += 1
			else:
				YassLine.append('Rev')
				FwRvInfo[YassLine[0]] -= 1
                        contig_name = YassLine[0]
                        contig_len = float(contig_name[contig_name.rfind('_')+1:])
                        contig_ref = float(YassLine[3])/contig_len
                        YassLine.append(round(contig_ref*100,2))
                        
                        SubjectID = split_SeqID(YassLine[1])
                        for hd in SubjectID_head:
                                if hd in SubjectID.keys():
                                        YassLine.append(SubjectID[hd])
                                else:
                                        YassLine.append("")
                        
			YassRslt.append({YassHead[i]:YassLine[i] for i in range(0,len(YassHead),1)})
			
		sorted_YassRslt = sorted(YassRslt, key=itemgetter('Error rate'))
		#sorted_YassRslt = sorted(YassRslt, key=itemgetter('Distance'))
		
		contigsAligned = []
		out_Head = ['Query id','Subject id','Ref name','gprId','% identity','Error rate','% ref cov','% contig cov','Contig orientation','q. start','q. end','s. start','s. end','RefLen','alignment length','mismatches','gap openings','e-value','bit score']
		out_Yass = ",".join(out_Head)+"\n"
		for l in sorted_YassRslt:
			if (('Paralog=1' not in l['Subject id']) or not args.excludeParalogs) and l['% ref cov']>config['cov_filter']:
				out_Yass += ",".join(["\"{}\"".format(str(l[v]).replace('.',config['float_sep'])) for v in out_Head])+"\n"
				contigsAligned.append(str(l['Query id']))
		
		out_yassFilename = os.path.join(outIndivFolder,"{iname}_{dip}ALign_Yass.csv".format(iname=IndivName,dip=dipStr))
		outYass = open(out_yassFilename,"w")
		outYass.write(out_Yass)
		outYass.close()
		
		orientedYassOutFastaContigs = ""
                truncatedYassOutFastaContigs = ""
		colWidth = 80
		
		for record in SeqIO.parse(outFastaContigsFilename, "fasta"):
			Fid = record.id
			seq = record.seq
			revseq = record.seq.reverse_complement()

			if Fid in FwRvInfo.keys():

                                qstart = min(QstartEnd[Fid]['start'])-1
                                qend = max(QstartEnd[Fid]['end'])
                                trucSeq = Seq.Seq(str(seq)[qstart:qend])
                                revTrucSeq = trucSeq.reverse_complement()

				if FwRvInfo[Fid] > 0:
			
					orientedYassOutFastaContigs += ">{}\n{}".format(Fid,fasta_SeqWidth(seq,colWidth))
                                        truncatedYassOutFastaContigs += ">{}\n{}".format(Fid,fasta_SeqWidth(trucSeq,colWidth))
				elif FwRvInfo[Fid] < 0:
			
					orientedYassOutFastaContigs += ">{}_C\n{}".format(Fid,fasta_SeqWidth(revseq,colWidth))
                                        truncatedYassOutFastaContigs += ">{}_C\n{}".format(Fid,fasta_SeqWidth(revTrucSeq,colWidth))
				
				else:
			
					orientedYassOutFastaContigs += ">{}\n{}".format(Fid,fasta_SeqWidth(seq,colWidth))
					orientedYassOutFastaContigs += ">{}_C\n{}".format(Fid,fasta_SeqWidth(revseq,colWidth))
                                        truncatedYassOutFastaContigs += ">{}\n{}".format(Fid,fasta_SeqWidth(trucSeq,colWidth))
                                        truncatedYassOutFastaContigs += ">{}_C\n{}".format(Fid,fasta_SeqWidth(revTrucSeq,colWidth))
					
		orientedOutFastaContigsFilename = os.path.join(outIndivFolder,"{iname}_{dip}yass_oriented_contigs.fasta".format(iname=IndivName,dip=dipStr))
		orientedOutFastaContigs = open(orientedOutFastaContigsFilename,"w")
		orientedOutFastaContigs.write(orientedYassOutFastaContigs)
		orientedOutFastaContigs.close()

                truncatedOutFastaContigsFilename = os.path.join(outIndivFolder,"{iname}_{dip}yass_truncated_contigs.fasta".format(iname=IndivName,dip=dipStr))
                truncatedOutFastaContigs = open(truncatedOutFastaContigsFilename,"w")
                truncatedOutFastaContigs.write(truncatedYassOutFastaContigs)
                truncatedOutFastaContigs.close()

                contig_quality(IndivName,outfolder,orientedOutFastaContigsFilename,ALL_filtered_fileName,FWD_filtered_fileName,REV_filtered_fileName,interlaced,paired)

                if args.phylotree:
                        t1 = threading.Thread(target=phylogenetic_tree,args=(IndivName,outfolder,truncatedOutFastaContigsFilename,list(set(contigsAligned)),ref_database,paralogs_List,))
                        t1.setDaemon(True)
                        return t1
	else:
	
		addTextToLogFile("\tWARNING: NO \'contigs.fasta\' file found for {}".format(IndivName))

        return False

def contig_quality(IndivName,outfolder,orientedOutFastaContigsFilename,ALL_filtered_fileName,FWD_filtered_fileName,REV_filtered_fileName,interlaced,paired):

        addTextToLogFile("Contig Quality for {}".format(IndivName),False)
        Indiv_folder = os.path.join(outfolder,IndivName)
        align_folder = os.path.join(Indiv_folder,"align_qual")
        ref_File = os.path.join(align_folder,"{}_oriented_Contigs_REF".format(IndivName))
        bwtStderr = os.path.join(align_folder,"bwt_stderr")
        bamFile = os.path.join(align_folder,"{}_VS_ContigRef.bam").format(IndivName)
        bamSortedFile = os.path.join(align_folder,"{}_VS_ContigRef_SORTED.bam").format(IndivName)

        statsFile = os.path.join(align_folder,"samtools_stats.txt")
        out_pdf = os.path.join(align_folder,"qual_report.pdf")
        out_csv = os.path.join(align_folder,"qual_report.csv")
        out_prot =  os.path.join(align_folder,"max_translated_seq.fa")

        if not os.path.exists(align_folder):
                os.mkdir(align_folder)

        contigsNames = list()
        contigsProtInfos={}
        
        contigsFile = SeqIO.parse(orientedOutFastaContigsFilename,"fasta")

        for rec in contigsFile:
                contigsNames.append(str(rec.id))
                contigsProtInfos[str(rec.id)]=get_maxProtInformations(str(rec.seq))

        #index ref
        if not os.path.exists("{}.1.bt2".format(ref_File)):
                cmd = "{path}bowtie2-build {ref} {outRef}".format(path=bowtie_folder,ref=orientedOutFastaContigsFilename,outRef=ref_File)
                os.system(cmd)
        else:
                addTextToLogFile("Bowtie2 reference index already exist -- PASS THIS STEP")
        
        bwtAlign = ""
        smtFlag = ""
        ISrange_STR = ""

        if interlaced or paired:
                ISrange = min_max_insertSizeFromAlignment(FWD_filtered_fileName,REV_filtered_fileName,ref_File,align_folder)
                ISrange_STR = " (minIS: {},maxIS: {})".format(ISrange['min'],ISrange['max'])
                bwtAlign = "-I {minIS} -X {maxIS} -1 {FWDfq} -2 {REVfq}".format(FWDfq=FWD_filtered_fileName,REVfq=REV_filtered_fileName,minIS=ISrange['min'],maxIS=ISrange['max'])
                #smtFlag = "-f 0x2"
                smtFlag = "-F 0x4"
        else:
                bwtAlign = "-U {fq}".format(fq=ALL_filtered_fileName)
                smtFlag = "-F 0x4"
        if not os.path.exists(out_csv):
                cmd = "{path}bowtie2 --phred33 -x {ref} {bwtAln} 2> {stderrOut} | {pth}samtools view -S -b {flg} - > {outBam}".format(path=bowtie_folder,ref=ref_File,bwtAln=bwtAlign,stderrOut=bwtStderr,outBam=bamFile,flg=smtFlag,pth=samtools_folder)
                os.system(cmd)

                cmd = "{pth}samtools sort {inBam} > {sortBam}".format(inBam=bamFile,sortBam=bamSortedFile,pth=samtools_folder)
                os.system(cmd)

                cmd = "{pth}samtools index {inBam}".format(inBam=bamSortedFile,pth=samtools_folder)
                os.system(cmd)
        
                cmd = "{pth}samtools stats {inBam} | grep ^SN | cut -f2-3 > {outStats}".format(inBam=bamSortedFile,outStats=statsFile,pth=samtools_folder)
                os.system(cmd)
        
                covDatas = {}
                statsDatas = {}
                statsPPDatas = {}
                for refName in contigsNames:
                        covDatas[refName] = bamDepthCoverage(bamSortedFile,refName,interlaced,paired)
                        statsDatas[refName] = bam_stats(bamSortedFile,refName)
                        if interlaced or paired:
                                statsPPDatas[refName] = bam_stats(bamSortedFile,refName,"-f 0x2")
                        else:
                                statsPPDatas[refName] = statsDatas[refName]

                plot_QALreports(out_pdf,out_csv,out_prot,covDatas,statsDatas,statsPPDatas,interlaced,paired,ISrange_STR,contigsProtInfos)

def get_maxProtInformations(inSeq):
    bestSeqPhase=[]
    tmpSizeFragm=0
    protStart=0
    protEnd=0
    protSeq=""
    return_Phase=0
    strand="+"
    seqSize=len(inSeq)
    protSize=0
    for Allphase in range(0,6):
        phase = Allphase%3
        if Allphase>2:
            seq = str(Seq.Seq(inSeq).reverse_complement())
        else:
            seq=inSeq
        NbrOfNMod=abs(len(seq[phase:])%-3)
        seqPhase = seq[phase:]+("N"*NbrOfNMod)
        prot = Seq.translate(seqPhase)
        prot = prot.replace("*","*_")
        protSplited=prot.split("_")
        tabLenProts = [len(p) for p in protSplited]
        tabSeqProts = [p for p in protSplited]
        nbrTransSeq = len(protSplited)
        maxFragSize = max(tabLenProts)
        if maxFragSize>tmpSizeFragm:
            tmpSizeFragm = maxFragSize
            bestSeqPhase=tabLenProts
            return_Phase=phase
            protStart=sum(bestSeqPhase[:bestSeqPhase.index(max(bestSeqPhase))])*3+phase
            protEnd=protStart+max(bestSeqPhase)*3-NbrOfNMod
            if Allphase>2:
                strand="-"
                protStart = seqSize - protStart
                protEnd = seqSize - protEnd
            protSeq = tabSeqProts[bestSeqPhase.index(max(bestSeqPhase))].strip("X")
            protSize = float(len(protSeq))

    percentOfTranslate=round((protSize*3)/float(seqSize)*100,2)
    return {'pStart':protStart,'pEnd':protEnd,'pSeq':protSeq,'phase':return_Phase,'strand':strand,'TransPct':percentOfTranslate,'dnaSize':seqSize}

def min_max_insertSizeFromAlignment(FQFwd,FQRev,Reference,outResultsFolder):

        Rslt = {}
        bamFile = os.path.join(outResultsFolder,"pre_align.bam")

        cmd = "{path}bowtie2 --phred33 -x {ref} -1 {FWDfq} -2 {REVfq} | {pth}samtools view -S -b -F 0x4 - > {outBam}".format(path=bowtie_folder,ref=Reference,outBam=bamFile,FWDfq=FQFwd,REVfq=FQRev,pth=samtools_folder)
        os.system(cmd)
        if os.path.exists(bamFile):
                stats = bam_stats(bamFile,"",flag="")
                insertAvg = float(stats['insertSizeAVG'])
                insertStdDev = float(stats['insertStdDev'])
                stdDevFactor = config['qualInsertSizeDevFactor']
                minIS = int(insertAvg - (insertStdDev * stdDevFactor))
                maxIS = int(insertAvg + (insertStdDev * stdDevFactor))
                Rslt = {'min':minIS,'max':maxIS}
                addTextToLogFile("Contigs Quality - insert Size informations for {},{}: minimum = {} -- maximum = {}".format(FQFwd,FQRev,minIS,maxIS),False)
                os.remove(bamFile)
                return Rslt
        else:
                return False

def get_zeroDepthLength(ydatas):
        start=False
        rslt=list()
        sval=0
        for i in range(0,len(ydatas)-1):
                if ydatas[i]!=0 and ydatas[i+1]==0 and start==False:
                        sval=i+1
                        start=True
                if ydatas[i]==0 and ydatas[i+1]!=0 and start:
                        start=False
                        rslt.append(i-sval+1)
        return rslt

def plot_QALreports(out_pdf,out_csv,out_prot,covDatas,statsDatas,statsPPDatas,interlaced,paired,ISrange_STR,contigsProtInfos):

        if interlaced or paired:
                modifTitle="Properly paired"
        else:
                modifTitle="Unpaired"
        
        protFasta=""
        csvRslt = "flag,Contig,Coverage,mean depth,insert size average,bases mapped,mismatchs,error rate,raw sequences,reads mapped,reads properly paired,inward oriented pairs,outward oriented pairs,pairs with other orientation,pairs on different contigs,inZero depth count,inZero depth total size, Max Translate Region percent,Prot Size (AA)\n"
        with PdfPages(out_pdf) as pdf:
                for refName in sorted(covDatas.keys()):
                        flag = ""
                        cov = covDatas[refName]
                        contigSize=contigsProtInfos[refName]['dnaSize']
                        basesMapped = int(statsPPDatas[refName]['basesMapped'])
                        mismatchs = int(statsPPDatas[refName]['mismatchs'])

                        protSeq = contigsProtInfos[refName]['pSeq']
                        protFasta +=">{rf}_{plen}_{phase}_{strand}\n{pseq}\n".format(rf=refName,plen=len(protSeq),phase=contigsProtInfos[refName]['phase'],strand=contigsProtInfos[refName]['strand'],pseq=protSeq)
                        if basesMapped>0:
                                errRate = round(100*float(mismatchs)/float(basesMapped),1)

                                PPratio = round(100*float(statsDatas[refName]['properlyPaired'])/float(statsDatas[refName]['mapped']),1)

                                xdatas = cov['xdepth']
                                ydatas = cov['ydepth']
                        
                                xzero = [xdatas[i] for i in range(0,len(xdatas)) if ydatas[i]==0]
                                revCumSum = cov['revCumSum']
                                inZeroD = get_zeroDepthLength(ydatas)
                                inZeroD_count = len(inZeroD)
                                inZeroD_sum = sum(inZeroD)

                                pairsDifferentChr = int(statsDatas[refName]['pairsDifferentChr'])
                                pairsOtherOrient = int(statsDatas[refName]['pairsOtherOrient'])
                                insertSizeAvg_STR = "{}{}".format(statsPPDatas[refName]['insertSizeAVG'],ISrange_STR)

                                if inZeroD_count>0 or errRate>5 or pairsOtherOrient>0 or pairsDifferentChr>0:
                                        flag="W"

                                fig = plt.figure()

                                plt.clf()
                                ax = fig.add_subplot(111)
                                ax.set_ylim([-0.05,int(max(ydatas)*1.1)])
                                ax.set_xlim([-0.05,contigSize])
                                ax.set_title("{} ({})".format(refName,modifTitle), fontsize=8)
                                ax.annotate('', xytext=(contigsProtInfos[refName]['pStart'],1), xy=(contigsProtInfos[refName]['pEnd'],1),arrowprops={'arrowstyle': '->', 'lw': 4, 'color': 'blue', 'alpha': 0.5})
                                plt.plot(xdatas,ydatas,color='silver')
                                plt.scatter(xzero,[cov['mean']]*len(xzero),color='salmon')
                                plt.plot([0,max(xdatas)],[cov['mean'],cov['mean']],color='lightgreen',linestyle="--")
                                plt.text(min(xdatas),int(max(ydatas)*0.5),"{title} stats:\n coverage: {cvr} %\n mean depth: {mean}\n insert size average: {ins}\n bases mapped: {bm}\n mismatchs: {mis}\n error rate: {er} %\n\nAlignment stats:\n raw seq.: {rws}\n mapped: {mapp}\n properly paired: {pp} ({ppr}%)\n pairs with other orientation: {oop}\n pairs on diff. contig: {dcp}\nmax translate region: {MTP}".format(title=modifTitle,ins=insertSizeAvg_STR,rws=statsDatas[refName]['rawSeq'],mapp=statsDatas[refName]['mapped'],pp=statsDatas[refName]['properlyPaired'],ppr=PPratio,oop=pairsOtherOrient,dcp=pairsDifferentChr,bm=basesMapped,mis=mismatchs,er=errRate,cvr=round(cov['pSeqCov']*100,1),mean=round(cov['mean'],2),MTP="{} % ({} AA)".format(contigsProtInfos[refName]['TransPct'],len(protSeq))),color='navy',size=7)
                                pdf.savefig(plt.gcf())

                        else:
                                PPratio = 0.0
                                errRate = 1.0
                                inZeroD_count = 0
                                inZeroD_sum = 0
                                flag = "E"

                        csvRslt += "{flg},{ref},{cvr},{mean},{ins},{bm},{mis},{er},{rws},{mapp},{pp},{iw},{ow},{oop},{dcp},{inZc},{inZs},{MTP},{pSizeAA}\n".format(ref=refName,ins=statsPPDatas[refName]['insertSizeAVG'],rws=statsDatas[refName]['rawSeq'],mapp=statsDatas[refName]['mapped'],pp=statsDatas[refName]['properlyPaired'],ppr=PPratio,oop=statsDatas[refName]['pairsOtherOrient'],dcp=statsDatas[refName]['pairsDifferentChr'],bm=basesMapped,mis=mismatchs,er=errRate,cvr=round(cov['pSeqCov']*100,1),mean=round(cov['mean'],2),iw=statsDatas[refName]['inwardOriented'],ow=statsDatas[refName]['outwardOriented'],inZc=inZeroD_count,inZs=inZeroD_sum, flg=flag,MTP=contigsProtInfos[refName]['TransPct'],pSizeAA=len(protSeq))

                csvFile = open(out_csv,"w")
                csvFile.write(csvRslt)
                csvFile.close()

                fastaProt = open(out_prot,"w")
                fastaProt.write(protFasta)
                fastaProt.close()

def bam_stats(bamSortedFile,refName,flag=""):
        cmd = "{pth}samtools stats {fl} {bam} {rfName} | grep ^SN | cut -f 3".format(bam=bamSortedFile,rfName=refName,fl=flag,pth=samtools_folder)
        job = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        cmdRslt = job.communicate()
        statsRslt = cmdRslt[0].split('\n')
        statsDatas = {'rawSeq':statsRslt[0],'mapped':statsRslt[6],'mappedPaired':statsRslt[7],'properlyPaired':statsRslt[9],'basesMapped':statsRslt[17],'mismatchs':statsRslt[20],'pairsOtherOrient':statsRslt[29],'pairsDifferentChr':statsRslt[30],'insertSizeAVG':statsRslt[25],'insertStdDev':statsRslt[26],'inwardOriented':statsRslt[27],'outwardOriented':statsRslt[28]}
        return statsDatas

def bamDepthCoverage(bamfile,refName,interlaced,paired):

        if interlaced or paired:
                flag = "-f 0x2 "
        else:
                flag = ""

        cmd = "{pth}samtools view -b {fl}{bam} | {pth}samtools depth -a - | grep {rfName} | cut -f 2-".format(fl=flag,bam=bamfile,rfName=refName,pth=samtools_folder)
        job = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        depthrstl = job.communicate()

	xdepth = []
	ydepth = []

	for line in depthrstl[0].split("\n"):
		if len(line.split("\t"))>1:
			xdepth.append(int(line.split("\t")[0]))
			ydepth.append(int(line.split("\t")[1]))
	
	fxdepth = xdepth
	fydepth = ydepth

	if (len(ydepth)>0 and sum(ydepth)>0):
	
		mean = sum(ydepth)/float(len(ydepth))
		gmean = sum(fydepth)/float(len(fydepth))
		pSeqCov = sum([1 for v in fydepth if v>0])/float(len(fydepth))
		
		revCumSum = []
                revCumSum_Pos = []
		for i in range(0,int(max(ydepth))):
			pos = 0
			tmp = 0
			for val in ydepth:
				if val >= i:
					tmp += 1
				pos += 1
                        revCumSum.append(float(tmp)/float(len(ydepth)))
                        revCumSum_Pos.append(tmp)

		covfreq = getFrequence(ydepth)
		
		maxdepth = int(max(covfreq[0]))

		med = np.median(ydepth)

		q75, q25 = np.percentile(ydepth, [75 ,25])

		iqr = q75 - q25

		riqr=iqr/maxdepth
		rmed=med/maxdepth
	
	else:
		mean = 0
		gmean = 0
		pSeqCov = 0
		xdepth = []
		xdepth = []
		revCumSum = []
                revCumSum_Pos = []
		fxdepth = []
		fydepth = []
		covfreq = []
		med = 0
		q25 = 0
		q75 = 0
		maxdepth = 0
		iqr = 0
		riqr = 0
		rmed = 0
	
	return {'mean':mean,'gmean':gmean,'pSeqCov':pSeqCov,'xdepth':xdepth,'ydepth':ydepth,'revCumSum':revCumSum,'fxdepth':fxdepth,'fydepth':fydepth,'covfreq':covfreq,'median':med,'q25':q25,'q75':q75,'maxdepth':maxdepth,'iqr':iqr, 'riqr':riqr, 'rmed':rmed,'revCumSum_Pos':revCumSum_Pos}

def getFrequence(data,excludeZero=False):

	if excludeZero:
		data = remove_values_from_list(data,0)
	
	x = np.array(data)
	
	unique, counts = np.unique(x, return_counts=True)
	
	return [unique, counts]
		
def phylogenetic_tree(IndivName,outfolder,orientedOutFastaContigsFilename,contigsAligned,ref_database,paralogs_List):
        addTextToLogFile("Phylogenetic tree for {}".format(IndivName))
        IndivFloder = os.path.join(outfolder,IndivName)
	phylOutFolder = os.path.join(IndivFloder,"{}_phylogeny".format(IndivName))
        phyloFastaFilename = os.path.join(phylOutFolder,"{}_ContigsRef.fasta".format(IndivName))
        MuscleFilename = os.path.join(phylOutFolder,"MUSCLE_{}_ContigsRef.fasta".format(IndivName))
        PhylipFilename = os.path.join(phylOutFolder,"MUSCLE_{}_ContigsRef.phy".format(IndivName))

	if not os.path.exists(phylOutFolder):
		os.mkdir(phylOutFolder)
	
	phyloFasta = ""
	
	refDbFasta = SeqIO.parse(ref_database,"fasta")
	
	for rec in refDbFasta:
                tmp = str(rec.id)

                if tmp.find('|')>0:
                        seqId = tmp[:tmp.find('|')]
                else:
                        seqId = tmp

                seq = str(rec.seq)

                if args.excludeParaPhylo:
                        if seqId not in paralogs_List:
                                phyloFasta += ">{}\n{}\n".format(seqId,seq)
                else:
                        phyloFasta += ">{}\n{}\n".format(seqId,seq)

        contigFasta = SeqIO.parse(orientedOutFastaContigsFilename,"fasta")
        
        for rec in contigFasta:
                tmp = str(rec.id)
                seqId = tmp[tmp.find("Contig"):]
                seq = str(rec.seq)

                if str(rec.id).rstrip("_C") in contigsAligned:
                        addTextToLogFile("Contig add for phylo: {}".format(rec.id))
                        phyloFasta += ">{}\n{}\n".format(seqId,seq)

        if not os.path.exists(phyloFastaFilename) or args.force:
                outphyloFasta = open(phyloFastaFilename,"w")
                outphyloFasta.write(phyloFasta)
                outphyloFasta.close()
        else:
                addTextToLogFile("{} already exist -- PASS THIS STEP".format(phyloFastaFilename))
        
        if not os.path.exists(MuscleFilename) or args.force:
                cmd = "{path}muscle -in {infasta} -out {outMsl}".format(path=muscle_folder,infasta=phyloFastaFilename,outMsl=MuscleFilename)
                addTextToLogFile("muscle command: {}".format(cmd))
                os.system(cmd)
        else:
                addTextToLogFile("Muscle file already exist -- PASS THIS STEP")

        if fastaToPhylip(MuscleFilename,PhylipFilename):
                if not os.path.exists("{}_phyml_tree.txt".format(PhylipFilename)) or os.path.getsize("{}_phyml_tree.txt".format(PhylipFilename))==0 or args.force:
                        cmd = "{path}phyml --quiet -i {phy} -d nt".format(path=phyml_folder,phy=PhylipFilename)
                        addTextToLogFile("phyml command: {}".format(cmd))
                        os.system(cmd)
                else:
                        addTextToLogFile("Phyml tree file already exist -- PASS THIS STEP")

def fastaToPhylip(inFasta,outPhy):
        fastaFile = SeqIO.parse(inFasta,"fasta")
        phyOut = ""
        seqLen=-1
        N = 0
        for rec in fastaFile:
                N += 1
                seq = str(rec.seq)
                phyOut += "{}\t{}\n".format(rec.id,seq)
                if seqLen>=0 and len(seq)!=seqLen:
                        addTextToLogFile("ERROR fasta sequences lentgh should be identical...")
                        return False
                seqLen = len(seq)

        of = open(outPhy,"w")
        of.write("{} {}\n".format(N,seqLen))
        of.write(phyOut)
        of.close()
        return True

def split_SeqID(seqID):
        tmp = seqID.split('|')
        rslt = {}
        rslt['Ref name']=tmp[0]
        for val in tmp[1:]:
                sval = val.split('=')
                rslt[sval[0]]=sval[1]
        return rslt

def apply_cap3(contigsFastaFile,IndivName,delCapFiles=False):

        ContigsCapFile = contigsFastaFile + ".cap.contigs"
        SingletsCapFile = contigsFastaFile + ".cap.singlets"
                
        if not os.path.exists(ContigsCapFile) and not args.force:
                cmd = "{path}cap3 {cfile}".format(path=cap3_folder,cfile=contigsFastaFile)
                os.system(cmd)
        else:
                addTextToLogFile("cap3 files already exist -- PASS THIS STEP")

        fastaAllContigs,fastaUnfilteredAllContigs = merge_fasta([ContigsCapFile,SingletsCapFile],IndivName)

        if delCapFiles:
                cmd = "rm {}.cap.*".format(contigsFastaFile)
                os.system(cmd)

        return fastaAllContigs,fastaUnfilteredAllContigs

def check_num(val):
	if val.replace('.','',1).isdigit():
		if val.find('.')>0:
			return float(val)
		else:
			return int(val)
	return val
	
def merge_fasta(fastaFiles,IndivName):
        addTextToLogFile("MERGE Fasta files for {}: {}".format(IndivName,",".join(fastaFiles)),False)
	AllContigs = ""
	Unfiltered_AllContigs = ""
	colWidth=80
	contigNbr = 1
	
	if args.contigminlen:
		contigminlen = args.contigminlen
	else:
		contigminlen = 0

	if args.contigmaxlen:
		contigmaxlen = args.contigmaxlen
	else:
		contigmaxlen = 0
        
	
	fastaList = []
	for faFile in fastaFiles:
                if os.path.exists(faFile):
                        for record in SeqIO.parse(faFile, "fasta"):
                                Unfiltered_AllContigs += ">{}\n".format(record.id)
                                Unfiltered_AllContigs += fasta_SeqWidth(record.seq,colWidth)
                                if len(record.seq)>=contigminlen:
                                        if contigmaxlen==0 or len(record.seq)<=contigmaxlen:
                                                fastaList.append({'seq':record.seq, 'seqlen':len(record.seq)})
        if len(fastaList)>0:
                sorted_fastaList = sorted(fastaList, key=itemgetter('seqlen'),reverse=True)
	
                for v in sorted_fastaList:
                        AllContigs += ">{}_Contig_N{}_length_{}\n".format(IndivName,contigNbr,v['seqlen'])
                        AllContigs += fasta_SeqWidth(v['seq'],colWidth)
                        contigNbr += 1
				
	return AllContigs,Unfiltered_AllContigs
	
def fasta_SeqWidth(seq,colWidth):
	rslt = ""
	for i in range(0,len(seq),colWidth):
		rslt += "{}\n".format(seq[i:i+colWidth])
	return rslt
		
def split_fastq(fastqTXT):

	rslt = ""
	fwd = ""
        rev = ""

	arr_fq = fastqTXT.split('\n')
	arr_fq.remove('')
	for i in range(0,len(arr_fq),4):
		h1 = arr_fq[i]
		seq = arr_fq[i+1]
		h2 = arr_fq[i+2]
		q = arr_fq[i+3]
		
		m = len(seq)/2
		
		if h1.rfind(' ')>0:
			h1_bis = h1[:h1.rfind(' ')]+"_2"+h1[h1.rfind(' '):]
		else:
			h1_bis = h1+"_2"
		
		if len(h2)>1:
			if h2.rfind(' ')>0:
				h2_bis = h2[:h2.rfind(' ')]+"_2"+h2[h2.rfind(' '):]
			else:
				h2_bis = h2+"_2"
		else:
			h2_bis = h2
			
		fwdLine = "{}\n{}\n{}\n{}\n".format(h1,seq[:m],h2,q[:m])
		revLine = "{}\n{}\n{}\n{}\n".format(h1_bis,seq[m:],h2_bis,q[m:])
                fwd += fwdLine
                rev += revLine
                rslt += fwdLine + revLine
		
	return rslt,fwd,rev
	
#inputfastq: fastq query STRING
#searchfastqArray: fastq where search reads (array with multiple files)
def recup_fastq(inputfastq,searchfastqArray,paired):
	
	rslt = ""
	inputDict = get_fastq_seqdict(inputfastq)
	
	if paired and len(searchfastqArray)%2==0:
		for i in range(0,len(searchfastqArray),2):
			fq1 = openFastq(searchfastqArray[i])
			fq2 = openFastq(searchfastqArray[i+1])
			
			while True:
				h1_1 = fq1.readline().strip()
				seq_1 = fq1.readline().strip()
				h2_1 = fq1.readline().strip()
				q_1 = fq1.readline().strip()
				
				h1_2 = fq2.readline().strip()
				seq_2 = fq2.readline().strip()
				h2_2 = fq2.readline().strip()
				q_2 = fq2.readline().strip()
				
				h1_1head = h1_1[:h1_1.rfind(' ')]
				h1_2head = h1_2[:h1_2.rfind(' ')]
				
				if not h1_1: break
				
				if h1_1head!=h1_2head:
					addTextToLogFile("ERROR PAIRED FILES {} -- {}".format(searchfastqArray[i],searchfastqArray[i+1]))
					break
					
				if (seq_1 in inputDict) or (seq_2 in inputDict):
					rslt += "{}\n{}\n{}\n{}\n".format(h1_1,seq_1,h2_1,q_1)
					rslt += "{}\n{}\n{}\n{}\n".format(h1_2,seq_2,h2_2,q_2)
			
			fq1.close()
			fq2.close()
	else:
	
		for fastq in searchfastqArray:
		
			fq = openFastq(fastq)
		
			while True:
				h1 = fq.readline().strip()
				seq = fq.readline().strip()
				h2 = fq.readline().strip()
				q = fq.readline().strip()
			
				if not h1: break
			
				m = len(seq)/2
				if (seq[:m] in inputDict) or (seq[m:] in inputDict):
					rslt += "{}\n{}\n{}\n{}\n".format(h1,seq,h2,q)
			fq.close()
		
	return rslt

def openFastq(fname):

	if is_gzip(fname):
	
		fq = gzip.open(fname, "r")
	else:
		fq = open(fname,"r")
		
	return fq
	
def get_fastq_seqdict(fastq):
	rslt = {}
	c=1
	arr_fastq = fastq.split('\n')
	arr_fastq.remove('')
	for val in arr_fastq:
		if c%4==2:
			rslt[val.strip()]=1
		c+=1
		
	return rslt
	
def check_FileExtensions(filename, extList):
	fileExt = filename[filename.rfind('.')+1:]
	
	if fileExt.upper() in [ext.upper() for ext in extList]:
		return True
	
	return False
	
def load_readsList(readsList):
	rslt = {}
	
	rl = open(readsList,'r')
	
	rslt['reads_dir'] = rl.readline().strip()
	
	rslt['reads_files'] = {}
	
	i=0
	for read in rl:
		read_info = read.strip().split(',')
		rslt['reads_files']["reads_{}".format(str(i))] = read_info
		i+=1
	
	rl.close()
	
	return rslt

def find_Magic_Bytes(filename,magic_bytes):
	
	with open(filename) as infile:
		file_start = infile.read(len(magic_bytes))
		
	if file_start.startswith(magic_bytes):
		return True
		
	return False

def is_gzip(filename):

	magic_bytes = "\x1f\x8b\x08"
	
	return find_Magic_Bytes(filename,magic_bytes)

def TitleFrame(title):
	cnt = "* - {} - *".format(title)
	updn = "*"*len(cnt)
	
	return "\n{}\n{}\n{}".format(updn,cnt,updn)

def exitProg(msg=''):
	if msg!='':
		addTextToLogFile("PROGRAM EXIT\t{}".format(msg))
	sys.exit(msg)

def dirFile_exist(path,creatdir = 0):

	if not os.path.exists(path):
		if (creatdir == 0):
			sys.exit("ERROR: {} not found".format(path))
		else:
			dirmsg = "-- create {} directory".format(path)
			print dirmsg
			os.makedirs(path)

def addTextToLogFile(logtxt,printScreen=True,printTime=True):

        if printScreen:
                print logtxt

	if printTime:
		ttime = "{}\t".format(datetime.now().strftime("[%d/%m/%Y %H:%M:%S]"))
	else:
		ttime=""
	
	logFolder = config['logFolder']

	dirFile_exist(logFolder,1)

	logFile = os.path.join(logFolder,config['log_file'])
	lf = open(logFile,'a+')
	if logtxt.count('\n')>0:
		for l in logtxt.split('\n'):
			lf.write("{tm}{log}\n".format(tm=" "*len(ttime),log=l))
	else:
		lf.write("{tm}{log}\n".format(tm=ttime,log=logtxt))
	lf.close()

	
if	__name__ == '__main__':

	ArgsVal = sys.argv[1:]
	run(ArgsVal)
