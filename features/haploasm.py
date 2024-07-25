#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Mathieu Genete

@NGS-feature
@NGS-description: Assembly pipeline from genotyping results
"""
#=== imports libraries ===
import sys
import os
import argparse
import threading
import time
from datetime import datetime
from ete3 import Tree
from Bio import Phylo
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

#home libraries
import utils
import load_dump as ld

utils.set_paramFileName('NGSgenotyp2_haploasm')

__version__ = "2.0.1"
AppName = "NGSgenotyp haploAsm"

args = None
config = None

def main(ArgsVal,main_configs):

    global args
    global config

    asm_config_file = os.path.join(main_configs['configs_folder'],'assembly_config.yaml')
    config = ld.yaml_loadDic(asm_config_file)

    assemblers_list=config['ASMtools'].keys()
    
    description = """ NGSgenotyp Haplotyp Assembly """

    parser = argparse.ArgumentParser(prog="haploAsm",description=description)
	
    parser.add_argument("-V", "--verbose", help="full verbose", action="store_const", const="-V", default="")
    parser.add_argument("-v", "--tinyverbose", help="verbose",action="store_const", const="-v", default="")
    parser.add_argument("-f", "--force", action="store_const", const="-f", default="")
    parser.add_argument("-c", "--checks", help=argparse.SUPPRESS, action="store_const", const="-c", default="")
    parser.add_argument("-t", "--test", help=argparse.SUPPRESS, action="store_const", const="-t", default="")
    parser.add_argument("-p","--paired",help="orignals reads files with paired-end reads forward and reverse files distincts",action="store_const", const=True, default=False)
    parser.add_argument("-y","--includeParalogsyass",help="include paralogs for yass analyse",action="store_const", const=True, default=False)
    parser.add_argument("-x","--includeParalogstree",help="include paralogs for phylogeny",action="store_const", const=True, default=False)

    parser.add_argument("-T","--threaded",help="define threads number used for phylogeny. By default 8 cpus",type=int,default=8)	
    parser.add_argument("-st","--spadesThreads",help="define threads number used for spades. By default 50%% of cpu numbers",type=int)	
    parser.add_argument("-l","--indivlist", help="list of haplotype name for assembly")
    parser.add_argument("-m","--contigminlen", help="keep contigs up from minimum lenght (default=500)",type=int,default=500)
    parser.add_argument("-M","--contigmaxlen", help="keep contigs down to maximum lenght",type=int)
    parser.add_argument("-g","--genotypfolder",help="genotyping folder", required = True)
    parser.add_argument("-s","--asmsuffix",help="assembly folder suffix", required = True)
    parser.add_argument("-a","--assembler",help="Program used to assemble alleles ({}) defaut: spades".format(",".join(assemblers_list)), default="spades")
	
    args = parser.parse_args(ArgsVal)

    startTime=time.time()

    genotypfolder_sublen=len(str(args.genotypfolder).strip("/").split("/"))

    if genotypfolder_sublen>1:
        genotypfolder_sub="/".join(str(args.genotypfolder).split("/")[:-1])
    else:
        genotypfolder_sub=args.genotypfolder

    #Test if asm_config.yaml exists in the genotyping configs folder
    assembly_config_file=os.path.join(args.genotypfolder,config['assembly_config'])
    if not os.path.exists(assembly_config_file):
        print("{} file not found".format(assembly_config_file))
        sys.exit(1)

    out_asembly=os.path.join(args.genotypfolder,"Assembly_{}".format(args.asmsuffix))
    
    #update config variable
    config['outfolder']=out_asembly
    config['log_folder']=os.path.join(config['outfolder'],'logs')

    config['reffolder']=os.path.join(config['outfolder'],"refs")
    dir_create(config['outfolder'])
    dir_create(config['reffolder'])
    
    #haploasm log informations
    cmdline = " ".join(sys.argv)
    currentDateTime = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    config['file_time']=datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
    config['log_file'] = "{appn}_Log_{dt}.txt".format(appn=AppName,dt=config['file_time'])
    startText = "{appn} {appv} -- {cdate}".format(appn=AppName,appv=__version__,cdate=currentDateTime)
    starDateline = "*"*len(startText)

    addTextToLogFile("{}\n{}\n{}".format(starDateline,startText,starDateline),False)
    addTextToLogFile("user command: {}".format(cmdline))


    assembly_config = ld.yaml_loadDic(assembly_config_file)

    #=======================
    #==== MAIN WORKFLOW ====
    #============================>>>
    reads_data = get_reads_data(assembly_config,genotypfolder_sub)
    pkfilename,reference_dic = get_reference(assembly_config,genotypfolder_sub)
    ref_fasta=os.path.join(config['reffolder'],"{}.fasta".format(pkfilename))
    write_reference_fasta(ref_fasta,reference_dic,add_paralogs=args.includeParalogsyass)
    asm_contigs=assembly(reads_data,args.assembler,config,main_configs)

    asm_contigs=clean_contigs_file(asm_contigs)

    contigs_infos=parse_contigs_infos(asm_contigs,args.assembler)
    yass_results=do_yass(asm_contigs,ref_fasta,main_configs)
    raw_yass_file=parse_yass_results(yass_results,contigs_infos,reference_dic,args.assembler)
    oriented_contigs=get_oriented_contigs(raw_yass_file,contigs_infos,asm_contigs)
    all_and_ref,sample_phylo=get_phylo_fasta(oriented_contigs,reference_dic,asm_contigs,add_paralogs=args.includeParalogstree)

    phylip_merged,phylip_samples=do_muscle_alignments(all_and_ref,sample_phylo,main_configs)

    phyml_merged,phyml_samples=do_phyml(main_configs,phylip_merged,phylip_samples)

    trees_to_pdf(phyml_merged,phyml_samples)
    #END of main analysis
    #produce documentation

    list_contigs_nearest_node(phyml_samples)

    #delete utils param file
    utils.delete_paramFileName()

    totalTime = time.time() - startTime
    endMessage = "-- ENDED IN {} --".format(time.strftime('%Hh%Mm%Ss', time.gmtime(totalTime)))
    addTextToLogFile(endMessage,True)
    #<<<<===========================
    #    ==== END MAIN WORKFLOW ====
    #    ===========================        

def trees_to_pdf(phyml_merged,phyml_samples):
    #merged tree
    addTextToLogFile("create tree pdf files",True)
    merged_dirname = os.path.dirname(phyml_merged)
    merged_outpdf = os.path.join(merged_dirname,"All_samples_phylogeny.pdf")
    draw_tree(phyml_merged,merged_outpdf)

    #samples trees
    for sample,treefile in phyml_samples.items():
        tree_dirname = os.path.dirname(treefile)
        out_pdf = os.path.join(tree_dirname,"{}_phylogeny.pdf".format(sample))
        draw_tree(treefile,out_pdf)

def draw_tree(treefile,outpdf):
    tree = Phylo.read(treefile,'newick')
    Phylo.draw(tree)
    plt.axis('off')
    fig=plt.gcf()
    fig.set_size_inches(tree.total_branch_length()*3, len(tree.get_terminals())*0.3)
    plt.savefig(outpdf)
    
def list_contigs_nearest_node(phyml_samples):
    out_nearest={}
    for sample,phyml_tree in phyml_samples.items():
        tree = Tree(phyml_tree)
        out_nearest[sample]={}
        for leaf in tree:
            if 'contig' in leaf.name:
                prev_nod=leaf
                end_while=False
                out_nearest[sample][leaf.name]=[]
                while True:
                    prev_nod=prev_nod.up
                    for child_leaf in prev_nod:
                        if 'contig' not in child_leaf.name:
                            out_nearest[sample][leaf.name].append((child_leaf.name,tree.get_distance(leaf,child_leaf)))
                            end_while=True

                    if prev_nod.is_root() or end_while:
                        break

    out_txt="#Contigs closest Alleles\n"
    out_txt+="#\tContig\tAllele\tdistance\n"
    for sample,contigs in out_nearest.items():
        out_txt+="sample: {}\n".format(sample)
        for contig,alleles in contigs.items():
            if len(alleles)>0:
                closest_allele=sorted(alleles, key = lambda x:x[1])[0]
                out_txt+="\t{}\t{}\t{}\n".format(contig,closest_allele[0],closest_allele[1])

    out_closest_alleles=os.path.join(config['outfolder'],"contigs_closest_alleles.txt")
    stdout_print("Write contigs closests alleles file: {}".format(out_closest_alleles),printInLog=True)
    with open(out_closest_alleles,'w') as outc:
        outc.write(out_txt)
    
def return_leaf_dist(node,outlist=None):
    if outlist is None:
        outlist=[]
    for child_node in prev_nod.get_children():
        if child_node.is_leaf:
            outlist

def do_phyml(main_configs,phylip_merged,phylip_samples,phyml_cpu=0):
    out_phyml={}
    phyml_configs=main_configs['tools_avail']['phyml']
    if phyml_configs['uselocal']:
        phyml_path=phyml_configs['folder']
    else:
        phyml_path=""

    phyml_jblst = utils.Jobslist("do phyml for all")
    
    #Prepare phyml for each sample
    for sample,files in phylip_samples.items():
        PhylipFilename=files['phylip']
        phyml_tree_file=PhylipFilename+"_phyml_tree.txt"
        if os.path.exists(PhylipFilename):
            out_phyml[sample]=phyml_tree_file
            cmd = "{path}phyml --quiet --no_memory_check -i {phy} -d nt > {logphyml} 2>&1".format(path=phyml_path,phy=PhylipFilename,logphyml=os.path.join(os.path.dirname(PhylipFilename),"phyml_log"))
            phyml_jblst.add_a_job(cmd,"phyml {}".format(sample),phyml_tree_file)
        else:
            stdout_print("WARNING: file '{}' NOT EXISTS".format(PhylipFilename),printInLog=True)

    #Prepare phyml for the merged samples file
    PhylipFilename_merged=phylip_merged['phylip']
    phyml_tree_file_merged=PhylipFilename_merged+"_phyml_tree.txt"
    cmd = "{path}phyml --quiet --no_memory_check -i {phy} -d nt > {logphyml} 2>&1".format(path=phyml_path,phy=PhylipFilename_merged,logphyml=os.path.join(os.path.dirname(PhylipFilename_merged),"phyml_log"))
    phyml_jblst.add_a_job(cmd,"phyml {}".format(sample),phyml_tree_file_merged)
        
    if args.tinyverbose:
        t1 = threading.Thread(target=progress_bar,args=['Phyml progress',[v[2] for v in phyml_jblst.get_joblist()]])
        t1.start()
        
    stdout_print("Launch phyml",printInLog=True)
    utils.trun(args, phyml_jblst)

    return phyml_tree_file_merged,out_phyml

def do_muscle_alignments(all_and_ref,sample_phylo,main_configs):
    out_muscle={}
    muscle_configs=main_configs['tools_avail']['muscle']
    if muscle_configs['uselocal']:
        muscle_path=muscle_configs['folder']
    else:
        muscle_path=""

    muscle_jblst = utils.Jobslist("do muscle for all")

    #Prepare Muscle for each sample
    for sample, phylo_fasta in sample_phylo.items():
        bname_phylo_fasta=os.path.basename(phylo_fasta)
        MuscleFilename=os.path.join(os.path.dirname(phylo_fasta),"MUSCLE_"+bname_phylo_fasta)
        phylipFilename=os.path.join(os.path.dirname(phylo_fasta),"PHYLIP_"+bname_phylo_fasta[:bname_phylo_fasta.rfind('.')]+".phylip")
        if os.path.exists(phylo_fasta):
            out_muscle[sample]={'fasta':MuscleFilename,'phylip':phylipFilename}
            cmd = "{path}muscle -in {infasta} -out {outMsl} > {logmsl} 2>&1".format(path=muscle_path,infasta=phylo_fasta,outMsl=MuscleFilename,logmsl=os.path.join(os.path.dirname(phylo_fasta),"muscle_log.txt"))
            muscle_jblst.add_a_job(cmd,"muscle {}".format(sample),MuscleFilename)
        else:
            stdout_print("WARNING: file '{}' NOT EXISTS".format(phylo_fasta),printInLog=True)

    #Prepare Muscle for the merged samples file
    bname_phylo_merged=os.path.basename(all_and_ref)
    MuscleMergedFilename=os.path.join(os.path.dirname(all_and_ref),"MUSCLE_"+bname_phylo_merged)
    phylipMergedFilename=os.path.join(os.path.dirname(all_and_ref),"PHYLIP_"+bname_phylo_merged[:bname_phylo_merged.rfind('.')]+".phylip")
    cmd = "{path}muscle -in {infasta} -out {outMsl} > {logmsl} 2>&1".format(path=muscle_path,infasta=all_and_ref,outMsl=MuscleMergedFilename,logmsl=os.path.join(os.path.dirname(all_and_ref),"muscle_log.txt"))
    muscle_jblst.add_a_job(cmd,"muscle merged",MuscleMergedFilename)
        
    if args.tinyverbose:
        t1 = threading.Thread(target=progress_bar,args=['Muscle progress',[v[2] for v in muscle_jblst.get_joblist()]])
        t1.start()
        
    stdout_print("Launch muscle",printInLog=True)
    utils.trun(args, muscle_jblst)

    for sample,datas in out_muscle.items():
        fastaToPhylip(datas['fasta'],datas['phylip'])
        
    fastaToPhylip(MuscleMergedFilename,phylipMergedFilename)

    return {'fasta':MuscleMergedFilename,'phylip':phylipMergedFilename},out_muscle

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

        with open(outPhy,"w") as of:
            of.write("{} {}\n".format(N,seqLen))
            of.write(phyOut)
            
        return True
    
def get_phylo_fasta(oriented_contigs,reference,asm_contigs,add_paralogs=True):
    reffasta=""
    outall_and_ref=os.path.join(config['outfolder'],"phylo_all_samples.fa")
    for allele,datas in reference.items():
            if (not datas['Paralog'] or add_paralogs) and not datas['combined']:
                allele+="|Group={}|haploID={}".format(datas['grpID'],datas['haploID'])
                if datas['Paralog']: allele+="_Paralog"
                reffasta+=">{}\n{}\n".format(allele,datas['seq'])

    all_and_ref=""
    out_sample_phylo={}
    for sample,datas in oriented_contigs.items():
        outfasta=""
        if datas['nb_contigs']>0:
            fasta = SeqIO.parse(datas['c_phylo'],'fasta')
            for rec in fasta:
                outfasta+=">{}\n{}\n".format(sample+"_"+rec.description,rec.seq)
            outfasta_all=os.path.join(asm_contigs[sample]['asmfolder'],"phylo_all_contigs.fa")
            out_sample_phylo[sample]=outfasta_all
            all_and_ref+=outfasta
            outfasta+=reffasta
            with open(outfasta_all,"w") as outall:
                outall.write(outfasta)
            
    all_and_ref+=reffasta
    with open(outall_and_ref,"w") as outall_r:
        outall_r.write(all_and_ref)
        
    return outall_and_ref,out_sample_phylo

def get_oriented_contigs(raw_yass_file,contigs_infos,asm_contigs):
    yass_datas=ld.pickle_loadDic(raw_yass_file)
    out_yass_contigs={}
    for sample,all_datas in yass_datas.items():
        contigs= {rec.description:rec.seq for rec in SeqIO.parse(asm_contigs[sample]['cl_contigs'],'fasta')}
        stdout_print("\n====> Extract contigs from yass -- sample: {} -- total number of contigs: {}".format(sample,len(contigs)),printInLog=True)
        rc_count=0
        unique_contigs_list=set([(str(datas['Query_id']),datas['strand'],datas['contig_depth']) for datas in all_datas])
        phylo_cut={}
        for datas in all_datas:
            q_id=str(datas['Query_id'])
            if q_id not in phylo_cut.keys():
                phylo_cut[q_id]={'start':-1,'end':0}
            q_start=datas['q_start']-1
            q_end=datas['q_end']

            if q_start<phylo_cut[q_id]['start'] or phylo_cut[q_id]['start']==-1:
                phylo_cut[q_id]['start']=q_start

            if q_end>phylo_cut[q_id]['end']:
                phylo_cut[q_id]['end']=q_end
        extracted_contigs={}
        for contig_id,strand,contig_depth in unique_contigs_list:
            contig_description=contigs_infos[sample][contig_id]['description']
            contig_seq=str(contigs[contig_description])
            phylo_seq=contig_seq[phylo_cut[contig_id]['start']:phylo_cut[contig_id]['end']]
            if strand=='-':
                rc_count+=1
                stdout_print("\t{} reverse complement: {}".format(rc_count,contig_description),printInLog=True)
                contig_seq=reverse_complement(contig_seq)
                phylo_seq=reverse_complement(phylo_seq)
            extracted_contigs_id="contig={}|l={}|d={}|s={}".format(contig_id,len(contig_seq),contig_depth,strand)
            if strand=="-":
                extracted_contigs_id+="_RC"
            extracted_contigs[extracted_contigs_id]={'c_seq':contig_seq,'phylo_seq':phylo_seq}
        stdout_print("{} contigs extracted <====".format(len(extracted_contigs)),printInLog=True)
        outfasta=os.path.join(asm_contigs[sample]['asmfolder'],"yass_oriented_contigs.fa")
        outfasta_phylo=os.path.join(asm_contigs[sample]['asmfolder'],"phylo_truncated_contigs.fa")
        out_yass_contigs[sample]={'contigs':outfasta,'c_phylo':outfasta_phylo,'nb_contigs':len(extracted_contigs.keys())}
        with open(outfasta_phylo,'w') as outphylo:
            with open(outfasta,'w') as outf:
                for c_id,seq in extracted_contigs.items():
                    outf.write(">{}\n{}\n".format(c_id,seq['c_seq']))
                    outphylo.write(">{}\n{}\n".format(c_id,seq['phylo_seq']))
    return out_yass_contigs
    
def parse_yass_results(yass_results,contigs_infos,reference_dic,assembler,filter_paralogs=True):
    out_yass={}
    outfile_yass=os.path.join(config['outfolder'],"raw_yass_datas.pk")
    strand_val={1:"+",-1:"-"}
    for sample,yass_stdout in yass_results.items():
        yass_rslt=parse_yass_stdout(yass_stdout,assembler)
        out_rslt=[]
        for datas in yass_rslt:
            ref_len=reference_dic[datas['Subject_id']]['len_seq']
            datas['ref_length']=ref_len
            ident_cov=datas['percent_identity']*(datas['alignment_length']/ref_len)
            datas['coverage']=(datas['alignment_length']/ref_len)
            datas['ident_cov']=ident_cov
            datas['contig_depth']=contigs_infos[sample][str(datas['Query_id'])]['cov']
            datas['strand']=strand_val[int((datas['s_end']-datas['s_start'])/abs(datas['s_end']-datas['s_start']))]
                
            #Apply Thresholds
            if all([eval(str(datas[param])+thrld) for param,thrld in config['YASSstep']['Thresholds'].items()]):
                out_rslt.append(datas)
        out_yass[sample]=out_rslt

    ld.pickle_dumpDic(outfile_yass,out_yass)
    return outfile_yass
            

def parse_yass_stdout(yass_stdout,assembler):
    #query => contigs
    #subject => ref allele
    out_rslt=[]
    fields=['Query_id','Subject_id', 'percent_identity', 'alignment_length', 'mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'e-value', 'bit_score']
    contig_nbr_regex=config['ASMtools'][assembler]['regex_contig_nbr']
    with open(yass_stdout,'r') as yassfile:
        for line in yassfile:
            if not line.startswith("#"):
                tmp=line.strip().split('\t')
                tmp[0]=re.search(contig_nbr_regex,tmp[0]).group(0)
                out_rslt.append({fields[i]:string_to_number(val) for i,val in enumerate(tmp)})
    return out_rslt

def string_to_number(text):
    if bool(re.search("^[+-]?([0-9]*[.])?[0-9]+",text)):
        if "." not in text and "e" not in text:
            return int(text)
        return float(text)
    return text

def parse_contigs_infos(asm_contigs,assembler):
    out_infos={}
    cov_regex=config['ASMtools'][assembler]['regex_cov']
    contig_nbr_regex=config['ASMtools'][assembler]['regex_contig_nbr']
    for sample,datas in asm_contigs.items():
        out_infos[sample]={}
        contigs_file=datas['contigs']
        fasta=SeqIO.parse(contigs_file,'fasta')
        for rec in fasta:
            seq_id=str(rec.description)
            COV_rslt=re.search(cov_regex,seq_id)
            contig_nbr=re.search(contig_nbr_regex,seq_id).group(0)
            if COV_rslt is not None:
                cov_val=float(COV_rslt.group(0))
            else:
                cov_val=0.0
                
            out_infos[sample][contig_nbr]={'description':seq_id,'length':len(rec.seq),'cov':cov_val}
    return out_infos
        
    
def clean_contigs_file(asm_contigs):
    out_asm_contigs={}
    for sample,datas in asm_contigs.items():
        contigs_file=datas['contigs']
        cleaned_contigs=os.path.join(datas['asmfolder'],"cleaned_contigs.fa")
        datas['cl_contigs']=cleaned_contigs
        if os.path.exists(contigs_file):
            out_asm_contigs[sample]=datas
            if not os.path.exists(cleaned_contigs):
                fasta = SeqIO.parse(contigs_file,'fasta')
                with open(cleaned_contigs,'w') as outcontigs:
                    for rec in fasta:
                        if len(rec.seq)>=args.contigminlen and (not args.contigmaxlen or len(rec.seq)<=contigmaxlen):
                            outcontigs.write(rec.format('fasta'))
    return out_asm_contigs
                
def do_yass(asm_contigs,ref_fasta,main_configs):
    yass_out={}
    yass_jblst = utils.Jobslist("do yass for all")
    
    yass_configs=main_configs['tools_avail']['yass']
    if yass_configs['uselocal']:
        yass_folder=yass_configs['folder']
    else:
        yass_folder=""
        
    for sample,datas in asm_contigs.items():
        indiv_folder=datas['asmfolder']
        yass_stdout = os.path.join(indiv_folder,"yass_stdout")
        yass_stderr = os.path.join(indiv_folder,"yass_stderr")
        yass_out[sample]=yass_stdout
        cmd="{path}yass -d 2 {params}{contigFile} {refFile} > {ystdout} 2> {ystderr}".format(ystderr=yass_stderr,ystdout=yass_stdout,path=yass_folder,params=config['YASSstep']['yass_params'],contigFile=datas['cl_contigs'],refFile=ref_fasta)
        yass_jblst.add_a_job(cmd,"yass {}".format(sample),yass_stdout)
        
    if args.tinyverbose:
        t1 = threading.Thread(target=progress_bar,args=['Yass progress',[v[2] for v in yass_jblst.get_joblist()]])
        t1.start()
        
    stdout_print("Launch yass",printInLog=True)
    utils.trun(args, yass_jblst)

    if args.tinyverbose:
        wait_progressbar_finished(t1)
        
    return yass_out

def assembly(reads_data,assembler,config,main_configs):
    #usedKmerSize = "21,41,81"
    #cmd = "{path}spades.py -t {thrds} --careful -k {kmerS} -o {outf} {fastq} > {logf}".format(path=spades_folder,outf=spadesOutFolder,kmerS=usedKmerSize,fastq=inputSpadesFQ,thrds=spThrd,logf=spadesLogFile)
    stdout_print("Launch Assembly using '{}'".format(assembler),printInLog=True)
    out_asm={}
    cur_dir=os.getcwd()
    asm_config=config['ASMtools'][assembler]
    tools_folder=main_configs['tools_folder']
    asm_path=os.path.join(tools_folder,asm_config['path'])
    if os.path.exists(asm_path):
        c=0
        for indiv,datas in reads_data.items():
            out_assembly=os.path.join(config['outfolder'],indiv,"Assembly_{}".format(assembler))
            #log_assembly=os.path.join(out_assembly,"log_assembly")
            log_assembly=os.path.join("log_assembly")
            dir_create(out_assembly)
            if args.paired and len(datas['reads'])%2==0:
                fwd_reads=asm_config['fastq_sep'].join(["{}{}".format(asm_config['options']['fwd'],os.path.join(cur_dir,fq)) for c,fq in enumerate(datas['reads']) if c%2==0])
                rev_reads=asm_config['fastq_sep'].join(["{}{}".format(asm_config['options']['rev'],os.path.join(cur_dir,fq)) for c,fq in enumerate(datas['reads']) if c%2==1])
                fastq_str=fwd_reads+asm_config['fastq_sep']+rev_reads
            else:
                fastq_str=asm_config['fastq_sep'].join(["{}{}".format(asm_config['options']['single'],os.path.join(cur_dir,fq)) for fq in datas['reads']])
                
            options = asm_config['cmdline'].format(thrds=args.threaded,kmerS=asm_config['kmers'],outf=out_assembly,fastq=fastq_str,logf=log_assembly)
            cmd = asm_path + " " + options
            #print(cmd)
            
            out_contigs=os.path.join(out_assembly.format(assembler),asm_config['out_contigs'])
            
            #Create contig simlink:
            out_contigs_local=os.path.join("Assembly_{}".format(assembler),asm_config['out_contigs'])
            out_contigs_link="contigs.fa"
            os.chdir(os.path.join(cur_dir,config['outfolder'],indiv))
            if os.path.islink(out_contigs_link):
                os.remove(out_contigs_link)
            os.symlink(out_contigs_local,out_contigs_link)
            os.chdir(cur_dir)
            out_contigs_link_abs=os.path.join(config['outfolder'],indiv,out_contigs_link)
            
            if not os.path.exists(out_contigs):
                os.chdir(os.path.join(cur_dir,out_assembly))
                os.system(cmd)
                os.chdir(cur_dir)

            asm_OK=True
            if not os.path.exists(out_contigs_link_abs):
                asm_OK=False

            out_asm[indiv]={'contigs':out_contigs_link_abs,'asmfolder':os.path.join(config['outfolder'],indiv),'asm_OK':asm_OK}
            
            c+=1
            prog_pct=100*float(c)/float(len(reads_data.keys()))
            progress_bar_pct("Assembly progress",prog_pct)
    else:
        exitProg("Tool not found: {}".format(asm_path))
    print("")
    return out_asm

def write_reference_fasta(out_fasta,reference,add_paralogs=False):
    #'seq', 'len_seq', 'combined', 'Paralog', 'specie', 'grpRef'
    #'allelePart', 'grpRefPart', 'grpIDPart', 'haploIDPart', 'grpID', 'haploID'
    with open(out_fasta,'w') as f_alleles:
        for allele,datas in reference.items():
            if (not datas['Paralog'] or add_paralogs) and not datas['combined']:
                f_alleles.write(">{}\n{}\n".format(allele,datas['seq']))
                    
def get_reference(assembly_config,genotypfolder_sub):
    return os.path.basename(assembly_config['ref']),ld.pickle_loadDic(os.path.join(genotypfolder_sub,assembly_config['ref']))

def get_reads_data(assembly_config,genotypfolder_sub):
    reads_data = ld.yaml_loadDic(os.path.join(genotypfolder_sub,assembly_config['readsconfig']))
    if args.indivlist:
        indiv_to_assemble=[]
        with open(args.indivlist,'r') as indivfile:
            for line in indivfile:
                if line!="":
                    indiv_to_assemble.append(line.strip())
        ret_reads_datas={indiv:datas for indiv,datas in reads_data.items() if indiv in indiv_to_assemble}

        if len(ret_reads_datas.keys())>0:
            return ret_reads_datas
    return reads_data

def reverse_complement(seq):
    rev_seq=seq[::-1].upper()
    bases={'A':'T','T':'A','C':'G','G':'C','N':'N'}
    rev_comp=""
    for b in rev_seq:
        if b in bases.keys():
            rev_comp+=bases[b]
        else:
            rev_comp+="N"
    return rev_comp

def progress_bar(title,file_list):
    finished=0
    nbr=len(file_list)
    while finished<nbr:
        time.sleep(2)
        finished=sum(1 for f in file_list if os.path.exists(f) and os.path.getsize(f)>0 )
        percent=round(100*float(finished)/float(nbr),1)

        progstr="="*(int(percent/2)-1)+">"+" "*(50-int(percent/2))
        if percent>=100:
            progstr="="*int(percent/2)

        printProgress("{}: [{}] {}%".format(title,progstr,round(percent,1)))
    print("")

def progress_bar_pct(title,percent):
    progstr="="*(int(percent/2)-1)+">"+" "*(50-int(percent/2))
    if percent>=100:
        progstr="="*int(percent/2)
    printProgress("{}: [{}] {}%".format(title,progstr,round(percent,1)))
        
def wait_progressbar_finished(pbThrd):
    while(pbThrd.is_alive()):
        time.sleep(2)

def printProgress(txt):
    sys.stdout.write("\r"+str(txt))
    sys.stdout.flush()

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
