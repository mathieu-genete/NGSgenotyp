#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import math
import numpy as np
import os
import gzip
import pickle
from scipy.integrate import quad
import load_dump as ld

def determine_genotypes(outstats,pk_refdatabase,lbda,mu,sigma,genotyp_alleleProb_THRLD):

    refinfo=ld.pickle_loadDic(pk_refdatabase)
    
    for sample,sample_alleles in outstats.items():
        for allele,datas in sample_alleles.items():
            prob_rslt=allele_prob(datas['stats']['error rate'],datas['coverage']['covered_pct'],mu,sigma,lbda)
            #print(sample,refinfo[allele],prob_rslt)

            IsPositive_val=(prob_rslt['Prob_allele']>=genotyp_alleleProb_THRLD)
            
            outstats[sample][allele].update({'IsPositive':IsPositive_val,'gscore':prob_rslt['Prob_allele'],'Prob_err':prob_rslt['Prob_err'],'Prob_cov':prob_rslt['Prob_cov']})

            #Is Allele positive ? => genotypescore >= genotyp_alleleProb_THRLD
            if IsPositive_val:
                print(sample,allele,datas['stats']['error rate'],datas['coverage']['covered_pct'],refinfo[allele]['grpRef'],prob_rslt)

def get_combined_allels(pk_refdatabase):
    DByaml=ld.pickle_loadDic(pk_refdatabase)
    AllelesCombined={}
    for allele,adatas in DByaml.items():
        if adatas['allelePart'] and adatas['haploIDPart'] and adatas['grpRefPart']:
            for subAllele,subGrp,subHaplo in zip(adatas['allelePart'],adatas['grpIDPart'],adatas['haploIDPart']):
                if subAllele not in AllelesCombined:
                    AllelesCombined[subAllele]={'allelesPart':[],'haploID':0,'Grp':0}
                AllelesCombined[subAllele]['allelesPart'].append(allele)
                AllelesCombined[subAllele]['haploID']=subHaplo
                AllelesCombined[subAllele]['GrpID']=subGrp
    return AllelesCombined

def combine_alleles_stats(outstats,pk_refdatabase,removeAllelesPart=True):

    AllelesCombined=get_combined_allels(pk_refdatabase)

    print("combined:",AllelesCombined)
    
    for sample,sample_alleles in outstats.items():
        MainAllele_stats={}
        for MainAllele,madatas in AllelesCombined.items():
            MainAllele_stats[MainAllele]={}
            for allele in madatas['allelesPart']:
                if allele in sample_alleles.keys():
                    MainAllele_stats[MainAllele][allele]=outstats[sample][allele]
                    if removeAllelesPart:
                        del outstats[sample][allele]
            if len(MainAllele_stats[MainAllele].keys())==0:
                del MainAllele_stats[MainAllele]
        merge_result=merge_alleles_stats(MainAllele_stats)
        #print(sample,[(v,len(d))for v,d in merge_result.items()])
        outstats[sample].update(merge_result)
        
    return outstats
                
def merge_alleles_stats(MainAllele_stats):
    merged_stats={}
    for MainAllele,allelesPart in MainAllele_stats.items():
        stats_val={'stats': {'raw sequences': 0.0, 'bases mapped': 0, 'mismatches': 0, 'reads mapped': 0, 'average length': 0.0, 'error rate': None}, 'reduced_stats': {'bases mapped': 0, 'mismatches': 0,'reads mapped': 0, 'average length': 0, 'error rate': None}, 'MQ_prob': 0,'readsNb': 0, 'ratio_readsNoInDel': 0.0, 'readsNoInDelNbr': 0, 'mismatch_mean': 0, 'distrib_XM_values': None, 'ref_length': 0, 'coverage': {'depth_mean': 0.0, 'depth_covered_mean': 0.0, 'covered_pos': 0, 'covered_pct': 0.0, 'min_depth': None, 'max_depth': None}, 'reduced_coverage': {'depth_mean': 0.0, 'depth_covered_mean': 0.0, 'covered_pos': 0, 'covered_pct': 0.0, 'min_depth': None, 'max_depth': None}}
        for allele,stats in allelesPart.items():
            for param  in stats_val['stats'].keys():
                if stats_val['stats'][param]!=None:
                    stats_val['stats'][param]+=stats['stats'][param]
            
            for param  in stats_val['reduced_stats'].keys():
                if stats_val['reduced_stats'][param]!=None:
                    stats_val['reduced_stats'][param]+=stats['reduced_stats'][param]

            for param  in stats_val['coverage'].keys():
                if stats_val['coverage'][param]!=None:
                    stats_val['coverage'][param]+=stats['coverage'][param]

            for param  in stats_val['reduced_coverage'].keys():
                if stats_val['reduced_coverage'][param]!=None:
                    stats_val['reduced_coverage'][param]+=stats['reduced_coverage'][param]

            for param  in [k for k,vals in stats_val.items() if type(vals)!=dict]:
                if stats_val[param]!=None:
                    stats_val[param]+=stats[param]

        params_to_mean=['average length','ratio_readsNoInDel','mismatch_mean','depth_mean','depth_covered_mean','covered_pct','MQ_prob']

        #print("MainAllele =>",MainAllele, allelesPart.keys())
        for param  in [k for k,vals in stats_val.items() if type(vals)!=dict]:
            if stats_val[param]!=None and param in params_to_mean:
                stats_val[param]=stats_val[param]/float(len(allelesPart))

        for param  in [k for k,vals in stats_val.items() if type(vals)==dict]:
            for sub_param in stats_val[param].keys():
                if stats_val[param][sub_param]!=None and sub_param in params_to_mean:
                    stats_val[param][sub_param]=stats_val[param][sub_param]/float(len(allelesPart))
                elif sub_param=='error rate':
                    if stats_val[param]['bases mapped']>0:
                        stats_val[param][sub_param]=stats_val[param]['mismatches']/stats_val[param]['bases mapped']
                    else:
                        stats_val[param][sub_param]=None
                
        merged_stats[MainAllele]=stats_val
    return merged_stats
            

#                 Score calculation
#------------------------------------------------------
def loi_Exponentielle(x,l):
    if(x>0):
        return l*math.exp(-l*(1-x))
    else:
        return l

def loi_LogNormal(x,mu,sigma):
        if(x>0):
                p1 = 1/(x*sigma*math.sqrt(2*math.pi))
                p2 = math.exp(-((math.log(x)-mu)**2/(2*sigma**2)))
                return (p1*p2)
        else:
                return 0.0

def loi_Normal(x,mu,sigma):
        p1 = 1/(sigma*math.sqrt(2*math.pi))
        p2 = math.exp(-0.5*((x-mu)/sigma)**2)
        return p1*p2
        
def get_ErrProb(err,mu,sigma):
    res, err = quad(loi_Normal,err,1,args=(mu,sigma))
    return res

def get_CovProb(cov,l):
    res, err= quad(loi_Exponentielle,0,cov,args=(l))
    return res

def allele_prob(err,cov,mu,sigma,l):
        errProb = get_ErrProb(err,mu,sigma)
        covProb = get_CovProb(cov,l)
        
        # use errProb for error probability and covProb for coverage probability
        allProb = errProb*covProb**0.5
        
        return {'Prob_allele':allProb,'Prob_err':errProb,'Prob_cov':covProb}
