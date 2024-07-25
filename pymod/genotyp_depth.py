#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import load_dump as ld
import numpy as np

def normalize_depth(outstats,pk_refdatabase):
    refinfo=ld.pickle_loadDic(pk_refdatabase)

    paralogs_ref=set([ref for ref in refinfo.keys() if refinfo[ref]['Paralog']])
    
    for sample,sample_alleles in outstats.items():
        PosAlleles_list=list([PosAllele for PosAllele in sample_alleles.keys() if sample_alleles[PosAllele]['IsPositive']])
        PosParalogs_list=set(PosAlleles_list).intersection(paralogs_ref)
        mean_paralog_covered_depth=float(np.median([sample_alleles[paralog]['coverage']['depth_covered_mean'] for paralog in PosParalogs_list]))
        print(sample,"==>",mean_paralog_covered_depth)
        for allele in sample_alleles.keys():
            depth_covered_mean=outstats[sample][allele]['coverage']['depth_covered_mean']
            NORM_depth_covered_mean = depth_covered_mean/mean_paralog_covered_depth
            outstats[sample][allele]['coverage'].update({'NORM_depth_covered_mean':NORM_depth_covered_mean})
