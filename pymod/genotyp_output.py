#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import xlwt
import load_dump as ld
import os

config = None

def write_xls_output(out_folder,results_folder,outstats,reads_nbr,pk_refdatabase,mismatchthrld,inconfig):

    global config
    config = inconfig
    
    refinfo=ld.pickle_loadDic(pk_refdatabase)
    
    samples_per_file=config['MaxSheetNbr']
    projectname=os.path.basename(out_folder)
    sample_list=sorted(outstats.keys())
    sample_list_byFile=[sample_list[i:i+samples_per_file] for i in range(0,len(sample_list),samples_per_file)]
    xls_file_nbr=len(sample_list_byFile)
    xls_generic_name="Genotyp_stats_{project}_{fileNbr}.xls"
    if samples_per_file>1:
        xls_file_list=[os.path.join(results_folder,xls_generic_name.format(project=projectname,fileNbr=i)) for i in range(1,xls_file_nbr+1)]
    else:
        xls_file_list=[os.path.join(results_folder,xls_generic_name.format(project=projectname,fileNbr=sampname)) for sampname in sample_list]

    xls_samplename_txt=""

    #XLS formatting
    style = xlwt.easyxf('alignment: horizontal center, vertical center;')
    hstyle = xlwt.easyxf('font: bold on; alignment: horizontal center, vertical center;')
    greenStyle = xlwt.easyxf('alignment: horizontal center, vertical center;pattern: pattern solid, fore_colour light_green;')
    orangeStyle = xlwt.easyxf('alignment: horizontal center, vertical center;pattern: pattern solid, fore_colour light_orange;')
    
    for xls_filename,samples_names in zip(xls_file_list,sample_list_byFile):
        workbook = xlwt.Workbook()
        headers=['Genotyp Score','Error rate','Mean covered Depth','Normalized covered Depth','Homolog IDs','Reference Length (bp)','bases mapped','mismatches','Region Coverage','reads <{} mismatchs\nratio'.format(mismatchthrld),'Mean covered Depth\nmismatch ratio corrected','Normalized covered Depth\nmismatch ratio corrected','Warnings']
        xls_samplename_txt+="file: {}\n".format(xls_filename)
        
        for sample in samples_names:
            xls_samplename_txt+="\t{}\n".format(sample)
            sheet = workbook.add_sheet(sample[:31])
            readsNbr_cell="Reads count: {}".format(reads_nbr[sample])
            sheet.write(0,0,readsNbr_cell,style)

            row=0
            col=1
            for h in headers:
                sheet.write(row,col,h,hstyle)
                col += 1


            raw=1
            col=0
            
            items = sorted(outstats[sample].items(), key = SortKey, reverse=config['sortOrderReverse'])
            for allele,stats in items:
                cstyle = style
                print(sample,allele)
                reads_ratio=float(stats['reduced_stats']['reads mapped'])/float(stats['stats']['reads mapped'])
                values=[allele,stats['gscore'],stats['stats']['error rate'],stats['coverage']['depth_covered_mean'],stats['coverage']['NORM_depth_covered_mean'],refinfo[allele]['grpRef'],refinfo[allele]['len_seq'],stats['stats']['bases mapped'],stats['stats']['mismatches'],stats['coverage']['covered_pct'],reads_ratio,stats['coverage']['depth_covered_mean']*reads_ratio,stats['coverage']['NORM_depth_covered_mean']*reads_ratio,""]
                if stats['IsPositive']:
                    if refinfo[allele]['Paralog']:
                        cstyle=orangeStyle
                    else:
                        cstyle=greenStyle
                col=0
                for v in values:
                    sheet.write(raw,col,v,cstyle)
                    col+=1
                raw+=1

        workbook.save(xls_filename)

    with open(os.path.join(results_folder,"xls_samples_index.txt"),"w") as outindex:
        outindex.write(xls_samplename_txt)
        
def write_genotypTXT(out_folder,results_folder,outstats,pk_refdatabase,sep="\t"):
    sample_list=sorted(outstats.keys())
    refinfo=ld.pickle_loadDic(pk_refdatabase)
    
    posAllele=[(sample,[allele for allele,stats in outstats[sample].items() if (stats['IsPositive'] and not refinfo[allele]['Paralog'])]) for sample in sample_list]

    postGroups=[(p[0],set([refinfo[allele]['grpRef'] for allele in p[1] if refinfo[allele]['grpRef'] is not None])) for p in posAllele]

    all_Alleles=set()
    for al in posAllele:
        all_Alleles.update(al[1])
    all_Alleles=sorted(all_Alleles)

    allele_counts=[0]*len(all_Alleles)
    out_allele=sep.join(['Sample']+all_Alleles+['Nb Alleles'])+"\n"
    for sample, allele_list in posAllele:
        binary_list=[int((allele in allele_list)) for allele in all_Alleles]
        for i,v in enumerate(binary_list):
            allele_counts[i]+=v
        out_allele+=sample+sep+sep.join([str(v) for v in binary_list])+sep+"{}".format(sum(binary_list))+"\n"
    out_allele+="Nb Samples"+sep+sep.join([str(v) for v in allele_counts])+"\n"

    with open(os.path.join(results_folder,"putative_alleles.txt"),"w") as outpalleles:
        outpalleles.write(out_allele)
    
    all_Groups=set()
    for gr in postGroups:
        print(gr)
        all_Groups.update(gr[1])
    all_Groups=sorted(all_Groups)

    group_counts=[0]*len(all_Groups)
    out_group=sep.join(['Sample']+all_Groups+['Nb Groups'])+"\n"
    for sample, group_list in postGroups:
        binary_list=[int((group in group_list)) for group in all_Groups]
        for i,v in enumerate(binary_list):
            group_counts[i]+=v
        out_group+=sample+sep+sep.join([str(v) for v in binary_list])+sep+"{}".format(sum(binary_list))+"\n"
    out_group+="Nb Samples"+sep+sep.join([str(v) for v in group_counts])+"\n"

    with open(os.path.join(results_folder,"putative_groups.txt"),"w") as outpgroups:
        outpgroups.write(out_group)
    
def SortKey(tup):
    key, d = tup
    return d[config['SortKey']]
