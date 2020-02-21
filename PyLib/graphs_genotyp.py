#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import pickle
import operator
import numpy as np
import math
import random
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import lognorm
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.switch_backend('agg')

__version__ = "1.0.1"

args = None

def run(ArgsVal):

    global args

    description = """ graphic for genotyp """

    parser = argparse.ArgumentParser(prog="graphs_genotyp",description=description)
        
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument("-g","--sortbygroup",help="sort by group",action="store_const", const=True, default=False)
    parser.add_argument("-S","--ShowScore",help="print genotyp scores on plots",action="store_const", const=True, default=False)

    parser.add_argument("-s","--samtviewstats", help="samview_stats.p file to analyse", required = True)
    parser.add_argument("-d","--graphsPerDocument", help="number of graphs per documents (15 by default)", type=int, default=15)
    parser.add_argument("-o","--outpdf", help="out pdf filename", required = True)
    parser.add_argument("-f","--refFilterName", help="filter for reference", nargs="*", default=[])

    args = parser.parse_args(ArgsVal)

    refFilterName = args.refFilterName
    outpdf = args.outpdf.strip(".pdf")
    graphs_per_page = args.graphsPerDocument

    stats = pickle.load(open(args.samtviewstats,'rb'))
    indivNbr=len(stats.keys())
    n=0
    nbrPdf = int(math.ceil(float(indivNbr)/graphs_per_page))

    print "{} - version {}".format(os.path.basename(__file__),__version__)
    print " ".join(sys.argv)
    print "*** {} pdf files will be generated ***".format(nbrPdf)

    pdfTab=[]
    pdfFileNameList=[]
    for i in range(0,nbrPdf):
        pdfFileName = "{}_{}.pdf".format(outpdf,i+1)
        pdfFileNameList.append(os.path.basename(pdfFileName))
        pdfTab.append(PdfPages(pdfFileName))

    idpdf=0
    distribDatas=[]
    groupsList=set()
    for v in stats.values():
        for w in v.values():
            groupsList.add(w['Group ID'])

    colorsByGroupID = colorsByKeys(list(groupsList),colormap='gist_ncar')

    for indiv in sorted(stats.keys()):
        paralogs={'scores':[],'labels':[],'mean_cover':[]}
        references={'scores':[],'labels':[],'Grp':[]}
        refsGrp = {r:stats[indiv][r]['Group ID'] for r in sorted(stats[indiv].keys())}
        if args.sortbygroup:
            refsSorted = [(k,v) for k,v in sorted(refsGrp.items(), key=operator.itemgetter(1))]
        else:
            refsSorted = [(k,v) for k,v in sorted(refsGrp.items(), key=operator.itemgetter(0))]
        for Tref in refsSorted:
            ref=Tref[0]
            if len(refFilterName)>0 and not any([True for v in refFilterName if v.upper() in ref.upper()]) and not stats[indiv][ref]['IsParalog']:
                continue
            score = stats[indiv][ref]['NormScore']
            mean_cover=stats[indiv][ref]['mean cover']
            GS_Threshold = stats[indiv][ref]['GS_Threshold']

            if stats[indiv][ref]['IsParalog'] and stats[indiv][ref]['IsPositiv'] and mean_cover>0:
                paralogs['labels'].append(ref)
                paralogs['scores'].append(score)
                paralogs['mean_cover'].append(mean_cover)
            elif mean_cover>0 and not stats[indiv][ref]['IsParalog']:
                references['labels'].append(ref)
                references['scores'].append(score)
                references['Grp'].append(Tref[1])
        if n%graphs_per_page==0 and n>0:
            pdfTab[idpdf].close()
            idpdf+=1
        n+=1

        pdfTab[idpdf].savefig(plot_graph(indiv,paralogs,references,GS_Threshold,colorsByGroupID))
        print "{} => plot {} - {}/{}".format(pdfFileNameList[idpdf],indiv,n,indivNbr)

    pdfTab[idpdf].close()

def get_cmap(n, name):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def colorsByKeys(inList,colormap='gist_rainbow'):
    AllKeys = sorted(list(set(inList)))
    colorSet = get_cmap(len(AllKeys),colormap)
    typesColor = {AllKeys[i]:colorSet(i) for i in range(0,len(AllKeys))}
    return typesColor

def plot_hist(indiv,datas):
    x= np.array(datas)
    unique, counts = np.unique(x,return_counts=True)
    shape, loc, scale = lognorm.fit(x, floc=0)
    x2 = np.linspace(min(unique),max(unique))
    p = lognorm.pdf(x2,shape, loc=loc, scale=scale)
    plt.clf()
    plt.hist(x,normed=True,bins=50)
    plt.plot(x2,p,'k')
    plt.savefig("hist_{}.pdf".format(indiv))

def get_coordsFromListGrp(GroupList):
    tmp=''
    coord=[]
    listCoord={}
    for i,c in enumerate(GroupList):
        if c!=tmp:
            if len(coord)<1:
                coord.append(i)
            elif len(coord)==1:
                coord.append(i-1)
                listCoord[tmp]=tuple(coord)
                coord=[]
                coord.append(i)
            tmp=c
    return listCoord

def plot_graph(indiv,paralogs,references,GS_Threshold,colorsByGroupID):
    plt.clf()
    fig, ax = plt.subplots()
    datas = paralogs['scores'] + references['scores']
    labels = paralogs['labels'] + references['labels']

    parScores = np.array(paralogs['scores'])
    parMed = np.median(parScores)
    parMean = np.mean(parScores)
    paralogsColor='white'
    referencesColor='gray'
    if args.sortbygroup:
        colors = [paralogsColor]*len(paralogs['labels']) + [colorsByGroupID[g] for g in references['Grp']]
        bracketsCoords = get_coordsFromListGrp(['Paralogs']*len(paralogs['labels'])+references['Grp'])
    else:
        colors = [paralogsColor]*len(paralogs['labels']) + [referencesColor]*len(references['labels'])

    patterns = ["///"]*len(paralogs['labels']) + [""]*len(references['labels'])
    for i in range(len(datas)):
        ax.bar(i,datas[i],color=colors[i],edgecolor='black',hatch=patterns[i])
        if args.ShowScore:
            ax.text(i, 1.05*datas[i],"{:.5E}".format(datas[i]), ha='center', va='bottom',rotation=90,fontsize=6)

    ax.set_xticks(range(0,len(datas),1))
    ax.set_xticklabels(labels,rotation=90, fontsize=5)

    if args.sortbygroup:
        ToRoman = {0:'',1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V',6: 'VI', 7: 'VII', 8: 'VIII', 9: 'IX'}
        ymax = max(datas)*(1+0.03)
        for k,v in bracketsCoords.items():
            x1=v[0]-0.25
            x2=v[1]+0.25
            if k!='Paralogs' and k!='':
                BraColor = colorsByGroupID[k]
                braText=ToRoman[int(k[1])]
            else:
                BraColor = 'black'
                braText=k
            plt.text(x1+(x2-x1)/2.0,ymax*(1+0.01),braText,horizontalalignment='center',verticalalignment='bottom',fontname='DejaVu Sans Mono',fontsize=6,color='grey')
            plt.plot((x1,x2),(ymax,ymax),color=BraColor,linewidth=1)
            plt.plot((x1,x1),(ymax,ymax*(1-0.02)),color=BraColor,linewidth=1)
            plt.plot((x2,x2),(ymax,ymax*(1-0.02)),color=BraColor,linewidth=1)
        plt.text(len(datas),ymax,"Groups of transpesifics alleles\n with allelic class",horizontalalignment='left',verticalalignment='bottom',fontsize=7,color='black')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.title("Genotyp Scores for {}".format(indiv))
    ax.set_ylabel("Genotyp score")

    if not math.isnan(parMed):
        ax.axhline(y=parMed/2,linestyle='--',color='black', linewidth=0.5)
        ax.axhline(y=parMed,color='black', linewidth=0.5)
        plt.text(len(datas),parMed,"paralogs\nscore median",horizontalalignment='left',verticalalignment='bottom',fontsize=10)
        plt.text(len(datas),parMed/2,"half score med.",horizontalalignment='left',verticalalignment='bottom',fontsize=10)

    #genotyp score threshold
    if not math.isnan(GS_Threshold):
        plt.axhspan(0, GS_Threshold, facecolor='0.5', alpha=0.2)
        plt.text(len(datas),GS_Threshold/2,"Genotyp Score\nThreshold",horizontalalignment='left',verticalalignment='center',fontsize=10)

    fig.set_size_inches(15,10)
    return fig

if __name__=='__main__':    
    ArgsVal = sys.argv[1:]
    run(ArgsVal)
