#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker

plt.switch_backend('agg')

def plot_distrib(datas,title,xtitle,ytitle,outfile,n_bins=80):
    plt.hist(datas,bins=n_bins)
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title(title)
    plt.savefig(outfile)
