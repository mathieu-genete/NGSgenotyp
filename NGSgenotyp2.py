#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NGSgenotyp2

@author: Mathieu GENETE
22/03/2021
"""

import os
import sys
import importlib

AppPath = os.path.realpath(__file__)
App_Folder = os.path.dirname(AppPath)
Features_Path = os.path.join(App_Folder,"features/")
pymod_path = os.path.join(App_Folder,"pymod/")
sys.path.insert(0,Features_Path)
sys.path.insert(0,pymod_path)

import yaml

__version__="2.0.1"

def run(ArgVal):
    configs={}

    configs['Features_Path']=Features_Path
    configs['pymod_path']=pymod_path
    configs['AppPath']=AppPath
    configs['App_Folder']=App_Folder
    
    configs['tools_folder']=os.path.join(App_Folder,'tools')

    configs['configs_folder']=os.path.join(App_Folder,"configs")
    configs['tools_config']=os.path.join(configs['configs_folder'],"tools_config.yaml")

    configs['tools_avail']=get_tools_list(configs['tools_config'],configs['tools_folder'])
        
    feature_list=[f[:f.rfind('.')] for f in os.listdir(Features_Path) if (f[-3:]==".py" and f[0] not in ['#','.'] and check_NGS_feature(os.path.join(Features_Path,f))) ]
    
    if len(ArgVal)>1:
        feature=ArgVal[1]
        feature_args=ArgVal[2:]
        if feature == "help":
            show_help(feature_list)

        if feature == "version":
            print("version",__version__)
            sys.exit()
            
        if feature in feature_list:
            NGSgenotyp_module = importlib.import_module(feature)
            NGSgenotyp_module.main(feature_args,configs)
        else:
            print("{} : unknown feature\n".format(feature))
            show_help(feature_list)
    else:
        show_help(feature_list)
    
def show_help(feature_list):
    print("usage: NGSgenotyp2.py [feature]\n")
    print("features list:")
    print("\thelp\t\t--\tshow this help message")
    print("\tversion\t\t--\tprogram version")
    for feature in sorted(feature_list):
        description=get_feature_description(os.path.join(Features_Path,"{}.py".format(feature)))
        feature+=" "*(max([len(s) for s in feature_list])-len(feature))
        print("\t{}\t--\t{}".format(feature,description))
    sys.exit()

def get_tools_list(tools_conf_file,tools_folder):
    tools_config=yaml.load(open(tools_conf_file),Loader=yaml.FullLoader)
    for tname,tconfig in tools_config.items():
        tconfig['folder']=os.path.join(tools_folder,tconfig['bin_folder'])
    return tools_config
    
def get_feature_description(feature_file):
    description="NO description"
    with open(feature_file,'r') as ft:
        for line in ft:
            if "@NGS-description" in line:
                tmp=line.strip().split(':')
                description=tmp[1].strip()
    return description

def check_NGS_feature(feature_py):
    with open(feature_py,'r') as ft:
        for line in ft:
            if "@NGS-feature" in line:
                return True
    return False


if __name__=="__main__":
    run(sys.argv)
