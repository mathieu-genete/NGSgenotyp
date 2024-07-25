#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gzip
import yaml
import pickle

def yaml_dumpDic(out_yamlfile,dic):
    with open(out_yamlfile,'w') as out_yaml:
        yaml.dump(dic, out_yaml, default_flow_style=False)

def yaml_loadDic(yamlfile):
    return yaml.load(open(yamlfile),Loader=yaml.FullLoader)
    
def pickle_dumpDic(outpickle,dic):
    pf = gzip.open(outpickle,'wb')
    pickle.dump(dic,pf,protocol=pickle.HIGHEST_PROTOCOL)
    pf.close()

def pickle_loadDic(picklefile):
    pf = gzip.open(picklefile,'rb')
    retval = pickle.load(pf)
    pf.close()
    return retval
