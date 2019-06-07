#!/bin/bash
#installation script for Ubuntu image
sudo apt-get update
sudo apt-get install -y python python-tk python-pip python-qt4* python-lxml python-six
sudo python2 -m pip install --upgrade pip
sudo python2 -m pip install pyyaml biopython matplotlib numpy scipy ipython jupyter pandas sympy nose pysam
sudo python2 -m pip install --upgrade ete3
echo "export PATH=$PATH:/ifb/NGSgenotyp" > /etc/profile.d/NGSgenotyp_PATH.sh 
