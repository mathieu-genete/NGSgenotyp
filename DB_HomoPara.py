#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os
import sys

#Library path
sys.path.insert(0, "PyLib/")

import xlrd
import argparse
from Bio import SeqIO

def Load_Paralogs(xlsfile):
	
	sh_paralogs = 1
	
	wb = xlrd.open_workbook(xlsfile)	
	
	sh = wb.sheet_by_index(sh_paralogs)
	
	row_nbr = sh.nrows
	
	list_paralogs = [ v.lstrip().rstrip() for v in sh.col_values(0)]
				
	return list_paralogs
	
def normalize_Str(text):
	rslt = text.replace(' ','_')
	return rslt
	
def Load_Homo_Group_Color(xlsfile):
	
	sh_homologues = 0
	
	wb = xlrd.open_workbook(xlsfile)	
	
	sh = wb.sheet_by_index(sh_homologues)
	
	col_nbr = sh.ncols
	row_nbr = sh.nrows
	
	col_grpRef = [ v.lstrip().rstrip() for v in sh.col_values(0)]
	
	col_species = []
	
	for index in range(1,col_nbr):
		col_species.append([ v.lstrip().rstrip() for v in sh.col_values(index)])
	
	returnDict = {}
	
	for species in col_species:
		index = 1
		specie_name = normalize_Str(species[0])
		for val in species[1:]:
			if val != '':
				grp_Ref = normalize_Str(col_grpRef[index])
				grp_Id = int(grp_Ref[1:-3])
				Haplo_Id = int(grp_Ref[-3:])
				returnDict[normalize_Str(val)] = {'specie':specie_name,'grpRef':grp_Ref, 'gprId':grp_Id, 'HaploId':Haplo_Id}
			index += 1
				
	return returnDict

def get_paramsFromSeqID(seqID):
	tabID = seqID.split("|")
	params = {}
	
	if len(tabID)>1:
		for p in tabID[1:]:
			kv = p.split("=")
			if len(kv)==2:
				params[kv[0]]=kv[1]
		
	return tabID[0],params

def main():
	
	xlsfile = args.xlsHomoPara
	dbFile = args.databaseFastaFile
	outFasta = args.outfasta
	
	if dbFile==outFasta:
		print "ERROR: input and output file name should be different"
		sys.exit(1)
	
	Para_List=Load_Paralogs(xlsfile)
	Homo_List=Load_Homo_Group_Color(xlsfile)
	
	dbFasta = SeqIO.parse(open(dbFile,'r'),'fasta')
	with open(outFasta, 'w') as output_handle:
		for seq in dbFasta:
			SeqID,params=get_paramsFromSeqID(seq.id)
			params = {}
			if SeqID in Para_List:
				params['Paralog']=1
			if SeqID in Homo_List:
				params.update(Homo_List[SeqID])
			TabP=[str(k)+"="+str(v) for k,v in params.items()]
			strParams = "|".join(TabP)
		
			fastaID = SeqID
		
			if strParams!="":
				fastaID += "|"+normalize_Str(strParams)
		
			seq.id = fastaID
			seq.description = ''
			SeqIO.write(seq, output_handle, "fasta")
	
if	__name__ == '__main__':
	
	description = """ Database completion with Homologues and paralogues informations """
	
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument("-v", "--verbose", action="store_const", const="-v", default="")	
	parser.add_argument("-o","--outfasta", help="destination file for database", required=True)
	parser.add_argument("-x","--xlsHomoPara", help="xls file with homolog ans paralog informations", required=True)
	parser.add_argument("-d","--databaseFastaFile", help="database fasta file", required=True)
	args = parser.parse_args()

	main()
