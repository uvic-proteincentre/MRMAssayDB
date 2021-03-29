import urllib,urllib2
from bioservices.kegg import KEGG
import os,subprocess,psutil,re,shutil,datetime,sys,glob
from operator import itemgetter
import numpy as np
import random, time
from itertools import count, groupby
import pandas as pd
import csv
def split_lines(sentence, step=4):
	c = count()
	chunks = sentence.split()
	return [' '.join(g) for k, g in groupby(chunks, lambda i: c.next() // step)]


statKeggDic={}
speciesdic={}
homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
filename='ReportBook_mother_file.csv'
filepath = os.path.join(homedir, 'src/mappermotherfile', filename)
calfilename='calculationprog.py'
calmovefilepath=os.path.join(homedir, 'updatefile', calfilename)
calfilepath = os.path.join(homedir, 'src/pepmapperapp', calfilename)

with open(filepath) as repcsvfile:
	reppepreader = csv.DictReader(repcsvfile, delimiter='\t')
	for reppeprow in reppepreader:
		statuniID=str(reppeprow['UniProtKB Accession']).strip()
		statPathway=str(reppeprow['Kegg Pathway Name']).strip()
		statSpecies=str(reppeprow['Organism']).strip()
		statSpeciesID=str(reppeprow['Organism ID']).strip()
		speciesdic[statSpeciesID]=statSpecies
		if len(statPathway)>0:
			statPathwayList=statPathway.split('|')
			for stpItem in statPathwayList:
				if stpItem.lower() !='na':
					keggid=stpItem.strip()+'|'+statSpeciesID
					if statKeggDic.has_key(keggid):
						statKeggDic[keggid].append(statuniID.strip())
					else:
						statKeggDic[keggid] =[statuniID.strip()]

#kegg pathway name extraction


sys.path.append(os.path.join(homedir, 'src/pepmapperapp'))
from calculationprog import *
pepfinalresult=finalresult['prodataseries']
#kegg pathway name extraction

statCoverageDic={}
for kskey in statKeggDic:
	keggpathwayname=(kskey.strip()).split('|')[0]
	tempUniqKeggUniIDList=list(set(statKeggDic[kskey]))
	cptac=[]
	panweb=[]
	passel=[]
	srmatlas=[]
	peptrack=[]
	for ckey in pepfinalresult:
		if ckey == "PeptideTracker":
			peptrack=list(set(pepfinalresult[ckey]).intersection(tempUniqKeggUniIDList))
		if ckey == "PASSEL":
			passel=list(set(pepfinalresult[ckey]).intersection(tempUniqKeggUniIDList))
		if ckey == "SRMAtlas":
			srmatlas=list(set(pepfinalresult[ckey]).intersection(tempUniqKeggUniIDList))
		if ckey == "CPTAC":
			cptac=list(set(pepfinalresult[ckey]).intersection(tempUniqKeggUniIDList))
		if ckey == "PanoramaWeb":
			panweb=list(set(pepfinalresult[ckey]).intersection(tempUniqKeggUniIDList))
	tempcptac=len(list(set(cptac)))
	temppanweb=len(list(set(panweb)))
	temppassel=len(list(set(passel)))
	temppeptrack=len(list(set(peptrack)))
	tempsrmatlas=len(list(set(srmatlas)))
	tempTotal=len(list(set(tempUniqKeggUniIDList)))
	templist=[keggpathwayname,tempTotal,temppeptrack,tempcptac,temppassel,tempsrmatlas,temppanweb]
	for tempItem in tempUniqKeggUniIDList:
		if statCoverageDic.has_key(tempItem):
			statCoverageDic[tempItem].append(templist)
		else:
			statCoverageDic[tempItem] =[templist]

keggstatfileCoveragefile="keggcoverage.csv"
koutput= open(keggstatfileCoveragefile,'w')
with open(filepath,'r') as f:
	for line in f:
		info=line.rstrip().split('\t')
		if 'UniProtKB Accession' in info:
			info.append("Kegg Coverage")
			koutput.write(('\t'.join(info))+'\n')
		else:
			if info[0].split('-')[0] in statCoverageDic:
				tempdata=str(statCoverageDic[info[0].split('-')[0]])
				info.append(tempdata)
				koutput.write(('\t'.join(info))+'\n')
			else:
				tempdata='NA'
				info.append(tempdata)
				koutput.write(('\t'.join(info))+'\n')
koutput.close()
movefilepath=os.path.join(homedir, 'updatefile', filename)
os.rename(keggstatfileCoveragefile,filename)
shutil.move(movefilepath,filepath)