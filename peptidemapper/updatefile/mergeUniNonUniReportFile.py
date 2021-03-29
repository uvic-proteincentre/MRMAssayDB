import os,subprocess,psutil,re,shutil,datetime,sys,glob
from socket import error as SocketError
import errno
import random, time
import csv
import json
import pandas as pd
from collections import Counter
from itertools import combinations
import operator
import collections


homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
filename='ReportBook_mother_file.csv'
pepfilepath = os.path.join(homedir, 'src/mappermotherfile', filename)
uniprotpepresult=[]
with open(pepfilepath) as pepresult:
	for line in pepresult:
		data=line.strip()
		info=data.split('\t')
		if 'UniProtKB Accession' in info:
			info.append('UniprotKb entry status')
			uniprotpepresult.append(info)
		else:
			info.append('Yes')
			uniprotpepresult.append(info)

nonunipepfilepath = os.path.join(homedir, 'updatefile','NonUniprot_ReportBook_mother_file.csv')
with open(nonunipepfilepath) as nonunipepresult:
	for nonuniline in nonunipepresult:
		nonunidata=nonuniline.strip()
		nonuniinfo=nonunidata.split('\t')
		uniprotpepresult.append(nonuniinfo)

mergerepfile="mergereprot.csv"
moutput= open(mergerepfile,'w')
totalpepseq=[]
for item in uniprotpepresult:
	totalpepseq.append(item[5].strip().upper())
	if 'UniprotID' not in item[0]:
		if len(item)==34:
			item.insert(18,'NA')
			item.insert(18,'NA')
			item.insert(18,'NA')
		moutput.write(('\t'.join(map(str, item)))+'\n')
moutput.close()

movefilepath=os.path.join(homedir, 'updatefile', filename)
os.rename(mergerepfile,filename)
shutil.move(movefilepath,pepfilepath)

unqtotalpepseq=list(set(totalpepseq))

calfilename='totalpepassay.py'
calmovefilepath=os.path.join(homedir, 'updatefile', calfilename)
calfilepath = os.path.join(homedir, 'src/pepmapperapp', calfilename)

totalassayresult={}
totalassayresult['totalpepassay']=unqtotalpepseq

calfileoutput=open(calfilename,'w')
calfileoutput.write("totalpepassay=")
calfileoutput.write(json.dumps(totalassayresult))
calfileoutput.close()
shutil.move(calmovefilepath,calfilepath)