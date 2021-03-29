import os,subprocess,psutil,re,shutil,datetime,sys,glob
from operator import itemgetter
import numpy as np
import random, time
from itertools import count, groupby
import pandas as pd
import csv
import pickle
import cPickle
import itertools
import collections

countProt=0
homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
filename='ReportBook_mother_file.csv'
filepath = os.path.join(homedir, 'src/mappermotherfile', filename)


countLine=0
modrepfile="moduserreprot.csv"
moutput= open(modrepfile,'w')
with open(filepath,'r') as repfile:
	for line in repfile:
		info=line.rstrip().split('\t')
		if 'UniProtKB Accession' in info:
			moutput.write(('\t'.join(info))+'\n')
		else:
			if len(info[12].strip())>0 and (str(info[12]).strip()).lower() !='na':
				#peptidetracker fragment ion always in lower case so skip this step
				info[12]=info[12]
			if len(info[14].strip())>0 and (str(info[14]).strip()).lower() !='na':
				transinfo=(info[14].strip()).split(',')
				if len(transinfo)>1:
					tempT=[]
					tempT.append(transinfo[0])
					for t in transinfo[1:]:
						subtransinfo=t.split('|')
						subtransinfo[4]=subtransinfo[4].replace(" ","")
						subtransinfo[4]=subtransinfo[4].lower()
						tempT.append('|'.join(subtransinfo))
					info[14]=','.join(tempT)
			if len(info[16].strip())>0 and (str(info[16]).strip()).lower() !='na':
				transinfo=(info[16].strip()).split(',')
				if len(transinfo)>1:
					tempT=[]
					tempT.append(transinfo[0])
					for t in transinfo[1:]:
						subtransinfo=t.split('|')
						subtransinfo[-4]=subtransinfo[-4].replace(" ","")
						subtransinfo[-4]=subtransinfo[-4].lower()
						tempT.append('|'.join(subtransinfo))
					info[16]=','.join(tempT)
			if len(info[18].strip())>0 and (str(info[18]).strip()).lower() !='na':
				info[18]=info[18].replace("||","|")
				transinfo=(info[18].strip()).split(',')
				if len(transinfo)>1:
					tempT=[]
					tempT.append(transinfo[0])
					for t in transinfo[1:]:
						subtransinfo=t.split('|')
						subtransinfo[-3]=subtransinfo[-3].replace(" ","")
						subtransinfo[-3]=subtransinfo[-3].lower()
						tempT.append('|'.join(subtransinfo))
					info[18]=','.join(tempT)
			if len(info[20].strip())>0 and (str(info[20]).strip()).lower() !='na':
				info[20]=info[20].replace("||","|")
				transinfo=(info[20].strip()).split(',')
				if len(transinfo)>1:
					tempT=[]
					tempT.append(transinfo[0])
					for t in transinfo[1:]:
						subtransinfo=t.split('|')
						try:
							subtransinfo[-3]=subtransinfo[-3].replace(" ","")
							subtransinfo[-3]=subtransinfo[-3].lower()
							tempT.append('|'.join(subtransinfo))
						except IndexError:
							pass
					info[20]=','.join(tempT)

			sumtrancol=[29,31,33,35,37]
			for sc in sumtrancol:
				try:
					if len(info[sc].strip())>0 and (str(info[sc]).strip()).lower() !='na':
						sumtransinfo=(info[sc].strip()).split('<br/>')
						sumtransinfo[1]=sumtransinfo[1].split(':')
						sumtransinfo[1][1]=sumtransinfo[1][1].lower()
						sumtransinfo[1]=':'.join(sumtransinfo[1])
						info[sc]='<br/>'.join(sumtransinfo)
				except IndexError:
					pass

			if len(info[39].strip())>0 and (str(info[39]).strip()).lower() !='na':
				info[39]=info[39].lower()
			if len(info[40].strip())>0 and (str(info[40]).strip()).lower() !='na':
				sumtran=(info[40].strip()).split('.')
				tempL=[]
				for x in sumtran:
					l=x.split('|')
					l[-1]=l[-1].lower()
					tempL.append('|'.join(l))
				info[40]='.'.join(tempL)
			moutput.write(('\t'.join(map(str, info)))+'\n')
moutput.close()
movefilepath=os.path.join(homedir, 'updatefile', filename)
os.rename(modrepfile,filename)
shutil.move(movefilepath,filepath)