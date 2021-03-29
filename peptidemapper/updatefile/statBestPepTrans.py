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

def statBestPepTrans(filename):
	countProt=0
	homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
	filepath = os.path.join(homedir, 'src/mappermotherfile', filename)


	countLine=0
	modrepfile="bestreprot.csv"
	moutput= open(modrepfile,'w')
	with open(filepath,'r') as repfile:
		for line in repfile:
			info=line.rstrip().split('\t')
			if 'UniProtKB Accession' in info:
				info.append("Peptide Occurrence")
				info.append("Best Transition")
				info.append("Summary Transition")
				moutput.write(('\t'.join(info))+'\n')
			else:
				info.extend(['NA']*3)
				countLine+=1
				countOcurancePep=0
				transitioninfopeptrack=[]
				besttransionpeptrack=[]
				transitioninfopassel=[]
				besttransionpassel=[]
				transitioninfosrm=[]
				besttransionsrm=[]
				transitioninfocptac=[]
				besttransioncptac=[]
				transitioninfopan=[]
				besttransionpan=[]
				if len(info[10].strip())>0 and (str(info[10]).strip()).lower() !='na':
					countOcurancePep+=1
				if len(info[11].strip())>0  and (str(info[11]).strip()).lower() !='na':
					transinfo=(info[11].strip()).split(',')
					if len(transinfo)>1:
						for t in transinfo[1:]:
							if t.upper() !="NA":
								subtransinfo=t.split('|')
								intrtument=subtransinfo[0]
								fragionlist=''
								if 'agilent' in subtransinfo[0].lower():
									fragionlist=subtransinfo[6].split('+')
								if 'qtrap' in subtransinfo[0].lower():
									fragionlist=subtransinfo[7].split('+')
								if 'thermo' in subtransinfo[0].lower():
									fragionlist=subtransinfo[5].split('+')
								charge=0
								try:
									if 'Q3' in str(transinfo[0].split('|')[-1]).strip().upper():
										try:
											charge=int(float(str(subtransinfo[-1]).strip()))
										except ValueError:
											pass
									else:
										if len(fragionlist[-1])==0:
											charge=fragionlist.count('')
										else:
											if len(fragionlist)==1:
												charge=1
											else:
												charge=int(fragionlist[-1])
									fragion=fragionlist[0].strip()
									if 'y' in fragion.lower():
										besttransionpeptrack.append([charge,fragion])
									transitioninfopeptrack.append([intrtument,"PeptideTracker",charge,fragion])
								except IndexError:
									pass
					besttransionpeptrack=map(list,set(map(tuple,besttransionpeptrack)))
					transitioninfopeptrack=map(list,set(map(tuple,transitioninfopeptrack)))
				if len(info[12].strip())>0 and (str(info[12]).strip()).lower() !='na':
					countOcurancePep+=1
				if len(info[13].strip())>0  and (str(info[13]).strip()).lower() !='na':
					transinfo=(info[13].strip()).split(',')
					if len(transinfo)>1:
						for t in transinfo[1:]:
							subtransinfo=t.split('|')
							intrtument=subtransinfo[0]
							charge=0
							fragionlist=subtransinfo[4].split('+')
							if 'Q3' in str(transinfo[0].split('|')[-1]).strip().upper():
								try:
									charge=int(float(str(subtransinfo[-1]).strip()))
								except ValueError:
									pass
							else:
								if len(fragionlist[-1])==0:
									charge=fragionlist.count('')
								else:
									if len(fragionlist)==1:
										charge=1
									else:
										charge=int(fragionlist[-1])
							fragion=fragionlist[0].strip()
							if 'y' in fragion.lower():
								besttransionpassel.append([charge,fragion])
							transitioninfopassel.append([intrtument,"Passel",charge,fragion])
					besttransionpassel=map(list,set(map(tuple,besttransionpassel)))
					transitioninfopassel=map(list,set(map(tuple,transitioninfopassel)))
				if len(info[14].strip())>0 and (str(info[14]).strip()).lower() !='na':
					countOcurancePep+=1
				if len(info[15].strip())>0  and (str(info[15]).strip()).lower() !='na':
					transinfo=(info[15].strip()).split(',')
					if len(transinfo)>1:
						for t in transinfo[1:]:
							subtransinfo=t.split('|')
							intrtument=subtransinfo[0]
							charge=0
							fragionlist=subtransinfo[-4].split('+')
							if 'Q3' in str(transinfo[0].split('|')[-1]).strip().upper():
								try:
									charge=int(float(str(subtransinfo[-1]).strip()))
								except ValueError:
									pass
							else:
								if len(fragionlist[-1])==0:
									charge=fragionlist.count('')
								else:
									if len(fragionlist)==1:
										charge=1
									else:
										charge=int(fragionlist[-1])
							fragion=fragionlist[0].strip()
							if 'y' in fragion.lower():
								besttransionsrm.append([charge,fragion])
							transitioninfosrm.append([intrtument,"SRMAtlas",charge,fragion])
					besttransionsrm=map(list,set(map(tuple,besttransionsrm)))
					transitioninfosrm=map(list,set(map(tuple,transitioninfosrm)))
				if len(info[16].strip())>0 and (str(info[16]).strip()).lower() !='na':
					countOcurancePep+=1
				if len(info[17].strip())>0  and (str(info[17]).strip()).lower() !='na':
					transinfo=(info[17].strip()).split(',')
					if len(transinfo)>1:
						for t in transinfo[1:]:
							subtransinfo=t.split('|')
							intrtument=subtransinfo[0]
							charge=0
							fragionlist=subtransinfo[4].split('+')
							if 'Q3' in str(transinfo[0].split('|')[-1]).strip().upper():
								try:
									charge=int(float(str(subtransinfo[-1]).strip()))
								except ValueError:
									pass
							else:
								if len(fragionlist[-1])==0:
									charge=fragionlist.count('')
								else:
									if len(fragionlist)==1:
										charge=1
									else:
										charge=int(fragionlist[-1])
							fragion=fragionlist[0].strip()
							if 'y' in fragion.lower():
								besttransioncptac.append([charge,fragion])
							transitioninfocptac.append([intrtument,"CPTAC",charge,fragion])
					besttransioncptac=map(list,set(map(tuple,besttransioncptac)))
					transitioninfocptac=map(list,set(map(tuple,transitioninfocptac)))
				if len(info[18].strip())>0 and (str(info[18]).strip()).lower() !='na':
					countOcurancePep+=1
				if len(info[19].strip())>0  and (str(info[19]).strip()).lower() !='na':
					transinfo=(info[19].strip()).split(',')
					if len(transinfo)>1:
						for t in transinfo[1:]:
							subtransinfo=t.split('|')
							intrtument="NA"
							charge=0
							try:
								fragionlist=subtransinfo[-3].split('+')
								intrtument=subtransinfo[0]
								if 'Q3' in str(transinfo[0].split('|')[-1]).strip().upper():
									try:
										charge=int(float(str(subtransinfo[-1]).strip()))
									except ValueError:
										pass
								else:
									if len(fragionlist[-1])==0:
										charge=fragionlist.count('')
									else:
										if len(fragionlist)==1:
											charge=1
										else:
											charge=int(fragionlist[-1])
								fragion=fragionlist[0].strip()
								if 'y' in fragion.lower():
									besttransionpan.append([charge,fragion])
								transitioninfopan.append([intrtument,"PanoramaWeb",charge,fragion])
							except IndexError:
								pass
					besttransionpan=map(list,set(map(tuple,besttransionpan)))
					transitioninfopan=map(list,set(map(tuple,transitioninfopan)))
				info[-3]=countOcurancePep

				finalbesttransion=besttransionpeptrack+besttransionpassel+besttransionsrm+besttransioncptac+besttransionpan
				finalbesttransiondata=[str(i[1])+"+"*i[0] for i in finalbesttransion]
				
				if len(finalbesttransiondata)>0:
					top5best=list(collections.Counter(finalbesttransiondata).items())
					top5best=sorted(top5best, key=lambda l:l[1], reverse=True)
					top5best=top5best[:5]
					info[-2]='|'.join([i[0] for i in top5best])
				finaltransitioninfo=transitioninfopeptrack+transitioninfopassel+transitioninfosrm+transitioninfocptac+transitioninfopan
				finaltransitiondata=['|'.join(map(str, i)) for i in finaltransitioninfo]
				if len(finaltransitiondata)>0:
					allFragIons=[]
					for fl in finaltransitiondata:
						allFragIons.append(str(fl.split('|')[-1]).strip().lower())
					fragionCounts=list(collections.Counter(allFragIons).items())
					tempFargDic={k[0]:k[1] for k in fragionCounts}

					tempfinaltransitiondata=[]
					for tfl in finaltransitiondata:
						tempTransData=tfl.split('|')
						tempTransData.append(str(tempFargDic[str(tempTransData[-1]).strip().lower()]))
						tempfinaltransitiondata.append('|'.join(tempTransData))
					info[-1]='.'.join(tempfinaltransitiondata)
				moutput.write(('\t'.join(map(str, info)))+'\n')

	moutput.close()
	movefilepath=os.path.join(homedir, 'updatefile', filename)
	os.rename(modrepfile,filename)
	shutil.move(movefilepath,filepath)
