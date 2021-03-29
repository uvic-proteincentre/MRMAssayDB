import os,subprocess,psutil,re,shutil,datetime,sys,glob
import csv

import json
import pandas as pd

class GenerateStatIDUpdateDiff():
	def __init__(self):
		self.listOfResourceFile={'PeptideTracker':'final_report_peptrack_data.csv',\
		'PASSEL':'final_report_passel_data.csv','CPTAC':'final_report_cptac_data.csv',\
		'PanoramaWEB':'final_report_panorama_data.csv','SRM':'final_report_srm_data.csv'}
		self.mergedictdata={}
		#get home directory path
		self.protPresencefilename='uniprotpresence'
		self.unquniidlist=[]


	def makeMergeDict(self):
		for fitem in self.listOfResourceFile:
			with open(self.listOfResourceFile[fitem],'r') as pepfile:
				for line in pepfile:
					data=line.strip()
					if not data.startswith('UniprotID'):
						info=data.split('\t')
						if '?' not in info[0].strip() or '=' not in info[0].strip() or ':' not in info[0].strip():
							if fitem in self.mergedictdata:
								self.mergedictdata[fitem].append(info[0].split('-')[0])
							else:
								self.mergedictdata[fitem]=[info[0].split('-')[0]]
			print(str(fitem),"data dictionay job done",str(datetime.datetime.now()))
		self.mergedictdata={a:list(set(b)) for a, b in self.mergedictdata.items()}

	def checkTotalUniqProtein(self):
		df=pd.read_csv(self.protPresencefilename+'.csv',  sep='\t',keep_default_na=False)
		self.unquniidlist=list(set(df['ActualUniID']))
		self.unquniidlist=[i.split('-')[0] for i in self.unquniidlist]

	def getStat(self):
		self.makeMergeDict()
		self.checkTotalUniqProtein()
		for key in self.mergedictdata:
			uniqIDs=list(set(self.mergedictdata[key]))
			diff=list(set(uniqIDs)-set(self.unquniidlist))
			print('Difference ID mapping for '+key+' is '+str(len(diff)))

getdiff=GenerateStatIDUpdateDiff()
getdiff.getStat()