import urllib,urllib2
from bioservices.kegg import KEGG
import os,subprocess,psutil,re,shutil,datetime,sys,glob
from operator import itemgetter
import numpy as np
import random, time
from itertools import count, groupby
import pandas as pd
import csv
import itertools
import json

class SummaryStat():
	def __init__(self):
		self.curr_dir=os.getcwd()
		self.filename='ReportBook_mother_file.csv'
		self.homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
		self.filepath = os.path.join(self.homedir, 'src/mappermotherfile', self.filename)
		self.statfilename='overallstat.py'
		self.statmovefilepath=os.path.join(self.homedir, 'updatefile', self.statfilename)
		self.statfilepath = os.path.join(self.homedir, 'src/pepmapperapp', self.statfilename)
		self.orglist=[]
		self.orglistID=[]
		self.unqorglist=[]
		self.statKeggDic={}
		self.speciesList=[]
		self.speciesProt={}
		self.speciesPep={}
		self.godic={}
		self.speciesdic={}
		self.keggpathwaycoverage=[]
		self.statsepcies=[]
		self.golist=[]
		self.sortedstatsepcies=[]
		self.sortedgolist=[]
		self.sortedkeggpathwaycoverage=[]

	def orgKeggStat(self):
		csv.field_size_limit(sys.maxsize)
		with open(self.filepath) as pepcsvfile:
			pepreader=csv.DictReader(pepcsvfile, delimiter='\t')
			for peprow in pepreader:
				if  str(peprow['UniprotKb entry status']).strip().upper()=='YES':
					self.orglist.append(' '.join(str(peprow['Organism']).strip().split(' ')[:2]))
					self.orglistID.append(str(peprow['Organism ID']).strip())
					statuniID=str(peprow['UniProtKB Accession']).strip().split('-')[0]
					statPathway=str(peprow['Kegg Pathway Name']).strip()
					statSpeciesID=str(peprow['Organism ID']).strip()
					if len(statPathway)>0:
						statPathwayList=statPathway.split('|')
						for stpItem in statPathwayList:
							if stpItem.lower() !='na':
								keggid=stpItem.strip()
								if keggid in self.statKeggDic:
									self.statKeggDic[keggid].append(statuniID.strip())
								else:
									self.statKeggDic[keggid] =[statuniID.strip()]

		self.unqorglist=list(set(self.orglist))
		self.unqorglist.sort(key=str.lower)
		self.speciesList=self.unqorglist

	def speciesGoStat(self):
		#print speciesList
		csv.field_size_limit(sys.maxsize)
		with open(self.filepath) as pepcsvfile:
			pepreader=csv.DictReader(pepcsvfile, delimiter='\t')
			for frow in pepreader:
				if str(frow['UniprotKb entry status']).strip().upper()=='YES':
					pepseq=frow['Peptide Sequence'].strip()
					Calorg=str(frow['Organism']).strip()
					CalorgID=str(frow['Organism ID']).strip()
					self.speciesdic[CalorgID]=Calorg
					acccode=None
					acccode=str(frow['UniProtKB Accession']).split('-')[0]
					if acccode != None:
						for spitem in self.speciesList:
							if spitem in Calorg:
								if spitem in self.speciesProt:
									self.speciesProt[spitem].append(acccode)
								else:
									self.speciesProt[spitem]=[acccode]

								if spitem in self.speciesPep:
									self.speciesPep[spitem].append(pepseq)
								else:
									self.speciesPep[spitem]=[pepseq]
						if frow["Go Name"].upper() !='NA' and len(str(frow["Go Name"]).strip())>0:
							goname=(str(frow["Go Name"]).strip()).split('|')
							for gitem in goname:
								if str(gitem).strip() in self.godic:
									self.godic[str(gitem).strip()].append(acccode)
								else:
									self.godic[str(gitem).strip()]=[acccode]

	def keggCoverageStat(self):
		sys.path.append(os.path.join(self.homedir, 'src/pepmapperapp'))
		from calculationprog import *
		pepfinalresult=finalresult['prodataseries']

		for kskey in self.statKeggDic:
			keggpathwayname=(kskey.strip()).split('|')[0]
			tempUniqKeggUniIDList=list(set(self.statKeggDic[kskey]))
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
			self.keggpathwaycoverage.append(templist)

	def sortData(self):
		self.speciesProt={k:len(list(set(j))) for k,j in self.speciesProt.items()}
		self.speciesPep={k:len(list(set(j))) for k,j in self.speciesPep.items()}

		self.godic={k:len(set(v)) for k, v in self.godic.items()}
		 
		for skey in self.speciesProt:
			if skey in self.speciesPep:
				self.statsepcies.append([skey,self.speciesProt[skey],self.speciesPep[skey]])

		self.sortedstatsepcies=sorted(self.statsepcies, key= itemgetter(1), reverse=True)


		for gkey in self.godic:
			self.golist.append([gkey,self.godic[gkey]])

		self.sortedgolist=sorted(self.golist, key= itemgetter(1), reverse=True)

		self.keggpathwaycoverage.sort()
		unqkeggpathwaycoverage=list(self.keggpathwaycoverage for self.keggpathwaycoverage,_ in itertools.groupby(self.keggpathwaycoverage))
		self.sortedkeggpathwaycoverage=sorted(unqkeggpathwaycoverage, key= itemgetter(1), reverse=True)

	def finalStat(self):
		self.orgKeggStat()
		self.speciesGoStat()
		self.keggCoverageStat()
		self.sortData()

		overallSumresult={}
		overallSumresult['organism']=self.unqorglist
		self.unqorglist.insert(0,"")
		self.unqorglist.remove('NA')
		overallSumresult['species']=self.unqorglist
		overallSumresult['speciesstat']=self.sortedstatsepcies
		overallSumresult['gostat']=self.sortedgolist
		overallSumresult['keggstat']=self.sortedkeggpathwaycoverage

		statfileoutput=open(self.statfilename,'w')
		statfileoutput.write("overallSumresult=")
		statfileoutput.write(json.dumps(overallSumresult))
		statfileoutput.close()
		shutil.move(self.statmovefilepath,self.statfilepath)