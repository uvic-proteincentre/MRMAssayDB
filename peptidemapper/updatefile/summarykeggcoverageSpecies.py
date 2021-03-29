import os,sys
from operator import itemgetter
import numpy as np
from itertools import count, groupby
import pandas as pd
import csv
import ast
from operator import itemgetter

def summarykeggcal(searchtype,matchhit,pepfinalresult):
	homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
	filename='ReportBook_mother_file100V3.csv'
	filepath = os.path.join(homedir, 'src/mappermotherfile', filename)
	keggpathwaycoverage=[]
	statKeggDic={}
	for searchterm in matchhit:
		with open(filepath) as repcsvfile:
			pepresult = csv.DictReader(repcsvfile, delimiter='\t')
			for reppeprow in pepresult:
				statuniID=str(reppeprow['UniProtKB Accession']).strip()
				statPathway=str(reppeprow['Kegg Pathway Name']).strip()
				statSpecies=str(reppeprow['Organism']).strip()
				if len(statPathway)>0:
					statPathwayList=statPathway.split('|')
					for stpItem in statPathwayList:
						if stpItem.lower() !='na':
							keggid=stpItem.strip()+'|'+statSpecies
							if statKeggDic.has_key(keggid):
								statKeggDic[keggid].append(statuniID.strip())
							else:
								statKeggDic[keggid] =[statuniID.strip()]

		#kegg pathway name extraction

		#kegg pathway name extraction

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
			keggpathwaycoverage.append(templist)

		unqkeggpathwaycoverage=[list(tupl) for tupl in {tuple(item) for item in keggpathwaycoverage }]
		keggchart=[]
		if len(unqkeggpathwaycoverage)>0:
			sortedkeggpathwaycoverage=sorted(unqkeggpathwaycoverage, key= itemgetter(1), reverse=True)
			keggchart=[['Pathway Name', 'Total', 'PeptideTracker','CPTAC', 'PASSEL','SRMAtlas', 'PanoramaWeb']]
			keggchart=keggchart+sortedkeggpathwaycoverage
		return keggchart