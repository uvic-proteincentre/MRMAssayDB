import os,sys
from operator import itemgetter
import numpy as np
from itertools import count, groupby
import pandas as pd
import ast
from operator import itemgetter
import gc
def summaryPathwaycal(top10Pathwaydict,prodataseries):
	pathwaycoverage=[]

	for pskey in top10Pathwaydict:
		pathwayname=(pskey.strip()).split('|')[0]
		tempUniqPathwayUniIDList=list(set(top10Pathwaydict[pskey]))
		cptac=[]
		panweb=[]
		passel=[]
		srmatlas=[]
		peptrack=[]
		for ckey in prodataseries:
			if "peptidetracker" in ckey.lower():
				peptrack=list(set(prodataseries[ckey]).intersection(tempUniqPathwayUniIDList))
			if "passel" in ckey.lower():
				passel=list(set(prodataseries[ckey]).intersection(tempUniqPathwayUniIDList))
			if "srmatlas" in ckey.lower():
				srmatlas=list(set(prodataseries[ckey]).intersection(tempUniqPathwayUniIDList))
			if "cptac" in ckey.lower():
				cptac=list(set(prodataseries[ckey]).intersection(tempUniqPathwayUniIDList))
			if "panoramaweb" in ckey.lower():
				panweb=list(set(prodataseries[ckey]).intersection(tempUniqPathwayUniIDList))
		tempcptac=len(list(set(cptac)))
		temppanweb=len(list(set(panweb)))
		temppassel=len(list(set(passel)))
		temppeptrack=len(list(set(peptrack)))
		tempsrmatlas=len(list(set(srmatlas)))
		tempTotal=len(list(set(tempUniqPathwayUniIDList)))
		templist=[pathwayname,tempTotal,temppeptrack,tempcptac,temppassel,tempsrmatlas,temppanweb]
		pathwaycoverage.append(templist)

	unqpathwaycoverage=[list(tupl) for tupl in {tuple(item) for item in pathwaycoverage }]
	pathwaychart=[]
	if len(unqpathwaycoverage)>0:
		sortedpathwaycoverage=sorted(unqpathwaycoverage, key= itemgetter(1), reverse=True)
		pathwaychart=[['Pathway Name', 'Total unique proteins', 'PeptideTracker','CPTAC', 'PASSEL','SRMAtlas', 'PanoramaWeb']]
		pathwaychart=pathwaychart+sortedpathwaycoverage
	gc.collect()
	return pathwaychart