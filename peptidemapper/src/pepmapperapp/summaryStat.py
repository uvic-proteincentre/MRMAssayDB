from collections import Counter
from itertools import combinations
import ast
from operator import itemgetter
from summaryPathwayCoverage import *
import operator
import json
import re,sys
import itertools
from collections import OrderedDict
import datetime
import gc
from jvenn import *
# function to generate stat based on user search result
def summaryStatcal(pepresult):
	specieslist=[]
	unqisostat=[]
	subcell={}
	godic={}
	statPathwayDic={}
	prodataseries={}
	pepseqdataseries={}
	databaseNameList=['PeptideTracker ID','Passel ID','SRMAtlas ID','Cptac ID','Panoramaweb ID']
	pattern = re.compile(r"-\d+")

	listOfSubcell=list((object['SubCellular'] for object in pepresult)) #reading subcellular data from result
	listOfSubcell=filter(None, listOfSubcell) #filter out empty value
	if len(listOfSubcell)>0:
		listOfSubcell=filter(lambda k: 'NA' != k.upper(), listOfSubcell) #filter out NA value
		if len(listOfSubcell)>0:
			subcellcount=Counter((('|'.join(listOfSubcell)).split('|'))) #split subcellular location then count the occurance
			tempsortedsubcell=OrderedDict(sorted(dict(subcellcount).items(), key=lambda x: x[1],reverse=True)) #sorting dictionary based on number of unique items in dict and reverse order
			listOfSubcell=(dict(itertools.islice(tempsortedsubcell.items(), 0, 50))).keys() # only top 50 restult will be displayed

	listOfSpecies=list((object['Organism'] for object in pepresult)) #reading organism data from result
	listOfSpecies=list(set(filter(None, listOfSpecies))) #filter out empty value

	listOfGoName=list((object['Go Name'] for object in pepresult)) #reading GO name data from result
	listOfGoName=filter(None, listOfGoName) #filter out empty value
	if len(listOfGoName)>0:
		listOfGoName=filter(lambda k: 'NA' != k.upper(), listOfGoName) #filter out NA value
		if len(listOfGoName)>0:
			listOfGoNamecount=Counter((('|'.join(listOfGoName)).split('|'))) #split data then count the occurance
			tempsortedgo=OrderedDict(sorted(dict(listOfGoNamecount).items(), key=lambda x: x[1],reverse=True)) #sorting dictionary based on number of unique items in dict and reverse order
			listOfGoName=(dict(itertools.islice(tempsortedgo.items(), 0, 50))).keys() # only top 50 restult will be displayed

	listOfPathwayName=list((object['Pathway Name'] for object in pepresult)) #reading Pathway data from result
	listOfPathwayName=filter(None, listOfPathwayName) #filter out empty value
	if len(listOfPathwayName)>0:
		listOfPathwayName=filter(lambda k: 'NA' != k.upper(), listOfPathwayName) #filter out NA value
		if len(listOfPathwayName)>0:
			listOfPathwayNamecount=Counter((('|'.join(listOfPathwayName)).split('|'))) #split data then count the occurance
			tempsortedPathway=OrderedDict(sorted(dict(listOfPathwayNamecount).items(), key=lambda x: x[1],reverse=True)) #sorting dictionary based on number of unique items in dict and reverse order
			listOfPathwayName=(dict(itertools.islice(tempsortedPathway.items(), 0, 20))).keys()# only top 20 restult will be displayed

	for i in listOfSubcell:
		subdatafilter=filter(lambda pepdata: (str(i).strip()).lower() in (str(pepdata['SubCellular']).strip()).lower(), pepresult)
		subdata=list((object['UniProtKB Accession'] for object in subdatafilter))
		subcell[str(i).strip()]=len(list(set([pattern.sub("", item) for item in subdata])))

	for i in listOfSpecies:
		subdatafilter=filter(lambda pepdata: str(i).strip().lower() in (str(pepdata["Organism"]).strip()).lower(), pepresult)
		subdataprot=list((object['UniProtKB Accession'] for object in subdatafilter))
		subdatapep=list((object['Peptide Sequence'] for object in subdatafilter))
		subdataprot=list(set([pattern.sub("", item) for item in subdataprot]))
		subdatapep=list(set([pattern.sub("", item) for item in subdatapep]))
		specieslist.append([i,len(subdataprot),len(subdatapep)])

	for i in listOfPathwayName:
		subdatafilter=filter(lambda pepdata: str(i).strip().lower() in (str(pepdata["Pathway Name"]).strip()).lower(), pepresult)
		subdata=list((object['UniProtKB Accession'] for object in subdatafilter))
		statPathwayDic[str(i).strip()]=list(set([pattern.sub("", item) for item in subdata]))

	for i in listOfGoName:
		subdatafilter=filter(lambda pepdata: str(i).strip().lower() in (str(pepdata["Go Name"]).strip()).lower(), pepresult)
		subdata=list((object['UniProtKB Accession'] for object in subdatafilter))
		godic[str(i).strip()]=len(list(set([pattern.sub("", item) for item in subdata])))

	for i in databaseNameList:
		subdatafilter=filter(lambda pepdata: ((str(pepdata[i]).strip()).lower() !='na' and len(str(pepdata[i]).strip()) >0), pepresult)
		subdataprot=list((object['UniProtKB Accession'] for object in subdatafilter))
		subdatapep=list((object['Peptide Sequence'] for object in subdatafilter))
		prodataseries[(str(i).strip()).split(' ')[0]]=list(set([pattern.sub("", item) for item in subdataprot]))
		pepseqdataseries[(str(i).strip()).split(' ')[0]]=list(set([pattern.sub("", item) for item in subdatapep]))

		isostatus=list((object['Present in isoforms'] for object in subdatafilter))

		unqstatus=list((object['Unique in protein'] for object in subdatafilter))
		tempcalunq='0.0%'
		tempcaliso='0.0%'
		if len(isostatus)>0:
			isostatus = map(lambda x: 'NA' if x.strip() == '' else x, isostatus)
			tempcaliso=str(100-(round((100.00-(float(isostatus.count('NA'))/float(len(isostatus)))*100),2)))+'%'
		if len(unqstatus)>0:
			tempcalunq=str(round(((float(unqstatus.count('Yes'))/float(len(unqstatus)))*100),2))+'%'
		unqisostat.append([(str(i).strip()).split(' ')[0],tempcalunq,tempcaliso])

	sortedPathwayStatdic=OrderedDict(sorted(statPathwayDic.items(), key=lambda x: len(x[1]), reverse=True))
	top10Pathwaydict=dict(itertools.islice(sortedPathwayStatdic.items(), 0, 10))
	pathwaychart=summaryPathwaycal(top10Pathwaydict,prodataseries)
	prodatavalues,pepdatavalues,modmrmdatabase=jvennstat(prodataseries,pepseqdataseries)
	statfinalresult={}
	statfinalresult['total']=[len(set(reduce(operator.concat, prodataseries.values()))),len(set(reduce(operator.concat, pepseqdataseries.values())))]
	statfinalresult['pepseqdataseries']=pepseqdataseries
	statfinalresult['prodataseries']=prodataseries
	statfinalresult['specieslist']=specieslist
	statfinalresult['unqisostat']=unqisostat
	statfinalresult['pathwaychart']=pathwaychart
	sortedGodic=OrderedDict(sorted(godic.items(), key=lambda x: x[1],reverse=True))

	sortedSubcelldic=OrderedDict(sorted(subcell.items(), key=lambda x: x[1],reverse=True))

	statfinalresult['jevennstat']=[prodatavalues,pepdatavalues,modmrmdatabase]
	statfinalresult['subcell']=sortedSubcelldic
	statfinalresult['godic']=sortedGodic
	gc.collect()
	return statfinalresult