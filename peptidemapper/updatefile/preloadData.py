#!/usr/bin/env python
# -*- coding: utf-8 -*-
import unicodedata
import datetime,glob
import time
import os,subprocess,psutil,re,sys,shutil,datetime
import csv
import datetime
import urllib,urllib2,urllib3
import json
import numpy as np
import ctypes
from itertools import combinations
import operator

def preLoadDataJson(homedir,filepath):
	csv.field_size_limit(sys.maxsize)
	calfilename='calculationprog.py'
	calmovefilepath=os.path.join(homedir, 'updatefile', calfilename)
	calfilepath = os.path.join(homedir, 'src/pepmapperapp', calfilename)
	speciesdic={}
	pepseqdic={}
	prodic={}
	speciesProt={}
	speciesPep={}
	unquniprotdic={}
	isodic={}
	godic={}
	subcell={}
	totalpeptideseq=''
	aacountlist=[]
	aalist=['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
	finalreader = csv.DictReader(open(filepath),delimiter='\t')
	for frow in finalreader:
		pepseq=frow['Peptide Sequence'].strip()
		Calorg=str(frow['Organism']).strip()
		CalorgID=str(frow['Organism ID']).strip()
		speciesdic[CalorgID]=Calorg
		acccode=None
		acccode=str(frow['UniProtKB Accession']).split('-')[0]
		if acccode !=None and str(frow['UniprotKb entry status']).strip().upper()=='YES':
			if (CalorgID).lower()!='na' and len(CalorgID) >0:
				totalpeptideseq +=pepseq.upper()
				if speciesProt.has_key(CalorgID):
					speciesProt[CalorgID].append(acccode)
				else:
					speciesProt[CalorgID]=[acccode]

				if speciesPep.has_key(CalorgID):
					speciesPep[CalorgID].append(pepseq)
				else:
					speciesPep[CalorgID]=[pepseq]

			if frow["Go Name"].upper() !='NA' and len(str(frow["Go Name"]).strip())>0:
				goname=(str(frow["Go Name"]).strip()).split('|')
				for gitem in goname:
					if godic.has_key(str(gitem).strip()):
						godic[str(gitem).strip()].append(acccode)
					else:
						godic[str(gitem).strip()]=[acccode]

			if frow["SubCellular"].upper() !='NA' and len(str(frow["SubCellular"]).strip())>0:
				subcelllist=str(frow["SubCellular"]).strip().split('|')
				for subcitem in subcelllist:
					if subcell.has_key(str(subcitem).strip()):
						subcell[str(subcitem).strip()].append(acccode)
					else:
						subcell[str(subcitem).strip()]=[acccode]

			if (frow['PeptideTracker ID']).lower()!='na' and len(str(frow['PeptideTracker ID']).strip()) >0:
				if pepseqdic.has_key(pepseq):
					pepseqdic[pepseq].append('PeptideTracker')
				else:
					pepseqdic[pepseq]=['PeptideTracker']
				if prodic.has_key(acccode):
					prodic[acccode].append('PeptideTracker')
				else:
					prodic[acccode]=['PeptideTracker']

			if (frow['Passel ID']).lower()!='na' and len(str(frow['Passel ID']).strip()) >0:
				if pepseqdic.has_key(pepseq):
					pepseqdic[pepseq].append('PASSEL')
				else:
					pepseqdic[pepseq]=['PASSEL']
				if prodic.has_key(acccode):
					prodic[acccode].append('PASSEL')
				else:
					prodic[acccode]=['PASSEL']


			if (frow['SRMAtlas ID']).lower()!='na' and len(str(frow['SRMAtlas ID']).strip()) >0:
				if pepseqdic.has_key(pepseq):
					pepseqdic[pepseq].append('SRMAtlas')
				else:
					pepseqdic[pepseq]=['SRMAtlas']
				if prodic.has_key(acccode):
					prodic[acccode].append('SRMAtlas')
				else:
					prodic[acccode]=['SRMAtlas']


			if (frow['Cptac ID']).lower()!='na' and len(str(frow['Cptac ID']).strip()) >0:
				if pepseqdic.has_key(pepseq):
					pepseqdic[pepseq].append('CPTAC')
				else:
					pepseqdic[pepseq]=['CPTAC']
				if prodic.has_key(acccode):
					prodic[acccode].append('CPTAC')
				else:
					prodic[acccode]=['CPTAC']

			if (frow['Panoramaweb ID']).lower() !='na' and len(str(frow['Panoramaweb ID']).strip()) >0:
				if pepseqdic.has_key(pepseq):
					pepseqdic[pepseq].append('PanoramaWeb')
				else:
					pepseqdic[pepseq]=['PanoramaWeb']
				if prodic.has_key(acccode):
					prodic[acccode].append('PanoramaWeb')
				else:
					prodic[acccode]=['PanoramaWeb']


	speciesProt={k:list(set(j)) for k,j in speciesProt.items()}
	speciesPep={k:list(set(j)) for k,j in speciesPep.items()}
	modspeciesProt={}
	modspeciesPep={}
	for sprok in speciesProt:
		if sprok in speciesdic:
			modspeciesProt[speciesdic[sprok]]=speciesProt[sprok]

	for spepk in speciesPep:
		if spepk in speciesdic:
			modspeciesPep[speciesdic[spepk]]=speciesPep[spepk]

	pepseqdic={k:list(set(j)) for k,j in pepseqdic.items()}
	prodic={k:list(set(j)) for k,j in prodic.items()}
	subcell={k:list(set(j)) for k,j in subcell.items()}

	list_of_pepdb=pepseqdic.values()
	flattenedpepdb  = [val for sublist in list_of_pepdb for val in sublist]
	uniqueflattenedpepdb=list(set(flattenedpepdb))
	pepseqdataseries={}
	for pepitem in uniqueflattenedpepdb:
		subdata=[]
		for pepkey in pepseqdic:
			if pepitem in pepseqdic[pepkey]:
				subdata.append(pepkey)
		pepseqdataseries[pepitem]=subdata

	prodataseries={}
	for proitem in uniqueflattenedpepdb:
		subdata=[]
		for prokey in prodic:
			if proitem in prodic[prokey]:
				subdata.append(prokey)
		prodataseries[proitem]=subdata

	for aa in aalist:
		aacountlist.append([aa,round((((float(totalpeptideseq.count(aa)))/(float(len(totalpeptideseq))))*100),2)])

	mrmdatabase={'CPTAC':'A','PanoramaWeb':'B','PeptideTracker':'C','SRMAtlas':'D','PASSEL':'E'}
	prodatavalues={}
	pepdatavalues={}
	modmrmdatabase={mrmdatabase[k]:k for k in prodataseries.keys()}

	for i in range(len(prodataseries)):
		for v in combinations(prodataseries.keys(),i+1):
			vsets = [set(prodataseries[x]) for x in v ]
			sets=tuple(sorted(v))
			sets= list(sets)
			mrmdbname=[mrmdatabase[k] for k in sets]
			mrmdbname.sort()
			label=list(reduce(lambda x,y: x.intersection(y), vsets))
			prodatavalues[ ''.join(mrmdbname)]=label

	for j in range(len(pepseqdataseries)):
		for v in combinations(pepseqdataseries.keys(),j+1):
			vsets = [set(pepseqdataseries[x]) for x in v ]
			sets=tuple(sorted(v))
			sets= list(sets)
			mrmdbname=[mrmdatabase[k] for k in sets]
			mrmdbname.sort()
			label=list(reduce(lambda x,y: x.intersection(y), vsets))
			pepdatavalues[ ''.join(mrmdbname)]=label

	valueskey=prodatavalues.keys()
	sortedvalueskey=sorted(valueskey, key=len,reverse=True)
	newprodatavalues={}
	newpepdatavalues={}
	for vk in sortedvalueskey:
		try:
			tempvalprot=prodatavalues[vk]
			newprodatavalues[vk]=len(tempvalprot)
			tempvalpep=pepdatavalues[vk]
			newpepdatavalues[vk]=len(tempvalpep)
			prodatavalues.pop(vk, None)
			pepdatavalues.pop(vk, None)
			tmpprodatavalues = dict((k, list(set(v)-set(tempvalprot))) for k, v in prodatavalues.items())
			prodatavalues.update(tmpprodatavalues)

			tmppepdatavalues = dict((k, list(set(v)-set(tempvalpep))) for k, v in pepdatavalues.items())
			pepdatavalues.update(tmppepdatavalues)
		except KeyError:
			pass


	finalresult={}
	finalresult['pepseqdic']=pepseqdic
	finalresult['prodic']=prodic
	finalresult['pepseqdataseries']=pepseqdataseries
	finalresult['prodataseries']=prodataseries
	finalresult['prodatavalues']=newprodatavalues
	finalresult['pepdatavalues']=newpepdatavalues
	finalresult['mrmdatabase']=modmrmdatabase
	finalresult['speciesProt']=modspeciesProt
	finalresult['speciesPep']=modspeciesPep
	finalresult['aacountlist']=aacountlist
	finalresult['godic']={k:len(set(v)) for k, v in godic.items()}
	finalresult['subcell']={k:len(v) for k, v in subcell.items()}
	finalresult['subcell']=sorted(finalresult['subcell'].items(), key=operator.itemgetter(1))
	finalresult['subcell']=map(list,finalresult['subcell'])
	finalresult['godic']=sorted(finalresult['godic'].items(), key=operator.itemgetter(1))
	finalresult['godic']=map(list,finalresult['godic'])

	calfileoutput=open(calfilename,'w')
	calfileoutput.write("finalresult=")
	calfileoutput.write(json.dumps(finalresult))
	calfileoutput.close()
	speciesdic.clear()
	pepseqdic.clear()
	prodic.clear()
	speciesProt.clear()
	speciesPep.clear()
	modspeciesProt.clear()
	modspeciesPep.clear()
	unquniprotdic.clear()
	isodic.clear()
	godic.clear()
	subcell.clear()
	finalresult.clear()
	newprodatavalues.clear()
	newpepdatavalues.clear()
	prodatavalues.clear()
	pepdatavalues.clear()
	shutil.move(calmovefilepath,calfilepath)
