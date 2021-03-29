#!/usr/bin/env.python
# -*- coding: utf-8 -*-
# encoding: utf-8
from __future__ import unicode_literals
from django.utils.encoding import force_text
from django.shortcuts import render
from django.http import HttpResponse,HttpResponseRedirect
import sys,re,os,glob,shutil,subprocess,socket
import urllib,urllib3
import fileinput
from .models import IpAddressInformation
from ipware.ip import get_ip
from time import gmtime, strftime,sleep
import csv
import hashlib, random
import datetime
from django.conf import settings
from django.core.mail import send_mail
from django.contrib import messages
# from django.core.context_processors import csrf
from django.core.mail import EmailMultiAlternatives
from django.template import RequestContext
import json
import json as simplejson
import calendar
from django.contrib import auth
from bioservices.kegg import KEGG
from xml.etree import cElementTree as ET
import xmltodict
from xml.dom import minidom
from xml.parsers.expat import ExpatError
from Bio.PDB.Polypeptide import *
from .forms import ContactForm
from Bio import SeqIO
from requests.exceptions import ConnectionError
import requests
from django.utils.datastructures import MultiValueDictKeyError
from Bio.SeqUtils import seq1
from goatools import obo_parser
from elasticsearch import Elasticsearch,helpers,RequestsHttpConnection
import pandas as pd
from .colName import *
from summaryStat import summaryStatcal
from .calculationprog import *
from .totalpepassay import *
from .overallstat import *
from updatedstat import *
import random
import names
from operator import itemgetter
import ast
from collections import OrderedDict
from itertools import combinations
import pandas as pd
import operator
from pyteomics import parser
# Import Biopython modules to interact with KEGG
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
#from Bio.Graphics.KGML_vis import KGMLCanvas,get_temp_imagefilename
#from Bio.Graphics.ColorSpiral import ColorSpiral
# pdf to image
#from pdf2image import convert_from_path
#from PIL import Image
#delete after publication
from django.contrib.auth.decorators import login_required
from django.contrib.auth.views import login,logout

searchFields=["UniProtKB Accession.ngram","Protein.ngram","Gene.ngram",\
"Organism.ngram","Organism ID.ngram","SubCellular.ngram","Peptide Sequence.ngram",\
"Pathway Name.ngram","Disease Name.ngram",\
"Go ID.ngram","Go Name.ngram","Go Term.ngram","Drug Bank.ngram",\
"Assays for FDA approved Marker.ngram"]

es = Elasticsearch(
	['http://172.16.1.60:9200/'],
	connection_class=RequestsHttpConnection
)
# Create your views here.

def search_form(request):
	"""
	This is basic search page for peptidemapper.
	It also store information for user ip address
	""" 
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	title= "Please search by:"
	# opening files for plotting stat
	organismName=overallSumresult['organism']
	speciesName=overallSumresult['species']
	speciesstat=overallSumresult['speciesstat'][0:10]
	speciesName=list(set(speciesName))
	speciesName=sorted(speciesName)
	speciesstat.insert(0,['Species','Unique protein','Unique peptide'])
	gostat=overallSumresult['gostat'][:10]
	gostat.insert(0,['Go Term','Unique proteins in various species'])
	keggstat=overallSumresult['keggstat'][:10]
	keggstat.insert(0,['Pathway Name', 'Unique proteins in various species', 'PeptideTracker', 'CPTAC', 'PASSEL', 'SRMAtlas', 'PanoramaWeb'])
	pepseqdic=finalresult['pepseqdic']
	prodic=finalresult['prodic']
	pepdatavalues=finalresult['pepdatavalues']
	prodatavalues=finalresult['prodatavalues']
	mrmdatabase=finalresult['mrmdatabase']
	allpepassay=totalpepassay['totalassayNonValid']
	allvalidpepassay=totalpepassay['totalassayValid']
	allunqStripPep=totalpepassay['totalstripPep']
	uqpep=len(pepseqdic)
	uqpro=len(prodic)
	keggstat=[i[:2] for i in keggstat]
	speciesstat=[i[:2] for i in speciesstat]
	contextindex ={"title": title,"uqpro":uqpro, "uqpep":uqpep,\
					"speciesName":speciesName,"speciesnumber":len(speciesName)-1,\
					"speciesstat":json.dumps(speciesstat),\
					"gostat":json.dumps(gostat),"keggstat":json.dumps(keggstat),\
					'allpepassay':allpepassay,\
					'allvalidpepassay':allvalidpepassay,\
					'allunqStripPep':len(allunqStripPep),\
					'jvennpep':json.dumps(pepdatavalues),\
					'jvennprot':json.dumps(prodatavalues),\
					'jvennmrmdb':json.dumps(mrmdatabase)\
					}
	return render(request, 'index.html', contextindex)


def search(request):
	"""
	This is basic search function, based on given search parameters it will generate result datatable along with stat
	from database.
	"""
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if request.method=='POST':
		searchterm=request.POST.getlist('searchterm')# user input for searching result
		searchterm=map(str, searchterm)
		searchterm=searchterm[0]
		searchterm= searchterm.strip()

		if len(searchterm)>0:
			contextres =[]
			#build elasticsearch query to search data
			query={
				"query":{
					"bool":{
						"should":[{
								"multi_match":{
									"query":searchterm,
									"type":"best_fields",
									"fields":searchFields,
									"minimum_should_match":"100%"
								}
						}]
					}
				}
			}
			#generate random file name to store search result in json format
			nameFIle=names.get_first_name()
			jsonfilename=nameFIle+'_basic_search.json'
			jsonfilepath=os.path.join(settings.BASE_DIR, 'resultFile', 'jsonData','resultJson', 'basicsearch', 'results', jsonfilename)
			jsonfileoutput= open(jsonfilepath,'w')
			jfinaldata=[]
			es.indices.refresh(index="mrmassaydb-index")
			#elasticsearch will search data
			res=helpers.scan(client=es,scroll='2m',index="mrmassaydb-index", doc_type="mrmassaydb-type",query=query,request_timeout=30)
			jfinaldata=[]
			#if data is valid based on uniprotkb release then it will display
			for i in res:
				jdic=i['_source']
				jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
				if jdic["UniprotKb entry status"] =="Yes" and jdic['UniProtKB Accession'] !='502':
					jdic["PPI"] ="View"
					jdic["sel"] =""
					jdic["Drug Bank"]=jdic["Drug Bank"].replace('\\','')
					jdic["Drug Bank"]=jdic["Drug Bank"].replace('<br>','|')
					jdic["SRMAtlas URL"]=jdic["SRMAtlas URL"].replace('\\','')
					jdic["Passel URL"]=jdic["Passel URL"].replace('\\','')
					jdic["CPTAC URL"]=jdic["CPTAC URL"].replace('\\','')
					jdic["Panoramaweb URL"]=jdic["Panoramaweb URL"].replace('\\','')
					jdic["PeptideTracker URL"]=jdic["PeptideTracker URL"].replace('\\','')
					#if jdic["Pathway Name"].lower() !='na':
					#jdic["Pathway Name"]=re.sub(r"(\w)([A-Z])",r"\1|\2",jdic["Pathway Name"])
					jfinaldata.append(jdic)
			es.indices.refresh(index="mrmassaydb-index")
			#checking any result generated by database
			foundHits=len(jfinaldata)
			#storing only 10000 rows in json format
			json.dump(jfinaldata[:10000],jsonfileoutput)
			jsonfileoutput.close()
			# if result found then do other job
			if foundHits >0:
				statsummary=summaryStatcal(jfinaldata) # sent data to this funcation for generating stat
				pathwaychart=statsummary['pathwaychart']
				pathwaychart=[i[:2] for i in pathwaychart]
				specieslist=statsummary['specieslist']
				totallist=statsummary['total']
				subcell=statsummary['subcell']
				godic=statsummary['godic']
				jvennprot=statsummary['jevennstat'][0]
				jvennpep=statsummary['jevennstat'][1]
				mrmdatabase=statsummary['jevennstat'][2]
				sortedgodic=OrderedDict(sorted(godic.items(), key=lambda t: t[1])) # sorting GO data
				updatedgodic=dict(list(sortedgodic.items()))
				pepseqdataseries=ast.literal_eval(json.dumps(statsummary['pepseqdataseries'])) #dumping data into json format
				prodataseries=statsummary['prodataseries']
				unqisostat=statsummary['unqisostat']
				jsonfilepathStat=os.path.join(settings.BASE_DIR, 'resultFile', 'jsonData','resultJson', 'basicsearch', 'statsummary', jsonfilename) #storing stat result in json format
				jsonfileoutputStat= open(jsonfilepathStat,'w')
				json.dump(statsummary,jsonfileoutputStat)
				jsonfileoutputStat.close()
				urlname="'/resultFile/jsonData/resultJson/basicsearch/results/"+jsonfilename+"'"

				contextindex={
					"filename":urlname,"colname":json.dumps(colname),
					'query': searchterm,'foundHits':foundHits,
					'pathwaychart':pathwaychart[:11],'specieslist':specieslist,
					'totallist':totallist,'subcell':subcell,
					'updatedgodic':updatedgodic,'pepseqdataseries':pepseqdataseries,
					'prodataseries':prodataseries,'unqisostat':unqisostat,
					'jvennprot':json.dumps(jvennprot),'jvennpep':json.dumps(jvennpep),'jvennmrmdb':json.dumps(mrmdatabase)
					}
				return render(request,'resultform.html',contextindex)
			else:
				return render(request,'resultform.html',{'foundHits':foundHits})
		else:
			return render(request,'resultform.html',{'foundHits':0})

def advanced_search(request):
	"""
	This is advance search function, based on given search parameters it will generate result datatable along with stat
	from database.
	"""
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if request.method=='POST':
		searchterm =request.POST.getlist('searchterm') # list of search term value associated with searchtype
		searchtype =request.POST.getlist('searchtype') # list of search parameter
		searchtermorg=request.POST.getlist('searchtermorg') # search term value for organism
		searchtermfda=request.POST.getlist('searchtermfda') # search term value for FDA
		searchtermlist=[]
		nameFIle=names.get_first_name() # generate random file name to store user search result
		fastaseq=[]
		finalsearhdata=''
		unique_peptides = set()
		tryptic_peptide={}
		userSuppliedPepSeqStatus=0
		try:
			fastafile = request.FILES["fileupload"].read()
			finalsearhdata+='File'+':'+'Fasta Sequence'+' '
			currdate=str(datetime.datetime.now())
			currdate=currdate.replace('-','_')
			currdate=currdate.replace(' ','_')
			currdate=currdate.replace(':','_')
			currdate=currdate.replace('.','_')
			nameFIle=currdate+'_'+str(request.FILES["fileupload"]).split('.')[0] # if user upload fasta file then file name will be replaced with user provided file name along with current data and time
			fastafilename=nameFIle+'.fasta'
			#storing user provided fasta file
			fastafilepath=os.path.join(settings.BASE_DIR, 'resultFile', 'fastaFIle', fastafilename)
			fastafilewrite=open(fastafilepath,"w")
			fastafilewrite.write(fastafile)
			fastafilewrite.close()

			#reading fasta file
			seqCounter=0
			for useq_record in SeqIO.parse(fastafilepath, 'fasta'):
				seqCounter+=1
				seqheader = useq_record.id
				sequniID = seqheader.split(' ')[0]
				sequniID=sequniID.replace('>','')
				tempseqs = str(useq_record.seq).strip()
				new_peptides = parser.cleave(tempseqs, 'trypsin')
				new_peptides=[pep for pep in new_peptides if len(pep.strip()) > 3 and len(pep.strip()) <50]
				tryptic_peptide[seqCounter]=list(new_peptides)
				new_peptides=list(set(new_peptides))
				unique_peptides.update(new_peptides)
				fastaseq.append(str(sequniID)+'_'+tempseqs.upper())
		except MultiValueDictKeyError:
			pass

		try:
			fastafileindex=searchtype.index("FastaFile")
			#delete data based on index from list
			del searchtype[fastafileindex]
			del searchterm[fastafileindex]
		except ValueError:
			pass

		try:
			orgindex=searchtype.index("Organism")
			#delete data based on index from list
			del searchtype[orgindex]
			del searchterm[orgindex]
		except ValueError:
			pass
		if len(fastaseq)>0:
			unique_peptides=list(unique_peptides)
			unique_peptides=list(map(lambda x:x.lower(),unique_peptides))
		searchtermorg=map(str, searchtermorg) # convert data into string
		searchtermorg=map(lambda j: j.strip(), searchtermorg) # remove space
		searchtermorg=filter(None, searchtermorg) # remove empty value
		unqsearchtermorg=list(set(searchtermorg))
		if len(unqsearchtermorg)>0:
			finalsearhdata+='Organism'+':'+unqsearchtermorg[0].strip()+' '
			#build elasticsearch query for organism to search data

			orgquery={"should":[
							{
								"multi_match":{
									"query":unqsearchtermorg[0].strip(),
									"type":"best_fields",
									"fields":["Organism.ngram"],
									"minimum_should_match":"100%"
								}
							}
						]
					}
			booldic={}
			booldic["bool"]=orgquery
			searchtermlist.append(booldic)

		try:
			fdaindex=searchtype.index("Assays for FDA approved Marker")
			#delete data based on index from list
			del searchtype[fdaindex]
			del searchterm[fdaindex]
		except ValueError:
			pass

		searchtermfda=map(str, searchtermfda) # convert data into string
		searchtermfda=map(lambda j: j.strip(), searchtermfda) # remove space
		searchtermfda=filter(None, searchtermfda) # remove empty value
		unqsearchtermfda=list(set(searchtermfda))
		if len(unqsearchtermfda)>0:
			finalsearhdata+='Assays for FDA approved Marker'+':'+unqsearchtermfda[0].strip()+' '
			#build elasticsearch query for FDA to search data

			fdaquery={"should":[
							{
								"multi_match":{
									"query":unqsearchtermfda[0].strip(),
									"type":"best_fields",
									"fields":["Assays for FDA approved Marker.ngram"],
									"minimum_should_match":"100%"
								}
							}
						]
					}
			booldic={}
			booldic["bool"]=fdaquery
			searchtermlist.append(booldic)
		if 'Peptide Sequence' in searchtype:
			userSuppliedPepSeqStatus=1
		for i in range(0,len(searchtype)):
			subsearchtype=searchtype[i]
			subsearchterm=searchterm[i]
			#build elasticsearch query for all except organism and FDA to search data
			if '|' in subsearchterm:
				subsearchterm=(subsearchterm.strip()).split('|')
			else:
				subsearchterm=(subsearchterm.strip()).split('\n')
			subsearchterm=map(str, subsearchterm)
			subsearchterm=map(lambda j: j.strip(), subsearchterm)
			subsearchterm=filter(None, subsearchterm)
			if subsearchtype == 'Peptide Sequence':
				if userSuppliedPepSeqStatus==1:
					finalsearhdata+=''.join(subsearchtype)+':'+';'.join(subsearchterm)+' '
					if len(unique_peptides)>0:
						subsearchterm=[(item.strip()).lower() for item in subsearchterm]
						subsearchterm=list(set(subsearchterm) & set(unique_peptides))
			else:
				finalsearhdata+=''.join(subsearchtype)+':'+';'.join(subsearchterm)+' '
			if len(subsearchterm)>0:
				subsearchterm=[(item.strip()).lower() for item in subsearchterm] #converting into lower case
				subsearchterm=list(set(subsearchterm))
				shouldlist=[]
				
				for x in subsearchterm:
					tempquery={
								"multi_match":{
									"query":x.strip(),
									"type":"best_fields",
									"fields":[str(subsearchtype)+".ngram"],
									"minimum_should_match":"100%"
								}
							}
					shouldlist.append(tempquery)
				booldic={}
				booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
				searchtermlist.append(booldic)

		if userSuppliedPepSeqStatus==0 and len(unique_peptides)>0:
			shouldlist=[]
			for x in unique_peptides:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Peptide Sequence.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		unqfastaseq=list(set(fastaseq))
		if len(searchtermlist)>0 or len(unqfastaseq)>0:
			es.indices.refresh(index="mrmassaydb-index")

			query=""
			#if len(searchtermlist)>0:
			query={
				"query": {
					"bool": {
						"must":searchtermlist
					}
				}
			}
			# if len(searchtermlist)==0:
			# 	query={
			# 		"query": {
			# 			"match_all": {}
			# 		}
			# 	}
			#storing user search result into json format
			jsonfilename=nameFIle+'_advance_search.json'
			jsonfilepath=os.path.join(settings.BASE_DIR, 'resultFile', 'jsonData','resultJson', 'adavancesearch', 'results', jsonfilename)
			jsonfileoutput= open(jsonfilepath,'w')
			jfinaldata=[]
			res=helpers.scan(client=es,size=1000,scroll='2m',index="mrmassaydb-index", doc_type="mrmassaydb-type",query=query,request_timeout=60)
			#res=helpers.scan(client=es,size=1000,scroll='2m',index="my-index", doc_type="my-type",query=query,request_timeout=30)
			jfinaldata=[]
			usersequnq=[]
			for i in res:
				jdic=i['_source']
				jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
				if jdic["UniprotKb entry status"] =="Yes" and jdic['UniProtKB Accession'] !='502':
					jdic["PPI"] ="View"
					jdic["sel"] =""
					jdic["Drug Bank"]=jdic["Drug Bank"].replace('\\','')
					jdic["Drug Bank"]=jdic["Drug Bank"].replace('<br>','|')
					jdic["SRMAtlas URL"]=jdic["SRMAtlas URL"].replace('\\','')
					jdic["Passel URL"]=jdic["Passel URL"].replace('\\','')
					jdic["CPTAC URL"]=jdic["CPTAC URL"].replace('\\','')
					jdic["Panoramaweb URL"]=jdic["Panoramaweb URL"].replace('\\','')
					jdic["PeptideTracker URL"]=jdic["PeptideTracker URL"].replace('\\','')
					#if jdic["Pathway Name"].lower() !='na':
					#	jdic["Pathway Name"]=re.sub(r"(\w)([A-Z])",r"\1|\2",jdic["Pathway Name"])
					seqhit=0
					# checking any peptide present in user provided fasta sequence
					# classified into 3 catagories
					if len(unqfastaseq)>0:
						pepseq=str(jdic['Peptide Sequence']).strip()
						#if 
						#matchCount = tryptic_peptide.count(pepseq.upper())
						indices = [k for k in tryptic_peptide if pepseq.upper() in tryptic_peptide[k]]
						if len(indices)>0:
							tempuserseqheadermatch='NA'
							tempmatchlist=[]
							for i in indices:
								tempmatchlist.append('_'.join(fastaseq[i-1].split('_')[:-1]))
							if len(tempmatchlist)>0:
								tempuserseqheadermatch='<br/>'.join(tempmatchlist)
						
							if len(indices) > 1:
								seqhit=len(indices)
								jdic["Peptide in user's database"] =tempuserseqheadermatch
								jdic["Peptide unique in user's database"] ="Present but not unique"
							if len(indices) == 1:
								seqhit=len(indices)
								jdic["Peptide in user's database"] =tempuserseqheadermatch
								jdic["Peptide unique in user's database"] ="Present and unique"
								usersequnq.append("Present and unique")
							jfinaldata.append(jdic)
							
					else:
						jfinaldata.append(jdic)
			es.indices.refresh(index="mrmassaydb-index")
			#checking any result generated by database
			foundHits=len(jfinaldata)
			#storing only 10000 rows in json format
			json.dump(jfinaldata[:10000],jsonfileoutput)
			jsonfileoutput.close()
			# if result found then do other job
			if foundHits >0:
				statsummary=summaryStatcal(jfinaldata) # sent data to this funcation for generating stat
				pathwaychart=statsummary['pathwaychart']
				pathwaychart=[i[:2] for i in pathwaychart]
				specieslist=statsummary['specieslist']
				totallist=statsummary['total']
				subcell=statsummary['subcell']
				godic=statsummary['godic']
				jvennprot=statsummary['jevennstat'][0]
				jvennpep=statsummary['jevennstat'][1]
				mrmdatabase=statsummary['jevennstat'][2]
				sortedgodic=OrderedDict(sorted(godic.items(), key=lambda t: t[1])) # sorting GO data
				updatedgodic=dict(list(sortedgodic.items()))
				pepseqdataseries=ast.literal_eval(json.dumps(statsummary['pepseqdataseries'])) #dumping data into json format
				prodataseries=statsummary['prodataseries']
				unqisostat=statsummary['unqisostat']
				jsonfilepathStat=os.path.join(settings.BASE_DIR, 'resultFile', 'jsonData','resultJson', 'adavancesearch', 'statsummary', jsonfilename)
				jsonfileoutputStat= open(jsonfilepathStat,'w')
				json.dump(statsummary,jsonfileoutputStat)
				jsonfileoutputStat.close()
				urlname="'/resultFile/jsonData/resultJson/adavancesearch/results/"+jsonfilename+"'"
				if len(unqfastaseq)>0:
					tempcalunq=str(round(((float(usersequnq.count('Present and unique'))/float(len(jfinaldata)))*100),2))+'%'
					unqisostat.append(["User data",tempcalunq,"NA"])
					contextindex={
						"filename":urlname,"fastacolname":json.dumps(fastacolname),
						'query': finalsearhdata,'foundHits':foundHits,
						'pathwaychart':pathwaychart[:11],'specieslist':specieslist,
						'totallist':totallist,'subcell':subcell,
						'updatedgodic':updatedgodic,'pepseqdataseries':pepseqdataseries,
						'prodataseries':prodataseries,'unqisostat':unqisostat,
						'jvennprot':json.dumps(jvennprot),'jvennpep':json.dumps(jvennpep),'jvennmrmdb':json.dumps(mrmdatabase),'fastafilename':json.dumps(nameFIle)
						}
					return render(request,'resultformuserseq.html',contextindex)
				else:
					contextindex={
						"filename":urlname,"colname":json.dumps(colname),
						'query': finalsearhdata,'foundHits':foundHits,
						'pathwaychart':pathwaychart[:11],'specieslist':specieslist,
						'totallist':totallist,'subcell':subcell,
						'updatedgodic':updatedgodic,'pepseqdataseries':pepseqdataseries,
						'prodataseries':prodataseries,'unqisostat':unqisostat,
						'jvennprot':json.dumps(jvennprot),'jvennpep':json.dumps(jvennpep),'jvennmrmdb':json.dumps(mrmdatabase)
						}
					return render(request,'resultform.html',contextindex)
			else:
				return render(request,'resultform.html',{'foundHits':foundHits})
		else:
			return render(request,'resultform.html',{'foundHits':0})

# def fastaseq(request):
# 	ip = get_ip(request, right_most_proxy=True)
# 	IpAddressInformation.objects.create(ip_address=ip)
# 	if 'Uniprotkb' in request.GET and request.GET['Uniprotkb'] and 'pepseq' in request.GET and request.GET['pepseq']:
# 		uniprotkb=request.GET['Uniprotkb']
# 		uniprotkb=uniprotkb.strip()
# 		pepseq=request.GET['pepseq']
# 		pepseq=pepseq.strip()
# 		pepseqlist=[]
# 		pepseqlist.append(pepseq)
# 		reachable=True
# 		presentunidpepseqstat=False
# 		es.indices.refresh(index="mrmassaydb-index")
# 		query={"query": {
# 			"bool": {
# 				"must": [
# 					{"match": {"UniProtKB Accession": uniprotkb}},
# 					{"match": {"Peptide Sequence": pepseq}},
# 					{"match": {"UniprotKb entry status": "Yes"}}
# 				]
# 			}
# 		}
# 		}
# 		res = es.search(index="mrmassaydb-index", doc_type="mrmassaydb-type", body=query)

# 		foundHits=res["hits"]["total"]
# 		if foundHits >0:
# 			presentunidpepseqstat=True
# 		es.indices.refresh(index="mrmassaydb-index")
# 		if presentunidpepseqstat:
# 			try:
# 				sleep(random.randint(5,10))
# 				requests.get("https://www.uniprot.org/", timeout=5)
# 				proseq=''
# 				match_info=[]
# 				fasthead=''
# 				fastasq=''
# 				pepstartend=[]
# 				pdbstrucexist=False
# 				disorderstatlist=[]
# 				defaultPDBId=''
# 				unidatafasta = urllib.urlopen("https://www.uniprot.org/uniprot/" + uniprotkb + ".fasta")
# 				for mseq in SeqIO.parse(unidatafasta, "fasta"):
# 					fasthead=str((mseq.id).strip())
# 					proseq=str((mseq.seq).strip())
# 					fastasq=str((mseq.seq).strip())
# 				unidatafasta.close()

# 				peptide_pattern = re.compile(pepseq,re.IGNORECASE)
# 				for match in re.finditer(peptide_pattern,proseq):
# 					pepstartend=[(int(match.start())+1),match.end(),'peptide']
# 					disorderstatlist.append([(int(match.start())+1),'null','null',0,'null'])
# 					disorderstatlist.append([(int(match.start())+1),'null','null',1,'null'])
# 					disorderstatlist.append([int(match.end()),'null','null',1,'null'])
# 					disorderstatlist.append([int(match.end()),'null','null',0,'null'])
# 					match_info.append([uniprotkb,(int(match.start())+1),match.end(),match.group()])
# 				pepseqjava='["'+pepseq+'"]'
# 				try:
# 					sleep(random.randint(5,10))
# 					requests.get("https://iupred2a.elte.hu/", timeout=5)
# 					disfile = urllib.urlopen("https://iupred2a.elte.hu/iupred2a/short/" + uniprotkb)
# 					for disline in disfile:
# 						disdata = disline.strip()
# 						if not disdata.startswith('#') and '<' not in disdata and len(disdata)>1:
# 							disinfo = disdata.split()
# 							disdes="Residue:"+str(disinfo[1])+" Position:"+str(disinfo[0])+" &Score:"+str(disinfo[2])
# 							tempdisinfo=[int(disinfo[0]),float(disinfo[2]),str(disdes),'null','null']
# 							disorderstatlist.append(tempdisinfo)
# 					disorderstatlist.append([0,'null','null','null',0.5])
# 					disorderstatlist.append([int(len(fastasq)),'null','null','null',0.5])
# 					disfile.close()
# 				except ConnectionError as e:
# 					reachable=False
# 				fastalen=len(fastasq)
# 				if '-' not in uniprotkb:
# 					pdbuniprotlist=[]

# 					try:
# 						sleep(random.randint(5,10))
# 						requests.get("https://www.uniprot.org/", timeout=5)
# 						try:
# 							requestURL="https://www.uniprot.org/uniprot/"+str(uniprotkb)+".xml"
# 							unifile=urllib.urlopen(requestURL)
# 							unidata= unifile.read()
# 							unifile.close()
# 							try:
# 								unidata=minidom.parseString(unidata)
# 								try:
# 									pdbUnidata=(unidata.getElementsByTagName('dbReference'))
# 									for item in pdbUnidata:
# 										if (item.attributes['type'].value).upper() == 'PDB':
# 											try:
# 												pdbResolution=''
# 												pdbchaininfo=''
# 												pdbid=str(item.attributes['id'].value).strip()
# 												pdbMethod=(str(item.getElementsByTagName('property')[0].attributes['value'].value).strip())
# 												if pdbMethod.upper() == 'NMR':
# 													pdbchaininfo=(str(item.getElementsByTagName('property')[1].attributes['value'].value).strip()).split('=')
# 												else:
# 													pdbResolution=(str(item.getElementsByTagName('property')[1].attributes['value'].value).strip())
# 													pdbchaininfo=(str(item.getElementsByTagName('property')[2].attributes['value'].value).strip()).split('=')
# 												pdbchain=pdbchaininfo[0].strip()
# 												pdbStartEnd=pdbchaininfo[1].strip()
# 												pdbStart=int(pdbStartEnd.split('-')[0])
# 												pdbEnd=int(pdbStartEnd.split('-')[1])
# 												peptidecoverage=[]
# 												for mI in match_info:
# 													peptidecoverage.append(float(len(set(range(pdbStart,pdbEnd+1))&set(range(mI[1],mI[2]+1))))/float(len(pepseq)))
# 												if len(pdbResolution.strip())>0:
# 													pdbuniprotlist.append([pdbid,pdbMethod,pdbResolution,pdbchain,pdbStartEnd,round(max(peptidecoverage)*100,2)])
# 												else:
# 													if pdbMethod.upper() == 'NMR':
# 														pdbuniprotlist.append([pdbid,pdbMethod,pdbResolution,pdbchain,pdbStartEnd,round(max(peptidecoverage)*100,2)])
# 											except:
# 												pass
# 								except IndexError:
# 									pass

# 							except ExpatError:
# 								pass
# 						except IOError:
# 							pass

# 					except ConnectionError as e:
# 						reachable=False
# 					if len(pdbuniprotlist)>0:
# 						pdbstrucexist=True
# 					print(pdbuniprotlist)

# 				if reachable:
# 					defaultPDBId=pdbuniprotlist[0][0]
# 					#defaultPDBId="CYR1"
# 					return render(request, 'fastasequence.html', {
# 						'match_info':match_info,'fasthead':fasthead,\
# 						'fastasq':fastasq,'pepseqjava':pepseqjava,\
# 						'fastalen':fastalen,'pdbstrucexist':pdbstrucexist,\
# 						'disorderstatlist':json.dumps(disorderstatlist),\
# 						'pdbuniprotlist':pdbuniprotlist,'defaultPDBId':json.dumps(defaultPDBId),\
# 						'reachable':reachable} )
# 				else:
# 					return render(request, 'fastasequence.html', {'reachable':reachable})
# 			except ConnectionError as e:
# 				reachable=False
# 				return render(request, 'fastasequence.html', {'reachable':reachable})

def pathway(request):
	'''
	This function will display result for pathways.
	'''
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if 'Uniprotkb' in request.GET and request.GET['Uniprotkb'] and 'organismid' in request.GET and request.GET['organismid']:
		uniprotkb=request.GET['Uniprotkb']
		uniprotkb=uniprotkb.strip()
		listOfPathways=[]
		listOfPathwaySources=[]
		listOfPathwayLinks=[]
		listOfPCPathwayStatus=[]
		pathwaysInfo=[]
		presentunidstat=False
		OSid=request.GET['organismid']
		OSid=OSid.strip()

		homeURL=str((request.build_absolute_uri()).split('pathway')[0]).strip()
		if len(uniprotkb)>0 and len(OSid)>0:
			es.indices.refresh(index="mrmassaydb-index")
			query={"query": {
				"bool": {
					"must": [
						{"match": {"UniProtKB Accession": uniprotkb}},
						{"match": {"Organism ID": OSid}},
						{"match": {"UniprotKb entry status": "Yes"}}
					]
				}
			}
			}
			res = es.search(index="mrmassaydb-index", doc_type="mrmassaydb-type", body=query)
			foundHits=res["hits"]["total"]
			es.indices.refresh(index="mrmassaydb-index")
			if foundHits >0 :
				presentunidstat=True
				for hit in res['hits']['hits']:
					jdic=hit["_source"]
					jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
					jdic["sel"] =""
					if len(str(jdic['PC pathway URL']).strip())>0 and str(jdic['PC pathway URL']).strip().upper() !='NA':
						listOfPathwayLinks.extend(str(jdic['PC pathway URL']).strip().split('|'))
						listOfPathways.extend(str(jdic['PC pathway name']).strip().split('|'))
						listOfPathwaySources.extend(str(jdic['PC pathway source']).strip().split('|'))
						listOfPCPathwayStatus.extend([True]*len(str(jdic['PC pathway source']).strip().split('|')))
						break
				code=None
				unikeggid='NA'
				if '-' in uniprotkb:
					code=(str(uniprotkb).split('-'))[0]
				else:
					code=str(uniprotkb.strip())

				try:
					sleep(random.randint(5,10))
					requests.get("https://www.uniprot.org/", timeout=5)

					try:
						requestURL="https://www.uniprot.org/uniprot/"+str(code)+".xml"
						unifile=urllib.urlopen(requestURL)
						unidata= unifile.read()
						unifile.close()
						try:
							unidata=minidom.parseString(unidata)
							try:
								keggUnidata=(unidata.getElementsByTagName('dbReference'))
								for item in keggUnidata:
									if (item.attributes['type'].value).upper() == 'KEGG':
										try:
											unikeggid=str(item.attributes['id'].value).strip()
											break
										except:
											pass
							except IndexError:
								pass


						except ExpatError:
							pass
					except IOError:
						pass

				except ConnectionError as e:
					pass
				if len(unikeggid) >0:
					k = KEGG()
					kegg=k.get(unikeggid)
					dict_data = k.parse(kegg)
					try:
						keggpathwayid= (dict_data['PATHWAY'].keys())
						for kegpathidietm in keggpathwayid:
							keggentryid = (unikeggid.split(':'))[1].strip()
							subkeggmapid=str(kegpathidietm)+'+'+str(keggentryid)
							if len(subkeggmapid) >0:
								kegggeneiddicpeptrack={}
								temppathwayname=str(dict_data['PATHWAY'][kegpathidietm]).strip()
								if str(kegpathidietm) not in '\t'.join(listOfPathwayLinks):
									listOfPathways.append(temppathwayname)
									listOfPathwayLinks.append("https://www.kegg.jp/kegg-bin/show_pathway?"+str(kegpathidietm))
									listOfPathwaySources.append('kegg')
									listOfPCPathwayStatus.append(False)

					except KeyError:
						pass

			if presentunidstat:
				for i in zip(listOfPathways,listOfPathwayLinks,listOfPathwaySources,listOfPCPathwayStatus):
						tempList=[]
						tempList.append(i[0])
						tempPathId=i[1].split('/')[-1]
						if i[3] == False:
							tempPathId=i[1].split('?')[-1]
						tempList.append(tempPathId)
						tempList.append(i[1])
						tempList.append(i[2])
						tempList.append(i[3])
						pathwaysInfo.append(tempList)
				es.indices.refresh(index="mrmassaydb-index")

		return render(request, 'pathway.html', {'uniprotid':uniprotkb,'pathwaysInfo':pathwaysInfo,'OSid':OSid} )


def ppi(request):
	'''
	This function will display result for STRING PPI.
	'''
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if 'Uniprotkb' in request.GET and request.GET['Uniprotkb'] and 'organismid' in request.GET and request.GET['organismid']:
		uniprotkb=request.GET['Uniprotkb']
		uniprotkb=uniprotkb.strip()
		stringdbaid=''
		uniprotname=''
		nodesscript=[]
		edgesscript=[]
		nodeslist=[]
		reachable=True
		presentunidstat=False
		OSid=request.GET['organismid']
		OSid=OSid.strip()
		stringjsondata=[]
		searchtermlist=[]
		uniGeneDic={}
		stringScoreDic={
			'score':'combined',
			'nscore':'gene neighborhood',
			'fscore':'gene fusion',
			'pscore':'phylogenetic profile',
			'ascore':'coexpression',
			'escore':'experimental',
			'dscore':'database',
			'tscore':'textmining'
		}
		pepfilepath = os.path.join(settings.BASE_DIR, 'mappermotherfile', 'ReportBook_mother_file.csv')
		homeURL=str((request.build_absolute_uri()).split('ppi')[0]).strip()
		if len(uniprotkb)>0 and len(OSid)>0:
			pepfilegenidlistpeptrack=[]
			pepfilegenidlistremain=[]
			es.indices.refresh(index="mrmassaydb-index")
			query={"query": {
				"bool": {
					"must": [
						{"match": {"UniProtKB Accession": uniprotkb}},
						{"match": {"Organism ID": OSid}},
						{"match": {"UniprotKb entry status": "Yes"}}
					]
				}
			}
			}
			res = es.search(index="mrmassaydb-index", doc_type="mrmassaydb-type", body=query)
			foundHits=res["hits"]["total"]
			es.indices.refresh(index="mrmassaydb-index")
			if foundHits >0 :
				presentunidstat=True
			if presentunidstat:

				code=None
				if '-' in uniprotkb:
					code=(str(uniprotkb).split('-'))[0]
				else:
					code=str(uniprotkb.strip())

				try:
					sleep(random.randint(5,10))
					requests.get("https://www.uniprot.org/", timeout=5)
					try:
						requestURL="https://www.uniprot.org/uniprot/"+str(code)+".xml"
						unifile=urllib.urlopen(requestURL)
						unidata= unifile.read()
						unifile.close()
						try:
							unidata=minidom.parseString(unidata)
							try:
								stringdata=(unidata.getElementsByTagName('dbReference'))
								for item in stringdata:
									if (item.attributes['type'].value).upper() == 'STRING':
										try:
											stringdbaid=str(item.attributes['id'].value).strip()
										except:
											pass
							except IndexError:
								pass

							try:
								try:
									uniprotid=unidata.getElementsByTagName('name')[0].firstChild.nodeValue
									uniprotid=str(uniprotid)
								except:
									uniprotid='NA'
							except IndexError:
								pass

							try:
								try:
									uniprotname=((unidata.getElementsByTagName('gene')[0]).getElementsByTagName('name')[0]).firstChild.nodeValue
									uniprotname=str(uniprotname)
								except:
									uniprotname='NA'
							except IndexError:
								pass


						except ExpatError:
							pass
					except IOError:
						pass

				except ConnectionError as e:
					reachable=False

			if len(stringdbaid.strip()) >0 and reachable:
				OSid=str(int(float(OSid)))
				try:
					sleep(random.randint(5,10))
					requests.get("https://string-db.org/", timeout=5)
					stringUrl="https://string-db.org/api/json/interactionsList?identifiers="+stringdbaid+"&limit=50"
					if not requests.get(stringUrl).ok:
					  requests.get(stringUrl).raise_for_status()
					  sys.exit()
					stringdbData = requests.get(stringUrl).text
					# convert 'str' to Json
					stringdbData = json.loads(stringdbData)
					for sdbwline in stringdbData:
						if len(sdbwline['stringId_A'].strip()) >0  and len(sdbwline['stringId_B'].strip()):
							nodeslist.append(str(sdbwline['preferredName_A']).strip())
							nodeslist.append(str(sdbwline['preferredName_B']).strip())
							
							for stringScoreKey in stringScoreDic:
								if float(sdbwline[stringScoreKey]) >0:
									if (str(uniprotname).upper()).strip() == (str(sdbwline['preferredName_A']).upper()).strip() or (str(uniprotname).upper()).strip() ==(str(sdbwline['preferredName_B']).upper()).strip():
										edgesscript.append({"data": { "source": (str(sdbwline['preferredName_A']).upper()).strip(), "target": (str(sdbwline['preferredName_B']).upper()).strip(),"weight": float(sdbwline[stringScoreKey]),"group":stringScoreDic[stringScoreKey]},"group": "edges","grabbable": True })
									else:
										edgesscript.append({"data": { "source": (str(sdbwline['preferredName_A']).upper()).strip(), "target": (str(sdbwline['preferredName_B']).upper()).strip(),"weight": float(sdbwline[stringScoreKey]),"group":stringScoreDic[stringScoreKey]},"group": "edges","grabbable": True })

				except ConnectionError as e:
					reachable=False
			if len(nodeslist) >0:
				uniqnodelist=list(set(nodeslist))
				shouldlist=[]
				for x in uniqnodelist:
					tempquery={
								"multi_match":{
									"query":x.strip(),
									"type":"best_fields",
									"fields":["Gene.ngram"],
									"minimum_should_match":"100%"
								}
							}
					shouldlist.append(tempquery)
				booldic={}
				booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
				searchtermlist.append(booldic)
				if len(OSid)>0:
					ogshouldlist=[]
					ogtempquery={
								"multi_match":{
									"query":OSid.strip(),
									"type":"best_fields",
									"fields":["Organism ID.ngram"],
									"minimum_should_match":"100%"
								}
							}
					ogshouldlist.append(ogtempquery)
					ogbooldic={}
					ogbooldic["bool"]={"should":ogshouldlist,"minimum_should_match": 1}
					searchtermlist.append(ogbooldic)					

			if len(searchtermlist)>0:
				es.indices.refresh(index="mrmassaydb-index")

				query={
					"query": {
						"bool": {
							"must":searchtermlist
						}
					}
				}

				resOrg=helpers.scan(client=es,scroll='2m',index="mrmassaydb-index", doc_type="mrmassaydb-type",query=query,request_timeout=30)
				jfinaldata=[]
				for i in resOrg:
					jdic=i['_source']
					jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
					uniGeneDic[jdic['Gene'].lower()]=jdic["UniProtKB Accession"]
					jfinaldata.append(jdic)
				es.indices.refresh(index="mrmassaydb-index")
				databaseNameList=['PeptideTracker ID','Passel ID','SRMAtlas ID','Cptac ID','Panoramaweb ID']
				for i in databaseNameList:
					subdatafilter=filter(lambda pepdata: ((str(pepdata[i]).strip()).lower() !='na' and len(str(pepdata[i]).strip()) >0), jfinaldata)
					subdatagene=list((object['Gene'] for object in subdatafilter))
					subdatagene=filter(None, subdatagene)
					if i =="PeptideTracker ID":
						pepfilegenidlistpeptrack=list(set(subdatagene))
					else:
						pepfilegenidlistremain.append(list(set(subdatagene)))

				pepfilegenidlistremain=list(set(reduce(operator.concat, pepfilegenidlistremain)))
				pepfilegenidlistpeptrack=[x.lower() for x in pepfilegenidlistpeptrack]
				pepfilegenidlistremain=[x.lower() for x in pepfilegenidlistremain]
				jfinaldata=[]

		if len(stringdbaid.strip()) >0 and reachable:
			jsonfilestatus=False
			jsonfilename=str((str(uniprotkb).split('-'))[0]).strip()+'.json'
			jsonfilepath=os.path.join(settings.BASE_DIR, 'resultFile','jsonData', 'stringJson', jsonfilename)
			if os.path.exists(jsonfilepath):
				if datetime.datetime.fromtimestamp(os.path.getmtime(pepfilepath))< datetime.datetime.fromtimestamp(os.path.getmtime(jsonfilepath)):
					jsonfilestatus=True
					jsonprotvistastatus=True
					stringjsondata=json.load(open(jsonfilepath))
				else:
					os.remove(jsonfilepath)
					jsonfilestatus=False
			else:
				jsonfilestatus=False
			if not jsonfilestatus:
				if len(nodeslist) >0:
					uniqnodelist=list(set(nodeslist))
					if (str(uniprotname).upper()).strip() in list(map(lambda x:x.upper(),uniqnodelist)):
						for nitem in uniqnodelist:
							if (str(uniprotname).upper()).strip() == nitem:
								nodesscript.append({"data": { "id": nitem.upper().strip(), "name": nitem.upper().strip(),"Uniprotkb":uniGeneDic[(str(nitem)).lower()], "organismid": OSid, "href": homeURL+"search/hyperlink/?", "score": 0.5, "gene": True,"nodecolor": "violet"},"group": "nodes","grabbable": True, "classes": "fn10273 fn6944 fn9471 fn10569 fn8023 fn6956 fn6935 fn8147 fn6939 fn6936 fn6629 fn7928 fn6947 fn8612 fn6957 fn8786 fn6246 fn9367 fn6945 fn6946 fn10024 fn10022 fn6811 fn9361 fn6279 fn6278 fn8569 fn7641 fn8568 fn6943"})
							else:
								if (str(nitem)).lower() in pepfilegenidlistpeptrack:
									nodesscript.append({"data": { "id": nitem.upper().strip(), "name": nitem.upper().strip(),"Uniprotkb":uniGeneDic[(str(nitem)).lower()], "organismid":OSid, "href": homeURL+"search/hyperlink/?", "score": 0.35, "gene": True,"nodecolor": "#C64646"},"group": "nodes","grabbable": True, "classes": "fn10273 fn6944 fn9471 fn10569 fn8023 fn6956 fn6935 fn8147 fn6939 fn6936 fn6629 fn7928 fn6947 fn8612 fn6957 fn8786 fn6246 fn9367 fn6945 fn6946 fn10024 fn10022 fn6811 fn9361 fn6279 fn6278 fn8569 fn7641 fn8568 fn6943"})
								elif (str(nitem)).lower() in pepfilegenidlistremain:
									nodesscript.append({"data": { "id": nitem.upper().strip(), "name": nitem.upper().strip(),"Uniprotkb":uniGeneDic[(str(nitem)).lower()], "organismid":OSid, "href": homeURL+"search/hyperlink/?", "score": 0.35, "gene": True,"nodecolor": "orange"},"group": "nodes","grabbable": True, "classes": "fn10273 fn6944 fn9471 fn10569 fn8023 fn6956 fn6935 fn8147 fn6939 fn6936 fn6629 fn7928 fn6947 fn8612 fn6957 fn8786 fn6246 fn9367 fn6945 fn6946 fn10024 fn10022 fn6811 fn9361 fn6279 fn6278 fn8569 fn7641 fn8568 fn6943"})
								else:
									nodesscript.append({"data": { "id": nitem.upper().strip(), "name": nitem.upper().strip(),"Uniprotkb":'NA', "organismid":OSid, "href": homeURL+"search/hyperlink/?","score": 0.25, "gene": True,"nodecolor": "skyblue"},"group": "nodes","grabbable": True, "classes": "fn10273 fn6944 fn9471 fn10569 fn8023 fn6956 fn6935 fn8147 fn6939 fn6936 fn6629 fn7928 fn6947 fn8612 fn6957 fn8786 fn6246 fn9367 fn6945 fn6946 fn10024 fn10022 fn6811 fn9361 fn6279 fn6278 fn8569 fn7641 fn8568 fn6943"})
						stringjsondata=nodesscript+edgesscript
						jsonfileoutput= open(jsonfilepath,'w')
						json.dump(stringjsondata,jsonfileoutput)
						jsonfileoutput.close()
					else:
						edgesscript=[]
		if presentunidstat:
			cy_style_path=os.path.join(settings.BASE_DIR, 'resultFile','jsonData', 'cytoscapeStyleJson', 'cy-style.json')
			cystylejson=json.load(open(cy_style_path))
			if reachable:
				return render(request, 'ppi.html', {'cystylejson':json.dumps(cystylejson),'stringjsondata':json.dumps(stringjsondata),'stringdbaid':stringdbaid,'uniprotid':uniprotid,'uniprotname':uniprotname,'reachable':reachable} )
			else:
				return render(request, 'ppi.html', {'reachable':reachable})


def pathwayview(request):
	'''
	This function will display result for KEGG pathways.
	'''
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if 'Uniprotkb' in request.GET and request.GET['Uniprotkb'] and 'organismid' in request.GET and request.GET['organismid'] and 'pathwayid' in request.GET and request.GET['pathwayid'] and 'pathwayname' in request.GET and request.GET['pathwayname']:
		uniprotkb=request.GET['Uniprotkb']
		uniprotkb=uniprotkb.strip()
		uniprotname=''
		unikeggid=''
		nodesscript=[]
		edgesscript=[]
		uniGeneDic={}
		keggurl=''
		reachable=True
		presentunidstat=False
		OSid=request.GET['organismid']
		OSid=OSid.strip()
		pathwayid=request.GET['pathwayid']
		pathwayid=pathwayid.strip()
		pathwayname=request.GET['pathwayname']
		pathwayname=pathwayname.strip()
		homeURL=str((request.build_absolute_uri()).split('pathwayview')[0]).strip()
		if len(uniprotkb)>0:
			pepfilegenidlistpeptrack=[]
			pepfilegenidlistremain=[]

			es.indices.refresh(index="mrmassaydb-index")
			query={"query": {
				"bool": {
					"must": [
						{"match": {"UniProtKB Accession": uniprotkb}},
						{"match": {"Organism ID": OSid}},
						{"match": {"UniprotKb entry status": "Yes"}}
					]
				}
			}
			}
			res = es.search(index="mrmassaydb-index", doc_type="mrmassaydb-type", body=query)

			foundHits=res["hits"]["total"]
			es.indices.refresh(index="mrmassaydb-index")
			if foundHits >0 :
				presentunidstat=True
			if presentunidstat:
				viewquery={"query": {
					"bool": {
						"must": [
							{"match": {"Pathway Name": pathwayname}},
							{"match": {"Organism ID": OSid}},
							{"match": {"UniprotKb entry status": "Yes"}}
						]
					}
				}
				}
				resOrg = helpers.scan(client=es,scroll='2m',index="mrmassaydb-index", doc_type="mrmassaydb-type",query=viewquery)
				
				jfinaldata=[]
				for i in resOrg:
					jdic=i['_source']
					jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
					jdic["sel"] =""
					uniGeneDic[jdic['Gene'].lower()]=jdic["UniProtKB Accession"]
					jfinaldata.append(jdic)
				es.indices.refresh(index="mrmassaydb-index")
				databaseNameList=['PeptideTracker ID','Passel ID','SRMAtlas ID','Cptac ID','Panoramaweb ID']
				for i in databaseNameList:
					subdatafilter=filter(lambda pepdata: ((str(pepdata[i]).strip()).lower() !='na' and len(str(pepdata[i]).strip()) >0), jfinaldata)
					subdatagene=list((object['Gene'] for object in subdatafilter))
					subdatagene=filter(None, subdatagene)
					if i =="PeptideTracker ID":
						pepfilegenidlistpeptrack=list(set(subdatagene))
					else:
						pepfilegenidlistremain.append(list(set(subdatagene)))

				pepfilegenidlistremain=list(set(reduce(operator.concat, pepfilegenidlistremain)))
				pepfilegenidlistpeptrack=[x.lower() for x in pepfilegenidlistpeptrack]
				pepfilegenidlistremain=[x.lower() for x in pepfilegenidlistremain]
				code=None
				if '-' in uniprotkb:
					code=(str(uniprotkb).split('-'))[0]
				else:
					code=str(uniprotkb.strip())

				try:
					sleep(random.randint(5,10))
					requests.get("https://www.uniprot.org/", timeout=5)

					try:
						requestURL="https://www.uniprot.org/uniprot/"+str(code)+".xml"
						unifile=urllib.urlopen(requestURL)
						unidata= unifile.read()
						unifile.close()
						try:
							unidata=minidom.parseString(unidata)
							try:
								keggUnidata=(unidata.getElementsByTagName('dbReference'))
								for item in keggUnidata:
									if (item.attributes['type'].value).upper() == 'KEGG':
										try:
											unikeggid=str(item.attributes['id'].value).strip()
										except:
											pass
							except IndexError:
								pass

							try:
								try:
									uniprotid=unidata.getElementsByTagName('name')[0].firstChild.nodeValue
									uniprotid=str(uniprotid)
								except:
									uniprotid='NA'
							except IndexError:
								pass

							try:
								try:
									uniprotname=((unidata.getElementsByTagName('gene')[0]).getElementsByTagName('name')[0]).firstChild.nodeValue
									uniprotname=str(uniprotname)
								except:
									uniprotname='NA'
							except IndexError:
								pass


						except ExpatError:
							pass
					except IOError:
						pass

				except ConnectionError as e:
					reachable=False

		keggimagedict={}
		if len(unikeggid) >0 and reachable:
			k = KEGG()
			kegg=k.get(unikeggid)
			dict_data = k.parse(kegg)
			try:
				keggpathwayid= (dict_data['PATHWAY'].keys())
				for kegpathidietm in keggpathwayid:
					keggentryid = (unikeggid.split(':'))[1].strip()
					subkeggmapid=str(kegpathidietm)+'+'+str(keggentryid)
					if len(subkeggmapid) >0 and str(kegpathidietm).lower() == pathwayid.lower():
						kegggeneiddicpeptrack=[]
						kegggeneiddicremain=[]
						notpresentkegggeneid=[]
						updatedcorddatalist=[]
						keggimagepath=''
						pathway = KGML_parser.read(REST.kegg_get(kegpathidietm, "kgml"))
						pathwayName=(str(pathway)).split('\n')[0].split(':')[1].strip()
						# colorspiral = ColorSpiral()
						# colorlist = colorspiral.get_colors(len(pathway.genes))
						listUnqGene=[]
						keggOrgInitial=''
						getGeneList=[]
						for g in pathway.genes:
							tempGInfo=(g.name).split( )
							for tg in tempGInfo:
								geneID=tg.split(':')[1].strip()
								keggOrgInitial=tg.split(':')[0].strip()
								getGeneList.append(tg.strip())
								listUnqGene.append(geneID)
						listUnqGene=list(set(listUnqGene))
						getGeneList=list(set(getGeneList))
						geneKeggIdDics={}
						for gI in range(0,len(getGeneList),100):
							keggRESTURL='http://rest.kegg.jp/list/'+'+'.join(getGeneList[gI:gI+100])
							keggResponseREST = requests.get(keggRESTURL,verify=False)
							if not keggResponseREST.ok:
								keggResponseREST.raise_for_status()
								sys.exit()
							keggResponseBodyREST = keggResponseREST.text
							for gN in keggResponseBodyREST.split('\n')[:-1]:
								gnInfo=gN.split('\t')
								tempGeneList=[i.strip().lower() for i in gnInfo[1].split(';')[0].split(',')]
								tempKeggGeneID=gnInfo[0].split(':')[1]
								geneKeggIdDics[tempKeggGeneID]=[i.strip() for i in gnInfo[1].split(';')[0].split(',')]
								if gnInfo[0] != unikeggid:
									for tg in tempGeneList:
										if tg in pepfilegenidlistpeptrack:
											kegggeneiddicpeptrack.append(tempKeggGeneID)
										elif tg in pepfilegenidlistremain:
											kegggeneiddicremain.append(tempKeggGeneID)

						kegggeneiddicpeptrack=list(set(kegggeneiddicpeptrack))
						kegggeneiddicremain=list(set(kegggeneiddicremain))
						notpresentkegggeneid=set(geneKeggIdDics.keys())-set(kegggeneiddicpeptrack+kegggeneiddicremain+[unikeggid.split(':')[1]])
						notpresentkegggeneid=list(notpresentkegggeneid)
												
						presentKeggGeneIDList=set(kegggeneiddicpeptrack+kegggeneiddicremain+[unikeggid.split(':')[1]])
						presentKeggGeneIDList=list(presentKeggGeneIDList)
						presentGeneList=list(set(pepfilegenidlistpeptrack+pepfilegenidlistremain+[uniprotname.lower()]))
						presentGeneIDListWithoutUserGeneName=presentKeggGeneIDList
						presentGeneIDListWithoutUserGeneName.remove(unikeggid.split(':')[1])

						# for color, element in zip(colorlist, pathway.genes):
						# 	x=element.graphics[0].name.split(',')[0].strip()
						# 	element.graphics[0].name=x.split(',')[0].strip()
						# 	element.graphics[0].y=str(int(element.graphics[0].y)-1)
						# 	tempArray=element.name.split(' ')
						# 	updatedtempArray=[i.split(':')[1].split(' ')[0] for i in tempArray]
						# 	matchResArray=list(set(updatedtempArray) & set(presentGeneIDListWithoutUserGeneName))
						# 	matchResArrayWithUserGene=list(set(updatedtempArray) & set([unikeggid.split(':')[1]]))
						# 	kshape='rect'
						# 	x1=element.graphics[0].x
						# 	y1=element.graphics[0].y
						# 	x2=element.graphics[0].x+element.graphics[0].width
						# 	y2=element.graphics[0].y+element.graphics[0].height

						# 	#
						# 	#
						# 	tempTitle=[]
						# 	commonKeggGeneList=[]
						# 	if len(matchResArray)>0:
						# 		for m1 in matchResArray:
						# 			tempMatchGene1=list(set(map(lambda i:i.lower(),geneKeggIdDics[m1])) & set(presentGeneList))
						# 			tempTitle.append(m1+' ('+tempMatchGene1[0]+')')
						# 			commonKeggGeneList.append(tempMatchGene1[0])
						# 		element.graphics[0].bgcolor = 'orange'
						# 	if len(matchResArrayWithUserGene)>0:
						# 		for m2 in matchResArray:
						# 			tempMatchGene2=list(set(map(lambda i:i.lower(),geneKeggIdDics[m2])) & set(presentGeneList))
						# 			tempTitle.append(m2+' ('+tempMatchGene2[0]+')')
						# 			commonKeggGeneList.append(tempMatchGene2[0])
						# 		element.graphics[0].bgcolor = 'violet'
						# 	if len(tempTitle)>0:
						# 		ktitleData=','.join(tempTitle)
						# 		khref='/search/hyperlink/?Gene='+str('|'.join(commonKeggGeneList))+'&Organismid='+OSid
						# 		kcoords=str(x1)+','+str(y1)+','+str(x2)+','+str(y2)
						# 		updatedcorddatalist.append([kshape,kcoords,str(khref),ktitleData])
						# canvas = KGMLCanvas(pathway, import_imagemap=True,show_genes=True,fontsize=11)
						# originalkeggIMGURL='https://www.kegg.jp/kegg/pathway/'+str(unikeggid.split(':')[0])+'/'+kegpathidietm+'.png'
						# kwidth,kheight = Image.open(requests.get(originalkeggIMGURL,stream=True).raw).size
						# keggimgfilenamePDF=kegpathidietm+'_'+unikeggid.split(':')[1]+'.pdf'
						# keggimgfilenamePNG=kegpathidietm+'_'+unikeggid.split(':')[1]+'.png'
						# keggimgfilepathPDF=os.path.join(settings.BASE_DIR, 'resultFile', 'image','kegg', keggimgfilenamePDF)
						# keggimgfilepathPNG=os.path.join(settings.BASE_DIR, 'resultFile', 'image','kegg', keggimgfilenamePNG)
						# canvas.draw(keggimgfilepathPDF)
						# pathwayImage = convert_from_path(keggimgfilepathPDF, dpi=300, thread_count=4)
						# pathwayImage[0].save(keggimgfilepathPNG, 'PNG')
						# pimg=Image.open(keggimgfilepathPNG)
						# pimg=pimg.resize((kwidth,kheight),Image.ANTIALIAS)
						# pimg.save(keggimgfilepathPNG)

						keggurl = "http://www.kegg.jp/kegg-bin/show_pathway?" + str(kegpathidietm)
						# if len(kegggeneiddicpeptrack)>0:
						# 	for pgc in kegggeneiddicpeptrack:
						# 		keggurl += "/%s%%20%s/" % (pgc,'maroon')
						# if len(kegggeneiddicremain)>0:
						# 	for rgc in kegggeneiddicremain:
						# 		keggurl += "/%s%%20%s/" % (rgc,'orange')
						# if len(notpresentkegggeneid)>0:
						# 	for ngc in notpresentkegggeneid:
						# 		keggurl += "/%s%%20%s/" % (ngc,'#ffffff')
						try:
							sleep(random.randint(5,10))
							requests.get("https://www.kegg.jp/", timeout=1)
							kegghttp=urllib3.PoolManager()
							urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
							keggurlreq = kegghttp.request('GET',keggurl)
							keggurldata =keggurlreq.data
							keggurldata=keggurldata.decode('utf-8')
							keggurldata=keggurldata.replace('\t','')
							for kline in keggurldata.split('\n'):
								kline=kline.lstrip()
								try:
									if '<area' in kline.lower() and 'href='.lower() in kline.lower() and 'title=' in kline.lower() and 'shape=' in kline.lower():
										kline=kline.replace('"','')
										kshape=kline.split('shape=')[1].split(' ')[0]
										kcoords=kline.split('coords=')[1].split(' ')[0]
										khref=''
										ktitle=kline.split('title=')[1].replace('/>','').replace('class=module','')
										ktitleData=ktitle.lstrip()
										kclass=''
										if kshape=='rect':
											ktitleInfo=ktitleData.split(',')
											ktitleInfoID=[k.strip().split(' ')[0] for k in ktitleInfo]
											commonKeggGeneIDListWithoutUserGeneName=list(set(ktitleInfoID)&set(presentGeneIDListWithoutUserGeneName))
											matchResArrayWithUserGene=list(set(ktitleInfoID) & set([unikeggid.split(':')[1]]))
											if len(commonKeggGeneIDListWithoutUserGeneName)>0:
												tempTitle=[]
												commonKeggGeneList=[]
												for ki in ktitleInfoID:
													if ki in geneKeggIdDics:
														keggGeneList=geneKeggIdDics[ki]
														tempcommonKeggGeneList=[str(cg) for cg in keggGeneList if cg.lower() in presentGeneList]
														tempTitle.extend([ki+' ('+str(x)+')' for x in tempcommonKeggGeneList])
														commonKeggGeneList.extend(tempcommonKeggGeneList)
													
												ktitleData=','.join(tempTitle)
												commonUniIdList=[uniGeneDic[cg.lower()] for cg in commonKeggGeneList]
												khref='/search/hyperlink/?UniProtKB Accession='+str('|'.join(commonUniIdList))
												kclass='originalMRMAssayDB'
											elif len(matchResArrayWithUserGene)>0:
												tempTitle=[]
												commonKeggGeneList=[]
												for ki in ktitleInfoID:
													if ki in geneKeggIdDics:
														keggGeneList=geneKeggIdDics[ki]
														tempcommonKeggGeneList=[str(cg) for cg in keggGeneList if cg.lower() in presentGeneList]
														tempTitle.extend([ki+' ('+str(x)+')' for x in tempcommonKeggGeneList])
														commonKeggGeneList.extend(tempcommonKeggGeneList)
													
												ktitleData=','.join(tempTitle)
												commonUniIdList=[uniGeneDic[cg.lower()] for cg in commonKeggGeneList]
												khref='/search/hyperlink/?UniProtKB Accession='+str('|'.join(commonUniIdList))
												kclass='originalUser'
											else:
													
												ktitleData=''
												khref=''
												kclass='original'	
											updatedcorddatalist.append([kshape,kcoords,str(khref),ktitleData,kclass])											
									
									if kline.startswith('<img src="/kegg/pathway'):
										kinfo=kline.split('name=')
										keggimagepath=(kinfo[0].split('"'))[1].strip()
								except UnicodeDecodeError:
									pass
						except ConnectionError as e:
							reachable=False
						except urllib3.exceptions.MaxRetryError:
							reachable=False
						#keggimagedict[subkeggmapid]=[keggimgfilenamePNG,pathwayName,updatedcorddatalist]
						keggimagedict[subkeggmapid]=[keggimagepath,pathwayName,updatedcorddatalist]

			except KeyError:
				pass
		if presentunidstat:
			if reachable:
				contextindex= {
					'uniprotid':uniprotid,'uniprotname':uniprotname,
					'keggimagedict':keggimagedict,'keggurl':keggurl,
					'reachable':reachable
				}
				return render(request, 'pathwayview.html',contextindex)
			else:
				return render(request, 'pathwayview.html', {'reachable':reachable})


def mutptmdom(request):
	'''
	This function will display result for Mutation,PTM & Domain.
	'''
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if 'Uniprotkb' in request.GET and request.GET['Uniprotkb'] and 'pepseq' in request.GET and request.GET['pepseq']:
		uniprotkb=request.GET['Uniprotkb']
		uniprotkb=uniprotkb.strip()
		molartuniID="'"+str(uniprotkb)+"'"
		pepseq=request.GET['pepseq']
		pepseq=pepseq.strip()
		domainplotscript=[]
		ptmmutplotscript=[]
		valid=False
		contextpep={}
		proseq=''
		match_info=[]
		fasthead=''
		fastasq=''
		pepstartend=[]
		fastalen=0
		peptidematchstartend=[]
		reachable=True
		pepfilepath = os.path.join(settings.BASE_DIR, 'mappermotherfile', 'ReportBook_mother_file.csv')
		presentunidpepseqstat=False
		listOfPeptide=[]
		pepstart=0
		pepend=0
		jsonmolartstatus=False
		protname=None
		es.indices.refresh(index="mrmassaydb-index")
		query={"query": {
			"bool": {
				"must": [
					{"match": {"UniProtKB Accession": uniprotkb}},
					{"match": {"UniprotKb entry status": "Yes"}}
				]
			}
		}
		}

		res=helpers.scan(client=es,scroll='2m',index="mrmassaydb-index", doc_type="mrmassaydb-type",query=query,request_timeout=30)
		for hit in res:
			jdic=hit["_source"]
			jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
			jdic["sel"] =""
			if (str(jdic["Peptide Sequence"]).strip()).upper() == (str(pepseq)).upper():
				presentunidpepseqstat=True
				protname=str(jdic["Gene"]).strip()
			listOfPeptide.append(str(jdic["Peptide Sequence"]).strip())
		es.indices.refresh(index="mrmassaydb-index")
		foundHits=len(listOfPeptide)
		if '-' not in uniprotkb and presentunidpepseqstat:
			code=(str(uniprotkb).strip()).upper()
			try:
				sleep(random.randint(5,10))
				requests.get("https://www.uniprot.org/", timeout=1)
				unidatafasta = urllib.urlopen("https://www.uniprot.org/uniprot/" + code + ".fasta")
				for mseq in SeqIO.parse(unidatafasta, "fasta"):
					fasthead=str((mseq.id).strip())
					proseq=str((mseq.seq).strip())
					fastasq=str((mseq.seq).strip())
				unidatafasta.close()
			except ConnectionError as e:
				reachable=False


			if reachable:
				jsonfilestatus=False
				jsonfilename='externalLabeledFeatures_'+uniprotkb+'.json'
				jsonfilepath=os.path.join(settings.BASE_DIR, 'resultFile','jsonData', 'protvistadataJson', jsonfilename)
				if os.path.exists(jsonfilepath):
					if datetime.datetime.fromtimestamp(os.path.getmtime(pepfilepath))< datetime.datetime.fromtimestamp(os.path.getmtime(jsonfilepath)):
						jsonfilestatus=True
						jsonmolartstatus=True
					else:
						os.remove(jsonfilepath)
						jsonfilestatus=False
				else:
					jsonfilestatus=False
				jsonpepdata={}
				for pepseqitem in listOfPeptide:
					peptide_pattern = re.compile(pepseqitem,re.IGNORECASE)
					for match in re.finditer(peptide_pattern,proseq):
						pepstartend=[(int(match.start())+1),match.end(),'peptide']
						peptidematchstartend.append((str((int(match.start())+1)).strip())+"-"+(str(match.end()).strip()))
						jsonpepdata[pepseqitem.upper()] =[(int(match.start())+1),match.end()]
						if pepseqitem.lower() == pepseq.lower():
							pepstart=int(match.start())+1
							pepend=int(match.end())
							match_info.append([(int(match.start())+1),match.group(),match.end()])
				contextpep[uniprotkb]=[list(x) for x in set(tuple(x) for x in match_info)]
				pepseqjava='["'+pepseq+'"]'
				fastalen=len(fastasq)

				if not jsonfilestatus:
					jsonfileoutput= open(jsonfilepath,'w')
					jsonformatdatamolart={}
					jsonformatdatamolart["sequence"]=str(proseq)
					tempfearures=[]
					for jsonkey in jsonpepdata.keys():
						tempdicjsondic={}
						tempdicjsondic["type"]="MRM"
						tempdicjsondic["category"]="Targeted_Proteomics_Assay "
						tempdicjsondic["description"]="Suitable MRM Assay"
						tempdicjsondic["begin"]=str(jsonpepdata[jsonkey][0])
						tempdicjsondic["end"]=str(jsonpepdata[jsonkey][1])
						tempdicjsondic["color"]="#00B88A"
						tempdicjsondic["accession"]=uniprotkb
						tempfearures.append(tempdicjsondic)
					jsonformatdatamolart["features"]=tempfearures
					json.dump(jsonformatdatamolart,jsonfileoutput)
					jsonfileoutput.close()
					jsonmolartstatus=True

		for pepmatchitem in peptidematchstartend:
			domainplotscript.append({"name":"Peptide","coord":pepmatchitem})
		if reachable:
			return render(request, 'mutptmdom.html',{'contextpep':contextpep,\
				'fastalen':fastalen,'uniprotkb':uniprotkb,\
				'molartuniID':molartuniID,'pepstart':pepstart,\
				'pepend':pepend,'jsonmolartstatus':jsonmolartstatus,'protname':protname,'reachable':reachable})
		else:
			return render(request, 'mutptmdom.html',{'reachable':reachable})

def disease(request):
	'''
	This function will display result for Disease.
	'''
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if 'Uniprotkb' in request.GET and request.GET['Uniprotkb']:
		uniprotkb=request.GET['Uniprotkb']
		uniprotkb=uniprotkb.strip()
		diseasePresence=0
		diseaseUniProt=['NA']
		diseaseDisGeNet=['NA']
		protname=None
		geneName=None
		es.indices.refresh(index="mrmassaydb-index")
		query={"query": {
			"bool": {
				"must": [
					{"match": {"UniProtKB Accession": uniprotkb}},
					{"match": {"UniprotKb entry status": "Yes"}}
				]
			}
		}
		}

		res=helpers.scan(client=es,scroll='2m',index="mrmassaydb-index", doc_type="mrmassaydb-type",query=query,request_timeout=30)
		for hit in res:
			jdic=hit["_source"]
			jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
			protname=str(jdic["Protein"]).strip()
			genename=str(jdic["Gene"]).strip()
			jdic["sel"] =""
			if len(jdic["Disease Name"].strip()) >0:
				diseasePresence=1

			if len(jdic["DisGen DiseaseData"].strip()) > 0 and jdic["DisGen DiseaseData"].strip() != "NA":
				jdic["DisGen DiseaseData URL"]=jdic["DisGen DiseaseData URL"].replace('</a>"','</a>')
				jdic["DisGen DiseaseData URL"]=jdic["DisGen DiseaseData URL"].replace('href=\\','href=')
				jdic["DisGen DiseaseData URL"]=jdic["DisGen DiseaseData URL"].replace('\\">','">')
				diseaseDisGeNet=jdic["DisGen DiseaseData URL"].strip().split('|')
			else:
				diseaseDisGeNet=['NA']
			if len(jdic["UniProt DiseaseData"].strip()) > 0 and jdic["UniProt DiseaseData"].strip() != "NA":
				jdic["UniProt DiseaseData URL"]=jdic["UniProt DiseaseData URL"].replace('</a>"','</a>')
				jdic["UniProt DiseaseData URL"]=jdic["UniProt DiseaseData URL"].replace('href=\\','href=')
				jdic["UniProt DiseaseData URL"]=jdic["UniProt DiseaseData URL"].replace('\\">','">')
				diseaseUniProt=jdic["UniProt DiseaseData URL"].strip().split('|')				
			else:
				diseaseUniProt=['NA']
			break
		es.indices.refresh(index="mrmassaydb-index")

		return render(request, 'disease.html',{'uniprotkb':uniprotkb,\
			'diseasePresence':diseasePresence,'diseaseUniProtRaw':json.dumps(diseaseUniProt),\
			'diseaseDisGeNetRaw':json.dumps(diseaseDisGeNet),'protname':protname,\
			'gene':genename})

def goterm(request):
	'''
	This function will display result for Go term.
	'''
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if 'Uniprotkb' in request.GET and request.GET['Uniprotkb']:
		uniprotkb=request.GET['Uniprotkb']
		uniprotkb=uniprotkb.strip()
		contextgoterm =[]
		reachable=True
		if len(uniprotkb)>0:
			code=None
			gene=None
			protname=None
			es.indices.refresh(index="mrmassaydb-index")
			#build elasticsearch query based provided uniprotkb acc
			query={"query": {
				"bool": {
					"must": [
						{"match": {"UniProtKB Accession": uniprotkb}},
						{"match": {"UniprotKb entry status": "Yes"}}
					]
				}
			}
			}
			res = es.search(index="mrmassaydb-index", doc_type="mrmassaydb-type", body=query)

			foundHits=res["hits"]["total"]
			for hit in res['hits']['hits']:
				jdic=hit["_source"]
				jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
				jdic["sel"] =""
				gene=str(jdic["Gene"]).strip()
				protname=str(jdic["Protein"]).strip()
			es.indices.refresh(index="mrmassaydb-index")
			if foundHits >0:
				if '-' in uniprotkb:
					code=str(uniprotkb).split('-')[0]
				else:
					code=str(uniprotkb.strip())

				try:
					sleep(random.randint(5,10)) # delay between 5 and 10 sec
					requests.get("https://www.uniprot.org/", timeout=5) # check uniprot website live
					try:
						GrequestURL="https://www.uniprot.org/uniprot/"+str(code)+".xml"
						Gunifile=urllib.urlopen(GrequestURL)
						Gunidata= Gunifile.read()
						Gunifile.close()
						#reading uniprot XML data
						try:
							Gunidata=minidom.parseString(Gunidata)
							try:
								godata=(Gunidata.getElementsByTagName('dbReference'))
								for gItem in godata:
									if (gItem.attributes['type'].value).upper() == 'GO':
										try:
											gonamedetails=(str(gItem.getElementsByTagName('property')[0].attributes['value'].value).strip()).split(':')[1] #get go term name
											gotermdetails=(str(gItem.getElementsByTagName('property')[0].attributes['value'].value).strip()).split(':')[0] #get go term
											goid=str(gItem.attributes['id'].value).strip() # get GO ID
											tempGolist=[]
											tempGolist.append(uniprotkb.upper())
											tempGolist.append(protname)
											tempGolist.append(gene)
											tempGolist.append(goid)
											tempGolist.append(gonamedetails)
											if gotermdetails.lower()=='p':
												tempGolist.append('Biological Process')
											if gotermdetails.lower()=='f':
												tempGolist.append('Molecular Function')
											if gotermdetails.lower()=='c':
												tempGolist.append('Cellular Component')

											if tempGolist not in contextgoterm:
												contextgoterm.append(tempGolist)

										except:
											pass
							except IndexError:
								pass
						except ExpatError:
							pass
					except IOError:
						pass
				except ConnectionError as e:
					reachable=False

				if reachable:
					return render(request, 'goterm.html',  {'foundHits': foundHits,'contextgoterminfo': contextgoterm,'uniprotkb':uniprotkb,'reachable':reachable})
				else:
					return render(request, 'goterm.html',  {'reachable':reachable,'foundHits': foundHits})
			else:
				return render(request, 'goterm.html',  {'foundHits': foundHits})

def transition(request):
	"""
	This function is searching information for tranisition, based on intrument type
	in database.
	"""
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if 'uniprotkb' in request.GET and request.GET['uniprotkb'] and 'pepseq' in request.GET and request.GET['pepseq'] and 'modpepseq' in request.GET and request.GET['modpepseq']:
		uniprotkb=request.GET['uniprotkb']
		uniprotkb= uniprotkb.strip()
		pepseq= request.GET['pepseq']
		pepseq= pepseq.strip()
		modpepseq= request.GET['modpepseq']
		modpepseq= modpepseq.strip()
		modpepseq=modpepseq.replace(" ","")
		es.indices.refresh(index="mrmassaydb-index")
		#build elasticsearch query based on peptide seq and uniprotkb acc
		query={"query": {
			"bool": {
				"must": [
					{"match": {"UniProtKB Accession": uniprotkb}},
					{"match": {"Peptide Sequence": pepseq}},
					{"match": {"UniprotKb entry status": "Yes"}}
				]
			}
		}
		}
		res = es.search(index="mrmassaydb-index", doc_type="mrmassaydb-type", body=query)

		transdic={}
		listOfTransCol=["PeptideTracker Transition","Passel Transition","SRMAtlas Transition","CPTAC Transitions","Panoramaweb Transition"]
		foundHits=res["hits"]["total"]
		protList=[]
		summarydata=[]
		for hit in res['hits']['hits']:
			jdic=hit["_source"]
			jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
			jdic["sel"] =""
			if modpepseq.upper() == (str(jdic["Modified Peptide Sequence"]).strip().upper()).replace("+",""):
				tempsummarydata=str(jdic["Summary Transition"]).strip()
				for i in tempsummarydata.split('.'):
					if 'instrument' not in i.lower():
						summarydata.append(i.split('|'))
				if str(jdic["UniProtKB Accession"]).strip() == uniprotkb:
					protList.append(str(jdic["UniProtKB Accession"]).strip())
					protList.append(str(jdic["Protein"]).strip())
					protList.append(str(jdic["Peptide Sequence"]).strip())
					protList.append(str(jdic["Modified Peptide Sequence"]).strip().upper())
				for trancol in listOfTransCol:
					if (str(jdic[trancol]).strip()) >0 and (str(jdic[trancol]).strip()).lower() !='na':
						transdata=str(jdic[trancol]).strip()
						transinfo=transdata.split(';')
						temptransdic={}
						for titem in transinfo:
							subtransinfo=titem.split(',')
							for stitem in subtransinfo[1:]:
								if 'instrument' not in stitem.lower():
									tinfo =stitem.split('|')
									#formatting SRMAtlas transition data
									if trancol =="SRMAtlas Transition":
										subintrument=tinfo[:-5]
										for y in subintrument:
											subtinfo=tinfo[-5:]
											subtinfo.insert(0,y)
											if temptransdic.has_key(subtinfo[0].strip()):
												temptransdic[subtinfo[0].strip()].append(subtinfo)
											else:
												temptransdic[subtinfo[0].strip()]=[subtinfo]
									else:
										if temptransdic.has_key(tinfo[0].strip()):
											temptransdic[tinfo[0].strip()].append(tinfo)
										else:
											temptransdic[tinfo[0].strip()]=[tinfo]
						if len(temptransdic)>0:
							transdic[trancol]=temptransdic
		es.indices.refresh(index="mrmassaydb-index")
		if foundHits >0:
			contextindex={'transdic':transdic,'foundHits':foundHits,'protList':protList,'summarydata':summarydata}
			return render(request, 'transition.html', contextindex)
		else:
			return render(request,'transition.html',{'foundHits':foundHits})

def concentration(request):
	'''
	This function will display result for KEGG pathways and STRING PPI.
	'''
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if 'uniprotkb' in request.GET and request.GET['uniprotkb'] and \
	'pepseq' in request.GET and request.GET['pepseq'] and 'peptrackid' in request.GET and request.GET['peptrackid']:

		uniprotkb=request.GET['uniprotkb']
		uniprotkb=uniprotkb.strip()
		pepseq= request.GET['pepseq']
		pepseq= pepseq.strip()
		peptrackid=request.GET['peptrackid']
		peptrackid=peptrackid.strip()
		es.indices.refresh(index="mrmassaydb-index")
		query={"query": {
			"bool": {
				"must": [
					{"match": {"UniProtKB Accession": uniprotkb}},
					{"match": {"Peptide Sequence": pepseq}},
					{"match": {"UniprotKb entry status": "Yes"}}
				]
			}
		}
		}
		res = es.search(index="mrmassaydb-index", doc_type="mrmassaydb-type", body=query)
		conclist=[]
		foundHits=res["hits"]["total"]
		protList=[]
		for hit in res['hits']['hits']:
			jdic=hit["_source"]
			jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
			jdic["sel"] =""
			if str(jdic['UniProtKB Accession']).strip() == uniprotkb:
				protList.append(str(jdic['UniProtKB Accession']).strip())
				protList.append(str(jdic["Protein"]).strip())
				protList.append(str(jdic['Peptide Sequence']).strip())
			if (str(jdic["Concentration"]).strip()) >0 and (str(jdic["Concentration"]).strip()).lower() !='no' and str(jdic['PeptideTracker ID']).strip() == peptrackid:
				try:
					coninfo=(jdic['Concentration'].strip()).split(';')
					tempconinfo=[]
					if len(coninfo)>0:
						for c in coninfo:
							subconinfo=c.split('|')
							tempmatrix=str(subconinfo[4]).split()[-1]
							subconinfo[0]=subconinfo[0].replace('tissue',tempmatrix)
							c='|'.join(subconinfo)
							tempconinfo.append(c)
					jdic['Concentration Range']=';'.join(tempconinfo)
				except KeyError:
					pass
				concdata=str(jdic["Concentration"]).strip()
				concdata=concdata.replace('fmol/','fmol/µ')
				concinfo=concdata.split(';')
				for citem in concinfo:
					cinfo =citem.split('|')
					conclist.append(cinfo)
		es.indices.refresh(index="mrmassaydb-index")
		if foundHits >0:
			contextindex={'conclist':conclist,'foundHits':foundHits,'protList':protList}
			return render(request, 'concentration.html',contextindex)
		else:
			return render(request, 'concentration.html',{'foundHits':foundHits})

def help(request):
	"""
	This is function will generate data for help page inlcuding updating data
	"""
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	updatedstatresult=updatedstat
	countaxis=len(updatedstatresult)/10
	if countaxis==0:
		countaxis=len(updatedstatresult)
	return render(request, 'help.html', {'updatedstatresult':updatedstatresult,'countaxis':countaxis})

def contact(request):
	"""
	This is contact form where 
	Admin will recieve email from contact page.
	"""
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	title="For any inquiry, please contact us using following contact form or email us at bioinformatics@proteincentre.com!"
	form = ContactForm(request.POST or None)
	contextcontact={"title":title, "form":form}
	if form.is_valid():
		form_email = form.cleaned_data.get("email") # get email address
		form_email=form_email.strip()
		form_full_name = form.cleaned_data.get("full_name") # get name
		form_full_name=form_full_name.strip()
		form_message = form.cleaned_data.get("message") # get message
		form_message=form_message.strip()
		subject='Site contact form for MRMAssayDB'
		from_email=settings.EMAIL_HOST_USER
		to_email=[from_email, 'Pallab@proteincentre.com','Yassene@proteincentre.com','bioinformatics@proteincentre.com']
		contact_message="%s: %s via %s"%(
			form_full_name,
			form_message,
			form_email)
		send_mail(subject,contact_message,from_email,to_email,fail_silently=True) # sent email

		contextcontact={"title":"Thank you for your enquiry","hide":True}
	return render(request, 'contact.html', contextcontact)

@login_required
def adminsite(request):
	"""
	This is admin site where only djano superuser
	can access informations.
	User ip details and statistics
	"""
	if request.user.is_superuser and request.user.is_active:
		contextadminip={}
		if 'searchtype' in request.GET and request.GET['searchtype']:
			searchtype= request.GET['searchtype'] # adminsite for accessing informations
			searchtype= searchtype.strip()
			#generate datatable for user ip address
			if searchtype=='userip':
				ipinfo=IpAddressInformation.objects.all().values() # get all information from ip database
				for ipitem in ipinfo:
					ipaddess=ipitem['ip_address']
					accesstime=ipitem['access_date']
					if contextadminip.has_key(str(ipaddess)):
						contextadminip[str(ipaddess)].append(accesstime)
					else:
						contextadminip[str(ipaddess)]= [accesstime]



		return render(request, 'adminsite.html', {'contextadminip':contextadminip})
	else:
		return HttpResponse("You don't have permission to access this page.Please contact site administrator for more details")

def restfulapisearch(request):
	"""
	This function is searching results, based on given multi search parameters
	in database.
	"""
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if 'UniProtKB Accession' in request.GET and request.GET['UniProtKB Accession'] or \
	'Protein' in request.GET and request.GET['Protein'] or \
	'Gene' in request.GET and request.GET['Gene'] or \
	'Organism' in request.GET and request.GET['Organism'] or \
	'Organism ID' in request.GET and request.GET['Organism ID'] or \
	'SubCellular' in request.GET and request.GET['SubCellular']  or \
	'Peptide Sequence' in request.GET and request.GET['Peptide Sequence']  or \
	'Pathway Name' in request.GET and request.GET['Pathway Name']  or \
	'Disease Name' in request.GET and request.GET['Disease Name']  or \
	'Go ID' in request.GET and request.GET['Go ID']  or \
	'Go Name' in request.GET and request.GET['Go Name'] or \
	'Go Term' in request.GET and request.GET['Go Term'] or 'Drug Bank' in request.GET and request.GET['Drug Bank'] or \
	'AssayFdaApproveMark' in request.GET and request.GET['AssayFdaApproveMark']:
		useruniprotkb =""
		userprotein =""
		usergeneid =""
		userorg=""
		userorgid=""
		usersubcell =""
		userpepseq =""
		userpathway =""
		userdis =""
		usergoid =""
		usergotn =""
		usergot=""
		userdrug=""
		userassayfdaapprovemark=""
		try:
			useruniprotkb = request.GET["UniProtKB Accession"]
		except MultiValueDictKeyError:
			pass
		if '|' in useruniprotkb:
			useruniprotkb=(useruniprotkb.strip()).split('|')
		else:
			useruniprotkb=(useruniprotkb.strip()).split('\n')
		useruniprotkb=[(item.strip()).lower() for item in useruniprotkb]
		useruniprotkb=map(str, useruniprotkb)
		useruniprotkb=filter(None, useruniprotkb)

		try:
			userprotein = request.GET["Protein"]
		except MultiValueDictKeyError:
			pass
		if '|' in userprotein:
			userprotein=(userprotein.strip()).split('|')
		else:
			userprotein=(userprotein.strip()).split('\\n')
		userprotein=[(item.strip()).lower() for item in userprotein]
		userprotein=map(str, userprotein)
		userprotein=filter(None, userprotein)

		try:
			usergeneid = request.GET["Gene"]
		except MultiValueDictKeyError:
			pass
		if '|' in usergeneid:
			usergeneid=(usergeneid.strip()).split('|')
		else:
			usergeneid=(usergeneid.strip()).split('\\n')
		usergeneid=[(item.strip()).lower() for item in usergeneid]
		usergeneid=map(str, usergeneid)
		usergeneid=filter(None, usergeneid)

		try:
			userorg = request.GET["Organism"]
		except MultiValueDictKeyError:
			pass
		if '|' in userorg:
			userorg=(userorg.strip()).split('|')
		else:
			userorg=(userorg.strip()).split('\\n')
		userorg=[(item.strip()).lower() for item in userorg]
		userorg=map(str, userorg)
		userorg=filter(None, userorg)

		try:
			userorgid = request.GET["Organism ID"]
		except MultiValueDictKeyError:
			pass
		if '|' in userorgid:
			userorgid=(userorgid.strip()).split('|')
		else:
			userorgid=(userorgid.strip()).split('\\n')
		userorgid=[(item.strip()).lower() for item in userorgid]
		userorgid=map(str, userorgid)
		userorgid=filter(None, userorgid)

		try:
			usersubcell = request.GET["SubCellular"]
		except MultiValueDictKeyError:
			pass
		if '|' in usersubcell:
			usersubcell=(usersubcell.strip()).split('|')
		else:
			usersubcell=(usersubcell.strip()).split('\\n')
		usersubcell=[(item.strip()).lower() for item in usersubcell]
		usersubcell=map(str, usersubcell)
		usersubcell=filter(None, usersubcell)

		try:
			userpepseq = request.GET["Peptide Sequence"]
		except MultiValueDictKeyError:
			pass
		if '|' in userpepseq:
			userpepseq=(userpepseq.strip()).split('|')
		else:
			userpepseq=(userpepseq.strip()).split('\\n')
		userpepseq=[(item.strip()).lower() for item in userpepseq]
		userpepseq=map(str, userpepseq)
		userpepseq=filter(None, userpepseq)

		try:
			userpathway = request.GET["Pathway Name"]
		except MultiValueDictKeyError:
			pass
		if '|' in userpathway:
			userpathway=(userpathway.strip()).split('|')
		else:
			userpathway=(userpathway.strip()).split('\\n')
		userpathway=[(item.strip()).lower() for item in userpathway]
		userpathway=map(str, userpathway)
		userpathway=filter(None, userpathway)

		try:
			userdis = request.GET["Disease Name"]
		except MultiValueDictKeyError:
			pass
		if '|' in userdis:
			userdis=(userdis.strip()).split('|')
		else:
			userdis=(userdis.strip()).split('\\n')
		userdis=[(item.strip()).lower() for item in userdis]
		userdis=map(str, userdis)
		userdis=filter(None, userdis)

		try:
			usergoid = request.GET["Go ID"]
		except MultiValueDictKeyError:
			pass
		if '|' in usergoid:
			usergoid=(usergoid.strip()).split('|')
		else:
			usergoid=(usergoid.strip()).split('\\n')
		usergoid=[(item.strip()).lower() for item in usergoid]
		usergoid=map(str, usergoid)
		usergoid=filter(None, usergoid)

		try:
			usergotn = request.GET["Go Name"]
		except MultiValueDictKeyError:
			pass
		if '|' in usergotn:
			usergotn=(usergotn.strip()).split('|')
		else:
			usergotn=(usergotn.strip()).split('\\n')
		usergotn=[(item.strip()).lower() for item in usergotn]
		usergotn=map(str, usergotn)
		usergotn=filter(None, usergotn)

		try:
			usergot = request.GET["Go Term"]
		except MultiValueDictKeyError:
			pass
		if '|' in usergot:
			usergot=(usergot.strip()).split('|')
		else:
			usergot=(usergot.strip()).split('\\n')
		usergot=[(item.strip()).lower() for item in usergot]
		usergot=map(str, usergot)
		usergot=filter(None, usergot)

		try:
			userdrug = request.GET["Drug Bank"]
		except MultiValueDictKeyError:
			pass
		if '|' in userdrug:
			userdrug=(userdrug.strip()).split('|')
		else:
			userdrug=(userdrug.strip()).split('\\n')
		userdrug=[(item.strip()).lower() for item in userdrug]
		userdrug=map(str, userdrug)
		userdrug=filter(None, userdrug)

		try:
			userassayfdaapprovemark = request.GET["AssayFdaApproveMark"]
		except MultiValueDictKeyError:
			pass
		if '|' in userassayfdaapprovemark:
			userassayfdaapprovemark=(userassayfdaapprovemark.strip()).split('|')
			userassayfdaapprovemark=list(set(userassayfdaapprovemark))
		else:
			userassayfdaapprovemark=(userassayfdaapprovemark.strip()).split('\\n')
			userassayfdaapprovemark=list(set(userassayfdaapprovemark))
		userassayfdaapprovemark=[(item.strip()).lower() for item in userassayfdaapprovemark]
		userassayfdaapprovemark=map(str, userassayfdaapprovemark)
		userassayfdaapprovemark=filter(None, userassayfdaapprovemark)

		spquerylist =[]
		searchtermlist=[]

		if len(useruniprotkb) >0:
			shouldlist=[]
			for x in useruniprotkb:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["UniProtKB Accession.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userprotein)> 0:
			shouldlist=[]
			for x in userprotein:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Protein.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(usergeneid) >0:
			shouldlist=[]
			for x in usergeneid:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Gene.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userorg) > 0:
			shouldlist=[]
			for x in userorg:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Organism.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userorgid) > 0:
			shouldlist=[]
			for x in userorgid:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Organism ID.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(usersubcell) >0:
			shouldlist=[]
			for x in usersubcell:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["SubCellular.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userpepseq) >0:
			shouldlist=[]
			for x in userpepseq:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Peptide Sequence.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userpathway) >0:
			shouldlist=[]
			for x in userpathway:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Pathway Name.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userdis) >0:
			shouldlist=[]
			for x in userdis:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Disease Name.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(usergoid) >0:
			sdict={}
			sdict["Go ID.ngram"]=[i.split(' ')[0] for i in usergoid]
			tdict={}
			tdict["terms"]=sdict
			searchtermlist.append(tdict)
			shouldlist=[]
			for x in usergoid:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Go ID.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(usergotn) >0:
			shouldlist=[]
			for x in usergotn:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Go Name.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(usergot) > 0:
			shouldlist=[]
			for x in usergot:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Go Term.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userdrug) > 0:
			shouldlist=[]
			for x in userdrug:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Drug Bank.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)

		if len(userassayfdaapprovemark) > 0:
			shouldlist=[]
			for x in userassayfdaapprovemark:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Assays for FDA approved Marker.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)

		if len(searchtermlist)>0:
			es.indices.refresh(index="mrmassaydb-index")

			query={
				"query": {
					"bool": {
						"must":searchtermlist
					}
				}
			}
			jsonfilename=names.get_first_name()+'_restapi_search.json'
			jsonfilepath=os.path.join(settings.BASE_DIR, 'resultFile', 'jsonData','resultJson', 'restapisearch', jsonfilename)
			jsonfileoutput= open(jsonfilepath,'w')
			jfinaldata=[]
			res=helpers.scan(client=es,scroll='2m',index="mrmassaydb-index", doc_type="mrmassaydb-type",query=query,request_timeout=30)
			jfinaldata=[]
			for i in res:
				jdic=i['_source']
				jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
				if jdic["UniprotKb entry status"] =="Yes" and jdic['UniProtKB Accession'] !='502':
					jdic["sel"] =""
					jdic["Drug Bank"]=jdic["Drug Bank"].replace('\\','')
					jdic["Drug Bank"]=jdic["Drug Bank"].replace('<br>','|')
					jdic["SRMAtlas URL"]=jdic["SRMAtlas URL"].replace('\\','')
					jdic["Passel URL"]=jdic["Passel URL"].replace('\\','')
					jdic["CPTAC URL"]=jdic["CPTAC URL"].replace('\\','')
					jdic["Panoramaweb URL"]=jdic["Panoramaweb URL"].replace('\\','')
					jdic["PeptideTracker URL"]=jdic["PeptideTracker URL"].replace('\\','')
					#if jdic["Pathway Name"].lower() !='na':
					#	jdic["Pathway Name"]=re.sub(r"(\w)([A-Z])",r"\1|\2",jdic["Pathway Name"])
					jdic["Mean Concentration"] =jdic["Mean Concentration"].replace('fmol/','fmol/µ')
					jdic["Concentration"] =str(jdic["Concentration"].replace('fmol/','fmol/µ'))
					jfinaldata.append(jdic)

			foundHits=len(jfinaldata)
			es.indices.refresh(index="mrmassaydb-index")

			json.dump(jfinaldata,jsonfileoutput)
			jsonfileoutput.close()
			df=pd.read_json(jsonfilepath)
			donwloadfilename ="MRMAssayDBResult.csv"
			response = HttpResponse(content_type='text/csv')
			response['Content-Disposition'] = 'attachment; filename=%s'%donwloadfilename
			writer = csv.writer(response)
			if foundHits >0:
				df.to_csv(path_or_buf=response,index=False, columns=downloadcolname, encoding='utf-8')
			return response

		else:
			return render(request, 'index.html', {'error': True})
	else:
		return render(request, 'index.html', {'error': True})


def hyperlink_search(request):
	"""
	This function is searching results, based on given multi search parameters
	in database.
	"""
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if 'UniProtKB Accession' in request.GET and request.GET['UniProtKB Accession'] or \
	'Protein' in request.GET and request.GET['Protein'] or \
	'Gene' in request.GET and request.GET['Gene'] or \
	'Organism' in request.GET and request.GET['Organism'] or \
	'Organismid' in request.GET and request.GET['Organismid'] or \
	'SubCellular' in request.GET and request.GET['SubCellular']  or \
	'Peptide Sequence' in request.GET and request.GET['Peptide Sequence']  or \
	'Pathway Name' in request.GET and request.GET['Pathway Name']  or \
	'Disease Name' in request.GET and request.GET['Disease Name']  or \
	'Go ID' in request.GET and request.GET['Go ID']  or \
	'Go Name' in request.GET and request.GET['Go Name'] or \
	'Go Term' in request.GET and request.GET['Go Term'] or \
	'AssayFdaApproveMark' in request.GET and request.GET['AssayFdaApproveMark']:
		useruniprotkb =""
		userprotein =""
		usergeneid =""
		userorg=""
		userorgid=""
		usersubcell =""
		userpepseq =""
		userpathway =""
		userdis =""
		usergoid =""
		usergotn =""
		usergot=""
		userassayfdaapprovemark=""
		finalsearhdata=''
		try:
			useruniprotkb = request.GET["UniProtKB Accession"]
		except MultiValueDictKeyError:
			pass
		if '|' in useruniprotkb:
			useruniprotkb=(useruniprotkb.strip()).split('|')
		else:
			useruniprotkb=(useruniprotkb.strip()).split('\\n')
		useruniprotkb=[(item.strip()).lower() for item in useruniprotkb]
		useruniprotkb=map(str, useruniprotkb)
		useruniprotkb=filter(None, useruniprotkb)

		try:
			userprotein = request.GET["Protein"]
		except MultiValueDictKeyError:
			pass
		if '|' in userprotein:
			userprotein=(userprotein.strip()).split('|')
		else:
			userprotein=(userprotein.strip()).split('\\n')
		userprotein=[(item.strip()).lower() for item in userprotein]
		userprotein=map(str, userprotein)
		userprotein=filter(None, userprotein)

		try:
			usergeneid = request.GET["Gene"]
		except MultiValueDictKeyError:
			pass
		if '|' in usergeneid:
			usergeneid=(usergeneid.strip()).split('|')
		else:
			usergeneid=(usergeneid.strip()).split('\\n')
		usergeneid=[(item.strip()).lower() for item in usergeneid]
		usergeneid=map(str, usergeneid)
		usergeneid=filter(None, usergeneid)

		try:
			userorg = request.GET["Organism"]
		except MultiValueDictKeyError:
			pass
		if '|' in userorg:
			userorg=(userorg.strip()).split('|')
		else:
			userorg=(userorg.strip()).split('\\n')
		userorg=[(item.strip()).lower() for item in userorg]
		userorg=map(str, userorg)
		userorg=filter(None, userorg)

		try:
			userorgid = request.GET["Organismid"]
		except MultiValueDictKeyError:
			pass
		if '|' in userorgid:
			userorgid=(userorgid.strip()).split('|')
		else:
			userorgid=(userorgid.strip()).split('\\n')
		userorgid=[(item.strip()).lower() for item in userorgid]
		userorgid=map(str, userorgid)
		userorgid=filter(None, userorgid)

		try:
			usersubcell = request.GET["SubCellular"]
		except MultiValueDictKeyError:
			pass
		if '|' in usersubcell:
			usersubcell=(usersubcell.strip()).split('|')
		else:
			usersubcell=(usersubcell.strip()).split('\\n')
		usersubcell=[(item.strip()).lower() for item in usersubcell]
		usersubcell=map(str, usersubcell)
		usersubcell=filter(None, usersubcell)

		try:
			userpepseq = request.GET["Peptide Sequence"]
		except MultiValueDictKeyError:
			pass
		if '|' in userpepseq:
			userpepseq=(userpepseq.strip()).split('|')
		else:
			userpepseq=(userpepseq.strip()).split('\\n')
		userpepseq=[(item.strip()).lower() for item in userpepseq]
		userpepseq=map(str, userpepseq)
		userpepseq=filter(None, userpepseq)

		try:
			userpathway = request.GET["Pathway Name"]
		except MultiValueDictKeyError:
			pass
		if '|' in userpathway:
			userpathway=(userpathway.strip()).split('|')
		else:
			userpathway=(userpathway.strip()).split('\\n')
		userpathway=[(item.strip()).lower() for item in userpathway]
		userpathway=map(str, userpathway)
		userpathway=filter(None, userpathway)

		try:
			userdis = request.GET["Disease Name"]
		except MultiValueDictKeyError:
			pass
		if '|' in userdis:
			userdis=(userdis.strip()).split('|')
		else:
			userdis=(userdis.strip()).split('\\n')
		userdis=[(item.strip()).lower() for item in userdis]
		userdis=map(str, userdis)
		userdis=filter(None, userdis)

		try:
			usergoid = request.GET["Go ID"]
		except MultiValueDictKeyError:
			pass
		if '|' in usergoid:
			usergoid=(usergoid.strip()).split('|')
		else:
			usergoid=(usergoid.strip()).split('\\n')
		usergoid=[(item.strip()).lower() for item in usergoid]
		usergoid=map(str, usergoid)
		usergoid=filter(None, usergoid)

		try:
			usergotn = request.GET["Go Name"]
		except MultiValueDictKeyError:
			pass
		if '|' in usergotn:
			usergotn=(usergotn.strip()).split('|')
		else:
			usergotn=(usergotn.strip()).split('\\n')
		usergotn=[(item.strip()).lower() for item in usergotn]
		usergotn=map(str, usergotn)
		usergotn=filter(None, usergotn)

		try:
			usergot = request.GET["Go Term"]
		except MultiValueDictKeyError:
			pass
		if '|' in usergot:
			usergot=(usergot.strip()).split('|')
		else:
			usergot=(usergot.strip()).split('\\n')
		usergot=[(item.strip()).lower() for item in usergot]
		usergot=map(str, usergot)
		usergot=filter(None, usergot)

		try:
			userassayfdaapprovemark = request.GET["AssayFdaApproveMark"]
		except MultiValueDictKeyError:
			pass
		if '|' in userassayfdaapprovemark:
			userassayfdaapprovemark=(userassayfdaapprovemark.strip()).split('|')
			userassayfdaapprovemark=list(set(userassayfdaapprovemark))
		else:
			userassayfdaapprovemark=(userassayfdaapprovemark.strip()).split('\\n')
			userassayfdaapprovemark=list(set(userassayfdaapprovemark))
		userassayfdaapprovemark=[(item.strip()).lower() for item in userassayfdaapprovemark]
		userassayfdaapprovemark=map(str, userassayfdaapprovemark)
		userassayfdaapprovemark=filter(None, userassayfdaapprovemark)

		spquerylist =[]
		searchtermlist=[]

		if len(useruniprotkb) >0:
			finalsearhdata+='UniProtKB Accession:'+';'.join(useruniprotkb)+' '
			shouldlist=[]
			for x in useruniprotkb:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["UniProtKB Accession.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userprotein)> 0:
			finalsearhdata+='Protein:'+';'.join(userprotein)+' '
			shouldlist=[]
			for x in userprotein:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Protein.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(usergeneid) >0:
			finalsearhdata+='Gene:'+';'.join(usergeneid)+' '
			shouldlist=[]
			for x in usergeneid:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Gene.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userorg) > 0:
			finalsearhdata+='Organism:'+';'.join(userorg)+' '
			shouldlist=[]
			for x in userorg:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Organism.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userorgid) > 0:
			finalsearhdata+='Organism ID:'+';'.join(userorgid)+' '
			shouldlist=[]
			for x in userorgid:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Organism ID.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(usersubcell) >0:
			finalsearhdata+='SubCellular:'+';'.join(usersubcell)+' '
			shouldlist=[]
			for x in usersubcell:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["SubCellular.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userpepseq) >0:
			finalsearhdata+='Peptide Sequence:'+';'.join(userpepseq)+' '
			shouldlist=[]
			for x in userpepseq:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Peptide Sequence.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userpathway) >0:
			finalsearhdata+='Pathway Name:'+';'.join(userpathway)+' '
			shouldlist=[]
			for x in userpathway:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Pathway Name.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(userdis) >0:
			finalsearhdata+='Disease Name:'+';'.join(userdis)+' '
			shouldlist=[]
			for x in userdis:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Disease Name.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(usergoid) >0:
			finalsearhdata+='Go ID:'+';'.join(usergoid)+' '
			sdict={}
			sdict["Go ID.ngram"]=[i.split(' ')[0] for i in usergoid]
			tdict={}
			tdict["terms"]=sdict
			searchtermlist.append(tdict)
			shouldlist=[]
			for x in usergoid:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Go ID.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]+={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(usergotn) >0:
			finalsearhdata+='Go Name:'+';'.join(usergotn)+' '
			shouldlist=[]
			for x in usergotn:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Go Name.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)
		if len(usergot) > 0:
			finalsearhdata+='Go Term:'+';'.join(usergot)+' '
			shouldlist=[]
			for x in usergot:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Go Term.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)

		if len(userassayfdaapprovemark) > 0:
			finalsearhdata+='Assays for FDA approved Marker::'+';'.join(userassayfdaapprovemark)+' '
			shouldlist=[]
			for x in userassayfdaapprovemark:
				tempquery={
							"multi_match":{
								"query":x.strip(),
								"type":"best_fields",
								"fields":["Assays for FDA approved Marker.ngram"],
								"minimum_should_match":"100%"
							}
						}
				shouldlist.append(tempquery)
			booldic={}
			booldic["bool"]={"should":shouldlist,"minimum_should_match": 1}
			searchtermlist.append(booldic)

		if len(searchtermlist)>0:
			es.indices.refresh(index="mrmassaydb-index")

			query={
				"query": {
					"bool": {
						"must":searchtermlist
					}
				}
			}
			nameFIle=names.get_first_name()
			jsonfilename=nameFIle+'_advance_search.json'
			jsonfilepath=os.path.join(settings.BASE_DIR, 'resultFile', 'jsonData','resultJson', 'adavancesearch', 'results', jsonfilename)
			jsonfileoutput= open(jsonfilepath,'w')
			jfinaldata=[]
			res=helpers.scan(client=es,scroll='2m',index="mrmassaydb-index", doc_type="mrmassaydb-type",query=query,request_timeout=30)
			jfinaldata=[]
			for i in res:
				jdic=i['_source']
				jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
				if jdic["UniprotKb entry status"] =="Yes" and jdic['UniProtKB Accession'] !='502':
					jdic["sel"] =""
					jdic["PPI"] ="View"
					jdic["Drug Bank"]=jdic["Drug Bank"].replace('\\','')
					jdic["Drug Bank"]=jdic["Drug Bank"].replace('<br>','|')
					jdic["SRMAtlas URL"]=jdic["SRMAtlas URL"].replace('\\','')
					jdic["Passel URL"]=jdic["Passel URL"].replace('\\','')
					jdic["CPTAC URL"]=jdic["CPTAC URL"].replace('\\','')
					jdic["Panoramaweb URL"]=jdic["Panoramaweb URL"].replace('\\','')
					jdic["PeptideTracker URL"]=jdic["PeptideTracker URL"].replace('\\','')
					#if jdic["Pathway Name"].lower() !='na':
					#	jdic["Pathway Name"]=re.sub(r"(\w)([A-Z])",r"\1|\2",jdic["Pathway Name"])
					jdic["Mean Concentration"] =jdic["Mean Concentration"].replace('fmol/','fmol/µ')
					jdic["Concentration"] =jdic["Concentration"].replace('fmol/','fmol/µ')					
					jfinaldata.append(jdic)

			foundHits=len(jfinaldata)
			json.dump(jfinaldata[:10000],jsonfileoutput)
			jsonfileoutput.close()

			if foundHits >0:
				statsummary=summaryStatcal(jfinaldata)
				pathwaychart=statsummary['pathwaychart']
				pathwaychart=[i[:2] for i in pathwaychart]
				specieslist=statsummary['specieslist']
				totallist=statsummary['total']
				subcell=statsummary['subcell']
				godic=statsummary['godic']
				jvennprot=statsummary['jevennstat'][0]
				jvennpep=statsummary['jevennstat'][1]
				mrmdatabase=statsummary['jevennstat'][2]
				sortedgodic=OrderedDict(sorted(godic.items(), key=lambda t: t[1]))
				updatedgodic=dict(list(sortedgodic.items())[:10])
				pepseqdataseries=ast.literal_eval(json.dumps(statsummary['pepseqdataseries']))
				prodataseries=statsummary['prodataseries']
				unqisostat=statsummary['unqisostat']
				jsonfilepathStat=os.path.join(settings.BASE_DIR, 'resultFile', 'jsonData','resultJson', 'adavancesearch', 'statsummary', jsonfilename)
				jsonfileoutputStat= open(jsonfilepathStat,'w')
				json.dumps(statsummary,jsonfileoutputStat)
				jsonfileoutputStat.close()
				urlname="'/resultFile/jsonData/resultJson/adavancesearch/results/"+jsonfilename+"'"
				contextindex={
					"filename":urlname,"colname":json.dumps(colname),
					'query': finalsearhdata,'foundHits':foundHits,
					'pathwaychart':pathwaychart[:11],'specieslist':specieslist,
					'totallist':totallist,'subcell':subcell,
					'updatedgodic':updatedgodic,'pepseqdataseries':pepseqdataseries,
					'prodataseries':prodataseries,'unqisostat':unqisostat,
					'jvennprot':json.dumps(jvennprot),'jvennpep':json.dumps(jvennpep),'jvennmrmdb':json.dumps(mrmdatabase)
					}
				return render(request,'resultform.html',contextindex)
			else:
				return render(request,'resultform.html',{'foundHits':foundHits})

def peptide_uniqueness(request):
	"""
	This function is searching results, based on given multi search parameters
	in database.
	"""
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)
	if 'Uniprotkb' in request.GET and request.GET['Uniprotkb'] and \
		'pepseq' in request.GET and request.GET['pepseq']  and \
		'fastafile' in request.GET and request.GET['fastafile']:

		useruniprotkb = request.GET["Uniprotkb"]
		userpepseq = request.GET["pepseq"]
		userfastafilename = request.GET["fastafile"]
		useruniprotkb=str(useruniprotkb).strip()
		userpepseq=str(userpepseq).strip()
		userpepseq=userpepseq.upper()
		userfastafilename=str(userfastafilename).strip()
		reachable=True
		valid=False
		pepjsondatalist=[]
		es.indices.refresh(index="mrmassaydb-index")
		query={"query": {
			"bool": {
				"must": [
					{"match": {"UniProtKB Accession": useruniprotkb}},
					{"match": {"Peptide Sequence": userpepseq}}				]
			}
		}
		}
		res = es.search(index="mrmassaydb-index", doc_type="mrmassaydb-type", body=query)
		foundHits=res["hits"]["total"]
		originalpepunqstatus="Not unique"
		uniprotkbstatus="NA"
		for hit in res['hits']['hits']:
			jdic=hit["_source"]
			jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
			jdic["sel"] =""
			uniprotkbstatus=str(jdic["UniprotKb entry status"]).strip()
			if len(str(jdic["Unique in protein"]).strip())>0:
				originalpepunqstatus=str(jdic["Unique in protein"]).strip()
		es.indices.refresh(index="mrmassaydb-index")
		fastafilepath=os.path.join(settings.BASE_DIR, 'resultFile', 'fastaFIle', userfastafilename+'.fasta')
		if os.path.exists(fastafilepath) and foundHits>0:
			proseq=''
			try:
				sleep(random.randint(5,10))
				requests.get("https://www.uniprot.org/", timeout=5)
				unidatafasta = urllib.urlopen("https://www.uniprot.org/uniprot/" + useruniprotkb + ".fasta")
				for i in range(1):
					header=unidatafasta.next()
					if '>' in header[0]:
						valid=True
				if valid:
					for line in unidatafasta:
						proseq+=line.strip()
				unidatafasta.close()
			except ConnectionError as e:
				reachable=False

			if reachable:
				fastaseq=[]
				contextpepuser=[]
				peptide_pattern = re.compile(userpepseq,re.IGNORECASE)
				
				for match in re.finditer(peptide_pattern,proseq):
					matchpdbseqpos=range(int(match.start()),match.end())
					proseq=proseq.upper()
					tseq=proseq
					seqlen=len(proseq)
					proseqlist=list(proseq)
					for y in range(0,len(proseqlist)):
						if y in matchpdbseqpos:
							proseqlist[y]='<b><font color="red">'+proseqlist[y]+'</font></b>'
						if y%60 ==0 and y !=0:
							proseqlist[y]=proseqlist[y]+'<br/>'
					proseq="".join(proseqlist)
					contextpepuser.append(['<a target="_blank" href="https://www.uniprot.org/uniprot/'+useruniprotkb+'">'+useruniprotkb+'</a>',(int(match.start())+1),match.end(),matchpdbseqpos,tseq,seqlen,proseq,'MRMAssayDB'])

				for useq_record in SeqIO.parse(fastafilepath, 'fasta'):
					seqheader = useq_record.id
					sequniID = seqheader.split(' ')[0]
					sequniID=sequniID.replace('>','')
					tempseqs = str(useq_record.seq).strip()
					fastaseq.append(tempseqs.upper())
					for fmatch in re.finditer(peptide_pattern,tempseqs.upper()):
						fmatchpdbseqpos=range(int(fmatch.start()),fmatch.end())
						tempseqs=tempseqs.upper()
						tseq=tempseqs
						seqlen=len(tempseqs)
						tempseqslist=list(tempseqs)
						for y in range(0,len(tempseqslist)):
							if y in fmatchpdbseqpos:
								tempseqslist[y]='<b><font color="red">'+tempseqslist[y]+'</font></b>'
							if y%60 ==0 and y !=0:
								tempseqslist[y]=tempseqslist[y]+'<br/>'
						tempseqs="".join(tempseqslist)
						contextpepuser.append([sequniID,(int(fmatch.start())+1),fmatch.end(),fmatchpdbseqpos,tseq,seqlen,tempseqs,'User data'])

				useruniqstatus=''
				unqfastaseq=list(set(fastaseq))
				indices = [i for i, s in enumerate(fastaseq) if userpepseq.upper() in s]
				if len(indices)==0:
					useruniqstatus ="Not present"
				elif len(indices) > 1:
					useruniqstatus="Present but not unique"
				else:
					useruniqstatus="Present and unique"
				for x in contextpepuser:
					tempdatadic={}
					tempdatadic["proteinID"]=x[0]
					tempdatadic["peptideSequence"]=userpepseq
					tempdatadic["start"]=str(x[1])
					tempdatadic["end"]=str(x[2])
					tempdatadic["seqlen"]=str(x[-3])
					tempdatadic["seq"]=str(x[-4])
					if 'target' in x[0]:
						tempdatadic["peptideuniqueinprotein"]=originalpepunqstatus
					else:
						tempdatadic["peptideuniqueinprotein"]=useruniqstatus
					tempdatadic["datasource"]=str(x[-1])
					tempdatadic["highlightedpepseq"]=str(x[-2])
					pepjsondatalist.append(tempdatadic)
				pepunqdata={"data":pepjsondatalist}
				finalcontextpep={
								'pepunqdata':json.dumps(pepunqdata),'reachable':reachable
				}
				return render(request, 'peptideuniqueness.html', finalcontextpep)
			else:
				return render(request, 'peptideuniqueness.html', {'reachable':reachable})
		return render(request, 'peptideuniqueness.html', {'reachable':reachable})

def fdaAssay(request):
	"""
	This is for FDA Assay Tab.
	"""
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)

	contextres =[]
	#build elasticsearch query to search data
	query={"query": {
		"bool": {
			"must": [
				{"match": {"Assays for FDA approved Marker": "Yes"}},
				{"match": {"UniprotKb entry status": "Yes"}}
			]
		}
	}
	}
	#generate random file name to store search result in json format
	nameFIle=names.get_first_name()
	jsonfilename=nameFIle+'_basic_search_fda.json'
	jsonfilepath=os.path.join(settings.BASE_DIR, 'resultFile', 'jsonData','resultJson', 'basicsearch', 'results', jsonfilename)
	jsonfileoutput= open(jsonfilepath,'w')
	jfinaldata=[]
	es.indices.refresh(index="mrmassaydb-index")
	#elasticsearch will search data
	res=helpers.scan(client=es,scroll='2m',index="mrmassaydb-index", doc_type="mrmassaydb-type",query=query,request_timeout=30)
	jfinaldata=[]
	#if data is valid based on uniprotkb release then it will display
	for i in res:
		jdic=i['_source']
		jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
		if jdic["UniprotKb entry status"] =="Yes" and jdic['UniProtKB Accession'] !='502':
			jdic["PPI"] ="View"
			jdic["sel"] =""
			jdic["Drug Bank"]=jdic["Drug Bank"].replace('\\','')
			jdic["Drug Bank"]=jdic["Drug Bank"].replace('<br>','|')
			jdic["SRMAtlas URL"]=jdic["SRMAtlas URL"].replace('\\','')
			jdic["Passel URL"]=jdic["Passel URL"].replace('\\','')
			jdic["CPTAC URL"]=jdic["CPTAC URL"].replace('\\','')
			jdic["Panoramaweb URL"]=jdic["Panoramaweb URL"].replace('\\','')
			jdic["PeptideTracker URL"]=jdic["PeptideTracker URL"].replace('\\','')
			#if jdic["Pathway Name"].lower() !='na':
			#	jdic["Pathway Name"]=re.sub(r"(\w)([A-Z])",r"\1|\2",jdic["Pathway Name"])
			jdic["Mean Concentration"] =jdic["Mean Concentration"].replace('fmol/','fmol/µ')
			jdic["Concentration"] =jdic["Concentration"].replace('fmol/','fmol/µ')
			jfinaldata.append(jdic)
	es.indices.refresh(index="mrmassaydb-index")
	#checking any result generated by database
	foundHits=len(jfinaldata)
	#storing only 10000 rows in json format
	json.dump(jfinaldata[:10000],jsonfileoutput)
	jsonfileoutput.close()
	# if result found then do other job
	if foundHits >0:
		statsummary=summaryStatcal(jfinaldata) # sent data to this funcation for generating stat
		pathwaychart=statsummary['pathwaychart']
		pathwaychart=[i[:2] for i in pathwaychart]
		specieslist=statsummary['specieslist']
		totallist=statsummary['total']
		subcell=statsummary['subcell']
		godic=statsummary['godic']
		jvennprot=statsummary['jevennstat'][0]
		jvennpep=statsummary['jevennstat'][1]
		mrmdatabase=statsummary['jevennstat'][2]
		sortedgodic=OrderedDict(sorted(godic.items(), key=lambda t: t[1])) # sorting GO data
		updatedgodic=dict(list(sortedgodic.items()))
		pepseqdataseries=ast.literal_eval(json.dumps(statsummary['pepseqdataseries'])) #dumping data into json format
		prodataseries=statsummary['prodataseries']
		unqisostat=statsummary['unqisostat']
		jsonfilepathStat=os.path.join(settings.BASE_DIR, 'resultFile', 'jsonData','resultJson', 'basicsearch', 'statsummary', jsonfilename) #storing stat result in json format
		jsonfileoutputStat= open(jsonfilepathStat,'w')
		json.dump(statsummary,jsonfileoutputStat)
		jsonfileoutputStat.close()
		urlname="'/resultFile/jsonData/resultJson/basicsearch/results/"+jsonfilename+"'"

		contextindex={
			"filename":urlname,"colname":json.dumps(colname),'foundHits':foundHits,
			'pathwaychart':pathwaychart[:11],'specieslist':specieslist,
			'totallist':totallist,'subcell':subcell,
			'updatedgodic':updatedgodic,'pepseqdataseries':pepseqdataseries,
			'prodataseries':prodataseries,'unqisostat':unqisostat,
			'jvennprot':json.dumps(jvennprot),'jvennpep':json.dumps(jvennpep),'jvennmrmdb':json.dumps(mrmdatabase)
			}
		return render(request,'fdaAssay.html',contextindex)
	else:
		return render(request,'fdaAssay.html',{'foundHits':foundHits})

def covid19(request):
	"""
	This is for COVID 19 Tab.
	"""
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)

	contextres =[]
	#build elasticsearch query to search data
	query={"query": {
		"bool": {
			"must": [
				{"match": {"Associated with COVID-19": "Yes"}},
				{"match": {"UniprotKb entry status": "Yes"}}
			]
		}
	}
	}
	#generate random file name to store search result in json format
	nameFIle=names.get_first_name()
	jsonfilename=nameFIle+'_basic_search_covid19.json'
	jsonfilepath=os.path.join(settings.BASE_DIR, 'resultFile', 'jsonData','resultJson', 'basicsearch', 'results', jsonfilename)
	jsonfileoutput= open(jsonfilepath,'w')
	jfinaldata=[]
	es.indices.refresh(index="mrmassaydb-index")
	#elasticsearch will search data
	res=helpers.scan(client=es,scroll='2m',index="mrmassaydb-index", doc_type="mrmassaydb-type",query=query,request_timeout=30)
	jfinaldata=[]
	pepSeqList=[]
	proteinList=[]
	#if data is valid based on uniprotkb release then it will display
	for i in res:
		jdic=i['_source']
		jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
		if jdic["UniprotKb entry status"] =="Yes" and jdic['UniProtKB Accession'] !='502':
			jdic["PPI"] ="View"
			jdic["sel"] =""
			jdic["Drug Bank"]=jdic["Drug Bank"].replace('\\','')
			jdic["Drug Bank"]=jdic["Drug Bank"].replace('<br>','|')
			jdic["SRMAtlas URL"]=jdic["SRMAtlas URL"].replace('\\','')
			jdic["Passel URL"]=jdic["Passel URL"].replace('\\','')
			jdic["CPTAC URL"]=jdic["CPTAC URL"].replace('\\','')
			jdic["Panoramaweb URL"]=jdic["Panoramaweb URL"].replace('\\','')
			jdic["PeptideTracker URL"]=jdic["PeptideTracker URL"].replace('\\','')
			#if jdic["Pathway Name"].lower() !='na':
			#	jdic["Pathway Name"]=re.sub(r"(\w)([A-Z])",r"\1|\2",jdic["Pathway Name"])
			jdic["Mean Concentration"] =jdic["Mean Concentration"].replace('fmol/','fmol/µ')
			jdic["Concentration"] =jdic["Concentration"].replace('fmol/','fmol/µ')
			if str(jdic["Associated with COVID-19"]).strip().upper() =='YES':
				pepSeqList.append(jdic["Peptide Sequence"].strip())
				proteinList.append(jdic["UniProtKB Accession"].strip().split('-')[0])
			jfinaldata.append(jdic)
	es.indices.refresh(index="mrmassaydb-index")
	#checking any result generated by database
	foundHits=len(jfinaldata)
	#storing only 10000 rows in json format
	json.dump(jfinaldata[:10000],jsonfileoutput)
	jsonfileoutput.close()
	# if result found then do other job
	if foundHits >0:
		statsummary=summaryStatcal(jfinaldata) # sent data to this funcation for generating stat
		pathwaychart=statsummary['pathwaychart']
		pathwaychart=[i[:2] for i in pathwaychart]
		specieslist=statsummary['specieslist']
		totallist=statsummary['total']
		subcell=statsummary['subcell']
		godic=statsummary['godic']
		jvennprot=statsummary['jevennstat'][0]
		jvennpep=statsummary['jevennstat'][1]
		mrmdatabase=statsummary['jevennstat'][2]
		sortedgodic=OrderedDict(sorted(godic.items(), key=lambda t: t[1])) # sorting GO data
		updatedgodic=dict(list(sortedgodic.items()))
		pepseqdataseries=ast.literal_eval(json.dumps(statsummary['pepseqdataseries'])) #dumping data into json format
		prodataseries=statsummary['prodataseries']
		unqisostat=statsummary['unqisostat']
		jsonfilepathStat=os.path.join(settings.BASE_DIR, 'resultFile', 'jsonData','resultJson', 'basicsearch', 'statsummary', jsonfilename) #storing stat result in json format
		jsonfileoutputStat= open(jsonfilepathStat,'w')
		json.dump(statsummary,jsonfileoutputStat)
		jsonfileoutputStat.close()
		urlname="'/resultFile/jsonData/resultJson/basicsearch/results/"+jsonfilename+"'"

		contextindex={
			"filename":urlname,"colname":json.dumps(colname),'foundHits':foundHits,
			'pathwaychart':pathwaychart[:11],'specieslist':specieslist,
			'totallist':totallist,'subcell':subcell,
			'updatedgodic':updatedgodic,'pepseqdataseries':pepseqdataseries,
			'prodataseries':prodataseries,'unqisostat':unqisostat,
			'uniquePepSeq':len(set(pepSeqList)),'uniqueProtein':len(set(proteinList)),
			'jvennprot':json.dumps(jvennprot),'jvennpep':json.dumps(jvennpep),'jvennmrmdb':json.dumps(mrmdatabase)
			}
		return render(request,'covid19.html',contextindex)
	else:
		return render(request,'covid19.html',{'foundHits':foundHits})

def drugData(request):
	"""
	This is for FDA Assay Tab.
	"""
	ip = get_ip(request, right_most_proxy=True)
	IpAddressInformation.objects.create(ip_address=ip)

	if 'Uniprotkb' in request.GET and request.GET['Uniprotkb']:
		uniprotkb=request.GET['Uniprotkb']
		uniprotkb=uniprotkb.strip()
		drugPresence=0
		drugData=['NA']
		protname=None
		geneName=None
		es.indices.refresh(index="mrmassaydb-index")
		query={"query": {
			"bool": {
				"must": [
					{"match": {"UniProtKB Accession": uniprotkb}},
					{"match": {"UniprotKb entry status": "Yes"}}
				]
			}
		}
		}

		res=helpers.scan(client=es,scroll='2m',index="mrmassaydb-index", doc_type="mrmassaydb-type",query=query,request_timeout=30)
		for hit in res:
			jdic=hit["_source"]
			if jdic["UniProtKB Accession"] ==uniprotkb.upper():
				jdic={str(tkey):force_text(tvalue) for tkey,tvalue in jdic.items()}
				protname=str(jdic["Protein"]).strip()
				genename=str(jdic["Gene"]).strip()
				jdic["sel"] =""
				jdic["Drug Bank"]=jdic["Drug Bank"].replace('<br>','|')
				jdic["Drug Bank"]=jdic["Drug Bank"].replace('</a>"','</a>')
				jdic["Drug Bank"]=jdic["Drug Bank"].replace('href=\\','href=')
				jdic["Drug Bank"]=jdic["Drug Bank"].replace('\\">','">')
				if len(jdic["Drug Bank"].strip()) > 0 and jdic["Drug Bank"].strip() != "NA":
					#tempLst =resitem["Human Drug Bank"].strip().split('|')
					drugData=jdic["Drug Bank"].strip().split('|')
					drugPresence=len(drugData)
				else:
					drugData=['NA']
				break
		es.indices.refresh(index="mrmassaydb-index")


		return render(request, 'drug.html',{'uniprotkb':uniprotkb,\
			'drugPresence':drugPresence,'drugDataRaw':json.dumps(drugData),\
			'protname':protname,'gene':genename})