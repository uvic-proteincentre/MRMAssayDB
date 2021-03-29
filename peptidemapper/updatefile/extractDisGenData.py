#!/usr/bin/env python
# -*- coding: utf-8 -*-
#!/usr/bin/env.python
# -*- coding: utf-8 -*-
# encoding: utf-8
import os,subprocess,psutil,re,shutil,datetime,sys,glob
import urllib,urllib2,urllib3
from socket import error as SocketError
import errno
import random, time
import requests
from sh import gunzip
import pandas as pd
import pickle
import cPickle
import chardet
import unicodedata
def disGenData():
	disGenDataDic={}
	disgenFileName='all_gene_disease_associations.tsv'
	disGenURL="https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz"
	filepath = os.getcwd()

	print("Extracting DisGen data, job starts",str(datetime.datetime.now()))
	try:
		urllib.urlretrieve(disGenURL,filepath+'/all_gene_disease_associations.tsv.gz')
		urllib.urlcleanup()
		print("Extracting DisGen data, job done",str(datetime.datetime.now()))
	except:
		print ("Can't able to download all_gene_disease_associations.tsv.gz file!!")

	if os.path.exists(filepath+'/all_gene_disease_associations.tsv'):
		os.remove(filepath+'/all_gene_disease_associations.tsv')
	print("Extracting .gz data, job starts",str(datetime.datetime.now()))
	gunzip(filepath+'/all_gene_disease_associations.tsv.gz')
	print("Extracting .gz data, job done",str(datetime.datetime.now()))
	disgendf= pd.read_csv(disgenFileName, delimiter='\t')
	disGenList=list(disgendf['geneSymbol'].unique())
	for gene in disGenList:
		#tempDisGenID=list(disgendf['geneId'][disgendf['geneSymbol']==gene].unique())[0]
		tempDF=disgendf[['diseaseName','diseaseId']][(disgendf['geneSymbol']==gene) & (disgendf['diseaseType'] =='disease')]
		tempDisNamesInfo=list(zip(tempDF['diseaseName'],tempDF['diseaseId']))
		tempDisNames=[i[0] for i in tempDisNamesInfo]
		tempDisNamesURL=['<a target="_blank" href="https://www.disgenet.org/search/0/'+i[1]+'">'+i[0]+'</a>' for i in tempDisNamesInfo]
		tempDisNames=map(str,tempDisNames)
		if len(tempDisNamesInfo)>0:
			disGenDataDic[gene]=[tempDisNames,tempDisNamesURL]
	dicfile='disGen.obj'
	dicf = open(dicfile, 'wb')
	pickle.dump(disGenDataDic, dicf , pickle.HIGHEST_PROTOCOL)
	dicf.close()
	return dicfile
