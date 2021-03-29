#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
import ctypes
from xml.etree import cElementTree as ET
import xmltodict
from xml.dom import minidom
from xml.parsers.expat import ExpatError
from extractDisGenData import disGenData
def updateDisData():
	filename='ReportBook_mother_file.csv'
	homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
	filepath = os.path.join(homedir, 'src/mappermotherfile', filename)
	#disGenDataDicName=disGenData()
	disGenDataDicName='disGen.obj'
	disGenDataDic = cPickle.load(open(disGenDataDicName, 'rb'))
	df= pd.read_csv(filepath, delimiter='\t',keep_default_na=False)
	newDF=df.dropna(subset=['Disease Name'])
	disAssocUniID=list(newDF['UniProtKB Accession'][(newDF['Disease Name'] !='NA')].unique())
	disAssocUniID=list(set(disAssocUniID))
	df['Disease Name']='NA'
	df['UniProt DiseaseData']='NA'
	df['UniProt DiseaseData URL']='NA'
	df['DisGen DiseaseData']='NA'
	df['DisGen DiseaseData URL']='NA'
	countProt=0
	RETRY_TIME=20
	print("total proteins needs to update their unirpot disease URL:",len(disAssocUniID))
	for unID in disAssocUniID:
		GN='NA'
		dislist=[]
		unidislist=[]
		unidisURLlist=[]
		disgendislist=[]
		disgendisURLlist=[]
		subunID=unID.split('-')[0]
		while True:
			try:
				countProt+=1
				if countProt%100 ==0:
					print str(countProt), "th protein Gene, ,disease data job starts",str(datetime.datetime.now())

				SGrequestURL="https://www.uniprot.org/uniprot/"+str(subunID)+".xml"
				SGunifile=urllib.urlopen(SGrequestURL)
				SGunidata= SGunifile.read()
				SGunifile.close()

				try:
					SGunidata=minidom.parseString(SGunidata)

					try:
						try:
							GN=((SGunidata.getElementsByTagName('gene')[0]).getElementsByTagName('name')[0]).firstChild.nodeValue
						except:
							GN='NA'
					except IndexError:
						pass

					try:
						disdata=SGunidata.getElementsByTagName('disease')
						for dItem in disdata:
							disname=''
							disshort=''
							disURL=''
							disID=''
							try:
								disname=(dItem.getElementsByTagName('name')[0]).firstChild.nodeValue
								disID=(dItem.attributes['id'].value).upper()
							except:
								pass
							try:
								disshort=(dItem.getElementsByTagName('acronym')[0]).firstChild.nodeValue
							except:
								pass
							if len(disname.strip())>0:
								disURL='<a target="_blank" href="https://www.uniprot.org/diseases/'+disID+'">'+str(disname.strip())+'('+str(disshort)+')'+'</a>'
								dislist.append(str(disname.strip())+'('+str(disshort)+')')
								unidislist.append(str(disname.strip())+'('+str(disshort)+')')
								unidisURLlist.append(disURL)
					except IndexError:
						pass

				except ExpatError:
					pass
				break
			except IOError:
				time.sleep(RETRY_TIME)
				print('Hey, I am trying again until succeeds to get funcational data from UniProt!',str(datetime.datetime.now()))
				pass
		
		disdata='NA'
		uniDisData='NA'
		uniDisURLData='NA'
		disgenDisData='NA'
		disgenDisURLData='NA'
		if GN != 'NA' and GN in disGenDataDic:
			disgendislist=disGenDataDic[GN][0]
			disgendisURLlist=disGenDataDic[GN][1]
			if len(dislist)>0:
				dislist=dislist+disGenDataDic[GN][0]
			else:
				dislist=disGenDataDic[GN][0]
		if len(dislist)>0:
			disdata='|'.join(list(set(dislist)))
		if len(unidislist)>0:
			uniDisData='|'.join(list(set(unidislist)))
		if len(unidisURLlist)>0:
			uniDisURLData='|'.join(list(set(unidisURLlist)))
		if len(disgendislist)>0:
			disgenDisData='|'.join(list(set(disgendislist)))
		if len(disgendisURLlist)>0:
			disgenDisURLData='|'.join(list(set(disgendisURLlist)))

		df.loc[df['UniProtKB Accession']==unID, 'Disease Name'] = disdata
		df.loc[df['UniProtKB Accession']==unID, 'UniProt DiseaseData'] = uniDisData
		df.loc[df['UniProtKB Accession']==unID, 'UniProt DiseaseData URL'] = uniDisURLData
		df.loc[df['UniProtKB Accession']==unID, 'DisGen DiseaseData'] = disgenDisData
		df.loc[df['UniProtKB Accession']==unID, 'DisGen DiseaseData URL'] = disgenDisURLData
	df.to_csv(filepath, index=False, sep="\t")
