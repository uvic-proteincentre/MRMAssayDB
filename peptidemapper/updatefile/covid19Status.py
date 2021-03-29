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
from Bio import SeqIO

class Covid19Status():
	def updateJob(self):
		filename='ReportBook_mother_file.csv'
		homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
		filepath = os.path.join(homedir, 'src/mappermotherfile', filename)
		covidUniProtAcc=self.extractCovidData()
		df= pd.read_csv(filepath, delimiter='\t',keep_default_na=False)
		df.fillna('NA',inplace=True)
		df = df.replace(r'^\s*$', 'NA', regex=True)
		df['Associated with COVID-19']=df['UniProtKB Accession'].isin(covidUniProtAcc)
		df['Associated with COVID-19'] = df['Associated with COVID-19'].replace({True: 'Yes', False: 'No'})
		df=df.apply(lambda row: self.modifyUniDisName(row), axis=1)
		df.to_csv(filepath, index=False, sep="\t")


	def extractCovidData(self):
		covidFastaFileName='covid-19.fasta'
		filepath = os.getcwd()
		covidFastaURL="ftp://ftp.uniprot.org/pub/databases/uniprot/pre_release/covid-19.fasta"
		try:
			urllib.urlretrieve(covidFastaURL,filepath+'/'+covidFastaFileName)
			urllib.urlcleanup()
		except:
			print ("Can't able to download covid-19.fasta file!!")

		print "Extracting Covid-19 disease data, job done",str(datetime.datetime.now())
		covidUniProtAccList=[]
		for seqrecord in SeqIO.parse(filepath+'/'+covidFastaFileName, 'fasta'):
			seqheader = seqrecord.id
			covidUniProtAccList.append(seqheader.split('|')[1].strip())
		return list(set(covidUniProtAccList))

	def modifyUniDisName(self,row):
		if len(row["UniProt DiseaseData"].strip()) > 0 and row["UniProt DiseaseData"].strip() != "NA":
			row["UniProt DiseaseData URL"]=row["UniProt DiseaseData URL"].replace('</a>"','</a>')
			row["UniProt DiseaseData URL"]=row["UniProt DiseaseData URL"].replace('href=\\','href=')
			row["UniProt DiseaseData URL"]=row["UniProt DiseaseData URL"].replace('\\">','">')
			diseaseUniProtURL=row["UniProt DiseaseData URL"].strip().split('|')
			diseaseUniProt=row["UniProt DiseaseData"].strip().split('|')
			if str(row["Associated with COVID-19"]).strip().upper() =='YES':
				covidURL='<a target="_blank" href="https://covid-19.uniprot.org/uniprotkb/'+str(row["UniProtKB Accession"]).strip()+'">'+'COVID-19'+'</a>'
				diseaseUniProtURL.append(covidURL)
				diseaseUniProt.append('COVID-19')
			row["UniProt DiseaseData URL"]='|'.join(diseaseUniProtURL)
			row["UniProt DiseaseData"]='|'.join(diseaseUniProt)
		else:
			if str(row["Associated with COVID-19"]).strip().upper() =='YES':
				covidURL='<a target="_blank" href="https://covid-19.uniprot.org/uniprotkb/'+str(row["UniProtKB Accession"]).strip()+'">'+'COVID-19'+'</a>'
				diseaseUniProtURL=[covidURL]
				diseaseUniProt=['COVID-19']
				row["UniProt DiseaseData URL"]='|'.join(diseaseUniProtURL)
				row["UniProt DiseaseData"]='|'.join(diseaseUniProt)
			else:
				diseaseUniProtURL=['NA']
				diseaseUniProt=['NA']
				row["UniProt DiseaseData URL"]='|'.join(diseaseUniProtURL)
				row["UniProt DiseaseData"]='|'.join(diseaseUniProt)
		if str(row["Associated with COVID-19"]).strip().upper() =='YES':
			if len(row["Disease Name"].strip()) >0 and row["Disease Name"].strip() !='NA':
				tempDiseaseData=row["Disease Name"].split('|')
				tempDiseaseData.append('COVID-19')
				row["Disease Name"]='|'.join(tempDiseaseData)
			else:
				row["Disease Name"]='COVID-19'
		return row