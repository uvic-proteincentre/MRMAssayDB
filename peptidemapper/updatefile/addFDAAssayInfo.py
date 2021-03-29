#!/usr/bin/env python
# -*- coding: utf-8 -*-
import urllib,urllib2
from bioservices.kegg import KEGG
import os,subprocess,psutil,re,shutil,datetime,sys,glob
from operator import itemgetter
import numpy as np
import random, time
from itertools import count, groupby
import pandas as pd
import csv

homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)

try:
	peptrackurl='http://peptidetracker.proteincentre.com/rest/api/?searchtype=Any&searchterm=PEP'
	u=urllib.URLopener()
	f=u.open(peptrackurl)
	output_peptrack_file= 'allPepTrackResult.csv'
	p=open(output_peptrack_file,'w')
	p.write(f.read())
	p.close()
except IOError:
	pass


peptrackerConc={}
fdaApproveFile='FDA_Approved_Assays.csv'
with open('allPepTrackResult.csv','r') as peptrackfile:
	peptrackcsvreader=csv.DictReader(peptrackfile)
	for prow in peptrackcsvreader:
		acc=prow['UniProtKB Accession'].strip().upper()
		pepid=prow['Peptide ID'].strip().upper()
		conc=prow['Concentration Range'].strip()
		if len(conc.strip()) >0 and conc.strip().upper() !='NA':
			coninfo=conc.strip().split(';')
			if len(coninfo) >0:
				subconinfo=coninfo[0].split('|')
				concdata='Mean Conc.:'+subconinfo[0].strip()+'<br/> Matrix:'+subconinfo[4].strip()
				peptrackerConc[acc]=[conc,pepid,concdata]

filename='ReportBook_mother_file.csv'
filepath = os.path.join(homedir, 'src/mappermotherfile', filename)

fdaApproveHumanAcc=[]
fdaApproveFile='FDA_Approved_Assays.csv'
with open(fdaApproveFile,'r') as fdafile:
	fdacsvreader=csv.DictReader(fdafile)
	for frow in fdacsvreader:
		humanAcc=frow['Accession(Human)'].strip().split(',')
		for i in humanAcc:
			if len(i.strip()) >0:
				fdaApproveHumanAcc.append(i.strip().upper())

fdaApproveHumanAcc=list(set(fdaApproveHumanAcc))


headers =[]
updatedResult=[]
with open(filepath,'r') as repfile:
	repcsvreader=csv.DictReader(repfile, delimiter='\t')
	headers = repcsvreader.fieldnames
	headers=headers+['Assays for FDA approved Marker','Mean Concentration','Concentration']
	updatedResult.append(headers)
	for reprow in repcsvreader:
		acc=reprow['UniProtKB Accession'].strip().upper()
		peptrackID=reprow['PeptideTracker ID'].strip().upper()
		if acc in fdaApproveHumanAcc:
			reprow['Assays for FDA approved Marker']='Yes'
			if acc in peptrackerConc and peptrackID in peptrackerConc[acc]:
				peptrackerConc[acc][0]=str(peptrackerConc[acc][0]).replace('µ','')
				reprow['Concentration']=peptrackerConc[acc][0]
				peptrackerConc[acc][2]=str(peptrackerConc[acc][2]).replace('µ','')
				reprow['Mean Concentration']=peptrackerConc[acc][2]
			else:
				reprow['Concentration']='No'
				reprow['Mean Concentration']='No'
		else:
			reprow['Assays for FDA approved Marker']='No'
			reprow['Concentration']='No'
			reprow['Mean Concentration']='No'
		tempupdatedResult=[]
		for c in headers:
			tempupdatedResult.append(reprow[c])
		updatedResult.append(tempupdatedResult)

modrepfile="modreprot.csv"
with open(modrepfile,'wb') as modpep:
	modpepwriter=csv.writer(modpep, delimiter='\t')
	modpepwriter.writerows(updatedResult)
	
movefilepath=os.path.join(homedir, 'updatefile', filename)
os.rename(modrepfile,filename)
shutil.move(movefilepath,filepath)