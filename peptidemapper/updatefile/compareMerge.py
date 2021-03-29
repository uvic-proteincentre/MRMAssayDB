#!/usr/bin/env python
# -*- coding: utf-8 -*-
import unicodedata
import datetime,glob
import time
import os,subprocess,psutil,re,sys,shutil
import csv

def compareMerge(reportFile,tempReportFile):
	csv.field_size_limit(sys.maxsize)
	cols=['UniProtKB Accession', 'Protein', 'Gene', 'Organism', 'Organism ID', \
	'SubCellular', 'Peptide Sequence', 'Modified Peptide Sequence', 'Unique in protein',\
	'Present in isoforms', 'PeptideTracker ID', 'PeptideTracker Transition',\
	'Passel ID', 'Passel Transition', 'SRMAtlas ID', 'SRMAtlas Transition', \
	'Cptac ID', 'CPTAC Transitions', 'Panoramaweb ID', 'Panoramaweb Transition', \
	'Kegg Pathway Name', 'Disease Name', 'Go ID', 'Go Name', 'Go Term', 'Drug Bank',\
	'UniprotKb entry status', 'PeptideTracker URL', 'PeptideTracker TransView', \
	'Passel URL', 'Passel TransView', 'SRMAtlas URL', 'SRMAtlas TransView', \
	'CPTAC URL', 'CPTAC TransView', 'Panoramaweb URL', 'Panoramaweb TransView', \
	'Peptide Occurrence', 'Best Transition', 'Summary Transition']
	homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
	movefilepathTempReport=os.path.join(homedir, 'updatefile', tempReportFile)
	filepathReport = os.path.join(homedir, 'src/mappermotherfile', reportFile)
	filepathTempReport = os.path.join(homedir, 'src/mappermotherfile', tempReportFile)

	tempUnIDPepList=[]
	tempReportData=[]
	tempReportData.append(cols)
	with open(filepathTempReport,'r') as tempRepfile:
		tempRepcsvreader = csv.DictReader(tempRepfile, delimiter='\t')
		for trRow in tempRepcsvreader:
			tempList=[]
			if len(trRow['Modified Peptide Sequence'].strip())==0:
				trRow['Modified Peptide Sequence']='NA'
			tempUnIDPepList.append(trRow['UniProtKB Accession'].strip()+trRow['Peptide Sequence'].strip()+trRow['Modified Peptide Sequence'].strip())
			for c in cols:
				tempList.append(trRow[c])
			tempReportData.append(tempList)

	print("Temp Done",str(datetime.datetime.now()))
	repUnIDPepList=[]
	with open(filepathReport,'r') as repfile:
		repcsvreader = csv.DictReader(repfile, delimiter='\t')
		for rRow in repcsvreader:
			if len(rRow['Modified Peptide Sequence'].strip())==0:
				rRow['Modified Peptide Sequence']='NA'
			repUnIDPepList.append(rRow['UniProtKB Accession'].strip()+rRow['Peptide Sequence'].strip()+rRow['Modified Peptide Sequence'].strip())
	
	print("Report Done and Length of old number of Unique UniID and new Unique UniID",len(set(repUnIDPepList)),len(set(tempUnIDPepList)),str(datetime.datetime.now()))
	diffRepUniIDPepList=list(set(repUnIDPepList) - set(tempUnIDPepList))
	print("Length of Difference",len(diffRepUniIDPepList),str(datetime.datetime.now()))
	getindex=[repUnIDPepList.index(i) for i in diffRepUniIDPepList]
	repUnIDPepList=[]
	tempUnIDPepList=[]
	diffRepUniIDPepList=[]
	print("Difference Done",str(datetime.datetime.now()))

	with open(filepathReport,'r') as arepfile:
		arepcsvreader = csv.DictReader(arepfile, delimiter='\t')
		for num, line in enumerate(arepcsvreader):
			if num in getindex:
				tempList=[]
				for c in cols:
					tempList.append((dict(line))[c])
				tempReportData.append(tempList)

	print("Difference Added",str(datetime.datetime.now()))

	with open(tempReportFile,'wb') as rf:
		rwriter =csv.writer(rf,delimiter='\t')
		rwriter.writerows(tempReportData)
	shutil.move(movefilepathTempReport,filepathTempReport)