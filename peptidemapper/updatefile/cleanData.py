#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os,subprocess,psutil,re,sys,shutil,datetime,time
import unicodedata
from collections import OrderedDict
import itertools

monthdic={"January":1,"February":2,"March":3,"April":4,"May":5,"June":6,"July":7,"August":8,"September":9,"October":10,"November":11,"December":12}
shortmonthdic={"Jan":1,"Feb":2,"Mar":3,"Apr":4,"May":5,"Jun":6,"Jul":7,"Aug":8,"Sep":9,"Oct":10,"Nov":11,"Dec":12}
currentdata=datetime.datetime.now().strftime("%B_%d_%Y_%H_%M_%S")
currentdataList=currentdata.split('_')
homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
backupDataPath=os.path.join(homedir, 'updatefile/backup/')
listOfbackUpdir=os.listdir(backupDataPath)
listOfresultDir=['src/resultFile/fastaFIle','src/resultFile/jsonData/protvistadataJson','src/resultFile/jsonData/stringJson','src/resultFile/jsonData/resultJson/restapisearch','src/resultFile/jsonData/resultJson/basicsearch/results','src/resultFile/jsonData/resultJson/basicsearch/statsummary','src/resultFile/jsonData/resultJson/adavancesearch/results','src/resultFile/jsonData/resultJson/adavancesearch/statsummary']
listOfFileExtn=['.fasta','.json']
for fname in listOfbackUpdir:
	fnameList=fname.split('_')
	temppath=backupDataPath=os.path.join(homedir, 'updatefile/backup/',fname)
	if abs(int(fnameList[3])-int(currentdataList[2])) ==0:
		if abs(monthdic[fnameList[1]]-monthdic[currentdataList[0]]) >6:
			shutil.rmtree(temppath)
	if abs(int(fnameList[3])-int(currentdataList[2])) !=0:
		if abs(monthdic[fnameList[1]]-monthdic[currentdataList[0]]) ==7:
			shutil.rmtree(temppath)	
	

for rname in listOfresultDir:
	resultDataPath=os.path.join(homedir, rname)
	for rfile in os.listdir(resultDataPath):
		for fexten in listOfFileExtn:
			if rfile.endswith(fexten):
				resultFilePath=os.path.join(homedir,rname,rfile)
				fileCreateDateList=str(time.ctime(os.path.getmtime(resultFilePath))).split(' ')
				if os.path.isfile(resultFilePath):
					if abs(int(fileCreateDateList[-1])-int(currentdataList[2])) ==0:
						if abs(shortmonthdic[fileCreateDateList[1]]-monthdic[currentdataList[0]]) >1:
							os.remove(resultFilePath)
					if abs(int(fileCreateDateList[-1])-int(currentdataList[2])) !=0:
						if abs(shortmonthdic[fileCreateDateList[1]]-monthdic[currentdataList[0]]) ==2:
							os.remove(resultFilePath)
