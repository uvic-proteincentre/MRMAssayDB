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
def reactomeData():
	reactomeDataDic={}
	reactomeFileName='UniProt2Reactome.txt'
	reactomeURL="https://reactome.org/download/current/UniProt2Reactome.txt"
	filepath = os.getcwd()

	print("Extracting Reactome data, job starts",str(datetime.datetime.now()))
	try:
		urllib.urlretrieve(reactomeURL,filepath+'/UniProt2Reactome.txt')
		urllib.urlcleanup()
		print("Extracting Reactome data, job done",str(datetime.datetime.now()))
	except:
		print ("Can't able to download UniProt2Reactome.txt file!!")

	reactomedf= pd.read_csv(reactomeFileName, delimiter='\t',header=None)

	reactomeUniList=list(reactomedf[0].unique())
	for rUni in reactomeUniList:
		tempReactomePathID=list(reactomedf[1][reactomedf[0]==rUni].unique())[0]
		tempReactomePathNameID='|'.join(list(reactomedf[[3,1]][reactomedf[0]==rUni].agg(':'.join, axis=1)))
		tempReactomePathNames='|'.join(list(reactomedf[3][reactomedf[0]==rUni].unique()))
		reactomeDataDic[rUni]=[tempReactomePathNames,tempReactomePathNameID]
		print(reactomeDataDic);sys.exit();

reactomeData()