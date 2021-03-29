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
class UpdatePathName():
	def updateJob(self):
		filename='ReportBook_mother_file.csv'
		homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
		filepath = os.path.join(homedir, 'src/mappermotherfile', filename)
		df= pd.read_csv(filepath, delimiter='\t',keep_default_na=False)
		#"Kegg Pathway Name"
		pathwaysCols=['Kegg Pathway Name','PC pathway URL','PC pathway name','PC pathway source']
		for pathCol in pathwaysCols:
			df[pathCol].fillna('NA',inplace=True)
		df['Pathway Name']='NA'
		df['Pathway Name']=df.apply(lambda row: self.modifyName(row), axis=1)
		df.fillna('NA',inplace=True)
		df.to_csv(filepath, index=False, sep="\t", na_rep='NA')

	def modifyName(self,row):
		if row['Kegg Pathway Name'].strip().lower() !='na':
			tempKeggData=row['Kegg Pathway Name'].strip().split('|')
			if row['PC pathway name'].strip().lower() !='na':
				tempPCData=row['PC pathway name'].strip().split('|')
				tempPathData='|'.join(list(set(tempKeggData+tempPCData)))
				return tempPathData
			else:
				return row['Kegg Pathway Name']
		else:
			if row['PC pathway name'].strip().lower() !='na':
				if len(row['PC pathway name'].strip())> 0:
					return row['PC pathway name']
				else:
					return 'NA'
			else:
				return 'NA'
