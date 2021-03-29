import unicodedata,datetime,time,glob
import os,subprocess,psutil,re,sys,shutil
import urllib,urllib2,urllib3,httplib
import html2text
from xml.etree import cElementTree as ET
import xmltodict
from xml.dom import minidom
from xml.parsers.expat import ExpatError
import itertools
import pandas as pd
import csv, sqlite3
import multiprocessing
import ctypes

class PasselIDMappingStat():
	def __init__(self):
		self.curr_dir=os.getcwd()
		self.passelpath=self.curr_dir+'/passel/'
		self.passeldic={}
		self.output_download_passel_file_id= 'unqPasselData_f10'
		self.RETRY_TIME = 20.0
		self.listExpId=[]
		self.uniProtIDMapFile='idmapping.dat.2015_03_col_1_3_uniref_uniq'
		self.passelProteinID=set()
		self.numberOfprocess=6
		self.sucessMappingID=set()
		self.failedMappingID=set()


	def getPasselProteinMapID(self,passProteinInfo):
		if '<I><FONT COLOR' not in passProteinInfo:
			passProteinID=''
			if len(passProteinInfo)>0:
				if 'DECOY_' in passProteinInfo:
					passProteinID=passProteinInfo.split('_')[1]
				elif 'gi|' == passProteinInfo[:3]:
					passProteinID=passProteinInfo.split('|')[1]
				elif 'NX_' == passProteinInfo[:3]:
					passProteinID=passProteinInfo.split('-')[0]
				elif len(passProteinInfo) > 6 and '-' == passProteinInfo[6]:
					passProteinID=passProteinInfo.split('-')[0]
				elif 'sp|' == passProteinInfo[:3] or 'tr|' == passProteinInfo[:3]:
					tempData=passProteinInfo.split('|')[1]
					if '_' in tempData:
						passProteinID=tempData.split('_')[0]
					else:
						passProteinID=tempData
				else:
					passProteinID=passProteinInfo
				return passProteinID

	def createPasselObject(self):
		if os.path.isfile(self.passelpath+self.output_download_passel_file_id):
			with open(self.passelpath+self.output_download_passel_file_id,'r') as filetoreadpassel:
				for pdline in filetoreadpassel:
					pddata=pdline.strip()
					try:
						tuid=self.getPasselProteinMapID(pddata)
						if tuid != None:
							self.passelProteinID.add(tuid)
					except IndexError:
						pass

	def findMatchID(self,pid):
		presult=''
		try:
			command="grep -F -m 1  '"+pid+"' "+self.uniProtIDMapFile
			out=subprocess.check_output(command,shell=True)
			out=out.strip()
			presult=out.split()[0]
		except subprocess.CalledProcessError:
			pass
		return presult
	def mapId(self,inputs):
		idmppaingfilepath = self.curr_dir
		print('ID mapping starts for Processer:',str(inputs[0]))
		countID=0
		for pkey in inputs[1]:
			countID+=1
			# if countID ==1:
			# 	print(str(countID), "th id map check job starts for processor",str(inputs[0]), ' at',str(datetime.datetime.now()))
			# if countID%100 ==0:
			# 	print(str(countID), "th id map job starts for processor",str(inputs[0]), ' at',str(datetime.datetime.now()))
			# if countID == len(inputs[1]):
			# 	print(str(countID), "th id map job starts for processor",str(inputs[0]), ' at',str(datetime.datetime.now()))

			accessid=self.findMatchID(pkey)
			if len(accessid.strip())>0:
				self.sucessMappingID.add(pkey)
			else:
				self.failedMappingID.add(pkey)
		print('ID mapping done for Processer:',str(inputs[0]),len(self.sucessMappingID),len(self.failedMappingID))
	def multiProcessMapData(self):
		print('Passel Data Extraction Done')
		self.createPasselObject()
		if len(self.passelProteinID) >0:
			processes=[]
			passelPrimaryProtID=list(self.passelProteinID)
			passelPrimaryProtID=sorted(passelPrimaryProtID)
			print('Number of Proteins:',len(passelPrimaryProtID))
			bins=int(len(passelPrimaryProtID)/self.numberOfprocess)
			counter=0
			start=time.time()
			for i in range(0,len(passelPrimaryProtID),bins):
				counter+=1
				if counter==self.numberOfprocess:
					print(i)
					p=multiprocessing.Process(target=self.mapId,args=[[counter,passelPrimaryProtID[i:]]])
					p.start()
					processes.append(p)
					break
				elif counter < self.numberOfprocess:
					p=multiprocessing.Process(target=self.mapId,args=[[counter,passelPrimaryProtID[i:i+bins]]])
					p.start()
					processes.append(p)
					print(i,i+bins)
			for process in processes:
				process.join()
			finish=time.time()
			print('Finished in '+str(finish-start)+' seconds for get Passel Mapping ID job')
			#increase the field size of CSV
			print('No of sucessfully Mapped ID:',len(self.sucessMappingID))
			print('No of failed Mapped ID:',len(self.failedMappingID))
			print('Passel Final Report Done')

passelIDMapping=PasselIDMappingStat()
passelIDMapping.multiProcessMapData()