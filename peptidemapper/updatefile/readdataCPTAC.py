import os,subprocess,psutil,re,shutil,sys,csv,datetime,time
from os import listdir

import urllib,urllib2
import unicodedata
import xmltodict
from xml.dom import minidom
from xml.parsers.expat import ExpatError
import labkey
import html2text
from bs4 import BeautifulSoup
from Bio import SeqIO
from urlparse import unquote
from xml.etree import cElementTree as ET
import itertools
from selenium import webdriver
from selenium.webdriver.support.ui import Select
import multiprocessing
import ctypes
import pandas as pd

class ExtractCPTAC():
	def __init__(self):
		self.curr_dir=os.getcwd()
		self.cptacpath=self.curr_dir+'/cptac/'
		self.output_download_cptac_file= 'cptac_data.csv'
		self.cptac_file_withPanoramaLink= 'cptac_data_withPanoramaLink'
		self.RETRY_TIME = 20.0
		self.finalcptacFile='final_report_cptac_data'
		self.numberOfprocess=4
		self.cptacPanoramaLinkData=[]
		self.cptacDataArray=[]

	def downloadCPTACdata(self):
		os.chdir(self.cptacpath)
		chrome_options = webdriver.ChromeOptions()
		chrome_options.add_argument("--allow-running-insecure-content")
		# chrome_options.add_argument('--headless')
		preferences = {"download.default_directory": self.cptacpath ,"directory_upgrade": True,"safebrowsing.enabled": True }
		chrome_options.add_experimental_option("prefs", preferences)
		driverpath=self.cptacpath +'driver/chrome/linux/chromedriver'
		driver = webdriver.Chrome(chrome_options=chrome_options,executable_path=driverpath)
		driver.get("https://assays.cancer.gov/available_assays")
		dropdowndata = Select(driver.find_element_by_name("browse_table_length"))
		#click analyze button
		#dropdowndata.select_by_visible_text("All") # no longer needed
		time.sleep(20)
		# #download file

		driver.find_element_by_xpath("//*[@id='table_csv']").click()
		time.sleep(20)
		driver.quit()
		filenames = listdir(self.cptacpath)
		
		try:
			cptacoutputfilename=[ filename for filename in filenames if filename.endswith( ".csv" ) and filename.startswith( "CPTAC_Assays_export" ) ][0]

			if os.path.isfile(cptacoutputfilename) and len(cptacoutputfilename.strip()) >0:
				outputresultcptac=open(self.output_download_cptac_file,'w')
				filetoread = open(cptacoutputfilename,'r')
				cptacdic={}
				with open(cptacoutputfilename) as csvfile:
					reader = csv.DictReader(csvfile)
					for row in reader:
						uid=(((row['Gene']).split('UniProt Accession ID:'))[1]).strip()
						pepseq=(row['Proteins and peptides with assays']).strip()
						cptacid=(row['CPTAC ID']).strip()
						cptacinfo=[pepseq,cptacid]
						if len(uid) >0 and cptacid>0:
							if cptacdic.has_key(uid):
								cptacdic[uid].append(cptacinfo)
							else:
								cptacdic[uid]=[cptacinfo]

				cptacheader=('UniprotID'+'\t'+'Peptide Sequence'+'\t'+'Cptac ID')
				outputresultcptac.write(cptacheader+'\n')
				for key in cptacdic.keys():
					uniq_set = set(map(tuple,cptacdic[key]))
					uniq_list = map(list,uniq_set)
					for citem in uniq_list:
						cptacdata='\t'.join(citem)
						outputresultcptac.write(key+'\t'+cptacdata+'\n')

				outputresultcptac.close()
			os.remove(cptacoutputfilename)
		except IndexError:
			pass  
		os.chdir(self.curr_dir)

	def getPanoramaLink(self,inputs):
		cptacdatalist=[]
		header=['UniprotID','Peptide Sequence','Peptide Modified Sequence','Cptac ID','PanoramaLink','Instrument']
		cptacdatalist.append(header)
		for info in inputs[1]:
			cptacstripseq=info[1]
			cptacID=info[-1]
			cptacURL="https://assays.cancer.gov/"+cptacID
			panromproturl="NA"
			modpepseq="NA"
			cpinstrument="NA"

			while True:
				try:
					cptacurlfile=urllib.urlopen(cptacURL)
					for cpline in cptacurlfile:
						cpurldata=cpline.strip()
						
						if '<div style="display:none;" class="sequence_table_shadow detail-modal all-details-'.lower() in cpurldata.lower():
							tempcptacstripseq=cpurldata.split('<div style="display:none;" class="sequence_table_shadow detail-modal all-details-')[1].split('"><')[0].strip()
							if tempcptacstripseq==cptacstripseq:
								tempmodpepseq=cpurldata.split('target="_blank">')[1].split('</a>')[0].strip()
								if '[' in tempmodpepseq:
									modpepseq=tempmodpepseq
						if '<dt>Instrument</dt>'.lower() in cpurldata.lower():
							cpinstrument=(next(cptacurlfile, '').strip()).split("<dd>")[-1].split("<")[0].strip()
						if 'showProtein.view?id='.lower() in cpurldata.lower():
							panromproturl=((cpurldata.split('href='))[-1]).split("class=")[0].strip()
							break
					panromproturl=((cpurldata.split('href='))[-1]).split("class=")[0].strip()
					panromproturl = panromproturl.replace('"', '')
					info.insert(2,modpepseq)
					info.append(panromproturl)
					info.append(cpinstrument)

					cptacurlfile.close()
					break
				except IOError:
					time.sleep(self.RETRY_TIME)
					print('Hey, I am trying again until succeeds to get Panorama Link from CPTAC!',str(datetime.datetime.now()))
					pass
			
			cptacdatalist.append(info)
		with open(self.cptacpath+self.cptac_file_withPanoramaLink+'_Part'+str(inputs[0])+'.csv','w') as cptacfilewithPanoramaLink:
			for citem in cptacdatalist:
				cptacfilewithPanoramaLink.write('\t'.join(citem)+'\n')

	def createCptacArray(self):
		if os.path.isfile(self.cptacpath+self.output_download_cptac_file):
			with open(self.cptacpath+self.output_download_cptac_file,'r') as f:
				reader= csv.reader(f,delimiter='\t')
				self.cptacDataArray=(list(reader))[1:]

	def multiProcssJobGetPanoramaLink(self):
		self.downloadCPTACdata()
		self.createCptacArray()
		if len(self.cptacDataArray)>0:
			print('Number of cptac jobs:',len(self.cptacDataArray))

			bins=int(len(self.cptacDataArray)/self.numberOfprocess)
			counter=0
			fileNames=[]
			processes=[]
			start=time.time()
			for i in range(0,len(self.cptacDataArray),bins):
				counter+=1
				if counter==self.numberOfprocess:
					print(i)
					fileNames.append(self.cptac_file_withPanoramaLink+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.getPanoramaLink,args=[[counter,self.cptacDataArray[i:]]])
					p.start()
					processes.append(p)
					break
				elif counter < self.numberOfprocess:
					fileNames.append(self.cptac_file_withPanoramaLink+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.getPanoramaLink,args=[[counter,self.cptacDataArray[i:i+bins]]])
					p.start()
					processes.append(p)
					print(i,i+bins)
			for process in processes:
				process.join()
			finish=time.time()
			print('Finished in '+str(finish-start)+' seconds for getting Panorama Link')
			#increase the field size of CSV
			csv.field_size_limit(int(ctypes.c_ulong(-1).value // 2))
			df = pd.concat((pd.read_csv(self.cptacpath+f, header = 0, sep='\t',keep_default_na=False) for f in fileNames))
			df_deduplicated = df.drop_duplicates()
			df_deduplicated.to_csv(self.cptacpath+self.cptac_file_withPanoramaLink+".csv",sep='\t', encoding='utf-8',index=False)
			for f in fileNames:
				if os.path.isfile(self.cptacpath+f):
					os.remove(self.cptacpath+f)

	def createCptacPanoramaArray(self):
		if os.path.isfile(self.cptacpath+self.cptac_file_withPanoramaLink+".csv"):
			with open(self.cptacpath+self.cptac_file_withPanoramaLink+".csv",'r') as f:
				reader= csv.reader(f,delimiter='\t')
				self.cptacPanoramaLinkData=(list(reader))[1:]

	def panoramamsrunlist(self,inputs):
		tempfinalreportfile=open(self.cptacpath+self.finalcptacFile+'_Part'+str(inputs[0])+'.csv','w')
		tempfinalreportfile.write('UniprotID'+'\t'+'Peptide Sequence'+'\t'+'Peptide Modified Sequence'+'\t'+'Cptac ID'+'\t'+'CPTAC Transitions'+'\n')
		for data in inputs[1]:
			cpinstrument=data[-1]
			panromproturl=data[-2]
			modpepseq='NA'
			if len(data[2].strip()) !=0:
				modpepseq=data[2].strip()
			cptacstripseq=data[1].strip()
			info=data[:4]
			transList=[]
			while True:
				try:
					proturldata=urllib.urlopen(panromproturl)
					for puline in proturldata:
						pudata=puline.strip()
						if "showPrecursorList.view?id=".lower() in pudata.lower():
							presurid=((pudata.split('</a></li><li><a'))[-1]).split('">')[0].split('href="')[-1].strip()
							msRundecodeURLName=presurid.split('/showPrecursorList.view?')[0].split('/targetedms/')[-1].strip()
							custommsRundecodeURLName=msRundecodeURLName.split('ResponseCurve')[0]+'ResponseCurve'
							msRunserver_context = labkey.utils.create_server_context('panoramaweb.org', custommsRundecodeURLName, use_ssl=True)
							while True:
								try:
									targetmsmy_results = labkey.query.select_rows(
										server_context=msRunserver_context,
										schema_name='targetedms',
										query_name='Transition'
									)

									for targetmsitem in targetmsmy_results.values()[0]:
										tempdic={}
										targetmsID=""
										targetmsmodpepseq=""
										targetmsurl=""
										targetmslabel=""
										targetmsQ1Charge=""
										targetmsQ3Charge=""
										targetmsQ1Mz=""
										targetmsQ3Mz=""
										targetmsFrag=""
										for targetmskey in targetmsitem.keys():
											if 'PeptideId/PeptideGroupId/Label'.lower() in targetmskey.lower() and '_labkeyurl_PrecursorId/PeptideId/PeptideGroupId/Label'.lower() not in targetmskey.lower():
												targetmsID=str(targetmsitem[targetmskey]).strip()
											
											if 'PeptideId/ModifiedPeptideDisplayColumn'.lower() in targetmskey.lower(): targetmsmodpepseq=str(targetmsitem[targetmskey]).strip()
											
											if '_labkeyurl_PrecursorId/PeptideId/PeptideGroupId/Label'.lower() in targetmskey.lower(): targetmsurl=str(targetmsitem[targetmskey]).strip()
											
											if 'PrecursorId/IsotopeLabelId/Name'.lower() in targetmskey.lower(): targetmslabel=str(targetmsitem[targetmskey]).strip()
											
											if 'PrecursorId/Charge'.lower() in targetmskey.lower(): targetmsQ1Charge=str(targetmsitem[targetmskey]).strip()
											if 'Charge'.lower() in targetmskey.lower() and 'PrecursorId/Charge'.lower() not in targetmskey.lower():
												targetmsQ3Charge=str(targetmsitem[targetmskey]).strip()

											if 'PrecursorId/Mz'.lower() in targetmskey.lower(): targetmsQ1Mz=str(targetmsitem[targetmskey]).strip()
											if 'Mz'.lower() in targetmskey.lower() and 'PrecursorId/Mz'.lower() not in targetmskey.lower():
												targetmsQ3Mz=str(targetmsitem[targetmskey]).strip()

											if 'Fragment'.lower() in targetmskey.lower(): targetmsFrag=str(targetmsitem[targetmskey]).strip()


										targetmsmodpepseq=str(targetmsmodpepseq).strip()
										tempdic["instrument"]=cpinstrument
										tempdic["targetmsID"]=targetmsID
										tempdic["targetmsmodpepseq"]=targetmsmodpepseq
										tempdic["targetmslabel"]=targetmslabel
										tempdic["targetmsFrag"]=targetmsFrag
										tempdic["targetmsQ1Charge"]=targetmsQ1Charge
										tempdic["targetmsQ1Mz"]=targetmsQ1Mz
										tempdic["targetmsQ3Charge"]=targetmsQ3Charge
										tempdic["targetmsQ3Mz"]=targetmsQ3Mz
										cppanprotid=panromproturl.split('showProtein.view?id=')[-1].strip()
										panprotid=targetmsurl.split('showProtein.view?id=')[-1].strip()
										if modpepseq !="NA":
											if cppanprotid ==panprotid and len(targetmslabel.strip()) >0 and targetmsmodpepseq.upper()==modpepseq.upper():
												transList.append(tempdic["instrument"]+"|"+tempdic["targetmslabel"]+"|"+tempdic["targetmsQ1Mz"]+"|"+tempdic["targetmsQ1Charge"]+tempdic["targetmsFrag"]+"|"+tempdic["targetmsQ3Mz"]+"|"+tempdic["targetmsQ3Charge"])
										else:
											if cppanprotid ==panprotid and len(targetmslabel.strip()) >0 and targetmsmodpepseq.upper()==cptacstripseq.upper():
												transList.append(tempdic["instrument"]+"|"+tempdic["targetmslabel"]+"|"+tempdic["targetmsQ1Mz"]+"|"+tempdic["targetmsQ1Charge"]+tempdic["targetmsFrag"]+"|"+tempdic["targetmsQ3Mz"]+"|"+tempdic["targetmsQ3Charge"])
									break
								except labkey.exceptions.ServerContextError:
									time.sleep(self.RETRY_TIME)
									print('Hey, I am trying again until succeeds to get data from Panorama by fixing labkey!',str(datetime.datetime.now()))
									pass
					proturldata.close()
					break
				except IOError:
					time.sleep(self.RETRY_TIME)
					print('Hey, I am trying again until succeeds to get data from panorama!',str(datetime.datetime.now()))
					pass
			if len(transList) >0:
				transData=",".join(list(set(transList)))
				transData="Instrument|Label|Q1 m/z|Q1 Z|Fragment|Q3 m/z|Q3 Z,"+transData
				info.append(transData)
			else:
				info.append("NA")
			tempfinalreportfile.write('\t'.join(info)+'\n')

		tempfinalreportfile.close()
	def multiProcssJobCPTAC(self):
		self.multiProcssJobGetPanoramaLink()
		self.createCptacPanoramaArray()
		if len(self.cptacPanoramaLinkData)>0:

			bins=int(len(self.cptacPanoramaLinkData)/self.numberOfprocess)
			counter=0
			fileNames=[]
			processes=[]
			start=time.time()
			for i in range(0,len(self.cptacPanoramaLinkData),bins):
				counter+=1
				if counter==self.numberOfprocess:
					print(i)
					fileNames.append(self.finalcptacFile+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.panoramamsrunlist,args=[[counter,self.cptacPanoramaLinkData[i:]]])
					p.start()
					processes.append(p)
					break
				elif counter < self.numberOfprocess:
					fileNames.append(self.finalcptacFile+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.panoramamsrunlist,args=[[counter,self.cptacPanoramaLinkData[i:i+bins]]])
					p.start()
					processes.append(p)
					print(i,i+bins)
			for process in processes:
				process.join()
			finish=time.time()
			print('Finished in '+str(finish-start)+' seconds for getting Panorama data')
			#increase the field size of CSV
			csv.field_size_limit(int(ctypes.c_ulong(-1).value // 2))
			df = pd.concat((pd.read_csv(self.cptacpath+f, header = 0, sep='\t',keep_default_na=False) for f in fileNames))
			df_deduplicated = df.drop_duplicates()
			df_deduplicated.to_csv(self.cptacpath+self.finalcptacFile+".csv",sep='\t', encoding='utf-8',index=False)
			for f in fileNames:
				if os.path.isfile(self.cptacpath+f):
					os.remove(self.cptacpath+f)
			movefilepath=os.path.join(self.cptacpath+self.finalcptacFile+".csv")
			filepath = os.path.join(self.curr_dir,self.finalcptacFile+".csv")
			shutil.move(movefilepath,filepath)
			print('CPTAC Final Report Done')