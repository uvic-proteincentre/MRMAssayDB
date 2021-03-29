#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os,subprocess,psutil,re,sys,shutil,csv,time,datetime,glob
import urllib,requests
import html2text
import unicodedata
from bs4 import BeautifulSoup
from Bio import SeqIO
import xmltodict
from xml.dom import minidom
from xml.parsers.expat import ExpatError
from collections import OrderedDict
import labkey
from urlparse import unquote
from xml.etree import cElementTree as ET
import itertools
import multiprocessing
import pandas as pd
import ctypes

class ExtractPanoramaData():
	def __init__(self):
		self.curr_dir=os.getcwd()
		self.tempdir='tempPanorama'
		if not os.path.exists(self.tempdir):
			os.makedirs(self.tempdir)
		self.finalreportfilename='final_report_panorama_data'
		self.panoramaurl='https://panoramaweb.org/project/Panorama%20Public/begin.view?Targeted%20MS%20Experiment%20List.Organism~isnonblank=&Targeted%20MS%20Experiment%20List.sort=Organism&Targeted MS Experiment List.offset='
		self.targetexpdetailsinfo=[]
		self.listOfTargetProjects=[0,100,200,300,400,500]
		self.RETRY_TIME = 20.0
		self.numberOfprocess=4
		self.processes=[]
		self.tempCSVfilepath=os.path.join(self.curr_dir,self.tempdir)

	def homepagedata(self):
		countProject=0
		for projRange in self.listOfTargetProjects:
			print('Project Range:',projRange)
			while True:
				try:
					panoramaurlhome=urllib.urlopen(self.panoramaurl+str(projRange)).read()
					phsoup = BeautifulSoup(panoramaurlhome,"lxml")
					phsoup.prettify()
					phtables = phsoup.find_all('table')
					for phtable in phtables:
						try:
							phtable_id = phtable['id']
							if phtable_id.startswith('lk-region-'):
								phtr_tags = phtable.find_all('tr')
								for phtr in phtr_tags:
									try:
										phtr_name = phtr['class']
										if phtr_name[0].startswith('labkey-alternate-row') or phtr_name[0].startswith('labkey-row'):
											phtd_tags = phtr.find_all('td')
											countTD=0
											tdinsntrument='NA'
											tempphurl='NA'
											for phtd in phtd_tags:
												if countTD==4:
													tdinsntrument=phtd.string
												phhrf=phtd.find_all('a', href=True)
												phtext = phtd.string
												try:
													if len(phhrf[0]) >0:
														if str(phhrf[0]).startswith('<a href="/Panorama%20Public/'):
															tempphurl="https://panoramaweb.org"+str(phhrf[0]['href'])
												except IndexError:
													pass
												countTD+=1
											if 'NA' !=tempphurl:
												countProject+=1
												self.targetexpdetailsinfo.append([countProject,tempphurl,tdinsntrument])
									except KeyError:
										pass
						except KeyError:
							pass
					break
				except IOError:
					time.sleep(self.RETRY_TIME)
					print('Hey, I am trying again until succeeds to get data from Panorama!',str(datetime.datetime.now()))
					pass

	def panoramamsrunlist(self,inputs):
		if os.path.exists(self.tempdir):
			os.chdir(self.tempdir)
			tempfinalreportfile=open(self.finalreportfilename+'_'+str(inputs[0])+'.csv','w')
			tempfinalreportfile.write('UniprotID'+'\t'+'PeptideSequence'+'\t'+'PeptideModifiedSequence'+'\t'+'URLaddress'+'\t'+'Transitions'+'\n')
			for hinfo in inputs[1]:
				countProt=0
				listofinstruments=[]
				for i in hinfo[2].split(','):
					if 'and' not in i:
						listofinstruments.append(i.encode('utf-8').strip())
					else:
						for x in i.split('and'):
							listofinstruments.append(x.encode('utf-8').strip())
				reload(sys)
				sys.setdefaultencoding("utf-8")
				urladdressexp=hinfo[1]
				msRundecodeURL=unquote(unquote(urladdressexp))
				print(hinfo[0],urladdressexp,str(datetime.datetime.now()))
				#msRundecodeURLName=(msRundecodeURL.split('project/')[1]).split('/begin.view?')[0]
				msRundecodeURLName=(msRundecodeURL.split('/project-begin.view?')[0]).split('panoramaweb.org/')[1]
				msRunserver_context = labkey.utils.create_server_context('panoramaweb.org', msRundecodeURLName, use_ssl=True)
				while True:
					try:
						targetmsmy_results = labkey.query.select_rows(
							server_context=msRunserver_context,
							schema_name='targetedms',
							query_name='Transition'
						)

						targetmsurldic={}
						targetpepmsdic={}
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


							tempdic["targetmsID"]=targetmsID
							tempdic["targetmsmodpepseq"]=targetmsmodpepseq
							tempdic["targetmslabel"]=targetmslabel
							tempdic["targetmsFrag"]=targetmsFrag
							tempdic["targetmsQ1Charge"]=targetmsQ1Charge
							tempdic["targetmsQ1Mz"]=targetmsQ1Mz
							tempdic["targetmsQ3Charge"]=targetmsQ3Charge
							tempdic["targetmsQ3Mz"]=targetmsQ3Mz
							tempmspepID=targetmsID+"@"+targetmsmodpepseq
							if len(tempmspepID.strip()) >0 and len(targetmslabel.strip()) >0:
								if targetpepmsdic.has_key(tempmspepID):
									if tempdic not in targetpepmsdic[tempmspepID]:
										targetpepmsdic[tempmspepID].append(tempdic)
								else:
									targetpepmsdic[tempmspepID]=[tempdic]
							targetmsurldic[targetmsID]=targetmsurl
						for protkey in targetmsurldic.keys():
							countProt+=1
							if 'showProtein.view?id' in targetmsurldic[protkey]:
								uniprotidexist=False
								while True:
									try:
										if 'uniprot.org/uniprot/' in urllib.urlopen("https://panoramaweb.org"+str(targetmsurldic[protkey])).read():
											uniprotidexist=True
										break
									except IOError:
										time.sleep(self.RETRY_TIME)
										print('Hey, I am trying again until succeeds to get data from Panorama!',str(datetime.datetime.now()))
										pass
								if uniprotidexist:
									uniprotlinklist=[]
									while True:
										try:
											proturldata=urllib.urlopen("https://panoramaweb.org"+str(targetmsurldic[protkey]))
											for puline in proturldata:
												pudata=puline.strip()
												if "href=" in pudata and \
												"https://www.uniprot.org/uniprot/" in pudata and \
												"target=" in pudata and \
												"protWindow" in pudata and \
												'<a' in pudata:
													subunid=(((pudata.split('href='))[-1]).split('>'))[0].split()[0]
													uniprotlinklist.append(subunid)
												if "href=" in pudata and \
												"http://www.uniprot.org/uniprot/" in pudata and \
												"target=" in pudata and \
												"protWindow" in pudata and \
												'<a' in pudata:
													subunid=(((pudata.split('href='))[-1]).split('>'))[0].split()[0]
													uniprotlinklist.append(subunid)
											proturldata.close()
											break
										except IOError:
											time.sleep(self.RETRY_TIME)
											print('Hey, I am trying again until succeeds to get data from Panorama!',str(datetime.datetime.now()))
											pass
									if len(uniprotlinklist)>0:
										if countProt %1000 ==0:
											print(str(countProt), "th protein job done",str(datetime.datetime.now()))
										for uniprotlink in uniprotlinklist:
											uniprotlink=uniprotlink.replace('http','https')
											uniprotlink = uniprotlink.replace("'","")
											uniprotlink = uniprotlink.replace('"','')
											try:
												r = requests.get(uniprotlink)
												AC='NA'
												AC=str(r.url).split('/')[-1].strip()
												if AC != 'NA':
													for pepkey in targetpepmsdic.keys():
														if protkey in pepkey.split('@')[0]:
															modpepseq='NA'
															pepseq=filter(str.isalpha, str(pepkey.split('@')[1]))
															protURL="https://panoramaweb.org"+str(targetmsurldic[protkey])
															transList=[]
															for peplabitem in targetpepmsdic[pepkey]:
																for insitem in listofinstruments:
																	transList.append(insitem+"|"+peplabitem["targetmslabel"]+"|"+peplabitem["targetmsQ1Mz"]+"|"+peplabitem["targetmsQ1Charge"]+peplabitem["targetmsFrag"]+"|"+peplabitem["targetmsQ3Mz"]+"|"+peplabitem["targetmsQ3Charge"])
															transData=",".join(list(set(transList)))
															transData="Instrument|Label|Q1 m/z|Q1 Z|Fragment|Q3 m/z|Q3 Z,"+transData
															if pepseq !=str(pepkey.split('@')[1]):
																modpepseq=str(pepkey.split('@')[1]).strip()
															rowdata=AC+'\t'+str(pepseq)+'\t'+str(modpepseq)+'\t'+str(protURL)+'\t'+"'"+transData+"'"
															tempfinalreportfile.write(rowdata+'\n')
											except requests.exceptions.ConnectionError:
												pass
											except requests.exceptions.ChunkedEncodingError:
												pass
						break
					except labkey.exceptions.ServerContextError:
						time.sleep(self.RETRY_TIME)
						print('Hey, I am trying again until succeeds to get data from Panorama by fixing labkey!',str(datetime.datetime.now()))
						pass
			tempfinalreportfile.close()
			os.chdir(self.curr_dir)

	def multiProcssJobPanoramaWeb(self):
		if os.path.exists(self.tempdir):
			self.homepagedata()
			print('Number of projects:',len(self.targetexpdetailsinfo))

			bins=int(len(self.targetexpdetailsinfo)/self.numberOfprocess)
			counter=0
			start=time.time()
			for i in range(0,len(self.targetexpdetailsinfo),bins):
				counter+=1
				if counter==self.numberOfprocess:
					print(i)
					p=multiprocessing.Process(target=self.panoramamsrunlist,args=[[counter,self.targetexpdetailsinfo[i:]]])
					p.start()
					self.processes.append(p)
					break
				elif counter < self.numberOfprocess:
					p=multiprocessing.Process(target=self.panoramamsrunlist,args=[[counter,self.targetexpdetailsinfo[i:i+bins]]])
					p.start()
					self.processes.append(p)
					print(i,i+bins)
			for process in self.processes:
				process.join()
			finish=time.time()
			print('Finished in '+str(finish-start)+' seconds')

			panoramaCSV_files = glob.glob(self.tempCSVfilepath+"/*.csv") 
			#increase the field size of CSV
			csv.field_size_limit(int(ctypes.c_ulong(-1).value // 2))
			df = pd.concat((pd.read_csv(f, header = 0, sep='\t',keep_default_na=False) for f in panoramaCSV_files))
			df_deduplicated = df.drop_duplicates()
			df_deduplicated.to_csv(self.finalreportfilename+".csv",sep='\t', encoding='utf-8',index=False)

			movefilepath=os.path.join(self.curr_dir,self.tempdir ,self.finalreportfilename+".csv")
			filepath = os.path.join(self.curr_dir, self.finalreportfilename+".csv")
			os.chdir(self.curr_dir)
			shutil.move(movefilepath,filepath)
			shutil.rmtree(self.tempdir)