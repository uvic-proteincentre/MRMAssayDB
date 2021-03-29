import os,subprocess,psutil,re,shutil,datetime,sys,glob
import urllib,urllib2,urllib3,httplib
from bioservices.kegg import KEGG
from socket import error as SocketError
import errno
from Bio import SeqIO
import xmltodict
from xml.dom import minidom
from xml.parsers.expat import ExpatError
import random, time
from goatools import obo_parser
import csv
from elasticsearch import Elasticsearch,helpers
import json
import pandas as pd
import requests
from collections import Counter
from itertools import combinations, chain
from xml.etree import cElementTree as ET
import xmltodict
from xml.dom import minidom
from xml.parsers.expat import ExpatError
import pickle
import cPickle
import operator
import multiprocessing
import ctypes
from sh import gunzip
from downloadUniProtIDMapping import downloadUniprotIDMapping
from readPassel import ExtractPassel
from readdataCPTAC import ExtractCPTAC
from panoramaDataExtraction import ExtractPanoramaData
from extractPeptideTracker import ExtractPeptideTrackerData

from preloadData import preLoadDataJson
from addModCol import addModCol
from statBestPepTrans import statBestPepTrans
from uploadDataElasticSearch import uploadData
from updateDiseaseData import updateDisData
from covid19Status import Covid19Status
from updatePathwayData import UpdatePathName
from generateSummaryReport import SummaryStat
from maketotalassaypep import TotalAssay
from addReviewed import AddReviewStatus

class GenerateMRMAssayDBFile():
	def __init__(self):
		self.curr_dir=os.getcwd()
		self.RETRY_TIME = 20.0
		self.numberOfprocess=6
		self.listOfResourceFile=['final_report_peptrack_data.csv',\
		'final_report_passel_data.csv','final_report_cptac_data.csv',\
		'final_report_panorama_data.csv','final_report_srm_data.csv']
		self.mergedictdata={}
		self.runcomplete=0
		#get home directory path
		self.homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
		self.outfilefilename ='outfilefile.csv'
		self.filename='ReportBook_mother_file.csv'
		self.filenameTempReport='Temp_ReportBook_mother_file.csv'
		self.keggfile='keggdata.csv'
		self.filepath = os.path.join(self.homedir, 'src/mappermotherfile', self.filename)
		self.filepathTempReport = os.path.join(self.homedir, 'src/mappermotherfile', self.filenameTempReport)
		self.protgnogscgofilename='uniprotfuncdata'
		self.protPresencefilename='uniprotpresence'
		self.below_curr_dir=os.path.normpath(self.curr_dir + os.sep + os.pardir)
		self.canonisounidic={}
		self.unqcanonicalUnId=[]
		self.unquniidlist=[]
		self.unqisocheckdic={}
		self.functionalData={}
		self.unifuncdic={}
		self.diseasefilepath = os.path.join(self.below_curr_dir, 'src/UniDiseaseInfo/' 'humsavar.txt')
		self.outfilefileUnqIsoname='UnqIsoresult'
		self.uniprotFastaFileName="uniprot-reviewed_yes.fasta"
		self.unikeggunqisofuncdic={}
		self.transdic={}
		self.movefilepath=os.path.join(self.homedir, 'updatefile', self.filename)
		self.movefilepathTempReport=os.path.join(self.homedir, 'updatefile', self.filenameTempReport)
		self.fastaDataDic={}

	def runResourceProg(self):
		# downloadUniprotIDMapping()
		# srmcmd='python readSRMAtlast.py'
		# subprocess.Popen(srmcmd, shell=True).wait()
		print("SRMAtlas data extraction job done")
		#cptac=ExtractCPTAC()
		#cptac.multiProcssJobCPTAC()
		print("CPTAC data extraction job done")
		#passel=ExtractPassel()
		#passel.multiProcessMapData()
		print("Passel data extraction job done")
		#panorama=ExtractPanoramaData()
		#panorama.multiProcssJobPanoramaWeb()
		print("Panorama data extraction job done")
		#peptrack=ExtractPeptideTrackerData()
		#peptrack.fetchDataFromPeptidetracker()
		print("PeptideTracker data extraction job done")
		self.runcomplete=1

	def copyMotherFile(self):
		# copy mother file from its source to working directory
		if os.path.exists(self.filepath):
			if os.path.exists(self.movefilepath):
				os.remove(self.movefilepath)
			if not os.path.exists(self.movefilepath):
				shutil.copy2(self.filepath, self.movefilepath)
				#create backup folder before update and then move that folder with old version mother file for backup
				mydate = datetime.datetime.now()
				folder_name="version_"+mydate.strftime("%B_%d_%Y_%H_%M_%S")
				backupPath='/home/bioinf/mnt/media/mrmassaydb-vm/backup/'
				if not os.path.exists(backupPath+folder_name):
					os.makedirs(backupPath+folder_name)
				if not os.path.exists(backupPath+folder_name+'/ReportBook_mother_file.csv'):
					shutil.copy2(self.movefilepath, backupPath+folder_name+'/ReportBook_mother_file.csv')

	def makeMergeDict(self):
		for fitem in self.listOfResourceFile:
			with open(fitem,'r') as pepfile:
				for line in pepfile:
					data=line.strip()
					if not data.startswith('UniprotID'):
						info=data.split('\t')
						if '?' not in info[0].strip() or '=' not in info[0].strip() or ':' not in info[0].strip():
							if info[0].strip() in self.mergedictdata:
								self.mergedictdata[info[0].strip()].append(info[1])
							else:
								self.mergedictdata[info[0].strip()]=[info[1]]
			print(str(fitem),"data dictionay job done",str(datetime.datetime.now()))
		self.mergedictdata={a:list(set(b)) for a, b in self.mergedictdata.items()}

	def checkTotalUniqProtein(self):
		totalUniId=list(set(self.mergedictdata.keys()))
		for tuid in totalUniId:
			tuempcode=(str(tuid).split('-'))[0]
			if '?' not in tuempcode or '=' not in tuempcode or ':' not in tuempcode:
				if tuempcode in self.canonisounidic:
					self.canonisounidic[tuempcode].append(tuid)
				else:
					self.canonisounidic[tuempcode]=[tuid]

		self.canonisounidic={a:list(set(b)) for a, b in self.canonisounidic.items()}
		self.unqcanonicalUnId=list(set(self.canonisounidic.keys()))
		self.unqcanonicalUnId=sorted(self.unqcanonicalUnId)
		print("Total Unique protein in this file: ",len(self.unqcanonicalUnId))

	def extractProteionFunctionalInformation(self,inputs):
		countProt=0
		protgnogscgofile=open(self.protgnogscgofilename+'_Part'+str(inputs[0])+'.csv','w')
		protgnogscgofile.write('ActualUniID'+'\t'+'PepSeq'+'\t'+'ProteinName'+'\t'+'Gene'+'\t'+'Organism'+'\t'+'OrganismID'+'\t'+'Subcellular'+'\t'+'GOID'+'\t'+'GOName'+'\t'+'GoTerm'+'\t'+'DrugBank'+'\t'+'DiseaseData'+'\n')
		for subcgcode in inputs[1]:
			time.sleep(2)
			ScAllLocList=[]
			GoIDList=[]
			GoNamList=[]
			GoTermList=[]
			drugbanklist=[]
			PN='NA'
			GN='NA'
			OG='NA'
			OGID='NA'
			dislist=[]
			while True:
				try:
					countProt+=1
					if countProt ==1:
						print(str(countProt), "th protein Protein Name, Gene, Organism Name, GO, sub cellular,drug bank data,disease data job starts for processor",str(inputs[0]), ' at',str(datetime.datetime.now()))
					if countProt%1000 ==0:
						print(str(countProt), "th protein Protein Name, Gene, Organism Name, GO, sub cellular,drug bank data,disease data job starts for processor",str(inputs[0]), ' at',str(datetime.datetime.now()))
					if countProt == len(inputs[1]):
						print(str(countProt), "th protein Protein Name, Gene, Organism Name, GO, sub cellular,drug bank data,disease data job starts for processor",str(inputs[0]), ' at',str(datetime.datetime.now()))

					SGrequestURL="https://www.uniprot.org/uniprot/"+str(subcgcode)+".xml"
					SGunifile=urllib.urlopen(SGrequestURL)
					SGunidata= SGunifile.read()
					SGunifile.close()

					try:
						SGunidata=minidom.parseString(SGunidata)
						try:
							subcelldata=(SGunidata.getElementsByTagName('subcellularLocation'))
							for subcItem in subcelldata:
								try:
									subloc=(subcItem.getElementsByTagName('location')[0]).firstChild.nodeValue
									if len(str(subloc).strip()) >0:
										ScAllLocList.append(str(subloc).strip())
								except:
									pass

						except IndexError:
							pass
						try:
							godata=(SGunidata.getElementsByTagName('dbReference'))
							for gItem in godata:
								if (gItem.attributes['type'].value).upper() == 'GO':
									try:
										gonamedetails=(str(gItem.getElementsByTagName('property')[0].attributes['value'].value).strip()).split(':')[1]
										gotermdetails=(str(gItem.getElementsByTagName('property')[0].attributes['value'].value).strip()).split(':')[0]
										GoNamList.append(gonamedetails)
										goid=str(gItem.attributes['id'].value).strip()
										GoIDList.append(goid)
										if gotermdetails.lower()=='p':
											GoTermList.append('Biological Process')
										if gotermdetails.lower()=='f':
											GoTermList.append('Molecular Function')
										if gotermdetails.lower()=='c':
											GoTermList.append('Cellular Component')
									except:
										pass
								if (gItem.attributes['type'].value).upper() == 'DRUGBANK':
									try:
										drugname=(str(gItem.getElementsByTagName('property')[0].attributes['value'].value).strip())
										drugid=str(gItem.attributes['id'].value).strip()
										durl='<a target="_blank" href="https://www.drugbank.ca/drugs/'+drugid+'">'+drugname+'</a>'
										drugbanklist.append(durl)
									except:
										pass
								if (gItem.attributes['type'].value).strip() == 'NCBI Taxonomy':
									try:
										OGID=str(gItem.attributes['id'].value).strip()
										OGID=str(int(float(OGID)))
									except:
										pass
						except IndexError:
							pass

						try:
							try:
								PN=(((SGunidata.getElementsByTagName('protein')[0]).getElementsByTagName('recommendedName')[0]).getElementsByTagName('fullName')[0]).firstChild.nodeValue

							except:
								PN=(((SGunidata.getElementsByTagName('protein')[0]).getElementsByTagName('submittedName')[0]).getElementsByTagName('fullName')[0]).firstChild.nodeValue

						except IndexError:
							pass

						try:
							try:
								GN=((SGunidata.getElementsByTagName('gene')[0]).getElementsByTagName('name')[0]).firstChild.nodeValue
							except:
								GN='NA'
						except IndexError:
							pass

						try:
							try:
								OG=((SGunidata.getElementsByTagName('organism')[0]).getElementsByTagName('name')[0]).firstChild.nodeValue
							except:
								OG='NA'
						except IndexError:
							pass

						try:
							disdata=SGunidata.getElementsByTagName('disease')
							for dItem in disdata:
								disname=''
								disshort=''
								try:
									disname=(dItem.getElementsByTagName('name')[0]).firstChild.nodeValue
								except:
									pass
								try:
									disshort=(dItem.getElementsByTagName('acronym')[0]).firstChild.nodeValue
								except:
									pass
								if len(disname.strip())>0:
									dislist.append(str(disname.strip())+'('+str(disshort)+')')
						except IndexError:
							pass

					except ExpatError:
						pass
					break
				except IOError:
					time.sleep(self.RETRY_TIME)
					print('Hey, I am trying again until succeeds to get funcational data from UniProt!',str(datetime.datetime.now()))
					pass
			subcelldata='NA'
			goiddata='NA'
			gonamedata='NA'
			gotermdata='NA'
			drugbankdata='NA'
			disdata='NA'
			if len(ScAllLocList)>0:
				subcelldata='|'.join(list(set(ScAllLocList)))
			if len(GoIDList)>0:
				goiddata='|'.join(list(set(GoIDList)))
			if len(GoNamList)>0:
				gonamedata='|'.join(list(set(GoNamList)))
			if len(GoTermList)>0:
				gotermdata='|'.join(list(set(GoTermList)))
			if len(drugbanklist)>0:
				drugbankdata='|'.join(list(set(drugbanklist)))
			if len(dislist)>0:
				disdata='|'.join(list(set(dislist)))

			for canisoitem in self.canonisounidic[subcgcode]:
				for temppgopepseq in self.mergedictdata[canisoitem]:
					protgnogscgofile.write(str(canisoitem)+'\t'+str(temppgopepseq)+'\t'+str(PN)+'\t'+str(GN)+'\t'+str(OG)+'\t'+str(OGID)+'\t'+str(subcelldata)+'\t'+str(goiddata)+'\t'+str(gonamedata)+'\t'+str(gotermdata)+'\t'+str(drugbankdata)+'\t'+str(disdata)+'\n')

		protgnogscgofile.close()
		
	def multiProcessGetFuncationalData(self):
		print("Extracting Protein Name, Gene, Organism Name,GO,Sub cellular data, drug bank data ,disease data and checking pep seq present in uniprot specified seq, job starts",str(datetime.datetime.now()))
		if len(self.unqcanonicalUnId)>0:
			fileNames=[]

			bins=int(len(self.unqcanonicalUnId)/self.numberOfprocess)
			counter=0
			processes=[]
			start=time.time()
			for i in range(0,len(self.unqcanonicalUnId),bins):
				counter+=1
				if counter==self.numberOfprocess:
					fileNames.append(self.protgnogscgofilename+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.extractProteionFunctionalInformation,args=[[counter,self.unqcanonicalUnId[i:]]])
					p.start()
					processes.append(p)
					break
				elif counter < self.numberOfprocess:
					fileNames.append(self.protgnogscgofilename+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.extractProteionFunctionalInformation,args=[[counter,self.unqcanonicalUnId[i:i+bins]]])
					p.start()
					processes.append(p)
			for process in processes:
				process.join()
			finish=time.time()
			print('Finished in '+str(finish-start)+' seconds for get Uniprot Funcational data job')
			#increase the field size of CSV
			csv.field_size_limit(int(ctypes.c_ulong(-1).value // 2))
			df = pd.concat((pd.read_csv(f, header = 0, sep='\t',keep_default_na=False) for f in fileNames))
			df_deduplicated = df.drop_duplicates()
			df_deduplicated.to_csv(self.protgnogscgofilename+".csv",sep='\t', encoding='utf-8',index=False)
			for f in fileNames:
				if os.path.isfile(f):
					os.remove(f)
		self.unqcanonicalUnId=[]
		print("Extracting Protein Name, Gene, Organism Name,GO,Sub cellular data, drug bank data, disease data and checking pep seq present in uniprot specified seq, job done",str(datetime.datetime.now()))

	def functionalDataArray(self):
		with open(self.protgnogscgofilename+".csv",'r') as f:
			for line in f:
				data=line.strip()
				if not data.startswith("ActualUniID"):
					info=data.split('\t')
					tempId=info[0].strip()
					tempSeq=info[1].strip()
					tempKey=tempId+'_'+tempSeq
					if '?' not in tempId or '=' not in tempId or ':' not in tempId:
						if tempId in self.functionalData:
							self.functionalData[tempId].append(info)
						else:
							self.functionalData[tempId]=[info]

	def checkPepSeqInUniprot(self,inputs):
		countRow=0
		protPresencefile=open(self.protPresencefilename+'_Part'+str(inputs[0])+'.csv','w')
		protPresencefile.write('ActualUniID'+'\t'+'UpdatedUniID'+'\t'+'PepSeq'+'\t'+'ProteinName'+'\t'+'Gene'+'\t'+'Organism'+'\t'+'OrganismID'+'\t'+'Subcellular'+'\t'+'GOID'+'\t'+'GOName'+'\t'+'GoTerm'+'\t'+'DrugBank'+'\t'+'DiseaseData'+'\t'+'PresentInSeq'+'\n')
		for item in inputs[1]:
			subcgcode=item

			time.sleep(1)
			while True:
				try:
					countRow+=1
					if countRow ==1:
						print(str(countRow), "th protein seq check job starts for processor",str(inputs[0]), ' at',str(datetime.datetime.now()))
					if countRow%1000 ==0:
						print(str(countRow), "th protein seq check job starts for processor",str(inputs[0]), ' at',str(datetime.datetime.now()))
					if countRow == len(inputs[1]):
						print(str(countRow), "th protein seq check job starts for processor",str(inputs[0]), ' at',str(datetime.datetime.now()))

					tempfastaseq=''
					unifastaurl="https://www.uniprot.org/uniprot/"+str(subcgcode)+".fasta"
					fr = requests.get(unifastaurl)
					fAC=(str(fr.url).split('/')[-1].strip()).split('.')[0].strip()
					fastaresponse = urllib.urlopen(unifastaurl)
					for seq in SeqIO.parse(fastaresponse, "fasta"):
						tempfastaseq=(seq.seq).strip()
					for pepSeqItem in self.functionalData[subcgcode]:
						tempPepSeq=pepSeqItem[1]
						tempData=pepSeqItem[2:]
						tempData=['NA' if len(i.strip()) ==0  else i for i in tempData]
						if len(tempfastaseq.strip()) >0:
							pepinfastapresent='No'
							if tempPepSeq in tempfastaseq:
								pepinfastapresent='Yes'
							if '-' in fAC:
								protPresencefile.write(str(subcgcode)+'\t'+str(fAC)+'\t'+str(tempPepSeq)+'\t'+str(tempData[0])+'\t'+str(tempData[1])+'\t'+str(tempData[2])+'\t'+str(tempData[3])+'\t'+str('NA')+'\t'+str('NA')+'\t'+str('NA')+'\t'+str('NA')+'\t'+str('NA')+'\t'+str('NA')+'\t'+str(pepinfastapresent)+'\n')
							else:
								protPresencefile.write(str(subcgcode)+'\t'+str(fAC)+'\t'+str(tempPepSeq)+'\t'+'\t'.join(tempData)+'\t'+str(pepinfastapresent)+'\n')
					
					break
				except IOError:
					time.sleep(self.RETRY_TIME)
					print('Hey, I am trying again until succeeds to get sequence data from UniProt!',str(datetime.datetime.now()))
					pass
		protPresencefile.close()

	def multiProcessPepSeqPresenceData(self):
		print("Data formatting and checking pep seq present in uniprot specified seq, job starts",str(datetime.datetime.now()))
		if len(self.functionalData)>0:
			fileNames=[]
			tempIDs=list(set(self.functionalData.keys()))
			tempIDs=sorted(tempIDs)
			print("Total Data with Canonical and Isoforms:",len(tempIDs))
			bins=int(len(tempIDs)/self.numberOfprocess)
			counter=0
			processes=[]
			start=time.time()
			for i in range(0,len(tempIDs),bins):
				counter+=1
				if counter==self.numberOfprocess:
					print(i)
					fileNames.append(self.protPresencefilename+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.checkPepSeqInUniprot,args=[[counter,tempIDs[i:]]])
					p.start()
					processes.append(p)
					break
				elif counter < self.numberOfprocess:
					fileNames.append(self.protPresencefilename+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.checkPepSeqInUniprot,args=[[counter,tempIDs[i:i+bins]]])
					p.start()
					processes.append(p)
					print(i,i+bins)
			for process in processes:
				process.join()
			finish=time.time()
			print('Finished in '+str(finish-start)+' seconds for get Uniprot Funcational data job')
			#increase the field size of CSV
			csv.field_size_limit(int(ctypes.c_ulong(-1).value // 2))
			df = pd.concat((pd.read_csv(f, header = 0, sep='\t',keep_default_na=False) for f in fileNames))
			df_deduplicated = df.drop_duplicates()
			df_deduplicated.to_csv(self.protPresencefilename+".csv",sep='\t', encoding='utf-8',index=False)
			for f in fileNames:
				if os.path.isfile(f):
					os.remove(f)
		self.mergedictdata.clear()
		print("Data formatting and checking pep seq present in uniprot specified seq, job done",str(datetime.datetime.now()))

	def createUniFuncObject(self):
		uniidlist=[]
		with open(self.protPresencefilename+".csv") as pgosgfile:
			preader = csv.DictReader(pgosgfile, delimiter='\t')
			for prow in preader:
				templist=[str(prow['ActualUniID']).strip(),str(prow['ProteinName']).strip(),str(prow['Gene']).strip(),str(prow['Organism']).strip(),str(prow['OrganismID']).strip(),str(prow['Subcellular']).strip(),str(prow['PepSeq']).strip(),str(prow['DiseaseData']).strip(),str(prow['GOID']).strip(),str(prow['GOName']).strip(),str(prow['GoTerm']).strip(),str(prow['DrugBank']).strip(),str(prow['PresentInSeq']).strip()]
				tempfuncid=str(prow['UpdatedUniID']).strip()+'_'+str(prow['PepSeq']).strip()
				if '?' not in str(prow['UpdatedUniID']).strip() or '=' not in str(prow['UpdatedUniID']).strip() or ':' not in str(prow['UpdatedUniID']).strip():
					self.unifuncdic[tempfuncid]=templist

					uniidlist.append((((prow['UpdatedUniID']).split('-'))[0]).strip())

					if str(prow['PresentInSeq']).strip() =='Yes':
						tempid=str(prow['UpdatedUniID']).strip()+'_'+str(prow['OrganismID']).strip()
						if tempid in self.unqisocheckdic:
							self.unqisocheckdic[tempid].append(str(prow['PepSeq']).strip())
						else:
							self.unqisocheckdic[tempid]=[str(prow['PepSeq']).strip()]

		self.unquniidlist=list(set(uniidlist))
		self.unquniidlist=sorted(self.unquniidlist)

	def getKeggData(self):
		print("Extracting KEGG pathway name, job starts",str(datetime.datetime.now()))
		keggf=open(self.keggfile,'w')
		keggf.write('UniprotID'+'\t'+'PathwayNames'+'\n')
		countProt=0
		uniproturl = 'https://www.uniprot.org/uploadlists/'
		k = KEGG()
		for kx in range(0,len(self.unquniidlist),2000):
			startEachKeggIteration=time.time()
			countProt+=2000
			print(str(countProt), "th protein kegg job starts",str(datetime.datetime.now()))

			uniprotcodes=' '.join(self.unquniidlist[kx:kx+2000])
			uniprotparams = {
			'from':'ACC',
			'to':'KEGG_ID',
			'format':'tab',
			'query':uniprotcodes
			}
			
			while True:
				try:
					kuniprotdata = urllib.urlencode(uniprotparams)
					kuniprotrequest = urllib2.Request(uniproturl, kuniprotdata)
					kuniprotresponse = urllib2.urlopen(kuniprotrequest)
					print(str(countProt), "th protein kegg reading starts",str(datetime.datetime.now()))
					for kuniprotline in kuniprotresponse:
						kudata=kuniprotline.strip()
						if not kudata.startswith("From"):
							kuinfo=kudata.split("\t")
							if len(kuinfo[1].strip())>0:
								kegg=k.get(kuinfo[1].strip())
								kudict_data = k.parse(kegg)
								try:
									try:
										if len(str(kuinfo[0]).strip()) >5:
											tempKeggData='|'.join(kudict_data['PATHWAY'].values())
											keggf.write(str(kuinfo[0]).strip()+'\t'+tempKeggData+'\n')
									except TypeError: 
										pass
								except KeyError:
									pass

					print(str(countProt), "th protein kegg reading ends",str(datetime.datetime.now()))

					break
				except urllib2.URLError:
					time.sleep(self.RETRY_TIME)
					print('Hey, I am trying again until succeeds to get data from KEGG! for processor',str(inputs[0]), ' at',str(datetime.datetime.now()))
					pass
				except urllib2.HTTPError:
					time.sleep(self.RETRY_TIME)
					print('Hey, I am trying again until succeeds to get data from KEGG! for processor',str(inputs[0]), ' at',str(datetime.datetime.now()))
					pass
				except httplib.BadStatusLine:
					time.sleep(self.RETRY_TIME)
					print('Hey, I am trying again until succeeds to get data from KEGG! for processor',str(inputs[0]), ' at',str(datetime.datetime.now()))
					pass
			endEachKeggIteration=time.time()
			print('Finished in '+str(endEachKeggIteration-startEachKeggIteration)+' seconds for each 2000 proteins to get KEGG')
		keggf.close()

		print("Extracting KEGG pathway name, job ends",str(datetime.datetime.now()))

	def extractUniProtDiseaseData(self):
		
		print("Extracting disease data, job starts",str(datetime.datetime.now()))
		try:
			urllib.urlretrieve('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt',self.diseasefilepath)
			urllib.urlcleanup()
		except:
			print("Can't able to download humsavar.txt file!!")

		print("Extracting Human disease data, job done",str(datetime.datetime.now()))

		print("Formatting disease data, job done",str(datetime.datetime.now()))

	def downloadUniprotFastaFile(self):
		print('Download UniProt Fasta File Starts',str(datetime.datetime.now()))
		filepathCanonicalGZ = os.path.join(os.getcwd(), 'uniprot_sprot.fasta.gz')
		filepathCanonical = os.path.join(os.getcwd(), 'uniprot_sprot.fasta')
		try:
			urllib.urlretrieve('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz',filepathCanonicalGZ)
			urllib.urlcleanup()
		except:
			print ("Can't able to download uniprot_sprot.fasta.gz file!")

		filepathIsoGZ = os.path.join(os.getcwd(), 'uniprot_sprot_varsplic.fasta.gz')
		filepathIso = os.path.join(os.getcwd(), 'uniprot_sprot_varsplic.fasta')
		try:
			urllib.urlretrieve('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz',filepathIsoGZ)
			urllib.urlcleanup()
		except:
			print ("Can't able to download uniprot_sprot_varsplic.fasta.gz file!")
		if os.path.exists(filepathCanonical):
			os.remove(filepathCanonical)
		if os.path.exists(filepathIso):
			os.remove(filepathIso)
		gunzip(filepathCanonicalGZ)
		gunzip(filepathIsoGZ)

		canonicalSeq = open(filepathCanonical)
		canonicalSeq_contents = canonicalSeq.read()
		canonicalSeq.close()

		isoSeq = open(filepathIso)
		isoSeq_contents = isoSeq.read()
		isoSeq.close()

		finalFastaFIle = open(self.uniprotFastaFileName, "w") # open in `w` mode to write
		finalFastaFIle.write(canonicalSeq_contents + isoSeq_contents) # concatenate the contents
		print('Download UniProt Fasta File Ends',str(datetime.datetime.now()))

	def createFastaDic(self):
		
		for seq in SeqIO.parse(self.uniprotFastaFileName, "fasta"):
			tempfastaseq=str((seq.seq).strip())
			headerID=seq.id.strip().split('|')[1]
			orgID=seq.description.strip().split('OX=')[1].split(' ')[0]
			if orgID in self.fastaDataDic:
				self.fastaDataDic[orgID].append(headerID+'_'+tempfastaseq)
			else:
				self.fastaDataDic[orgID]=[headerID+'_'+tempfastaseq]

	def protUnqCheck(self,inputs):
		countPepSeq=0
		outfilefileUnqIso = open(self.outfilefileUnqIsoname+'_Part'+str(inputs[0])+'.csv','w')
		outfilefileUnqIso.write('UniProtKB Accession'+'\t'+'Peptide Sequence'+'\t'+'Unique in protein'+'\t'+'Present in isoforms'+'\n')
		for mkey in inputs[1]:
			pepunid=mkey.split('_')[0]
			unqtemppepseqList=list(set(self.unqisocheckdic[mkey]))
			pepUnqDic={}
			pepIsodic={}
			nonprotuniqstatDic={}
			peppresentUniFastaDic={}
			canopepunid=''
			pepunidver=''
			if '-' in pepunid:
				pepunidinfo=pepunid.split('-')
				canopepunid=pepunidinfo[0]
				pepunidver=pepunidinfo[-1]
			else:
				canopepunid=pepunid
			uqorgid=mkey.split('_')[1]
			try:
				tempFastaDicData=self.fastaDataDic[uqorgid]

				countPepSeq+=1
				if countPepSeq ==1:
					print(str(countPepSeq), "th pepseq uniqueness job starts for processor",str(inputs[0]), ' at',str(datetime.datetime.now()))
				if countPepSeq%1000 ==0:
					print(str(countPepSeq), "th pepseq uniqueness job starts for processor",str(inputs[0]), ' at',str(datetime.datetime.now()))
				if countPepSeq == len(inputs[1]):
					print(str(countPepSeq), "th pepseq uniqueness job starts for processor",str(inputs[0]), ' at',str(datetime.datetime.now()))

				
				for matchpeptide in unqtemppepseqList:
					matched=[s for s in tempFastaDicData if matchpeptide in s]

					if len(matched)>0:
						for itemUnq in matched:
							uniID=itemUnq.split('_')[0]

							if pepunid.lower() == (str(uniID).lower()).strip():
								peppresentUniFastaDic[str(matchpeptide).strip()]=True

							canouniID=''
							uniIDver=''
							if '-' in uniID:
								uniIDinfo=uniID.split('-')
								canouniID=uniIDinfo[0]
								uniIDver=uniIDinfo[-1]
							else:
								canouniID=uniID

							if (canouniID.strip()).lower() == (canopepunid.strip()).lower():
								if len(uniIDver.strip()) ==0:
									pepUnqDic[str(matchpeptide).strip()]=True
								if len(uniIDver.strip()) !=0:
									if pepIsodic.has_key(str(matchpeptide).strip()):
										pepIsodic[str(matchpeptide).strip()].append(uniID)
									else:
										pepIsodic[str(matchpeptide).strip()]=[uniID]
							if canouniID.strip() !=canopepunid.strip():
								nonprotuniqstatDic[str(matchpeptide).strip()]=True
		
			except KeyError:
				pass
			for peptideseq in unqtemppepseqList:
				peptideunique='NA'
				pepisodata='No'
				if peptideseq not in nonprotuniqstatDic:
					if peptideseq in pepUnqDic:
						if pepUnqDic[peptideseq]:
							peptideunique='Yes'
						else:
							peptideunique='Not unique'
				else:
					peptideunique='NA'
				if peptideseq in pepIsodic:
					pepisodata=','.join(list(set(pepIsodic[peptideseq])))
				outfilefileUnqIso.write(str(pepunid)+'\t'+str(peptideseq)+'\t'+str(peptideunique)+'\t'+str(pepisodata)+'\n')
				
		outfilefileUnqIso.close()
		
	def multiProcessCheckUniquenes(self):
		print("Checking uniqueness of peptide sequence and presence in isoforms, job starts",str(datetime.datetime.now()))
		if len(self.unqisocheckdic)>0:
			fileNames=[]

			tempIds=self.unqisocheckdic.keys()
			tempIds=sorted(tempIds)
			print("Total Data:",len(tempIds))
			bins=int(len(tempIds)/self.numberOfprocess)
			counter=0
			processes=[]
			start=time.time()
			for i in range(0,len(tempIds),bins):
				counter+=1
				if counter==self.numberOfprocess:
					fileNames.append(self.outfilefileUnqIsoname+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.protUnqCheck,args=[[counter,tempIds[i:]]])
					p.start()
					processes.append(p)
					break
				elif counter < self.numberOfprocess:
					fileNames.append(self.outfilefileUnqIsoname+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.protUnqCheck,args=[[counter,tempIds[i:i+bins]]])
					p.start()
					processes.append(p)
			for process in processes:
				process.join()
			finish=time.time()
			print('Finished in '+str(finish-start)+' seconds for get Uniprot Funcational data job')
			#increase the field size of CSV
			csv.field_size_limit(int(ctypes.c_ulong(-1).value // 2))
			df = pd.concat((pd.read_csv(f, header = 0, sep='\t',keep_default_na=False) for f in fileNames))
			df_deduplicated = df.drop_duplicates()
			df_deduplicated.to_csv(self.outfilefileUnqIsoname+".csv",sep='\t', encoding='utf-8',index=False)
			for f in fileNames:
				if os.path.isfile(f):
					os.remove(f)
		self.fastaDataDic.clear()
		print("Checking uniqueness of peptide sequence and presence in isoforms, job done",str(datetime.datetime.now()))

	def createFunctionalDataDictionay(self):
		print("Functional data dictionay job starts",str(datetime.datetime.now()))
		tempunqisodic={}
		with open(self.outfilefileUnqIsoname+".csv") as unqisofile:
			uireader = csv.DictReader(unqisofile, delimiter='\t')
			for uirow in uireader:
				if '?' not in str(uirow['UniProtKB Accession']).strip() or '=' not in str(uirow['UniProtKB Accession']).strip() or ':' not in str(uirow['UniProtKB Accession']).strip():
					tempunqisodic[str(uirow['UniProtKB Accession']).strip()+'_'+str(uirow['Peptide Sequence']).strip()]=[str(uirow['Unique in protein']).strip(),str(uirow['Present in isoforms']).strip()]
		
		keggDF= pd.read_csv(self.keggfile, sep='\t',keep_default_na=False)
		keggDF_deduplicated = keggDF.drop_duplicates()
		keggdict=pd.Series(keggDF_deduplicated['PathwayNames'].values,index=keggDF_deduplicated['UniprotID']).to_dict()


		for tukey in self.unifuncdic:
			tempkeggdata='NA'
			tempunqdata='NA'
			tempisodata='NA'
			tuniID=tukey.split('_')[0]
			if tuniID in keggdict:
				tempkeggdata='|'.join(list(set(keggdict[tuniID].split('|'))))
			if tukey in tempunqisodic:
				tempunqdata=tempunqisodic[tukey][0]
				tempisodata=tempunqisodic[tukey][1]
			tuitem=self.unifuncdic[tukey]
			tuitem.insert(7,tempunqdata)
			tuitem.insert(8,tempisodata)
			tuitem.insert(9,tempkeggdata)
			tempolduniid=tuitem[0]
			tuitem[0]=tuniID
			modtukey=tempolduniid+'_'+tukey.split('_')[1]
			self.unikeggunqisofuncdic[modtukey]=tuitem
		self.unifuncdic.clear()
		print("Functional data dictionay job done",str(datetime.datetime.now()))

	def createTransitionDict(self):
		for fitem in self.listOfResourceFile:
			print(str(fitem),"transition data dictionay job starts",str(datetime.datetime.now()))
			with open(fitem,'r') as resfile:
				for resline in resfile:
					resdata=resline.strip()
					if not resdata.startswith('UniprotID'):
						resinfo=resdata.split('\t')
						restempid=resinfo[0].strip()+'_'+resinfo[1].strip()+'_'+resinfo[2].strip()
						if '?' not in restempid or '=' not in restempid or ':' not in restempid:
							if restempid in self.transdic:
								self.transdic[restempid].append(resinfo[3:])
							else:
								self.transdic[restempid]=[resinfo[3:]]
			print(str(fitem),"transition data dictionay job done",str(datetime.datetime.now()))

	def createInitialReport(self):
		countEntry =0
		sizeOfTransDic=len(self.transdic)
		time0=time.time()
		print("Initial report file creation, job starts",str(datetime.datetime.now()))
		outfilefile = open(self.outfilefilename,'w')
		outfilefile.write('UniProtKB Accession'+'\t'+'Protein'+'\t'+'Gene'+'\t'+'Organism'+'\t'+'Organism ID'+'\t'+'SubCellular'+'\t'+'Peptide Sequence'+'\t'+'Modified Peptide Sequence'+'\t'+'Unique in protein'+'\t'+'Present in isoforms'+'\t'+'PeptideTracker ID'+'\t'+'PeptideTracker Transition'+'\t'+'Passel ID'+'\t'+'Passel Transition'+'\t'+'SRMAtlas ID'+'\t'+'SRMAtlas Transition'+'\t'+'Cptac ID'+'\t'+'CPTAC Transitions'+'\t'+'Panoramaweb ID'+'\t'+'Panoramaweb Transition'+'\t'+'Kegg Pathway Name'+'\t'+'Disease Name'+'\t'+'Go ID'+'\t'+'Go Name'+'\t'+'Go Term'+'\t'+'Drug Bank'+'\t'+'UniprotKb entry status'+'\n')
		for key in self.transdic:
			countEntry+=1
			if countEntry%1000 ==0:
				print(str(countEntry), "entry out of ",str(sizeOfTransDic),str(datetime.datetime.now()))
				temptime=time.time()-time0
				calTime=((float(temptime)/float(countEntry))*(float(sizeOfTransDic)-float(countEntry)))/(60.0*60.0)
				print('Expected Remaining time:',str(calTime))

			tempsubid='_'.join(key.split('_')[0:2])
			tempmodpepseq=key.split('_')[-1]
			tpeptrackidlist=[]
			tpasselidlist=[]
			tsrmidlist=[]
			tcptacidlist=[]
			tpanoramidlist=[]

			tpeptrackTransList=[]
			tpasselTransList=[]
			tsrmTransList=[]
			tpanoramTransList=[]
			tcptacTransList=[]

			tpeptrackid='NA'
			tpasselid='NA'
			tsrmid='NA'
			tcptacid='NA'
			tpanoramid='NA'
			
			tpeptrackTrans='NA'
			tpasselTrans='NA'
			tsrmTrans='NA'
			tpanoramTrans='NA'
			tcptacTrans='NA'
			for subtempitem in self.transdic[key]:
				subtempitem[-1]=subtempitem[-1].replace("'","")

				if (subtempitem[0][0:3]).lower()=='pep':
					tpeptrackidlist.append(str(subtempitem[0]))
					tpeptrackTransList.append(str(subtempitem[-1]))
				elif (subtempitem[0][0:4]).lower()=='epep':
					 tpeptrackidlist.append(str(subtempitem[0]))
					 tpeptrackTransList.append(str(subtempitem[-1]))

				elif (subtempitem[0][0:4]).lower()=='pass':
					tpasselidlist.append(str(subtempitem[0]))
					tpasselTransList.append(str(subtempitem[-1]))

				elif 'systemsbiology.net' in (subtempitem[0]).lower():
					tsrmidlist.append(str(subtempitem[0]))
					tsrmTransList.append(str(subtempitem[-1]))

				elif (subtempitem[0][0:5]).lower()=='cptac':
					tcptacidlist.append(str(subtempitem[0]))
					tcptacTransList.append(str(subtempitem[-1]))
				elif (subtempitem[0][0:9]).lower()=='non-cptac':
					tcptacidlist.append(str(subtempitem[0]))
					tcptacTransList.append(str(subtempitem[-1]))

				else:
					tpanoramidlist.append(str(subtempitem[0]))
					tpanoramTransList.append(str(subtempitem[-1]))

			if len(tpeptrackidlist)>0:
				tpeptrackid=','.join(list(set(tpeptrackidlist)))
			if len(tpeptrackTransList)>0:
				tpeptrackTrans=','.join(list(set(tpeptrackTransList)))
			if len(tpasselidlist)>0:
				tpasselid=','.join(list(set(tpasselidlist)))
			if len(tpasselTransList)>0:
				tpasselTrans=','.join(list(set(tpasselTransList)))
			if len(tsrmidlist)>0:
				tsrmid=','.join(list(set(tsrmidlist)))
			if len(tsrmTransList)>0:
				tsrmTrans=','.join(list(set(tsrmTransList)))
			if len(tcptacidlist)>0:
				tcptacid=','.join(list(set(tcptacidlist)))
			if len(tcptacTransList)>0:
				tcptacTrans=','.join(list(set(tcptacTransList)))
			if len(tpanoramidlist)>0:
				tpanoramid=','.join(list(set(tpanoramidlist)))
			if len(tpanoramTransList)>0:
				tpanoramTrans=','.join(list(set(tpanoramTransList)))

			temprow=['NA']*27
			temprowpos=[0,1,2,3,4,5,6,8,9,20,21,22,23,24,25,26]
			if tempsubid in self.unikeggunqisofuncdic:
				for x,y in zip(self.unikeggunqisofuncdic[tempsubid],temprowpos):
					temprow[y]=x

			temprow[7]=tempmodpepseq

			temprow[10]=tpeptrackid
			temprow[11]=tpeptrackTrans

			temprow[12]=tpasselid
			temprow[13]=tpasselTrans

			temprow[14]=tsrmid
			temprow[15]=tsrmTrans

			temprow[16]=tcptacid
			temprow[17]=tcptacTrans

			temprow[18]=tpanoramid
			temprow[19]=tpanoramTrans
			if (temprow[0].strip()).upper() != 'NA' and (temprow[6].strip()).upper() != 'NA':
				temprow=['NA' if len(i.strip()) ==0  else i for i in temprow]
				finalreportdata='\t'.join(temprow)
				outfilefile.write(finalreportdata+'\n')
			temprow=[]

		outfilefile.close()
		print("Initial report file creation, job done",str(datetime.datetime.now()))
		self.transdic.clear()
		self.unikeggunqisofuncdic.clear()
		os.rename(self.outfilefilename,self.filenameTempReport)
		shutil.move(self.movefilepathTempReport,self.filepathTempReport)
		print("Initial report file transfer, job done",str(datetime.datetime.now()))


	def finalSteps(self):
		print("Update mother file job starts now",datetime.datetime.now())
		self.runResourceProg()
		if self.runcomplete==1:
			#self.copyMotherFile()
			#self.makeMergeDict()
			#self.checkTotalUniqProtein()
			#self.multiProcessGetFuncationalData()
			#self.functionalDataArray()
			#self.multiProcessPepSeqPresenceData()
			#self.createUniFuncObject()
			#self.getKeggData()
			#self.extractUniProtDiseaseData()
			#self.downloadUniprotFastaFile()
			#self.createFastaDic()
			#self.multiProcessCheckUniquenes()
			#self.createFunctionalDataDictionay()
			#self.createTransitionDict()
			#self.createInitialReport()
			# print ("Add new Columns, job Starts",str(datetime.datetime.now()))
			# addModCol(self.filenameTempReport)
			# print ("Add new Columns, job Done",str(datetime.datetime.now()))
			
			# print ("Best Pepetide Transition, job Starts",str(datetime.datetime.now()))
			# statBestPepTrans(self.filenameTempReport)
			# print ("Best Pepetide Transition, job Done",str(datetime.datetime.now()))
			

			# print ("PreLoad, job Starts",str(datetime.datetime.now()))
			# preLoadDataJson(self.homedir,self.filepathTempReport)
			# print ("PreLoad, job Done",str(datetime.datetime.now()))

			# print ("Rename, job Starts",str(datetime.datetime.now()))
			# os.rename(self.filepathTempReport,self.filepath)
			# print ("Rename, job Done",str(datetime.datetime.now()))

			# print ("Stat Kegg, job Starts",str(datetime.datetime.now()))
			# keggcmd='python statKEGGcoverage.py'
			# subprocess.Popen(keggcmd, shell=True).wait()
			# print ("Stat Kegg, job Done",str(datetime.datetime.now()))
			
			# print ("Generate Summary report, job Starts",str(datetime.datetime.now()))
			# sumStat=SummaryStat()
			# sumStat.finalStat()
			# print ("Generate Summary report, job Done",str(datetime.datetime.now()))

			# print ("Total Assay, job Starts",str(datetime.datetime.now()))
			# totalAssay=TotalAssay()
			# totalAssay.finalStep()
			# print ("Total Assay, job Done",str(datetime.datetime.now()))

			# print ("Update Help, job Starts",str(datetime.datetime.now()))		
			# updatestatjobcmd='python updatestatHelp.py'
			# subprocess.Popen(updatestatjobcmd, shell=True).wait()
			# print ("Update Help, job Done",str(datetime.datetime.now()))

			# print ("Add Sel Col, job Starts",str(datetime.datetime.now()))
			# addseljobcmd='python addSelCol.py'
			# subprocess.Popen(addseljobcmd, shell=True).wait()
			# print ("Add Sel Col, job Done",str(datetime.datetime.now()))
			# print ("Modification of File, job Starts",str(datetime.datetime.now()))
			# modfcolujobcmd='python modFileBasedUserRequirement.py'
			# subprocess.Popen(modfcolujobcmd, shell=True).wait()
			# print ("Modification of File, job Done",str(datetime.datetime.now()))

			# print ("Add FDA Assay, job Starts",str(datetime.datetime.now()))
			# fdaColjobcmd='python addFDAAssayInfo.py'
			# subprocess.Popen(fdaColjobcmd, shell=True).wait()
			# print ("Add FDA Assay, job Done",str(datetime.datetime.now()))

			# print ("Update disease, job Starts",str(datetime.datetime.now()))
			# updateDisData()
			# print ("Update disease, job Done",str(datetime.datetime.now()))

			# print ("Update COVID, job Starts",str(datetime.datetime.now()))
			# covid19=Covid19Status()
			# covid19.updateJob()
			# print ("Update COVID, job Done",str(datetime.datetime.now()))
			#Update pathway will run after yassene gives PC data
			# print ("Update Pathway after adding Pathway common data, job Starts",str(datetime.datetime.now()))
			# updatePath=UpdatePathName()
			# updatePath.updateJob()
			# print ("Update Pathway after adding Pathway common data, job Done",str(datetime.datetime.now()))

			# print ("Update Protein Review Status, job Starts",str(datetime.datetime.now()))
			# addreview=AddReviewStatus()
			# addreview.finalSteps()
			# print ("Update Protein Review Status, job Done",str(datetime.datetime.now()))

			print ("Upload, job Starts",str(datetime.datetime.now()))
			uploadData()
			print ("Upload, job Done",str(datetime.datetime.now()))


		print("Update mother file job done",datetime.datetime.now())


testG=GenerateMRMAssayDBFile()
testG.finalSteps()