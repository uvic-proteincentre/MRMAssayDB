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

class ExtractPassel():
	def __init__(self):
		self.curr_dir=os.getcwd()
		self.passelpath=self.curr_dir+'/passel/'
		self.passeldic={}
		self.output_download_passel_file= 'passeldownloaddata'
		self.RETRY_TIME = 20.0
		self.listExpId=[]
		self.uniProtIDMapFile='idmapping.dat.2015_03_col_1_3_uniref_uniq'
		self.passelProteinID=set()
		self.finalpasselFile='final_report_passel_data'
		self.numberOfprocess=4

	def getExpSRMId(self):
		passelhomeURL='https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetSELTransitions?_subtab=5'
		passelHomedata=urllib.urlopen(passelhomeURL)
		passelHomedataInfo=passelHomedata.readlines()
		passelHomedataInfo=map(lambda x:x.lower().strip(),passelHomedataInfo)
		passelHomedataInfo=filter(None,passelHomedataInfo)
		srmExpIndex=passelHomedataInfo.index('<td><b>srm experiments:</b></td>')

		for phind,phdata in enumerate(passelHomedataInfo[srmExpIndex:]):
			if '<option value=' in phdata:
				phvalue=phdata.split('<option value="')[1].split('">')[0]
				if len(phvalue)>0:
					self.listExpId.append(int(phvalue))
			if phdata =='<option value=""></option>':
				break
	def getPasselData(self,inputs):
		tempPasselDatafile=open(self.passelpath+self.output_download_passel_file+'_Part'+str(inputs[0])+'.csv','w')
		for passexpID in inputs[1]:
			print('Experiment ID:',passexpID)
			passelQueryUrl="https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetSELTransitions?row_limit=1000000&display_options=ShowEachTransition&SEL_experiments=&QUERY_NAME=AT_GetSELTransitions&action=QUERY&uploaded_file_not_saved=1&apply_action=QUERY"
			newexpinfo="SEL_experiments="+str(passexpID)
			passelQueryUrl=passelQueryUrl.replace("SEL_experiments=",newexpinfo)
			while True:
				try:
					passelurldata=urllib.urlopen(passelQueryUrl)
					for pline in passelurldata:
						pdata=pline.strip()
						if (pdata.lower()).startswith(('<A HREF="/sbeams/cgi/PeptideAtlas/GetSELTransitions/query_guest').lower()) and ('output_mode=tsv').lower() in pdata.lower():
							ptsvlink=(pdata.split('"'))[1]
							downloadlink='https://db.systemsbiology.net'+str(ptsvlink)
							try:
								passelu=urllib.URLopener()
								f=passelu.open(downloadlink)
								for l in f.readlines():
									tempPasselDatafile.write('\t'.join(l.strip().split('\t')[:29])+'\n')
							except IOError:
								pass
							break
					passelurldata.close()
					break
				except IOError:
					time.sleep(self.RETRY_TIME)
					print('Hey, I am trying again until succeeds to get data from Passel!',str(datetime.datetime.now()))
					pass
		tempPasselDatafile.close()


	def multiProcessGetPasselData(self):
		self.getExpSRMId()
		if len(self.listExpId)>0:
			fileNames=[]
			print('Number of projects:',len(self.listExpId))

			bins=int(len(self.listExpId)/self.numberOfprocess)
			counter=0
			processes=[]
			start=time.time()
			for i in range(0,len(self.listExpId),bins):
				counter+=1
				if counter==self.numberOfprocess:
					print(i)
					fileNames.append(self.output_download_passel_file+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.getPasselData,args=[[counter,self.listExpId[i:]]])
					p.start()
					processes.append(p)
					break
				elif counter < self.numberOfprocess:
					fileNames.append(self.output_download_passel_file+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.getPasselData,args=[[counter,self.listExpId[i:i+bins]]])
					p.start()
					processes.append(p)
					print(i,i+bins)
			for process in processes:
				process.join()
			finish=time.time()
			print('Finished in '+str(finish-start)+' seconds for get Passel data job')
			#increase the field size of CSV
			csv.field_size_limit(int(ctypes.c_ulong(-1).value // 2))
			df = pd.concat((pd.read_csv(self.passelpath+f, header = 0, sep='\t',keep_default_na=False) for f in fileNames))
			df_deduplicated = df.drop_duplicates()
			df_deduplicated.to_csv(self.passelpath+self.output_download_passel_file+".csv",sep='\t', encoding='utf-8',index=False)
			for f in fileNames:
				if os.path.isfile(self.passelpath+f):
					os.remove(self.passelpath+f)

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
		if os.path.isfile(self.passelpath+self.output_download_passel_file+".csv"):
			with open(self.passelpath+self.output_download_passel_file+".csv",'r') as filetoreadpassel:
				for pdline in filetoreadpassel:
					pddata=pdline.strip()
					if not (pddata.lower()).startswith(('modified_peptide_sequence').lower()) and not (pddata.lower()).startswith(('<I><FONT').lower()):
						pdinfo=pddata.split('\t')
						try:
							tuid=self.getPasselProteinMapID(pdinfo[9].strip())
							tlist=[pdinfo[0].strip(),pdinfo[25].strip(),pdinfo[11].strip(),pdinfo[24].strip(),pdinfo[12].strip(),pdinfo[16].strip(),pdinfo[15].strip(),pdinfo[1].strip(),pdinfo[14].strip()]
							if tuid != None:
								self.passelProteinID.add(tuid)
								if self.passeldic.has_key(tuid):
									self.passeldic[tuid].append(tlist)
								else:
									self.passeldic[tuid]= [tlist]
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
		outputresultpassel=open(self.passelpath+self.finalpasselFile+'_Part'+str(inputs[0])+'.csv','w')
		passelheader='UniprotID'+'\t'+'PeptideSeq'+'\t'+'Peptide Modified Sequence'+'\t'+'PasselID'+'\t'+'Transitions'
		outputresultpassel.write(passelheader+'\n')
		for pkey in inputs[1]:
			accver=''
			accessid=self.findMatchID(pkey)
			ptemplist=self.passeldic[pkey]
			p_set = set(map(tuple,ptemplist))
			uptemplist = map(list,p_set)
			if len(accessid.strip())>0:
				subpasseldic={}
				subpasselTransdic={}
				for pitem in uptemplist:
					#pepseq=re.sub("[\[].*?[\]]","",str(pitem[0]))
					modpepseq=pitem[0].strip()
					pepseq=filter(str.isalpha, str(pitem[0]))
					passelid=str(pitem[1])
					tempid=pepseq+'-'+passelid+'-'+modpepseq
					if len(str(pitem[3]).strip())>0:
						temptransData='|'.join(pitem[2:])
						templabellist=str(pitem[3]).split('+')
						if subpasseldic.has_key(tempid):
							subpasseldic[tempid].append(templabellist)
						else:
							subpasseldic[tempid]= [templabellist]

						if subpasselTransdic.has_key(tempid):
							subpasselTransdic[tempid].append(temptransData)
						else:
							subpasselTransdic[tempid]= [temptransData]
							
				if len(subpasseldic)>0:
					for spkey in subpasseldic.keys():
						if spkey in subpasselTransdic:
							labelList2d = subpasseldic[spkey]
							mergedlabel = list(itertools.chain.from_iterable(labelList2d))
							mergeddatalabel='+'.join(list(set(mergedlabel)))

							transList = subpasselTransdic[spkey]
							transdataFinal=','.join(transList)
							transdataFinal="Instrument|Label|PrecursorIon|CollisionEnergy|FragmentIon|ParentChargeState|ProductIon,"+transdataFinal
							spkeyinfo=spkey.split('-')
							tempmodpepseq='NA'
							if str(spkeyinfo[0]) != str(spkeyinfo[2]):
								tempmodpepseq=str(spkeyinfo[2])
							if len(accver.strip())>0:
								outputresultpassel.write((str(accessid)+'-'+str(accver))+'\t'+str(spkeyinfo[0])+'\t'+str(tempmodpepseq)+'\t'+str(spkeyinfo[1])+'\t'+str(transdataFinal)+'\n')
							else:
								outputresultpassel.write(str(accessid)+'\t'+str(spkeyinfo[0])+'\t'+str(tempmodpepseq)+'\t'+str(spkeyinfo[1])+'\t'+str(transdataFinal)+'\n')
		outputresultpassel.close()

	def multiProcessMapData(self):
		self.multiProcessGetPasselData()
		print('Passel Data Extraction Done')
		self.createPasselObject()
		if len(self.passeldic) >0:
			fileNames=[]
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
					fileNames.append(self.finalpasselFile+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.mapId,args=[[counter,passelPrimaryProtID[i:]]])
					p.start()
					processes.append(p)
					break
				elif counter < self.numberOfprocess:
					fileNames.append(self.finalpasselFile+'_Part'+str(counter)+'.csv')
					p=multiprocessing.Process(target=self.mapId,args=[[counter,passelPrimaryProtID[i:i+bins]]])
					p.start()
					processes.append(p)
					print(i,i+bins)
			for process in processes:
				process.join()
			finish=time.time()
			print('Finished in '+str(finish-start)+' seconds for get Passel Mapping ID job')
			#increase the field size of CSV
			csv.field_size_limit(int(ctypes.c_ulong(-1).value // 2))
			df = pd.concat((pd.read_csv(self.passelpath+f, header = 0, sep='\t',keep_default_na=False) for f in fileNames))
			df_deduplicated = df.drop_duplicates()
			df_deduplicated.to_csv(self.passelpath+self.finalpasselFile+".csv",sep='\t', encoding='utf-8',index=False)
			for f in fileNames:
				if os.path.isfile(self.passelpath+f):
					os.remove(self.passelpath+f)
			movefilepath=os.path.join(self.passelpath,self.finalpasselFile+".csv")
			filepath = os.path.join(self.curr_dir, self.finalpasselFile+".csv")
			shutil.move(movefilepath,filepath)
			print('Passel Final Report Done')