import os,subprocess,psutil,re,shutil,datetime,sys,glob
from socket import error as SocketError
import errno
import random, time
import csv
import json
import pandas as pd
from collections import Counter
from itertools import combinations
import operator
import collections

class TotalAssay():
	def __init__(self):
		self.homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
		self.filename='ReportBook_mother_file.csv'
		self.pepfilepath = os.path.join(self.homedir, 'src/mappermotherfile', self.filename)
		self.totalpepseq=[]
		self.totalassayresult={}
		self.calfilename='totalpepassay.py'
		self.unqtotalpepseq=[]
		self.calmovefilepath=os.path.join(self.homedir, 'updatefile', self.calfilename)
		self.calfilepath = os.path.join(self.homedir, 'src/pepmapperapp', self.calfilename)
		self.totalassayValid=0
		self.totalassayNonValid=0
	def totalSeq(self):
		csv.field_size_limit(sys.maxsize)
		pepresult = csv.DictReader(open(self.pepfilepath,'r'), delimiter='\t')
		for reppeprow in pepresult:
			self.totalpepseq.append(str(reppeprow['Peptide Sequence']).strip())
		self.unqtotalpepseq=list(set(self.totalpepseq))


	def assayStat(self):
		df =pd.read_csv(self.pepfilepath, delimiter='\t')
		portPepDF=pd.DataFrame(columns=["UniProtKB Accession","Peptide Sequence","UniprotKb entry status"])
		portPepDF=df[["UniProtKB Accession","Peptide Sequence","UniprotKb entry status"]]
		portPepDF=portPepDF.drop_duplicates()
		portPepDF["Canonical UniID"]= portPepDF["UniProtKB Accession"].apply(lambda x: x.split('-')[0])
		portPepDF["Peptide Sequence"]=portPepDF["Peptide Sequence"].apply(str)
		portPepDF["Peptide Sequence"]=portPepDF["Peptide Sequence"].str.lower()
		portPepDF["Canonical UniID"]=portPepDF["Canonical UniID"].apply(str)
		portPepDF["Canonical UniID"]=portPepDF["Canonical UniID"].str.lower()
		portPepDF=portPepDF[portPepDF["Canonical UniID"] !='502']
		portPepDF["UniIDPepSeq"]=portPepDF["Canonical UniID"]+portPepDF["Peptide Sequence"]
		portPepDF["UniprotKb entry status"]=portPepDF["UniprotKb entry status"].apply(str)
		portPepDF["UniprotKb entry status"]=portPepDF["UniprotKb entry status"].str.lower()
		portPepDF.drop(columns=["Canonical UniID","Peptide Sequence"],axis=1,inplace=True)
		portPepDF=portPepDF.drop_duplicates()
		self.totalassayValid=len(set(portPepDF["UniIDPepSeq"][portPepDF["UniprotKb entry status"] =='yes']))
		self.totalassayNonValid=len(set(portPepDF["UniIDPepSeq"]))

	def finalStep(self):
		self.totalSeq()
		self.assayStat()

		self.totalassayresult['totalstripPep']=self.unqtotalpepseq
		self.totalassayresult['totalassayValid']=self.totalassayValid
		self.totalassayresult['totalassayNonValid']=self.totalassayNonValid

		calfileoutput=open(self.calfilename,'w')
		calfileoutput.write("totalpepassay=")
		calfileoutput.write(json.dumps(self.totalassayresult))
		calfileoutput.close()
		shutil.move(self.calmovefilepath,self.calfilepath)
