import unicodedata,datetime,time,glob
import os,subprocess,psutil,re,sys,shutil
import itertools
import pandas as pd
import csv
import multiprocessing
from Bio import SeqIO

class AddReviewStatus():
	def __init__(self):
		self.homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
		self.filename='ReportBook_mother_file.csv'
		self.pepfilepath = os.path.join(self.homedir, 'src/mappermotherfile', self.filename)
		self.uniprotFastaFileName="uniprot-reviewed_yes.fasta"
		self.uniReviewStatusDataDic={}

	def createFastaDic(self):
		for seq in SeqIO.parse(self.uniprotFastaFileName, "fasta"):
			headerID=seq.id.strip().split('|')[1]
			review=seq.description.strip().split('|')[0]
			self.uniReviewStatusDataDic[headerID]=review.lower()

	def updateColumn(self):
		df= pd.read_csv(self.pepfilepath, delimiter='\t',keep_default_na=False)
		df['Reviewed']='No'
		df['Reviewed']=df.apply(lambda row: self.modifydata(row), axis=1)
		df.to_csv(self.pepfilepath, index=False, sep="\t", na_rep='NA')

	def modifydata(self,row):
		if row['UniProtKB Accession'] in self.uniReviewStatusDataDic:
			if self.uniReviewStatusDataDic[row['UniProtKB Accession']] =='sp':
				return 'Yes'
			else:
				return 'No'
		else:
			return 'No'


	def finalSteps(self):
		self.createFastaDic()
		self.updateColumn()
