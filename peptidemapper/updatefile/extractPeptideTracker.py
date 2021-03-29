#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pymssql
import unicodedata
import datetime,glob
import time
import os,subprocess,psutil,re,sys,shutil,time
import urllib
import csv

class ExtractPeptideTrackerData():
	def __init__(self):
		self.curr_dir=os.getcwd()
		self.finalPeptideTrackerResult= 'final_report_peptrack_data.csv'
		self.RETRY_TIME = 20.0
		self.output_peptrack_file= 'allPepTrackResult.csv'
		self.spresdic={}
		self.status=['active','optimization']
		self.modpepseqdic={}
		self.peptrackurl='http://peptidetracker.proteincentre.com/rest/api/?searchtype=Any&searchterm=PEP'


	def fectchTrackerData(self):
		connection = pymssql.connect(host='65.61.234.70', port=1433, user='proteincentre_db_admin', password='y2EjiriWEK', database='proteomics')
		cursor = connection.cursor()
		cursor.execute ("SELECT PeptideModifications.Position,Peptides.PeptideName,Peptides.PeptideSequence,\
					Peptides.Status,Peptides.IsDeleted,PeptideModifications.IsDeleted,modificationtype.TotalMassShift,modificationtype.IsDeleted FROM PeptideModifications,Peptides,modificationtype \
					WHERE PeptideModifications.PeptideID = Peptides.PeptideID AND modificationtype.ModificationTypeID =PeptideModifications.ModificationTypeID AND \
					modificationtype.ResidueTypeID =PeptideModifications.ResidueID ")

		for rowdata in  cursor.fetchall():
			tempsqllist=list(rowdata)
			sqllist=[]
			for item in tempsqllist:
				if isinstance(item,unicode):
					sqllist.append(unicodedata.normalize('NFKD', item).encode('ascii','ignore').strip())
				else:
					sqllist.append(str(item).strip())
			if sqllist[3].split('-')[0].strip().lower() in self.status and sqllist[4]=='False' and sqllist[5]=='False' and sqllist[7]=='False':
				pepid=str(sqllist[1]).strip()
				pepseq=str(sqllist[2]).strip()
				pos=str(sqllist[0]).strip()
				moddata=cursor.fetchone()
				massshift=str(sqllist[6]).strip()
				try:
					residue=pepseq[int(pos)-1]
					modresrep="[+"+str(massshift)+".0]"
					sepres=str(residue)+'+'+str(massshift)
					datalist=[pepseq,pos,modresrep]
					if pepid in self.spresdic:
						self.spresdic[pepid].append(datalist)
					else:
						self.spresdic[pepid]=[datalist]
				except IndexError:
					print("Specail residue error",'\t'.join(sqllist[:5]))
		connection.commit()
		connection.close()

	def createModSeqObject(self):
		self.fectchTrackerData()
		if len(self.spresdic)>0:
			for modkey in self.spresdic:
				temppepseq=self.spresdic[modkey][0][0]
				temppepseqlist=list(temppepseq)
				countinsert=0
				for moditem in self.spresdic[modkey]:
					temppepseqlist.insert(int(moditem[1])+countinsert,moditem[2])
					countinsert+=1
				if ''.join(temppepseqlist) == temppepseq:
					self.modpepseqdic[modkey]='NA'
				else:
					self.modpepseqdic[modkey]=''.join(temppepseqlist)


	def fixedFragIon(self):
		reToFind='[a-z]\\+.*?\\d+'
		re1='([a-z])'	# Any Single Word Character (Not Whitespace) 1
		re2='(\\++)'	# Any Single Character 1
		re3='.*?'	# Non-greedy match on filler
		re4='(\\d+)'	# Integer Number 1
		rg = re.compile(reToFind,re.IGNORECASE|re.DOTALL)
		rg2 = re.compile(re1+re2+re3+re4,re.IGNORECASE|re.DOTALL)


		colName=[]
		data=[]
		tempFragion={}
		finalResults=[]
		with open(self.finalPeptideTrackerResult) as csvfile:
			reader = csv.DictReader(csvfile,delimiter='\t')
			colName=reader.fieldnames
			for row in reader:
				templist=[]
				for c in colName:
					templist.append(row[c])
				data.append(templist)
				if str(row['PeptideID'])[0:3] =='PEP':
					transData=str(row['Transitions'])
					m = rg.findall(transData)
					if len(m) >0:
						m=map(str,list(set(m)))
						for i in m:
							mg = rg2.search(i)
							if mg:
								w1=mg.group(1)
								c1=mg.group(2)
								int1=mg.group(3)
								tempFragion[i]=w1+int1+c1
		finalResults.append(colName)

		for d in data:
			for f in tempFragion:
				d[-1]=d[-1].replace(f,tempFragion[f])
			finalResults.append(d)

		with open(self.finalPeptideTrackerResult,'wb') as trackf:
			twriter=csv.writer(trackf,delimiter='\t')
			twriter.writerows(finalResults)


	def fetchDataFromPeptidetracker(self):
		self.createModSeqObject()
		while True:
			try:
				u=urllib.URLopener()
				f=u.open(self.peptrackurl)
				p=open(self.output_peptrack_file,'w')
				p.write(f.read())
				p.close()
				if os.path.isfile(self.output_peptrack_file):
					outputresultpeptrack=open(self.finalPeptideTrackerResult,'w')
					outputresultpeptrack.write('UniprotID'+'\t'+'PeptideSeq'+'\t'+'Peptide Modified Sequence'+'\t'+'PeptideID'+'\t'+'Transitions'+'\n')
					with open(self.output_peptrack_file) as csvfile:
						reader = csv.DictReader(csvfile)
						for row in reader:
							if str(row['Peptide ID'])[0:3] =='PEP':
								if str(row['Peptide ID']) in self.modpepseqdic:
									outputresultpeptrack.write(str(row['UniProtKB Accession'])+'\t'+str(row['Peptide Sequence'])+'\t'+str(self.modpepseqdic[str(row['Peptide ID'])])+'\t'+str(row['Peptide ID'])+'\t'+str(row['Transitions'])+'\n')
								else:
									outputresultpeptrack.write(str(row['UniProtKB Accession'])+'\t'+str(row['Peptide Sequence'])+'\t'+'NA'+'\t'+str(row['Peptide ID'])+'\t'+str(row['Transitions'])+'\n')
							else:
								if str(row['Special Residues']).strip().upper() != 'NA' and len(str(row['Special Residues']).strip())>0:
									epeptempspres=str(row['Special Residues'])
									epeptempres=epeptempspres.split('+')[0].strip()
									epeptempmassshift=epeptempspres.split('+')[1].strip()
									epeptempmodresrep="[+"+str(epeptempmassshift)+".0]"
									epepcountresoccur=str(row['Peptide Sequence']).count(epeptempres)
									epeptemppepseq=str(row['Peptide Sequence'])
									epeptemppepseqlist=list(epeptemppepseq)
									if epepcountresoccur >1:
										outputresultpeptrack.write(str(row['UniProtKB Accession'])+'\t'+str(row['Peptide Sequence'])+'\t'+str(row['Special Residues'])+'\t'+str(row['Peptide ID'])+'\t'+str(row['Transitions'])+'\n')
									else:
										epeplocres=str(row['Peptide Sequence']).index(epeptempres)
										epeptemppepseqlist.insert(int(epeplocres),epeptempmodresrep)
										outputresultpeptrack.write(str(row['UniProtKB Accession'])+'\t'+str(row['Peptide Sequence'])+'\t'+str(''.join(epeptemppepseqlist))+'\t'+str(row['Peptide ID'])+'\t'+str(row['Transitions'])+'\n')
								else:
									outputresultpeptrack.write(str(row['UniProtKB Accession'])+'\t'+str(row['Peptide Sequence'])+'\t'+'NA'+'\t'+str(row['Peptide ID'])+'\t'+str(row['Transitions'])+'\n')
					outputresultpeptrack.close()
					self.fixedFragIon()
					os.remove(self.output_peptrack_file)
				break
			except IOError:
				time.sleep(self.RETRY_TIME)
				print('Hey, I am trying again until succeeds to get data from PeptideTracker!',str(datetime.datetime.now()))
				pass

