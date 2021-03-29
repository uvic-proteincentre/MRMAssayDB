import os,subprocess,psutil,re,shutil,datetime,sys,glob
import urllib,urllib2,urllib3
import errno
from Bio import SeqIO
import random, time
import csv
import pandas as pd
import requests
from collections import Counter
from itertools import combinations
import pickle
import cPickle
import operator
from sh import gunzip

def protUnqCheck(unqisocheckdic,outfilefileUnqIsoname):

	# filepathCanonicalGZ = os.path.join(os.getcwd(), 'uniprot_sprot.fasta.gz')
	# filepathCanonical = os.path.join(os.getcwd(), 'uniprot_sprot.fasta')
	# try:
	# 	urllib.urlretrieve('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz',filepathCanonicalGZ)
	# 	urllib.urlcleanup()
	# except:
	# 	print ("Can't able to download uniprot_sprot.fasta.gz file!")

	# filepathIsoGZ = os.path.join(os.getcwd(), 'uniprot_sprot_varsplic.fasta.gz')
	# filepathIso = os.path.join(os.getcwd(), 'uniprot_sprot_varsplic.fasta')
	# try:
	# 	urllib.urlretrieve('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz',filepathIsoGZ)
	# 	urllib.urlcleanup()
	# except:
	# 	print ("Can't able to download uniprot_sprot_varsplic.fasta.gz file!")
	# if os.path.exists(filepathCanonical):
	# 	os.remove(filepathCanonical)
	# if os.path.exists(filepathIso):
	# 	os.remove(filepathIso)
	# gunzip(filepathCanonicalGZ)
	# gunzip(filepathIsoGZ)

	# canonicalSeq = open(filepathCanonical)
	# canonicalSeq_contents = canonicalSeq.read()
	# canonicalSeq.close()

	# isoSeq = open(filepathIso)
	# isoSeq_contents = isoSeq.read()
	# isoSeq.close()

	# finalFastaFIle = open("uniprot-reviewed_yes.fasta", "w") # open in `w` mode to write
	# finalFastaFIle.write(canonicalSeq_contents + isoSeq_contents) # concatenate the contents


	fastaDataDic={}
	for seq in SeqIO.parse("uniprot-reviewed_yes.fasta", "fasta"):
		tempfastaseq=str((seq.seq).strip())
		headerID=seq.id.strip().split('|')[1]
		orgID=seq.description.strip().split('OX=')[1].split(' ')[0]
		if orgID in fastaDataDic:
			fastaDataDic[orgID].append(headerID+'_'+tempfastaseq)
		else:
			fastaDataDic[orgID]=[headerID+'_'+tempfastaseq]

	countProt=0
	
	outfilefileUnqIso = open(outfilefileUnqIsoname,'w')
	outfilefileUnqIso.write('UniProtKB Accession'+'\t'+'Peptide Sequence'+'\t'+'Unique in protein'+'\t'+'Present in isoforms'+'\n')
	for mkey in unqisocheckdic.keys():
		pepunid=mkey.split('_')[0]
		unqtemppepseqList=list(set(unqisocheckdic[mkey]))
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
			tempFastaDicData=fastaDataDic[uqorgid]
			countProt+=1
			if countProt%1000 ==0:
				print str(countProt), "th protein peptide uniqueness job starts",str(datetime.datetime.now())
				time.sleep(10)
			
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
	print "Checking uniqueness of peptide sequence and presence in isoforms, job done",str(datetime.datetime.now())


