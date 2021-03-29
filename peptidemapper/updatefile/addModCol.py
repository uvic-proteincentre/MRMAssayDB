import urllib,urllib2
from bioservices.kegg import KEGG
import os,subprocess,psutil,re,shutil,datetime,sys,glob
from operator import itemgetter
import numpy as np
import random, time
from itertools import count, groupby
import pandas as pd
import csv


def addModCol(filename):
	homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
	filepath = os.path.join(homedir, 'src/mappermotherfile', filename)

	modrepfile="modreprot.csv"
	moutput= open(modrepfile,'w')
	with open(filepath,'r') as repfile:
		for line in repfile:
			info=line.rstrip().split('\t')
			if 'UniProtKB Accession' in info:
				info.append("PeptideTracker URL")
				info.append("PeptideTracker TransView")
				info.append("Passel URL")
				info.append("Passel TransView")
				info.append("SRMAtlas URL")
				info.append("SRMAtlas TransView")
				info.append("CPTAC URL")
				info.append("CPTAC TransView")
				info.append("Panoramaweb URL")
				info.append("Panoramaweb TransView")
				moutput.write(('\t'.join(info))+'\n')
			else:
				info.extend(['NA']*10)
				if len(info[0].strip())>0:
					uid=(str(info[0].strip()).split('-'))[0]
					if len(info[10].strip())>0 and (str(info[10]).strip()).lower() !='na':
						subinfo=(info[10].strip()).split(',')
						urlList=[]
						for sitem in subinfo:
							urlList.append('<a target="_blank" href="http://peptidetracker.proteincentre.com/search/advanced/?protein=&uniprotkb=&gene=&pepid=' + str(sitem).strip()+ '&pepseq=&org=">' + str(sitem).strip() + '</a>')
						if len(urlList)>0:
							urldata='|'.join(urlList)
							info[27]=urldata

					if len(info[11].strip())>0  and (str(info[11]).strip()).lower() !='na':
						transinfo=(info[11].strip()).split(',')
						if len(transinfo)>1:
							subtransinfo=transinfo[1].split('|')
							transdata='NA'
							if 'agilent' in subtransinfo[0].lower():
								transdata="Precursor Ion:"+str(subtransinfo[3])+",Fragment Ion:"+str(subtransinfo[6])+"<br/>Product Ion:"+str(subtransinfo[10])
							if 'qtrap' in subtransinfo[0].lower():
								transdata="Precursor Ion:"+str(subtransinfo[3])+"<br/>Fragment Ion:"+str(subtransinfo[7])+"<br/>Product Ion:"+str(subtransinfo[4])
							if 'thermo' in subtransinfo[0].lower():
								transdata="Precursor Ion:"+str(subtransinfo[3])+"<br/>Fragment Ion:"+str(subtransinfo[5])+"<br/>Product Ion:"+str(subtransinfo[8])
							info[28]=transdata

					if len(info[12].strip())>0 and (str(info[12]).strip()).lower() !='na':
						subinfo=(info[12].strip()).split(',')
						urlList=[]
						for sitem in subinfo:
							urlList.append('<a target="_blank" href="https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetSELTransitions?stripped_peptide_sequence_constraint=&pass_identifier_constraint=' + str(sitem).strip()+ '&action=QUERY">' + str(sitem).strip() + '</a>')
						if len(urlList)>0:
							urldata='|'.join(urlList)
							info[29]=urldata

					if len(info[13].strip())>0  and (str(info[13]).strip()).lower() !='na':
						transinfo=(info[13].strip()).split(',')
						if len(transinfo)>1:
							subtransinfo=transinfo[1].split('|')
							transdata='NA'
							transdata="Precursor Ion:"+str(subtransinfo[2])+"<br/>Fragment Ion:"+str(subtransinfo[4])+"<br/>Product Ion:"+str(subtransinfo[6])
							info[30]=transdata

					if len(info[14].strip())>0 and (str(info[14]).strip()).lower() !='na':
						subinfo=(info[14].strip()).split(',')
						urlList=[]
						for sitem in subinfo:
							try:
								srmID=sitem.split('protein_name_constraint=')[1].split(';')[0]
								if 'pabst_build_id=146' in sitem.lower():
									urlList.append('<a target="_blank" href="'+sitem.strip()+'">' + str(srmID).strip() + '</a>')
								elif 'pabst_build_id=161' in sitem.lower():
									urlList.append('<a target="_blank" href="'+str(sitem).strip()+'">' + str(srmID).strip() + '</a>')
								else:
									urlList.append('<a target="_blank" href="'+str(sitem).strip()+'">' + str(srmID).strip() + '</a>')
							except IndexError:
								pass
						if len(urlList)>0:
							urldata='|'.join(urlList)
							info[31]=urldata

					if len(info[15].strip())>0  and (str(info[15]).strip()).lower() !='na':
						transinfo=(info[15].strip()).split(',')
						if len(transinfo)>1:
							subtransinfo=transinfo[1].split('|')
							transdata='NA'
							transdata="Precursor Ion:"+str(subtransinfo[1])+"<br/>Fragment Ion:"+str(subtransinfo[2])+"<br/>Product Ion:"+str(subtransinfo[4])
							info[32]=transdata

					if len(info[16].strip())>0 and (str(info[16]).strip()).lower() !='na':
						subinfo=(info[16].strip()).split(',')
						urlList=[]
						for sitem in subinfo:
							urlList.append('<a target="_blank" href="https://assays.cancer.gov/'+str(sitem).strip()+'">' + str(sitem).strip() + '</a>')
						if len(urlList)>0:
							urldata='|'.join(urlList)
							info[33]=urldata

					if len(info[17].strip())>0  and (str(info[17]).strip()).lower() !='na':
						transinfo=(info[17].strip()).split(',')
						if len(transinfo)>1:
							subtransinfo=transinfo[1].split('|')
							transdata='NA'
							if len(subtransinfo)==7:
								transdata="Precursor Ion:"+str(subtransinfo[2])+"<br/>Fragment Ion:"+str(subtransinfo[4])+"<br/>Product Ion:"+str(subtransinfo[5])
							else:
								transdata="Precursor Ion:"+str(subtransinfo[2])+"<br/>Fragment Ion:"+str(re.findall('[A-Z][^A-Z]*', str(subtransinfo[3]).upper())[0])+"<br/>Product Ion:"+str(subtransinfo[4])

							subcptactranlist=[]
							for cptactan in transinfo:
								cptactraninfo=cptactan.split('|')
								if len(subtransinfo)==7:
									subcptactranlist.append(cptactan)
								else:
									fion=str(re.findall('[A-Z][^A-Z]*', str(cptactraninfo[3]).upper())[0])
									tempfdata=cptactraninfo[3].upper()
									q1z=tempfdata[:tempfdata.index(fion)]
									cptactraninfo[3]=str(q1z)
									cptactraninfo.insert(4,fion)
									subcptactranlist.append('|'.join(cptactraninfo))

							info[17]=','.join(subcptactranlist)
							info[34]=transdata


					if len(info[18].strip())>0 and (str(info[18]).strip()).lower() !='na':
						subinfo=(info[18].strip()).split(',')
						urlList=[]
						for sitem in subinfo:
							try:
								panid='PAN'+(str((str(sitem).strip()).split('showProtein.view?id=')[1])).strip()
								urlList.append('<a target="_blank" href="'+str(sitem).strip()+'">' + panid.strip() + '</a>')
							except IndexError:
								pass
						if len(urlList)>0:
							urldata='|'.join(urlList)
							info[35]=urldata

					if len(info[19].strip())>0  and (str(info[19]).strip()).lower() !='na':
						transinfo=(info[19].strip()).split(',')
						if len(transinfo)>1:
							subtransinfo=transinfo[1].split('|')
							transdata='NA'
							if len(subtransinfo)==7:
								transdata="Precursor Ion:"+str(subtransinfo[2])+"<br/>Fragment Ion:"+str(subtransinfo[4])+"<br/>Product Ion:"+str(subtransinfo[5])
							else:
								transdata="Precursor Ion:"+str(subtransinfo[2])+"<br/>Fragment Ion:"+str(re.findall('[A-Z][^A-Z]*', str(subtransinfo[3]).upper())[0])+"<br/>Product Ion:"+str(subtransinfo[4])

							subpantranlist=[]
							for pantan in transinfo:
								pantraninfo=pantan.split('|')
								if len(subtransinfo)==7:
									subpantranlist.append(pantan)
								else:
									try:
										fion=str(re.findall('[A-Z][^A-Z]*', str(pantraninfo[3]).upper())[0])
										tempfdata=pantraninfo[3].upper()
										q1z=tempfdata[:tempfdata.index(fion)]
										pantraninfo[3]=str(q1z)
										pantraninfo.insert(4,fion)
										subpantranlist.append('|'.join(pantraninfo))
									except IndexError:
										pass

							info[19]=','.join(subpantranlist)
							info[36]=transdata

					moutput.write(('\t'.join(info))+'\n')

	moutput.close()
	movefilepath=os.path.join(homedir, 'updatefile', filename)
	os.rename(modrepfile,filename)
	shutil.move(movefilepath,filepath)