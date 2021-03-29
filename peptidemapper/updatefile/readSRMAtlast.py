import time
from selenium import webdriver
import unicodedata
import os,subprocess,psutil,re,sys,shutil,csv
import urllib
import html2text
from xml.etree import cElementTree as ET
import xmltodict
from xml.dom import minidom
from xml.parsers.expat import ExpatError
import itertools
from urlparse import unquote
import requests
import datetime
from os import listdir
from selenium.webdriver.chrome.options import Options

curr_dir = os.getcwd()
srmatlaspath=curr_dir+'/SrmAtlasHuYeast/'

chrome_options = webdriver.ChromeOptions()
preferences = {"download.default_directory": srmatlaspath ,"directory_upgrade": True,"safebrowsing.enabled": True }
chrome_options.add_experimental_option("prefs", preferences)
driverpath=srmatlaspath+'driver/chrome/linux/chromedriver'
os.chdir(srmatlaspath)
driver = webdriver.Chrome(chrome_options=chrome_options,executable_path=driverpath)

srmfilelist=["HumanSRMAtlasPeptidesFinalAnnotated.csv","allbackground_nossrcalc_1Da1Da.csv"]

tempSRMidDichuman={}
tempSRMidDicyeast={}
for file in srmfilelist:
	with open(file, 'r') as outfile:
		for line in outfile:
			data=line.strip()
			if not data.startswith('protein') or not data.startswith('sequence'):
				info=data.split('\t')
				if 'human' in file.lower():
					pepseq=str(info[0]).strip()
					acclist=info[-1].split('.')
					for tempid in acclist:
						accid=""
						if 'd_' in tempid:
							accid=str(tempid[2:]).strip()
						elif '_' in tempid:
							accid=str(tempid.split('_')[0]).strip()
						else:
							accid=str(tempid).strip()

						if len(accid) >0 and len(pepseq.strip())>0:
							if tempSRMidDichuman.has_key(accid):
								tempSRMidDichuman[accid].append(pepseq)
							else:
								tempSRMidDichuman[accid]= [pepseq]


				if 'human' not in file.lower():
					pepseq=str(info[1]).strip()
					accid=str(info[0]).strip()
					if len(accid.strip())>0 and len(pepseq.strip())>0:
						if tempSRMidDicyeast.has_key(accid):
							tempSRMidDicyeast[accid].append(pepseq)
						else:
							tempSRMidDicyeast[accid]= [pepseq]


humanUnidic={}
yeastUnidic={}
with open("humanData.tab","r") as hf:
	for hline in hf:
		hdata=hline.strip()
		hinfo=hdata.split('\t')
		if not hdata.startswith('yourlist'):
			if (str(hinfo[2]).strip()).lower() != "deleted.":
				uidlist=(hinfo[0].strip()).split(',')
				for x in uidlist:
					if humanUnidic.has_key(x.strip()):
						humanUnidic[x.strip()].append(hinfo[1:])
					else:
						humanUnidic[x.strip()]= [hinfo[1:]]


TempyeastUnidic={}
with open("yeastData.tab","r") as yf:
	for yline in yf:
		ydata=yline.strip()
		yinfo=ydata.split('\t')
		if not ydata.startswith('Entry'):
			if (str(yinfo[1]).strip()).lower() != "deleted.":
				uidlist=(yinfo[0].strip()).split(',')
				for x in uidlist:
					if TempyeastUnidic.has_key(x.strip()):
						TempyeastUnidic[x.strip()].append(yinfo)
					else:
						TempyeastUnidic[x.strip()]= [yinfo]

with open("YeastmapData","r") as yf:
	for yline in yf:
		ydata=yline.strip()
		yinfo=ydata.split(' ')
		if yinfo[1].strip() in TempyeastUnidic:
			for yitem in TempyeastUnidic[yinfo[1].strip()]:
				if yeastUnidic.has_key(yinfo[0].strip()):
					yeastUnidic[yinfo[0].strip()].append(yitem)
				else:
					yeastUnidic[yinfo[0].strip()]= [yitem]

validatedSRMHumandic={}
validatedSRMYeastdic={}

for tykey in tempSRMidDicyeast.keys():
	if tykey in yeastUnidic:
		PN='NA'
		GN='NA'
		OG='NA'
		yeastAccessid=''
		for tyitem in yeastUnidic[tykey]:
			PN=str(tyitem[1]).strip()
			GN=str(tyitem[2]).strip()
			OG=str(tyitem[3]).strip()
			yeastAccessid=str(tyitem[0]).strip()
			if len(yeastAccessid.strip())>0:
				if OG.lower() !='na' and len(str(yeastAccessid).strip())>0:
					for tysrmitem in tempSRMidDicyeast[tykey]:
						tempyeastlist=[yeastAccessid,PN,GN,OG,tysrmitem]
						if validatedSRMYeastdic.has_key(tykey):
							validatedSRMYeastdic[tykey].append(tempyeastlist)
						else:
							validatedSRMYeastdic[tykey]= [tempyeastlist]

for thkey in tempSRMidDichuman.keys():
	uid=''
	if '_' in thkey:
		uid=thkey.split('_')[0]
	else:
		uid=thkey
	if uid in humanUnidic:
		PN='NA'
		GN='NA'
		OG='NA'
		accid=""
		for thitem in humanUnidic[uid]:
			PN=str(thitem[1]).strip()
			try:
				GN=str(thitem[2]).strip()
			except IndexError:
				pass
			try:
				OG=str(thitem[3]).strip()
			except IndexError:
				pass
			accid=str(thitem[0]).strip()
			if OG.lower() !='na' and accid.lower() in thkey.lower() and len(str(accid).strip())>0:
				for thsrmitem in tempSRMidDichuman[thkey]:
					temphumanlist=[thkey,PN,GN,OG,thsrmitem]
					if validatedSRMHumandic.has_key(uid):
						validatedSRMHumandic[uid].append(temphumanlist)
					else:
						validatedSRMHumandic[uid]= [temphumanlist]
						

tubUnidic={}
with open("APD_MTB_all_mapped.tab","r") as tf:
	for tline in tf:
		tdata=tline.strip()
		tinfo=tdata.split('\t')
		if not tdata.startswith('yourlist'):
			if (str(tinfo[-1]).strip()).lower() != "deleted.":
				uidlist=(tinfo[0].strip()).split(',')
				for x in uidlist:
					tubUnidic[x.strip()]= str(tinfo[-1]).strip()


validatedSRMTubIDList=tubUnidic.keys()
validatedSRMHumanIDList=validatedSRMHumandic.keys()
validatedSRMYeastIDList=validatedSRMYeastdic.keys()
output_download_srm_file= 'final_report_srm_data.csv'

outputresultsrm=open(output_download_srm_file,'w')
srmheader='UniprotID'+'\t'+'PeptideSeq'+'\t'+'PeptideModifiedSequence'+'\t'+'SRMLink'+'\t'+'Transitions'
outputresultsrm.write(srmheader+'\n')
movefilepath=os.path.join(curr_dir,'SrmAtlasHuYeast' ,output_download_srm_file)
filepath = os.path.join(curr_dir, output_download_srm_file)

countProt=0
for tx in range(0,len(validatedSRMTubIDList),50):
	print validatedSRMTubIDList[tx]
	txdata='%3B'.join(validatedSRMTubIDList[tx:tx+50])
	countProt+=len(validatedSRMTubIDList[tx:tx+50])
	if countProt%2000 ==0:
		print str(countProt), "th protein tuberculosis job starts",str(datetime.datetime.now())
	srmQueryUrl="https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetTransitions?pabst_build_id=163;protein_name_constraint=;protein_file=;peptide_sequence_constraint=;peptide_file=;peptide_length=;n_highest_intensity_fragment_ions=4;n_peptides_per_protein=1000000;target_instrument=2;exclusion_range=;multimap=results;4H=1;5H=1;Hper=1;ssr_p=0.5;C=1;D=1;M=1;P=1;R=1;S=1;W=1;nQ=1;NxST=1;nE=1;nM=1;Xc=1;nX=1;bAA=0;BA=1;EC2=1;obs=2;PATR=10;min_l=7;min_p=0.2;max_l=25;max_p=0.2;min_mz=;max_mz=;exclude_ions=;speclinks=on;rt_type=7;y_ions=on;b_ions=on;C%5B160%5D=on;K%5B136%5D=on;R%5B166%5D=on;N%5B115%5D=on;M%5B147%5D=on;C%5B143%5D=on;Q%5B111%5D=on;set_current_work_group=;set_current_project_id=;QUERY_NAME=AT_GetTransitions;apply_action_hidden=;action=QUERY"
	newprotidexp="protein_name_constraint="+str(txdata)
	srmQueryUrl=srmQueryUrl.replace("protein_name_constraint=",newprotidexp)
	driver.get(srmQueryUrl)
	time.sleep(5) # Let the user actually see something!

	driver.find_element_by_xpath('//*[@id="DWNLD"]/div/input[4]').click()
	time.sleep(5)
	tfilenames = listdir(srmatlaspath)
	tsrmatlasoutputfilename=[ tfilename for tfilename in tfilenames if tfilename.endswith( ".tsv" ) and tfilename.startswith( "SRMAtlasAssays" ) ][0]
	srmresultdic={}
	if os.path.isfile(tsrmatlasoutputfilename) and len(tsrmatlasoutputfilename.strip()) >0:
		with open(tsrmatlasoutputfilename) as tcsvfile:
			treader = csv.DictReader(tcsvfile,delimiter='\t')
			for trow in treader:
				precursorIon='NA'
				fragmentIon='NA'
				q1_charge='NA'
				productIon='NA'
				q3_charge='NA'
				instrument='NA'
				tmodpepseq='NA'
				tprot='NA'
				try:
					precursorIon=str(trow['Q1_mz']).strip()
				except KeyError:
					precursorIon='NA'
				try:
					fragmentIon=str(trow['Ion']).strip()
				except KeyError:
					fragmentIon='NA'
				try:
					q1_charge=str(trow['Q1_chg']).strip()
				except KeyError:
					q1_charge='NA'
				try:
					productIon=str(trow['Q3_mz']).strip()
				except KeyError:
					productIon='NA'
				try:
					q3_charge=str(trow['Q3_chg']).strip()
				except KeyError:
					q3_charge='NA'
				try:
					instrument=str(trow['Source']).strip()
				except KeyError:
					instrument='NA'
				try:
					tmodpepseq=str(trow['Sequence']).strip()
				except KeyError:
					tmodpepseq='NA'
				try:
					tprot=str(trow['Protein']).strip()
				except KeyError:
					tprot='NA'

				ttransitiondata=[instrument,precursorIon,fragmentIon,q1_charge,productIon,q3_charge]
				tstrippepseq=filter(str.isalpha, tmodpepseq)
				if tmodpepseq.upper() ==tstrippepseq.upper():
					tmodpepseq='NA'
				srmkey=tstrippepseq+'_'+tprot+'_'+tmodpepseq
				if tprot !='NA' and tstrippepseq !='NA':
					if srmresultdic.has_key(srmkey):
						srmresultdic[srmkey].append('|'.join(ttransitiondata))
					else:
						srmresultdic[srmkey]= ['|'.join(ttransitiondata)]

	for spkey in srmresultdic.keys():
		pID=spkey.split('_')[1]
		srmseq=spkey.split('_')[0]
		modsrmpepseq=spkey.split('_')[2]
		if pID in validatedSRMTubIDList:
			funid=tubUnidic[pID].upper()
			stransdata=','.join(srmresultdic[spkey])
			stransdata='Instrument|PrecursorIon|FragmentIon|Q1 charge|ProductIon|Q3 charge,'+stransdata
			SubsrmQueryUrl="https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetTransitions?pabst_build_id=163;protein_name_constraint=;protein_file=;peptide_sequence_constraint=;peptide_file=;peptide_length=;n_highest_intensity_fragment_ions=4;n_peptides_per_protein=1000000;target_instrument=2;exclusion_range=;multimap=results;4H=1;5H=1;Hper=1;ssr_p=0.5;C=1;D=1;M=1;P=1;R=1;S=1;W=1;nQ=1;NxST=1;nE=1;nM=1;Xc=1;nX=1;bAA=0;BA=1;EC2=1;obs=2;PATR=10;min_l=7;min_p=0.2;max_l=25;max_p=0.2;min_mz=;max_mz=;exclude_ions=;speclinks=on;rt_type=7;y_ions=on;b_ions=on;C%5B160%5D=on;K%5B136%5D=on;R%5B166%5D=on;N%5B115%5D=on;M%5B147%5D=on;C%5B143%5D=on;Q%5B111%5D=on;set_current_work_group=;set_current_project_id=;QUERY_NAME=AT_GetTransitions;apply_action_hidden=;action=QUERY"
			Subnewprotidexp="protein_name_constraint="+str(pID)
			Subnewpepseqexp="peptide_sequence_constraint="+str(srmseq)
			SubsrmQueryUrl=SubsrmQueryUrl.replace("protein_name_constraint=",Subnewprotidexp)
			SubsrmQueryUrl=SubsrmQueryUrl.replace("peptide_sequence_constraint=",Subnewpepseqexp)
			outputresultsrm.write(str(funid)+'\t'+str(srmseq)+'\t'+str(modsrmpepseq)+'\t'+str(SubsrmQueryUrl)+'\t'+str(stransdata)+'\n')

	os.remove(tsrmatlasoutputfilename)

print str(countProt), "th protein tuberculosis job ends",str(datetime.datetime.now())
print "Yeast job starts",str(datetime.datetime.now())

countProt=0
for yx in range(0,len(validatedSRMYeastIDList),50):
	print validatedSRMYeastIDList[yx]
	yxdata='%3B'.join(validatedSRMYeastIDList[yx:yx+50])
	countProt+=len(validatedSRMYeastIDList[yx:yx+50])
	if countProt%2000 ==0:
		print str(countProt), "th protein yeast job starts",str(datetime.datetime.now())
	srmQueryUrl="https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetTransitions?pabst_build_id=161;protein_name_constraint=;protein_file=;peptide_sequence_constraint=;peptide_file=;peptide_length=;n_highest_intensity_fragment_ions=4;n_peptides_per_protein=1000000;target_instrument=2;exclusion_range=;multimap=results;4H=1;5H=1;Hper=1;ssr_p=0.5;C=1;D=1;M=1;P=1;R=1;S=1;W=1;nQ=1;NxST=1;nE=1;nM=1;Xc=1;nX=1;bAA=0;BA=1;EC2=1;obs=2;PATR=10;min_l=7;min_p=0.2;max_l=25;max_p=0.2;cmod_opts=light;min_mz=;max_mz=;exclude_ions=;speclinks=on;rt_type=4;y_ions=on;b_ions=on;C%5B160%5D=on;K%5B136%5D=on;R%5B166%5D=on;N%5B115%5D=on;M%5B147%5D=on;C%5B143%5D=on;Q%5B111%5D=on;E%5B111%5D=on;set_current_work_group=;set_current_project_id=;QUERY_NAME=AT_GetTransitions;apply_action_hidden=;action=QUERY"
	newprotidexp="protein_name_constraint="+str(yxdata)
	srmQueryUrl=srmQueryUrl.replace("protein_name_constraint=",newprotidexp)
	driver.get(srmQueryUrl)
	time.sleep(5) # Let the user actually see something!

	driver.find_element_by_xpath('//*[@id="DWNLD"]/div/input[4]').click()
	time.sleep(5)
	yfilenames = listdir(srmatlaspath)
	ysrmatlasoutputfilename=[ yfilename for yfilename in yfilenames if yfilename.endswith( ".tsv" ) and yfilename.startswith( "SRMAtlasAssays" ) ][0]
	srmresultdic={}
	if os.path.isfile(ysrmatlasoutputfilename) and len(ysrmatlasoutputfilename.strip()) >0:
		with open(ysrmatlasoutputfilename) as ycsvfile:
			yreader = csv.DictReader(ycsvfile,delimiter='\t')
			for yrow in yreader:
				precursorIon='NA'
				fragmentIon='NA'
				q1_charge='NA'
				productIon='NA'
				q3_charge='NA'
				instrument='NA'
				ymodpepseq='NA'
				yprot='NA'
				try:
					precursorIon=str(yrow['Q1_mz']).strip()
				except KeyError:
					precursorIon='NA'
				try:
					fragmentIon=str(yrow['Ion']).strip()
				except KeyError:
					fragmentIon='NA'
				try:
					q1_charge=str(yrow['Q1_chg']).strip()
				except KeyError:
					q1_charge='NA'
				try:
					productIon=str(yrow['Q3_mz']).strip()
				except KeyError:
					productIon='NA'
				try:
					q3_charge=str(yrow['Q3_chg']).strip()
				except KeyError:
					q3_charge='NA'
				try:
					instrument=str(yrow['Source']).strip()
				except KeyError:
					instrument='NA'
				try:
					ymodpepseq=str(yrow['Sequence']).strip()
				except KeyError:
					ymodpepseq='NA'
				try:
					yprot=str(yrow['Protein']).strip()
				except KeyError:
					yprot='NA'
				ytransitiondata=[instrument,precursorIon,fragmentIon,q1_charge,productIon,q3_charge]
				ystrippepseq=filter(str.isalpha, ymodpepseq)
				if ymodpepseq.upper() ==ystrippepseq.upper():
					ymodpepseq='NA'
				srmkey=ystrippepseq+'_'+yprot+'_'+ymodpepseq
				if yprot !='NA' and ystrippepseq !='NA':
					if srmresultdic.has_key(srmkey):
						srmresultdic[srmkey].append('|'.join(ytransitiondata))
					else:
						srmresultdic[srmkey]= ['|'.join(ytransitiondata)]

	for spkey in srmresultdic.keys():
		pID=spkey.split('_')[1]
		srmseq=spkey.split('_')[0]
		modsrmpepseq=spkey.split('_')[2]
		if pID in validatedSRMYeastdic:
			for yitem in validatedSRMYeastdic[pID]:
				funid=yitem[0]
				stransdata=','.join(srmresultdic[spkey])
				stransdata='Instrument|PrecursorIon|FragmentIon|Q1 charge|ProductIon|Q3 charge,'+stransdata
				SubsrmQueryUrl="https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetTransitions?pabst_build_id=161;protein_name_constraint=;protein_file=;peptide_sequence_constraint=;peptide_file=;peptide_length=;n_highest_intensity_fragment_ions=4;n_peptides_per_protein=1000000;target_instrument=2;exclusion_range=;multimap=results;4H=1;5H=1;Hper=1;ssr_p=0.5;C=1;D=1;M=1;P=1;R=1;S=1;W=1;nQ=1;NxST=1;nE=1;nM=1;Xc=1;nX=1;bAA=0;BA=1;EC2=1;obs=2;PATR=10;min_l=7;min_p=0.2;max_l=25;max_p=0.2;cmod_opts=light;min_mz=;max_mz=;exclude_ions=;speclinks=on;rt_type=4;y_ions=on;b_ions=on;C%5B160%5D=on;K%5B136%5D=on;R%5B166%5D=on;N%5B115%5D=on;M%5B147%5D=on;C%5B143%5D=on;Q%5B111%5D=on;E%5B111%5D=on;set_current_work_group=;set_current_project_id=;QUERY_NAME=AT_GetTransitions;apply_action_hidden=;action=QUERY"
				Subnewprotidexp="protein_name_constraint="+str(pID)
				Subnewpepseqexp="peptide_sequence_constraint="+str(srmseq)
				SubsrmQueryUrl=SubsrmQueryUrl.replace("protein_name_constraint=",Subnewprotidexp)
				SubsrmQueryUrl=SubsrmQueryUrl.replace("peptide_sequence_constraint=",Subnewpepseqexp)
				outputresultsrm.write(str(funid)+'\t'+str(srmseq)+'\t'+str(modsrmpepseq)+'\t'+str(SubsrmQueryUrl)+'\t'+str(stransdata)+'\n')

	os.remove(ysrmatlasoutputfilename)

print str(countProt), "th protein yeast job ends",str(datetime.datetime.now())
print "human job starts",str(datetime.datetime.now())
countProt=0
for hx in range(0,len(validatedSRMHumanIDList),50):
	print validatedSRMHumanIDList[hx]
	hxdata='%3B'.join(validatedSRMHumanIDList[hx:hx+50])
	countProt+=len(validatedSRMHumanIDList[hx:hx+50])
	if countProt%2000 ==0:
		print str(countProt), "th protein human job starts",str(datetime.datetime.now())

	srmQueryUrl="https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetTransitions?pabst_build_id=146;protein_name_constraint=;protein_file=;peptide_sequence_constraint=;peptide_file=;peptide_length=;n_highest_intensity_fragment_ions=4;n_peptides_per_protein=1000000;target_instrument=2;exclusion_range=;SwissProt=on;multimap=results;4H=1;5H=1;Hper=1;ssr_p=0.5;C=1;D=1;M=1;P=1;R=1;S=1;W=1;nQ=1;NxST=1;nE=1;nM=1;Xc=1;nX=1;bAA=0;BA=1;EC2=1;obs=2;PATR=10;min_l=7;min_p=0.2;max_l=25;max_p=0.2;min_mz=;max_mz=;exclude_ions=;speclinks=on;rt_type=5;y_ions=on;b_ions=on;C[160]=on;K[136]=on;R[166]=on;N[115]=on;M[147]=on;C[143]=on;Q[111]=on;E[111]=on;set_current_work_group=;set_current_project_id=;QUERY_NAME=AT_GetTransitions;apply_action_hidden=;action=QUERY"
	newprotidexp="protein_name_constraint="+str(hxdata)
	srmQueryUrl=srmQueryUrl.replace("protein_name_constraint=",newprotidexp)
	driver.get(srmQueryUrl)
	time.sleep(5) # Let the user actually see something!

	driver.find_element_by_xpath('//*[@id="DWNLD"]/div/input[4]').click()
	time.sleep(5)
	srmresultdic={}
	hfilenames = listdir(srmatlaspath)
	hsrmatlasoutputfilename=[ hfilename for hfilename in hfilenames if hfilename.endswith( ".tsv" ) and hfilename.startswith( "SRMAtlasAssays" ) ][0]
	if os.path.isfile(hsrmatlasoutputfilename) and len(hsrmatlasoutputfilename.strip()) >0:
		with open(hsrmatlasoutputfilename) as hcsvfile:
			hreader = csv.DictReader(hcsvfile,delimiter='\t')
			for hrow in hreader:
				precursorIon='NA'
				fragmentIon='NA'
				q1_charge='NA'
				productIon='NA'
				q3_charge='NA'
				instrument='NA'
				hmodpepseq='NA'
				hprot='NA'
				try:
					precursorIon=str(hrow['Q1_mz']).strip()
				except KeyError:
					precursorIon='NA'
				try:
					fragmentIon=str(hrow['Ion']).strip()
				except KeyError:
					fragmentIon='NA'
				try:
					q1_charge=str(hrow['Q1_chg']).strip()
				except KeyError:
					q1_charge='NA'
				try:
					productIon=str(hrow['Q3_mz']).strip()
				except KeyError:
					productIon='NA'
				try:
					q3_charge=str(hrow['Q3_chg']).strip()
				except KeyError:
					q3_charge='NA'
				try:
					instrument=str(hrow['Source']).strip()
				except KeyError:
					instrument='NA'
				try:
					hmodpepseq=str(hrow['Sequence']).strip()
				except KeyError:
					hmodpepseq='NA'
				try:
					hprot=str(hrow['Protein']).strip()
				except KeyError:
					hprot='NA'
				htransitiondata=[instrument,precursorIon,fragmentIon,q1_charge,productIon,q3_charge]
				hstrippepseq=filter(str.isalpha, hmodpepseq)
				if hmodpepseq.upper() ==hstrippepseq.upper():
					hmodpepseq='NA'
				srmkey=hstrippepseq+'_'+hprot+'_'+hmodpepseq
				if hprot !='NA' and hstrippepseq !='NA':
					if srmresultdic.has_key(srmkey):
						srmresultdic[srmkey].append('|'.join(htransitiondata))
					else:
						srmresultdic[srmkey]= ['|'.join(htransitiondata)]

	for spkey in srmresultdic.keys():
		pID=spkey.split('_')[1]
		srmseq=spkey.split('_')[0]
		modsrmpepseq=spkey.split('_')[2]
		if pID in validatedSRMHumandic:
			for hitem in validatedSRMHumandic[pID]:
				funid=hitem[0]
				stransdata=','.join(srmresultdic[spkey])
				stransdata='Instrument|PrecursorIon|FragmentIon|Q1 charge|ProductIon|Q3 charge,'+stransdata
				SubsrmQueryUrl="https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetTransitions?pabst_build_id=146;protein_name_constraint=;protein_file=;peptide_sequence_constraint=;peptide_file=;peptide_length=;n_highest_intensity_fragment_ions=4;n_peptides_per_protein=1000000;target_instrument=2;exclusion_range=;SwissProt=on;multimap=results;4H=1;5H=1;Hper=1;ssr_p=0.5;C=1;D=1;M=1;P=1;R=1;S=1;W=1;nQ=1;NxST=1;nE=1;nM=1;Xc=1;nX=1;bAA=0;BA=1;EC2=1;obs=2;PATR=10;min_l=7;min_p=0.2;max_l=25;max_p=0.2;min_mz=;max_mz=;exclude_ions=;speclinks=on;rt_type=5;y_ions=on;b_ions=on;C[160]=on;K[136]=on;R[166]=on;N[115]=on;M[147]=on;C[143]=on;Q[111]=on;E[111]=on;set_current_work_group=;set_current_project_id=;QUERY_NAME=AT_GetTransitions;apply_action_hidden=;action=QUERY"
				Subnewprotidexp="protein_name_constraint="+str(pID)
				Subnewpepseqexp="peptide_sequence_constraint="+str(srmseq)
				SubsrmQueryUrl=SubsrmQueryUrl.replace("protein_name_constraint=",Subnewprotidexp)
				SubsrmQueryUrl=SubsrmQueryUrl.replace("peptide_sequence_constraint=",Subnewpepseqexp)
				outputresultsrm.write(str(funid)+'\t'+str(srmseq)+'\t'+str(modsrmpepseq)+'\t'+str(SubsrmQueryUrl)+'\t'+str(stransdata)+'\n')

	os.remove(hsrmatlasoutputfilename)

outputresultsrm.close()
shutil.move(movefilepath,filepath)
driver.quit()													 
os.chdir(curr_dir)

