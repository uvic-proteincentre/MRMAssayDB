import pandas as pd
from selenium import webdriver
import os,subprocess,psutil,re,shutil,datetime,sys,glob,time
import names
import csv

curr_dir = os.getcwd()
pepdic={}
finalresult=[]
userxl = pd.ExcelFile("201081024_PeptidePicker_Mouse_OctFails.xlsx")
userdframe = pd.read_excel(userxl)
userheader=list(userdframe.columns.values)
#userheader=userheader[:45]
userheader=[str(s) for s in userheader]
userheader.insert(3,'Gravy score')
userheader.insert(3,'SRM/MRM Compatibility')
userheader.insert(3,'Synthesis/ Purification')
finalresult.append(userheader)
for uindex, urow in userdframe.iterrows():
	templist=[]
	for i in range(0,len(urow)):
		templist.append(str(urow[i]))
	if len([s for s in map(str.lower,templist) if 'no peptides found' in s]):
		templist.insert(3,'NA')
		templist.insert(3,'NA')
		templist.insert(3,'NA')
		finalresult.append(templist)
	if pepdic.has_key(templist[2].strip()):
		pepdic[templist[2].strip()].append(templist)
	else:
		pepdic[templist[2].strip()]=[templist]
peptidelist=pepdic.keys()
peptidelist=list(set(peptidelist))
# peptidelist=peptidelist[:10]
if len(peptidelist)>0:
	mydate = datetime.datetime.now()
	for x in range(0,len(peptidelist),100):
		temppeplist=peptidelist[x:x+100]
		nameFolder=names.get_first_name()
		folder_name=nameFolder+"_"+mydate.strftime("%B_%d_%Y_%H_%M_%S")
		os.makedirs('./downloadfile/'+folder_name)
		downloadfilepath=curr_dir+'/downloadfile/'+folder_name
		chrome_options = webdriver.ChromeOptions()
		chrome_options.add_argument("--allow-running-insecure-content")
		# chrome_options.add_argument('--headless')
		preferences = {"download.default_directory": downloadfilepath ,"directory_upgrade": True,"safebrowsing.enabled": True }
		chrome_options.add_experimental_option("prefs", preferences)
		driverpath=curr_dir +'/driver/chrome/linux/chromedriver'
		driver = webdriver.Chrome(chrome_options=chrome_options,executable_path=driverpath)
		driver.get("https://www.thermofisher.com/ca/en/home/life-science/protein-biology/peptides-proteins/custom-peptide-synthesis-services/peptide-analyzing-tool.html")
		#enter peptide sequence
		driver.find_element_by_xpath(".//textarea").send_keys('\n'.join(temppeplist))
		#click analyze button
		driver.find_element_by_xpath('//*[@id="peptide-tool"]/div[1]/form/fieldset/button[1]').click()
		time.sleep(5)
		#download file

		if len(temppeplist)==1:
			driver.find_element_by_xpath('//*[@id="peptide-tool"]/div[3]/div/div[2]/button[1]').click()
		else:
			driver.find_element_by_xpath('//*[@id="peptide-tool"]/div[3]/div[1]/div[2]/button[1]').click()
		time.sleep(5)
		driver.quit()

		thermofile=os.path.join(downloadfilepath, 'peptide_analysis.xlsx')
		srmexpdic={1:"Difficult(hydrophilic)",2:"Moderate(hydrophilic)",3:"Good",4:"Moderate(hydrophobic)",5:"Difficult(hydrophobic)"}
		synthpudic={"A":"Good","B":"Moderate","C":"Difficult"}
		#Loading the workbook into python
		oxl = pd.ExcelFile(thermofile)
		odframe = pd.read_excel(oxl, sheetname='Results')
		for oindex, orow in odframe[2:].iterrows():
			if orow[0] in pepdic:
				for item in pepdic[orow[0]]:
					item.insert(3,str(orow[7]))
					item.insert(3,srmexpdic[int(str(orow[4]))])
					item.insert(3,synthpudic[str(orow[3])])
					finalresult.append(item)

		shutil.rmtree(downloadfilepath)

	outputfilename='result_201081024_PeptidePicker_Mouse_OctFails_'+mydate.strftime("%B_%d_%Y_%H_%M_%S")+'.csv'
	with open(outputfilename, "wb") as f:
		writer = csv.writer(f)
		writer.writerows(finalresult)
