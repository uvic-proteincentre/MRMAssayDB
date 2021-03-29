#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os,subprocess,psutil,re,sys,shutil,datetime,time
import unicodedata
import urllib,urllib2
from socket import error as SocketError
import errno
import xmltodict
import random
import csv,json
from elasticsearch import Elasticsearch,helpers,RequestsHttpConnection
import requests
import unicodedata
import pandas as pd
import ctypes

def uploadData():
	print (datetime.datetime.now())
	print ("Upload mother file job Starts")

	shortmonthdic={"Jan":"January","Feb":"February","Mar":"March","Apr":"April","May":"May","Jun":"June","Jul":"July","Aug":"August","Sep":"September","Oct":"October","Nov":"November","Dec":"December"}
	#increase the field size of CSV
	csv.field_size_limit(int(ctypes.c_ulong(-1).value // 2))

	uploadDateList=[]
	currentdate=datetime.datetime.now().strftime("%B-%d-%Y")

	uploadfilename='uploadElasticSearchInfo.txt'
	with open(uploadfilename,'r') as f:
		for l in f:
			d=l.strip()
			uploadDateList.append(d)




	homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
	filename='ReportBook_mother_file.csv'
	filepath = os.path.join(homedir, 'src/mappermotherfile', filename)

	fileCreateDateList=str(time.ctime(os.path.getmtime(filepath))).split(' ')
	fileCreateDateList=[f for f in fileCreateDateList if len(f.strip()) >0 ]
	
	fileCreateDate=str(shortmonthdic[fileCreateDateList[1]])+"-"+str(fileCreateDateList[2])+"-"+str(fileCreateDateList[-1])
	if fileCreateDate not in uploadDateList:
		uploadDateList.append(fileCreateDate)

		es = Elasticsearch(['http://localhost:9200/'],connection_class=RequestsHttpConnection)

		es.indices.delete(index='mrmassaydb-index', ignore=[400, 404])
		if not es.indices.exists(index="mrmassaydb-index"):
			customset={
				"settings": {
					"analysis": {
				  		"analyzer": {
				  			"ngram_analyzer": {
								"tokenizer": "ngram_tokenizer",
								"filter": [
									"lowercase",
									"asciifolding"
								]
				  			}
				  		},
				  		"tokenizer": {
				  			"ngram_tokenizer": {
								"type": "ngram",
								"min_gram": 1,
								"max_gram": 10000,
								"token_chars": [
									"letter",
									"digit",
									"punctuation"
								]
				  			}
				  		}
					}
				},
			  	"mappings": {
					"mrmassaydb-type": {
				  		"properties": {
				  			"UniProtKB Accession": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Protein": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Gene": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},				  			
				  			"Organism": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Organism ID": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"SubCellular": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Peptide Sequence": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Pathway Name": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Disease Name": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Go ID": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Go Name": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Go Term": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Kegg Coverage": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Assays for FDA approved Marker": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Associated with COVID-19": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			},
				  			"Drug Bank": {
								"type": "text",
								"fields": {
									"ngram": {
					  					"type": "text",
					  					"analyzer": "ngram_analyzer",
					  					"search_analyzer": "standard"
									}
								}
				  			}
				  		}
					}
			  	}
			}
			es.indices.create(index="mrmassaydb-index", body=customset, ignore=400)
			tempfile=open(filepath)
			tempreader=csv.reader(tempfile)
			totalrows=len(list(tempreader))-1 # exclude header
			tempreader=None
			tempfile.close()
			linecount=0
			jsonfileExt=0
			listOfJsonFilename=[]
			with open(filepath,'r') as f:
				jsonData=[]
				elasreader = csv.DictReader(f,delimiter="\t")
				for row in elasreader:
					linecount+=1
					jsonData.append(row)
					if linecount%100000 ==0:
						filenameJson='ReportBook_mother_file_'+str(jsonfileExt)+'.json'
						filepathJson = os.path.join(homedir, 'src/mappermotherfile', filenameJson)
						listOfJsonFilename.append(filepathJson)
						jsonfileoutput= open(filepathJson,'w')
						jsonfileoutput.write(json.dumps(jsonData))
						jsonfileoutput.close()
						jsonfileExt+=1
						jsonData=[]
					if totalrows == linecount:
						filenameJson='ReportBook_mother_file_'+str(jsonfileExt)+'.json'
						filepathJson = os.path.join(homedir, 'src/mappermotherfile', filenameJson)
						listOfJsonFilename.append(filepathJson)
						jsonfileoutput= open(filepathJson,'w')
						jsonfileoutput.write(json.dumps(jsonData))
						jsonfileoutput.close()

			count=1
			for tempJsonFile in listOfJsonFilename:
				print(tempJsonFile)
				json_file_data=open(tempJsonFile).read()
				json_read_data=json.loads(json_file_data)
				for d in json_read_data:
					print(count,tempJsonFile)
					res=es.index(index='mrmassaydb-index',doc_type='mrmassaydb-type', id=count, body=d, timeout="60s")
					res=es.get(index='mrmassaydb-index',doc_type='mrmassaydb-type', id=count)
					count=count+1
		es.indices.refresh(index="mrmassaydb-index")
		with open(uploadfilename,'w') as ufile:
			for i in uploadDateList:
				ufile.write(i+"\n")

		print (datetime.datetime.now())
		print ("Upload mother file job done")
