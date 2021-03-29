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
	filenameJson='ReportBook_mother_file.json'
	filepath = os.path.join(homedir, 'src/mappermotherfile', filename)
	filepathJson = os.path.join(homedir, 'src/mappermotherfile', filenameJson)


	fileCreateDateList=str(time.ctime(os.path.getmtime(filepath))).split(' ')
	fileCreateDateList=[f for f in fileCreateDateList if len(f.strip()) >0 ]


	fileCreateDate=str(shortmonthdic[fileCreateDateList[1]])+"-"+str(fileCreateDateList[2])+"-"+str(fileCreateDateList[-1])

	if fileCreateDate not in uploadDateList:
		uploadDateList.append(fileCreateDate)

		es = Elasticsearch(['http://localhost:9200/'],connection_class=RequestsHttpConnection)

		## after test change index from my-index to my-index
		es.indices.delete(index='my-index', ignore=[400, 404])
		if not es.indices.exists(index="my-index"):
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
					"my-type": {
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
				  			"Kegg Pathway Name": {
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
			es.indices.create(index="my-index", body=customset, ignore=400)

			jsonData=[]
			count=1
			with open(filepath,'r') as f:
				elasreader = csv.DictReader(f,delimiter="\t")
				for row in elasreader:
					#jsonData.append(row)
					if count%100 ==0:
						print str(count), "th row upload job starts",str(datetime.datetime.now())
					tempJsonData=json.dumps(row)
					res=es.index(index='my-index',doc_type='my-type', id=count, body=tempJsonData)
					res=es.get(index='my-index',doc_type='my-type', id=count)
					count=count+1
			'''
			finalresultJson={"data":jsonData}
			with open(filepathJson, 'w') as jsonfileoutput:
				json.dump(jsonData[:10], jsonfileoutput)

			print('Json Done')
			jsonData=[]
			count=1
			json_file_data=open(filepathJson).read()
			print('Check')
			json_read_data=json.loads(json_file_data)
			print('Load')
			for d in json_read_data:
				print(d);sys.exit();
				res=es.index(index='my-index',doc_type='my-type', id=count, body=d)
				res=es.get(index='my-index',doc_type='my-type', id=count)
				count=count+1
		'''
		es.indices.refresh(index="my-index")
		with open(uploadfilename,'w') as ufile:
			for i in uploadDateList:
				ufile.write(i+"\n")
		#os.remove(filepathJson)
		print (datetime.datetime.now())
		print ("Upload mother file job done")


uploadData()