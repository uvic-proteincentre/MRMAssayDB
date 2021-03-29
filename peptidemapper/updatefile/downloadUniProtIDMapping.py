import os,subprocess,psutil,re,shutil,datetime,sys,glob
import urllib,urllib2,urllib3
from socket import error as SocketError
import errno
import random, time
import requests
from sh import gunzip
def downloadUniprotIDMapping():
	idmppaingfilepath = os.getcwd()
	print("Extracting mapping data, job starts",str(datetime.datetime.now()))
	try:
		urllib.urlretrieve('ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.2015_03.gz',idmppaingfilepath+'/idmapping.dat.2015_03.gz')
		urllib.urlcleanup()
		print("Extracting mapping data, job done",str(datetime.datetime.now()))
	except:
		print ("Can't able to download idmapping.dat.2015_03.gz file!!")

	if os.path.exists(idmppaingfilepath+'/idmapping.dat.2015_03'):
		os.remove(idmppaingfilepath+'/idmapping.dat.2015_03')
	print("Extracting .gz data, job starts",str(datetime.datetime.now()))
	gunzip(idmppaingfilepath+'/idmapping.dat.2015_03.gz')
	print("Extracting .gz data, job done",str(datetime.datetime.now()))