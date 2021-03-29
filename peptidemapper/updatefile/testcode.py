#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os,subprocess,psutil,re,shutil,datetime,sys,glob
import urllib,urllib2,urllib3
from socket import error as SocketError
import errno
import random, time
import requests
from sh import gunzip
import pandas as pd
import pickle
import cPickle
import ctypes
from xml.etree import cElementTree as ET
import xmltodict
from xml.dom import minidom
from xml.parsers.expat import ExpatError
from extractDisGenData import disGenData


filenamePallab='ReportBook_mother_file.csv'
filenameYassene='ReportBook_mother_file_disgen_nows_20210128_v7_y.tsv'
homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
filepathPallab = os.path.join(homedir, 'src/mappermotherfile', filenamePallab)
filepathYassene = os.path.join(homedir, 'src/mappermotherfile', filenameYassene)

with open(filepathPallab) as fp:
	for lp in fp:
		dp=lp.strip()
		if not dp.startswith('sel'):
			ip=dp.split('\t')
			print(ip[1],ip[4],ip[5])
			break

# print('----------------------------------------')
# with open(filepathYassene) as fy:
# 	for ly in fy:
# 		dy=ly.strip()
# 		if not dy.startswith('sel'):
# 			iy=dy.split('\t')
# 			print(iy[1],iy[4],iy[5])
# 			break