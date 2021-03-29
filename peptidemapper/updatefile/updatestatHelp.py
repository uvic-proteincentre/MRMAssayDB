import os,subprocess,psutil,re,shutil,datetime,sys,glob
from operator import itemgetter
import numpy as np
import random, time
from itertools import count, groupby
import pandas as pd
import itertools
import calendar,json
homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
sys.path.append(os.path.join(homedir, 'src/pepmapperapp'))

from updatedstat import *
from totalpepassay import *
from overallstat import *
from calculationprog import *


updatedstatresult=updatedstat
pepseqdic=finalresult['pepseqdic']
prodic=finalresult['prodic']

speciesName=overallSumresult['species']

allpepassay=totalpepassay['totalpepassay']
uqpep=len(pepseqdic)
uqpro=len(prodic)
year=datetime.datetime.now().year
month=datetime.datetime.now().month
templist=[year,str(month),calendar.month_abbr[int(month)],uqpep,uqpro,len(speciesName),len(allpepassay)]
updatedstatresult.append(templist)
calfilename='updatedstat.py'
calmovefilepath=os.path.join(homedir, 'updatefile', calfilename)
calfilepath = os.path.join(homedir, 'src/pepmapperapp', calfilename)

calfileoutput=open(calfilename,'w')
calfileoutput.write("updatedstat=")
calfileoutput.write(json.dumps(updatedstatresult))
calfileoutput.close()
shutil.move(calmovefilepath,calfilepath)