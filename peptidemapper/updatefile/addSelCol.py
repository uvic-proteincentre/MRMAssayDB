import os,subprocess,psutil,re,shutil,datetime,sys,glob
import csv

homedir = os.path.normpath(os.getcwd() + os.sep + os.pardir)
filename='ReportBook_mother_file.csv'
filepath = os.path.join(homedir, 'src/mappermotherfile', filename)
addselfile="addsel.csv"
aoutput= open(addselfile,'w')
with open(filepath,'r') as f:
	for line in f:
		info=line.rstrip().split('\t')
		if 'UniProtKB Accession' in info:
			info.insert(0,"sel")
			aoutput.write(('\t'.join(info))+'\n')
		else:
			info.insert(0,"")
			aoutput.write(('\t'.join(info))+'\n')
aoutput.close()
movefilepath=os.path.join(homedir, 'updatefile', filename)
os.rename(addselfile,filename)
shutil.move(movefilepath,filepath)