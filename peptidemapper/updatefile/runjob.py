import schedule
import time
import subprocess, shutil,datetime,sys
def job():
	currentdate=datetime.datetime.now().strftime("%d")
	if currentdate =="01":
		commandjob = "python generate_final_data_report_pepmap.py"
		subprocess.Popen(commandjob, shell=True).wait()


		#elastucsearch2Restartcmd="echo 'bioinfbioinf' | sudo -S service elasticsearch start"
		#subprocess.Popen(elastucsearch2Restartcmd, shell=True).wait()


		scpcmdCalFile='sshpass -p "bioinfbioinf" scp -r /home/bioinf/Desktop/MRMAssayDB/Production/peptidemapper/src/pepmapperapp/calculationprog.py bioinf@172.16.2.212:/home/bioinf/Desktop/productionSoftware/MRMAssayDB/peptidemapper/src/pepmapperapp'
		subprocess.Popen(scpcmdCalFile, shell=True).wait()

		scpcmdUpdateStatFile='sshpass -p "bioinfbioinf" scp -r /home/bioinf/Desktop/MRMAssayDB/Production/peptidemapper/src/pepmapperapp/updatedstat.py bioinf@172.16.2.212:/home/bioinf/Desktop/productionSoftware/MRMAssayDB/peptidemapper/src/pepmapperapp'
		subprocess.Popen(scpcmdUpdateStatFile, shell=True).wait()

		scpcmdTotalAssayFile='sshpass -p "bioinfbioinf" scp -r /home/bioinf/Desktop/MRMAssayDB/Production/peptidemapper/src/pepmapperapp/totalpepassay.py bioinf@172.16.2.212:/home/bioinf/Desktop/productionSoftware/MRMAssayDB/peptidemapper/src/pepmapperapp'
		subprocess.Popen(scpcmdTotalAssayFile, shell=True).wait()

		scpcmdpepdataserFile='sshpass -p "bioinfbioinf" scp -r /home/bioinf/Desktop/MRMAssayDB/Production/peptidemapper/src/pepmapperapp/pepdataSeries.json bioinf@172.16.2.212:/home/bioinf/Desktop/productionSoftware/MRMAssayDB/peptidemapper/src/pepmapperapp'
		subprocess.Popen(scpcmdpepdataserFile, shell=True).wait()

		scpcmdoverstatFile='sshpass -p "bioinfbioinf" scp -r /home/bioinf/Desktop/MRMAssayDB/Production/peptidemapper/src/pepmapperapp/overallstat.py bioinf@172.16.2.212:/home/bioinf/Desktop/productionSoftware/MRMAssayDB/peptidemapper/src/pepmapperapp'
		subprocess.Popen(scpcmdoverstatFile, shell=True).wait()

		scpcmdhumdisFile='sshpass -p "bioinfbioinf" scp -r /home/bioinf/Desktop/MRMAssayDB/Production/peptidemapper/src/UniDiseaseInfo/humsavar.txt bioinf@172.16.2.212:/home/bioinf/Desktop/productionSoftware/MRMAssayDB/peptidemapper/src/UniDiseaseInfo'
		subprocess.Popen(scpcmdhumdisFile, shell=True).wait()

		scpcmdReportFile='sshpass -p "bioinfbioinf" scp -r /home/bioinf/Desktop/MRMAssayDB/Production/peptidemapper/src/mappermotherfile/ReportBook_mother_file.csv bioinf@172.16.2.212:/home/bioinf/Desktop/productionSoftware/MRMAssayDB/peptidemapper/src/mappermotherfile'
		subprocess.Popen(scpcmdReportFile, shell=True).wait()
		
		#commandcleanjob = "python cleanData.py"
		#subprocess.Popen(commandcleanjob, shell=True).wait()

schedule.every().day.at("01:00").do(job)

while 1:
	schedule.run_pending()
	time.sleep(1)