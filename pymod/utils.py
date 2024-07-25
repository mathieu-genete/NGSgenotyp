#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import bz2
import gzip
import sys
import time
import yaml
import signal
import threading
from threading import Thread, RLock
import multiprocessing
import subprocess

__version__="v2.2"
verrou = RLock()

#default values for FreeMemory Alert and CPU load Alert
gMemFreeAlert = 10000000
gCpuLoadAlert = 95.0
default_maxParallelJobs = 20

temporaryFolder="/tmp"
params_fileName = "{}_utils_params".format(os.getuid())
params_file = os.path.join(temporaryFolder,params_fileName)

#get term and int signal
def sig_handler(sig, frame):
        
        if sig==signal.SIGTERM or sig==signal.SIGINT:
                print("Terminated...")
                os.system('pkill -TERM -g {pgid}'.format(pgid=os.getpid()))
                sys.exit(0)

signal.signal(signal.SIGTERM, sig_handler)
signal.signal(signal.SIGINT, sig_handler)

class Jobslist:

	def __init__(self, name):
		self.jobs_list = []
		self.listname = name
		
	def add_a_job(self, cmd, jobtitle ,target=''):
		jobtab = [cmd, jobtitle, target]
		self.jobs_list.append(jobtab)
	
	def get_joblist(self):
		return self.jobs_list
		
	def get_SizeOfJoblist(self):
		return len(self.jobs_list)
	
	def get_ListName(self):
		return self.listname
	
	def set_joblist(self, joblist):
		self.jobs_list = joblist
	
	def split_joblist(self, factor):

		list_joblist = []

		if factor<len(self.jobs_list):
			size = len(self.jobs_list)/factor
		
			for i in range(0,len(self.jobs_list),size):
				tmp = Jobslist("{}_SPLITED".format(self.listname))
				tmp.set_joblist(self.jobs_list[i:i+size])
				list_joblist.append(tmp)
		
		else:
			sys.stderr.write("ERROR split_joblist: factor should be less than the list size ({} items)".format(len(self.jobs_list)))
		
		return list_joblist

class Trun(Thread):

	termsig = False
	
	def __init__(self, cmd, jobtitle, target='', force=False, verbose=False):

		Thread.__init__(self)

		self.cmd = cmd
		
		self.jobtitle = jobtitle
		
		self.target = target
		
		self.force = force
		
		self.verbose = verbose

	def run(self):

		msg = ''
		
		if self.target != '' and os.path.isfile(self.target):
			target_Notexist = False
			msg += "- File already exist "
		else:
			target_Notexist = True
			
		if self.force:
			msg += "- Force "
		
		if (target_Notexist or self.force) and not Trun.termsig:
		
			if self.verbose:
				with verrou:
					sys.stdout.write("{} Trun on {}\n".format(msg,self.jobtitle))		
			ret = os.system(self.cmd)
			if ret!=0:
				if ret==2 or ret==15:
					Trun.termsig = True
					sys.exit(0)
				with verrou:
					sys.stderr.write("\nerror ({}) from execution of {}".format(ret, self.cmd))
		else:
			if self.verbose:
				with verrou:
					sys.stderr.write("-- ABORD {msg} {jobtitle} --\n".format(msg=msg, jobtitle=self.jobtitle))

#return CPU load in percents
def getCPUload():

	ldavg = os.getloadavg()[1]
	nbrcpu = multiprocessing.cpu_count()
	return float((ldavg / nbrcpu)*100)
	
#return MEMORY load in percents
def getMEMload():
	meminfo = getMeminfo()
	memtotal = meminfo['MemTotal']
	memfree = meminfo['MemFree']
	return float(((float(memtotal) - float(memfree))/float(memtotal))*100)

#return Free Memory in kB
def getFreeMEM():
	meminfo = getMeminfo()
	memfree = meminfo['MemFree']
	return memfree

def getMeminfo():
	meminfo = {}
	for i in open('/proc/meminfo').readlines():
		meminfo[i.split()[0].rstrip(':')] = int(i.split()[1])
	return meminfo

def checkArgs(argument,args):
	
	argsDict = vars(args)
	if argument in argsDict:
		if argsDict[argument]:
			return True
	else:
		sys.stderr.write("DEV-WARNING: argument '{}' not present in parser\n".format(argument))
	
	return False

#Check memory and cpu load. Ask user to continue or not
def checkMemCpuLoad(args, msg, MemFreeAlert, CpuLoadAlert):
	
	rep ='Y'
	
	if checkArgs('checks',args):
	
		if (getFreeMEM()<MemFreeAlert):
			rep = str(raw_input("WARNING: High MEMORY Load {}% - {} gB Free - Continue '{}' ? (y/n) [n]:".format(round(getMEMload(),1),round(float(getMeminfo()['MemFree'])/1000000,1),msg)))
	
		if (getCPUload()>CpuLoadAlert):
			rep = str(raw_input("WARNING: High CPU Load {}% - Continue '{}'' ? (y/n) [n]:".format(round(getCPUload(),1),msg)))
	
	if rep.upper()=='Y':
		return True
	else:
		return False
		
def set_paramFileName(AppId):
	global params_file
	params_fileName = "{}_{}_utils_params".format(AppId,os.getuid())
	params_file = os.path.join(temporaryFolder,params_fileName)

def delete_paramFileName():
	if os.path.exists(params_file):
		os.remove(params_file)
	
	
def trun(args, joblstObj, waitTermin = True, MemFreeAlert = None, CpuLoadAlert = None,thrd = []):

        if MemFreeAlert==None:
                MemFreeAlert = gMemFreeAlert
        if CpuLoadAlert==None:
                CpuLoadAlert = gCpuLoadAlert

        jobslist = joblstObj.get_joblist()

        if checkMemCpuLoad(args, joblstObj.get_ListName(), MemFreeAlert, CpuLoadAlert):
                force = False
                verbose = False
                if checkArgs('force',args):
                        force = True
                if checkArgs('verbose',args):
                        verbose = True
                        print("-- Launch Jobslist {} --".format(joblstObj.get_ListName()))
                pstart = time.process_time()
                wstart = time.time()
                for jobtab in jobslist:
                        cmd = jobtab[0]
                        jobtitle = jobtab[1]
                        target = jobtab[2]
                        if checkArgs('test',args):
                                print("should make: {tgt}\nCommand: {cmd}\n\n".format(jbtitle=jobtitle, tgt=target, cmd=cmd))
                        else:
                                if checkArgs('verbose',args):
                                        print("Launch job {jbtitle}\ntarget: {tgt}\nforce : {force}\ncommand: {cmd}\n".format(jbtitle=jobtitle,tgt=target,force=force,cmd=cmd))
                                thrd.append(startJobList(cmd, jobtitle, target, force, verbose))
                                if waitTermin: waitForJobTermination(thrd)

                if not checkArgs('test',args):
                        if waitTermin: checkJobListTerminated(args, thrd)

                pend = time.process_time()
                wend = time.time()
                ptime = pend - pstart
                wtime = wend - wstart

                if checkArgs('verbose',args):
                        sys.stderr.write("\nJob list terminated in {ptime} sec (process time) -- {wtime} sec (wall time)\n\n".format(ptime=ptime, wtime=wtime))
        else:
                sys.stderr.write("MEMORY or CPU Load WARNING: Job canceled by user...\n")

        return thrd

def getMaxParallelJobsFromConfig():
	config = getConfigurationsFromFile(params_file)
	return config['MaxParallelJobs']

def setMaxParallelJobsToConfig(maxJobs):
	config = getConfigurationsFromFile(params_file)
	config['MaxParallelJobs'] = int(maxJobs)

	with open(params_file, 'w') as yaml_file:
		yaml.dump(config, yaml_file, default_flow_style=False)

def getConfigurationsFromFile(conf_file):
	if not os.path.exists(conf_file):
		config = {}
		config['MaxParallelJobs'] = default_maxParallelJobs
		with open(conf_file, 'w') as yaml_file:
			yaml.dump(config, yaml_file, default_flow_style=False)
			
	config = yaml.load(open(conf_file),Loader=yaml.FullLoader)
	return config
	
def getConfigModificationTime():
	return os.path.getmtime(params_file)
	
def startJobList(cmd, jobtitle, target, force, verbose):
	tmp = Trun(cmd, jobtitle, target, force, verbose)
	tmp.setName(jobtitle)
	tmp.start()
	return tmp
                                
def waitForJobTermination(thrdList):
	nbrmaxOfJobs = getMaxParallelJobsFromConfig()
	lastmod = getConfigModificationTime()
	
	startTime = time.time()
	
	if nbrmaxOfJobs>0:
		NbrJobsAlive = getNumberJobAlive(thrdList)
	
		while not(NbrJobsAlive<nbrmaxOfJobs):
			NbrJobsAlive = getNumberJobAlive(thrdList)
			
			if (time.time()-startTime)>5:
				startTime = time.time()
				if lastmod != getConfigModificationTime():
					nbrmaxOfJobs = getMaxParallelJobsFromConfig()
					print("--- MaxParallelJobs change to {} ---".format(nbrmaxOfJobs))

			time.sleep(0.1)

def getNumberJobAlive(thrdList):

	nbr = 0
	
	for obj in thrdList:
		if obj.is_alive():
			nbr += 1
			
	return nbr

def getNumberJobTerminated(thrdList):

	nbr = 0
	
	for obj in thrdList:
		if not obj.is_alive():
			nbr += 1
			
	return nbr

def checkJobListTerminated(args, thrdList):

	ended = False
	jobterminated = []
	
	while not ended:
		endc=0
		
		time.sleep(0.1)
		
		for obj in thrdList:
			if not obj.is_alive():
				if checkArgs('verbose',args) and (obj.getName() not in jobterminated):
					jobterminated.append(obj.getName())
					sys.stderr.write("\nOK: {} terminated\n".format(obj.getName()))
				endc+=1
				
		if endc==len(thrdList):
			ended = True
		else:
			ended = False
			

def run(args, msg, target, cmd, retOUT=False, MemFreeAlert = None, CpuLoadAlert = None):
	""" Execute a shell command that should create a target file depending on this target file presence and on arguments options
	:param args: command line arguments : verbose, test and force
	  verbose : print message
	  test : do not execute command (only print it)
	  force : execute dommand even if target already exists
	:param msg: a message to be printed if verose mode is on
	:param target: the file that the command should generate
	:param cmd: the shell comand that should generate the target file

	:type args: argparse object
	:type msg: str
	:type target: str
	:type cmd: str

	"""
	outputs = ['','']
	
	if MemFreeAlert==None:
		MemFreeAlert = gMemFreeAlert
	if CpuLoadAlert==None:
		CpuLoadAlert = gCpuLoadAlert
	
	if checkMemCpuLoad(args, msg, MemFreeAlert, CpuLoadAlert):
		
		if os.path.isfile(target):
			# target exists, look force and test
			if checkArgs('force',args):
				if checkArgs('test',args):
					# target exists -f -t
					if checkArgs('verbose',args): print("should re-make", target)
				else:
					# target exists -f
					outputs = do_it(args, msg, target, cmd, retOUT)
			else:
					if checkArgs('verbose',args): print("OK", target)
			
		else:
			# target does not exists, look only for test
			if checkArgs('test',args):
				# target does not exists -t
				if checkArgs('verbose',args): print("should make", target)
			else:
				# target does not exists
				outputs = do_it(args, msg, target, cmd, retOUT)
	else:
		sys.stderr.write("MEMORY or CPU Load WARNING: Job canceled by user...\n")
	
	return outputs

def do_it(args, msg, target, cmd, retOUT):

	if checkArgs('verbose',args):
		#print "create %s, cmd=[%s]" % (target, cmd)
		print("create", target)
		print(cmd)
		
	if not retOUT:
		ret = os.system(cmd)
		if ret!=0:
			print("erreur sur l'execution de", cmd)
		else:
			if checkArgs('verbose',args):
				print("ok", cmd)
	else:
		job = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		outputs = job.communicate()
		return outputs


def find_Magic_Bytes(filename,magic_bytes):
    with open(filename,'r',encoding="ISO-8859-1") as infile:
        file_start = infile.read(len(magic_bytes))
    if file_start.startswith(magic_bytes):
        return True
    return False

def is_gzip(filename):
    magic_bytes = "\x1f\x8b\x08"
    return find_Magic_Bytes(filename,magic_bytes)

def is_b2z(filename):
    magic_bytes = "\x42\x5a\x68"
    return find_Magic_Bytes(filename,magic_bytes)
