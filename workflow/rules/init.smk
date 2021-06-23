#########################################################
# IMPORT PYTHON LIBRARIES HERE
#########################################################
import sys
import os
import pandas as pd
import yaml
# import glob
# import shutil
#########################################################

#########################################################
# FILE-ACTION FUNCTIONS 
#########################################################
def check_existence(filename):
  if not os.path.exists(filename):
    exit("File: %s does not exists!"%(filename))

def check_readaccess(filename):
  check_existence(filename)
  if not os.access(filename,os.R_OK):
    exit("File: %s exists, but cannot be read!"%(filename))

def check_writeaccess(filename):
  check_existence(filename)
  if not os.access(filename,os.W_OK):
    exit("File: %s exists, but cannot be read!"%(filename))

def get_file_size(filename):
	filename=filename.strip()
	if check_readaccess(filename):
		return os.stat(filename).st_size
#########################################################

#########################################################
# DEFINE CONFIG FILE AND READ IT
#########################################################
CONFIGFILE = str(workflow.overwrite_configfiles[0])

# set memory limit 
# used for sambamba sort, etc
MEMORYG="100G"

# read in various dirs from config file
WORKDIR=config['workdir']
RESULTSDIR=join(WORKDIR,"results")
SCRIPTSDIR=config['scriptsdir']
RESOURCESDIR=config['resourcesdir']
if not os.path.exists(join(WORKDIR,"fastqs")):
	os.mkdir(join(WORKDIR,"fastqs"))
if not os.path.exists(join(WORKDIR,"results")):
	os.mkdir(join(WORKDIR,"results"))
for f in ["samplemanifest", "tools"]:
	check_readaccess(config[f])
#########################################################


#########################################################
# CREATE SAMPLE DATAFRAME
#########################################################
# each line in the samplemanifest is a replicate
# multiple replicates belong to a sample
# currently only 1,2,3 or 4 replicates per sample is supported
REPLICATESDF = pd.read_csv(config["samplemanifest"],sep="\t",header=0,index_col="replicateName")
REPLICATES = list(REPLICATESDF.index)
REPLICATESDF["R1"]=join(RESOURCESDIR,"dummy")
REPLICATESDF["R2"]=join(RESOURCESDIR,"dummy")
REPLICATESDF["PEorSE"]="PE"

for replicate in REPLICATES:
	R1file=REPLICATESDF["path_to_R1_fastq"][replicate]
	R2file=REPLICATESDF["path_to_R2_fastq"][replicate]
	# print(replicate,R1file,R2file)
	check_readaccess(R1file)
	R1filenewname=join(WORKDIR,"fastqs",replicate+".R1.fastq.gz")
	if not os.path.exists(R1filenewname):
		os.symlink(R1file,R1filenewname)
	REPLICATESDF.loc[[replicate],"R1"]=R1filenewname
	if str(R2file)!='nan':
		check_readaccess(R2file)
		R2filenewname=join(WORKDIR,"fastqs",replicate+".R2.fastq.gz")
		if not os.path.exists(R2filenewname):
			os.symlink(R2file,R2filenewname)
		REPLICATESDF.loc[[replicate],"R2"]=R2filenewname
	else:
# only PE samples are supported by the ATACseq pipeline at the moment
		print("Only Paired-end samples are supported by this pipeline!")
		print(config["samplemanifest"]+" is missing second fastq file for "+replicate)
		exit()
		REPLICATESDF.loc[[replicate],"PEorSE"]="SE"

SAMPLES=list(REPLICATESDF.sampleName.unique())
SAMPLE2REPLICATES=dict()
for g in SAMPLES:
	SAMPLE2REPLICATES[g]=list(REPLICATESDF[REPLICATESDF['sampleName']==g].index)

# print(REPLICATESDF.columns)
# print(REPLICATESDF.sampleName)
# print(SAMPLES[0])
# print(REPLICATESDF[REPLICATESDF['sampleName']==SAMPLES[0]].index)
# print(SAMPLE2REPLICATES)
# exit()

#########################################################
# READ IN TOOLS REQUIRED BY PIPELINE
# THESE INCLUDE LIST OF BIOWULF MODULES (AND THEIR VERSIONS)
# MAY BE EMPTY IF ALL TOOLS ARE DOCKERIZED
#########################################################

## Load tools from YAML file
check_readaccess(config["tools"])
with open(config["tools"]) as f:
	TOOLS = yaml.safe_load(f)
#########################################################


#########################################################
# READ CLUSTER PER-RULE REQUIREMENTS
#########################################################

## Load cluster.json
check_readaccess(config["clusterjson"])
with open(config["clusterjson"]) as json_file:
    CLUSTER = json.load(json_file)
## Create lambda functions to allow a way to insert read-in values
## as rule directives
getthreads=lambda rname:int(CLUSTER[rname]["threads"]) if rname in CLUSTER else int(CLUSTER["__default__"]["threads"])
getmemg=lambda rname:CLUSTER[rname]["mem"] if rname in CLUSTER else CLUSTER["__default__"]["mem"]
getmemG=lambda rname:getmemg(rname).replace("g","G")
#########################################################

#########################################################
# SET OTHER PIPELINE GLOBAL VARIABLES
#########################################################

GENOME=config["genome"]
INDEXDIR=config[GENOME]["indexdir"]
GENOMEFILE=join(INDEXDIR,GENOME+".genome") # genome file is required by macs2 peak calling
check_readaccess(GENOMEFILE)
QCDIR=join(RESULTSDIR,"QC")

#########################################################
