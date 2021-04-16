import sys
import os
import pandas as pd
import yaml
# import glob
# import shutil

# define file related functions
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

CONFIGFILE = str(workflow.overwrite_configfiles[0])

## set memory limit 
## used for sambamba sort, etc
MEMORYG="100G"

#resouce absolute path
WORKDIR=config['workdir']
SCRIPTSDIR=config['scriptsdir']
RESOURCESDIR=config['resourcesdir']
if not os.path.exists(join(WORKDIR,"fastqs")):
	os.mkdir(join(WORKDIR,"fastqs"))
if not os.path.exists(join(WORKDIR,"results")):
	os.mkdir(join(WORKDIR,"results"))
for f in ["samples", "tools"]:
	check_readaccess(config[f])

SAMPLESDF = pd.read_csv(config["samples"],sep="\t",header=0,index_col="sampleName")
SAMPLES = list(SAMPLESDF.index)
SAMPLESDF["R1"]=join(RESOURCESDIR,"dummy")
SAMPLESDF["R2"]=join(RESOURCESDIR,"dummy")
SAMPLESDF["PEorSE"]="PE"

for sample in SAMPLES:
	R1file=SAMPLESDF["path_to_R1_fastq"][sample]
	R2file=SAMPLESDF["path_to_R2_fastq"][sample]
	# print(sample,R1file,R2file)
	check_readaccess(R1file)
	R1filenewname=join(WORKDIR,"fastqs",sample+".R1.fastq.gz")
	if not os.path.exists(R1filenewname):
		os.symlink(R1file,R1filenewname)
	SAMPLESDF.loc[[sample],"R1"]=R1filenewname
	if str(R2file)!='nan':
		check_readaccess(R2file)
		R2filenewname=join(WORKDIR,"fastqs",sample+".R2.fastq.gz")
		if not os.path.exists(R2filenewname):
			os.symlink(R2file,R2filenewname)
		SAMPLESDF.loc[[sample],"R2"]=R2filenewname
	else:
		SAMPLESDF.loc[[sample],"PEorSE"]="SE"

## Load tools from YAML file
check_readaccess(config["tools"])
with open(config["tools"]) as f:
	TOOLS = yaml.safe_load(f)

## Load cluster.json
check_readaccess(config["clusterjson"])
with open(config["clusterjson"]) as json_file:
    CLUSTER = json.load(json_file)

getthreads=lambda rname:int(CLUSTER[rname]["threads"]) if rname in CLUSTER else int(CLUSTER["__default__"]["threads"])

GENOME=config["genome"]
INDEXDIR=config[GENOME]["indexdir"]
