import sys
import os
import Sniper
import argparse
import numpy as np

prs = argparse.ArgumentParser()
prs.add_argument('-input-list', '--input', help='Input esd file list')
prs.add_argument('-output', '--output', help='output file')
prs.add_argument('-corr-list', '--corr', help='Input list of rtraw files', default=None)

args = prs.parse_args()
cwd=os.getcwd()

outfilename = args.output
listname = args.input
corrlist = args.corr

Sniper.loadDll("/storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/BiPo212_analyzer/BiPo212_reader_cxx.so")
#Sniper.loadDll("libSimEvent.so")

task = Sniper.Task("task")
task.setLogLevel(1)

alg = task.createAlg("BiPo212_reader")

import BufferMemMgr
bufMgr = task.createSvc("BufferMemMgr")

File = "/cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/Jlatest/junosw/OEC/OECTutorial/share/DummyCommonConfig_1000t.json"
import OECComJSONSvc
oeccomjsonsvc = task.createSvc('OECComJSONSvc')
oeccomjsonsvc.property("OECReadComJSONFrom").set(0)
oeccomjsonsvc.property("OECComJSONFile").set(File)

import OECTagSvc
oectagsvc = task.createSvc('OECTagSvc')

import RootWriter
task.property("svcs").append("RootWriter")
rw = task.find("RootWriter")
rw.property("Output").set({"tree":outfilename})

import RootIOSvc
import RootIOTools
riSvc = task.createSvc("RootInputSvc/InputSvc")
inputFileNumpy = np.loadtxt(listname,usecols=(0),unpack=True,dtype=str)
inputFileList = inputFileNumpy.tolist()
riSvc.property("InputFile").set(inputFileList)

if (corrlist != None):
    corrFileNumpy = np.loadtxt(corrlist,usecols=(0),unpack=True,dtype=str)
    corrFileList = corrFileNumpy.tolist()
    riSvc.property("InputCorrelationFile").set(corrFileList)


task.setEvtMax(-1)
task.show()
task.run()
