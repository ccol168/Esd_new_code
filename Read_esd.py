import argparse
import os
def make_sh (infile, RunID, corrList) :

    #options_line = " --evtmax -1 --recMethod waterphase --recVersion OEC --waterPhase --method qctr --method wp-classifytrack --global-tag WaterPhase_J25.1"
    thispath = "/storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/"
    outfile = thispath + "root/" + RunID + ".root"
    #finalfile = thispath + "analysis/" + RunID + "_muons.root" 

    with open(thispath + "sh/" + RunID + ".sh" , "w") as file:
        file.write("#!/bin/bash \n")
        file.write("export LC_ALL=C \n")
        file.write("source /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/Jlatest/setup.sh \n")
        #file.write("python ${TUTORIALROOT}/share/tut_calib2rec.py --input-list " + infile + " --output " + outfile + options_line + "  >& " + thispath + "log/" + RunID + ".log \n" )
        #file.write("cd /storage/gpfs_data/juno/junofs/users/ccoletta/Muon_reco/Muon_analyzer \n")
        if (corrList == None) :
            file.write("python /storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/BiPo212_analyzer/Analyze_BiPo212.py -input-list " + infile + " -output " + outfile + "  >& " + thispath + "log/" + RunID + "_python.log")
        else:
            file.write("python /storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/BiPo212_analyzer/Analyze_BiPo212.py -input-list " + infile + " -output " + outfile + " -corr-list " + corrList +"  >& " + thispath + "log/" + RunID + "_python.log")
    os.chmod(thispath + "sh/" + RunID + ".sh",0o774 )
    return

def make_sub (RunID) :

    thispath = "/storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/"
    subpath = thispath + "sub/" + RunID + ".sub"

    with open(subpath,"w") as file:
        file.write("universe = vanilla \n")
        file.write("executable = " + thispath + "sh/" + RunID + ".sh \n")
        file.write("log = "+ thispath + "log/" + RunID + ".log \n")
        file.write("output = "+ thispath + "out/" + RunID + ".out \n")
        file.write("error = "+ thispath + "err/" + RunID + ".err \n")
        file.write("+MaxRuntime = 86400 \n")
        file.write("ShouldTransferFiles = YES \n")
        file.write("WhenToTransferOutput = ON_EXIT \n")
        file.write("+SingularityImage = false \n")
        file.write("queue 1 \n")

    return


parser = argparse.ArgumentParser(description="Make submission for a muon reco from a list of esd files")
parser.add_argument("-inList", help="Name of the input list",required=True)
parser.add_argument("-RunID",help="Identifying name of this run",required=True)
parser.add_argument("-corrList",help="Name of the rtraw files path",default=None)

args = parser.parse_args()

make_sh(args.inList,args.RunID,args.corrList)
make_sub(args.RunID)



