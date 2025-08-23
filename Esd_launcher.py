import os
import argparse
import re

def find_files_with_string(directory, search_string):
    matching_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if search_string in file:
                matching_files.append(os.path.join(root, file))
    return matching_files

def extract_name(filename):
    match = re.search(r'rawfile_(\w+)\.list', filename)
    return match.group(1) if match else None

def create_c_launch (RunName) :

    directory = "/storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/list/"

    names_string = find_files_with_string(directory,"RUN" + RunName)
    names_corr_string = find_files_with_string("/storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/list_rtraw/","RUN" + RunName)
    
    if len(names_string) == 0 :
        print(f"ERROR: No list file found for RUN{RunName}")
        print(f"Skipping RUN{RunName}")
        return
    elif len(names_corr_string) == 0 :
        print(f"ERROR: no rtraw file list found for RUN{RunName}")
        print(f"Skipping RUN{RunName}")
        return
  
    c_launch_file = open("/storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/c_launch/Launch_RUN" + RunName + ".sh","w")

    for listname,listname_rtraw in zip(names_string,names_corr_string) :
        RunID = extract_name(listname)
        os.system("python Read_esd.py -inList " + listname + " -corrList "+ listname_rtraw + " -RunID " + RunID )
        #print("python Read_esd.py -inList " + listname +  " -corrList "+ listname_rtraw + " -RunID " + RunID )
        c_launch_file.write("condor_submit -spool -name sn01-htc.cr.cnaf.infn.it -batch-name " + RunID + " /storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/sub/" + RunID + ".sub\n")

    os.chmod("/storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/c_launch/Launch_RUN" + RunName + ".sh", 0o775)

    return


parser = argparse.ArgumentParser(description="Launch a run from esd files, list in this dir")
parser.add_argument("-RunName", help="Name of the run")
parser.add_argument("-DatasetFile",help="Path of a file containing the list of runs to launch")
parser.add_argument("-DatasetName",help="Name of the resulting dataset")

args = parser.parse_args()

if (args.RunName == None and args.DatasetFile==None) :
    print("ERROR: you have not passed a run name or a list of run names to launch")
    raise FileNotFoundError

if (args.DatasetFile != None and args.DatasetName == None) :
    print("ERROR: -DatasetName missing, you need to give a name to the resulting c_launch file")
    raise ValueError

if (args.DatasetFile == None) :
    create_c_launch(args.RunName)

else :
    To_launch = []
    with open(args.DatasetFile, 'r') as f:
        for line in f:
            if line.strip():  # skip empty lines
                first_value = line.split()[0]  # split by whitespace and take first item
                To_launch.append(first_value)

    with open(f"c_launch/Launch_{args.DatasetName}.sh","w") as file :
        for element in To_launch :
            print(f"Launching RUN{element}")
            file.write(f"junosub Launch_RUN{element}.sh\n")
            create_c_launch(element)

    os.chmod(f"c_launch/Launch_{args.DatasetName}.sh", 0o775)


