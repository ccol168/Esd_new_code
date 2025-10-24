import argparse
from pathlib import Path
import re

FILES_PER_LIST = 100  # Number of files per output list
BASE_PATH_ESD = Path("/storage/gpfs_data/juno/junofs/production/storm/dirac/juno/juno-kup/")
BASE_PATH_RTRAW = Path("/storage/gpfs_data/juno/junofs/production/storm/dirac/juno/juno-rtraw/")

def get_valid_runs(run_list_file, version_dir_esd, version_dir_rtraw):
    """Return runs where both ESD and RTRAW directories exist and have the same number of files."""
    with open(run_list_file) as f:
        runs = [line.strip() for line in f if line.strip().isdigit()]

    valid_runs = []

    for run_str in runs:
        run = int(run_str)
        dir_lvl1 = f"{(run // 10000) * 10000:08d}"
        dir_lvl2 = f"{(run // 100) * 100:08d}"

        dir_esd = version_dir_esd / "global_trigger" / dir_lvl1 / dir_lvl2 / f"{run:05d}"
        dir_rtraw = version_dir_rtraw / "global_trigger" / dir_lvl1 / dir_lvl2 / f"{run:05d}"

        if not dir_esd.is_dir() or not dir_rtraw.is_dir():
            print(f"Skipping run {run}: missing directory")
            continue

        files_esd = list(dir_esd.glob(f"RUN.{run}.JUNODAQ.Physics.ds-2.global_trigger.*.esd"))
        files_rtraw = list(dir_rtraw.glob(f"RUN.{run}.JUNODAQ.Physics.ds-2.global_trigger.*.rtraw"))

        if len(files_esd) == 0 or len(files_rtraw) == 0:
            print(f"Skipping run {run}: no files found")
            continue

        if len(files_esd) != len(files_rtraw):
            print(f"Skipping run {run}: mismatch ({len(files_esd)} ESD vs {len(files_rtraw)} RTRAW)")
            continue

        valid_runs.append(run_str)

    print(f"\n {len(valid_runs)} valid runs found (out of {len(runs)})\n")
    return valid_runs

def generate_list_files(runs, version_dir, list_dir, file_ext):
    """Generate .list files for each run"""

    list_dir.mkdir(exist_ok=True)
    generated_lists = []

    for run_str in runs:
        run = int(run_str)

        # Correct directory naming
        dir_lvl1 = f"{(run // 10000) * 10000:08d}"  # e.g., 10554 -> 00010000
        dir_lvl2 = f"{(run // 100) * 100:08d}"      # e.g., 10554 -> 00010500

        run_dir = version_dir / "global_trigger" / dir_lvl1 / dir_lvl2 / f"{run:05d}"

        if not run_dir.is_dir():
            print(f"Warning: Directory not found for run {run}: {run_dir}")
            continue

        files = sorted(run_dir.glob(f"RUN.{run}.JUNODAQ.Physics.ds-2.global_trigger.*.{file_ext}"))
        if not files:
            print(f"Warning: No matching files for run {run}")
            continue

        for i in range(0, len(files), FILES_PER_LIST):
            chunk = files[i:i+FILES_PER_LIST]
            index = (i // FILES_PER_LIST)
            out_filename = list_dir / f"list_RUN{run}_{index:03d}.list"
            with open(out_filename, "w") as out:
                for f in chunk:
                    out.write(str(f.resolve()) + "\n")
            generated_lists.append(out_filename)
            print(f"Wrote {len(chunk)} paths to {out_filename}")

    return generated_lists

def generate_sh_file(list_file, corr_list_file, run, idx, sh_dir, output_root_dir, log_dir):
    """Generate a single .sh file for a given run/index."""
    sh_filename = sh_dir / f"run_rtraw2calib_RUN{run}_{idx}.sh"
    output_root = Path(output_root_dir) / f"RUN_{run}_{idx}.root"
    log_root = Path(log_dir) / f"RUN_{run}_{idx}.log"

    with open(sh_filename, "w") as sh:
        sh.write("#!/bin/bash\n")
        sh.write("export LC_ALL=C\n")
        sh.write("export CMTCONFIG=amd64_linux26\n\n")
        sh.write("source /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/Jlatest/setup.sh\n")

        sh.write("python /storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/BiPo212_analyzer/Analyze_BiPo212.py \\\n")
        sh.write(f"      -input-list {list_file.resolve()} \\\n")
        sh.write(f"      -corr-list {corr_list_file.resolve()} \\\n")
        sh.write(f"      -output {output_root} \\\n")
        sh.write(f"      >& {log_root}\n\n")

    sh_filename.chmod(0o755)
    print(f"Wrote {sh_filename}")
    return sh_filename

def write_sub_file(sh_file, run, idx, script_dir):
    sub_dir = script_dir / "sub"
    sub_dir.mkdir(exist_ok=True)

    sub_filename = sub_dir / f"RUN{run}_{idx}.sub"
    log_dir = script_dir / "log"
    out_dir = script_dir / "out"
    err_dir = script_dir / "err"
    log_dir.mkdir(exist_ok=True)
    out_dir.mkdir(exist_ok=True)
    err_dir.mkdir(exist_ok=True)

    log_file = log_dir / f"RUN{run}_{idx}.log"
    out_file = out_dir / f"RUN{run}_{idx}.out"
    err_file = err_dir / f"RUN{run}_{idx}.err"

    with open(sub_filename, "w") as sub:
        sub.write("universe = vanilla\n")
        sub.write(f"executable = {sh_file.resolve()}\n")
        sub.write(f"log = {log_file}\n")
        sub.write(f"output = {out_file}\n")
        sub.write(f"error = {err_file}\n")
        sub.write("+MaxRuntime = 86400\n")
        sub.write("ShouldTransferFiles = YES\n")
        sub.write("WhenToTransferOutput = ON_EXIT\n")
        sub.write("+SingularityImage = false\n")
        sub.write("queue 1\n")

    print(f"Wrote {sub_filename}")
    return sub_filename

def main():
    parser = argparse.ArgumentParser(description="Generate JUNO RTRAW .list, .sh, .sub, and master submit script.")
    parser.add_argument("--versionesd", type=str, default="J25.5.0.b", help="Software version (default J25.5.0.b)")
    parser.add_argument("--versionrtraw", type=str, default="J25.5.0", help="Software version (default J25.5.0)")
    parser.add_argument("--runlist", type=str, required=True, help="Text file with run numbers")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    list_dir_rtraw = script_dir / "list_rtraw"
    list_dir_esd = script_dir / "list"
    sh_dir = script_dir / "sh"
    sh_dir.mkdir(exist_ok=True)

    version_dir_esd = BASE_PATH_ESD / args.versionesd
    version_dir_rtraw = BASE_PATH_RTRAW / args.versionrtraw

    if not version_dir_esd.is_dir():
        print(f"Error: version directory not found: {version_dir_esd}")
        return
    if not version_dir_rtraw.is_dir():
        print(f"Error: version directory not found: {version_dir_rtraw}")
        return
    
    # Step 0: control for mismatches between rtraw and esd files
    valid_runs = get_valid_runs(Path(args.runlist), version_dir_esd, version_dir_rtraw)

    if not valid_runs:
        print("No valid runs found. Exiting.")
        return

    # Step 1: generate list files
    generated_lists_esd = generate_list_files(valid_runs, version_dir_esd, list_dir_esd, "esd")
    generated_lists_rtraw = generate_list_files(valid_runs, version_dir_rtraw, list_dir_rtraw, "rtraw")

    if not generated_lists_esd:
        print("No .list with esd files generated. Exiting.")
        return
    
    if not generated_lists_rtraw:
        print("No .list with rtraw files generated. Exiting.")
        return

    # Config paths
    OUTPUT_ROOT_DIR = "/storage/gpfs_data/juno/junofs/users/ccoletta/BiPo212/Esd_new_code/root"

    sub_files = []

    for list_file_rtraw,list_file_esd in zip(generated_lists_rtraw,generated_lists_esd) :
        m = re.search(r"RUN(\d+)_(\d+)\.list$", list_file_rtraw.name)
        if not m:
            continue
        run, idx = m.groups()

        sh_file = generate_sh_file(list_file_esd, list_file_rtraw, run, idx, sh_dir, OUTPUT_ROOT_DIR, script_dir / "log")
        sub_file = write_sub_file(sh_file, run, idx, script_dir)
        sub_files.append(sub_file)

    # Step 3: master submit script
    c_launch_dir = script_dir / "c_launch"
    c_launch_dir.mkdir(exist_ok=True)
    runlist_stem = Path(args.runlist).stem
    master_sh = c_launch_dir / f"submit_{runlist_stem}.sh"

    with open(master_sh, "w") as f:
        for sub_file in sub_files:
            batch_name = sub_file.stem
            f.write(f"condor_submit -spool -name sn01-htc.cr.cnaf.infn.it -batch-name {batch_name} {sub_file.resolve()}\n")

    master_sh.chmod(0o755)
    print(f"Wrote master submit script: {master_sh}")

if __name__ == "__main__":
    main()
