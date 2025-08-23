import os
import ROOT
import argparse

parser = argparse.ArgumentParser(description="Check veto rate and save results to file.")
parser.add_argument("--folder", default="/storage/gpfs_data/juno/junofs/users/caccianijuno/new_prod/results", help="Folder containing ROOT files.")
parser.add_argument("--output", required=True, help="Output filename for the table.")
parser.add_argument("--runlist", required=True, help="Text file with list of run numbers to analyze.")
args = parser.parse_args()

folder = args.folder
output_filename = args.output
runlist_filename = args.runlist

tree_name = "CdEvents"
npe_threshold = 30000
output_table = []
total_count = 0

# Read run numbers from runlist file
with open(runlist_filename) as f:
    run_numbers_to_analyze = set(line.strip() for line in f if line.strip())

for fname in os.listdir(folder):
    if fname.startswith("RUN") and fname.endswith(".root"):
        run_number = fname.split('_')[0][3:]  # Extract RUN number
        if run_number not in run_numbers_to_analyze:
            continue
        froot = ROOT.TFile(os.path.join(folder, fname))
        tree = froot.Get(tree_name)
        count = 0
        duration = 0.0
        rate = 0.0
        if tree:
            n_entries = tree.GetEntries()
            tree.SetBranchStatus("*", 0)
            tree.SetBranchStatus("npe", 1)
            tree.SetBranchStatus("fSec", 1)
            npe_arr = []
            fSec_arr = []
            for i in range(n_entries):
                tree.GetEntry(i)
                npe_arr.append(getattr(tree, "npe", 0))
                fSec_arr.append(getattr(tree, "fSec", None))
            indices = [i for i, npe in enumerate(npe_arr) if npe > npe_threshold]
            count = len(indices)
            fSec_vals = [fSec_arr[i] for i in range(n_entries) if fSec_arr[i] is not None]
            if fSec_vals:
                duration = max(fSec_vals) - min(fSec_vals)
            rate = count / duration if duration > 0 else 0.0
        output_table.append((run_number, count, duration, rate))
        total_count += count
        froot.Close()

# Save output table to file
with open(output_filename, "w") as fout:
    fout.write("RUNNumber\t#Events_npe_gt_30000\tDuration_s\tRate_Hz\n")
    for run_number, count, duration, rate in output_table:
        fout.write(f"{run_number}\t{count}\t{duration:.2f}\t{rate:.4f}\n")
    fout.write(f"Total\t{total_count}\n")