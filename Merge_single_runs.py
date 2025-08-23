import os
import re
import argparse
import ROOT

def is_valid_root_file(file_path, tree_name):
    """Check if the ROOT file is valid and contains a non-empty TTree."""
    try:
        f = ROOT.TFile(file_path)
        if not f or f.IsZombie():
            return False
        tree = f.Get(tree_name)
        if not tree or tree.GetEntries() == 0:
            return False
        f.Close()
        return True
    except Exception:
        return False

def merge_runs(input_folder, output_folder, run_min, run_max, tree_name="CdEvents"):
    # Regex pattern to extract run number
    pattern = re.compile(r"RUN(\d+)_.*\.root")

    runs = {}

    for filename in os.listdir(input_folder):
        match = pattern.match(filename)
        if match:
            run_number = int(match.group(1))
            if run_min <= run_number <= run_max:
                full_path = os.path.join(input_folder, filename)
                if is_valid_root_file(full_path, tree_name):
                    runs.setdefault(run_number, []).append(full_path)
                else:
                    print(f"Skipping invalid or empty file: {filename}")

    os.makedirs(output_folder, exist_ok=True)

    for run_number, file_list in runs.items():
        if not file_list:
            continue
        print(f"\nMerging {len(file_list)} files for RUN{run_number}...")

        chain = ROOT.TChain(tree_name)
        for file_path in sorted(file_list):
            chain.Add(file_path)

        output_path = os.path.join(output_folder, f"RUN{run_number}_merged.root")
        chain.Merge(output_path)
        print(f"Saved: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge ROOT files by run number.")
    parser.add_argument("--input", default="root", help="Input folder containing ROOT files.")
    parser.add_argument("--output", default="analysis", help="Output folder for merged ROOT files.")
    parser.add_argument("--run-start", type=int, required=True, help="Start run number (inclusive).")
    parser.add_argument("--run-end", type=int, required=True, help="End run number (inclusive).")

    args = parser.parse_args()
    merge_runs(args.input, args.output, args.run_start, args.run_end)

