"""
XNAT Incremental Downloader, DICOM Organizer, and Completeness Checker

This script connects to a specified XNAT project, identifies experiments that
have not yet been downloaded locally, retrieves them, reorganizes associated
resource files, standardizes DICOM folder structures, and generates
two reports: modified_data.csv and missing_sequences.csv.

Core functionalities:
- Delta download based on local_experiments.txt
- Automatic folder renaming based on Accession Number (ANCID)
- Behavioral resource classification (VFT, TRENDS, Examcards)
- DICOM completeness validation against expected counts
Author: Deepankan, ADC
Date: 29th Sept, 2025

"""

import xnat
import os
import urllib.parse
import shutil
from pydicom import dcmread as dr
import pandas as pd
from collections import defaultdict

# ============================================================================
# CONFIGURATION
# ============================================================================

PROJECT_ID = "YOUR_PROJECT_ID"
XNAT_URL = "http://YOUR_XNAT_HOST:PORT/xnat"
XNAT_USER = "YOUR_USERNAME"
XNAT_PASS = "YOUR_PASSWORD"

OUTPUT_DIR = "/path/to/output/DICOM_directory"
LOCAL_EXP_FILE = "local_experiments.txt"

# Expected number of files for each scan type
SCANS_TO_KEEP = {
    'dki': 4480,
    'ref_rest_ap': 88,
    'ref_dwi_pa': 112,
    'ref_dwi_ap': 112,
    'ref_rest_pa': 88,
    'ref_trends_pa': 84,
    'ref_trends_ap': 84,
    't2w': 136,
    'flair': 150,
    'dti_6dir': 448,
    'fieldmap': 176,
    'task-rest_bold': 12100,
    'task-trends_bold': 6720,
    't1w': 192,
    'ref_vft_ap': 88,
    'ref_vft_pa': 88,
    'task-vft_bold': 4752,
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def flatten_and_move(src_dir, dest_dir):
    """
    Move all contents from src_dir into dest_dir and remove src_dir.
    Ensures final scan folders do not contain nested directories.
    """
    if not os.path.exists(src_dir):
        return

    os.makedirs(dest_dir, exist_ok=True)

    for root, dirs, files in os.walk(src_dir):
        for f in files:
            src_file = os.path.join(root, f)
            dst_file = os.path.join(dest_dir, f)

            if not os.path.exists(dst_file):
                shutil.move(src_file, dst_file)

    shutil.rmtree(src_dir, ignore_errors=True)


def count_dicom_files(scan_path):
    """Return the number of DICOM (.dcm) files in a directory."""
    if not os.path.isdir(scan_path):
        return 0
    return len([f for f in os.listdir(scan_path) if f.endswith(".dcm")])


# ============================================================================
# LOAD PREVIOUSLY DOWNLOADED EXPERIMENTS
# ============================================================================

if os.path.exists(LOCAL_EXP_FILE):
    with open(LOCAL_EXP_FILE, "r") as f:
        local_experiments = set(line.strip() for line in f.readlines())
else:
    local_experiments = set()

print(f"Found {len(local_experiments)} experiments locally (tracked in text file).")

# ============================================================================
# XNAT CONNECTION AND DELTA DOWNLOAD
# ============================================================================

with xnat.connect(XNAT_URL, user=XNAT_USER, password=XNAT_PASS) as session:
    project = session.projects[PROJECT_ID]

    # Retrieve all experiments in the XNAT project
    xnat_experiments = {
        exp.label: exp
        for subj in project.subjects.values()
        for exp in subj.experiments.values()
    }

    print(f"\nFound {len(xnat_experiments)} experiments on XNAT.")

    # Determine which experiments are not yet downloaded
    missing_experiments = [
        label for label in xnat_experiments.keys()
        if label not in local_experiments
    ]

    print(f"\nMissing experiments to download: {len(missing_experiments)}")
    print(missing_experiments)

    # Download and process each missing experiment
    for exp_label in missing_experiments:
        exp = xnat_experiments[exp_label]
        exp_dir = os.path.join(OUTPUT_DIR, exp.label)
        os.makedirs(exp_dir, exist_ok=True)

        print(f"\n--- Downloading experiment: {exp.label} ---")

        scans_root = os.path.join(exp_dir, "scans")
        os.makedirs(scans_root, exist_ok=True)

        # Temporary folder for raw resource downloads
        tmp_save_path = os.path.join(exp_dir, "tmp_download")
        os.makedirs(tmp_save_path, exist_ok=True)

        # Resource download
        for res in exp.resources.values():
            res_label_raw = urllib.parse.unquote(res.label)
            tmp_resource_path = os.path.join(tmp_save_path, res_label_raw)

            os.makedirs(tmp_resource_path, exist_ok=True)
            print(f"  Resource -> {res.label} saved to {tmp_resource_path}")

            try:
                res.download_dir(tmp_resource_path)
            except Exception as e:
                print(f"  Failed to download resource {res.label}: {e}")

        # Reorganize resource files into meaningful categories
        res_root = os.path.join(exp_dir, "resource_files")
        os.makedirs(res_root, exist_ok=True)

        for root, dirs, files in os.walk(tmp_save_path):
            root_lower = root.lower()

            for file in files:
                src_file = os.path.join(root, file)

                if "trends" in root_lower:
                    dest_dir = os.path.join(res_root, "Behavioral_data-TRENDS")
                elif "vft" in root_lower:
                    dest_dir = os.path.join(res_root, "Behavioral_data-VFT")
                elif "exam" in root_lower:
                    dest_dir = os.path.join(res_root, "Examcards")
                else:
                    dest_dir = os.path.join(res_root, "Other_resources")

                os.makedirs(dest_dir, exist_ok=True)
                dest_file = os.path.join(dest_dir, file)

                if not os.path.exists(dest_file):
                    shutil.move(src_file, dest_file)

        shutil.rmtree(tmp_save_path, ignore_errors=True)

        # Download DICOM scans
        for scan in exp.scans.values():
            scan_label = scan.type.lower()

            if scan_label in SCANS_TO_KEEP:
                temp_out = os.path.join(exp_dir, scan_label)
                dicom_res = scan.resources.get("DICOM")

                if dicom_res:
                    os.makedirs(temp_out, exist_ok=True)
                    print(f"  Scan -> {scan_label}")
                    dicom_res.download_dir(temp_out)

                    final_out = os.path.join(scans_root, scan_label)
                    flatten_and_move(temp_out, final_out)

        # Rename downloaded folder by ANCID
        try:
            sample_dicom = None

            for root, _, files in os.walk(scans_root):
                for f in files:
                    if f.endswith(".dcm"):
                        sample_dicom = os.path.join(root, f)
                        break
                if sample_dicom:
                    break

            if sample_dicom:
                dcm_data = dr(sample_dicom, force=True)
                ancid = getattr(dcm_data, 'AccessionNumber', None)

                if ancid and exp.label != ancid:
                    new_path = os.path.join(OUTPUT_DIR, ancid)

                    if not os.path.exists(new_path):
                        print(f"Renaming folder {exp.label} to {ancid}")
                        os.rename(exp_dir, new_path)
                        exp_dir = new_path
                else:
                    print(f"Keeping original name for {exp.label}")

        except Exception as e:
            print(f"ANCID rename failed for {exp.label}: {e}")

        local_experiments.add(exp.label)

# Update tracking file
with open(LOCAL_EXP_FILE, "w") as f:
    for exp_label in sorted(local_experiments):
        f.write(exp_label + "\n")

print("local_experiments.txt updated.")
print("Delta download and folder renaming completed.")

# ============================================================================
# DICOM AND RESOURCE COMPLETENESS CHECKING
# ============================================================================

final_info = []
behavioral_status = {}

for subj in os.listdir(OUTPUT_DIR):

    subj_path = os.path.join(OUTPUT_DIR, subj)

    if not os.path.isdir(subj_path):
        continue

    behavioral_status[subj] = {"VFT": 0, "TRENDS": 0, "EXAMCARD": 0}

    scans_dir = os.path.join(subj_path, "scans")

    if os.path.isdir(scans_dir):
        for scan_label in os.listdir(scans_dir):

            scan_path = os.path.join(scans_dir, scan_label)

            if not os.path.isdir(scan_path):
                continue

            n_files = count_dicom_files(scan_path)

            sample_dicom = None

            if n_files > 0:
                sample_dicom = os.path.join(scan_path, os.listdir(scan_path)[0])

                try:
                    dcm_data = dr(sample_dicom, force=True)

                    date = getattr(dcm_data, 'PerformedProcedureStepStartDate', None)

                    final_info.append([
                        subj,
                        getattr(dcm_data, 'AccessionNumber', None),
                        getattr(dcm_data, 'PatientName', None),
                        getattr(dcm_data, 'StudyDescription', None),
                        getattr(dcm_data, 'StudyComments', None),
                        scan_label,
                        getattr(dcm_data, 'ProtocolName', None),
                        n_files,
                        date
                    ])

                except:
                    final_info.append([subj, None, None, None, None, scan_label, None, n_files, None])

    res_dir = os.path.join(subj_path, "resource_files")

    if os.path.isdir(res_dir):

        for folder in ["Behavioral_data-VFT", "Behavioral_data-TRENDS", "Examcards"]:

            resource_path = os.path.join(res_dir, folder)

            if not os.path.isdir(resource_path):
                continue

            files = os.listdir(resource_path)

            if folder == "Behavioral_data-VFT":
                has_wav = any(f.lower().endswith(".wav") for f in files)
                has_export = any("export.txt" in f.lower() for f in files)

                if has_wav and has_export:
                    behavioral_status[subj]["VFT"] = 1

            elif folder == "Behavioral_data-TRENDS":
                has_export = any("export.txt" in f.lower() for f in files)

                if has_export:
                    behavioral_status[subj]["TRENDS"] = 1

            elif folder == "Examcards":
                behavioral_status[subj]["EXAMCARD"] = 1

# ============================================================================
# CREATE DATAFRAME AND COMPLETENESS REPORTS
# ============================================================================

data = pd.DataFrame(final_info, columns=[
    'SUBJ', 'ANCID', 'Name', 'ADBS_ID', 'Assesment_ID',
    'Seq', 'Sequence_name', 'N_files', 'DATE'
])

data['Sequence_name_clean'] = data['Sequence_name'].str.replace(
    r'^WIP\s+', '', regex=True
).str.lower()

def is_complete(row):
    seq = row['Sequence_name_clean']
    if seq in SCANS_TO_KEEP:
        return 1 if row['N_files'] >= SCANS_TO_KEEP[seq] else 0
    return -1

data['N_files_Complete'] = data.apply(is_complete, axis=1)

modified_data_file = "/path/to/modified_data.csv"

if os.path.exists(modified_data_file):
    existing_data = pd.read_csv(modified_data_file)
    combined_data = pd.concat([existing_data, data], ignore_index=True)
else:
    combined_data = data

combined_data = combined_data.drop_duplicates(
    subset=['SUBJ', 'Seq']
).sort_values(by='ANCID')

combined_data.to_csv(modified_data_file, index=False)

print("modified_data.csv updated and sorted.")

missing_sequences = defaultdict(list)

for subj, group in combined_data.groupby('SUBJ'):

    group_clean = group.set_index('Sequence_name_clean')

    for seq, expected_count in SCANS_TO_KEEP.items():

        if (
            seq not in group_clean.index or
            group_clean.loc[seq, 'N_files'].max() < expected_count
        ):
            missing_sequences[subj].append(seq)

missing_data = []

for subj, sequences in missing_sequences.items():

    subj_data = combined_data[combined_data['SUBJ'] == subj]

    status = behavioral_status.get(subj, {"VFT": 0, "TRENDS": 0, "EXAMCARD": 0})

    row = [
        subj,
        subj_data['ANCID'].iloc[0] if not subj_data.empty else None,
        subj_data['Name'].iloc[0] if not subj_data.empty else None,
        subj_data['ADBS_ID'].iloc[0] if not subj_data.empty else None,
        subj_data['Assesment_ID'].iloc[0] if not subj_data.empty else None,
        subj_data['DATE'].iloc[0] if not subj_data.empty else None,
        ','.join(sorted(set(sequences))),
        status["VFT"],
        status["TRENDS"],
        status["EXAMCARD"]
    ]

    missing_data.append(row)

missing_df = pd.DataFrame(missing_data, columns=[
    'SUBJ', 'ANCID', 'Name', 'ADBS_ID', 'Assesment_ID', 'Date',
    'Missing Sequences', 'VFT', 'TRENDS', 'EXAMCARD'
])

missing_sequences_file = "/path/to/missing_sequences.csv"

if os.path.exists(missing_sequences_file):
    existing_missing = pd.read_csv(missing_sequences_file)
    combined_missing = pd.concat([existing_missing, missing_df], ignore_index=True)
else:
    combined_missing = missing_df

combined_missing = combined_missing.drop_duplicates(
    subset=['SUBJ']
).sort_values(by='ANCID')

combined_missing.to_csv(missing_sequences_file, index=False)

print("missing_sequences.csv updated and sorted.")
print("Processing complete.")
