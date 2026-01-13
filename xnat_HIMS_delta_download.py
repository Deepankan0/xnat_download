import xnat
import os
import urllib.parse
import shutil
from pydicom import dcmread as dr
import pandas as pd
from collections import defaultdict

# ----------------- CONFIG -----------------
PROJECT_ID = "PROJECT_ID"
XNAT_URL = "URL_TO_THE_XNAT_SERVER"
XNAT_USER = "XNAT_USER"
XNAT_PASS = "YOUR_PASSWORD"

OUTPUT_DIR = "PATH/TO/YOUR/OUTPUT_DIR"
LOCAL_EXP_FILE = "local_experiments.txt"

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

# ----------------- HELPER FUNCTIONS -----------------
def flatten_and_move(src_dir, dest_dir):
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
    if not os.path.isdir(scan_path):
        return 0
    return len([f for f in os.listdir(scan_path) if f.endswith(".dcm")])

def parse_study_comments(study_comments):
    """
    Expected format:
    <Assessment_ID>-<ANCID>
    Example: 111556677-ANC999
    """
    if not study_comments or '-' not in str(study_comments):
        return None, None, None
    try:
        assessment_id, ancid = str(study_comments).split('-', 1)
        assessment_id = assessment_id.strip()
        ancid = ancid.strip()
        if not assessment_id.isdigit():
            return None, None, None
        adbs_id = assessment_id[-6:]
        return assessment_id, ancid, adbs_id
    except Exception:
        return None, None, None

# ----------------- LOAD LOCAL EXPERIMENTS -----------------
if os.path.exists(LOCAL_EXP_FILE):
    with open(LOCAL_EXP_FILE, "r") as f:
        local_experiments = set(line.strip() for line in f.readlines())
else:
    local_experiments = set()

print(f"Found {len(local_experiments)} experiments locally (from txt).")

# ----------------- XNAT DOWNLOAD -----------------
with xnat.connect(XNAT_URL, user=XNAT_USER, password=XNAT_PASS) as session:
    project = session.projects[PROJECT_ID]

    xnat_experiments = {
        exp.label: exp
        for subj in project.subjects.values()
        for exp in subj.experiments.values()
    }

    print(f"\nFound {len(xnat_experiments)} experiments on XNAT.")

    missing_experiments = [
        label for label in xnat_experiments
        if label not in local_experiments
    ]

    print(f"\nMissing experiments (to download): {len(missing_experiments)}")
    print(missing_experiments)

    for exp_label in missing_experiments:
        exp = xnat_experiments[exp_label]
        exp_dir = os.path.join(OUTPUT_DIR, exp.label)
        os.makedirs(exp_dir, exist_ok=True)

        print(f"\n--- Downloading {exp.label} ---")

        scans_root = os.path.join(exp_dir, "scans")
        os.makedirs(scans_root, exist_ok=True)

        # ----------------- RESOURCES -----------------
        tmp_save_path = os.path.join(exp_dir, "tmp_download")
        os.makedirs(tmp_save_path, exist_ok=True)

        for res in exp.resources.values():
            res_label_raw = urllib.parse.unquote(res.label)
            tmp_resource_path = os.path.join(tmp_save_path, res_label_raw)
            os.makedirs(tmp_resource_path, exist_ok=True)
            try:
                res.download_dir(tmp_resource_path)
            except Exception as e:
                print(f"Failed to download {res.label}: {e}")

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

        # ----------------- SCANS -----------------
        for scan in exp.scans.values():
            scan_label = scan.type.lower()
            if scan_label in SCANS_TO_KEEP:
                temp_out = os.path.join(exp_dir, scan_label)
                dicom_res = scan.resources.get("DICOM")
                if dicom_res:
                    os.makedirs(temp_out, exist_ok=True)
                    dicom_res.download_dir(temp_out)
                    final_out = os.path.join(scans_root, scan_label)
                    flatten_and_move(temp_out, final_out)

        # ----------------- RENAME BY ANCID -----------------
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
                study_comments = getattr(dcm_data, 'StudyComments', None)
                assessment_id, ancid, adbs_id = parse_study_comments(study_comments)

                if ancid and exp.label != ancid:
                    new_path = os.path.join(OUTPUT_DIR, ancid)
                    if not os.path.exists(new_path):
                        print(f"Renaming {exp.label} -> {ancid}")
                        os.rename(exp_dir, new_path)
                        exp_dir = new_path
        except Exception as e:
            print(f"Failed to rename {exp.label}: {e}")

        local_experiments.add(exp.label)

# ----------------- UPDATE LOCAL EXP FILE -----------------
with open(LOCAL_EXP_FILE, "w") as f:
    for exp_label in sorted(local_experiments):
        f.write(exp_label + "\n")

print("\n--- Delta download complete ---")

# ----------------- DICOM & BEHAVIORAL CHECK -----------------
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
            if n_files == 0:
                continue

            sample_dicom = os.path.join(scan_path, os.listdir(scan_path)[0])
            try:
                dcm_data = dr(sample_dicom, force=True)
                study_comments = getattr(dcm_data, 'StudyComments', None)
                assessment_id, ancid, adbs_id = parse_study_comments(study_comments)
                date = getattr(dcm_data, 'PerformedProcedureStepStartDate', None)

                final_info.append([
                    subj,
                    ancid,
                    getattr(dcm_data, 'PatientName', None),
                    adbs_id,
                    assessment_id,
                    scan_label,
                    getattr(dcm_data, 'ProtocolName', None),
                    n_files,
                    date
                ])
            except Exception:
                final_info.append([subj, None, None, None, None, scan_label, None, n_files, None])

    res_dir = os.path.join(subj_path, "resource_files")
    if os.path.isdir(res_dir):
        for folder in ["Behavioral_data-VFT", "Behavioral_data-TRENDS", "Examcards"]:
            resource_path = os.path.join(res_dir, folder)
            if not os.path.isdir(resource_path):
                continue
            files = os.listdir(resource_path)
            if folder == "Behavioral_data-VFT":
                if any(f.lower().endswith(".wav") for f in files) and any("export.txt" in f.lower() for f in files):
                    behavioral_status[subj]["VFT"] = 1
            elif folder == "Behavioral_data-TRENDS":
                if any("export.txt" in f.lower() for f in files):
                    behavioral_status[subj]["TRENDS"] = 1
            elif folder == "Examcards":
                behavioral_status[subj]["EXAMCARD"] = 1

# ----------------- CREATE DATAFRAME -----------------
data = pd.DataFrame(final_info, columns=[
    'SUBJ', 'ANCID', 'Name', 'ADBS_ID', 'Assesment_ID',
    'Seq', 'Sequence_name', 'N_files', 'DATE'
])

data['Sequence_name_clean'] = (
    data['Sequence_name']
    .str.replace(r'^WIP\s+', '', regex=True)
    .str.lower()
)

def is_complete(row):
    seq = row['Sequence_name_clean']
    if seq in SCANS_TO_KEEP:
        return 1 if row['N_files'] >= SCANS_TO_KEEP[seq] else 0
    return -1

data['N_files_Complete'] = data.apply(is_complete, axis=1)

# ----------------- SAVE modified_data.csv -----------------
modified_data_file = "xnat/DICOM_information/modified_data.csv"
if os.path.exists(modified_data_file):
    existing_data = pd.read_csv(modified_data_file)
    combined_data = pd.concat([existing_data, data], ignore_index=True)
else:
    combined_data = data

combined_data = combined_data.drop_duplicates(subset=['SUBJ', 'Seq']).sort_values(by='ANCID')
combined_data.to_csv(modified_data_file, index=False)

# ----------------- CHECK MISSING SEQUENCES -----------------
missing_sequences = defaultdict(list)

for subj, group in combined_data.groupby('SUBJ'):
    group_clean = group.set_index('Sequence_name_clean')
    for seq, expected in SCANS_TO_KEEP.items():
        if seq not in group_clean.index or group_clean.loc[seq, 'N_files'].max() < expected:
            missing_sequences[subj].append(seq)
    if subj not in missing_sequences:
        missing_sequences[subj].append('NIL')

missing_data = []

for subj, sequences in missing_sequences.items():
    subj_data = combined_data[combined_data['SUBJ'] == subj]
    status = behavioral_status.get(subj, {"VFT": 0, "TRENDS": 0, "EXAMCARD": 0})

    missing_data.append([
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
    ])

missing_df = pd.DataFrame(missing_data, columns=[
    'SUBJ', 'ANCID', 'Name', 'ADBS_ID', 'Assesment_ID', 'Date',
    'Missing Sequences', 'VFT', 'TRENDS', 'EXAMCARD'
])

missing_sequences_file = "xnat/DICOM_information/missing_sequences.csv"
if os.path.exists(missing_sequences_file):
    existing_missing = pd.read_csv(missing_sequences_file)
    combined_missing = pd.concat([existing_missing, missing_df], ignore_index=True)
else:
    combined_missing = missing_df

combined_missing = combined_missing.drop_duplicates(subset=['SUBJ']).sort_values(by='ANCID')
combined_missing.to_csv(missing_sequences_file, index=False)

print("Modified data and missing sequences saved.")
