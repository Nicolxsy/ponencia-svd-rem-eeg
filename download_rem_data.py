import os
import re
import json
import urllib.request
import numpy as np
import pandas as pd

BASE_URL = "https://physionet.org/files/sleep-edfx/1.0.0/sleep-cassette/"
DATA_DIR = r"C:\Users\nicol\.gemini\antigravity\scratch\ponencia_svd_rem\data"
os.makedirs(DATA_DIR, exist_ok=True)

# Fetch directory listing to map PSG to Hypnogram filenames
req = urllib.request.Request(BASE_URL, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req, timeout=30) as resp:
    html = resp.read().decode("utf-8")

psg_files = sorted(set(re.findall(r'SC4\d{3}E0-PSG\.edf', html)))
hyp_files = sorted(set(re.findall(r'SC4\d{3}E[A-Z0-9]-Hypnogram\.edf', html)))

hyp_map = {h[:6]: h for h in hyp_files}
psg_map = {p[:6]: p for p in psg_files}

# Select first 10 distinct subjects (Night 1: SC4001 to SC4091)
subject_Prefixes = [f"SC4{i:02d}1" for i in range(10)]
for sp in subject_Prefixes:
    assert sp in psg_map and sp in hyp_map, f"Missing {sp}"

print("Selected 10 subjects (Night 1):", subject_Prefixes)
print("  Development set (5 subjects):", subject_Prefixes[:5])
print("  Evaluation set  (5 subjects):", subject_Prefixes[5:])

def fetch_bytes(url, start=None, end=None):
    headers = {"User-Agent": "Mozilla/5.0"}
    if start is not None and end is not None:
        headers["Range"] = f"bytes={start}-{end}"
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=30) as resp:
        return resp.read()

def parse_edf_header(hdr_bytes):
    header_nbytes = int(hdr_bytes[184:192].decode("ascii").strip())
    n_records = int(hdr_bytes[236:244].decode("ascii").strip())
    rec_dur = float(hdr_bytes[244:252].decode("ascii").strip())
    ns = int(hdr_bytes[252:256].decode("ascii").strip())
    pos = 256
    labels = [hdr_bytes[pos + i*16 : pos + (i+1)*16].decode("ascii").strip() for i in range(ns)]
    pos += ns * 16 + ns * 80
    phys_dim = [hdr_bytes[pos + i*8 : pos + (i+1)*8].decode("ascii").strip() for i in range(ns)]
    pos += ns * 8
    phys_min = [float(hdr_bytes[pos + i*8 : pos + (i+1)*8].decode("ascii").strip()) for i in range(ns)]
    pos += ns * 8
    phys_max = [float(hdr_bytes[pos + i*8 : pos + (i+1)*8].decode("ascii").strip()) for i in range(ns)]
    pos += ns * 8
    dig_min = [int(hdr_bytes[pos + i*8 : pos + (i+1)*8].decode("ascii").strip()) for i in range(ns)]
    pos += ns * 8
    dig_max = [int(hdr_bytes[pos + i*8 : pos + (i+1)*8].decode("ascii").strip()) for i in range(ns)]
    pos += ns * 8 + ns * 80
    n_samps = [int(hdr_bytes[pos + i*8 : pos + (i+1)*8].decode("ascii").strip()) for i in range(ns)]
    return {
        "header_nbytes": header_nbytes,
        "n_records": n_records,
        "rec_dur": rec_dur,
        "ns": ns,
        "labels": labels,
        "phys_dim": phys_dim,
        "phys_min": phys_min,
        "phys_max": phys_max,
        "dig_min": dig_min,
        "dig_max": dig_max,
        "n_samps": n_samps,
    }

def parse_hypnogram_rem(hyp_bytes):
    text = hyp_bytes[512:].decode("latin-1", errors="ignore")
    pattern = re.compile(r'\+(\d+(?:\.\d+)?)\x15(\d+(?:\.\d+)?)\x14([^\x14\x00]+)')
    rem_intervals = []
    stage_durations = {}
    for m in pattern.finditer(text):
        onset = float(m.group(1))
        dur = float(m.group(2))
        stage = m.group(3).strip()
        stage_durations[stage] = stage_durations.get(stage, 0.0) + dur
        if stage == "Sleep stage R":
            rem_intervals.append((onset, dur))
    return rem_intervals, stage_durations

all_epochs_data = [] # list of arrays (3000, 3)
all_epochs_meta = []
subject_summary = []

for idx, sp in enumerate(subject_Prefixes):
    split = "Desarrollo" if idx < 5 else "Evaluacion"
    psg_file = psg_map[sp]
    hyp_file = hyp_map[sp]
    
    # 1. Download and save Hypnogram
    hyp_bytes = fetch_bytes(BASE_URL + hyp_file)
    with open(os.path.join(DATA_DIR, hyp_file), "wb") as f:
        f.write(hyp_bytes)
    rem_intervals, stage_durations = parse_hypnogram_rem(hyp_bytes)
    
    # 2. Download EDF header
    hdr_raw = fetch_bytes(BASE_URL + psg_file, 0, 4095)
    info = parse_edf_header(hdr_raw)
    
    # Verify channels
    assert info["labels"][0] == "EEG Fpz-Cz"
    assert info["labels"][1] == "EEG Pz-Oz"
    assert info["labels"][2] == "EOG horizontal"
    assert info["rec_dur"] == 30.0
    assert info["n_samps"][0] == 3000 and info["n_samps"][1] == 3000 and info["n_samps"][2] == 3000
    
    bytes_per_record = sum(info["n_samps"]) * 2
    gains = [(info["phys_max"][c] - info["phys_min"][c]) / (info["dig_max"][c] - info["dig_min"][c]) for c in range(3)]
    offsets = [info["phys_min"][c] - gains[c] * info["dig_min"][c] for c in range(3)]
    
    n_rem_epochs_subj = 0
    for (onset, dur) in rem_intervals:
        start_rec = int(round(onset / 30.0))
        n_recs = int(round(dur / 30.0))
        if start_rec + n_recs > info["n_records"]:
            n_recs = max(0, info["n_records"] - start_rec)
        if n_recs <= 0:
            continue
        
        byte_start = info["header_nbytes"] + start_rec * bytes_per_record
        byte_end = byte_start + n_recs * bytes_per_record - 1
        chunk = fetch_bytes(BASE_URL + psg_file, byte_start, byte_end)
        raw_i16 = np.frombuffer(chunk, dtype="<i2").reshape(n_recs, sum(info["n_samps"]))
        
        for r in range(n_recs):
            rec_idx = start_rec + r
            fpz_cz = raw_i16[r, 0:3000].astype(np.float64) * gains[0] + offsets[0]
            pz_oz  = raw_i16[r, 3000:6000].astype(np.float64) * gains[1] + offsets[1]
            eog_h  = raw_i16[r, 6000:9000].astype(np.float64) * gains[2] + offsets[2]
            mat = np.column_stack([fpz_cz, pz_oz, eog_h]) # (3000, 3)
            all_epochs_data.append(mat.astype(np.float32))
            all_epochs_meta.append({
                "subject_id": sp,
                "psg_file": psg_file,
                "hyp_file": hyp_file,
                "split": split,
                "record_idx": rec_idx,
                "onset_sec": rec_idx * 30.0,
                "std_fpz": float(np.std(fpz_cz)),
                "std_pzoz": float(np.std(pz_oz)),
                "std_eog": float(np.std(eog_h)),
                "ptp_eog": float(np.ptp(eog_h)),
                "corr_fpz_eog": float(np.corrcoef(fpz_cz, eog_h)[0, 1]),
                "corr_pzoz_eog": float(np.corrcoef(pz_oz, eog_h)[0, 1]),
            })
            n_rem_epochs_subj += 1
            
    subject_summary.append({
        "subject_id": sp,
        "psg_file": psg_file,
        "hyp_file": hyp_file,
        "split": split,
        "total_records_30s": info["n_records"],
        "rem_bouts": len(rem_intervals),
        "rem_epochs_30s": n_rem_epochs_subj,
        "rem_duration_min": n_rem_epochs_subj * 0.5,
        "fs_hz": 100.0,
        "units": info["phys_dim"][0],
    })
    print(f"  [{split}] {sp} ({psg_file}): {n_rem_epochs_subj} REM epochs ({n_rem_epochs_subj*0.5:.1f} min) downloaded.")

X_all = np.stack(all_epochs_data, axis=0) # (M, 3000, 3)
df_meta = pd.DataFrame(all_epochs_meta)
df_subj = pd.DataFrame(subject_summary)

np.savez_compressed(os.path.join(DATA_DIR, "sleep_edf_rem_10subjects.npz"), X_rem=X_all)
df_meta.to_csv(os.path.join(DATA_DIR, "rem_epochs_metadata.csv"), index=False)
df_subj.to_csv(os.path.join(DATA_DIR, "subjects_summary.csv"), index=False)

print(f"\nSaved {X_all.shape[0]} real REM epochs of shape {X_all.shape[1:]} to sleep_edf_rem_10subjects.npz")
print(df_subj.to_string(index=False))
