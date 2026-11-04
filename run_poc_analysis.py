import os
import json
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, welch

DATA_DIR = r"C:\Users\nicol\.gemini\antigravity\scratch\ponencia_svd_rem\data"
FIG_DIR  = r"C:\Users\nicol\.gemini\antigravity\scratch\ponencia_svd_rem\figures"
os.makedirs(FIG_DIR, exist_ok=True)

# 1. Load real REM epochs from 10 Sleep-EDF subjects
npz = np.load(os.path.join(DATA_DIR, "sleep_edf_rem_10subjects.npz"))
X_all = npz["X_rem"].astype(np.float64) # (1594, 3000, 3) -> [Fpz-Cz, Pz-Oz, EOG]
df_meta = pd.read_csv(os.path.join(DATA_DIR, "rem_epochs_metadata.csv"))
df_subj = pd.read_csv(os.path.join(DATA_DIR, "subjects_summary.csv"))

FS = 100.0
N_SAMPLES = 3000

# 2. Criterion for selecting clean REM reference segments (x_ref) + active EOG
# Two-stage inspection criterion:
#   (a) Low baseline cross-correlation between EEG and EOG: |r(Fpz-Cz, EOG)| < 0.15 and |r(Pz-Oz, EOG)| < 0.15
#   (b) Absence of high-amplitude saturation/movement spikes in EEG: std(EEG) < 25 uV and max|EEG| < 120 uV
#   (c) Presence of genuine REM ocular activity in the EOG channel: std(EOG) > 15 uV
max_abs_fpz = np.max(np.abs(X_all[:, :, 0] - np.mean(X_all[:, :, 0], axis=1, keepdims=True)), axis=1)
max_abs_pzoz = np.max(np.abs(X_all[:, :, 1] - np.mean(X_all[:, :, 1], axis=1, keepdims=True)), axis=1)

mask_ref = (
    (df_meta["corr_fpz_eog"].abs() < 0.15) &
    (df_meta["corr_pzoz_eog"].abs() < 0.15) &
    (df_meta["std_fpz"] < 25.0) &
    (df_meta["std_pzoz"] < 25.0) &
    (max_abs_fpz < 120.0) &
    (max_abs_pzoz < 120.0) &
    (df_meta["std_eog"] > 15.0)
)

X_ref_raw = X_all[mask_ref]
df_ref = df_meta[mask_ref].copy().reset_index(drop=True)

# Zero-mean centering per epoch (removes DC offset)
X_ref = X_ref_raw - np.mean(X_ref_raw, axis=1, keepdims=True)

dev_mask = (df_ref["split"] == "Desarrollo").values
eval_mask = (df_ref["split"] == "Evaluacion").values

print(f"Total REM epochs analyzed: {len(df_meta)} (Dev: {(df_meta['split']=='Desarrollo').sum()}, Eval: {(df_meta['split']=='Evaluacion').sum()})")
print(f"Selected REM reference epochs (x_ref): {len(df_ref)} (Dev: {dev_mask.sum()}, Eval: {eval_mask.sum()})")
print("\nBreakdown by subject:")
subj_counts = df_ref.groupby(["split", "subject_id"]).size().reset_index(name="n_ref_epochs")
df_subj_merged = pd.merge(df_subj, subj_counts, on=["split", "subject_id"], how="left")
print(df_subj_merged.to_string(index=False))
df_subj_merged.to_csv(os.path.join(DATA_DIR, "subjects_table_complete.csv"), index=False)

# 3. Reference Butterworth bandpass filter: 4th-order, 0.5 - 30 Hz, zero-phase (filtfilt)
b_butt, a_butt = butter(4, [0.5 / (FS / 2.0), 30.0 / (FS / 2.0)], btype="bandpass")

# Helper metrics
def calc_snr_db(ref, proc):
    p_sig = np.sum(ref ** 2)
    p_err = np.sum((proc - ref) ** 2)
    return float(10.0 * np.log10(p_sig / max(p_err, 1e-15)))

def calc_rmse(ref, proc):
    return float(np.sqrt(np.mean((proc - ref) ** 2)))

def calc_corr(ref, proc):
    if np.std(proc) < 1e-12:
        return 0.0
    return float(np.corrcoef(ref, proc)[0, 1])

def calc_delta_psd_rel_error(ref, proc, fs=100.0):
    f, p_ref = welch(ref, fs=fs, nperseg=512)
    _, p_proc = welch(proc, fs=fs, nperseg=512)
    idx_delta = (f >= 0.5) & (f <= 4.0)
    e_ref = np.trapezoid(p_ref[idx_delta], f[idx_delta])
    e_proc = np.trapezoid(p_proc[idx_delta], f[idx_delta])
    return float(abs(e_proc - e_ref) / max(e_ref, 1e-15))

# 4. Synthetic mixture protocol:
# Following Waterman et al. (1992), ocular propagation is maximal at frontal leads (Fpz-Cz: factor alpha)
# and attenuated at parietal-occipital leads (Pz-Oz: factor beta * alpha, with beta = 0.30).
# We also evaluate the symmetric case (beta = 1.0) in sensitivity analysis.
BETA_PZOZ = 0.30
ALPHA_LEVELS = [0.25, 0.50, 0.80] # Low, Medium (primary benchmark), High

def build_synthetic_bank(X_ref_arr, alpha, beta=BETA_PZOZ):
    """
    Constructs X_cont of shape (M, 3000, 3):
      col 0: Fpz-Cz_cont = Fpz-Cz_ref + alpha * EOG
      col 1: Pz-Oz_cont  = Pz-Oz_ref  + (beta * alpha) * EOG
      col 2: EOG
    """
    X_cont = np.empty_like(X_ref_arr)
    eog = X_ref_arr[:, :, 2]
    X_cont[:, :, 0] = X_ref_arr[:, :, 0] + alpha * eog
    X_cont[:, :, 1] = X_ref_arr[:, :, 1] + (beta * alpha) * eog
    X_cont[:, :, 2] = eog
    return X_cont

# 5. Calibrate correlation threshold tau on the Development set (N_dev = 285 epochs across 5 subjects)
# We test tau in [0.15, 0.20, ..., 0.90] for both scaling modes:
#   - Mode A ("dev_zscore"): standardized by development-set channel standard deviations (as written in Obj. 3)
#   - Mode B ("physical_uV"): centered in physical units (uV), preserving relative channel energy
thresholds = np.round(np.arange(0.15, 0.91, 0.05), 2)
sensitivity_rows = []

for scale_mode in ["dev_zscore", "physical_uV"]:
    for alpha in ALPHA_LEVELS:
        X_cont = build_synthetic_bank(X_ref, alpha, BETA_PZOZ)
        dev_std = np.std(X_cont[dev_mask], axis=(0, 1)) if scale_mode == "dev_zscore" else np.ones(3)
        
        # Precompute SVD for all dev epochs
        U_list, S_list, Vt_list, r_eog_list = [], [], [], []
        for i in np.where(dev_mask)[0]:
            X_ep = X_cont[i] / dev_std
            U, S, Vt = np.linalg.svd(X_ep, full_matrices=False)
            eog_ep = X_cont[i, :, 2]
            r_k = np.array([abs(np.corrcoef(U[:, k], eog_ep)[0, 1]) for k in range(3)])
            U_list.append(U)
            S_list.append(S)
            Vt_list.append(Vt)
            r_eog_list.append(r_k)
            
        for tau in thresholds:
            snrs_fpz = []
            snrs_pzoz = []
            n_Zeroed = []
            for idx_local, i_global in enumerate(np.where(dev_mask)[0]):
                U = U_list[idx_local]
                S = S_list[idx_local].copy()
                Vt = Vt_list[idx_local]
                r_k = r_eog_list[idx_local]
                
                kill = (r_k >= tau)
                S[kill] = 0.0
                n_Zeroed.append(int(np.sum(kill)))
                X_hat = ((U * S) @ Vt) * dev_std
                snrs_fpz.append(calc_snr_db(X_ref[i_global, :, 0], X_hat[:, 0]))
                snrs_pzoz.append(calc_snr_db(X_ref[i_global, :, 1], X_hat[:, 1]))
                
            sensitivity_rows.append({
                "scale_mode": scale_mode,
                "alpha": alpha,
                "tau": float(tau),
                "mean_snr_fpz_db": float(np.mean(snrs_fpz)),
                "std_snr_fpz_db": float(np.std(snrs_fpz)),
                "mean_snr_pzoz_db": float(np.mean(snrs_pzoz)),
                "mean_components_removed": float(np.mean(n_Zeroed)),
            })

df_sens = pd.DataFrame(sensitivity_rows)
df_sens.to_csv(os.path.join(DATA_DIR, "threshold_sensitivity_dev.csv"), index=False)

# Find optimal tau on Development set (averaged across alpha levels)
for sm in ["dev_zscore", "physical_uV"]:
    sub = df_sens[df_sens["scale_mode"] == sm].groupby("tau")["mean_snr_fpz_db"].mean()
    best_tau = sub.idxmax()
    print(f"\n[{sm}] Optimal threshold on Development set: tau* = {best_tau:.2f} (mean Dev SNR Fpz-Cz = {sub.max():.2f} dB)")
    print(df_sens[(df_sens["scale_mode"] == sm) & (df_sens["alpha"] == 0.50)][["tau", "mean_snr_fpz_db", "mean_snr_pzoz_db", "mean_components_removed"]].to_string(index=False))
