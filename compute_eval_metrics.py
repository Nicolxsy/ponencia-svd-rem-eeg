import os
import json
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, welch

DATA_DIR = r"C:\Users\nicol\.gemini\antigravity\scratch\ponencia_svd_rem\data"
npz = np.load(os.path.join(DATA_DIR, "sleep_edf_rem_10subjects.npz"))
X_all = npz["X_rem"].astype(np.float64)
df_meta = pd.read_csv(os.path.join(DATA_DIR, "rem_epochs_metadata.csv"))

FS = 100.0
BETA_PZOZ = 0.30
ALPHA_LEVELS = [0.25, 0.50, 0.80]
TAU_OPT = 0.50

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

X_ref = X_all[mask_ref] - np.mean(X_all[mask_ref], axis=1, keepdims=True)
df_ref = df_meta[mask_ref].copy().reset_index(drop=True)

dev_mask = (df_ref["split"] == "Desarrollo").values
eval_mask = (df_ref["split"] == "Evaluacion").values
eval_indices = np.where(eval_mask)[0]

b_butt, a_butt = butter(4, [0.5 / (FS / 2.0), 30.0 / (FS / 2.0)], btype="bandpass")

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

def calc_frob_rel_err(X_eeg_ref, X_eeg_proc):
    num = np.linalg.norm(X_eeg_ref - X_eeg_proc, ord="fro")
    den = np.linalg.norm(X_eeg_ref, ord="fro")
    return float(num / max(den, 1e-15))

def calc_delta_rel_err(ref, proc):
    f, p_ref = welch(ref, fs=FS, nperseg=512)
    _, p_proc = welch(proc, fs=FS, nperseg=512)
    m = (f >= 0.5) & (f <= 4.0)
    er = np.trapezoid(p_ref[m], f[m])
    ep = np.trapezoid(p_proc[m], f[m])
    return float(abs(ep - er) / max(er, 1e-15))

epoch_records = []
svd_spectrum_records = []

for alpha in ALPHA_LEVELS:
    # Build synthetic mixture
    X_cont = np.empty_like(X_ref)
    eog = X_ref[:, :, 2]
    X_cont[:, :, 0] = X_ref[:, :, 0] + alpha * eog
    X_cont[:, :, 1] = X_ref[:, :, 1] + (BETA_PZOZ * alpha) * eog
    X_cont[:, :, 2] = eog
    
    dev_std = np.std(X_cont[dev_mask], axis=(0, 1))
    
    for i in eval_indices:
        subj = df_ref.loc[i, "subject_id"]
        rec_idx = int(df_ref.loc[i, "record_idx"])
        x_ref_eeg = X_ref[i, :, :2] # (3000, 2)
        x_raw_eeg = X_cont[i, :, :2]
        eog_i = X_cont[i, :, 2]
        
        # 1. Butterworth baseline (0.5 - 30 Hz, 4th order, zero phase)
        x_butt_eeg = np.column_stack([
            filtfilt(b_butt, a_butt, x_raw_eeg[:, 0]),
            filtfilt(b_butt, a_butt, x_raw_eeg[:, 1]),
        ])
        
        # 2. SVD with Dev Z-score scaling (Objective 3 baseline protocol)
        X_norm_z = X_cont[i] / dev_std
        U_z, S_z, Vt_z = np.linalg.svd(X_norm_z, full_matrices=False)
        r_z = np.array([abs(np.corrcoef(U_z[:, k], eog_i)[0, 1]) for k in range(3)])
        S_z_clean = S_z.copy()
        S_z_clean[r_z >= TAU_OPT] = 0.0
        X_hat_z = ((U_z * S_z_clean) @ Vt_z) * dev_std
        x_svd_z_eeg = X_hat_z[:, :2]
        
        # 3. SVD in physical units (uV, common scalar scale)
        U_p, S_p, Vt_p = np.linalg.svd(X_cont[i], full_matrices=False)
        r_p = np.array([abs(np.corrcoef(U_p[:, k], eog_i)[0, 1]) for k in range(3)])
        S_p_clean = S_p.copy()
        S_p_clean[r_p >= TAU_OPT] = 0.0
        X_hat_p = (U_p * S_p_clean) @ Vt_p
        x_svd_p_eeg = X_hat_p[:, :2]
        
        # Record SVD spectrum metrics (for both scaling variants)
        for mode_name, S_vec, r_vec, Vt_mat in [
            ("SVD_Zscore", S_z, r_z, Vt_z),
            ("SVD_Fisico_uV", S_p, r_p, Vt_p)
        ]:
            e_vec = (S_vec ** 2) / np.sum(S_vec ** 2) * 100.0
            svd_spectrum_records.append({
                "alpha": alpha,
                "subject_id": subj,
                "epoch_eval_idx": int(i),
                "record_idx": rec_idx,
                "mode": mode_name,
                "sigma_1": float(S_vec[0]),
                "sigma_2": float(S_vec[1]),
                "sigma_3": float(S_vec[2]),
                "energy_1_pct": float(e_vec[0]),
                "energy_2_pct": float(e_vec[1]),
                "energy_3_pct": float(e_vec[2]),
                "energy_top2_pct": float(e_vec[0] + e_vec[1]),
                "corr_u1_eog": float(r_vec[0]),
                "corr_u2_eog": float(r_vec[1]),
                "corr_u3_eog": float(r_vec[2]),
                "v1_fpz": float(abs(Vt_mat[0, 0])),
                "v1_pzoz": float(abs(Vt_mat[0, 1])),
                "v1_eog": float(abs(Vt_mat[0, 2])),
                "n_removed": int(np.sum(r_vec >= TAU_OPT)),
            })
            
        # Record reconstruction metrics for each method
        for method_name, x_proc_eeg in [
            ("Sin_procesar", x_raw_eeg),
            ("Butterworth_0.5_30Hz", x_butt_eeg),
            ("SVD_Zscore", x_svd_z_eeg),
            ("SVD_Fisico_uV", x_svd_p_eeg),
        ]:
            epoch_records.append({
                "alpha": alpha,
                "subject_id": subj,
                "epoch_eval_idx": int(i),
                "record_idx": rec_idx,
                "method": method_name,
                "snr_fpz_db": calc_snr_db(x_ref_eeg[:, 0], x_proc_eeg[:, 0]),
                "snr_pzoz_db": calc_snr_db(x_ref_eeg[:, 1], x_proc_eeg[:, 1]),
                "rmse_fpz_uv": calc_rmse(x_ref_eeg[:, 0], x_proc_eeg[:, 0]),
                "rmse_pzoz_uv": calc_rmse(x_ref_eeg[:, 1], x_proc_eeg[:, 1]),
                "corr_fpz": calc_corr(x_ref_eeg[:, 0], x_proc_eeg[:, 0]),
                "corr_pzoz": calc_corr(x_ref_eeg[:, 1], x_proc_eeg[:, 1]),
                "frob_rel_err_eeg": calc_frob_rel_err(x_ref_eeg, x_proc_eeg),
                "delta_psd_rel_err_fpz": calc_delta_rel_err(x_ref_eeg[:, 0], x_proc_eeg[:, 0]),
            })

df_eval = pd.DataFrame(epoch_records)
df_spec = pd.DataFrame(svd_spectrum_records)

df_eval.to_csv(os.path.join(DATA_DIR, "eval_metrics_all_epochs.csv"), index=False)
df_spec.to_csv(os.path.join(DATA_DIR, "eval_svd_spectrum.csv"), index=False)

print("=== 1. SVD SPECTRUM SUMMARY (Evaluation Set, N=183 epochs) ===")
spec_summary = df_spec.groupby(["mode", "alpha"])[
    ["sigma_1", "sigma_2", "sigma_3", "energy_1_pct", "energy_2_pct", "energy_3_pct", "corr_u1_eog", "corr_u2_eog", "corr_u3_eog", "n_removed"]
].agg(["mean", "std"]).round(2)
print(spec_summary.to_string())

print("\n=== 2. RECONSTRUCTION METRICS SUMMARY (Evaluation Set, N=183 epochs) ===")
eval_summary = df_eval.groupby(["alpha", "method"])[
    ["snr_fpz_db", "snr_pzoz_db", "rmse_fpz_uv", "corr_fpz", "frob_rel_err_eeg", "delta_psd_rel_err_fpz"]
].agg(["mean", "std"]).round(3)
print(eval_summary.to_string())

# 3. Paired Bootstrap Confidence Intervals (95% CI, B=10,000 resamples by subject and by epoch)
rng = np.random.default_rng(2026)
B = 10000
bootstrap_results = []

eval_subjects = sorted(df_eval["subject_id"].unique())

for alpha in ALPHA_LEVELS:
    sub_df = df_eval[df_eval["alpha"] == alpha]
    piv_fpz = sub_df.pivot(index=["subject_id", "epoch_eval_idx"], columns="method", values="snr_fpz_db")
    piv_pzoz = sub_df.pivot(index=["subject_id", "epoch_eval_idx"], columns="method", values="snr_pzoz_db")
    
    for svd_col in ["SVD_Zscore", "SVD_Fisico_uV"]:
        for ref_col in ["Sin_procesar", "Butterworth_0.5_30Hz"]:
            diff_fpz = piv_fpz[svd_col] - piv_fpz[ref_col]
            diff_pzoz = piv_pzoz[svd_col] - piv_pzoz[ref_col]
            
            # Subject-level means for cluster bootstrap by subject
            subj_means_fpz = diff_fpz.groupby("subject_id").mean().values
            subj_means_pzoz = diff_pzoz.groupby("subject_id").mean().values
            
            # Cluster bootstrap over the 5 evaluation subjects
            boot_idx_subj = rng.integers(0, len(eval_subjects), size=(B, len(eval_subjects)))
            boot_subj_fpz = subj_means_fpz[boot_idx_subj].mean(axis=1)
            boot_subj_pzoz = subj_means_pzoz[boot_idx_subj].mean(axis=1)
            
            # Stratified epoch bootstrap (resampling epochs within each subject)
            vals_fpz = diff_fpz.values
            boot_idx_ep = rng.integers(0, len(vals_fpz), size=(B, len(vals_fpz)))
            boot_ep_fpz = vals_fpz[boot_idx_ep].mean(axis=1)
            
            bootstrap_results.append({
                "alpha": alpha,
                "comparison": f"{svd_col} vs {ref_col}",
                "mean_diff_fpz_db": float(diff_fpz.mean()),
                "std_diff_fpz_db": float(diff_fpz.std()),
                "subj_ci95_low_fpz": float(np.percentile(boot_subj_fpz, 2.5)),
                "subj_ci95_high_fpz": float(np.percentile(boot_subj_fpz, 97.5)),
                "epoch_ci95_low_fpz": float(np.percentile(boot_ep_fpz, 2.5)),
                "epoch_ci95_high_fpz": float(np.percentile(boot_ep_fpz, 97.5)),
                "mean_diff_pzoz_db": float(diff_pzoz.mean()),
                "subj_ci95_low_pzoz": float(np.percentile(boot_subj_pzoz, 2.5)),
                "subj_ci95_high_pzoz": float(np.percentile(boot_subj_pzoz, 97.5)),
                "pct_epochs_improved_fpz": float((diff_fpz > 0).mean() * 100.0),
            })

df_boot = pd.DataFrame(bootstrap_results)
df_boot.to_csv(os.path.join(DATA_DIR, "bootstrap_ci95_results.csv"), index=False)
print("\n=== 3. PAIRED BOOTSTRAP 95% CI RESULTS ===")
print(df_boot.round(3).to_string(index=False))
