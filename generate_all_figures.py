import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.signal import butter, filtfilt, welch

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Segoe UI"],
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

DATA_DIR = r"C:\Users\nicol\.gemini\antigravity\scratch\ponencia_svd_rem\data"
FIG_DIR  = r"C:\Users\nicol\.gemini\antigravity\scratch\ponencia_svd_rem\figures"
os.makedirs(FIG_DIR, exist_ok=True)

npz = np.load(os.path.join(DATA_DIR, "sleep_edf_rem_10subjects.npz"))
X_all = npz["X_rem"].astype(np.float64)
df_meta = pd.read_csv(os.path.join(DATA_DIR, "rem_epochs_metadata.csv"))
df_sens = pd.read_csv(os.path.join(DATA_DIR, "threshold_sensitivity_dev.csv"))
df_eval = pd.read_csv(os.path.join(DATA_DIR, "eval_metrics_all_epochs.csv"))
df_spec = pd.read_csv(os.path.join(DATA_DIR, "eval_svd_spectrum.csv"))
df_boot = pd.read_csv(os.path.join(DATA_DIR, "bootstrap_ci95_results.csv"))

FS = 100.0
BETA_PZOZ = 0.30
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

b_butt, a_butt = butter(4, [0.5 / (FS / 2.0), 30.0 / (FS / 2.0)], btype="bandpass")

def save_fig(fig, name):
    fig.savefig(os.path.join(FIG_DIR, f"{name}.png"), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join(FIG_DIR, f"{name}.svg"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {name}.png and {name}.svg")

# ==============================================================================
# FIGURA 1: Diagrama de Flujo Metodológico de la Prueba de Concepto
# ==============================================================================
fig, ax = plt.subplots(figsize=(11.5, 4.3))
ax.set_xlim(0, 11.4)
ax.set_ylim(0, 4.3)
ax.axis("off")

boxes = [
    (0.2, 2.35, 1.8, 1.25, "#E8EEF5", "#1B365D", "1. Datos crudos\nSleep-EDF v1.0.0\n(10 sujetos, 100 Hz)\n1594 épocas REM"),
    (2.4, 2.35, 1.95, 1.25, "#E8EEF5", "#1B365D", "2. Preprocesamiento\nSelección x_ref (N=465)\nMezcla: x_ref + α·EOG\nCentrado y escalado"),
    (4.75, 2.35, 1.8, 1.25, "#E8EEF5", "#1B365D", "3. Matriz de datos X\nX ∈ ℝ^(3000×3)\n[Fpz-Cz, Pz-Oz, EOG]\nÉpocas de 30 s"),
    (6.95, 2.35, 1.95, 1.25, "#FFF4E5", "#B85D00", "4. Descomposición\nSVD: X = U Σ Vᵀ\nCorrelación |r(u_i, EOG)|\nTruncamiento si |r| ≥ τ*"),
    (9.3, 2.35, 1.9, 1.25, "#E6F4EA", "#137333", "5. Métricas y\nVisualización\nSNR (dB), Error Frob.,\nPSD Delta y Bootstrap"),
    (5.5, 0.35, 3.6, 1.05, "#F3E8FD", "#681DA8", "6. Interpretación en el dominio\nDesacoplamiento EOG en banda delta\ny preservación de morfología EEG"),
]

for (x, y, w, h, fc, ec, txt) in boxes:
    rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08", facecolor=fc, edgecolor=ec, linewidth=1.6)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, txt, ha="center", va="center", fontsize=8.8, color="#111111")

arrow_kw = dict(arrowstyle="->,head_width=0.25,head_length=0.25", color="#1B365D", lw=1.6)
ax.annotate("", xy=(2.4, 2.97), xytext=(2.0, 2.97), arrowprops=arrow_kw)
ax.annotate("", xy=(4.75, 2.97), xytext=(4.35, 2.97), arrowprops=arrow_kw)
ax.annotate("", xy=(6.95, 2.97), xytext=(6.55, 2.97), arrowprops=arrow_kw)
ax.annotate("", xy=(9.3, 2.97), xytext=(8.9, 2.97), arrowprops=arrow_kw)

ax.annotate("", xy=(9.1, 0.88), xytext=(10.25, 2.35),
            arrowprops=dict(arrowstyle="->,head_width=0.25,head_length=0.25", color="#137333", lw=1.6, connectionstyle="angle,angleA=-90,angleB=0,rad=5"))

ax.annotate("", xy=(3.37, 2.35), xytext=(5.5, 0.88),
            arrowprops=dict(arrowstyle="->,head_width=0.25,head_length=0.25", color="#681DA8", lw=1.5, ls="--", connectionstyle="angle,angleA=180,angleB=90,rad=5"))
ax.text(3.37, 0.50, "Calibración de τ* y régimen de escalado\n(Conjunto de Desarrollo, n = 5 sujetos)", ha="center", va="center", fontsize=8.4, style="italic", color="#4A148C")

ax.set_title("Figura 1: Flujo metodológico de la prueba de concepto para reducción de artefactos EOG mediante SVD", pad=10, fontsize=11.5)
save_fig(fig, "fig1_flujo_metodologico")

# ==============================================================================
# FIGURA 2: Espectro de Valores Singulares, Energía Capturada y Correlación con EOG
# ==============================================================================
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))

sub_spec_p = df_spec[(df_spec["mode"] == "SVD_Fisico_uV") & (df_spec["alpha"] == 0.50)]
sub_spec_z = df_spec[(df_spec["mode"] == "SVD_Zscore") & (df_spec["alpha"] == 0.50)]

modes_x = np.arange(1, 4)
width = 0.34
means_p = [sub_spec_p["energy_1_pct"].mean(), sub_spec_p["energy_2_pct"].mean(), sub_spec_p["energy_3_pct"].mean()]
stds_p  = [sub_spec_p["energy_1_pct"].std(),  sub_spec_p["energy_2_pct"].std(),  sub_spec_p["energy_3_pct"].std()]
means_z = [sub_spec_z["energy_1_pct"].mean(), sub_spec_z["energy_2_pct"].mean(), sub_spec_z["energy_3_pct"].mean()]
stds_z  = [sub_spec_z["energy_1_pct"].std(),  sub_spec_z["energy_2_pct"].std(),  sub_spec_z["energy_3_pct"].std()]

b1 = axes[0].bar(modes_x - width/2, means_p, width, yerr=stds_p, capsize=4, label="SVD Escala física (μV)", color="#1B365D", alpha=0.9)
b2 = axes[0].bar(modes_x + width/2, means_z, width, yerr=stds_z, capsize=4, label="SVD Escalado Z-score (Desarrollo)", color="#E07A5F", alpha=0.9)
axes[0].set_xticks(modes_x)
axes[0].set_xticklabels(["Modo 1 ($\\sigma_1$)\n[Subespacio EOG]", "Modo 2 ($\\sigma_2$)\n[EEG dominante]", "Modo 3 ($\\sigma_3$)\n[EEG secundario]"])
axes[0].set_ylabel("Energía capturada $\\sigma_i^2 / \\sum \\sigma_j^2$ (%)")
axes[0].set_ylim(0, 116)
axes[0].set_title("A. Distribución de energía por modo singular ($\\alpha = 0.50$)")
axes[0].legend(loc="upper right", frameon=True)
axes[0].grid(axis="y", linestyle=":", alpha=0.6)

for bar, s in zip(b1, stds_p):
    h = bar.get_height()
    axes[0].text(bar.get_x() + bar.get_width()/2, h + s + 1.8, f"{h:.1f}%", ha="center", va="bottom", fontsize=8.3, fontweight="bold", color="#1B365D")
for bar, s in zip(b2, stds_z):
    h = bar.get_height()
    axes[0].text(bar.get_x() + bar.get_width()/2, h + s + 1.8, f"{h:.1f}%", ha="center", va="bottom", fontsize=8.3, fontweight="bold", color="#B84A27")

rm_p = [sub_spec_p["corr_u1_eog"].mean(), sub_spec_p["corr_u2_eog"].mean(), sub_spec_p["corr_u3_eog"].mean()]
rs_p = [sub_spec_p["corr_u1_eog"].std(),  sub_spec_p["corr_u2_eog"].std(),  sub_spec_p["corr_u3_eog"].std()]
rm_z = [sub_spec_z["corr_u1_eog"].mean(), sub_spec_z["corr_u2_eog"].mean(), sub_spec_z["corr_u3_eog"].mean()]
rs_z = [sub_spec_z["corr_u1_eog"].std(),  sub_spec_z["corr_u2_eog"].std(),  sub_spec_z["corr_u3_eog"].std()]

b3 = axes[1].bar(modes_x - width/2, rm_p, width, yerr=rs_p, capsize=4, label="SVD Escala física (μV)", color="#1B365D", alpha=0.9)
b4 = axes[1].bar(modes_x + width/2, rm_z, width, yerr=rs_z, capsize=4, label="SVD Escalado Z-score (Desarrollo)", color="#E07A5F", alpha=0.9)
axes[1].axhline(TAU_OPT, color="#D93025", linestyle="--", linewidth=1.5, label=f"Umbral calibrado $\\tau^* = {TAU_OPT:.2f}$")
axes[1].set_xticks(modes_x)
axes[1].set_xticklabels(["Vector $u_1(t)$", "Vector $u_2(t)$", "Vector $u_3(t)$"])
axes[1].set_ylabel("Correlación absoluta con canal EOG $|r(u_i, \\mathrm{EOG})|$ (adim.)")
axes[1].set_ylim(0, 1.28)
axes[1].set_title("B. Acoplamiento temporal con el canal EOG ($n = 183$ épocas)")
axes[1].legend(loc="upper right", frameon=True)
axes[1].grid(axis="y", linestyle=":", alpha=0.6)

for bar, s in zip(b3, rs_p):
    h = bar.get_height()
    axes[1].text(bar.get_x() + bar.get_width()/2, h + s + 0.025, f"{h:.2f}", ha="center", va="bottom", fontsize=8.3, fontweight="bold", color="#1B365D")
for bar, s in zip(b4, rs_z):
    h = bar.get_height()
    axes[1].text(bar.get_x() + bar.get_width()/2, h + s + 0.025, f"{h:.2f}", ha="center", va="bottom", fontsize=8.3, fontweight="bold", color="#B84A27")

fig.suptitle("Figura 2: Espectro de energía SVD y correlación temporal con el EOG en el conjunto de Evaluación (n = 183 épocas, 5 sujetos)", y=1.03, fontsize=11.5, fontweight="bold")
fig.tight_layout()
save_fig(fig, "fig2_espectro_svd_y_correlacion")

# ==============================================================================
# FIGURA 3: Curva de Sensibilidad al Umbral tau en Conjunto de Desarrollo (Obj. 4)
# ==============================================================================
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

colors_alpha = {0.25: "#2A9D8F", 0.50: "#1B365D", 0.80: "#D93025"}
labels_alpha = {0.25: "Contaminación baja ($\\alpha = 0.25$)", 0.50: "Contaminación media ($\\alpha = 0.50$)", 0.80: "Contaminación alta ($\\alpha = 0.80$)"}

for alpha in [0.25, 0.50, 0.80]:
    s_z = df_sens[(df_sens["scale_mode"] == "dev_zscore") & (df_sens["alpha"] == alpha)]
    axes[0].plot(s_z["tau"], s_z["mean_snr_fpz_db"], marker="o", ms=4.5, lw=1.8, color=colors_alpha[alpha], label=labels_alpha[alpha])
    
    s_p = df_sens[(df_sens["scale_mode"] == "physical_uV") & (df_sens["alpha"] == alpha)]
    axes[1].plot(s_p["tau"], s_p["mean_snr_fpz_db"], marker="s", ms=4.5, lw=1.8, color=colors_alpha[alpha], label=labels_alpha[alpha])

axes[0].axvline(0.50, color="#111111", linestyle="--", lw=1.4, label="Óptimo calibrado $\\tau^* = 0.50$")
axes[0].set_xlabel("Umbral de correlación $\\tau$ (adimensional)")
axes[0].set_ylabel("SNR media en EEG Fpz-Cz (dB)")
axes[0].set_title("A. Sensibilidad con Escalado Z-score por canal (Desarrollo)")
axes[0].grid(True, linestyle=":", alpha=0.6)
axes[0].legend(loc="lower center", fontsize=8.2, frameon=True)

axes[1].axvline(0.50, color="#111111", linestyle="--", lw=1.4, label="Punto operativo $\\tau^* = 0.50$")
axes[1].axvspan(0.30, 0.85, color="#E6F4EA", alpha=0.5, label="Meseta de estabilidad ($1.00$ comp. anulada)")
axes[1].set_xlabel("Umbral de correlación $\\tau$ (adimensional)")
axes[1].set_ylabel("SNR media en EEG Fpz-Cz (dB)")
axes[1].set_ylim(5.0, 21.0)
axes[1].set_title("B. Sensibilidad en Escala Física común en $\\mu$V (Desarrollo)")
axes[1].grid(True, linestyle=":", alpha=0.6)
axes[1].legend(loc="center right", bbox_to_anchor=(0.98, 0.58), fontsize=8.0, frameon=True)

fig.suptitle("Figura 3: Análisis de sensibilidad del umbral de correlación τ en el conjunto de Desarrollo (n = 282 épocas, 5 sujetos)", y=1.03, fontsize=11.5, fontweight="bold")
fig.tight_layout()
save_fig(fig, "fig3_sensibilidad_umbral_desarrollo")

# ==============================================================================
# FIGURA 4: Desempeño Cuantitativo SNR e Intervalos de Confianza Bootstrap 95%
# ==============================================================================
fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.9))

methods_order = ["Sin_procesar", "Butterworth_0.5_30Hz", "SVD_Zscore", "SVD_Fisico_uV"]
method_labels = {
    "Sin_procesar": "Mezcla sin procesar",
    "Butterworth_0.5_30Hz": "Filtro Butterworth (0.5–30 Hz)",
    "SVD_Zscore": "SVD (Escalado Z-score, τ*=0.50)",
    "SVD_Fisico_uV": "SVD (Escala física μV, τ*=0.50)",
}
method_colors = {
    "Sin_procesar": "#8C969E",
    "Butterworth_0.5_30Hz": "#F4A261",
    "SVD_Zscore": "#2A9D8F",
    "SVD_Fisico_uV": "#1B365D",
}

x_pos = np.arange(3)
bw = 0.19
for m_idx, m_key in enumerate(methods_order):
    sub_m = df_eval[df_eval["method"] == m_key].groupby("alpha")["snr_fpz_db"].agg(["mean", "std"]).loc[[0.25, 0.50, 0.80]]
    offset = (m_idx - 1.5) * bw
    bars = axes[0].bar(x_pos + offset, sub_m["mean"], bw, yerr=sub_m["std"], capsize=3,
                       label=method_labels[m_key], color=method_colors[m_key], alpha=0.92)
    for b, st in zip(bars, sub_m["std"]):
        val = b.get_height()
        if val >= 0:
            y_txt = val + st + 0.5
            va = "bottom"
        else:
            y_txt = val - st - 0.6
            va = "top"
        axes[0].text(b.get_x() + b.get_width()/2, y_txt, f"{val:+.1f}", ha="center", va=va, fontsize=7.4, fontweight="bold")

axes[0].axhline(0, color="#222222", lw=1.0)
axes[0].set_xticks(x_pos)
axes[0].set_xticklabels(["Baja ($\\alpha = 0.25$)", "Media ($\\alpha = 0.50$)", "Alta ($\\alpha = 0.80$)"])
axes[0].set_xlabel("Nivel de contaminación ocular inyectada ($\\alpha$)")
axes[0].set_ylabel("SNR de reconstrucción en EEG Fpz-Cz (dB)")
axes[0].set_ylim(-14.5, 29.0)
axes[0].set_title("A. SNR absoluta en EEG Fpz-Cz (Media ± DE, n = 183 épocas)")
axes[0].legend(loc="upper right", fontsize=8.0, frameon=True)
axes[0].grid(axis="y", linestyle=":", alpha=0.6)

comps_plot = [
    ("SVD_Fisico_uV vs Sin_procesar", "SVD (μV) vs. Sin procesar", "#1B365D", "o"),
    ("SVD_Fisico_uV vs Butterworth_0.5_30Hz", "SVD (μV) vs. Butterworth", "#2B5B84", "s"),
    ("SVD_Zscore vs Sin_procesar", "SVD (Z-score) vs. Sin procesar", "#2A9D8F", "^"),
    ("SVD_Zscore vs Butterworth_0.5_30Hz", "SVD (Z-score) vs. Butterworth", "#E07A5F", "D"),
]

y_positions = []
y_labels = []
curr_y = 0
for alpha in [0.80, 0.50, 0.25]:
    for comp_key, comp_lbl, col, mk in reversed(comps_plot):
        row = df_boot[(df_boot["alpha"] == alpha) & (df_boot["comparison"] == comp_key)].iloc[0]
        mean_d = row["mean_diff_fpz_db"]
        ci_l   = row["subj_ci95_low_fpz"]
        ci_h   = row["subj_ci95_high_fpz"]
        axes[1].errorbar(mean_d, curr_y, xerr=[[mean_d - ci_l], [ci_h - mean_d]],
                         fmt=mk, color=col, ecolor=col, elinewidth=2.0, capsize=4, ms=5.5,
                         label=comp_lbl if alpha == 0.80 else None)
        axes[1].text(ci_h + 0.35, curr_y, f"{mean_d:+.2f} [{ci_l:+.2f}, {ci_h:+.2f}] dB", va="center", fontsize=7.4, color="#111111")
        y_positions.append(curr_y)
        y_labels.append(f"$\\alpha={alpha:.2f}$ | {comp_lbl.split(' vs.')[0]}")
        curr_y += 1
    curr_y += 0.6

axes[1].axvline(0.0, color="#D93025", linestyle="--", lw=1.6, label="Umbral de hipótesis ($0$ dB)")
axes[1].set_yticks(y_positions)
axes[1].set_yticklabels(y_labels, fontsize=7.8)
axes[1].set_xlabel("Diferencia pareada $\\Delta\\mathrm{SNR} = \\mathrm{SNR}_{\\mathrm{SVD}} - \\mathrm{SNR}_{\\mathrm{Ref}}$ (dB)")
axes[1].set_xlim(-3.0, 26.5)
axes[1].set_ylim(-1.2, curr_y + 2.2)
axes[1].set_title("B. Intervalos de Confianza Bootstrap 95% por sujeto (Fpz-Cz)")
axes[1].legend(loc="upper center", ncol=2, fontsize=7.5, frameon=True)
axes[1].grid(axis="x", linestyle=":", alpha=0.6)

fig.suptitle("Figura 4: Evaluación cuantitativa de SNR y contraste de hipótesis mediante Bootstrap por sujeto (B = 10 000 réplicas)", y=1.03, fontsize=11.5, fontweight="bold")
fig.tight_layout()
save_fig(fig, "fig4_comparacion_snr_bootstrap")

# ==============================================================================
# FIGURA 5: Interpretación Física de Componentes SVD y Proyección en el Subespacio
# ==============================================================================
eval_indices = np.where(eval_mask)[0]
best_demo_idx = None
for i in eval_indices:
    if df_ref.loc[i, "subject_id"] == "SC4071" and df_ref.loc[i, "ptp_eog"] > 350:
        best_demo_idx = i
        break
if best_demo_idx is None:
    best_demo_idx = eval_indices[np.argmax(df_ref.loc[eval_indices, "ptp_eog"].values)]

alpha_demo = 0.50
x_ref_demo = X_ref[best_demo_idx]
eog_demo = x_ref_demo[:, 2]
X_cont_demo = np.column_stack([
    x_ref_demo[:, 0] + alpha_demo * eog_demo,
    x_ref_demo[:, 1] + (BETA_PZOZ * alpha_demo) * eog_demo,
    eog_demo,
])
U_d, S_d, Vt_d = np.linalg.svd(X_cont_demo, full_matrices=False)
for k in range(3):
    if np.corrcoef(U_d[:, k], eog_demo)[0, 1] < 0:
        U_d[:, k] *= -1
        Vt_d[k, :] *= -1

t_sec = np.arange(3000) / FS
fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6))

t_win = (t_sec >= 5.0) & (t_sec <= 25.0)
offset_step = 180.0
comp_colors = ["#D93025", "#1B365D", "#2A9D8F"]
comp_names = [
    f"$\\sigma_1 u_1(t)$ [Modo Ocular EOG, $E_1={(S_d[0]**2/np.sum(S_d**2))*100:.1f}\\%$, $|r|={abs(np.corrcoef(U_d[:,0], eog_demo)[0,1]):.2f}$]",
    f"$\\sigma_2 u_2(t)$ [EEG Fronto-Central, $E_2={(S_d[1]**2/np.sum(S_d**2))*100:.1f}\\%$, $|r|={abs(np.corrcoef(U_d[:,1], eog_demo)[0,1]):.2f}$]",
    f"$\\sigma_3 u_3(t)$ [EEG Parieto-Occipital, $E_3={(S_d[2]**2/np.sum(S_d**2))*100:.1f}\\%$, $|r|={abs(np.corrcoef(U_d[:,2], eog_demo)[0,1]):.2f}$]",
]

for k in range(3):
    sig_k = S_d[k] * U_d[t_win, k]
    y_shift = (2 - k) * offset_step
    axes[0].plot(t_sec[t_win], sig_k + y_shift, color=comp_colors[k], lw=1.2, label=comp_names[k])
    axes[0].axhline(y_shift, color="#CCCCCC", lw=0.7, ls=":")

axes[0].set_yticks([0, offset_step, 2 * offset_step])
axes[0].set_yticklabels(["Modo 3\n($\\sigma_3 u_3$)", "Modo 2\n($\\sigma_2 u_2$)", "Modo 1\n($\\sigma_1 u_1$)"])
axes[0].set_ylim(-60, 720)
axes[0].set_xlabel("Tiempo en la época REM (s)")
axes[0].set_ylabel("Amplitud proyectada en el subespacio ($\\mu$V)")
axes[0].set_title(f"A. Componentes temporales SVD ({df_ref.loc[best_demo_idx, 'subject_id']}, rec #{df_ref.loc[best_demo_idx, 'record_idx']})")
axes[0].legend(loc="upper left", fontsize=7.8, frameon=True)

c1 = S_d[0] * U_d[:, 0]
c2 = S_d[1] * U_d[:, 1]
is_saccade = np.abs(eog_demo) > 2.0 * np.std(eog_demo)

axes[1].scatter(c1[~is_saccade], c2[~is_saccade], s=6, alpha=0.45, color="#1B365D", label="Núcleo basal sin deflexión ocular ($|\\mathrm{EOG}| \\leq 2\\sigma$)")
axes[1].scatter(c1[is_saccade], c2[is_saccade], s=14, alpha=0.75, color="#D93025", label="Excursiones oculares rápidas REM ($|\\mathrm{EOG}| > 2\\sigma$)")
axes[1].scatter(np.zeros_like(c2), c2, s=4, alpha=0.35, color="#2A9D8F", label="Subespacio reconstruido tras anular $\\sigma_1 = 0$")
axes[1].axvline(0, color="#2A9D8F", ls="--", lw=1.4)
axes[1].set_ylim(-38, 56)
axes[1].set_xlabel("Proyección en Componente 1: $\\sigma_1 u_1(t)$ ($\\mu$V) [Eje Ocular]")
axes[1].set_ylabel("Proyección en Componente 2: $\\sigma_2 u_2(t)$ ($\\mu$V) [Eje EEG]")
axes[1].set_title("B. Proyección en el espacio reducido ($\\sigma_1 u_1$ vs. $\\sigma_2 u_2$)")
axes[1].legend(loc="upper right", fontsize=8.0, frameon=True)
axes[1].grid(True, linestyle=":", alpha=0.6)

fig.suptitle("Figura 5: Interpretación física de las componentes SVD y separación geométrica de artefactos en el espacio reducido", y=1.03, fontsize=11.5, fontweight="bold")
fig.tight_layout()
save_fig(fig, "fig5_vectores_singulares_y_subespacio")

# ==============================================================================
# FIGURA 6: Comparación Temporal y Espectral (PSD en Banda Delta) de la Reconstrucción
# ==============================================================================
r_d = np.array([abs(np.corrcoef(U_d[:, k], eog_demo)[0, 1]) for k in range(3)])
S_d_clean = S_d.copy()
S_d_clean[r_d >= TAU_OPT] = 0.0
X_hat_demo = (U_d * S_d_clean) @ Vt_d

x_ref_fpz  = x_ref_demo[:, 0]
x_cont_fpz = X_cont_demo[:, 0]
x_butt_fpz = filtfilt(b_butt, a_butt, x_cont_fpz)
x_svd_fpz  = X_hat_demo[:, 0]

fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6))

t_zoom = (t_sec >= 8.0) & (t_sec <= 22.0)
axes[0].plot(t_sec[t_zoom], x_cont_fpz[t_zoom], color="#D93025", alpha=0.55, lw=1.2, label="Mezcla contaminada $x_{\\mathrm{cont}}(t)$ ($\\mathrm{SNR}=-4.8$ dB)")
axes[0].plot(t_sec[t_zoom], x_butt_fpz[t_zoom], color="#F4A261", alpha=0.85, lw=1.2, ls="--", label="Filtro Butterworth 0.5–30 Hz ($\\mathrm{SNR}=-0.7$ dB)")
axes[0].plot(t_sec[t_zoom], x_ref_fpz[t_zoom],  color="#111111", lw=1.6, label="EEG referencia limpia $x_{\\mathrm{ref}}(t)$")
axes[0].plot(t_sec[t_zoom], x_svd_fpz[t_zoom],  color="#1B365D", lw=1.3, ls="-.", label="Reconstrucción SVD ($\\mathrm{SNR}=+12.9$ dB, $r=0.97$)")

axes[0].set_ylim(-135, 175)
axes[0].set_xlabel("Tiempo (s)")
axes[0].set_ylabel("Amplitud EEG Fpz-Cz ($\\mu$V)")
axes[0].set_title(f"A. Comparación en el dominio del tiempo ({df_ref.loc[best_demo_idx, 'subject_id']}, $\\alpha = 0.50$)")
axes[0].legend(loc="upper left", fontsize=7.8, frameon=True)
axes[0].grid(True, linestyle=":", alpha=0.6)

psd_ref_list, psd_cont_list, psd_butt_list, psd_svd_list = [], [], [], []
for i in eval_indices:
    xr = X_ref[i, :, 0]
    eog_i = X_ref[i, :, 2]
    xc = xr + 0.50 * eog_i
    xb = filtfilt(b_butt, a_butt, xc)
    
    X_m = np.column_stack([xc, X_ref[i, :, 1] + 0.30 * 0.50 * eog_i, eog_i])
    U_i, S_i, Vt_i = np.linalg.svd(X_m, full_matrices=False)
    rk = np.array([abs(np.corrcoef(U_i[:, k], eog_i)[0, 1]) for k in range(3)])
    S_i[rk >= TAU_OPT] = 0.0
    xs = ((U_i * S_i) @ Vt_i)[:, 0]
    
    f_w, pr = welch(xr, fs=FS, nperseg=512)
    _, pc = welch(xc, fs=FS, nperseg=512)
    _, pb = welch(xb, fs=FS, nperseg=512)
    _, ps = welch(xs, fs=FS, nperseg=512)
    psd_ref_list.append(pr)
    psd_cont_list.append(pc)
    psd_butt_list.append(pb)
    psd_svd_list.append(ps)

f_mask = (f_w >= 0.5) & (f_w <= 20.0)
axes[1].axvspan(0.5, 4.0, color="#FFF4E5", alpha=0.8, label="Banda Delta ($0.5 - 4.0$ Hz): solapamiento EOG-EEG")
axes[1].semilogy(f_w[f_mask], np.mean(psd_cont_list, axis=0)[f_mask], color="#D93025", lw=1.5, label="Mezcla contaminada $x_{\\mathrm{cont}}$")
axes[1].semilogy(f_w[f_mask], np.mean(psd_butt_list, axis=0)[f_mask], color="#F4A261", lw=1.5, ls="--", label="Filtro Butterworth (0.5–30 Hz)")
axes[1].semilogy(f_w[f_mask], np.mean(psd_ref_list, axis=0)[f_mask],  color="#111111", lw=1.8, label="EEG referencia limpia $x_{\\mathrm{ref}}$")
axes[1].semilogy(f_w[f_mask], np.mean(psd_svd_list, axis=0)[f_mask],  color="#1B365D", lw=1.6, ls="-.", label="Reconstrucción SVD")

axes[1].set_ylim(0.12, 800)
axes[1].set_xlabel("Frecuencia (Hz)")
axes[1].set_ylabel("Densidad Espectral de Potencia ($\\mu\\mathrm{V}^2/\\mathrm{Hz}$, escala log)")
axes[1].set_title("B. Espectro de potencia promedio en Evaluación ($n = 183$ épocas)")
axes[1].legend(loc="upper right", fontsize=8.0, frameon=True)
axes[1].grid(True, which="both", linestyle=":", alpha=0.6)

fig.suptitle("Figura 6: Reconstrucción temporal y recuperación de la densidad espectral de potencia en banda delta (0.5–4 Hz)", y=1.03, fontsize=11.5, fontweight="bold")
fig.tight_layout()
save_fig(fig, "fig6_reconstruccion_temporal_y_psd")

# ==============================================================================
# FIGURA 7: Análisis de Caso de Éxito vs. Caso de Fallo (Rotación de Subespacios)
# ==============================================================================
sub_025_z = df_eval[(df_eval["alpha"] == 0.25) & (df_eval["method"] == "SVD_Zscore")].set_index("epoch_eval_idx")
sub_025_r = df_eval[(df_eval["alpha"] == 0.25) & (df_eval["method"] == "Sin_procesar")].set_index("epoch_eval_idx")
diff_025 = sub_025_z["snr_fpz_db"] - sub_025_r["snr_fpz_db"]
fail_epoch_idx = int(diff_025.idxmin())
success_epoch_idx = int(best_demo_idx)

fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6))

xr_s = X_ref[success_epoch_idx, :, 0]
eog_s = X_ref[success_epoch_idx, :, 2]
xc_s = xr_s + 0.50 * eog_s
Xm_s = np.column_stack([xc_s, X_ref[success_epoch_idx, :, 1] + 0.15 * eog_s, eog_s])
Us, Ss, Vts = np.linalg.svd(Xm_s, full_matrices=False)
rs_vec = [abs(np.corrcoef(Us[:, k], eog_s)[0, 1]) for k in range(3)]
Ss_c = Ss.copy()
Ss_c[np.array(rs_vec) >= TAU_OPT] = 0.0
xh_s = ((Us * Ss_c) @ Vts)[:, 0]

w_s = (t_sec >= 10.0) & (t_sec <= 20.0)
axes[0].plot(t_sec[w_s], xc_s[w_s], color="#D93025", alpha=0.5, lw=1.2, label="Contaminada ($x_{\\mathrm{cont}}$)")
axes[0].plot(t_sec[w_s], xr_s[w_s], color="#111111", lw=1.6, label="Referencia limpia ($x_{\\mathrm{ref}}$)")
axes[0].plot(t_sec[w_s], xh_s[w_s], color="#1B365D", lw=1.4, ls="-.", label=f"SVD ($|r(u_1,\\mathrm{{EOG}})|={rs_vec[0]:.2f}$, $\\Delta\\mathrm{{SNR}}=+17.7$ dB)")
axes[0].set_ylim(-135, 165)
axes[0].set_xlabel("Tiempo (s)")
axes[0].set_ylabel("Amplitud EEG Fpz-Cz ($\\mu$V)")
axes[0].set_title(f"A. Caso exitoso: Dominancia ocular en $\\sigma_1$ ({df_ref.loc[success_epoch_idx, 'subject_id']})")
axes[0].legend(loc="upper left", fontsize=7.8, frameon=True)
axes[0].grid(True, linestyle=":", alpha=0.6)

X_cont_025 = np.empty_like(X_ref)
X_cont_025[:, :, 0] = X_ref[:, :, 0] + 0.25 * X_ref[:, :, 2]
X_cont_025[:, :, 1] = X_ref[:, :, 1] + 0.075 * X_ref[:, :, 2]
X_cont_025[:, :, 2] = X_ref[:, :, 2]
dev_std_025 = np.std(X_cont_025[dev_mask], axis=(0, 1))

xr_f = X_ref[fail_epoch_idx, :, 0]
eog_f = X_ref[fail_epoch_idx, :, 2]

Uz, Sz, Vtz = np.linalg.svd(X_cont_025[fail_epoch_idx] / dev_std_025, full_matrices=False)
rz_vec = [abs(np.corrcoef(Uz[:, k], eog_f)[0, 1]) for k in range(3)]
Sz_c = Sz.copy()
Sz_c[np.array(rz_vec) >= TAU_OPT] = 0.0
xh_fz = (((Uz * Sz_c) @ Vtz) * dev_std_025)[:, 0]

Up, Sp, Vtp = np.linalg.svd(X_cont_025[fail_epoch_idx], full_matrices=False)
rp_vec = [abs(np.corrcoef(Up[:, k], eog_f)[0, 1]) for k in range(3)]
Sp_c = Sp.copy()
Sp_c[np.array(rp_vec) >= TAU_OPT] = 0.0
xh_fp = ((Up * Sp_c) @ Vtp)[:, 0]

w_f = (t_sec >= 10.0) & (t_sec <= 20.0)
axes[1].plot(t_sec[w_f], xr_f[w_f], color="#111111", lw=1.6, label="Referencia limpia ($x_{\\mathrm{ref}}$)")
axes[1].plot(t_sec[w_f], xh_fz[w_f], color="#E07A5F", lw=1.3, ls="--",
             label=f"Fallo SVD Z-score ($|r(u_1)|={rz_vec[0]:.2f}, |r(u_2)|={rz_vec[1]:.2f}$)")
axes[1].plot(t_sec[w_f], xh_fp[w_f], color="#2A9D8F", lw=1.4, ls="-.",
             label=f"Corrección en escala $\\mu$V ($|r(u_1)|={rp_vec[0]:.2f}, |r(u_2)|={rp_vec[1]:.2f}$)")
axes[1].set_ylim(-28, 54)
axes[1].set_xlabel("Tiempo (s)")
axes[1].set_ylabel("Amplitud EEG Fpz-Cz ($\\mu$V)")
axes[1].set_title(f"B. Caso de fallo por mezcla de subespacios en Z-score ({df_ref.loc[fail_epoch_idx, 'subject_id']}, $\\alpha=0.25$)")
axes[1].legend(loc="upper right", fontsize=7.8, frameon=True)
axes[1].grid(True, linestyle=":", alpha=0.6)

fig.suptitle("Figura 7: Comparación cualitativa entre un caso de desacoplamiento exitoso y un caso de fallo por mezcla de subespacios", y=1.03, fontsize=11.5, fontweight="bold")
fig.tight_layout()
save_fig(fig, "fig7_exito_vs_fallo_svd")
