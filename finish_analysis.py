# =============================================================================
# InhA Virtual Screening â€” Complete Analysis (no PLIP needed)
# Produces: master results table + publication figures
# =============================================================================
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
BASE       = r"C:\Users\user\Desktop\inha_virtual_screening"
PROJ       = os.path.join(BASE, "project")
OUT_TABLES = os.path.join(PROJ, "tables")
OUT_FIGS   = os.path.join(PROJ, "figures")
os.makedirs(OUT_TABLES, exist_ok=True)
os.makedirs(OUT_FIGS,   exist_ok=True)

# -----------------------------------------------------------------------------
# Hard-coded docking results (from final_docking_summary.csv, verified)
# -----------------------------------------------------------------------------
DOCKING = {
    "Lantadene_A":              {"docking": -10.142, "heavy_atoms": 40},
    "Rutin":                    {"docking":  -9.249, "heavy_atoms": 43},
    "Guaijaverin":              {"docking":  -8.908, "heavy_atoms": 30},
    "Quercitrin":               {"docking":  -8.762, "heavy_atoms": 32},
    "7-Methyljuglone":          {"docking":  -7.999, "heavy_atoms": 14},
    "PT70_Control_Inhibitor":   {"docking":  -7.514, "heavy_atoms": 20},
}

# Plant source mapping (for Table 1)
PLANT_SOURCE = {
    "Lantadene_A":              "Lantana camara (introduced)",
    "Rutin":                    "Multiple Zimbabwean spp.",
    "Guaijaverin":              "Psidium guajava (Mugawha)",
    "Quercitrin":               "Multiple Zimbabwean spp.",
    "7-Methyljuglone":          "Euclea / Diospyros spp. (Mushenje)",
    "PT70_Control_Inhibitor":   "Reference inhibitor (2X22)",
}

# InhA key active-site residues (for reference in text)
INHA_RESIDUES = ["TYR158", "SER94", "GLY96", "PHE149",
                 "MET155", "ALA198", "ILE215", "LEU218",
                 "PHE41", "ILE122", "VAL65", "GLY14"]

# =============================================================================
# STEP 1 â€” Load SwissADME
# =============================================================================
def load_swissadme():
    path = os.path.join(PROJ, "admet", "swissadme.csv")
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found")
        return pd.DataFrame()
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    df["Molecule"] = df["Molecule"].astype(str).str.strip()
    # Map PT70 â†’ PT70_Control_Inhibitor
    df["Molecule"] = df["Molecule"].replace({"PT70": "PT70_Control_Inhibitor"})
    rename = {
        "Molecule":                 "Compound",
        "MW":                       "MW",
        "Consensus Log P":          "LogP",
        "#H-bond donors":           "HBD",
        "#H-bond acceptors":        "HBA",
        "TPSA":                     "TPSA",
        "#Rotatable bonds":         "RotBonds",
        "Lipinski #violations":     "Lipinski_viol",
        "Bioavailability Score":    "Bioavailability",
        "GI absorption":            "GI_absorption",
        "BBB permeant":             "BBB_permeant",
        "Pgp substrate":            "Pgp_substrate",
        "CYP1A2 inhibitor":         "CYP1A2",
        "CYP2C19 inhibitor":        "CYP2C19",
        "CYP2C9 inhibitor":         "CYP2C9",
        "CYP2D6 inhibitor":         "CYP2D6",
        "CYP3A4 inhibitor":         "CYP3A4",
        "Synthetic Accessibility":  "SA_score",
        "PAINS #alerts":            "PAINS",
        "Brenk #alerts":            "Brenk",
    }
    keep = [k for k in rename if k in df.columns]
    return df[keep].rename(columns=rename)

# =============================================================================
# STEP 2 â€” Load pkCSM (tab-separated)
# =============================================================================
def load_pkcsm():
    path = os.path.join(PROJ, "admet", "adme_results.csv")
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found")
        return pd.DataFrame()
    # Try tab first, then comma
    try:
        df = pd.read_csv(path, sep="\t")
        if df.shape[1] < 3:
            raise ValueError
    except Exception:
        df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]

    # Match to compounds via SMILES (from your Compound list)
    smi_to_name = {
        "Cc1cc(O)c2c(c1)C(=O)C=CC2=O":                          "7-Methyljuglone",
        "Oc1cc(O)c2c(c1)oc(c(c2=O)OC1OC(C(C1O)O)O)c1ccc(c(c1)O)O": "Guaijaverin",
        "CC=C(C(=O)OC1CC(C)(C)CC2C1(CCC1(C2=CCC2C1(C)CCC1C2(C)CCC(=O)C1(C)C)C)C(=O)O)C":
                                                                "Lantadene_A",
        "Oc1cc(O)c2c(c1)oc(c(c2=O)OC1OC(C)C(C(C1O)O)O)c1ccc(c(c1)O)O": "Quercitrin",
        "Oc1cc(O)c2c(c1)oc(c(c2=O)OC1OC(COC2OC(C)C(C(C2O)O)O)C(C(C1O)O)O)c1ccc(c(c1)O)O":
                                                                "Rutin",
        "CCCCCc1ccc(cc1)Cc1ccc(o1)C(=O)O":                      "PT70_Control_Inhibitor",
    }
    df["Compound"] = df["SMILES"].map(smi_to_name)
    df = df.dropna(subset=["Compound"])

    keep_cols = {
        "Compound":                         "Compound",
        "AMES toxicity":                    "AMES",
        "hERG I inhibitor":                 "hERG_I",
        "hERG II inhibitor":                "hERG_II",
        "Hepatotoxicity":                   "Hepatotoxicity",
        "Carcinogenicity":                  "Carcinogenicity",
        "Skin Sensitisation":               "Skin_Sens",
        "Oral Rat Acute Toxicity (LD50)":   "LD50_log",
        "Oral Rat Chronic Toxicity (LOAEL)":"LOAEL_log",
        "Intestinal absorption (human)":    "Intestinal_abs",
        "BBB permeability":                 "BBB_perm",
        "CNS permeability":                 "CNS_perm",
        "Total Clearance":                  "Total_clearance",
        "Max. tolerated dose (human)":      "Max_tolerated_dose",
    }
    keep = [k for k in keep_cols if k in df.columns]
    df = df[keep].rename(columns=keep_cols)
    # Convert LD50 from log(mol/kg) to mg/kg
    if "LD50_log" in df.columns:
        df["LD50_mg_per_kg"] = (10 ** df["LD50_log"]).round(1)
    if "LOAEL_log" in df.columns:
        df["LOAEL_log_mg_per_kg"] = (10 ** df["LOAEL_log"]).round(1)
    return df

# =============================================================================
# STEP 3 â€” Merge everything
# =============================================================================
def build_master():
    print("\n" + "=" * 72)
    print("Building master results table")
    print("=" * 72)

    rows = []
    for name, info in DOCKING.items():
        rows.append({
            "Compound":          name,
            "Plant_Source":      PLANT_SOURCE[name],
            "Docking_kcal_mol":  info["docking"],
            "Heavy_Atoms":       info["heavy_atoms"],
            "Ligand_Efficiency": round(-info["docking"] / info["heavy_atoms"], 3),
        })
    master = pd.DataFrame(rows)

    swiss = load_swissadme()
    pkcsm = load_pkcsm()
    if not swiss.empty:
        master = master.merge(swiss, on="Compound", how="left")
    if not pkcsm.empty:
        master = master.merge(pkcsm, on="Compound", how="left")

    master = master.sort_values("Docking_kcal_mol").reset_index(drop=True)
    master.insert(0, "Rank", range(1, len(master) + 1))

    # Verdict column
    def verdict(row):
        if row["Compound"] == "PT70_Control_Inhibitor":
            return "Reference"
        concerns = []
        if pd.notna(row.get("Lipinski_viol")) and row["Lipinski_viol"] > 1:
            concerns.append("Lipinski")
        if row.get("GI_absorption") == "Low":
            concerns.append("Low GI abs")
        if row.get("hERG_II") == "Yes":
            concerns.append("hERG-II")
        if row.get("AMES") == "Yes":
            concerns.append("AMES+")
        if not concerns:
            return "â˜… Promising lead"
        return "âš  " + ", ".join(concerns)
    master["Verdict"] = master.apply(verdict, axis=1)

    out = os.path.join(OUT_TABLES, "master_results.csv")
    master.to_csv(out, index=False)
    print(f"  Saved: {out}\n")
    print(master.to_string(index=False))
    return master

# =============================================================================
# STEP 4 â€” Figures
# =============================================================================
def fig1_docking(master):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    m = master.sort_values("Docking_kcal_mol", ascending=True)
    colors = ["#d62728" if c == "PT70_Control_Inhibitor" else "#1f77b4"
              for c in m["Compound"]]
    bars = ax.barh(m["Compound"], m["Docking_kcal_mol"], color=colors,
                   edgecolor="black", linewidth=0.7)
    ax.axvline(-7.514, color="red", linestyle="--", alpha=0.7,
               label="PT70 control (âˆ’7.51 kcal/mol)")
    for b, v in zip(bars, m["Docking_kcal_mol"]):
        ax.text(v - 0.15, b.get_y() + b.get_height() / 2,
                f"{v:.2f}", va="center", ha="right",
                color="white", fontsize=9, fontweight="bold")
    ax.set_xlabel("Binding energy (kcal/mol)", fontsize=11)
    ax.set_title("Docking scores of top hits against Mtb InhA (PDB 2AQ8)",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="lower left", fontsize=9)
    ax.invert_yaxis()
    plt.tight_layout()
    out = os.path.join(OUT_FIGS, "fig1_docking_scores.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

def fig2_ligand_efficiency(master):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    m = master.sort_values("Ligand_Efficiency", ascending=False)
    colors = ["#d62728" if c == "PT70_Control_Inhibitor" else "#2ca02c"
              for c in m["Compound"]]
    bars = ax.bar(m["Compound"], m["Ligand_Efficiency"],
                  color=colors, edgecolor="black", linewidth=0.7)
    for b, v in zip(bars, m["Ligand_Efficiency"]):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.01,
                f"{v:.3f}", ha="center", fontsize=9, fontweight="bold")
    ax.axhline(0.3, color="gray", linestyle=":",
               label="LE = 0.30 (drug-like threshold)")
    ax.set_ylabel("Ligand efficiency (kcal/mol per heavy atom)", fontsize=11)
    ax.set_title("Ligand efficiency of screened compounds",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    out = os.path.join(OUT_FIGS, "fig2_ligand_efficiency.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

def fig3_size_vs_affinity(master):
    fig, ax = plt.subplots(figsize=(7, 5))
    for _, r in master.iterrows():
        color = "#d62728" if r["Compound"] == "PT70_Control_Inhibitor" else "#1f77b4"
        ax.scatter(r["Heavy_Atoms"], r["Docking_kcal_mol"],
                   s=140, color=color, alpha=0.8, edgecolor="black", linewidth=1)
        ax.annotate(r["Compound"],
                    (r["Heavy_Atoms"], r["Docking_kcal_mol"]),
                    xytext=(7, 5), textcoords="offset points", fontsize=9)
    ax.set_xlabel("Ligand size (heavy atoms)", fontsize=11)
    ax.set_ylabel("Binding energy (kcal/mol)", fontsize=11)
    ax.set_title("Affinity vs. ligand size â€” the ligand efficiency story",
                 fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(OUT_FIGS, "fig3_size_vs_affinity.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

def fig4_admet_heatmap(master):
    # Build a numerical matrix of ADMET risk flags
    rows = []
    for _, r in master.iterrows():
        rows.append({
            "Compound":          r["Compound"],
            "Lipinski viol.":    r.get("Lipinski_viol", 0),
            "Low GI abs.":       1 if r.get("GI_absorption") == "Low" else 0,
            "hERG-II":           1 if r.get("hERG_II") == "Yes" else 0,
            "AMES":              1 if r.get("AMES") == "Yes" else 0,
            "Hepatotoxic":       1 if r.get("Hepatotoxicity") == "Yes" else 0,
            "Carcinogenic":      1 if r.get("Carcinogenicity") == "Yes" else 0,
        })
    df = pd.DataFrame(rows).set_index("Compound")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    sns.heatmap(df.astype(float), annot=True, cmap="RdYlGn_r",
                cbar_kws={"label": "Risk (0 = safe, 1 = alert)"},
                linewidths=0.5, ax=ax)
    ax.set_title("ADMET risk profile of screened compounds",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(OUT_FIGS, "fig4_admet_heatmap.png")
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

# =============================================================================
# STEP 5 â€” Output manuscript-ready summary
# =============================================================================
def write_summary(master):
    path = os.path.join(OUT_TABLES, "results_summary.txt")
    with open(path, "w", encoding="utf-8") as f:
        f.write("=" * 72 + "\n")
        f.write("InhA Virtual Screening â€” Results Summary\n")
        f.write("=" * 72 + "\n\n")

        f.write("TOP 3 DOCKING HITS\n")
        f.write("-" * 40 + "\n")
        top3 = master.nsmallest(3, "Docking_kcal_mol")
        for _, r in top3.iterrows():
            f.write(f"  {r['Compound']:25s} {r['Docking_kcal_mol']:>7.2f} kcal/mol | "
                    f"LE={r['Ligand_Efficiency']:.3f} | {r['Plant_Source']}\n")

        f.write("\nBEST DRUG-LIKE LEAD (0 Lipinski violations, High GI)\n")
        f.write("-" * 40 + "\n")
        dl = master[
            (master["Lipinski_viol"] == 0) &
            (master["GI_absorption"] == "High")
        ]
        for _, r in dl.iterrows():
            f.write(f"  {r['Compound']:25s} MW={r['MW']:.1f}, LogP={r['LogP']:.2f}, "
                    f"LE={r['Ligand_Efficiency']:.3f}\n")

        f.write("\nTOXICITY SUMMARY\n")
        f.write("-" * 40 + "\n")
        for _, r in master.iterrows():
            f.write(f"  {r['Compound']:25s} LD50={r.get('LD50_mg_per_kg', 'n/a')} mg/kg | "
                    f"hERG-II={r.get('hERG_II','?')} | AMES={r.get('AMES','?')}\n")

    print(f"  Saved: {path}")

# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    master = build_master()
    print("\n" + "=" * 72)
    print("Generating figures...")
    print("=" * 72)
    fig1_docking(master)
    fig2_ligand_efficiency(master)
    fig3_size_vs_affinity(master)
    fig4_admet_heatmap(master)
    write_summary(master)
    print("\n" + "=" * 72)
    print("DONE")
    print(f"  Master table:  {os.path.join(OUT_TABLES, 'master_results.csv')}")
    print(f"  Summary:       {os.path.join(OUT_TABLES, 'results_summary.txt')}")
    print(f"  Figures:       {OUT_FIGS}")
    print("=" * 72)
