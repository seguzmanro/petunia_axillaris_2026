#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import sys

if len(sys.argv) < 3:
    print("Usage: python format_glmm_table.py <input_csv> <output_prefix>")
    sys.exit(1)

input_csv = sys.argv[1]
out_prefix = sys.argv[2]

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(input_csv)

# -----------------------------
# 1. Select top models
# -----------------------------
DELTA_THRESHOLD = 4  # you can change to 2 for stricter filtering

top_models = df[df["deltaDIC"] <= DELTA_THRESHOLD].copy()

# Sort by deltaDIC
top_models = top_models.sort_values("deltaDIC")

# -----------------------------
# 2. Create model comparison table
# -----------------------------
model_table = top_models[[
    "model", "deltaDIC", "DICweight", "R2_marginal"
]].copy()

model_table.columns = ["Model", "ΔDIC", "Weight", "R² (marginal)"]

model_table.to_csv(f"{out_prefix}_model_comparison.csv", index=False)

# -----------------------------
# 3. Extract effect sizes
# -----------------------------
# Identify columns
post_cols = [c for c in df.columns if c.startswith("post_mean_")]
low_cols  = [c for c in df.columns if c.startswith("lower_95_")]
up_cols   = [c for c in df.columns if c.startswith("upper_95_")]
p_cols    = [c for c in df.columns if c.startswith("pMCMC_")]

# Function to extract variable name
def clean_name(col):
    return re.sub(r"^(post_mean_|lower_95_|upper_95_|pMCMC_)", "", col)

# Build long-format table
rows = []

for _, row in top_models.iterrows():
    model_name = row["model"]
    
    for col in post_cols:
        var = clean_name(col)
        beta = row[col]
        
        if pd.isna(beta):
            continue
        
        lower = row.get(f"lower_95_{var}", np.nan)
        upper = row.get(f"upper_95_{var}", np.nan)
        pval  = row.get(f"pMCMC_{var}", np.nan)
        
        # Skip intercept (optional, comment out if you want it)
        if var == "(Intercept)":
            continue
        
        rows.append({
            "Model": model_name,
            "Predictor": var,
            "β": beta,
            "CI_lower": lower,
            "CI_upper": upper,
            "pMCMC": pval
        })

effects_df = pd.DataFrame(rows)

# -----------------------------
# 4. Clean formatting
# -----------------------------
# Round values for publication
effects_df["β"] = effects_df["β"].round(3)
effects_df["CI_lower"] = effects_df["CI_lower"].round(3)
effects_df["CI_upper"] = effects_df["CI_upper"].round(3)
effects_df["pMCMC"] = effects_df["pMCMC"].round(3)

# Create CI string
effects_df["95% CI"] = (
    effects_df["CI_lower"].astype(str) + ", " +
    effects_df["CI_upper"].astype(str)
)

# Reorder columns
effects_df = effects_df[[
    "Model", "Predictor", "β", "95% CI", "pMCMC"
]]

# Sort by model then p-value
effects_df = effects_df.sort_values(["Model", "pMCMC"])

effects_df.to_csv(f"{out_prefix}_effects_table.csv", index=False)

# -----------------------------
# 5. Optional: best model table
# -----------------------------
best_model_name = df.loc[df["deltaDIC"].idxmin(), "model"]

best_effects = effects_df[effects_df["Model"] == best_model_name]

best_effects.to_csv(f"{out_prefix}_best_model_effects.csv", index=False)

# -----------------------------
print("Done!")
print(f"Model table: {out_prefix}_model_comparison.csv")
print(f"Effects table: {out_prefix}_effects_table.csv")
print(f"Best model table: {out_prefix}_best_model_effects.csv")