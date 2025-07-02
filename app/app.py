import streamlit as st
import pandas as pd
import numpy as np
import joblib
import sys
import os
from rdkit import Chem
from rdkit.Chem import Draw

# Add the project root to Python path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)

from src.featurization.morgan_featurizer import MorganFeaturizer

st.set_page_config(page_title="Biomarker Predictor", layout="wide")
st.title("Chemical Compound Biomarker Predictor")
st.write("Upload a CSV with a 'Drug' or 'SMILES' column, or paste SMILES below.")

# Load model
model_path = os.path.join(project_root, "models", "best_model.pkl")
model = joblib.load(model_path)
featurizer = MorganFeaturizer(radius=2, n_bits=2048)

# Inputs
uploaded = st.file_uploader("Upload CSV", type="csv")
smiles_text = st.text_area("Or paste SMILES (one per line)")

input_smiles = []
if uploaded:
    df = pd.read_csv(uploaded)
    cols = df.columns.str.strip()
    if "SMILES" in cols:
        input_smiles = df["SMILES"].astype(str).tolist()
    elif "Drug" in cols:
        input_smiles = df["Drug"].astype(str).tolist()
    else:
        st.error("CSV must have 'SMILES' or 'Drug' column.")
elif smiles_text:
    input_smiles = [s.strip() for s in smiles_text.splitlines() if s.strip()]

if input_smiles:
    X = featurizer.transform(input_smiles).astype(np.float32)
    preds = model.predict(X)
    results = pd.DataFrame({"SMILES": input_smiles, "Predicted_BioMarker": preds})
    st.subheader("Results")
    st.dataframe(results)

    st.subheader("Structures")
    # Display structures in 3 columns
    cols_per_row = 3
    for i in range(0, len(input_smiles), cols_per_row):
        cols = st.columns(cols_per_row)
        for j in range(cols_per_row):
            if i + j < len(input_smiles):
                smi = input_smiles[i + j]
                p = preds[i + j]
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    img = Draw.MolToImage(mol, size=(300,300))
                    with cols[j]:
                        st.image(img, caption=f"{smi} → {p:.3f}", use_container_width=True)
