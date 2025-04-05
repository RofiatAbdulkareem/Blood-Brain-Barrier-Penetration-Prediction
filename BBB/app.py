import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
import numpy as np
import joblib

# Load your trained model
model = joblib.load("random_forest_bbb_model.pkl")

# Define descriptor function
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.TPSA(mol),
    ]

# Streamlit UI
st.title("🧠 Blood-Brain Barrier Penetration Predictor")
st.write("Enter a drug name and SMILES string to predict if the molecule can cross the blood-brain barrier (BBB).")

# Drug name input
drug_name = st.text_input("Drug Name (Optional):")

# SMILES input
user_smiles = st.text_input("SMILES:", "CCOC(=O)c1ccc(cc1)N")

if st.button("Predict"):
    descriptors = compute_descriptors(user_smiles)
    
    if descriptors is None:
        st.error("❌ Invalid SMILES. Please enter a valid molecule.")
    else:
        prediction = model.predict([descriptors])[0]
        result = "✅ Can cross BBB" if prediction == 1 else "❌ Cannot cross BBB"
        
        if drug_name.strip() == "":
            st.warning("Drug name not provided.")
            st.success(f"Prediction for SMILES '{user_smiles}': {result}")
        else:
            st.success(f"Prediction for '{drug_name}': {result}")
