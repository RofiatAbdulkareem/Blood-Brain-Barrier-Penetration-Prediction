import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
import numpy as np
import joblib

# Load your trained model
model = joblib.load("random_forest_bbb_model.pkl")

# Load known drugs
@st.cache_data
def load_known_drugs():
    df = pd.read_csv("known_drugs.csv")
    df["BBB Status"] = df["BBB_permeable"].map({1: "Crosses BBB", 0: "Doesn't Cross BBB"})
    return df

known_df = load_known_drugs()

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

# Title and description
st.title("Rofiat's Blood-Brain Barrier Penetration Predictor")
st.markdown("Enter a *SMILES string* to predict if a molecule can cross the blood-brain barrier (BBB).")

# Display known examples with "Try This" buttons
st.subheader("Try Known Drugs")
for i, row in known_df.iterrows():
    cols = st.columns([2, 4, 2])
    cols[0].markdown(f"{row['Drug']}")
    cols[1].markdown(f"{row['SMILES']}")
    if cols[2].button(f"Try This ({row['Drug']})", key=f"btn_{i}"):
        st.session_state["smiles_input"] = row["SMILES"]
        st.session_state["drug_name"] = row["Drug"]

# Get inputs from user
default_smiles = st.session_state.get("smiles_input", "")
default_name = st.session_state.get("drug_name", "")

drug_name = st.text_input("Drug Name (Optional):", value=default_name)
user_smiles = st.text_input("SMILES (Required):", value=default_smiles)

if st.button("Predict"):
    if not user_smiles.strip():
        st.warning("Please enter a SMILES string.")
    else:
        descriptors = compute_descriptors(user_smiles)
        if descriptors is None:
            st.error("❌ Invalid SMILES. Please enter a valid molecule.")
        else:
            prediction = model.predict([descriptors])[0]
            result = "✅ Can cross BBB" if prediction == 1 else "❌ Cannot cross BBB"

            if drug_name.strip() == "":
                st.success(f"Prediction for SMILES '{user_smiles}': {result}")
            else:
                st.success(f"Prediction for '{drug_name}': {result}")
