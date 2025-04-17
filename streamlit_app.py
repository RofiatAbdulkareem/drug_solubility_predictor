import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors
import pickle
import pandas as pd

# Load trained model
with open("solubility_model.pkl", "rb") as f:
    model = pickle.load(f)

def featurize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return pd.DataFrame([{
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol)
    }])

st.title("Molecular Solubility Predictor")
smiles_input = st.text_input("Enter SMILES string")

if st.button("Predict"):
    try:
        features = featurize(smiles_input)
        prediction = model.predict(features)[0]
        st.success(f"Predicted logS: {prediction:.2f}")
    except:
        st.error("Invalid SMILES string or error in prediction.")