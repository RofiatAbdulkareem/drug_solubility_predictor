# drug_solubility_predictor

This project is a lightweight machine learning application that predicts the *aqueous solubility (logS)* of molecules based on their chemical structure (SMILES). It leverages molecular descriptors extracted using RDKit and trains a regression model to estimate solubility — a crucial property in drug development and formulation.

## Project Goal

The primary goal was to build a simple, interpretable AI model for healthcare that demonstrates how computational tools can support early-phase drug discovery. By predicting how soluble a compound is in water, researchers can screen drug candidates more effectively and reduce experimental workload.

## What I Achieved

- *Built a regression model* to predict molecular solubility using Random Forest.
- *Extracted molecular descriptors* such as molecular weight, LogP, and hydrogen bond counts using RDKit.
- *Deployed an interactive Streamlit app*, allowing users to input any SMILES string and instantly receive a predicted log solubility.
- *Demonstrated the integration* of cheminformatics, machine learning, and web app deployment in a resource-efficient way — the entire project runs smoothly on a personal laptop.

## Use Cases

- Early drug discovery and screening
- Educational tool for AI in drug design
- Accessible showcase for combining AI with pharmacology

## Limitations

- *Small dataset (ESOL)*: The model was trained on a limited number of compounds, which may restrict generalization to more complex molecules.
- *Descriptor simplicity*: Only a few molecular descriptors were used for the sake of speed and simplicity. This limits model accuracy but keeps it lightweight.
- *No 3D structural information*: The model does not account for stereochemistry or molecular conformation, which can influence solubility.
- *Basic error handling*: Invalid or rare SMILES strings may produce errors or unreliable predictions.

## Future Improvements

- Incorporate larger, more diverse datasets (e.g., AqSolDB).
- Add 3D descriptors or graph-based features using graph neural networks.
- Improve model robustness with cross-validation and hyperparameter tuning.
- Extend to other ADMET properties like permeability or toxicity.

## Acknowledgements

- [ESOL dataset](https://deepchem.io/) by John Delaney
- [RDKit](https://www.rdkit.org/) for cheminformatics
- [Streamlit](https://streamlit.io/) for building the interactive app