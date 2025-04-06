# Blood-Brain Barrier Penetration Prediction using Machine Learning

## Overview

This project explores the use of machine learning to predict the ability of small molecules to cross the *blood-brain barrier (BBB)* — a critical factor in central nervous system (CNS) drug development. Using molecular descriptors derived from SMILES strings, I developed a lightweight and interpretable model that classifies compounds based on their potential to penetrate the BBB.

An interactive *Streamlit web app* accompanies the model, allowing users to input SMILES strings (or try known drug examples) and receive real-time predictions.


## Motivation

Crossing the BBB is a major hurdle for CNS-active drugs. Traditional lab experiments to assess BBB permeability are expensive and time-consuming. Advances in cheminformatics and machine learning now make it possible to build *fast, computational tools* to screen for BBB-permeable compounds early in the drug development pipeline.

This project aligns with my research interests and PhD aspirations — particularly the intersection of AI, molecular modeling, and drug discovery.


## What I Did

- *Data Selection:* Used a cleaned subset of the publicly available BBBP (Blood-Brain Barrier Penetration) dataset.
- *Descriptor Engineering:* Computed molecular descriptors using RDKit (e.g., molecular weight, LogP, TPSA, hydrogen bond donors/acceptors).
- *Modeling:* Trained a Random Forest classifier and tuned it using cross-validation and GridSearch.
- *Evaluation:* Compared performance across models and visualized feature importances.
- *Web App:* Deployed a user-friendly *Streamlit app* for live BBB predictions using SMILES input.
- *Interactive Features:*
  - Users can enter a custom SMILES string.
  - A curated list of known drugs is provided with a “Try This SMILES” button for fast testing.
  - Results are instantly displayed along with warnings for invalid inputs.


## Achievements

- Built an *end-to-end machine learning pipeline* from data preprocessing to web deployment.
- Designed a *lightweight, interpretable model* suitable for use in low-resource environments.
- Delivered an *interactive app* that enables real-time testing of molecules.
- Implemented *descriptor-based feature importance* to support model interpretability.
- Gained practical experience in cheminformatics, model evaluation, and web deployment.


## Limitations

- *Dataset Size & Quality:* The dataset is relatively small; generalizability may be limited.
- *Descriptor Simplicity:* The model uses a small set of physicochemical descriptors. More complex molecular embeddings or graph-based representations could improve performance.
- *2D Only:* The model does not currently consider stereochemistry or 3D structure.
- *No Probabilistic Output:* Predictions are binary; no confidence score is currently shown.


## Future Work

- Integrate *graph neural networks (GNNs)* for structure-aware modeling.
- Expand the app with *molecular visualization tools* (e.g., RDKit molecule rendering).
- Add *confidence scores* and SHAP-based interpretability features.
- Include *SMILES validation and canonicalization* for robustness.
- Extend the app to predict additional properties (e.g., toxicity, solubility, metabolism).


## Relevance

This project reflects my interests in *molecular neuropharmacology, **AI in healthcare, and **computational drug discovery. It demonstrates my ability to design useful tools that connect **domain knowledge* with *machine learning techniques*, and to deploy them in an accessible format.


## Demo

To run the app:

```bash
streamlit run app.py