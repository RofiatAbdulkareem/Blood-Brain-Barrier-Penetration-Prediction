
# Blood-Brain Barrier Penetration Prediction using Machine Learning

## Overview

This project explores the use of machine learning to predict the ability of small molecules to cross the *blood-brain barrier (BBB)* — a crucial factor in central nervous system (CNS) drug development. Using molecular descriptors derived from SMILES strings, I developed a predictive model that classifies compounds based on their potential to penetrate the BBB.

The goal was to build a lightweight, interpretable model that could serve as a prototype for more advanced AI-driven drug discovery tools.

---

## Motivation

Crossing the BBB is a major hurdle for CNS-active drugs. Traditional lab experiments to assess BBB permeability are expensive and time-consuming. With recent advances in cheminformatics and machine learning, there is a growing interest in building computational tools to screen for BBB-permeable compounds early in the drug development pipeline.

This project is also aligned with the core themes of the PhD position I’m applying for — particularly the intersection of AI, molecular interaction prediction, and drug discovery.

---

## What I Did

- *Data Selection:* I used a cleaned subset of the publicly available *BBBP (Blood-Brain Barrier Penetration)* dataset.
- *Descriptor Engineering:* Molecular descriptors (e.g., molecular weight, LogP, TPSA, hydrogen bond donors/acceptors) were calculated using RDKit.
- *Modeling:* Trained a simple yet effective machine learning model (Random Forest) to classify compounds.
- *Web App:* Deployed an interactive Streamlit app that takes a SMILES string and outputs a prediction of BBB penetration likelihood.

---

## Achievements

- Successfully trained a model with acceptable accuracy on a small, interpretable set of molecular descriptors.
- Built a functional Streamlit app with live prediction capabilities.
- Demonstrated an end-to-end pipeline from data preprocessing to deployment.
- Gained practical experience in cheminformatics, descriptor engineering, and model interpretability.

---

## Limitations

- *Dataset Size & Quality:* The dataset is relatively small, and some compounds may lack critical experimental annotations.
- *Model Generalizability:* The model may not generalize well to structurally novel compounds outside the training set.
- *Descriptor Simplicity:* Only a few descriptors were used; more sophisticated representations (e.g., graph embeddings) may improve performance.
- *Lack of 3D Information:* This model relies on 2D descriptors and does not capture stereochemistry or conformational flexibility.

---

## Future Work

- Incorporate *graph-based molecular representations* using Graph Neural Networks (GNNs).
- Explore *3D structure-based features* such as shape descriptors or docking scores.
- Extend the model to predict other pharmacokinetic properties (e.g., metabolic stability, toxicity).
- Fine-tune hyperparameters and experiment with ensemble or deep learning models.
- Integrate external APIs for real-time compound visualization or data enrichment.

---

## Relevance

This project reflects my core research interests at the intersection of AI and drug development. It highlights my ability to integrate domain knowledge in pharmacology and computational biology to build practical, lightweight tools for molecular prediction — which aligns with the objectives of the PhD program.

---

## Acknowledgment

The dataset used in this project is from the [Therapeutics Data Commons (TDC)](https://tdcommons.ai), which aggregates benchmark datasets for drug discovery.