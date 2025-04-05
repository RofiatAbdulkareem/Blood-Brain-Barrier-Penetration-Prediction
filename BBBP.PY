import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from rdkit import Chem
from rdkit.Chem import Descriptors

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix

# Load and clean dataset
df = pd.read_csv("BBBP.csv")
df = df.dropna(subset=['smiles'])

# Compute RDKit descriptors
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

df["features"] = df["smiles"].apply(compute_descriptors)
df = df.dropna(subset=["features"])

# Filter valid SMILES
df["valid_smiles"] = df["smiles"].apply(lambda x: Chem.MolFromSmiles(x) is not None)
df = df[df["valid_smiles"]].copy()
df.drop(columns=["valid_smiles"], inplace=True)

# Split feature columns
descriptor_names = ['MolWt', 'MolLogP', 'NumHDonors', 'NumHAcceptors', 'TPSA']
features_df = pd.DataFrame(df["features"].tolist(), columns=descriptor_names)
df = pd.concat([df.reset_index(drop=True), features_df], axis=1)

# Define X and y
X = df[descriptor_names]
y = df["p_np"]

# Class distribution check
print("Class Distribution:")
print(y.value_counts())
sns.countplot(x=y)
plt.title("Class Distribution")
plt.show()

# Stratified train-test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)

# Random Forest hyperparameter tuning
param_grid = {
    'n_estimators': [100, 200],
    'max_depth': [None, 10, 20],
    'min_samples_split': [2, 5],
    'min_samples_leaf': [1, 2],
    'max_features': ['sqrt', 'log2']
}

grid_search = GridSearchCV(
    RandomForestClassifier(random_state=42),
    param_grid,
    cv=5,
    scoring='accuracy',
    n_jobs=-1,
    verbose=1
)

grid_search.fit(X_train, y_train)
best_model = grid_search.best_estimator_
print("Best Parameters:", grid_search.best_params_)

# Evaluation
y_pred = best_model.predict(X_test)

print("\nClassification Report:")
print(classification_report(y_test, y_pred))
print("Accuracy:", accuracy_score(y_test, y_pred))

# Confusion Matrix
cm = confusion_matrix(y_test, y_pred)
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=['Non-BBB', 'BBB'], yticklabels=['Non-BBB', 'BBB'])
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.title("Confusion Matrix")
plt.show()

# Feature Importance
importances = best_model.feature_importances_
sns.barplot(x=importances, y=descriptor_names)
plt.title("Feature Importances (Random Forest)")
plt.xlabel("Importance")
plt.ylabel("Feature")
plt.show()

# Comparing with other models
models = {
    "RandomForest": best_model,
    "LogisticRegression": LogisticRegression(max_iter=1000),
    "GradientBoosting": GradientBoostingClassifier()
}

print("\nModel Comparison:")
for name, clf in models.items():
    clf.fit(X_train, y_train)
    pred = clf.predict(X_test)
    acc = accuracy_score(y_test, pred)
    print(f"{name}: {acc:.4f}")
