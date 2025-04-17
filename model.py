import pandas as pd
url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv"
df = pd.read_csv(url)

df.head()
df = df.rename(columns={"measured log solubility in mols per litre":"logS"})
from rdkit import Chem
from rdkit.Chem import Descriptors
def featurize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol)
    }
features = df["smiles"].apply(featurize)
features_df = pd.DataFrame(features.tolist())

data = pd.concat([features_df, df["logS"]],axis=1)
data.head()
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score

X = data.drop("logS", axis=1)
y = data["logS"]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
model = RandomForestRegressor()
model.fit(X_train, y_train)

y_pred = model.predict(X_test)

print("R2 Score:", r2_score(y_test, y_pred))
print("RMSE:", mean_squared_error(y_test, y_pred))
import pickle

with open("solubility_model.pkl", "wb") as f:
    pickle.dump(model,f)