import urllib.request
from io import StringIO
import pandas as pd
from scipy.stats import pearsonr
import numpy as np
from rdkit import Chem
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from featurize import Featurizer
import joblib
from urllib.parse import urlparse
import os

def is_url(string:str) -> bool:
    """Check if the string is a valid URL

    Parameters
    ----------
    string : str
        The input string to check

    Returns
    -------
    bool
        True if URL, False otherwise.
    """
    try:
        result = urlparse(string)
        return all([result.scheme, result.netloc])
    except:
        return False

def download_dataset(url):
    with urllib.request.urlopen(url) as f:
        content = f.read().decode('utf-8')
        csvstringio = StringIO(content)
        df = pd.read_csv(csvstringio, sep=",")
        return df


def prepare_model(name, source, y_transform = None, scaler = None):
    """
    Fetches dataset and build a predictive model.
    source could be a path to a csv file with the columns smiles and target;
    or a valid URL that return a csv file with the aforementioned columns.
    """

    my_featurizer = Featurizer()
    if is_url(source):
        df = download_dataset(source)
    elif os.path.isfile(source):
        df = pd.read_csv(source)


    scol = "smile" if "smile" in df.columns else "smiles"
    X = df[scol].apply(lambda x: my_featurizer.featurize(Chem.MolFromSmiles(x)))
    X = np.stack(X)

    y = df["target"].values
    if not (y_transform is None):
        y = y_transform(y)

    print(f"Original shape of X: {X.shape}")
    # Find rows and columns containing NaN or infinite values
    nan_rows = np.isnan(X).any(axis=1)
    inf_rows = np.isinf(X).any(axis=1)


    # Remove rows and columns containing NaN or infinite values
    X = X[~(nan_rows | inf_rows), :]
    y = y[~(nan_rows | inf_rows)]
    
    print(f"New shape of X: {X.shape}")
    if scaler:
        # Call the scaler
        scaler = scaler()
        X = scaler.fit_transform(X)


    X_tr, X_te, y_tr, y_te = train_test_split(X, y, train_size=0.9, random_state=1)

    model = RandomForestRegressor(n_estimators=1000, n_jobs=-1)
    model.fit(X_tr, y_tr)
    y_pred = model.predict(X_te)
    print("Pearson Correlation: {:.3f}".format(pearsonr(y_pred, y_te)[0]))

    # Change n_jobs to 1 in order to avoid warnings during MolDrug run.
    model.n_jobs = 1

    joblib.dump(
        {
            'model':model,
            'scaler':scaler,
        },
        
        "{}.jlib".format(name)
    )


if __name__ == '__main__':

    source1 = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/clearance.csv"
    name1 = "clearance"
    prepare_model(name1, source1)


    source2 = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/hppb.csv"
    name2 = "hppb"
    prepare_model(name2, source2, y_transform=lambda x: np.log10(100-x))

    source3 = 'EGFR_compounds.csv'
    name3 = "egfr"
    prepare_model(name3, source3, scaler=StandardScaler)