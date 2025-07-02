import pandas as pd

def load_csv(path: str) -> pd.DataFrame:
    """
    Load a CSV and strip whitespace from column names.
    """
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    return df
