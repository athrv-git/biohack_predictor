import pandas as pd
import tempfile
import os
from src.models.train import train_pipeline

def test_train_pipeline_creates_outputs(tmp_path):
    # Create toy train.csv
    train_df = pd.DataFrame({
        "Drug_ID": ["d1", "d2"],
        "Drug": ["CC", "OCC"],
        "Toxicity_Value": [1.0, 2.0]
    })
    valid_df = pd.DataFrame({
        "Drug_ID": ["d3"],
        "Drug": ["CCC"]
    })
    train_csv = tmp_path / "train.csv"
    valid_csv = tmp_path / "valid.csv"
    model_out = tmp_path / "model.pkl"
    submission_out = tmp_path / "sub.csv"
    train_df.to_csv(train_csv, index=False)
    valid_df.to_csv(valid_csv, index=False)

    # Run pipeline
    train_pipeline(str(train_csv), str(valid_csv),
                   str(model_out), str(submission_out))

    # Check outputs
    assert os.path.exists(model_out)
    sub = pd.read_csv(submission_out)
    assert "Bio_Marker_Value" in sub.columns
    assert len(sub) == 1
