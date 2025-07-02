import numpy as np
from src.featurization.morgan_featurizer import MorganFeaturizer

def test_transform_shape_and_type():
    smiles = ["CCO", "invalid_smiles"]
    featurizer = MorganFeaturizer(radius=1, n_bits=8)
    feats = featurizer.transform(smiles)
    assert feats.shape == (2, 8)
    assert feats.dtype == np.int8
