import numpy as np
import warnings
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

class MorganFeaturizer:
    """
    Featurizer to convert SMILES to Morgan fingerprint vectors.
    """
    def __init__(self, radius: int = 2, n_bits: int = 2048):
        self.radius = radius
        self.n_bits = n_bits

    def transform(self, smiles_list):
        features = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                # fallback zero vector if parsing fails
                features.append(np.zeros(self.n_bits, dtype=np.int8))
            else:
                # Suppress deprecation warnings temporarily
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                        mol, self.radius, nBits=self.n_bits
                    )
                arr = np.array(fp, dtype=np.int8)
                features.append(arr)
        return np.stack(features)
