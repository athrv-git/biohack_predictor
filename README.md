# 🧬 Biohack Predictor

This project aims to predict biomarker values from compound SMILES using classical machine learning models and RDKit-based molecular fingerprinting.

---

## 📂 Project Structure

```
biohack_predictor/
├── data/
│   └── raw/
│       ├── train.csv
│       └── valid.csv
├── models/
│   └── best_model.pkl
├── src/
│   ├── featurization/
│   │   └── morgan_featurizer.py
│   ├── models/
│   │   ├── train.py
│   │   └── predict.py
│   └── utils/
│       └── io.py
├── Submission.csv
└── app/
    └── app.py
```

---

## 🧹 1. Data Processing & Cleaning

Data loading and featurization are handled using:
- **`src/featurization/morgan_featurizer.py`**: Converts SMILES strings into Morgan fingerprints (2048-bit vectors).
- **`src/utils/io.py`**: Contains helper function `load_csv()` used for safe CSV reading.

---

## 🤖 2. Models Used

Model training and evaluation are implemented in:
- **`src/models/train.py`**:
  - Trains and evaluates 4 models:  
    - `XGBoost`
    - `RandomForest`
    - `Ridge Regression`
    - `MLP Regressor`
  - Selects best model using MAE on a hold-out validation set.
  - Saves the best model to `models/best_model.pkl`.

Model inference is handled by:
- **`src/models/predict.py`**:
  - Loads saved model and predicts biomarker values for validation/test data.

---

## 📊 3. Output

- The **best model** (selected based on MAE): _e.g., XGBoost_
- **Best MAE** on hold-out set: _reported in console during training (e.g., 0.2567)_
- Final predictions are saved in: `Submission.csv`

---

## ▶️ 4. How to Run

### 📌 Prerequisites

Install dependencies:
```bash
pip install -r requirements.txt
```

Or for the app:
```bash
pip install -r app/requirements-app.txt
```

---

### 🏋️ Train the Model

```bash
python src/models/train.py \
    --train-csv data/raw/train.csv \
    --valid-csv data/raw/valid.csv \
    --model-out models/best_model.pkl \
    --submission-out Submission.csv
```

---

### 📈 Make Predictions

```bash
python src/models/predict.py \
    --input-csv data/raw/valid.csv \
    --model-path models/best_model.pkl \
    --output-csv new_submission.csv
```

---

### 🌐 Run Streamlit App (Optional)

```bash
cd app
streamlit run app.py
```

---

## 💡 Notes

- Uses **RDKit** for molecule parsing and fingerprinting.
- Modularized code for reusability and extensibility.
- Predictions saved in `Submission.csv` follow the required format: `Drug_ID, Drug, Bio_Marker_Value`.
