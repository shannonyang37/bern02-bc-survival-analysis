## **Breast Cancer Survival Analysis**

### Table of Contents
- [1. Project Overview](#1-project-overview)
- [2. Data Availability & Citation](#2-data-availability--citation) 
- [3. Project Workflow](#3-project-workflow)
- [4. Key Results](#4-key-results)
- [5. Limitations and Future Work](#5-limitations-and-future-work)
- [6. Repository Structure](#6-repository-structure)
- [7. Installation](#7-installation) 
- [8. License](#8-license)
- [9. Version Pinning & Risk Assessment](#9-version-pinning-and-risk-assessment)
- [10. FAIR Data Principles Assessment](#10-fair-data-principles-assessment)


---
### 1. Project Overview
This project aims to investigate the clinical determinants of mortality in breast cancer patients. Utilizing Pandas, Seaborn, and Matplotlib, a custom reproducible statistical pipeline was engineered to conduct univariate exploratory data analysis across 13 clinical factors. By coupling visual distributions with rigorous mathematical validation (including Chi-Square, ANOVA, Welch's T-Test, and Multivariate Log-Rank tests), this analysis identifies the statistically significant prognostic drivers of patient survival.

---
### 2. Data Availability & Citation

The dataset used in this analysis is the de-identified **SEER Breast Cancer Dataset** originally curated and published by Jing Teng on IEEE Dataport, accessed via Kaggle.

- **Persistent Identifier (DOI):** [https://dx.doi.org/10.21227/a9qy-ph35](https://dx.doi.org/10.21227/a9qy-ph35)
- **Primary Source:** Jing Teng (2019), *SEER Breast Cancer Data*, IEEE Dataport.
- **Repository Mirror:** [Sujith Mandala on Kaggle](https://www.kaggle.com/datasets/sujithmandala/seer-breast-cancer-data)
- **Original Registry Authority:** National Cancer Institute (NCI), Surveillance, Epidemiology, and End Results (SEER) Program.

**Dataset Citation:**
> Jing Teng, (January 18, 2019). "SEER Breast Cancer Data", *IEEE Dataport*, doi: https://dx.doi.org/10.21227/a9qy-ph35.
---

### 3. Project Workflow
The project is broken down into three Jupyter Notebooks:
1. **Data Cleaning and Preprocessing (`1_data_cleaning.ipynb`):** Cleaned and preprocessed the SEER Breast Cancer Dataset.
2. **Defining Analysis Tools (`2_analysis_tools.ipynb`):** Serves as a backend graphical and statistical engine, designed functions to automate the exploratory data analysis. 
3. **Exploratory Data Analysis (`3_exploratory_analysis.ipynb`):** Deployed the analysis tools in Notebook 2 to perform univariate exploratory data analysis to analyse the 13 clinical factors of mortality in breast cancer patients. 

---
### 4. Key Results
The univariate exploratory data analysis revealed a crucial clinical distinction between the likelihood of mortality and the velocity of mortality among the 13 clinical factors evaluated.

- **Universal Prognostic Significance:** Every clinical factor analyzed—spanning demographic, anatomical, and biological variables—proved to be a highly statistically significant determinant of overall survival probability. Factors such as advanced T Stage, higher Tumor Grade, and negative hormone receptor status all drastically increased a patient's absolute risk of mortality.
- **Selective Acceleration of Terminal Decline:** While all 13 factors strongly influenced whether a patient would survive, only a specific subset of biological and systemic factors mathematically accelerated the timeline to death once the disease became terminal. Factors indicating systemic metastasis, specifically N Stage (lymphatic spread), A Stage (distant metastasis), and negative Estrogen or Progesterone receptor status, significantly shortened the median survival time of the deceased cohort.
---
### 5. Limitations and Future Work
The primary limitation of this project is the absence of multivariate analysis. While exploring factors in isolation successfully reveals individual statistical associations with patient survival, it does not account for potential confounding variables or multi-collinearity. Furthermore, because the clinical factors were evaluated individually, it is impossible to rank their prognostic strength or determine which factor is the single most powerful driver of mortality. 

In the future, we could introduce multivariate analysis such as Cox Proportional Hazards model for the simultaneous evaluation of all 13 clinical factors. We could also utilize machine learning techniques to rank the features to explicitly isolate the most dominant clinical determinants of breast cancer mortality.

---
### 6. Repository Structure
```
bc-survival-analysis/
├── .gitignore                      # Excludes cache files and virtual environments
├── LICENSE                         # MIT License
├── README.md                       # Documentation, risk assessment & FAIR analysis
├── environment.yml                 # Conda environment definition with version pins
├── requirements.txt                # Pip dependency specifications
├── data/                           # Clinical datasets 
│   ├── raw/
    │    └── SEER_Breast_Cancer_Dataset.csv
    └── processed/
          └── data_clean.csv
└── notebooks/                      # Reproducible workflows
    ├── 1_data_cleaning.ipynb
    ├── 2_analysis_tools.ipynb
    ├── 3_exploratory_analysis.ipynb
    └── analysis_tools.py 
```
---
### 7. Installation
**Prerequisites:** 
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/)
- Git

**1. Clone the repository:**
```
git clone https://github.com/shannonyang37/bc-survival-analysis.git
cd bc-survival-analysis
```

**2. Set up a virtual environment:**
```
conda env create -f environment.yml
conda activate bc-survival-env
```

Alternatively, if using standard virtual environments with pip:
```
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

**3. Launch Jupyter Notebook:**
```
jupyter notebook
```

Run the notebooks in sequential order to replicate the study:
* 1_data_cleaning.ipynb
* 2_analysis_tools.ipynb
* 3_exploratory_analysis.ipynb

---
### 8. License
Distributed under the MIT License. See `LICENSE` for more information.

---
### 9. Version Pinning and Risk Assessment
The environment for this project is managed using Conda to guarantee reproducibility and cross-platform compatibility across dependencies through the `environment.yml` file.
- **Semantic version range pinning:** Range pinning (>= and <) is used instead of exact pinning (==).The lower bounds guarantee that the required features are available, while upper bounds restricts updates to below the next major version to prevent silent code breakage when dependencies upgrade automatically.
- **Environment pruning:** Unused exploratory packages were audited and removed once the final pipeline ws established to keep the environment lightweight.
- **Modular code structure:** Rather than copying and pasting functions across multiple notebooks, the functions were extracted into a dedicated Python module (`.py`) 
- **Clean notebooks:** Redundant code and duplicate imports were removed so the notebooks stay clean and import directly from the helper script.

**Risk assessment:**
- **`lifelines` (High risk):** As a specialized domain library, minor version updates will alter the function arguments, plotting outputs, and summary table structures. Therefore, it is strictly pinned to `>=0.30.3,<0.31.0`.
- **Other packages (Low to medium risk):** Core libraries like `pandas`, `scipy`, `matplotlib`, and `seaborn` are more stable across minor releases. They are upper-bounded to the next major version (e.g., pandas <3.0.0) to avoid major breaking syntax shifts while still allowing safe bug and security fixes.
---
### 10. FAIR Data Principles Assessment
**1. Findable (F):**
This workflow is hosted on an open-source, publicly accessible GitHub repository with structured metadata, relevant topic keywords, descriptive file names, and a comprehensive `README.md` to make the project easy to find and navigate. The input dataset has a persistent, globally recognized DOI (DOI: 10.21227/a9qy-ph35) on IEEE Dataport and is mirrored on Kaggle with structured metadata.

**2. Accessible (A):**
All code, notebooks, and datasets can be retrieved via standard, open protocols (https, Git) without paywalls or proprietary client software. The reproduction instructions are provided in `README.md`.

**3.Interoperable (I):**
Accessible, shared, and broadly used data formats are used (eg: CSV, jupyter notebooks). To support different environments, both `environment.yml`(for Conda) and `requirements.txt` (for pip) are provided.

**4. Reusable (R):**
The repository is licensed under the MIT License, allowing unrestricted reuse, modification, and academic reproduction. Pinned dependency files ensure that anyone running the code will reproduce identical results. Furthermore, the clinical dataset is fully de-identified to ensure ethical and safe reuse. The dataset source is documented with citations, providing provenance and credit to the original authors.
