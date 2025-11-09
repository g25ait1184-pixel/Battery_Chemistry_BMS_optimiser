# Battery_Chemistery–BMS Optimizer

This repository implements a **multi-objective optimization framework** for selecting the most suitable **battery chemistry** and tuning **Battery Management System (BMS) parameters** across three representative energy-storage domains:

- 🚗 **Electric Vehicle (EV)**
- ⚡ **Battery Energy Storage System (BESS)**
- 🔋 **Stationary Energy Storage System (SESS)**

The framework evaluates **five classical optimization methods** on a unified objective function that captures **chemistry performance, BMS constraints, and chemistry–BMS coupling penalties**.


---

## 🔍 Objective Function
> \( J(x,p) = α^T(P^T x) + β^T G(p) + γ^T H(x,p) \)  
> for EV, BESS, and SESS cases respectively.

| Term | Description |
|:--|:--|
| **Chemistry term (αᵀ Pᵀ x)** | Evaluates each chemistry based on normalized metrics (Energy Density, Cost, Safety, Cycle Life, RTE, DoD, C-rate, Temperature performance). |
| **BMS term (βᵀ G(p))** | Evaluates BMS configuration parameters such as voltage, current, SOC, and temperature range. |
| **Penalty term (γᵀ H(x,p))** | Applies chemistry–BMS coupling constraints such as over-current, temperature, or imbalance penalties. |
See docs/docs.md for detailed mathematical derivation.
---

## 🧱 Project Structure

```
chem_bms_optimizer/
│
├── optimizer.py # Core objective, penalties, dataset loader
├── run_all_methods.py # Classical optimization runners + convergence logic
├── notebooks/
│ └── Run_All_Methods.ipynb # Jupyter demo + radar charts + summary tables
│
├── sample_data/ # Application and chemistry datasets
│ ├── alpha_weights.csv # α-weight priorities for EV/BESS/SESS
│ ├── beta_weights.csv # β-weight importance for BMS variables
│ ├── gamma_weights.csv # γ-weight scaling for coupling penalties
│ ├── chemistry_metrics.csv # Normalized performance metrics (Energy, Cost, etc.)
│ ├── chemistry_envelope.csv # Recommended envelopes (C-rate, temperature, etc.)
│ ├── bms_bounds.csv # Global parameter bounds
│ ├── bms_bounds_EV.csv # EV-specific bounds
│ ├── bms_bounds_BESS.csv # BESS-specific bounds
│ └── bms_bounds_SESS.csv # SESS-specific bounds
│
├── outputs/
│ ├── summary_methods.csv # Final results (per app × method)
│ ├── trace_EV_SQP_KKT.csv # Iteration trace logs
│ └── tables/ # Markdown + LaTeX formatted tables
│
└── README.md # ← You are here
```

---
🛠️ Installation
# Clone the repository
git clone https://github.com/g25ait1184-pixel/Battery_Chemistry_BMS_optimiser.git
cd Battery_Chemistry_BMS_optimiser

# Install dependencies
pip install -r requirements.txt

## ⚙️ How to Run

### ▶️ Run in terminal

```bash
python run_demo.py
# or
python run_all_methods.py
```

### 📓 Run in Jupyter Notebook

Open `Run_All_Methods.ipynb` and execute cells sequentially.  
It will generate convergence plots (`trace_*.csv`) and α-weight radar charts.

---
## 💻 Usage

### 🔹 Basic Usage
Run a full optimization for all three applications (**EV**, **BESS**, **SESS**):

bash
python run_all_methods.py

# --- Run optimizer for all applications ---
results = run_all_apps(PROJECT)

# --- Load the standardized summary output ---
summary = pd.read_csv(PROJECT / "outputs" / "summary_methods.csv")



## 🧮 Optimization Methods Implemented
- **Steepest Descent**
- **Quasi-Newton (BFGS)**
- **Sequential Quadratic Programming (Box)**
- **Newton’s Method (Diagonal Hessian)**
- **Line Search (Armijo)**

Each method minimizes the multi-objective cost \(J(x,p)\) across the three application types (EV, BESS, ESS).

---

## 🏆 Winners by Application (J ↓ = Lower Cost, S ↑ = Higher Suitability)

| Application | Method               | Best Chemistry | Total J | Suitability (S) | Iter | Time (s) |
|:-------------|:---------------------|:----------------|:-------:|:----------------:|:----:|:---------:|
| **BESS**     | Newton’s Method      | LFP             | 0.197   | 0.827            | 193  | 166.11 |
| **EV**       | Quasi-Newton (BFGS)  | SCiB            | 0.231   | 0.774            | 11   | 40.13 |
| **SESS**     | Newton’s Method      | LFP             | 0.225   | 0.799            | 193  | 195.23 |

> **Interpretation:**  
> - For **BESS**, Newton’s Method identifies **LFP** as the most cost-efficient and stable chemistry.  
> - For **EV**, Quasi-Newton (BFGS) favors **SCiB**, balancing energy density and high-rate performance.  
> - For **SESS**, Newton’s Method again prefers **LFP**, reflecting reliability and cycle-life dominance.

---

## 🔋 Battery Chemistry Performance Matrix

| Chemistry | Energy Density | Affordability | Cycle Life | Safety | RTE | DoD | C-rate | Low-Temp |
|------------|:--------------:|:-------------:|:-----------:|:------:|:---:|:---:|:------:|:--------:|
| **LFP**    | 0.70 | 0.70 | 0.85 | 0.95 | 0.94 | 0.85 | 0.70 | 0.65 |
| **NMC**    | 0.85 | 0.55 | 0.75 | 0.80 | 0.93 | 0.85 | 0.85 | 0.55 |
| **NCA**    | 0.88 | 0.50 | 0.70 | 0.75 | 0.93 | 0.85 | 0.90 | 0.50 |
| **LCO**    | 0.82 | 0.35 | 0.55 | 0.60 | 0.90 | 0.80 | 0.80 | 0.45 |
| **Na-ion** | 0.55 | 0.80 | 0.70 | 0.90 | 0.88 | 0.80 | 0.60 | 0.60 |
| **NiMH**   | 0.45 | 0.65 | 0.65 | 0.85 | 0.80 | 0.75 | 0.55 | 0.55 |
| **Li-S**   | 0.95 | 0.40 | 0.30 | 0.50 | 0.85 | 0.75 | 0.80 | 0.40 |
| **LTO**    | 0.60 | 0.45 | 0.95 | 0.98 | 0.95 | 0.85 | 0.90 | 0.70 |

> **Note:** Higher values (closer to 1.0) indicate superior performance for the corresponding property.



## 📊 Quick Visual Interpretation

| Metric (Key Property)     | EV Focus | BESS Focus | SESS Focus |
|----------------------------|:--------:|:-----------:|:-----------:|
| **E – Energy Density**     | ★★★★★    | ★★☆☆☆      | ★★★☆☆      |
| **S – Safety**             | ★★☆☆☆    | ★★★★☆      | ★★★★☆      |
| **L – Cycle Life**         | ★★☆☆☆    | ★★★★★      | ★★★★☆      |
| **C_inv – Cost (1/Cost)**  | ★★☆☆☆    | ★★★★☆      | ★★★☆☆      |
| **RTE – Efficiency**       | ★★★☆☆    | ★★★☆☆      | ★★★☆☆      |
| **DoD – Depth of Discharge** | ★★☆☆☆  | ★★☆☆☆      | ★★☆☆☆      |
| **T – Low-Temperature**    | ★★☆☆☆    | ★★☆☆☆      | ★★☆☆☆      |
| **R – C-rate (Power Capability)** | ★★★★★ | ★☆☆☆☆ | ★★☆☆☆ |

## 🧾 Requirements

Install Python ≥ 3.10 and dependencies:

```bash
pip install -r requirements.txt
```

---


## 🪪 License

MIT License © 2025 svk  
Use freely with attribution for academic and research purposes.

---

## 📚 Citation (if used in papers)

> g25ait1184@iitj.ac.in, *Multi-Objective Optimization of Battery Performance Parameters for EV, BESS and ESS Using Classical Optimization Techniques*, 2025.  
> GitHub: [https://github.com/g25ait1184-pixel/Battery_Chemistry_BMS_optimiser)

---

**Maintainer:**  g25ait1184@iitj.ac.in 
**Language:** Python 3.13  
**Platform:** Windows / Jupyter Notebook
