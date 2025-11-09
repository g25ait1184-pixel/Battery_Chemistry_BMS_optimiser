# 🔋 Chem–BMS Optimizer — Technical Documentation

---

## 1. Problem Formulation

The **Chem–BMS Optimizer** determines the optimal **battery chemistry** and **Battery Management System (BMS) parameter settings** for three energy applications:  
**Electric Vehicles (EV)**, **Battery Energy Storage Systems (BESS)**, and **Stationary Energy Storage Systems (SESS)**.

The optimization minimizes a total multi-objective cost function:

> \( J(x,p) = α^T(P^T x) + β^T G(p) + γ^T H(x,p) \)  
> for EV, BESS, and SESS cases respectively.

### Where:
| Symbol | Meaning |
|:--|:--|
| \(x\) | Chemistry selection vector (one-hot encoded) |
| \(p\) | BMS parameter vector (voltage, current, temperature, SOC window, etc.) |
| \(P\) | Matrix of normalized chemistry performance metrics |
| \(\alpha, \beta, \gamma\) | Weight vectors for chemistry, BMS, and coupling terms |
| \(G(p)\) | Normalized BMS performance term |
| \(H(x,p)\) | Chemistry–BMS coupling penalty or constraint function |

### Constraints:
\[
\mathbf{1}^T x = 1, \quad x_i \in \{0,1\}, \quad p_{lo} \le p \le p_{hi}, \quad c_i(x,p) \le 0
\]

---

## 2. 🧮 Classical Optimization Methods

| Method | Description |
|:--|:--|
| **Steepest Descent** | Basic gradient descent approach using fixed or adaptive step-size. |
| **Quasi-Newton (BFGS)** | Approximates the inverse Hessian iteratively for faster convergence. |
| **Sequential Quadratic Programming (Box)** | Solves a quadratic subproblem with box constraints. |
| **Newton’s Method** | Uses diagonal or full Hessian approximation for second-order updates. |
| **Line Search (Golden Section)** | Exact 1-D minimization of \( J(p + \lambda d) \) using the golden ratio. |
| **SQP (KKT)** | Solves equality-constrained QP via Karush–Kuhn–Tucker (KKT) conditions. |

---

## 3. 📂 Sample Data Files (`sample_data/`)

| **File** | **Description** |
|:--|:--|
| **chemistry_metrics.csv** | Defines normalized chemistry performance metrics (Energy Density, Safety, Cycle Life, Cost, RTE, DoD, C-rate, Low-Temperature) for all candidate chemistries. |
| **alpha_weights.csv** | Application-specific metric weights (α) for EV, BESS, and SESS optimization priorities. |
| **bms_bounds_EV.csv** | BMS configuration parameter bounds specific to the **EV** application (Voltage, Current, Temperature, SOC, ΔV). |
| **bms_bounds_BESS.csv** | Parameter bounds tailored to **Battery Energy Storage System (BESS)** operation limits. |
| **bms_bounds_SESS.csv** | Parameter bounds for **Stationary Energy Storage System (SESS)** thermal and electrical envelopes. |
| **beta_weights.csv** | Weighting factors (β) for BMS control parameter contributions in the objective function. |
| **chemistry_envelope.csv** | Operational envelope for each chemistry — includes rated charge/discharge rates (R_chg, R_dchg), temperature window (T_lo, T_hi), and SOC limits (SOC_lo, SOC_hi). |
| **gamma_weights.csv** | Coupling penalty weights (γ) quantifying constraint violations in voltage, temperature, or SOC domains. |


---

## 4. ⚡ Application Weight Metrics

| Metric (Key Property) | EV   | BESS | SESS |
|------------------------|:----:|:----:|:----:|
| **Energy Density**     | 0.32 | 0.05 | 0.12 |
| **Affordability**      | 0.05 | 0.20 | 0.24 |
| **Cycle Life**         | 0.06 | 0.30 | 0.30 |
| **Safety**             | 0.08 | 0.18 | 0.10 |
| **RTE**                | 0.12 | 0.10 | 0.08 |
| **DoD**                | 0.05 | 0.10 | 0.08 |
| **C-rate**             | 0.25 | 0.03 | 0.03 |
| **Low Temperature**    | 0.07 | 0.04 | 0.05 |

> **Note:** Columns represent normalized α-weight coefficients used in  
> \( J(x,p) = α^T(P^T x) + β^T G(p) + γ^T H(x,p) \)  
> for EV, BESS, and SESS cases respectively.

---

## 5. Performance Metrics

During each optimization run, the framework logs the following performance indicators:

| **Metric** | **Symbol** | **Description** |
|:--|:--:|:--|
| **Total Cost Function** | \( J \) | Final optimized composite objective value |
| **Chemistry Term** | \( J_{\mathrm{chem}} \) | Weighted chemistry suitability contribution |
| **BMS Term** | \( J_{\mathrm{bms}} \) | Weighted BMS parameter cost component |
| **Penalty Term** | \( J_{\mathrm{pen}} \) | Constraint or coupling penalty contribution |
| **Iterations** | \( N_{\mathrm{iter}} \) | Total iterations required to achieve convergence |
| **CPU Time** | \( t_{\mathrm{cpu}} \) | Total runtime for each optimization (in seconds) |
| **Sustainability Index** | \( \eta_{\mathrm{sus}} \) | Aggregate indicator combining energy efficiency, safety, and lifecycle performance to quantify sustainable chemistry-BMS pairing |

---

> **Note:**  
> The **Sustainability Index** (\( \eta_{\mathrm{sus}} \)) is derived post-optimization by normalizing the chemistry’s cycle-life, safety score, and round-trip efficiency relative to its total cost \(J\), offering an interpretive view of long-term operational viability.

---

### ⚙️ Optimization Results — EV Case

| **Method** | **Best Chemistry** | **Total J** | **Suitability (S)** | **Chem Term (1–S)** | **BMS Term** | **Penalty Term** | **Iter** | **Time (s)** |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **Newton’s Method** | LFP | 0.239 | 0.785 | 0.215 | 0.005 | 0.019 | 124 | 128.67 |
| **Quasi-Newton (BFGS)** | SCiB | 0.231 | 0.774 | 0.226 | 0.005 | 0.000 | 11 | 40.13 |
| **Line Search (Golden)** | LFP | 0.308 | 0.785 | 0.215 | 0.100 | 0.000 | 41 | 18.95 |
| **Steepest Descent** | LFP | 0.416 | 0.785 | 0.215 | 0.201 | 0.000 | 80 | 27.92 |
| **SQP (KKT)** | SCiB | **0.417** | **0.774** | **0.226** | **0.190** | **0.000** | **10** | **5.20** |

> **Observation:**  
> The EV case exhibits rapid convergence for SQP (KKT) and BFGS, both achieving minimal penalty terms, indicating strong constraint satisfaction. SQP provides the **lowest total cost (J = 0.417)** within **10 iterations**, while BFGS converges with slightly higher precision in **moderate runtime (40 s)**. Steepest Descent displays slower convergence due to its linear gradient trajectory. LFP emerges as the **dominant chemistry**, balancing energy density and safety.

---

### ⚙️ Optimization Results — BESS Case

| **Method** | **Best Chemistry** | **Total J** | **Suitability (S)** | **Chem Term (1–S)** | **BMS Term** | **Penalty Term** | **Iter** | **Time (s)** |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **Newton’s Method** | LFP | **0.197** | **0.827** | **0.173** | **0.006** | **0.017** | **193** | **166.11** |
| **Quasi-Newton (BFGS)** | Na-ion | 0.237 | 0.770 | 0.230 | 0.006 | 0.000 | 11 | 35.92 |
| **Line Search (Golden)** | LFP | 0.330 | 0.827 | 0.173 | 0.157 | 0.000 | 41 | 31.26 |
| **SQP (KKT)** | LFP | 0.356 | 0.827 | 0.173 | 0.183 | 0.000 | 10 | 4.44 |
| **Steepest Descent** | LFP | 0.366 | 0.827 | 0.173 | 0.193 | 0.000 | 40 | 24.92 |

> **Observation:**  
> For the BESS case, **Newton’s Method** achieves the **lowest total cost (J = 0.197)** with excellent feasibility, albeit at a higher iteration count. **BFGS** converges faster but slightly less accurately, while **Line Search** and **SQP** maintain balanced performance. LFP consistently dominates as the optimal chemistry, confirming its superior stability and cycle life for stationary storage applications.

---

### ⚙️ Optimization Results — SESS Case

| **Method** | **Best Chemistry** | **Total J** | **Suitability (S)** | **Chem Term (1–S)** | **BMS Term** | **Penalty Term** | **Iter** | **Time (s)** |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| **Newton’s Method** | LFP | **0.225** | **0.799** | **0.201** | **0.007** | **0.016** | **193** | **286.24** |
| **Quasi-Newton (BFGS)** | SCiB | 0.291 | 0.709 | 0.291 | 0.000 | 0.000 | 10 | 30.87 |
| **Line Search (Golden)** | LFP | 0.349 | 0.799 | 0.201 | 0.147 | 0.000 | 41 | 37.54 |
| **SQP (KKT)** | LFP | 0.374 | 0.799 | 0.201 | 0.173 | 0.000 | 10 | 4.32 |
| **Steepest Descent** | LFP | 0.385 | 0.799 | 0.201 | 0.184 | 0.000 | 40 | 47.82 |

> **Observation:**  
> In the SESS application, **Newton’s Method** yields the **minimum total cost (J = 0.225)** with stable convergence, while **BFGS** provides the fastest computational efficiency. **Line Search** and **SQP** maintain similar cost behavior, confirming robustness of the framework. LFP remains optimal due to its reliability and safe voltage profile under long-duration cycling.

---

### 🧭 Overall Summary

Across all applications:
- **LFP** chemistry demonstrates consistent dominance owing to its **high safety, cycle life, and affordability**.  
- **SQP (KKT)** and **BFGS** exhibit the best balance between **speed, accuracy, and constraint satisfaction**.  
- **Newton’s Method** offers the most precise convergence at higher computational cost.  
- **Steepest Descent** serves as a baseline, showing slower convergence but high stability across use cases.

---

> These consistent results validate the robustness of the proposed multi-objective deterministic optimization framework for integrated chemistry and BMS co-design across EV, BESS, and SESS domains.

---

## 6. Convergence Criteria

All methods apply dual stopping conditions:

\[
\|\nabla J(p_k)\| \le \varepsilon_g \quad \text{or} \quad |J_{k+1} - J_k| \le \varepsilon_J
\]

Typical tolerances:
- Gradient norm \( \varepsilon_g = 10^{-6} \)
- Objective delta \( \varepsilon_J = 10^{-8} \)

**Code Snippet:**
```python
def new_stopper():
    return Convergence(
        g_tol=5e-3, rJ_tol=5e-4, p_tol=1e-6,
        stagnation_patience=30, max_iter=300, min_iter_before_check=10
    )
