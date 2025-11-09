# optimizer.py
# ============================================================
# Data loading + objective with per-app bounds/weights and
# normalized, unit-safe coupling penalties.
# Now: chemistry suitability S is converted to a cost (1 - S).
# ============================================================

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple
import numpy as np
import pandas as pd

@dataclass
class DataRepo:
    P: pd.DataFrame                 # (n_chem x m_metrics) chemistry metrics, higher better
    alpha: pd.DataFrame             # columns: metric, EV, BESS, SESS
    beta: pd.DataFrame              # either (var, weight) OR wide (var, EV, BESS, SESS)
    gamma: pd.DataFrame             # columns: penalty, weight
    bms_bounds: pd.DataFrame        # default bounds (var, lo, hi)
    envelope: pd.DataFrame          # per-chem recommended envelope & C-rates
    bms_bounds_map: Dict[str, pd.DataFrame]  # optional per-app bounds

def load_repo(root: Path) -> DataRepo:
    sd = root / "sample_data"
    P = pd.read_csv(sd / "chemistry_metrics.csv", index_col=0)
    alpha = pd.read_csv(sd / "alpha_weights.csv")
    beta = pd.read_csv(sd / "beta_weights.csv")
    gamma = pd.read_csv(sd / "gamma_weights.csv")

    # default (fallback) bounds
    bms_bounds_path = sd / "bms_bounds.csv"
    if bms_bounds_path.exists():
        bms_bounds = pd.read_csv(bms_bounds_path)
    else:
        bms_bounds = pd.DataFrame(columns=["var", "lo", "hi"])

    # per-app bounds (preferred)
    bms_bounds_map: Dict[str, pd.DataFrame] = {}
    for app, fname in [("EV", "bms_bounds_EV.csv"),
                       ("BESS", "bms_bounds_BESS.csv"),
                       ("SESS", "bms_bounds_SESS.csv")]:
        f = sd / fname
        if f.exists():
            bms_bounds_map[app] = pd.read_csv(f)

    envelope = pd.read_csv(sd / "chemistry_envelope.csv")
    return DataRepo(P, alpha, beta, gamma, bms_bounds, envelope, bms_bounds_map)

# -----------------------------
# Helpers
# -----------------------------

def normalize01(x, lo, hi):
    return (x - lo) / np.maximum(1e-9, (hi - lo))

def _bounds_for_app(repo: DataRepo, app: str) -> pd.DataFrame:
    return repo.bms_bounds_map.get(app, repo.bms_bounds)

def _beta_for_app(repo: DataRepo, app: str) -> np.ndarray:
    beta = repo.beta.set_index("var")
    if app in beta.columns:
        return beta[app].to_numpy()
    if "weight" in beta.columns:
        return beta["weight"].to_numpy()
    raise ValueError("beta_weights.csv must have either columns [var, weight] or [var, EV, BESS, SESS].")

# -----------------------------
# Objective parts
# -----------------------------

def chem_term(P: pd.DataFrame, alpha_vec: np.ndarray, x_onehot: np.ndarray) -> float:
    """
    Suitability S = alpha^T (P^T x), higher is better (benefit).
    x_onehot selects the chemistry row from P.
    """
    i = int(np.argmax(x_onehot))
    return float(alpha_vec @ P.iloc[i, :].to_numpy())

def G_of_p(p_norm: np.ndarray) -> np.ndarray:
    # BMS cost mapping; identity keeps "larger normalized values cost more"
    return p_norm

def coupling_penalties(x_onehot: np.ndarray, p: np.ndarray, repo: DataRepo) -> Dict[str, float]:
    """Dimensionless, normalized penalties (each ~O(1) at tolerance)."""
    i = int(np.argmax(x_onehot))
    chem = repo.envelope.iloc[i]

    # p order: [Vchg, Vdis_min, Ichg, Idchg, Tmin, Tmax, Tcool, SOCmin, SOCmax, dVbal(mV)]
    Vchg, Vdis_min, Ichg, Idchg, Tmin, Tmax, Tcool, SOCmin, SOCmax, dVbal_mV = p

    # env ΔV may be in V — convert to mV when ≤1.0
    dVmax_env = chem.get("DeltaV_bal_max", np.nan)
    if pd.notna(dVmax_env):
        dVmax_mV = dVmax_env * 1000.0 if dVmax_env <= 1.0 else dVmax_env
    else:
        dVmax_mV = np.inf
    Vdis_min_max = chem.get("V_dis_min_max", np.inf)

    tol = {
        "C_rate":      0.5,   # ±0.5 C
        "Temp":        5.0,   # ±5 °C
        "SOC_window": 10.0,   # ±10 %
        "V":           0.10,  # ±0.10 V
        "dVbal_mV":    5.0,   # ±5 mV
    }

    pen = {}
    pen["C_rate_chg"]  = (max(0.0, float(Ichg)  - float(chem["R_chg"]))  / tol["C_rate"])**2
    pen["C_rate_dchg"] = (max(0.0, float(Idchg) - float(chem["R_dchg"])) / tol["C_rate"])**2
    pen["Temp_low"]    = (max(0.0, float(chem["T_lo"]) - float(Tmin)) / tol["Temp"])**2
    pen["Temp_high"]   = (max(0.0, float(Tmax) - float(chem["T_hi"]))  / tol["Temp"])**2
    env_window = float(chem["SOC_hi"]) - float(chem["SOC_lo"])
    p_window   = float(SOCmax) - float(SOCmin)
    pen["SOC_window"] = (max(0.0, p_window - env_window) / tol["SOC_window"])**2
    pen["V_dis_min"]  = (max(0.0, float(Vdis_min) - float(Vdis_min_max)) / tol["V"])**2 if np.isfinite(Vdis_min_max) else 0.0
    pen["DeltaV_bal"] = (max(0.0, float(dVbal_mV) - float(dVmax_mV)) / tol["dVbal_mV"])**2 if np.isfinite(dVmax_mV) else 0.0

    # clip extreme outliers to keep scale stable
    for k in pen:
        if pen[k] > 25.0:
            pen[k] = 25.0
    return pen

def J_objective(app: str, x_onehot: np.ndarray, p: np.ndarray, repo: DataRepo) -> Tuple[float, Dict[str, float]]:
    """
    Composite cost to MINIMIZE:
        J_total = J_chem + J_bms + J_pen
    with:
        Suitability S = alpha^T (P^T x)   (benefit, higher is better)
        Chemistry cost J_chem = 1 - S     (lower is better)
        J_bms = beta^T G(p_norm)          (lower is better)
        J_pen = sum(gamma[k] * penalty_k) (lower is better; zero if feasible)
    """
    alpha_vec = repo.alpha.set_index("metric")[app].to_numpy()
    beta_w    = _beta_for_app(repo, app)
    gamma_w   = repo.gamma.set_index("penalty")["weight"].to_dict()

    # ---- Suitability (benefit) and Chemistry cost (1 - S)
    S = chem_term(repo.P, alpha_vec, x_onehot)     # suitability S ∈ [0,1], higher = better
    J_chem = 1.0 - S                               # convert to cost

    # ---- BMS term (cost)
    btab = _bounds_for_app(repo, app)
    lo = btab["lo"].to_numpy(); hi = btab["hi"].to_numpy()
    p_norm = normalize01(p, lo, hi)
    J_bms = float(beta_w @ G_of_p(p_norm))

    # ---- Coupling penalties (cost)
    pens = coupling_penalties(x_onehot, p, repo)
    J_pen = sum(gamma_w.get(name, 0.0) * val for name, val in pens.items())

    # ---- Total cost
    J_total = J_chem + J_bms + J_pen

    parts = {
        "suitability": S,   # benefit (S↑ good)
        "chem_term": J_chem,
        "bms_term": J_bms,
        "pen_term": J_pen,
    }
    parts |= {f"pen_{k}": v for k, v in pens.items()}
    return J_total, parts

# --------- small derivative-free tuner (for quick checks; unused by 5-method runner) ---------

def tune_bms_for_chem(app: str, chem_index: int, repo: DataRepo, steps: int = 200, seed: int = 0):
    btab = _bounds_for_app(repo, app)
    lo = btab["lo"].to_numpy(); hi = btab["hi"].to_numpy()
    p = (lo + hi) / 2.0
    x = np.zeros(len(repo.P)); x[chem_index] = 1.0

    best_p = p.copy()
    best_J, best_parts = J_objective(app, x, p, repo)

    rng = np.random.default_rng(seed)
    for _ in range(steps):
        k = int(rng.integers(0, len(p)))
        span = hi[k] - lo[k]
        delta = (rng.random() - 0.5) * 0.2 * span
        cand = p.copy(); cand[k] = np.clip(cand[k] + delta, lo[k], hi[k])
        J_new, parts = J_objective(app, x, cand, repo)
        if J_new < best_J:
            best_J, best_parts = J_new, parts
            best_p = cand
            p = cand
    return best_p, {"J": best_J, **best_parts}
