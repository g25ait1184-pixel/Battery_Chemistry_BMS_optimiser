# run_all_methods.py
# ============================================================
# 5 classical optimizers × 3 applications with per-app bounds,
# robust convergence, and anti-stagnation safeguards.
# Produces:
#   outputs/summary_methods.csv
#   outputs/trace_<APP>_<METHOD>_<CHEM>.csv
# ============================================================

from pathlib import Path
import time
from copy import deepcopy
import numpy as np
import pandas as pd

from optimizer import load_repo, J_objective

# ---------- project & base repo ----------
root = Path(__file__).resolve().parent
_base_repo = load_repo(root)
CHEMS = _base_repo.P.index.tolist()

# ============================================================
# SD behavior knobs (make SD slower/longer; others unchanged)
# ============================================================
SD_ARMIJO_C      = 1e-2      # stricter sufficient-decrease
SD_ARMIJO_TAU    = 0.50      # heavier backtracking
SD_ARMIJO_A0     = 0.20      # smaller initial step
SD_BACKTRACKS    = 25
SD_FALLBACK_STEP = 0.002     # fallback step when line-search fails (was 0.01)
SD_EPS_GRAD      = 1e-5      # sharper finite-diff gradient (was 1e-4)
SD_STOPPER = dict(           # tighter termination for SD only
    g_tol=1e-4,
    rJ_tol=1e-5,
    p_tol=1e-8,
    stagnation_patience=60,
    max_iter=500,
    min_iter_before_check=40,
)

# ============================================================
# Per-app repo & globals
# ============================================================
def make_repo_for_app(base_repo, app: str):
    repo = deepcopy(base_repo)

    # per-app bounds if present
    per_app = {
        "EV":   root / "sample_data" / "bms_bounds_EV.csv",
        "BESS": root / "sample_data" / "bms_bounds_BESS.csv",
        "SESS": root / "sample_data" / "bms_bounds_SESS.csv",
    }
    bp = per_app.get(app)
    if bp and bp.exists():
        repo.bms_bounds = pd.read_csv(bp)

    # narrow beta to this app if wide
    if app in base_repo.beta.columns:
        repo.beta = base_repo.beta[["var", app]].rename(columns={app: "weight"}).copy()

    return repo

# active repo + cached bounds for speed
repo = None
b_lo = None
b_hi = None

def set_active_repo(_repo):
    global repo, b_lo, b_hi
    repo = _repo
    b_lo = repo.bms_bounds["lo"].to_numpy()
    b_hi = repo.bms_bounds["hi"].to_numpy()

# ============================================================
# Convergence controller (smarter stagnation)
# ============================================================
class Convergence:
    def __init__(self,
                 g_tol=5e-3, rJ_tol=5e-4, p_tol=1e-6,
                 stagnation_patience=30, max_iter=300, min_iter_before_check=10):
        self.g_tol = g_tol
        self.rJ_tol = rJ_tol
        self.p_tol  = p_tol
        self.stag_pat = stagnation_patience
        self.max_iter = max_iter
        self.min_chk  = min_iter_before_check
        self.k = 0
        self.prev_J = None
        self.prev_p = None
        self.stag = 0
        self.reason = None

    def update(self, J, p, g=None, step=None):
        self.k += 1
        if self.k >= self.max_iter:
            self.reason = "max_iter"; return True, self.reason

        if self.k < self.min_chk:
            self.prev_J = J; self.prev_p = p.copy(); return False, None

        J_prev = self.prev_J if self.prev_J is not None else J
        denom = max(1.0, abs(J), abs(J_prev))
        rel = abs(J - J_prev) / denom
        self.prev_J = J

        # true convergence checks
        if rel <= self.rJ_tol:
            self.reason = f"rel_delta_J<={self.rJ_tol:g}"; return True, self.reason
        if g is not None and float(np.linalg.norm(g)) <= self.g_tol:
            self.reason = f"grad_norm<={self.g_tol:g}"; return True, self.reason
        if step is not None and float(np.linalg.norm(step)) <= self.p_tol:
            self.reason = f"step_norm<={self.p_tol:g}"; return True, self.reason

        # stagnation only if BOTH improvement and step are tiny
        small_step = (step is not None) and (float(np.linalg.norm(step)) <= 5*self.p_tol)
        if rel <= 0.5*self.rJ_tol and small_step:
            self.stag += 1
        else:
            self.stag = 0

        if self.stag >= self.stag_pat:
            self.reason = f"stagnation({self.stag_pat})"; return True, self.reason

        self.prev_p = p.copy()
        return False, None

def new_stopper():
    # default stopper for non-SD methods
    return Convergence(
        g_tol=5e-3, rJ_tol=5e-4, p_tol=1e-6,
        stagnation_patience=30, max_iter=300, min_iter_before_check=10
    )

def sd_stopper():
    # stricter stopper for SD (more iterations)
    return Convergence(**SD_STOPPER)

# ============================================================
# Utilities
# ============================================================
def proj(p):
    return np.minimum(np.maximum(p, b_lo), b_hi)

def J_only(app, chem_idx, p):
    x = np.zeros(len(CHEMS)); x[chem_idx] = 1.0
    val, _ = J_objective(app, x, proj(p), repo)
    return val

def num_grad_forward(app, chem_idx, p, eps=1e-4):
    g = np.zeros_like(p, dtype=float)
    J0 = J_only(app, chem_idx, p)
    for k in range(len(p)):
        d = np.zeros_like(p); d[k] = eps
        g[k] = (J_only(app, chem_idx, proj(p + d)) - J0) / eps
    return g

def diag_hess_approx(app, chem_idx, p, eps=1e-3):
    Hdiag = np.zeros_like(p, dtype=float)
    J0 = J_only(app, chem_idx, p)
    for k in range(len(p)):
        e = np.zeros_like(p); e[k] = eps
        J1 = J_only(app, chem_idx, proj(p + e))
        J2 = J_only(app, chem_idx, proj(p + 2 * e))
        Hdiag[k] = (J2 - 2 * J1 + J0) / (eps ** 2)
        if Hdiag[k] <= 0: Hdiag[k] = 1e-3
    return Hdiag

# ============================================================
# Line searches
# ============================================================
def armijo_backtracking(J_at, grad_at, p, d,
                        c=1e-3, tau=0.80, alpha0=0.50, max_backtracks=25):
    J0 = J_at(p); g0 = grad_at(p)
    a = alpha0
    for _ in range(max_backtracks):
        if J_at(p + a * d) <= J0 + c * a * float(np.dot(g0, d)):
            return a
        a *= tau
    return 0.0

def armijo_backtracking_sd(J_at, grad_at, p, d,
                           c=SD_ARMIJO_C, tau=SD_ARMIJO_TAU,
                           alpha0=SD_ARMIJO_A0, max_backtracks=SD_BACKTRACKS):
    J0 = J_at(p); g0 = grad_at(p)
    a = alpha0
    for _ in range(max_backtracks):
        if J_at(p + a * d) <= J0 + c * a * float(np.dot(g0, d)):
            return a
        a *= tau
    return 0.0

def golden_section(J_at, p, d, a=0.0, b=1.0, tol=1e-4, max_iters=60):
    phi = (1.0 + np.sqrt(5.0)) / 2.0; rho = phi - 1.0
    x1 = b - rho * (b - a); x2 = a + rho * (b - a)
    f1 = J_at(p + x1 * d); f2 = J_at(p + x2 * d)
    it = 0
    while (b - a) > tol and it < max_iters:
        if f1 <= f2:
            b, x2, f2 = x2, x1, f1; x1 = b - rho*(b - a); f1 = J_at(p + x1*d)
        else:
            a, x1, f1 = x1, x2, f2; x2 = a + rho*(b - a); f2 = J_at(p + x2*d)
        it += 1
    return 0.5 * (a + b)

# ============================================================
# Optimizers
# ============================================================
def opt_sd(app, chem_idx, p0, stopper: Convergence):
    p = p0.copy(); trace = []
    for _ in range(stopper.max_iter):
        # sharper finite-diff gradient ONLY for SD
        g = num_grad_forward(app, chem_idx, p, eps=SD_EPS_GRAD)
        d = -g
        # SD-specific Armijo
        a = armijo_backtracking_sd(lambda q: J_only(app, chem_idx, q),
                                   lambda q: num_grad_forward(app, chem_idx, q, eps=SD_EPS_GRAD),
                                   p, d)
        # smaller fallback step for SD
        s = a * d if a > 0 else SD_FALLBACK_STEP * d
        pn = proj(p + s); Jt = J_only(app, chem_idx, pn)
        trace.append((len(trace)+1, Jt))
        done, _ = stopper.update(Jt, pn, g=g, step=s)
        p = pn
        if done: break
    return p, trace, stopper.reason or "done"

def opt_bfgs(app, chem_idx, p0, stopper: Convergence):
    p = p0.copy(); n = len(p); H = np.eye(n); trace = []
    g = num_grad_forward(app, chem_idx, p)
    for _ in range(stopper.max_iter):
        d = -H @ g
        a = armijo_backtracking(lambda q: J_only(app, chem_idx, q),
                                lambda q: num_grad_forward(app, chem_idx, q), p, d)
        s = a * d if a > 0 else 0.01 * d
        pn = proj(p + s); g_new = num_grad_forward(app, chem_idx, pn)
        y = g_new - g
        # Powell damping to keep H positive definite
        ys = float(y @ s)
        if ys <= 1e-12:
            theta = 0.8
            y = theta*y + (1-theta)*(H @ s)
            ys = float(y @ s)
        if ys > 1e-12:
            rho = 1.0 / ys
            I = np.eye(n)
            H = (I - rho*np.outer(s, y)) @ H @ (I - rho*np.outer(y, s)) + rho*np.outer(s, s)
        p, g = pn, g_new
        Jt = J_only(app, chem_idx, p); trace.append((len(trace)+1, Jt))
        done, _ = stopper.update(Jt, p, g=g, step=s)
        if done: break
    return p, trace, stopper.reason or "done"

def opt_newton_diag(app, chem_idx, p0, stopper: Convergence):
    p = p0.copy(); trace = []
    for _ in range(stopper.max_iter):
        g = num_grad_forward(app, chem_idx, p); Hdiag = diag_hess_approx(app, chem_idx, p)
        d = -(g / Hdiag)
        a = armijo_backtracking(lambda q: J_only(app, chem_idx, q),
                                lambda q: num_grad_forward(app, chem_idx, q), p, d)
        s = a * d if a > 0 else 0.01 * d
        pn = proj(p + s); Jt = J_only(app, chem_idx, pn)
        trace.append((len(trace)+1, Jt))
        done, _ = stopper.update(Jt, pn, g=g, step=s)
        p = pn
        if done: break
    return p, trace, stopper.reason or "done"

def opt_linesearch_golden(app, chem_idx, p0, stopper: Convergence):
    p = p0.copy(); trace = []
    for _ in range(stopper.max_iter):
        g = num_grad_forward(app, chem_idx, p); d = -g
        a = golden_section(lambda q: J_only(app, chem_idx, q), p, d)
        s = a * d; pn = proj(p + s); Jt = J_only(app, chem_idx, pn)
        trace.append((len(trace)+1, Jt))
        done, _ = stopper.update(Jt, pn, g=g, step=s)
        p = pn
        if done: break
    return p, trace, stopper.reason or "done"

# ---- KKT-SQP (active-set on linearized inequalities) ----
def _unpack_p(p):
    Vchg, Vdis_min, Ichg, Idchg, Tmin, Tmax, Tcool, SOCmin, SOCmax, dVbal = p
    return Vchg, Vdis_min, Ichg, Idchg, Tmin, Tmax, Tcool, SOCmin, SOCmax, dVbal

def constraints_and_jac(chem_row, p):
    Vchg, Vdis_min, Ichg, Idchg, Tmin, Tmax, Tcool, SOCmin, SOCmax, dVbal = _unpack_p(p)
    Rchg   = float(chem_row["R_chg"]); Rdchg  = float(chem_row["R_dchg"])
    Tlo    = float(chem_row["T_lo"]);  Thi    = float(chem_row["T_hi"])
    SOCl   = float(chem_row["SOC_lo"]); SOCh  = float(chem_row["SOC_hi"])
    VdisMx = float(chem_row.get("V_dis_min_max", np.inf))
    dVmax  = float(chem_row.get("DeltaV_bal_max", np.inf))

    c = np.array([
        Ichg - Rchg,
        Idchg - Rdchg,
        Tlo - Tmin,
        Tmax - Thi,
        (SOCmax - SOCmin) - (SOCh - SOCl),
        Vdis_min - VdisMx,
        dVbal - (dVmax if dVmax > 1 else dVmax*1000.0)
    ], float)

    Jc = np.zeros((7,10), float)
    Jc[0,2] = 1.0; Jc[1,3] = 1.0; Jc[2,4] = -1.0; Jc[3,5] = 1.0
    Jc[4,8] = 1.0; Jc[4,7] = -1.0; Jc[5,1] = 1.0; Jc[6,9] = 1.0
    return c, Jc

def solve_kkt_eq_qp(B, g, A, b):
    n = B.shape[0]; m = A.shape[0]
    K = np.zeros((n+m, n+m), float)
    K[:n,:n] = B + 1e-6*np.eye(n); K[:n,n:] = A.T; K[n:,:n] = A
    rhs = np.concatenate([-g, b])
    sol = np.linalg.lstsq(K, rhs, rcond=None)[0]
    return sol[:n], sol[n:]

def opt_sqp_kkt(app, chem_idx, p0, stopper: Convergence):
    chem_row = repo.envelope.iloc[chem_idx]
    p = p0.copy(); n = len(p); H = np.eye(n); trace = []

    def J_at(pp):
        x = np.zeros(len(CHEMS)); x[chem_idx] = 1.0
        val, _ = J_objective(app, x, proj(pp), repo)
        return val
    def grad_at(pp): return num_grad_forward(app, chem_idx, proj(pp), eps=1e-4)

    g = grad_at(p)
    for _ in range(stopper.max_iter):
        c, Jc = constraints_and_jac(chem_row, p)
        active = np.where(c > -1e-8)[0]
        if active.size > 0:
            Aeq = Jc[active, :]; beq = -c[active]
            d, _ = solve_kkt_eq_qp(H, g, Aeq, beq)
        else:
            d = -np.linalg.solve(H + 1e-8*np.eye(n), g)

        alpha = 1.0
        for _ in range(20):
            pn = proj(p + alpha*d)
            c_new, _ = constraints_and_jac(chem_row, pn)
            if np.all(c_new <= 0.0 + 1e-10) and J_at(pn) <= J_at(p) + 1e-4*alpha*np.dot(g, d):
                break
            alpha *= 0.5

        s = alpha * d; p = proj(p + s)
        g_new = grad_at(p); y = g_new - g; ys = float(y @ s)
        if ys <= 1e-12:
            theta = 0.8
            y = theta*y + (1-theta)*(H @ s)
            ys = float(y @ s)
        if ys > 1e-12:
            I = np.eye(n)
            H = (I - np.outer(s,y)/ys) @ H @ (I - np.outer(y,s)/ys) + np.outer(s,s)/ys
        g = g_new

        Jt = J_at(p); trace.append((len(trace)+1, Jt))
        done, _ = stopper.update(Jt, p, g=g, step=s)
        if done: break
    return p, trace, stopper.reason or "done"

# ============================================================
# Soft restarts (gentle, deterministic)
# ============================================================
_rng = np.random.default_rng(0)

def run_with_soft_restarts(mfunc, app, chem_idx, p0, lo, hi, proj_fn, J_fn,
                           max_restarts=2, sigma=0.01):
    p_star, trace, reason = mfunc(app, chem_idx, p0)
    best_J = J_fn(app, chem_idx, p_star)
    best = (best_J, p_star, trace, reason)

    if not str(reason).startswith("stagnation"):
        return best

    span = (hi - lo)
    for _ in range(max_restarts):
        p0_new = proj_fn(p0 + sigma * span * _rng.normal(size=p0.shape))
        p2, tr2, r2 = mfunc(app, chem_idx, p0_new)
        J2 = J_fn(app, chem_idx, p2)
        if J2 < best[0]:
            best = (J2, p2, tr2, r2)
        if not str(r2).startswith("stagnation"):
            break
    return best

# ============================================================
# Registry + driver
# ============================================================
def METHODS_factory():
    return {
        "Steepest Descent":     lambda app, idx, p0: opt_sd(app, idx, p0, sd_stopper()),
        "Quasi-Newton (BFGS)":  lambda app, idx, p0: opt_bfgs(app, idx, p0, new_stopper()),
        "Newton's Method":      lambda app, idx, p0: opt_newton_diag(app, idx, p0, new_stopper()),
        "Line Search (Golden)": lambda app, idx, p0: opt_linesearch_golden(app, idx, p0, new_stopper()),
        "SQP (KKT)":            lambda app, idx, p0: opt_sqp_kkt(app, idx, p0, new_stopper()),
    }

def run_all():
    METHODS = METHODS_factory()
    summaries = []; traces_to_save = []
    print(f"[INFO] Enabled methods: {list(METHODS.keys())}")

    for app in ["EV", "BESS", "SESS"]:
        repo_app = make_repo_for_app(_base_repo, app)
        set_active_repo(repo_app)

        # base start at mid-bounds
        p0_base_nominal = (b_lo + b_hi) / 2.0

        for mname, mfunc in METHODS.items():
            # SD: start slightly off mid-bound to lengthen the path
            if mname == "Steepest Descent":
                span = (b_hi - b_lo)
                p0_base = proj(p0_base_nominal + 0.02 * span)
            else:
                p0_base = p0_base_nominal

            t0 = time.time()
            best = {"J": float("inf")}
            for i, chem in enumerate(CHEMS):
                # Disable soft-restarts for SD to avoid shortcuts
                restarts = 0 if mname == "Steepest Descent" else 2
                J_star, p_star, trace, reason = run_with_soft_restarts(
                    mfunc, app, i, p0_base, b_lo, b_hi, proj, J_only,
                    max_restarts=restarts, sigma=0.01
                )
                if J_star < best["J"]:
                    best = {"chem": chem, "p": p_star, "J": J_star,
                            "iters": len(trace), "trace": trace, "reason": reason,
                            "chem_idx": i}
            elapsed = time.time() - t0

            # parts (including Suitability) for the best
            x = np.zeros(len(CHEMS)); x[best["chem_idx"]] = 1.0
            J_total, parts = J_objective(app, x, proj(best["p"]), repo)

            summaries.append({
                "Application": app,
                "Method": mname,
                "Best Chemistry": best["chem"],
                "Total J": float(J_total),                 # cost (lower is better)
                "Suitability (S)": float(parts.get("suitability", 1.0 - parts["chem_term"])),
                "Chem Cost (1-S)": float(parts["chem_term"]),
                "BMS Term": float(parts["bms_term"]),
                "Penalty Term": float(parts["pen_term"]),
                "Iter": int(best["iters"]),
                "Converged": True,
                "Stop Reason": best["reason"],
                "Time (s)": elapsed,
            })

            df_trace = pd.DataFrame(best["trace"], columns=["Iter", "J"])
            traces_to_save.append((app, mname, best["chem"], df_trace))

    return pd.DataFrame(summaries), traces_to_save

# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    summary_df, traces = run_all()
    outdir = root / "outputs"; outdir.mkdir(exist_ok=True)

    for app, mname, chem, df in traces:
        fn = outdir / f"trace_{app}_{mname.replace(' ', '_').replace('(', '').replace(')', '').replace('-', '')}_{chem}.csv"
        df.to_csv(fn, index=False)

    summary_sorted = summary_df.sort_values(by=["Application", "Total J"])
    summary_sorted.to_csv(outdir / "summary_methods.csv", index=False, encoding="utf-8")
    pd.set_option("display.max_columns", None)
    print(summary_sorted.to_string(index=False))
