"""
alpha_theory.py
Computes the theoretical thermal-expansion coefficient α for each LJ glass
inherent state using the microscopic formula (eq. 24 of the PDF):

    α = 1/(2·K·V·d̄) × [ −∂³U/∂η∂x² : H⁻¹  +  ∂²U/∂η∂x · H⁻¹ · U‴ : H⁻¹ ]

Expanded in the normal-mode basis (H⁻¹ = Σ_ω 1/ω² Ψω⊗Ψω):

    α = −1/(2KVd̄) × Σ_ω (1/ω²) [ u_eta_x_x(Ψω,Ψω) + tess3(vna_file,Ψω,Ψω) ]

SIGN CONVENTION:
    compression_nonaffine_velocities() in C returns vna_eta_file = −H⁻¹·ξ_η
    (sign flip at line 836 of the C file).  The anharmonic term becomes:
        ξ_η · H⁻¹ · tess = (−vna_file) · Σ_ω 1/ω² tess3(Ψω,Ψω)
    giving a combined −1/(2KVd̄) prefactor for BOTH terms.

INPUTS  (written by 3dlj_elasticity_analysis.c):
    hessian_g{g}.txt     sparse COO Hessian        [write_for_python]
    vna_eta_g{g}.txt     nonaffine compression vels [write_for_python]
    bonds_g{g}.txt       per-bond geometry+derivs   [write_bond_network]
    meta_g{g}.txt        scalar globals              [write_bond_network]
    elasticity_results.txt  G, K per glass           [main]

CACHING:
    Expensive diagonalisation (scipy.linalg.eigh on 12000x12000) is done once
    per glass and saved to  alpha_cache/alpha_cache_g{g}.npz.
    All subsequent runs skip diagonalisation entirely.

OUTPUT:
    alpha_results.txt        alpha_theory + harmonic/anharmonic split per glass
    alpha_comparison.png     bar chart of contributions + scatter vs numerical
"""

import os
import sys
import argparse
import numpy as np
import scipy.sparse as sp
import scipy.linalg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def _parse_args():
    p = argparse.ArgumentParser(description="Compute alpha_theory from C output files")
    p.add_argument("--data-dir",  default=".",           help="directory with C output files (default: .)")
    p.add_argument("--cache-dir", default="alpha_cache", help="directory for .npz caches (default: alpha_cache)")
    p.add_argument("--out-dir",   default=".",           help="directory for results and plot (default: .)")
    p.add_argument("--glasses",   default="0-9",         help="glasses to process e.g. 0-9 or 0,2,5 (default: 0-9)")
    return p.parse_args()


def _parse_glass_range(spec):
    if "-" in spec and "," not in spec:
        lo, hi = spec.split("-")
        return list(range(int(lo), int(hi) + 1))
    return [int(x) for x in spec.split(",")]

# Constants

N_ATOMS          = 4000
DIM              = 3
DOFS             = N_ATOMS * DIM       # 12 000
ZERO_MODE_CUTOFF = 1e-4                # omega^2 threshold for zero-mode removal
BATCH_SIZE       = 64                  # modes processed per numpy batch

# File loaders

def _p(data_dir, fname):
    return os.path.join(data_dir, fname)


def load_hessian(g, data_dir):
    print(f"    loading hessian_g{g}.txt", flush=True)
    raw  = np.loadtxt(_p(data_dir, f"hessian_g{g}.txt"))
    rows = raw[:, 0].astype(np.int32)
    cols = raw[:, 1].astype(np.int32)
    vals = raw[:, 2]
    return sp.csr_matrix((vals, (rows, cols)), shape=(DOFS, DOFS))


def load_vna_eta(g, data_dir):
    print(f"    loading vna_eta_g{g}.txt", flush=True)
    return np.loadtxt(_p(data_dir, f"vna_eta_g{g}.txt")).reshape(N_ATOMS, DIM)


def load_bonds(g, data_dir):
    print(f"    loading bonds_g{g}.txt", flush=True)
    data     = np.loadtxt(_p(data_dir, f"bonds_g{g}.txt"), comments="#")
    ii       = data[:, 0].astype(np.int32)
    jj       = data[:, 1].astype(np.int32)
    rij_arr  = data[:, 2:5]
    first_b  = data[:, 5]
    second_b = data[:, 6]
    third_b  = data[:, 7]
    return ii, jj, rij_arr, first_b, second_b, third_b


def load_meta(g, data_dir):
    meta = {}
    with open(_p(data_dir, f"meta_g{g}.txt")) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            key, val = line.split()
            try:
                meta[key] = int(val)
            except ValueError:
                meta[key] = float(val)
    return meta


def load_elasticity_results(data_dir):
    raw = np.loadtxt(_p(data_dir, "elasticity_results.txt"), comments="#")
    out = {}
    for row in raw:
        g = int(row[0])
        out[g] = {"pressure": row[1], "u_per_N": row[2], "G": row[3], "K": row[4]}
    return out


def load_numerical_alpha(data_dir):
    """Read alpha_numerical.txt written by the Jupyter notebook."""
    fname = _p(data_dir, "alpha_numerical.txt")
    if not os.path.exists(fname):
        return None
    meta = {}
    with open(fname) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                k, v = line.split()
                meta[k] = float(v)
    if "dPdT" in meta and "K" in meta and "alpha" not in meta:
        meta["alpha"] = meta["dPdT"] / meta["K"]
    return meta

# Bond-loop functions — numpy mirrors of the C bond loops
# All loops are batched over BATCH_SIZE modes at once.
# Naming mirrors the C code exactly.
#   Psi_diff : (batch, contacts, DIM)  =  evec[jj] - evec[ii]  per mode

def _bond_diffs(evecs_batch, ii, jj):
    """evecs_batch (batch, N, DIM) -> (batch, contacts, DIM)"""
    return evecs_batch[:, jj, :] - evecs_batch[:, ii, :]


def u_eta_x_x_batch(Psi_diff, rij_arr, r2, second_b, third_b):
    """
    Mirrors C: u_eta_x_x_on_two_vectors(Psi, Psi)

    = sum_bonds  third*r^2*(Psi.r)^2  +  second*(2*(Psi.r)^2 + r^2*|Psi|^2)

    Returns (batch,) array.
    """
    psi_dot_r = np.einsum("bcd,cd->bc", Psi_diff, rij_arr)    # (batch, contacts)
    psi_sq    = np.einsum("bcd,bcd->bc", Psi_diff, Psi_diff)  # (batch, contacts)
    contrib   = (third_b * r2 * psi_dot_r**2
                 + second_b * (2.0 * psi_dot_r**2 + r2 * psi_sq))
    return contrib.sum(axis=1)   # (batch,)


def tess3_batch(Psi_diff, rij_arr, vna_diff, vna_dot_r, r2, second_b, third_b):
    """
    Mirrors C: tessianOnThreeVectors(vna_eta_file, Psi, Psi)

    = sum_bonds  third*(Psi.r)^2*(vna.r)
               + second*(2*(vna.Psi)*(Psi.r) + |Psi|^2*(vna.r))

    vna_diff  : (contacts, DIM)  constant across modes
    vna_dot_r : (contacts,)      constant across modes

    Returns (batch,) array.
    """
    psi_dot_r   = np.einsum("bcd,cd->bc", Psi_diff, rij_arr)   # (batch, contacts)
    psi_sq      = np.einsum("bcd,bcd->bc", Psi_diff, Psi_diff)  # (batch, contacts)
    vna_dot_psi = np.einsum("cd,bcd->bc", vna_diff, Psi_diff)   # (batch, contacts)
    contrib = (third_b * psi_dot_r**2 * vna_dot_r
               + second_b * (2.0 * vna_dot_psi * psi_dot_r
                              + psi_sq * vna_dot_r))
    return contrib.sum(axis=1)   # (batch,)

# Core: diagonalise H and accumulate per-mode sums

def compute_mode_contributions(g, H, vna_eta, ii, jj, rij_arr,
                                first_b, second_b, third_b):
    """
    Diagonalise H then accumulate per-mode sums for both terms in eq. 24.

    Returns
    -------
    valid_evals     (n_valid,)  omega^2 for each physical mode
    term1_per_mode  (n_valid,)  u_eta_x_x(Psi,Psi)            (no 1/omega^2 yet)
    term2_per_mode  (n_valid,)  tess3(vna_file, Psi,Psi)       (no 1/omega^2 yet)
    """
    # 1. Dense diagonalisation via LAPACK DSYEVD
    print(f"    [g={g}] converting sparse H to dense {DOFS}x{DOFS}", flush=True)
    H_dense = H.toarray()
    print(f"    [g={g}] scipy.linalg.eigh",
          flush=True)
    evals, evecs = scipy.linalg.eigh(H_dense)
    del H_dense

    # 2. Strip zero modes (translations)
    valid        = evals > ZERO_MODE_CUTOFF
    valid_evals  = evals[valid]
    valid_evecs  = evecs[:, valid]      # (DOFS, n_valid)
    del evecs
    n_valid = valid_evals.size
    print(f"    [g={g}] {n_valid} physical modes  ({DOFS - n_valid} near-zero removed)",
          flush=True)

    # 3. Precompute vna quantities that are constant across all modes
    vna_diff  = vna_eta[jj] - vna_eta[ii]                       # (contacts, DIM)
    vna_dot_r = np.einsum("cd,cd->c", vna_diff, rij_arr)         # (contacts,)
    r2        = np.einsum("cd,cd->c", rij_arr, rij_arr)          # (contacts,)

    # 4. Batched mode-sum loop
    term1_per_mode = np.zeros(n_valid)
    term2_per_mode = np.zeros(n_valid)
    n_batches      = (n_valid + BATCH_SIZE - 1) // BATCH_SIZE

    for b in range(n_batches):
        if b % 50 == 0:
            print(f"      mode batch {b+1}/{n_batches} ...", flush=True)
        lo, hi = b * BATCH_SIZE, min((b + 1) * BATCH_SIZE, n_valid)
        bsz    = hi - lo

        Psi_mat  = valid_evecs[:, lo:hi].T.reshape(bsz, N_ATOMS, DIM)
        Psi_diff = _bond_diffs(Psi_mat, ii, jj)

        term1_per_mode[lo:hi] = u_eta_x_x_batch(
            Psi_diff, rij_arr, r2, second_b, third_b)
        term2_per_mode[lo:hi] = tess3_batch(
            Psi_diff, rij_arr, vna_diff, vna_dot_r, r2, second_b, third_b)
    print(f"    [g={g}] term1 unweighted sum = {term1_per_mode.sum():+.6e}")
    print(f"    [g={g}] term2 unweighted sum = {term2_per_mode.sum():+.6e}")

    return valid_evals, valid_evecs, term1_per_mode, term2_per_mode


def alpha_from_mode_sums(valid_evals, term1_per_mode, term2_per_mode, K, V):
    """
    Apply eq. 24 from the accumulated per-mode sums.

    alpha = -1/(2KVd) * sum_omega (1/omega^2) [ term1(omega) + term2(omega) ]

    Returns alpha_total, alpha_harmonic, alpha_anharmonic.
    """
    w   = 1.0 / valid_evals
    pre = -1.0 / (2.0 * K * V * DIM)
    ah  = pre * np.dot(w, term1_per_mode)
    aa  = pre * np.dot(w, term2_per_mode)
    return ah + aa, ah, aa

# Per-glass driver with load/save caching

def process_glass(g, K, data_dir, cache_dir):
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = os.path.join(cache_dir, f"alpha_cache_g{g}.npz")

    if os.path.exists(cache_file):
        print(f"  [g={g}] cache hit -- loading {cache_file}", flush=True)
        cached         = np.load(cache_file)
        valid_evals    = cached["valid_evals"]
        valid_evecs    = cached["valid_evecs"]
        term1_per_mode = cached["term1_per_mode"]
        term2_per_mode = cached["term2_per_mode"]
        V              = float(cached["V"])
        print(f"    [g={g}] term1 unweighted sum = {term1_per_mode.sum():+.6e}")   # add
        print(f"    [g={g}] term2 unweighted sum = {term2_per_mode.sum():+.6e}")   # add
    else:
        print(f"  [g={g}] no cache -- running full computation ...", flush=True)
        meta = load_meta(g, data_dir)
        V    = meta["V"]
        H       = load_hessian(g, data_dir)
        vna_eta = load_vna_eta(g, data_dir)
        ii, jj, rij_arr, first_b, second_b, third_b = load_bonds(g, data_dir)
        valid_evals, valid_evecs, term1_per_mode, term2_per_mode = compute_mode_contributions(
            g, H, vna_eta, ii, jj, rij_arr, first_b, second_b, third_b)
        np.savez(cache_file,
                 valid_evals    = valid_evals,
                 valid_evecs    = valid_evecs,
                 term1_per_mode = term1_per_mode,
                 term2_per_mode = term2_per_mode,
                 V              = np.array([V]))
        print(f"  [g={g}] saved cache -> {cache_file}", flush=True)

    alpha, ah, aa = alpha_from_mode_sums(valid_evals, term1_per_mode, term2_per_mode, K, V)
    return alpha, ah, aa, V

# Plot

def make_plot(results, num_alpha, out_dir):
    arr    = np.array(results)
    gs     = arr[:, 0].astype(int)
    alphas = arr[:, 2]
    harms  = arr[:, 3]
    anhs   = arr[:, 4]
    x      = np.arange(len(gs))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: per-glass contribution bar chart
    ax = axes[0]
    ax.bar(x - 0.25, harms,   0.25, label="harmonic term",   color="steelblue",  alpha=0.85)
    ax.bar(x,        anhs,    0.25, label="anharmonic term", color="darkorange", alpha=0.85)
    ax.bar(x + 0.25, alphas,  0.25, label="alpha_theory",    color="darkgreen",  alpha=0.85)
    if num_alpha:
        ax.axhline(num_alpha["alpha"], color="red", ls="--", lw=1.8,
                   label=f"alpha_numerical = {num_alpha['alpha']:.4e}")
    ax.axhline(alphas.mean(), color="darkgreen", ls=":", lw=1.5,
               label=f"alpha_theory mean = {alphas.mean():.4e}")
    ax.set_xticks(x)
    ax.set_xticklabels([f"g{g}" for g in gs], fontsize=9)
    ax.set_xlabel("Glass"); ax.set_ylabel("alpha")
    ax.set_title("alpha_theory per glass  (eq. 24) -- harmonic vs anharmonic")
    ax.legend(fontsize=8); ax.grid(axis="y", alpha=0.3)

    # Right: theory vs numerical
    ax2 = axes[1]
    if num_alpha:
        alpha_num = num_alpha["alpha"]
        ax2.scatter(np.full(len(alphas), alpha_num), alphas,
                    color="navy", s=70, zorder=3, label="individual glasses")
        ax2.errorbar(alpha_num, alphas.mean(), yerr=alphas.std(),
                     fmt="r*", ms=14, capsize=5, zorder=5,
                     label=f"mean +/- std = {alphas.mean():.3e} +/- {alphas.std():.3e}")
        pad  = max(abs(alpha_num), abs(alphas.mean())) * 0.15
        lims = (min(alpha_num, alphas.min()) - pad,
                max(alpha_num, alphas.max()) + pad)
        ax2.plot(lims, lims, "k--", lw=1, label="1:1 line")
        ax2.set_xlim(lims); ax2.set_ylim(lims)
        ax2.set_xlabel("alpha_numerical  (NVT dP/dT / K)", fontsize=11)
        ax2.set_ylabel("alpha_theory  (eq. 24)", fontsize=11)
        ax2.set_title("Theory vs. Numerical")
        ax2.legend(fontsize=9); ax2.grid(alpha=0.3)
    else:
        ax2.hist(alphas, bins=max(3, len(gs)//2),
                 color="steelblue", edgecolor="white")
        ax2.set_xlabel("alpha_theory"); ax2.set_ylabel("Count")
        ax2.set_title("Distribution of alpha_theory across glasses")
        ax2.grid(alpha=0.3)

    fig.suptitle("Thermal Expansion -- alpha_theory (eq. 24)", fontsize=13, y=1.01)
    fig.tight_layout()
    out = os.path.join(out_dir, "alpha_comparison.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved plot -> {out}")
    return out

# Main

def main():
    args    = _parse_args()
    glasses = _parse_glass_range(args.glasses)

    print("=" * 65)
    print("  Alpha Theory (microscopic, inherent-state)")
    print(f"  data-dir  : {args.data_dir}")
    print(f"  cache-dir : {args.cache_dir}")
    print(f"  glasses   : {glasses}")
    print("=" * 65)

    try:
        elas = load_elasticity_results(args.data_dir)
    except FileNotFoundError:
        sys.exit("ERROR: elasticity_results.txt not found in data-dir. "
                 "Run the C binary first.")

    results = []

    for g in glasses:
        if g not in elas:
            print(f"\n[g={g}] not in elasticity_results.txt -- skipping")
            continue
        K = elas[g]["K"]
        print(f"\n[g={g}]  K = {K:.6f}")
        alpha, ah, aa, V = process_glass(g, K, args.data_dir, args.cache_dir)
        print(f"  alpha_theory     = {alpha:+.6e}")
        print(f"    harmonic       = {ah:+.6e}")
        print(f"    anharmonic     = {aa:+.6e}")
        results.append((g, K, alpha, ah, aa))

    if not results:
        sys.exit("No glasses processed -- check data-dir and glass indices.")

    # Write alpha_results.txt
    out_path = os.path.join(args.out_dir, "alpha_results.txt")
    with open(out_path, "w") as f:
        f.write("# glass   K              alpha_theory   "
                "alpha_harm     alpha_anharm\n")
        for row in results:
            f.write(f"{int(row[0]):6d}  {row[1]:14.8g}  {row[2]:+14.8g}  "
                    f"{row[3]:+14.8g}  {row[4]:+14.8g}\n")
    print(f"\nWrote {out_path}")

    arr    = np.array(results)
    alphas = arr[:, 2]
    num    = load_numerical_alpha(args.data_dir)

    print("\n" + "=" * 65)
    print(f"  Mean alpha_theory = {alphas.mean():+.6e}  +/-  {alphas.std():.6e}")
    if num:
        ratio = alphas.mean() / num["alpha"]
        diff  = abs(alphas.mean() - num["alpha"]) / abs(num["alpha"]) * 100
        print(f"  alpha_numerical   = {num['alpha']:+.6e}")
        print(f"  Ratio  theory / numerical = {ratio:.4f}")
        print(f"  Discrepancy               = {diff:.2f} %")
    else:
        print("  (alpha_numerical.txt not found -- run notebook Cell A first)")
    print("=" * 65)

    make_plot(results, num, args.out_dir)


if __name__ == "__main__":
    main()
