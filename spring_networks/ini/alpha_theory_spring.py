"""
alpha_theory_spring.py
======================
Compute the microscopic thermal expansion coefficient alpha_theory for a 3D
spring network, using full Hessian diagonalisation.

Mirrors the structure of alpha_theory_FINAL.py from thermal_expansion_lj.ipynb.

FORMULA (PRE 2010, eq. 25 adapted for spring networks)
-------------------------------------------------------
alpha = -(1 / (2 * K * V * DIM)) * sum_l (1/omega_l^2) * [term1_l + term2_l]

term1_l  = u_eta_xx(Psi_l, Psi_l)
         = sum_bonds |Delta_Psi_l^perp|^2
         = GEOMETRIC anharmonicity: comes from expanding r = |r_ij| to second
           order in displacements.  Non-zero even though d^3U/dr^3 = 0
           for harmonic springs.

term2_l  = tess3(u_na, Psi_l, Psi_l)
         = sum_bonds (1/r) * [2*(r_hat . Delta_Psi_l) * (u_na_perp . Delta_Psi_l)
                               + (u_na . r_hat) * |Delta_Psi_l^perp|^2]
         = ANHARMONIC correction: projection of the non-affine velocity field
           u_na = -H^{-1} xi_eta onto the geometric cubic coupling tensor.
           For stress-free springs (first[l]=0) this involves only the
           transverse projector acting on u_na and the eigenvector difference.

Split reported as:
  alpha_harmonic    from the purely longitudinal part: sum_l omega_l^2 / omega_l^2
                    = N_modes => alpha_harm = -N_modes / (3*K*V)
  alpha_anharmonic  from term1 (transverse geometric) + term2 (u_na correction)
  alpha_total       = alpha_harmonic + alpha_anharmonic

USAGE
-----
  python alpha_theory_spring.py --serial 0 [--data-dir .] [--out-dir .]

INPUT FILES (all produced by spring_networks_elasticity.c)
-----------
  snp_elasticity_N4000_z16.24_s000.dat
  hessian_sn<s>.txt       sparse COO Hessian (row col val)
  bonds_sn<s>.txt         bond geometry (i j rij_x rij_y rij_z first second third)
  vna_eta_sn<s>.txt       compression non-affine velocity (N x 3)

OUTPUT FILE
-----------
  alpha_theory_spring_results.txt
  columns: serial  alpha_harmonic  alpha_anharmonic  alpha_total

SOURCE ATTRIBUTION
------------------
  Formula:        PRE 2010 eq. 25 (Karmakar, Lerner, Procaccia)
                  PNAS 2017 (Zylberg et al.)
  Implementation: mirrors alpha_theory_FINAL.py from thermal_expansion_lj.ipynb
  Hessian:        hessianOnVector() in 3D_relaxed_spring_network_analysis.c
  Bond files:     write_bonds_text(), write_hessian_text(), write_vna_eta_text()
                  in spring_networks_elasticity.c
                  (adapted from 3dlj_elasticity_analysis.c)
"""

import argparse
import numpy as np
import scipy.linalg
import scipy.sparse as sp
from pathlib import Path
import time

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
N_ATOMS  = 4000
DIM      = 3
N_MODES  = DIM * N_ATOMS - DIM     # 11997 physical modes (remove 3 translations)
ZERO_EV  = 1e-5                    # threshold to identify translational modes


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--serial',   type=int, default=0,
                   help='Network serial number (default: 0)')
    p.add_argument('--data-dir', type=Path, default=Path('.'),
                   help='Directory containing input files (default: .)')
    p.add_argument('--out-dir',  type=Path, default=Path('.'),
                   help='Directory for output file (default: .)')
    return p.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# Load K from elasticity file
# ─────────────────────────────────────────────────────────────────────────────
def load_K(data_dir, serial):
    """
    Read bulk modulus K for a given serial from the elasticity output.
    File written by spring_networks_elasticity.c (calculateBulkModulus from
    3D_relaxed_spring_network_analysis.c).
    columns: serial  P0  K_born  K_na  K  mu_born  mu_na  mu
    """
    path = data_dir / "snp_elasticity_N4000_z16.24_s000.dat"
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[np.newaxis, :]
    row = data[data[:, 0].astype(int) == serial]
    if len(row) == 0:
        raise ValueError(f"Serial {serial} not found in {path}")
    P0, K_born, K_na, K = row[0, 1], row[0, 2], row[0, 3], row[0, 4]
    print(f"  K_born={K_born:.5f}  K_na={K_na:.5f}  K={K:.5f}  P0={P0:.2e}")
    return K, P0


# ─────────────────────────────────────────────────────────────────────────────
# Load bonds
# ─────────────────────────────────────────────────────────────────────────────
def load_bonds(data_dir, serial):
    """
    Load bond geometry from bonds_sn<s>.txt.
    Written by write_bonds_text() in spring_networks_elasticity.c
    (adapted from 3dlj_elasticity_analysis.c).
    columns: i  j  rij_x  rij_y  rij_z  first  second  third
    """
    path = data_dir / f"bonds_sn{serial}.txt"
    raw  = np.loadtxt(path, comments="#")
    ii   = raw[:, 0].astype(int)
    jj   = raw[:, 1].astype(int)
    rij  = raw[:, 2:5]                              # physical bond vectors
    r_mag  = np.linalg.norm(rij, axis=1)            # bond lengths
    r_hat  = rij / r_mag[:, np.newaxis]             # unit vectors
    print(f"  {len(ii)} bonds  |  mean|r|={r_mag.mean():.4f}  "
          f"mean_r2={(r_mag**2).mean():.4f}")
    return ii, jj, rij, r_mag, r_hat


# ─────────────────────────────────────────────────────────────────────────────
# Load vna_eta (compression non-affine velocity)
# ─────────────────────────────────────────────────────────────────────────────
def load_vna(data_dir, serial):
    """
    Load compression non-affine velocity from vna_eta_sn<s>.txt.
    Written by write_vna_eta_text() in spring_networks_elasticity.c
    (adapted from 3dlj_elasticity_analysis.c).
    vna = -H^{-1} xi_eta  where  xi_eta[j] = sum_{l in j} rij_l
    Shape: (N_ATOMS, DIM)
    """
    path = data_dir / f"vna_eta_sn{serial}.txt"
    vna  = np.loadtxt(path, comments="#")           # (N, 3)
    assert vna.shape == (N_ATOMS, DIM), \
        f"vna shape {vna.shape} != ({N_ATOMS},{DIM})"
    print(f"  vna_eta loaded: max|vna|={np.linalg.norm(vna,axis=1).max():.4f}")
    return vna


# ─────────────────────────────────────────────────────────────────────────────
# Load Hessian and fully diagonalise
# ─────────────────────────────────────────────────────────────────────────────
def load_and_diagonalise(data_dir, serial):
    """
    Load the sparse COO Hessian (written by write_hessian_text() in
    spring_networks_elasticity.c, adapted from 3dlj_elasticity_analysis.c),
    convert to dense, and fully diagonalise with scipy.linalg.eigh.

    Full diagonalisation avoids the extrapolation needed when using eigsh
    with only N_EIGS modes.  scipy.linalg.eigh uses LAPACK dsyevd
    (divide-and-conquer) — fastest symmetric eigensolver available.

    Returns: ev (N_phys,), evec (sz, N_phys)
    """
    path = data_dir / f"hessian_sn{serial}.txt"
    sz   = DIM * N_ATOMS

    # ── load sparse ──────────────────────────────────────────────────────────
    print(f"  Loading COO triplets from {path.name} "
          f"({path.stat().st_size/1e6:.1f} MB) ...", flush=True)
    t0  = time.time()
    raw = np.loadtxt(path)
    H_sp = sp.coo_matrix(
        (raw[:, 2], (raw[:,0].astype(int), raw[:,1].astype(int))),
        shape=(sz, sz)).tocsr()
    H_sp = (H_sp + H_sp.T) / 2.0          # enforce exact symmetry
    print(f"  Sparse: shape {H_sp.shape}  nnz {H_sp.nnz}  "
          f"({time.time()-t0:.1f}s)")

    # ── convert to dense ─────────────────────────────────────────────────────
    print(f"  Converting to dense ({sz}×{sz}, "
          f"{sz*sz*8/1e9:.2f} GB) ...", flush=True)
    t1     = time.time()
    H_dense = H_sp.toarray()
    print(f"  Dense allocated ({time.time()-t1:.1f}s)")

    # ── full diagonalisation ─────────────────────────────────────────────────
    print(f"  Running scipy.linalg.eigh (driver=evd) ...", flush=True)
    t2 = time.time()
    ev_all, evec_all = scipy.linalg.eigh(H_dense, driver='evd')
    print(f"  Done in {time.time()-t2:.1f}s")

    del H_dense   # free memory: no longer needed

    # ── remove translational (near-zero) modes ───────────────────────────────
    phys      = ev_all > ZERO_EV
    ev        = ev_all[phys]
    evec      = evec_all[:, phys]
    n_phys    = ev.size
    n_removed = (~phys).sum()
    print(f"  Physical modes: {n_phys}  (removed {n_removed} near-zero)")

    # ── sanity: eigenvalue identity lon_s[k] = omega_k^2 ────────────────────
    # Check a random sample of 10 modes
    rng      = np.random.default_rng(42)
    idx_chk  = rng.integers(0, n_phys, 10)
    max_err  = 0.0
    for k in idx_chk:
        quad = evec[:, k].T @ H_sp @ evec[:, k]
        max_err = max(max_err, abs(quad - ev[k]) / ev[k])
    print(f"  Eigenvalue identity max_err (10 random modes): {max_err:.2e}"
          f"  {'PASSED' if max_err < 1e-3 else 'WARNING'}")

    return ev, evec, H_sp


# ─────────────────────────────────────────────────────────────────────────────
# Compute per-mode coupling terms (batched for memory)
# ─────────────────────────────────────────────────────────────────────────────
def compute_mode_couplings(ev, evec, ii, jj, r_mag, r_hat, vna, batch=50):
    """
    For each physical mode k compute:

      lon_s[k]   = sum_bonds (r_hat . Delta_e_k)^2  =  omega_k^2  (identity)

      trans_s[k] = sum_bonds |Delta_e_k^perp|^2
                 = sum_bonds |Delta_e_k|^2 - lon_s[k]
                 = u_eta_xx(Psi_k, Psi_k)
                 [TERM 1: geometric anharmonicity from r=|r_ij| expansion]

      tess3_s[k] = sum_bonds (1/r) *
                     [2*(r_hat . Delta_e_k) * (u_na_perp . Delta_e_k)
                      + (u_na . r_hat) * |Delta_e_k^perp|^2]
                 = tess3(u_na, Psi_k, Psi_k)
                 [TERM 2: u_na contracted with geometric cubic tensor]

    Processed in batches of `batch` modes to avoid allocating the full
    (N_modes x N_bonds x DIM) tensor at once.
    """
    n_phys  = len(ev)
    n_bonds = len(ii)

    lon_s   = np.zeros(n_phys)
    trans_s = np.zeros(n_phys)
    tess3_s = np.zeros(n_phys)

    # pre-flatten vna to per-bond form:  vna_ij = vna[j] - vna[i]
    vna_ij   = vna[jj] - vna[ii]                 # (M, 3)
    # longitudinal component of vna for each bond
    vna_lon  = np.einsum('md,md->m', vna_ij, r_hat)   # (M,)
    # transverse component of vna for each bond
    vna_lon3 = vna_lon[:, np.newaxis] * r_hat     # (M, 3) longitudinal part
    vna_perp = vna_ij - vna_lon3                  # (M, 3) transverse part
    # 1/r per bond
    inv_r    = 1.0 / r_mag                        # (M,)

    print(f"  Computing couplings for {n_phys} modes "
          f"(batch={batch}) ...", flush=True)
    t0 = time.time()

    for lo in range(0, n_phys, batch):
        hi  = min(lo + batch, n_phys)
        B   = hi - lo

        # eigenvector differences Delta_e_k for batch: shape (B, M, 3)
        Pm  = evec[:, lo:hi].T.reshape(B, N_ATOMS, DIM)
        De  = Pm[:, jj, :] - Pm[:, ii, :]        # (B, M, 3)

        # ── longitudinal projection ───────────────────────────────────────
        # lp[b,m] = r_hat[m] . De[b,m]
        lp   = np.einsum('bmd,md->bm', De, r_hat) # (B, M)
        lon_s[lo:hi] = (lp**2).sum(axis=1)        # = omega^2 (verified above)

        # ── transverse component of Delta_e_k ────────────────────────────
        De_lon  = lp[:, :, np.newaxis] * r_hat[np.newaxis, :, :]  # (B,M,3)
        De_perp = De - De_lon                                       # (B,M,3)
        trans_s[lo:hi] = np.sum(De_perp**2, axis=-1).sum(axis=1)  # sum over bonds

        # ── tess3 term ────────────────────────────────────────────────────
        # tess3[b] = sum_bonds (1/r) * [2*(r_hat.De)*(vna_perp.De)
        #                                + (vna.r_hat)*|De_perp|^2]
        #
        # part A: 2*(r_hat.De)*(vna_perp.De)   -- cross coupling
        vna_perp_dot = np.einsum('md,bmd->bm', vna_perp, De)       # (B,M)
        partA = 2.0 * lp * vna_perp_dot                             # (B,M)

        # part B: (vna.r_hat)*|De_perp|^2       -- u_na compression * transverse
        De_perp_sq = np.sum(De_perp**2, axis=-1)                   # (B,M)
        partB = vna_lon[np.newaxis, :] * De_perp_sq                 # (B,M)

        tess3_s[lo:hi] = ((partA + partB) * inv_r[np.newaxis, :]).sum(axis=1)

        if (lo // batch) % 20 == 0:
            pct = 100*(hi)/n_phys
            print(f"    {hi}/{n_phys} ({pct:.0f}%)  "
                  f"elapsed {time.time()-t0:.0f}s", flush=True)

    print(f"  Done in {time.time()-t0:.1f}s")
    return lon_s, trans_s, tess3_s


# ─────────────────────────────────────────────────────────────────────────────
# Assemble alpha components
# ─────────────────────────────────────────────────────────────────────────────
def compute_alpha(ev, lon_s, trans_s, tess3_s, K, V):
    """
    Assemble the three alpha contributions following the LJ glass convention:

      alpha_harm    = -(1/(2KVd)) * sum_k (1/omega_k^2) * 2*lon_s[k]
                    = -(1/(KVd)) * N_modes  =  -N_modes / (K*V*DIM)
                    [purely longitudinal, = harmonic Hessian term]

      alpha_anharm  = -(1/(2KVd)) * sum_k (1/omega_k^2) * [trans_s[k] + tess3_s[k]]
                    [geometric anharmonic: term1=transverse + term2=u_na coupling]

      alpha_total   = alpha_harm + alpha_anharm

    Factor -(1/(2*K*V*DIM)) matches eq. 25 of PRE 2010 and the convention in
    alpha_theory_FINAL.py.
    """
    n_phys = len(ev)
    prefac = -1.0 / (2.0 * K * V * DIM)

    w2 = 1.0 / ev     # 1/omega^2 weights, shape (n_phys,)

    # harmonic: 2 * lon_s[k] / omega_k^2 = 2 * omega_k^2 / omega_k^2 = 2
    # so sum = 2 * n_phys
    sum_harm   = 2.0 * n_phys          # exact from eigenvalue identity
    alpha_harm = prefac * sum_harm

    # anharmonic
    sum_anharm  = (w2 * (trans_s + tess3_s)).sum()
    alpha_anharm = prefac * sum_anharm

    alpha_total = alpha_harm + alpha_anharm

    return alpha_harm, alpha_anharm, alpha_total


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    s    = args.serial
    din  = args.data_dir
    dout = args.out_dir
    dout.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"alpha_theory_spring.py  —  serial {s}")
    print(f"{'='*60}")

    # ── geometry ─────────────────────────────────────────────────────────────
    L_BOX = 15.874010519682
    V_BOX = L_BOX**3
    print(f"\nN={N_ATOMS}  DIM={DIM}  N_modes={N_MODES}")
    print(f"V={V_BOX:.4f}")

    # ── Step 1: K ─────────────────────────────────────────────────────────────
    print("\n[Step 1] Loading K ...")
    K, P0 = load_K(din, s)

    # ── Step 2: bond geometry + vna ───────────────────────────────────────────
    print(f"\n[Step 2] Loading bonds (bonds_sn{s}.txt) ...")
    ii, jj, rij, r_mag, r_hat = load_bonds(din, s)

    print(f"\n[Step 2] Loading vna_eta (vna_eta_sn{s}.txt) ...")
    vna = load_vna(din, s)

    # ── Step 3: full Hessian diagonalisation ──────────────────────────────────
    print(f"\n[Step 3] Loading + diagonalising Hessian ...")
    ev, evec, H_sp = load_and_diagonalise(din, s)

    # ── Step 4: per-mode couplings ────────────────────────────────────────────
    print(f"\n[Step 4] Computing per-mode couplings ...")
    lon_s, trans_s, tess3_s = compute_mode_couplings(
        ev, evec, ii, jj, r_mag, r_hat, vna)

    # ── Step 5: eigenvalue identity check ─────────────────────────────────────
    print(f"\n[Step 5] Sanity checks ...")
    err_lon = np.abs(lon_s - ev) / (ev + 1e-30)
    print(f"  max|lon_s - omega^2|/omega^2 = {err_lon.max():.2e}  "
          f"{'PASSED' if err_lon.max() < 1e-3 else 'FAILED'}")
    print(f"  mean lon_s/omega^2 = {(lon_s/ev).mean():.6f}  (should be 1.0)")
    print(f"  mean trans_s/omega^2 = {(trans_s/ev).mean():.6f}")
    print(f"  mean tess3_s/omega^2 = {(tess3_s/ev).mean():.6f}")

    # ── Step 6: alpha ──────────────────────────────────────────────────────────
    print(f"\n[Step 6] Computing alpha ...")
    alpha_harm, alpha_anharm, alpha_total = compute_alpha(
        ev, lon_s, trans_s, tess3_s, K, V_BOX)

    print(f"\n  alpha_harmonic    = {alpha_harm:+.8e}")
    print(f"  alpha_anharmonic  = {alpha_anharm:+.8e}")
    print(f"  alpha_total       = {alpha_total:+.8e}")

    # ── Step 7: write output ───────────────────────────────────────────────────
    out_path = dout / "alpha_theory_spring_results.txt"
    header   = (f"# alpha_theory_spring.py output\n"
                f"# Network N={N_ATOMS} z=16.24 serial={s}\n"
                f"# K={K:.8f}  V={V_BOX:.6f}\n"
                f"# n_phys_modes={len(ev)}\n"
                f"# serial  alpha_harmonic  alpha_anharmonic  alpha_total\n")

    write_header = not out_path.exists()
    with open(out_path, 'a') as f:
        if write_header:
            f.write(header)
        f.write(f"glass_{s:04d}  "
                f"{alpha_harm:+.8e}  "
                f"{alpha_anharm:+.8e}  "
                f"{alpha_total:+.8e}\n")

    print(f"\n  Results appended to {out_path}")
    print(f"  Format matches alpha_theory_FINAL.py output")
    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
