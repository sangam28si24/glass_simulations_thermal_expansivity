# Alpha Theory Calculation 

## 1. The Formula 

$$\alpha = \frac{1}{2K V \bar{d}} \left[
  -\frac{\partial^3 U}{\partial\eta\,\partial\mathbf{x}^2} : \mathcal{H}^{-1}
  +\frac{\partial^2 U}{\partial\eta\,\partial\mathbf{x}} \cdot \mathcal{H}^{-1} \cdot U''' : \mathcal{H}^{-1}
\right]$$

where:
- $K$ = athermal bulk modulus
- $V$ = box volume
- $\bar{d} = 3$ (space dimension)
- $\eta$ = volumetric strain (isotropic deformation parameter, $\mathbf{H} = e^\eta \mathbf{I}$)
- $\mathcal{H}$ = Hessian (dynamical matrix), $\mathcal{H}_{k\ell} = \partial^2 U/\partial x_k\partial x_\ell$
- $U'''$ = third-derivative tensor, $U'''_{k\ell m} = \partial^3 U/\partial x_k\partial x_\ell\partial x_m$

### Term names used throughout

| Symbol | Name | Sign in α formula |
|--------|------|-------------------|
| $-\partial^3U/\partial\eta\partial\mathbf{x}^2:\mathcal{H}^{-1}$ | **Harmonic term** | negative |
| $+\partial^2U/\partial\eta\partial\mathbf{x}\cdot\mathcal{H}^{-1}\cdot U''':\mathcal{H}^{-1}$ | **Anharmonic term** | positive |

---

## 2. Normal-Mode Expansion

Substituting $\mathcal{H}^{-1} = \sum_\ell \frac{1}{\omega_\ell^2}\,\Psi_\ell \otimes \Psi_\ell$ (eigenvectors $\Psi_\ell$, eigenvalues $\omega_\ell^2$):

### Harmonic term expansion

$$
-\frac{\partial^3 U}{\partial \eta\,\partial \mathbf{x}^2} : \mathcal{H}^{-1}
= -\sum_\ell \frac{1}{\omega_\ell^2}\;
\underbrace{\sum_{k,m} \frac{\partial^3 U}{\partial \eta\,\partial x_k\,\partial x_m}
(\Psi_\ell)_k(\Psi_\ell)_m}_{\operatorname{uEtaXXonTwoVectors}(\Psi_\ell,\Psi_\ell)}
$$
### Anharmonic term expansion

Define the **compression mismatch force**:
$$\xi_\eta = \frac{\partial^2 U}{\partial\eta\,\partial\mathbf{x}} \quad \text{(3N-vector)}$$

The non-affine compression response is $\mathbf{u}_\text{na} = -\mathcal{H}^{-1}\xi_\eta$,
i.e. the CG solution of $\mathcal{H}\,\mathbf{u}_\text{na} = -\xi_\eta$.

Using $\mathcal{H}^{-1}\xi_\eta = -\mathbf{u}_\text{na}$ (note the sign):

$$+\xi_\eta\cdot\mathcal{H}^{-1}\cdot U''':\mathcal{H}^{-1}
= -\sum_\ell\frac{1}{\omega_\ell^2}\;
  \underbrace{\sum_{k,m,p}(\mathbf{u}_\text{na})_k\,U'''_{kmp}\,(\Psi_\ell)_m\,(\Psi_\ell)_p}_{
  \texttt{tessianOnThreeVectors}(\mathbf{u}_\text{na},\,\Psi_\ell,\,\Psi_\ell)}$$

### Combined formula used in code

$$
\boxed{
\alpha = -\frac{1}{2 K V \bar{d}} \sum_\ell \frac{1}{\omega_\ell^2}
\Big[
\operatorname{uEtaXX}(\Psi_\ell,\Psi_\ell)
+
\operatorname{tessThree}(\mathbf{u}_{\mathrm{na}},\Psi_\ell,\Psi_\ell)
\Big]
}
$$

Both terms enter with the **same sign** and the **same $-1/(2KV\bar{d})$ prefactor**.
This is the single expression implemented in `alpha_theory.py`.

---

## 3. Potential and Per-Bond Derivatives

The LJ glass potential (shifted and smoothed, cutoff $r/\sigma = 2.5$):

$$\varphi(r) = 4\!\left(\frac{\sigma^{12}}{r^{12}}-\frac{\sigma^6}{r^6}
  + C_6 r^6 + C_4 r^4 + C_2 r^2 + C_0\right)$$

Constants hard-coded in `3dlj_elasticity_analysis.c`:
```c
#define C_0   ( 0.075421526038318)
#define C_2   (-0.035416041040039)
#define C_4   ( 0.005708815153688)
#define C_6   (-0.000318029410118)
#define N_REP  12.0
#define N_ATT   6.0
#define CUTOFF  2.5
```

### Shorthand notation (PDF page 2)

$$[1\text{st}] = \frac{\varphi'}{r}, \quad
[2\text{nd}] = \frac{\varphi''}{r^2}-\frac{\varphi'}{r^3}, \quad
[3\text{rd}] = \frac{\varphi'''}{r^3}-\frac{3\varphi''}{r^4}+\frac{3\varphi'}{r^5}$$

### How C stores them — inside `calculateEverything()`

Per contact bond $\ell$ between atoms $i, j$:

```c
sigma2OverR2 = 1.0 / r2OverSigma2;
sigma6OverR6 = sigma2OverR2 * sigma2OverR2 * sigma2OverR2;
sigma8OverR8 = sigma6OverR6 * sigma2OverR2;

// goo = φ'(r)/r in actual (not scaled) distance  →  this equals −[1st]·σ²
goo = 4.0*(sigma8OverR8*(N_REP*sigma6OverR6 - N_ATT)
           - 2.0*C_2 - r2OverSigma2*(4.0*C_4 + 6.0*r2OverSigma2*C_6)) * invSigmaSqrd;

third[l]  = 4.0*(-2688.0*sigma8OverR8*sigma8OverR8*sigma2OverR2
                 + 480.0*sigma6OverR6*sigma6OverR6
                 + 48.0*C_6) * (invSigmaSqrd*invSigmaSqrd*invSigmaSqrd);  // = [3rd]

second[l] = 4.0*(168.0*sigma8OverR8*sigma8OverR8
                - 48.0*sigma8OverR8*sigma2OverR2
                + 24.0*C_6*r2OverSigma2
                + 8.0*C_4) * (invSigmaSqrd*invSigmaSqrd);                  // = [2nd]

first[l]  = -goo;                                                            // = [1st] = +φ'(r)/r

// actual displacement vector x_j - x_i (real units, not fractional)
rij[DIM*l + X_COMP] = L*dx;
rij[DIM*l + Y_COMP] = L*dy;
rij[DIM*l + Z_COMP] = L*dz;
```

**Key sign:** `first[l] = -goo = +φ'(r)/r` = $+[1\text{st}]$.  
`second[l]` = $[2\text{nd}]$, `third[l]` = $[3\text{rd}]$.

---

## 4. C Functions Used by `alpha_theory.py`

### 4.1 `compression_nonaffine_velocities()` → writes `vna_eta_g{g}.txt`

```c
int compression_nonaffine_velocities(double *res){
    double *xi = malloc(DIM*N * sizeof(double));

    // Build ξ_η = ∂²U/∂η∂x  (PDF eq. 14)
    // For each contact l:  (second[l]*r² + first[l]) acts as the prefactor
    for (l=0; l<contacts; l++){
        r2 = sum over d of rij[DIM*l+d]^2;
        for (d=0; d<DIM; d++){
            xi[DIM*j + d] += (second[l]*r2 + first[l]) * rij[DIM*l + d];
            xi[DIM*i + d] -= (second[l]*r2 + first[l]) * rij[DIM*l + d];
        }
    }

    cgSolver(res, xi);       // solves H·res = ξ_η

    // Negate to get u_na = −H⁻¹ξ_η
    for (a=0; a<DIM*N; a++)
        res[a] = -res[a];
}
```

**Output `vna_eta_g{g}.txt`:** N rows × 3 columns. Each row $i$ is
$\bigl(u_\text{na}\bigr)_{i,x},\; \bigl(u_\text{na}\bigr)_{i,y},\; \bigl(u_\text{na}\bigr)_{i,z}$.

> **Sign to verify:** The last line negates the CG solution.
> The written file therefore contains $\mathbf{u}_\text{na} = -\mathcal{H}^{-1}\xi_\eta$.
> `alpha_theory.py` uses this directly — it does NOT negate again.

> **Sanity check:** Print `typicalGrad / typicalInteractionStrength` (C stdout).
> Should be $\lesssim 10^{-8}$. If larger, the snapshot is not well-minimised
> and all athermal expressions are unreliable.

### 4.2 `u_eta_x_x_on_two_vectors(vec1, vec2)` — harmonic integrand

Implements $\sum_{k,m}(\partial^3U/\partial\eta\partial x_k\partial x_m)\,v^{(1)}_k\,v^{(2)}_m$ (PDF eq. 15):

```c
double u_eta_x_x_on_two_vectors(double *vec1, double *vec2){
    double res = 0.0;
    for (l=0; l<contacts; l++){
        r2        = sum_d rij[d]^2
        v1_dot_r  = sum_d rij[d] * (vec1[j,d] - vec1[i,d])
        v2_dot_r  = sum_d rij[d] * (vec2[j,d] - vec2[i,d])
        v1_dot_v2 = sum_d (vec1[j,d] - vec1[i,d]) * (vec2[j,d] - vec2[i,d])

        res += third[l]*r2*v1_dot_r*v2_dot_r
             + second[l]*(2.0*v1_dot_r*v2_dot_r + r2*v1_dot_v2);
    }
    return res;
}
```

Called with `vec1 = vec2 = Ψ_ℓ`.

### 4.3 `tessianOnThreeVectors(v, u, w)` — anharmonic integrand

Implements $\sum_{k,m,p} U'''_{kmp}\,v_k\,u_m\,w_p$:

```c
double tessianOnThreeVectors(double *v, double *u, double *w){
    double res = 0.0;
    for (l=0; l<contacts; l++){
        u_dot_r = sum_d rij[d]*(u[j,d]-u[i,d])
        v_dot_r = sum_d rij[d]*(v[j,d]-v[i,d])
        w_dot_r = sum_d rij[d]*(w[j,d]-w[i,d])
        v_dot_u = sum_d (v[j,d]-v[i,d])*(u[j,d]-u[i,d])
        u_dot_w = sum_d (u[j,d]-u[i,d])*(w[j,d]-w[i,d])
        w_dot_v = sum_d (w[j,d]-w[i,d])*(v[j,d]-v[i,d])

        res += third[l]*u_dot_r*v_dot_r*w_dot_r
             + second[l]*(v_dot_u*w_dot_r + u_dot_w*v_dot_r + w_dot_v*u_dot_r);
    }
    return res;
}
```

Called as `tessianOnThreeVectors(u_na, Ψ_ℓ, Ψ_ℓ)`, i.e. `v = u_na`, `u = w = Ψ_ℓ`.

### 4.4 `get_bulk_modulus()` - K used in α prefactor

```c
int get_bulk_modulus(double *k_born, double *k_na){
    *k_born = 0.0;
    for (l=0; l<contacts; l++){
        r2 = sum_q rij[q]^2
        *k_born += (second[l]*r2 + first[l])*r2 + first[l]*r2;
        // (also builds ξ_η in xi[] simultaneously)
    }
    *k_born = *k_born / V;

    cgSolver(nav, xi);            // nav = H⁻¹ξ_η  (no negation here)
    *k_na = cdot(nav, xi) / V;
}
```

Then in `main()`:
```c
bulk_modulus = (k_born - k_na) / (double)(DIM*DIM) + pressure;
```

This implements $K = (K_\text{Born} - K_\text{NA})/d^2 + p$ (PDF eq. 6).  
Column 4 of `3dlj_elasticity_N4000_s000.dat` contains this value.

---

## 5. Data Flow: C to Python

### Files written by C, per glass `g`

| File | Content | Column/shape |
|------|---------|--------------|
| `3dlj_elasticity_N4000_s000.dat` | Elasticity results | see §5.1 |
| `bonds_g{g}.txt` | Per-contact geometry + derivatives | 8 columns |
| `meta_g{g}.txt` | Scalar globals | key-value pairs |
| `hessian_g{g}.txt` | Sparse COO Hessian (text) | `row  col  val` |
| `vna_eta_g{g}.txt` | Non-affine compression velocities | N rows × 3 cols |

### 5.1 Column mapping for `3dlj_elasticity_N4000_s000.dat`

```
col 0 : glass index (q)
col 1 : u/N   (potential energy per atom)
col 2 : P     (virial pressure at T=0)
col 3 : G     (shear modulus = g_born - g_na)
col 4 : K     (bulk modulus  = (k_born - k_na)/d² + P)
col 5 : first_order_dilatancy
col 6 : second_order_dilatancy
```
---

## 6. Python Implementation (`alpha_theory.py`)

### 6.1 Load bond data

```python
data     = np.loadtxt(f"bonds_g{g}.txt", comments="#")
ii       = data[:, 0].astype(np.int32)   # atom i per bond
jj       = data[:, 1].astype(np.int32)   # atom j per bond  (always j > i, as in C)
rij_arr  = data[:, 2:5]                  # x_j - x_i  in real units (L*dx, L*dy, L*dz)
first_b  = data[:, 5]                    # [1st] = φ'(r)/r
second_b = data[:, 6]                    # [2nd]
third_b  = data[:, 7]                    # [3rd]
```

### 6.2 Load non-affine compression velocities

```python
vna_eta = np.loadtxt(f"vna_eta_g{g}.txt").reshape(N_ATOMS, DIM)
# shape: (4000, 3)
# vna_eta[i, d]  =  (u_na)_{i,d}  =  −(H⁻¹ξ_η)_{i,d}
```

### 6.3 Diagonalise Hessian

```python
H_dense = H.toarray()                            # (12000, 12000) dense symmetric
evals, evecs = scipy.linalg.eigh(H_dense)        # LAPACK DSYEVD

# Remove translation zero modes
valid        = evals > ZERO_MODE_CUTOFF           # threshold = 1e-4
valid_evals  = evals[valid]                       # ω² values, shape (n_valid,)
valid_evecs  = evecs[:, valid]                    # shape (12000, n_valid)
```

### 6.4 Bond-projected eigenvectors (batched over modes)

```python
# bsz = batch size (up to BATCH_SIZE=64 modes at once)
Psi_mat  = valid_evecs[:, lo:hi].T.reshape(bsz, N_ATOMS, DIM)
Psi_diff = Psi_mat[:, jj, :] - Psi_mat[:, ii, :]
# shape: (bsz, contacts, DIM)
# mirrors C:  uij[q] = vec[DIM*j+q] - vec[DIM*i+q]
```

### 6.5 Harmonic term: `u_eta_x_x_batch`

Mirrors `u_eta_x_x_on_two_vectors(Ψ, Ψ)`:

```python
psi_dot_r = np.einsum("bcd,cd->bc", Psi_diff, rij_arr)    # (bsz, contacts)  = v_dot_r
psi_sq    = np.einsum("bcd,bcd->bc", Psi_diff, Psi_diff)   # (bsz, contacts)  = v1_dot_v2

contrib = (third_b * r2 * psi_dot_r**2
           + second_b * (2.0 * psi_dot_r**2 + r2 * psi_sq))

term1_per_mode[lo:hi] = contrib.sum(axis=1)    # sum over contacts → (bsz,)
```

**Variable-by-variable match with C:**

| C variable | Python variable | Formula |
|-----------|----------------|---------|
| `v1_dot_r` | `psi_dot_r` | $\Psi_{ij}\cdot \mathbf{r}_{ij}$ |
| `v2_dot_r` | `psi_dot_r` | same (vec1=vec2=Ψ) |
| `v1_dot_v2` | `psi_sq` | $\|\Psi_{ij}\|^2$ |
| `third[l]*r2*v1_dot_r*v2_dot_r` | `third_b*r2*psi_dot_r**2` | identical |
| `second[l]*(2*v1_dot_r*v2_dot_r + r2*v1_dot_v2)` | `second_b*(2*psi_dot_r**2 + r2*psi_sq)` | identical |

### 6.6 Anharmonic term: `tess3_batch`

Mirrors `tessianOnThreeVectors(u_na, Ψ, Ψ)` with `v = u_na`, `u = w = Ψ`:

```python
# Precomputed outside the mode loop (constant for all modes):
vna_diff  = vna_eta[jj] - vna_eta[ii]                          # (contacts, DIM)
vna_dot_r = np.einsum("cd,cd->c", vna_diff, rij_arr)            # (contacts,)

# Inside mode loop:
psi_dot_r   = np.einsum("bcd,cd->bc", Psi_diff, rij_arr)        # (bsz, contacts)
psi_sq      = np.einsum("bcd,bcd->bc", Psi_diff, Psi_diff)      # (bsz, contacts)
vna_dot_psi = np.einsum("cd,bcd->bc", vna_diff, Psi_diff)       # (bsz, contacts)

contrib = (third_b * psi_dot_r**2 * vna_dot_r
           + second_b * (2.0 * vna_dot_psi * psi_dot_r
                          + psi_sq * vna_dot_r))

term2_per_mode[lo:hi] = contrib.sum(axis=1)
```


| C variable | Python variable | Note |
|-----------|----------------|------|
| `u_dot_r` | `psi_dot_r` | Ψ·r |
| `w_dot_r` | `psi_dot_r` | same, u=w=Ψ |
| `v_dot_r` | `vna_dot_r` | u_na·r |
| `v_dot_u` | `vna_dot_psi` | u_na·Ψ |
| `u_dot_w` | `psi_sq` | Ψ·Ψ = |Ψ|² |
| `w_dot_v` | `vna_dot_psi` | same as v_dot_u |
| C formula: `third*u_dot_r*v_dot_r*w_dot_r` | `third_b*psi_dot_r**2*vna_dot_r` | ✓ |
| C formula: `second*(v_dot_u*w_dot_r + u_dot_w*v_dot_r + w_dot_v*u_dot_r)` | `second_b*(2*vna_dot_psi*psi_dot_r + psi_sq*vna_dot_r)` | ✓ |

### 6.7 Final assembly: `alpha_from_mode_sums`

```python
w   = 1.0 / valid_evals           # 1/ω²  per mode,  shape (n_valid,)
pre = -1.0 / (2.0 * K * V * DIM)  # = −1/(2KVd̄)

alpha_harmonic   = pre * np.dot(w, term1_per_mode)
alpha_anharmonic = pre * np.dot(w, term2_per_mode)
alpha_total      = alpha_harmonic + alpha_anharmonic
```

This directly implements:
$$\alpha = -\frac{1}{2KV\bar{d}}\sum_\ell\frac{1}{\omega_\ell^2}
\bigl[\text{term1}_\ell + \text{term2}_\ell\bigr]$$

---

## 7. Sign Convention — Complete Chain

This is the most likely source of errors. Full derivation:

```
PDF eq. 24 (anharmonic piece):  + ξ_η · H⁻¹ · U''' : H⁻¹

Expand H⁻¹ = Σ_ℓ (1/ω²) Ψ_ℓ⊗Ψ_ℓ :
           = + Σ_ℓ (1/ω²) (H⁻¹ξ_η) · U'''(·, Ψ_ℓ, Ψ_ℓ)

C code:  u_na = −H⁻¹ξ_η   (negation at end of compression_nonaffine_velocities)
⟹  H⁻¹ξ_η = −u_na

           = − Σ_ℓ (1/ω²) u_na · U'''(·, Ψ_ℓ, Ψ_ℓ)
           = − Σ_ℓ (1/ω²) tessianOnThreeVectors(u_na, Ψ_ℓ, Ψ_ℓ)
           = − Σ_ℓ (1/ω²) term2_ℓ

So anharmonic contribution to α:
  +1/(2KVd) × (anharmonic piece) = −1/(2KVd) × Σ_ℓ (1/ω²) term2_ℓ

Harmonic contribution:
  −1/(2KVd) × Σ_ℓ (1/ω²) term1_ℓ

Combined:
  pre = −1/(2KVd)
  α = pre × Σ_ℓ (1/ω²) [term1_ℓ + term2_ℓ]       exact formula in code
```

**Conclusion:** Both terms enter with the same negative prefactor `pre`.
The code is correct **only if** the file `vna_eta_g{g}.txt` carries
the sign $\mathbf{u}_\text{na} = -\mathcal{H}^{-1}\xi_\eta$.

**Confirm in C:**
```c
// End of compression_nonaffine_velocities():
if ( cgSolver(res, xi) ){
    for (a=0; a<DIM*N; a++)
        res[a] = -res[a];     // ← this negation produces u_na = −H⁻¹ξ_η
    free(xi);
    return 1;
}
```
And `write_vna_eta_text()` writes this negated result as-is.

---

## 8. Complete Call Graph

```
main()  [C]
 ├── readSnapShot()
 │    └── initializeSystem()
 │         ├── updateNebzLists()
 │         └── calculateEverything()         ← populates first[], second[], third[], rij[]
 │                                               contacts, lookupBond, pressure, V, L
 ├── get_shear_modulus()  → G  (col 3 of .dat)
 ├── get_bulk_modulus()   → K  (col 4 of .dat, includes pressure correction)
 ├── get_first_order_dilatancy()   (col 5)
 ├── get_second_order_dilatancy()  (col 6)
 ├── fprintf(outFile, ...)         → 3dlj_elasticity_N4000_s000.dat
 │
 ├── write_bonds_text(q)           → bonds_g{q}.txt
 ├── write_meta_text(q)            → meta_g{q}.txt
 ├── write_hessian_text(q)         → hessian_g{q}.txt
 └── write_vna_eta_text(q)
      └── compression_nonaffine_velocities()
           ├── builds ξ_η from second[], first[], rij[]
           ├── cgSolver(res, xi)   → res = H⁻¹ξ_η
           └── negates res         → u_na = −H⁻¹ξ_η  → vna_eta_g{q}.txt

alpha_theory.py
 ├── load_elasticity_results()    → K per glass (col 4)
 ├── load_bonds()                 → ii, jj, rij_arr, first_b, second_b, third_b
 ├── load_meta()                  → V
 ├── load_hessian()               → sparse H (12000×12000)
 ├── load_vna_eta()               → u_na, shape (N, 3)
 ├── scipy.linalg.eigh(H_dense)   → ω², Ψ  (full diagonalisation)
 │
 ├── [mode loop, batched over BATCH_SIZE=64 modes]
 │    ├── Psi_diff = Psi[jj] - Psi[ii]
 │    ├── u_eta_x_x_batch(Psi_diff, rij_arr, r2, second_b, third_b)
 │    │    → term1_per_mode[ℓ]  =  [3rd]·r²·(Ψ·r)² + [2nd]·(2(Ψ·r)² + r²|Ψ|²)
 │    └── tess3_batch(Psi_diff, rij_arr, vna_diff, vna_dot_r, r2, second_b, third_b)
 │         → term2_per_mode[ℓ]  =  [3rd]·(Ψ·r)²·(u_na·r)
 │                                + [2nd]·(2·(u_na·Ψ)·(Ψ·r) + |Ψ|²·(u_na·r))
 │
 └── alpha_from_mode_sums(ω², term1, term2, K, V)
      ├── pre = −1/(2KVd)
      ├── α_harm    = pre · Σ_ℓ (1/ω_ℓ²) · term1_ℓ
      ├── α_anharm  = pre · Σ_ℓ (1/ω_ℓ²) · term2_ℓ
      └── α = α_harm + α_anharm
```
