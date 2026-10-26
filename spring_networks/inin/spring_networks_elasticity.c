/* =============================================================================
 * spring_networks_elasticity.c
 *
 * PURPOSE
 * -------
 * Compute the zero-temperature elastic moduli and auxiliary data files for a
 * 3D spring network.  Outputs are read by the Python notebook to compute the
 * thermal expansion coefficient via  alpha = (1/K) * dP_vir/dT.
 *
 * OUTPUTS
 *   snp_elasticity_N4000_z16.24_s<s>.dat  —  scalar results (one row per serial)
 *     columns: serial  P0  K_born  K_na  K  mu_born  mu_na  mu
 *
 *   bonds_sn<s>.txt      —  per-bond data for Python mode-sum analysis
 *   hessian_sn<s>.txt    —  sparse COO Hessian (row col val) for eigsh()
 *   vna_eta_sn<s>.txt    —  compression non-affine velocity (N×3)
 *
 * USAGE:   ./snp_elast <serial>
 * EXAMPLE: ./snp_elast 0
 *
 * COMPILE: gcc -O2 -o snp_elast spring_networks_elasticity.c -lm
 *
 * =============================================================================
 * SOURCE ATTRIBUTION
 *
 *  readNetwork()          — 3D_relaxed_spring_network_analysis.c
 *  cdot()                 — 3D_relaxed_spring_network_analysis.c
 *  hessianOnVector()      — 3D_relaxed_spring_network_analysis.c
 *  cgSolver()             — 3D_relaxed_spring_network_analysis.c
 *                           (uses maxNorm stopping criterion as in that file)
 *  u_epsilon()            — 3D_relaxed_spring_network_analysis.c
 *  u_epsilon_epsilon()    — 3D_relaxed_spring_network_analysis.c
 *  u_epsilon_x()          — 3D_relaxed_spring_network_analysis.c
 *  calculateBulkModulus() — 3D_relaxed_spring_network_analysis.c
 *  calculateShearModulus()— 3D_relaxed_spring_network_analysis.c
 *  test_network()         — 3D_relaxed_spring_network_analysis.c
 *  writeHessian()         — 3D_relaxed_spring_network_analysis.c (adapted:
 *                           output file name now parameterised by serial)
 *  write_bonds_text()     — 3dlj_elasticity_analysis.c (adapted for springs:
 *                           third=0 since d^3U/dr^3=0 for harmonic springs)
 *  write_hessian_text()   — 3dlj_elasticity_analysis.c (adapted, plain-text
 *                           sparse COO; replaces binary writeHessian)
 *  write_vna_eta_text()   — 3dlj_elasticity_analysis.c (adapted)
 *  compression_nonaffine_velocities() — 3dlj_elasticity_analysis.c (adapted)
 *  free_everything()      — 3D_relaxed_spring_network_analysis.c
 *  main()                 — new
 * =============================================================================
 */

#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* ── system size ─────────────────────────────────────────────────────────── */
#define  N              4000
#define  DOUBLE_N       ((double)N)
#define  DIM            3
#define  CGSOLVER_TOL   1e-10

/* ── index aliases ───────────────────────────────────────────────────────── */
#define  DELTA(a,b)  ((a)==(b) ? 1.0 : 0.0)
#define  X_COMP  0
#define  Y_COMP  1
#define  Z_COMP  2
#define  I_INDEX 0
#define  J_INDEX 1

/* ── PBC geometry ────────────────────────────────────────────────────────── */
#define  LENGTH         1.0
#define  HALF_LENGTH    0.5

/* ── misc ────────────────────────────────────────────────────────────────── */
#define  HUGE_INT  10000000

/* function pointer for CG solver */
void (*linearOperator)(double *actOn, double *res);

/* ── globals ─────────────────────────────────────────────────────────────── */
double L, invL, V;
double pressure = 0.0;   /* T=0 virial pressure (0 for stress-free network)  */
int    contacts;

double *rx, *ry, *rz;

/* ── bond network ────────────────────────────────────────────────────────── */
int    *lookupBond;
double *second;
double *first;
double *rij;
double *inv_r;

/* ─────────────────────────────────────────────────────────────────────────────
 * readNetwork()  —  SOURCE: 3D_relaxed_spring_network_analysis.c
 *
 * Reads the spring network file.  For a stress-free network the springs are
 * at their natural length, so:
 *   first[l]  = 0        (no pre-stress force)
 *   second[l] = 1/r²     (spring stiffness = 1, normalised by r²)
 * These values feed into hessianOnVector() and the modulus functions.
 * ─────────────────────────────────────────────────────────────────────────── */
int readNetwork(char *inFileName) {
    int i, j, l, foo;
    double dummy, dx, dy, dz, r2;
    FILE *inFile;

    rx = (double *)malloc(sizeof(double) * N);
    ry = (double *)malloc(sizeof(double) * N);
    rz = (double *)malloc(sizeof(double) * N);

    inFile = fopen(inFileName, "rb");
    if (!inFile) { printf("readNetwork: cannot open %s\n", inFileName); return 0; }

    foo = fscanf(inFile, "%lf %lf %lf", &L, &dummy, &dummy);
    invL = 1.0 / L;
    V    = L * L * L;

    for (i = 0; i < N; i++)
        foo = fscanf(inFile, "%lf %lf %lf", &rx[i], &ry[i], &rz[i]);

    foo = fscanf(inFile, "%d %lf %lf", &contacts, &dummy, &dummy);

    lookupBond = (int    *)malloc(sizeof(int)    * 2 * contacts);
    rij        = (double *)malloc(sizeof(double) * DIM * contacts);
    inv_r      = (double *)malloc(sizeof(double) *       contacts);
    first      = (double *)malloc(sizeof(double) *       contacts);
    second     = (double *)malloc(sizeof(double) *       contacts);

    for (l = 0; l < contacts; l++) {
        foo = fscanf(inFile, "%d %d %lf",
                     &lookupBond[2*l + I_INDEX],
                     &lookupBond[2*l + J_INDEX], &dummy);
        i  = lookupBond[2*l + I_INDEX];
        j  = lookupBond[2*l + J_INDEX];
        dx = rx[j] - rx[i];
        dy = ry[j] - ry[i];
        dz = rz[j] - rz[i];
        if      (dz >=  HALF_LENGTH) dz -= LENGTH;
        else if (dz < -HALF_LENGTH)  dz += LENGTH;
        if      (dy >=  HALF_LENGTH) dy -= LENGTH;
        else if (dy < -HALF_LENGTH)  dy += LENGTH;
        if      (dx >=  HALF_LENGTH) dx -= LENGTH;
        else if (dx < -HALF_LENGTH)  dx += LENGTH;

        r2              = L*L * (dx*dx + dy*dy + dz*dz);
        rij[DIM*l+X_COMP] = L * dx;
        rij[DIM*l+Y_COMP] = L * dy;
        rij[DIM*l+Z_COMP] = L * dz;
        inv_r[l]  = 1.0 / sqrt(r2);
        first[l]  = 0.0;       /* stress-free: no pre-stress                  */
        second[l] = 1.0 / r2; /* Hessian coefficient = k/r² with k=1         */
    }
    fclose(inFile);
    return 1;
}

/* ─────────────────────────────────────────────────────────────────────────────
 * cdot()  —  SOURCE: 3D_relaxed_spring_network_analysis.c
 * ─────────────────────────────────────────────────────────────────────────── */
double cdot(double *vec1, double *vec2) {
    double res = 0.0;
    int i;
    for (i = 0; i < DIM*N; i++) res += vec1[i] * vec2[i];
    return res;
}

/* ─────────────────────────────────────────────────────────────────────────────
 * hessianOnVector()  —  SOURCE: 3D_relaxed_spring_network_analysis.c
 *
 * Applies the dynamical matrix H to a trial vector vec and stores the result
 * in res.  For spring networks (first[l]=0):
 *   H_{ij}^{αβ} = second[l] * rij_α * rij_β = (1/r²) * rij_α * rij_β
 * This is the standard harmonic-spring Hessian.
 * ─────────────────────────────────────────────────────────────────────────── */
void hessianOnVector(double *vec, double *res) {
    int q, i, j, l;
    double uij[DIM], u_dot_r;

    for (q = 0; q < DIM; q++)
        for (i = 0; i < N; i++)
            res[DIM*i + q] = 0.0;

    for (l = 0; l < contacts; l++) {
        i = lookupBond[2*l + I_INDEX];
        j = lookupBond[2*l + J_INDEX];
        for (q = 0; q < DIM; q++)
            uij[q] = vec[DIM*j + q] - vec[DIM*i + q];

        for (u_dot_r = 0.0, q = 0; q < DIM; q++)
            u_dot_r += rij[DIM*l + q] * uij[q];

        for (q = 0; q < DIM; q++) {
            res[DIM*i + q] -= second[l]*u_dot_r*rij[DIM*l + q] + first[l]*uij[q];
            res[DIM*j + q] += second[l]*u_dot_r*rij[DIM*l + q] + first[l]*uij[q];
        }
    }
}

/* ─────────────────────────────────────────────────────────────────────────────
 * cgSolver()  —  SOURCE: 3D_relaxed_spring_network_analysis.c
 *
 * Conjugate-gradient solver for H x = b (H = linearOperator).
 * Solves for the non-affine response vector x.
 * Returns 1 on convergence, 0 on failure.
 * ─────────────────────────────────────────────────────────────────────────── */
int cgSolver(double *x, double *b) {
    double *r, *w, *z;
    double a, B, sum, scale, maxNorm, typicalNorm;
    int i, k, q;

    scale = sqrt(cdot(b, b) / ((double)DIM * DOUBLE_N));

    r = (double *)malloc(sizeof(double) * DIM * N);
    w = (double *)malloc(sizeof(double) * DIM * N);
    z = (double *)malloc(sizeof(double) * DIM * N);

    linearOperator(x, r);
    for (i = 0; i < DIM*N; i++) { r[i] = b[i] - r[i]; w[i] = -r[i]; }
    linearOperator(w, z);
    a = cdot(r, w) / cdot(w, z);
    for (i = 0; i < DIM*N; i++) x[i] += a * w[i];
    B = 0.0;

    for (k = 0; k < HUGE_INT; k++) {
        maxNorm = 0.0; typicalNorm = 0.0;
        for (i = 0; i < DIM*N; i++) {
            r[i] -= a * z[i];
            if (r[i]*r[i] > maxNorm) maxNorm = r[i]*r[i];
            typicalNorm += r[i]*r[i];
        }
        maxNorm = sqrt(maxNorm);
        typicalNorm = sqrt(typicalNorm / ((double)DIM * DOUBLE_N));
        if (maxNorm < CGSOLVER_TOL) break;
        if (k > 0 && !(k % 16384))
            printf("  CG it%d  typNorm=%g  maxNorm=%g\n", k, typicalNorm, maxNorm);
        B = cdot(r, z) / cdot(w, z);
        for (i = 0; i < DIM*N; i++) w[i] = -r[i] + B * w[i];
        linearOperator(w, z);
        a = cdot(r, w) / cdot(w, z);
        for (i = 0; i < DIM*N; i++) x[i] += a * w[i];
    }
    /* remove residual drift from solution */
    for (q = 0; q < DIM; q++) {
        sum = 0.0;
        for (i = 0; i < N; i++) sum += x[DIM*i + q];
        sum /= DOUBLE_N;
        for (i = 0; i < N; i++) x[DIM*i + q] -= sum;
    }
    free(r); free(w); free(z);
    return (k < HUGE_INT) ? 1 : 0;
}

/* ─────────────────────────────────────────────────────────────────────────────
 * u_epsilon()  —  SOURCE: 3D_relaxed_spring_network_analysis.c
 * First strain derivative of the potential: ∂U/∂ε_{αβ}
 * = Σ_l first[l] * rij_α * rij_β
 * For stress-free network (first[l]=0): always 0.
 * Kept for generality and to document the formula.
 * ─────────────────────────────────────────────────────────────────────────── */
double u_epsilon(int alpha, int beta) {
    int l;
    double res = 0.0;
    for (l = 0; l < contacts; l++)
        res += first[l] * rij[DIM*l + alpha] * rij[DIM*l + beta];
    return res;
}

/* ─────────────────────────────────────────────────────────────────────────────
 * u_epsilon_epsilon()  —  SOURCE: 3D_relaxed_spring_network_analysis.c
 * Second strain derivative: ∂²U/∂ε_{αβ}∂ε_{κχ}
 * = Σ_l second[l] * rij_α * rij_β * rij_κ * rij_χ
 * This is the Born (affine) contribution to the elastic tensor.
 * For spring networks: second[l]=1/r², so this picks out the 4th moment of
 * bond-vector directions weighted by 1/r².
 * ─────────────────────────────────────────────────────────────────────────── */
double u_epsilon_epsilon(int alpha, int beta, int kappa, int chi) {
    int l;
    double res = 0.0;
    for (l = 0; l < contacts; l++)
        res += second[l]
             * rij[DIM*l + alpha] * rij[DIM*l + beta]
             * rij[DIM*l + kappa] * rij[DIM*l + chi];
    return res;
}

/* ─────────────────────────────────────────────────────────────────────────────
 * u_epsilon_x()  —  SOURCE: 3D_relaxed_spring_network_analysis.c
 * Mixed strain–displacement derivative: ∂²U/∂ε_{αβ}∂u_i^q
 * Stored as a 3N-vector in res[].
 * Used to build the right-hand side for the CG solver (non-affine driving
 * force) in both calculateBulkModulus and calculateShearModulus.
 * ─────────────────────────────────────────────────────────────────────────── */
void u_epsilon_x(int alpha, int beta, double *res) {
    int q, i, j, l;

    for (q = 0; q < DIM; q++)
        for (i = 0; i < N; i++)
            res[DIM*i + q] = 0.0;

    for (l = 0; l < contacts; l++) {
        i = lookupBond[2*l + I_INDEX];
        j = lookupBond[2*l + J_INDEX];
        for (q = 0; q < DIM; q++) {
            double contrib = second[l] * rij[DIM*l+alpha] * rij[DIM*l+beta] * rij[DIM*l+q]
                           + first[l]  * (DELTA(alpha,q)*rij[DIM*l+beta]
                                        + DELTA(beta, q)*rij[DIM*l+alpha]);
            res[DIM*i + q] -= contrib;
            res[DIM*j + q] += contrib;
        }
    }
}

/* ─────────────────────────────────────────────────────────────────────────────
 * calculateBulkModulus()  —  SOURCE: 3D_relaxed_spring_network_analysis.c
 *
 * Computes the isothermal bulk modulus K = K_born − K_na via:
 *
 *   K_born = (1/V) Σ_l [(second[l]*r² + first[l])*r² + first[l]*r²]
 *          = (1/V) Σ_l r²      (for first=0, second=1/r²)
 *          = contacts / V      (= z*ρ/2 in terms of coordination and density)
 *
 *   K_na   = (1/V) ξ · H⁻¹ · ξ    where ξ[j] = Σ_{l∋j}(second[l]*r²)*rij_l
 *
 *   K      = (K_born − K_na) / (DIM*DIM)  +  P_0
 *           = (K_born − K_na) / 9          (P_0 = 0 for stress-free net)
 *
 * The factor of DIM*DIM comes from the isotropic average over all 9 strain
 * components in 3D (three diagonal components each counted once).
 * ─────────────────────────────────────────────────────────────────────────── */
int calculateBulkModulus(double *k_born, double *k_na) {
    int a, q, i, j, l;
    double r2;
    double *nav, *xi;

    nav = (double *)malloc(sizeof(double) * N * DIM);
    xi  = (double *)malloc(sizeof(double) * N * DIM);
    for (a = 0; a < DIM*N; a++) { nav[a] = 0.0; xi[a] = 0.0; }
    *k_born = 0.0;

    for (l = 0; l < contacts; l++) {
        i  = lookupBond[2*l + I_INDEX];
        j  = lookupBond[2*l + J_INDEX];
        r2 = 0.0;
        for (q = 0; q < DIM; q++) r2 += rij[DIM*l+q] * rij[DIM*l+q];

        /* Born term: (second*r² + first)*r² + first*r² */
        *k_born += (second[l]*r2 + first[l])*r2 + first[l]*r2;

        /* driving force for non-affine displacement */
        for (q = 0; q < DIM; q++) {
            xi[DIM*j + q] += (second[l]*r2 + first[l]) * rij[DIM*l+q];
            xi[DIM*i + q] -= (second[l]*r2 + first[l]) * rij[DIM*l+q];
        }
    }
    *k_born /= V;

    if (cgSolver(nav, xi)) {
        *k_na = cdot(nav, xi) / V;
        free(nav); free(xi);
        return 1;
    }
    free(nav); free(xi);
    return 0;
}

/* ─────────────────────────────────────────────────────────────────────────────
 * calculateShearModulus()  —  SOURCE: 3D_relaxed_spring_network_analysis.c
 *
 * Computes the shear modulus μ for the xy shear strain ε_xy:
 *
 *   μ_born = [u_εε(x,y,x,y) + u_ε(y,y)] / V
 *          = Σ_l second[l]*rij_x²*rij_y² / V    (first=0 term vanishes)
 *
 *   μ_na   = (1/V) ξ · H⁻¹ · ξ
 *            where ξ = u_ε_x(x, y, ·) (shear driving force on each node)
 *
 *   μ      = μ_born − μ_na
 * ─────────────────────────────────────────────────────────────────────────── */
int calculateShearModulus(double *mu_born, double *mu_na) {
    int i, q;
    double *nav, *xi;

    nav = (double *)malloc(sizeof(double) * N * DIM);
    xi  = (double *)malloc(sizeof(double) * N * DIM);

    u_epsilon_x(X_COMP, Y_COMP, xi);
    for (i = 0; i < N; i++)
        for (q = 0; q < DIM; q++) nav[DIM*i+q] = 0.0;

    if (cgSolver(nav, xi)) {
        *mu_na   = cdot(nav, xi) / V;
        *mu_born = (u_epsilon_epsilon(X_COMP,Y_COMP,X_COMP,Y_COMP)
                   + u_epsilon(Y_COMP,Y_COMP)) / V;
        free(nav); free(xi);
        return 1;
    }
    printf("calculateShearModulus: CG failed\n");
    free(nav); free(xi);
    return 0;
}

/* ─────────────────────────────────────────────────────────────────────────────
 * test_network()  —  SOURCE: 3D_relaxed_spring_network_analysis.c
 * Prints min/max per-node coordination numbers as a sanity check.
 * ─────────────────────────────────────────────────────────────────────────── */
void test_network() {
    int i, j, l, *num_of_nebz;
    int mn = 100000, mx = 0;
    num_of_nebz = (int *)calloc(N, sizeof(int));
    for (l = 0; l < contacts; l++) {
        num_of_nebz[lookupBond[2*l+I_INDEX]]++;
        num_of_nebz[lookupBond[2*l+J_INDEX]]++;
    }
    for (i = 0; i < N; i++) {
        if (num_of_nebz[i] < mn) mn = num_of_nebz[i];
        if (num_of_nebz[i] > mx) mx = num_of_nebz[i];
    }
    printf("Network: N=%d  contacts=%d  z_mean=%.4f  z_min=%d  z_max=%d\n",
           N, contacts, 2.0*(double)contacts/DOUBLE_N, mn, mx);
    free(num_of_nebz);
}

/* ─────────────────────────────────────────────────────────────────────────────
 * compression_nonaffine_velocities()
 *   SOURCE: 3dlj_elasticity_analysis.c::compression_nonaffine_velocities()
 *   (adapted for spring networks: first[l]=0 so xi simplifies to Σ rij)
 *
 * Solves H · v = −ξ_η  where ξ_η is the affine driving force under uniform
 * volumetric compression.  The result v is the non-affine velocity field.
 * Returned in res[].
 * ─────────────────────────────────────────────────────────────────────────── */
int compression_nonaffine_velocities(double *res) {
    int a, i, j, l, q;
    double r2;
    double *xi = (double *)malloc(sizeof(double) * N * DIM);

    for (a = 0; a < DIM*N; a++) { res[a] = 0.0; xi[a] = 0.0; }

    for (l = 0; l < contacts; l++) {
        i  = lookupBond[2*l + I_INDEX];
        j  = lookupBond[2*l + J_INDEX];
        r2 = 0.0;
        for (q = 0; q < DIM; q++) r2 += rij[DIM*l+q]*rij[DIM*l+q];
        for (q = 0; q < DIM; q++) {
            xi[DIM*j + q] += (second[l]*r2 + first[l]) * rij[DIM*l+q];
            xi[DIM*i + q] -= (second[l]*r2 + first[l]) * rij[DIM*l+q];
        }
    }

    if (cgSolver(res, xi)) {
        for (a = 0; a < DIM*N; a++) res[a] = -res[a];
        free(xi);
        return 1;
    }
    free(xi);
    return 0;
}

/* ─────────────────────────────────────────────────────────────────────────────
 * write_bonds_text()
 *   SOURCE: 3dlj_elasticity_analysis.c::write_bonds_text()
 *   (adapted: writes i j rij_x rij_y rij_z first second third
 *    for spring networks: first=0 and third=0 since d³U/dr³=0 for k(r-r₀)²)
 *
 * Python notebook uses this to build the sparse Hessian matrix and to
 * compute the mode-sum prediction for (∂P/∂T)_harm.
 * ─────────────────────────────────────────────────────────────────────────── */
void write_bonds_text(int serial_idx) {
    int l;
    char fname[256];
    FILE *f;

    sprintf(fname, "bonds_g%d.txt", serial_idx); /* CHANGED to _g */
    f = fopen(fname, "w");
    fprintf(f, "# i j rij_x rij_y rij_z r first second third\n"); /* CHANGED: added r */
    for (l = 0; l < contacts; l++) {
        fprintf(f, "%d %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
                lookupBond[2*l + I_INDEX],
                lookupBond[2*l + J_INDEX],
                rij[DIM*l + X_COMP],
                rij[DIM*l + Y_COMP],
                rij[DIM*l + Z_COMP],
                1.0 / inv_r[l],  
                first[l],   
                second[l],  
                -3.0 * inv_r[l] * inv_r[l] * inv_r[l] * inv_r[l]);  /* CHANGED: -3/r^4 Cartesian correction */
    }
    fclose(f);
    printf("  wrote %s (%d bonds)\n", fname, contacts);
}

/* ─────────────────────────────────────────────────────────────────────────────
 * write_hessian_text()
 *   SOURCE: 3dlj_elasticity_analysis.c::write_hessian_text()
 *   (adapted for spring networks: first[l]=0 eliminates delta-function terms)
 *
 * Writes the full sparse Hessian in coordinate-list (COO) format:
 *   row   col   value
 * one entry per line.  scipy.sparse.coo_matrix() reads this directly.
 *
 * H_{(i,α),(j,β)} = −second[l] * rij_α * rij_β
 * H_{(i,α),(i,α)} = +Σ_{l∋i} second[l] * rij_α * rij_α  (diagonal)
 * ─────────────────────────────────────────────────────────────────────────── */
void write_hessian_text(int serial_idx) {
    int i, j, l, m, alpha, beta, a, b;
    double thisHes;
    double *diagonal;
    char fname[256];
    FILE *f;

    
    sprintf(fname, "hessian_g%d.txt", serial_idx);
    f = fopen(fname, "w");

    diagonal = (double *)calloc(N * DIM * DIM, sizeof(double));

    /* off-diagonal blocks */
    for (l = 0; l < contacts; l++) {
        i = lookupBond[2*l + I_INDEX];
        j = lookupBond[2*l + J_INDEX];
        m = 0;
        for (alpha = 0; alpha < DIM; alpha++) {
            for (beta = 0; beta < DIM; beta++) {
                /* spring network: first[l]=0  →  no delta(α,β) term          */
                thisHes = -second[l] * rij[DIM*l+alpha] * rij[DIM*l+beta]
                          -first[l]  * DELTA(alpha, beta);
                a = DIM*i + alpha;
                b = DIM*j + beta;
                fprintf(f, "%d\t%d\t%.15g\n%d\t%d\t%.15g\n", a, b, thisHes, b, a, thisHes);
                diagonal[DIM*DIM*i + m] -= thisHes;
                diagonal[DIM*DIM*j + m] -= thisHes;
                m++;
            }
        }
    }
    /* diagonal blocks */
    for (i = 0; i < N; i++) {
        m = 0;
        for (alpha = 0; alpha < DIM; alpha++) {
            for (beta = 0; beta < DIM; beta++) {
                a = DIM*i + alpha;
                b = DIM*i + beta;
                fprintf(f, "%d\t%d\t%.15g\n", a, b, diagonal[DIM*DIM*i + m]);
                m++;
            }
        }
    }
    free(diagonal);
    fclose(f);
    printf("  wrote %s (%d×%d sparse matrix)\n", fname, DIM*N, DIM*N);
}

/* ─────────────────────────────────────────────────────────────────────────────
 * write_vna_eta_text()
 *   SOURCE: 3dlj_elasticity_analysis.c::write_vna_eta_text()
 *   (adapted for spring networks)
 *
 * Solves H · v_na = −ξ_η  and writes the N×3 non-affine velocity field.
 * The Python notebook uses this for the mode-sum correction to (∂P/∂T).
 * ─────────────────────────────────────────────────────────────────────────── */
void write_vna_eta_text(int serial_idx) {
    int i;
    double *vna;
    char fname[256];
    FILE *f;

    sprintf(fname, "vna_eta_g%d.txt", serial_idx);
    vna = (double *)calloc(DIM * N, sizeof(double));

    if (!compression_nonaffine_velocities(vna)) {
        printf("write_vna_eta_text: CG solver failed\n");
        free(vna); return;
    }

    f = fopen(fname, "w");
    fprintf(f, "# vna_eta: compression non-affine velocity (N rows × 3 cols)\n");
    for (i = 0; i < N; i++) {
        fprintf(f, "%.15g\t%.15g\t%.15g\n",
                vna[DIM*i + X_COMP],
                vna[DIM*i + Y_COMP],
                vna[DIM*i + Z_COMP]);
    }
    fclose(f);
    printf("  wrote %s\n", fname);
    free(vna);
}

/* ─────────────────────────────────────────────────────────────────────────────
 * free_everything()  —  SOURCE: 3D_relaxed_spring_network_analysis.c
 * ─────────────────────────────────────────────────────────────────────────── */
void free_everything() {
    free(rx); free(ry); free(rz);
    free(lookupBond);
    free(second); free(first);
    free(rij); free(inv_r);
}

/* =============================================================================
 * main()
 *
 * Workflow for each serial number:
 *  1. Read the spring network
 *  2. Network diagnostic (test_network)
 *  3. Compute static (T=0) pressure P0  — should be 0 for stress-free net
 *  4. Compute bulk modulus K
 *  5. Compute shear modulus μ
 *  6. Write scalar results to snp_elasticity_N4000_z16.24_s<s>.dat
 *  7. Write per-bond data for Python (bonds, hessian, vna_eta text files)
 * =============================================================================
 */
int main(int argc, char **argv) {
    int serial_idx;
    double k_born, k_na, K;
    double mu_born, mu_na, mu;
    char networkFileName[256];
    char outFileName[256];
    FILE *outFile;
    int t0 = (int)time(0);

    if (argc < 2) {
        printf("Usage: ./snp_elast <serial>\n");
        return 1;
    }
    sscanf(argv[1], "%d", &serial_idx);

    linearOperator = &hessianOnVector;

    sprintf(networkFileName, "network_N%d_z16.24_s%.4d.dat", N, serial_idx);
    if (!readNetwork(networkFileName)) {
        printf("ERROR: could not read %s\n", networkFileName);
        return 1;
    }
    printf("Read %s\n", networkFileName);
    test_network();

    /* ── T=0 virial pressure ─────────────────────────────────────────────── */
    /* pressure global = 0 for stress-free spring network (all forces = 0)    */
    printf("T=0 pressure P0 = %.10e  (expected ~0)\n", pressure);

    /* ── bulk modulus ────────────────────────────────────────────────────── */
    printf("Computing bulk modulus K...\n");
    if (!calculateBulkModulus(&k_born, &k_na)) {
        printf("ERROR: calculateBulkModulus failed\n");
        return 1;
    }
    K = (k_born - k_na) / ((double)(DIM * DIM)) + pressure;
    printf("  K_born = %.8f\n  K_na   = %.8f\n  K      = %.8f\n",
           k_born, k_na, K);

    /* ── shear modulus ───────────────────────────────────────────────────── */
    printf("Computing shear modulus mu...\n");
    if (!calculateShearModulus(&mu_born, &mu_na)) {
        printf("ERROR: calculateShearModulus failed\n");
        return 1;
    }
    mu = mu_born - mu_na;
    printf("  mu_born = %.8f\n  mu_na   = %.8f\n  mu      = %.8f\n",
           mu_born, mu_na, mu);

    /* ── write scalar results ────────────────────────────────────────────── */
    sprintf(outFileName, "snp_elasticity_N%d_z16.24_s%.3d.dat", N, serial_idx);
    outFile = fopen(outFileName, "a");   /* append so multiple serials go in one file */
    if (serial_idx == 0)
        fprintf(outFile, "# serial  P0  K_born  K_na  K  mu_born  mu_na  mu\n");
    fprintf(outFile, "%d\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
            serial_idx, pressure, k_born, k_na, K, mu_born, mu_na, mu);
    fclose(outFile);
    printf("Scalar results appended to %s\n", outFileName);

    /* ── write Python-readable auxiliary files ───────────────────────────── */
    printf("Writing auxiliary files for Python notebook...\n");
    write_bonds_text(serial_idx);
    write_hessian_text(serial_idx);
    write_vna_eta_text(serial_idx);

/* --- REPLACE THOSE 5 LINES WITH THIS --- */
    char metaName[256];
    sprintf(metaName, "meta_g%d.txt", serial_idx);
    FILE *fmeta = fopen(metaName, "w");
    fprintf(fmeta, "V\t%.15g\n", V);  /* <-- CHANGED: explicitly labels it 'V' */
    fclose(fmeta);
    /* --------------------------------------- */

    printf("Done. Elapsed: %d s\n", (int)time(0) - t0);
    free_everything();
    return 0;
}
