/* =============================================================================
 * spring_networks_pressure_NVT.c
 *
 * PURPOSE
 * -------
 * Run Berendsen-thermostatted NVT molecular dynamics on a 3D spring network
 * and record the time-averaged virial pressure P_vir(T) as a function of
 * temperature.  The output feeds directly into the Python notebook that
 * computes the thermal expansion coefficient
 *
 *       alpha = (1/K) * dP_vir/dT |_V
 *
 * OUTPUT FILE:  snp_NVT_T<T>_s<serial>.dat
 * COLUMNS:      time   KE/N   PE/N   P_virial
 * (identical format to 3dlj_NVT_NH.c so the same read_nvt() helper works)
 *
 * USAGE:   ./snp_nvt <temperature>  <serial>
 * EXAMPLE: ./snp_nvt 0.05 0
 *
 * COMPILE: gcc -O2 -o snp_nvt spring_networks_pressure_NVT.c -lm
 *
 * =============================================================================
 * SOURCE ATTRIBUTION (every function is labelled)
 *
 *  normal()               — 3dlj_NVT_NH.c  (Box-Muller Maxwell-Boltzmann)
 *  fixDrift()             — spring_networks_finite_T_dynamics.c
 *  readNetwork()          — 3D_relaxed_spring_network_analysis.c  (base)
 *                           + restlength[] addition from
 *                             spring_networks_finite_T_dynamics.c
 *                           + px/py/pz/fx/fy/fz allocation added here
 *  calculateForces()      — spring_networks_finite_T_dynamics.c
 *                           BUG FIX 1: reset u = 0, pressure = 0 at entry
 *                           BUG FIX 2: normalize pressure /= (DIM * V) at exit
 *  advanceTimeBerendsen() — spring_networks_finite_T_dynamics.c
 *                           (correctness now guaranteed by calculateForces fix)
 *  run_equilibrate()      — structure from 3dlj_NVT_NH.c::run_berendsen()
 *  run_production()       — structure from 3dlj_NVT_NH.c::run_NH()
 *  free_everything()      — spring_networks_finite_T_dynamics.c (extended)
 *  main() velocity init   — 3dlj_NVT_NH.c  (MB draw + fixDrift + run)
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

/* ── integrator ──────────────────────────────────────────────────────────── */
#define  TIME_STEP      0.1          /* same as spring_networks_finite_T_dynamics.c */
#define  HALF_TIME_STEP (0.5*TIME_STEP)
#define  TAU_B          1.0          /* Berendsen coupling time (in time-step units) */

/* ── run lengths ─────────────────────────────────────────────────────────── */
#define  EQUIL_STEPS    10000        /* 1000 time units  — thermalisation       */
#define  PROD_STEPS     50000        /* 5000 time units  — production           */
#define  PRINT_EVERY    50           /* print every 50 steps = 5 time units     */

/* ── PBC geometry ────────────────────────────────────────────────────────── */
#define  LENGTH         1.0          /* reduced-coordinate box side = 1         */
#define  HALF_LENGTH    0.5

/* ── index aliases ───────────────────────────────────────────────────────── */
#define  X_COMP  0
#define  Y_COMP  1
#define  Z_COMP  2
#define  I_INDEX 0
#define  J_INDEX 1

/* ── RNG ──────────────────────────────────────────────────────────────────── */
#define  UNI  ((double)rand() / ((double)RAND_MAX + 1.0))

/* ── globals ─────────────────────────────────────────────────────────────── */
double L, invL, V;      /* box geometry                                       */
double T;               /* target temperature                                  */
double u;               /* total potential energy (reset each force call)      */
double pressure;        /* virial pressure P_vir (reset each force call)       */
double kinetic;         /* kinetic energy (updated each step)                  */
int    contacts;        /* number of bonds                                     */

/* ── position / momentum / force arrays ─────────────────────────────────── */
double *rx, *ry, *rz;   /* reduced coordinates ∈ [0,1)                        */
double *px, *py, *pz;   /* momenta (= velocities, mass = 1)                   */
double *fx, *fy, *fz;   /* forces on each particle                             */

/* ── bond network ────────────────────────────────────────────────────────── */
int    *lookupBond;     /* 2*contacts ints: pairs (i,j)                        */
double *restlength;     /* natural lengths r0_l (physical units)               */

/* needed by readNetwork but not by NVT forces; kept for completeness          */
double *rij;            /* physical bond vectors at read time                  */
double *inv_r;
double *first;
double *second;

/* ─────────────────────────────────────────────────────────────────────────────
 * normal()  —  SOURCE: 3dlj_NVT_NH.c
 * Box–Muller transform: returns a Gaussian deviate with given mean/variance.
 * Used to initialise velocities from the Maxwell–Boltzmann distribution.
 * ─────────────────────────────────────────────────────────────────────────── */
double normal(double mean, double variance) {
    double v1, v2, rsq;
    do {
        v1  = 2.0 * UNI - 1.0;
        v2  = 2.0 * UNI - 1.0;
        rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    return sqrt(variance) * (v1 * sqrt(-2.0 * log(rsq) / rsq)) + mean;
}

/* ─────────────────────────────────────────────────────────────────────────────
 * fixDrift()  —  SOURCE: spring_networks_finite_T_dynamics.c
 * Remove the centre-of-mass momentum drift that accumulates over long runs.
 * Called every 16384 steps in the production loop.
 * ─────────────────────────────────────────────────────────────────────────── */
void fixDrift() {
    int i;
    double Px = 0.0, Py = 0.0, Pz = 0.0;
    for (i = 0; i < N; i++) { Px += px[i]; Py += py[i]; Pz += pz[i]; }
    Px /= DOUBLE_N; Py /= DOUBLE_N; Pz /= DOUBLE_N;
    for (i = 0; i < N; i++) { px[i] -= Px; py[i] -= Py; pz[i] -= Pz; }
}

/* ─────────────────────────────────────────────────────────────────────────────
 * readNetwork()
 *   SOURCE (base):      3D_relaxed_spring_network_analysis.c::readNetwork()
 *   SOURCE (additions): spring_networks_finite_T_dynamics.c
 *     — allocates restlength[] so calculateForces() can use actual spring law
 *     — allocates px/py/pz/fx/fy/fz for the dynamics
 *
 * File format (network_N4000_z16_24_s0000.dat):
 *   Line 1:        L   0   0
 *   Lines 2..N+1:  rx[i]  ry[i]  rz[i]   (reduced coords in [0,1))
 *   Line N+2:      contacts  0  0
 *   Lines N+3..:   i  j  spring_const(ignored)
 *
 * restlength[l] = sqrt(r2) is the natural length of each spring (i.e. the
 * spring is at rest at the configuration encoded in the file, so first[l]=0).
 * ─────────────────────────────────────────────────────────────────────────── */
int readNetwork(char *inFileName) {
    int i, j, l, foo;
    double dummy, dx, dy, dz, r2;
    FILE *inFile;

    rx = (double *)malloc(sizeof(double) * N);
    ry = (double *)malloc(sizeof(double) * N);
    rz = (double *)malloc(sizeof(double) * N);
    px = (double *)malloc(sizeof(double) * N);
    py = (double *)malloc(sizeof(double) * N);
    pz = (double *)malloc(sizeof(double) * N);
    fx = (double *)malloc(sizeof(double) * N);
    fy = (double *)malloc(sizeof(double) * N);
    fz = (double *)malloc(sizeof(double) * N);

    inFile = fopen(inFileName, "rb");
    if (!inFile) {
        printf("readNetwork: could not open %s\n", inFileName);
        return 0;
    }

    foo = fscanf(inFile, "%lf %lf %lf", &L, &dummy, &dummy);
    invL = 1.0 / L;
    V    = L * L * L;

    for (i = 0; i < N; i++)
        foo = fscanf(inFile, "%lf %lf %lf", &rx[i], &ry[i], &rz[i]);

    foo = fscanf(inFile, "%d %lf %lf", &contacts, &dummy, &dummy);

    lookupBond  = (int    *)malloc(sizeof(int)    * 2 * contacts);
    restlength  = (double *)malloc(sizeof(double) *     contacts);
    rij         = (double *)malloc(sizeof(double) * DIM * contacts);
    inv_r       = (double *)malloc(sizeof(double) *       contacts);
    first       = (double *)malloc(sizeof(double) *       contacts);
    second      = (double *)malloc(sizeof(double) *       contacts);

    for (l = 0; l < contacts; l++) {
        foo = fscanf(inFile, "%d %d %lf",
                     &lookupBond[2*l + I_INDEX],
                     &lookupBond[2*l + J_INDEX],
                     &dummy);
        i  = lookupBond[2*l + I_INDEX];
        j  = lookupBond[2*l + J_INDEX];
        dx = rx[j] - rx[i];
        dy = ry[j] - ry[i];
        dz = rz[j] - rz[i];
        /* minimum image convention (reduced coordinates) */
        if      (dz >=  HALF_LENGTH) dz -= LENGTH;
        else if (dz < -HALF_LENGTH)  dz += LENGTH;
        if      (dy >=  HALF_LENGTH) dy -= LENGTH;
        else if (dy < -HALF_LENGTH)  dy += LENGTH;
        if      (dx >=  HALF_LENGTH) dx -= LENGTH;
        else if (dx < -HALF_LENGTH)  dx += LENGTH;

        r2 = L*L * (dx*dx + dy*dy + dz*dz);   /* (physical distance)^2 */

        rij[DIM*l + X_COMP] = L * dx;
        rij[DIM*l + Y_COMP] = L * dy;
        rij[DIM*l + Z_COMP] = L * dz;
        inv_r[l]      = 1.0 / sqrt(r2);
        restlength[l] = sqrt(r2);   /* spring at rest → natural length = current dist */
        first[l]      = 0.0;        /* no pre-stress at rest configuration            */
        second[l]     = 1.0 / r2;  /* spring constant / r^2 (for Hessian use)        */
    }
    fclose(inFile);
    return 1;
}

/* ─────────────────────────────────────────────────────────────────────────────
 * calculateForces()
 *   SOURCE: spring_networks_finite_T_dynamics.c
 *   BUG FIX 1: reset u = 0 and pressure = 0 at the start of every call
 *              (the original code omitted this, causing both to accumulate
 *               across calls instead of holding the current-step value)
 *   BUG FIX 2: normalize pressure by (DIM * V) at the end so that
 *              pressure = P_vir = (1/(3V)) Σ_bonds F_l * r_l
 *              (the original code left the raw virial sum un-normalised)
 *
 * Force law (harmonic spring, unit spring constant):
 *   F_l = (r0_l - r_l)   [restoring, in physical units]
 *   P_vir = (1/3V) Σ_bonds (r0_l - r_l) * r_l
 * ─────────────────────────────────────────────────────────────────────────── */
void calculateForces() {
    int l, i, j;
    double dx, dy, dz, r, g, temp;

    u        = 0.0;   /* BUG FIX 1 */
    pressure = 0.0;   /* BUG FIX 1 */
    for (i = 0; i < N; i++) { fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0; }

    for (l = 0; l < contacts; l++) {
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

        r  = L * sqrt(dx*dx + dy*dy + dz*dz);   /* physical distance */
        g  = (restlength[l] - r) / r;            /* g = F / r        */

        temp   = L * dx * g;
        fx[j] += temp;  fx[i] -= temp;
        temp   = L * dy * g;
        fy[j] += temp;  fy[i] -= temp;
        temp   = L * dz * g;
        fz[j] += temp;  fz[i] -= temp;

        u        += 0.5 * (restlength[l] - r) * (restlength[l] - r);
        pressure += g * r * r;   /* contributes (r0-r)*r to virial sum */
    }
    pressure /= ((double)DIM * V);   /* BUG FIX 2: P_vir = Σ/(3V) */
}

/* ─────────────────────────────────────────────────────────────────────────────
 * advanceTimeBerendsen()
 *   SOURCE: spring_networks_finite_T_dynamics.c  (unchanged logic;
 *           correctness depends on calculateForces() being fixed above)
 *
 * Velocity-Verlet + Berendsen thermostat:
 *   1. Half-kick momenta with current forces
 *   2. Advance positions (PBC wrapped in reduced coords)
 *   3. Recalculate forces at new positions
 *   4. Complete momentum half-kick
 *   5. Rescale momenta toward target temperature T (Berendsen coupling TAU_B)
 * ─────────────────────────────────────────────────────────────────────────── */
void advanceTimeBerendsen() {
    int i;
    double temp, dx, dy, dz, inst_T, xi_T;

    /* half-kick + position update */
    for (i = 0; i < N; i++) {
        px[i] += HALF_TIME_STEP * fx[i];
        py[i] += HALF_TIME_STEP * fy[i];
        pz[i] += HALF_TIME_STEP * fz[i];

        dx   = TIME_STEP * px[i];           /* real-space displacement         */
        temp = rx[i] + dx * invL;           /* convert back to reduced coords  */
        if      (temp >= LENGTH) rx[i] = temp - LENGTH;
        else if (temp <  0.0)   rx[i] = temp + LENGTH;
        else                    rx[i] = temp;

        dy   = TIME_STEP * py[i];
        temp = ry[i] + dy * invL;
        if      (temp >= LENGTH) ry[i] = temp - LENGTH;
        else if (temp <  0.0)   ry[i] = temp + LENGTH;
        else                    ry[i] = temp;

        dz   = TIME_STEP * pz[i];
        temp = rz[i] + dz * invL;
        if      (temp >= LENGTH) rz[i] = temp - LENGTH;
        else if (temp <  0.0)   rz[i] = temp + LENGTH;
        else                    rz[i] = temp;
    }

    calculateForces();   /* forces at new positions; also updates u, pressure  */

    /* complete velocity step + measure kinetic energy */
    kinetic = 0.0;
    for (i = 0; i < N; i++) {
        px[i] += HALF_TIME_STEP * fx[i];
        py[i] += HALF_TIME_STEP * fy[i];
        pz[i] += HALF_TIME_STEP * fz[i];
        kinetic += px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
    }
    kinetic = 0.5 * kinetic;

    /* Berendsen rescaling — drives <T_inst> toward T */
    inst_T = 2.0 * kinetic / ((double)DIM * DOUBLE_N);
    xi_T   = sqrt(1.0 + TIME_STEP * (T / inst_T - 1.0) / TAU_B);
    for (i = 0; i < N; i++) { px[i] *= xi_T; py[i] *= xi_T; pz[i] *= xi_T; }
}

/* ─────────────────────────────────────────────────────────────────────────────
 * run_equilibrate()
 *   SOURCE: structure from 3dlj_NVT_NH.c::run_berendsen()
 * Silent thermalisation run (no file output).  Drives the network to T before
 * the production measurement.
 * ─────────────────────────────────────────────────────────────────────────── */
void run_equilibrate() {
    int k;
    printf("  [equilibration: %d steps]\n", EQUIL_STEPS);
    for (k = 0; k < EQUIL_STEPS; k++) {
        advanceTimeBerendsen();
        if (!(k % 2000))
            printf("    step %6d | KE/N = %.6f | PE/N = %.6f | P_vir = %.8f\n",
                   k, kinetic / DOUBLE_N, u / DOUBLE_N, pressure);
        if (!(k % 16384)) fixDrift();
    }
}

/* ─────────────────────────────────────────────────────────────────────────────
 * run_production()
 *   SOURCE: structure from 3dlj_NVT_NH.c::run_NH()
 * Production run: writes (time  KE/N  PE/N  P_virial) every PRINT_EVERY steps.
 * Column format is IDENTICAL to 3dlj_NVT_NH.c so that the same Python
 * read_nvt() helper works without modification.
 * ─────────────────────────────────────────────────────────────────────────── */
void run_production(char *outFileName) {
    int k;
    FILE *outFile;

    outFile = fopen(outFileName, "w");
    if (!outFile) {
        printf("ERROR: cannot open %s for writing\n", outFileName);
        return;
    }
    fprintf(outFile, "# time   KE_per_N   PE_per_N   P_virial\n");

    printf("  [production: %d steps → %s]\n", PROD_STEPS, outFileName);
    for (k = 0; k < PROD_STEPS; k++) {
        advanceTimeBerendsen();
        if (!(k % PRINT_EVERY))
            fprintf(outFile, "%.8g\t%.8g\t%.8g\t%.8g\n",
                    (double)k * TIME_STEP,
                    kinetic / DOUBLE_N,
                    u       / DOUBLE_N,
                    pressure);
        if (!(k % 16384)) fixDrift();
    }
    fclose(outFile);
}

/* ─────────────────────────────────────────────────────────────────────────────
 * free_everything()  —  SOURCE: spring_networks_finite_T_dynamics.c (extended)
 * ─────────────────────────────────────────────────────────────────────────── */
void free_everything() {
    free(rx);  free(ry);  free(rz);
    free(px);  free(py);  free(pz);
    free(fx);  free(fy);  free(fz);
    free(lookupBond);
    free(restlength);
    free(rij); free(inv_r);
    free(first); free(second);
}

/* ─────────────────────────────────────────────────────────────────────────────
 * main()
 *   SOURCE: 3dlj_NVT_NH.c  (argument parsing, MB velocity init, run sequence)
 *
 * Workflow:
 *   1. Read network file
 *   2. Initialise velocities from Maxwell-Boltzmann distribution at T
 *      (uses normal() from 3dlj_NVT_NH.c)
 *   3. Remove centre-of-mass drift
 *   4. Run Berendsen equilibration (EQUIL_STEPS)
 *   5. Run Berendsen production (PROD_STEPS), writing output file
 * ─────────────────────────────────────────────────────────────────────────── */
int main(int argc, char **argv) {
    int i;
    char networkFileName[256];
    char outFileName[256];

    if (argc < 3) {
        printf("Usage: ./snp_nvt <temperature> <serial>\n");
        printf("Example: ./snp_nvt 0.05 0\n");
        return 1;
    }

    sscanf(argv[1], "%lf", &T);
    int serial;
    sscanf(argv[2], "%d",  &serial);

    srand((unsigned int)(time(0)) + (unsigned int)serial);

    sprintf(networkFileName, "network_N%d_z16.24_s%.4d.dat", N, serial);
    if (!readNetwork(networkFileName)) {
        printf("ERROR: could not read %s\n", networkFileName);
        return 1;
    }
    printf("Loaded %s  (N=%d, contacts=%d, L=%.6f, V=%.6f)\n",
           networkFileName, N, contacts, L, V);

    /* First force evaluation to fill fx/fy/fz before the integrator starts  */
    calculateForces();
    printf("T=0 static pressure: P_vir = %.8e  (should be ~0 for stress-free net)\n",
           pressure);

    /* ── Maxwell-Boltzmann velocity initialisation  (from 3dlj_NVT_NH.c) ── */
    for (i = 0; i < N; i++) {
        px[i] = normal(0.0, T);
        py[i] = normal(0.0, T);
        pz[i] = normal(0.0, T);
    }
    fixDrift();
    printf("Velocities initialised at T = %.6f\n", T);

    /* ── equilibration ───────────────────────────────────────────────────── */
    run_equilibrate();
    /* ── production ──────────────────────────────────────────────────────── */   
    /* Change this line to match the LJ naming convention: */
    sprintf(outFileName, "3dlj_NVT_T%.6f_s%d.dat", T, serial); 
    

    run_production(outFileName);

    printf("Done. Output: %s\n", outFileName);
    free_everything();
    return 0;
}
