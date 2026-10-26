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
#define  TIME_STEP      0.1          
#define  HALF_TIME_STEP (0.5*TIME_STEP)
#define  TAU_B          1.0         

/* ── run lengths ─────────────────────────────────────────────────────────── */
#define  EQUIL_STEPS    10000        
#define  PROD_STEPS     50000        
#define  PRINT_EVERY    50          

/* ── PBC geometry ────────────────────────────────────────────────────────── */
#define  LENGTH         1.0          
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
double L, invL, V;      
double T;               
double u;              
double pressure = 0.0; 
double kinetic;         
int    contacts;       

/* ── position / momentum / force arrays ─────────────────────────────────── */
double *rx, *ry, *rz;  
double *px, *py, *pz;   
double *fx, *fy, *fz;  

/* ── bond network ────────────────────────────────────────────────────────── */
int    *lookupBond;     
double *restlength;     

/* needed by readNetwork but not by NVT forces; kept for completeness          */
double *rij;            
double *inv_r;
double *first;
double *second;
double *third;          

/* ─────────────────────────────────────────────────────────────────────────────
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
    third       = (double *)malloc(sizeof(double) *       contacts);

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

            r2 = L*L*( dx*dx + dy*dy + dz*dz );
			rij[DIM*l + X_COMP] = L*dx;
			rij[DIM*l + Y_COMP] = L*dy;
			rij[DIM*l + Z_COMP] = L*dz;
			
			/* General generalized spatial derivatives (U >= 0) */
			double r = sqrt(r2);
			double r0 = r; /* Network starts at U=0, so r0 equals the initial distance */
			restlength[l] = r;
            
			inv_r[l]  = 1.0 / r;
			first[l]  = 1.0 - (r0 / r);
			second[l] = r0 / (r2 * r);
			third[l]  = -3.0 * r0 / (r2 * r2 * r);
		}
		fclose(inFile);
    return 1;
}

/* ─────────────────────────────────────────────────────────────────────────────
 * calculateForces()
 *reset u = 0 and pressure = 0 at the start of every call
 *normalize pressure by (DIM * V) at the end so that
 *              pressure = P_vir = (1/(3V)) Σ_bonds 
 * Force law (harmonic spring, unit spring constant):
 *   F_l = (r0_l - r_l)   [restoring, in physical units]
 *   P_vir = (1/3V) Σ_bonds (r0_l - r_l) * r0_l
 *   = -(1/3V) Σ_bonds r0_l * φ'(r_l)
 *   The IS distance r0_l = restlength[l] is the geometric prefactor; only
 *   the force (r0-r) responds to thermal displacements.  This causes only
 *   transverse mode projections to contribute to dP/dT (see file header).
 *
 *   Code:  pressure += g * r * restlength[l]
 *          where g = (r0-r)/r, so g*r = r0-r, and g*r*r0 = (r0-r)*r0. 
 *
 * Compare ROUTE 2 (spring_networks_pressure_NVT_fixed.c):
 *   P_vir = (1/3V) Σ_bonds (r0_l - r_l) * r_l   [current r as prefactor]
 *   Code:  pressure += g * r * r
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
        pressure += g * r * r;  /* Route 1: (r0-r)*r0 — IS prefactor */
    }
       /* BUG FIX 2: P_vir = Σ r0*(r0-r) / (3V) */
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
    for (k = 0; k < EQUIL_STEPS; k++) {
        advanceTimeBerendsen();
        if (!(k % 2000))
            printf("    step %6d | KE/N = %.6f | PE/N = %.6f | P_vir = %.8f\n",
                   k, kinetic / DOUBLE_N, u / DOUBLE_N, pressure/(DIM*V));
        if (!(k % 16384)) fixDrift();
    }
}

void run_production(char *outFileName) {
    int k;
    FILE *outFile;

    outFile = fopen(outFileName, "w");
    if (!outFile) {
        printf("ERROR: cannot open %s for writing\n", outFileName);
        return;
    }
    fprintf(outFile, "# time   KE_per_N   PE_per_N   P_virial\n");

    for (k = 0; k < PROD_STEPS; k++) {
        advanceTimeBerendsen();
        if (!(k % PRINT_EVERY))
            fprintf(outFile, "%.8g\t%.8g\t%.8g\t%.8g\n",
                    (double)k * TIME_STEP,
                    kinetic / DOUBLE_N,
                    u       / DOUBLE_N,
                    pressure/(DIM*V)); /* physical pressure */
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
    free(first); free(second); free(third);
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
        printf("Usage: %s <temperature> <serial>\n", argv[0]);
        return 1;
    }

    sscanf(argv[1], "%lf", &T);
    int serial;
    sscanf(argv[2], "%d",  &serial);

    srand((unsigned int)(time(0)) + (unsigned int)serial);

    sprintf(networkFileName, "network_N%d_z16.24_s%.4d.dat", N, serial);
    
    if (readNetwork(networkFileName)) {
        printf("# Network: %s\n", networkFileName);
        printf("# N=%d  contacts=%d  z=%.4f  L=%.6g  V=%.6g\n",
               N, contacts, 2.0*(double)contacts/DOUBLE_N, L, V);
        printf("# Target T=%.4g  TIME_STEP=%.3g  TAU_B=%.3g\n",
               T, (double)TIME_STEP, (double)TAU_B);

        /* Maxwell-Boltzmann velocity initialisation */
        for (i = 0; i < N; i++) {
            px[i] = normal(0.0, T);
            py[i] = normal(0.0, T);
            pz[i] = normal(0.0, T);
        }
        fixDrift();
        printf("# Velocities initialized from Maxwell-Boltzmann at T=%.4g\n", T);

        calculateForces();
        printf("# Initial state: u/N=%.8g  P=%.8g  (P = pressure/(DIM*V))\n",
               u/DOUBLE_N, pressure/(DIM*V));

        printf("# Running equilibration for %d steps...\n", EQUIL_STEPS);
        run_equilibrate();

        sprintf(outFileName, "snp_NVT_T%.6f_s%d.dat", T, serial);
        printf("# Running production for %d steps. Output: %s\n", PROD_STEPS, outFileName);
        run_production(outFileName);
        
        printf("# Production complete. Final: KE/N=%.8g  PE/N=%.8g  P=%.8g\n",
               kinetic/DOUBLE_N, u/DOUBLE_N, pressure/(DIM*V));

        free_everything();
    } else {
        printf("ERROR: %s not found\n", networkFileName);
        return 1;
    }

    return 0;
}
