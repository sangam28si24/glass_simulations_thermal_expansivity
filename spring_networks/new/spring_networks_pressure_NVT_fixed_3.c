
/*
 * spring_networks_finite_T_dynamics.c
 *
 * Berendsen-NVT molecular dynamics for a 3D harmonic spring network.
 *
 * =======================================================================
 * PHYSICAL CONTEXT  (see thermal_expansion_of_glasses_continued.pdf)
 * =======================================================================
 *
 * Model: a disordered amorphous solid (glass) as a random harmonic spring
 * network at finite temperature T.  Every bond (i,j) carries the harmonic
 * spring potential
 *
 *   phi(r) = 0.5 * (r - r0)^2       (spring constant k = 1, mass m = 1)
 *
 * Goal: compute the isobaric thermal expansion coefficient
 *
 *   alpha = (1/V)(dV/dT)_p  =  (1/K)(dp/dT)_V
 *
 * where K is the athermal bulk modulus.  From the PNAS-2017 inherent-state
 * low-temperature expansion (eq. 8 / eq. 21 of the PDF):
 *
 *   2(dp/dT)_V = -(1/(V*d)) * [d^3U/deta dxx] : H^{-1}
 *               + (1/(V*d)) * [d^2U/deta dx]  . H^{-1} . U''' : H^{-1}
 *
 * where d = DIM, eta is the volumetric strain (V ~ e^{d*eta}),
 * H = d^2U/dx^2 is the Hessian (dynamical matrix), and U''' = d^3U/dx^3
 * is the tessian (third-order tensor).  The inverse Hessian is
 *
 *   H^{-1} = sum_omega  (1/omega^2)  Psi_omega x Psi_omega
 *
 * in the normal-mode basis.
 *
 * TWO CONTRIBUTIONS TO alpha  (both scale as omega^{-4} for soft modes):
 *
 *   Harmonic term:   -(1/(2KVd)) * [d^3U/deta dxx] : H^{-1}
 *      -> pressure-field curvature projected on soft modes
 *
 *   Anharmonic term: +(1/(2KVd)) * [d^2U/deta dx]  . H^{-1} . U''' : H^{-1}
 *      -> mismatch forces (xi vector) coupled via the tessian
 *
 * SHORTHAND BOND COEFFICIENTS  (PDF notation, eq. 2nd page):
 *
 *   [1st] = phi'/r
 *   [2nd] = phi''/r^2 - phi'/r^3
 *   [3rd] = phi'''/r^3 - 3*phi''/r^4 + 3*phi'/r^5
 *
 * They satisfy d[1st]/dr = [2nd]*r,  d[2nd]/dr = [3rd]*r.
 *
 * For the harmonic spring phi(r) = 0.5*(r-r0)^2  ->  phi'=r-r0, phi''=1, phi'''=0:
 *
 *   [1st]  = (r-r0)/r                 ->  0        at inherent state (r=r0)
 *   [2nd]  = 1/r^2 - (r-r0)/r^3      ->  1/r0^2   at inherent state
 *   [3rd]  = -3/r^4 + 3*(r-r0)/r^5   ->  -3/r0^4  at inherent state
 *
 * The arrays first[], second[], third[] store these at the inherent-state
 * (equilibrium) configuration for Hessian and tessian operations.
 * Dynamic forces in calculateForces() are recomputed from current positions.
 *
 * WORKFLOW:
 *   1. readNetwork()        -> load bond geometry, init first/second/third
 *   2. initializeVelocities()  -> Maxwell-Boltzmann at temperature T
 *   3. calculateForces()    -> initial forces at inherent state
 *   4. run_berendsen()      -> NVT equilibration
 *   (5. tessian/xi routines available for subsequent dP/dT computation)
 * =======================================================================
 */

#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define N                   4000
#define DOUBLE_N            ((double)N)
#define DIM                 3
#define CGSOLVER_TOL        1e-10

/* --- index / component aliases --- */
#define DELTA(a,b)          ((a)==(b)?1.0:0.0)
#define X_COMP              0
#define Y_COMP              1
#define Z_COMP              2
#define I_INDEX             0
#define J_INDEX             1

/* --- time-stepping --- */
#define TIME_STEP           0.1
#define HALF_TIME_STEP      (0.5*TIME_STEP)
#define TAU_B               1.0         /* Berendsen coupling time constant */

/* --- box geometry (positions stored in [0,1) box units) --- */
#define LENGTH              1.0
#define HALF_LENGTH         0.5
#define TWO_PI              6.28318530717958

/* --- misc --- */
#define UNI                 ((double)rand()/((double)RAND_MAX + 1.0))
#define HUGE_INT            10000000

/* --- run control --- */
#define EQUIL_DURATION      50000.0     /* equilibration time (time units)   */
#define PRINT_INTERVAL      128         /* steps between terminal output      */


/* ---- global function pointer used by the CG solver ---- */
void (*linearOperator)(double *actOn, double *res);

int serial;  /* network file serial number */

/* ---- thermodynamic scalars ---- */
double L;           /* box side length (real units)                      */
double invL;        /* 1/L                                               */
double V;           /* box volume = L^3                                  */
double u;           /* instantaneous potential energy (reset each call)  */
double pressure;    /* unnormalised virial; physical P = pressure/(DIM*V)*/
double kinetic;     /* instantaneous kinetic energy                      */
double T;           /* target temperature (k_B=1, mass=1 units)          */

int contacts;       /* number of bonds in the network                    */

/* ---- per-particle arrays ---- */
double *rx, *ry, *rz;   /* positions (box units, range [0,1))            */
double *px, *py, *pz;   /* momenta   (real units)                        */
double *fx, *fy, *fz;   /* forces    (real units)                        */

/* ---- per-bond arrays (all evaluated at the inherent-state geometry) ---- */
int    *lookupBond;  /* endpoint indices: lookupBond[2*l+I_INDEX/J_INDEX] */
double *rij;         /* bond vectors xij = xj-xi in real units (size contacts*DIM) */
double *inv_r;       /* 1/|rij|                                           */
double *restlength;  /* rest lengths r0 (real units)                      */
double *first;       /* [1st] = phi'/r             (= 0 for springs at eq.)*/
double *second;      /* [2nd] = phi''/r^2-phi'/r^3 (= 1/r0^2 at eq.)     */
double *third;       /* [3rd] = phi'''/r^3-3phi''/r^4+3phi'/r^5 (=-3/r0^4)*/


/* ======================================================================= */
/*  Basic utilities                                                          */
/* ======================================================================= */

/* Remove net center-of-mass momentum. */
void fixDrift(){
    int i;
    double Px=0.0, Py=0.0, Pz=0.0;
    for (i=0; i<N; i++){ Px += px[i]; Py += py[i]; Pz += pz[i]; }
    Px /= DOUBLE_N; Py /= DOUBLE_N; Pz /= DOUBLE_N;
    for (i=0; i<N; i++){ px[i] -= Px; py[i] -= Py; pz[i] -= Pz; }
}

/* Project out the mean displacement from each Cartesian component of a 3N-vector. */
void removeTranslation(double *v){
    int q, i;
    double D;
    for (q=0; q<DIM; q++){
        D = 0.0;
        for (i=0; i<N; i++) D += v[i*DIM + q];
        D /= DOUBLE_N;
        for (i=0; i<N; i++) v[i*DIM + q] -= D;
    }
}

/* Inverse participation ratio of a 3N-vector (mode localisation measure). */
double getParticipationRatio(double *z){
    int i, q;
    double sum=0.0, temp;
    for (i=0; i<N; i++){
        temp = 0.0;
        for (q=0; q<DIM; q++) temp += z[DIM*i+q]*z[DIM*i+q];
        sum += temp*temp;
    }
    return 1.0/(DOUBLE_N * sum);
}

/* Subtract translations and normalize to unit length. */
void normalizeVector(double *v){
    int q, i, a;
    double D, moo=0.0;
    for (q=0; q<DIM; q++){
        D = 0.0;
        for (i=0; i<N; i++) D += v[i*DIM + q];
        D /= DOUBLE_N;
        for (i=0; i<N; i++) v[i*DIM + q] -= D;
    }
    for (a=0; a<DIM*N; a++) moo += v[a]*v[a];
    moo = 1.0/sqrt(moo);
    for (a=0; a<DIM*N; a++) v[a] *= moo;
}

/* Draw a uniformly random unit vector in 3N-space (with translations removed). */
void setRandomDirection(double *v){
    int a;
    for (a=0; a<DIM*N; a++) v[a] = 2.0*UNI - 1.0;
    normalizeVector(v);
}

/* Write a 3N-vector alongside positions to a text file. */
void writeVector(double *vec, char *outFileName){
    FILE *outFile;
    int i;
    outFile = fopen(outFileName,"wb");
    for (i=0; i<N; i++)
        fprintf(outFile,"%.8g\t%.8g\t%.8g\t%.15g\t%.15g\t%.15g\n",
                rx[i], ry[i], rz[i],
                vec[DIM*i+X_COMP], vec[DIM*i+Y_COMP], vec[DIM*i+Z_COMP]);
    fclose(outFile);
}

/*
 * gaussianRandom — Box-Muller transform.
 * Returns one sample from the standard normal N(0,1) using two U(0,1) draws.
 * The guard against u1=0 prevents log(0).
 */
double gaussianRandom(){
    double u1, u2;
    do { u1 = UNI; } while (u1 <= 0.0);
    u2 = UNI;
    return sqrt(-2.0*log(u1)) * cos(TWO_PI*u2);
}


/* ======================================================================= */
/*  Network I/O                                                             */
/* ======================================================================= */

int readNetwork(char *inFileName){
    int i, j, l, foo;
    double dummy, dx, dy, dz, r2;
    FILE *inFile;

    /* --- allocate per-particle arrays --- */
    rx = (double *)malloc(sizeof(double)*N);
    ry = (double *)malloc(sizeof(double)*N);
    rz = (double *)malloc(sizeof(double)*N);
    px = (double *)malloc(sizeof(double)*N);
    py = (double *)malloc(sizeof(double)*N);
    pz = (double *)malloc(sizeof(double)*N);
    fx = (double *)malloc(sizeof(double)*N);
    fy = (double *)malloc(sizeof(double)*N);
    fz = (double *)malloc(sizeof(double)*N);

    inFile = fopen(inFileName,"rb");
    if (inFile){
        /* box side length L; positions are stored in box units [0,1) */
        foo = fscanf(inFile,"%lf %lf %lf",&L, &dummy, &dummy);
        invL = 1.0/L;
        V    = L*L*L;

        /* node positions in box units */
        for (i=0; i<N; i++)
            foo = fscanf(inFile,"%lf %lf %lf",&rx[i],&ry[i],&rz[i]);

        /* number of bonds */
        foo = fscanf(inFile,"%d %lf %lf",&contacts, &dummy, &dummy);

        /* --- allocate per-bond arrays --- */
        lookupBond = (int    *)malloc(sizeof(int   )*2*contacts);
        rij        = (double *)malloc(sizeof(double)*contacts*DIM);
        inv_r      = (double *)malloc(sizeof(double)*contacts);
        restlength = (double *)malloc(sizeof(double)*contacts);
        first      = (double *)malloc(sizeof(double)*contacts);
        second     = (double *)malloc(sizeof(double)*contacts);
        third      = (double *)malloc(sizeof(double)*contacts); /* [3rd] NEW */

        for (l=0; l<contacts; l++){
            foo = fscanf(inFile,"%d %d %lf",
                         &lookupBond[2*l + I_INDEX],
                         &lookupBond[2*l + J_INDEX],
                         &dummy);
            i = lookupBond[2*l + I_INDEX];
            j = lookupBond[2*l + J_INDEX];

            /* minimum-image bond vector (box units -> real via *L) */
            dx = rx[j] - rx[i];
            dy = ry[j] - ry[i];
            dz = rz[j] - rz[i];
            if      (dz >=  HALF_LENGTH) dz -= LENGTH;
            else if (dz < -HALF_LENGTH ) dz += LENGTH;
            if      (dy >=  HALF_LENGTH) dy -= LENGTH;
            else if (dy < -HALF_LENGTH ) dy += LENGTH;
            if      (dx >=  HALF_LENGTH) dx -= LENGTH;
            else if (dx < -HALF_LENGTH ) dx += LENGTH;

            r2 = L*L*(dx*dx + dy*dy + dz*dz);   /* |rij|^2 in real units^2 */
            rij[DIM*l + X_COMP] = L*dx;
            rij[DIM*l + Y_COMP] = L*dy;
            rij[DIM*l + Z_COMP] = L*dz;
            inv_r[l]      = 1.0/sqrt(r2);        /* 1/r0                    */
            restlength[l] = sqrt(r2);            /* r0 (real units)         */

            /*
             * Inherent-state bond coefficients for phi(r) = 0.5*(r-r0)^2
             * at r = r0  (phi' = 0, phi'' = 1, phi''' = 0):
             *
             *   first  = [1st] = phi'/r         = 0
             *   second = [2nd] = phi''/r^2       = 1/r0^2
             *   third  = [3rd] = -3*phi''/r^4    = -3/r0^4  = -3*second^2
             *
             * These are fixed at the inherent-state geometry and used by
             * hessianOnVector / tessianOnTwoVectors / tessianOnThreeVectors.
             */
            first[l]  = 0.0;
            second[l] = 1.0/r2;
            third[l]  = -3.0 * second[l] * second[l]; /* = -3/r0^4 */
        }
        fclose(inFile);
        return 1;
    }
    printf("couldn't read %s\n", inFileName);
    return 0;
}


/* ======================================================================= */
/*  Velocity (momentum) initialisation — Maxwell-Boltzmann                  */
/* ======================================================================= */

/*
 * initializeVelocities
 *
 * Draw momenta from the Maxwell-Boltzmann distribution at temperature T.
 * With mass = k_B = 1, each Cartesian momentum component satisfies
 *
 *   P(p_{i,alpha}) ~ exp( -p^2 / (2T) )
 *
 * so we draw p ~ N(0, sqrt(T)) using Box-Muller.
 *
 * After drawing we:
 *   1. Remove center-of-mass drift  (zero total momentum).
 *   2. Uniformly rescale to achieve exactly the target temperature
 *      via the equipartition relation  T = 2*KE / (DIM*N).
 *
 * This gives a clean starting point for the Berendsen thermostat.
 */
void initializeVelocities(double target_T){
    int i;
    double sigma, ke, actual_T, scale;

    sigma = sqrt(target_T);

    /* sample each component from N(0, sigma) */
    for (i=0; i<N; i++){
        px[i] = sigma * gaussianRandom();
        py[i] = sigma * gaussianRandom();
        pz[i] = sigma * gaussianRandom();
    }

    /* remove center-of-mass drift */
    fixDrift();

    /* rescale so instantaneous kinetic temperature equals target_T exactly */
    ke = 0.0;
    for (i=0; i<N; i++)
        ke += px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
    ke *= 0.5;
    actual_T = 2.0*ke / (DIM * DOUBLE_N);   /* T_inst = 2KE/(d*N) */
    scale    = sqrt(target_T / actual_T);
    for (i=0; i<N; i++){ px[i] *= scale; py[i] *= scale; pz[i] *= scale; }
}


/* ======================================================================= */
/*  Inner product and Hessian  (d^2U/dx^2)                                 */
/* ======================================================================= */

/* Standard 3N inner product. */
double cdot(double *vec1, double *vec2){
    double res = 0.0;
    int i;
    for (i=0; i<DIM*N; i++) res += vec1[i]*vec2[i];
    return res;
}

/*
 * hessianOnVector:  res = H * vec
 *
 * The Hessian is (PDF eq. 11, inherent-state geometry):
 *   H_{(i,a),(j,b)} = -[2nd]*rij_a*rij_b - [1st]*delta_{ab}   (off-diag)
 *   H_{(i,a),(i,b)} = +sum_{bonds of i}([2nd]*rij_a*rij_b + [1st]*delta_{ab})
 *
 * For springs at equilibrium first[] = 0, so only the [2nd] term survives.
 * The result res is the matrix-vector product H|vec>.
 *
 * This is passed as linearOperator to the CG solver to obtain
 *   v_na = H^{-1} . xi_eta    (nonaffine velocity field).
 */
void hessianOnVector(double *vec, double *res){
    int q, i, j, l;
    double uij[DIM], u_dot_r;

    for (q=0; q<DIM; q++)
        for (i=0; i<N; i++)
            res[DIM*i + q] = 0.0;

    for (l=0; l<contacts; l++){
        i = lookupBond[2*l + I_INDEX];
        j = lookupBond[2*l + J_INDEX];
        for (q=0; q<DIM; q++)
            uij[q] = vec[DIM*j+q] - vec[DIM*i+q];
        for (u_dot_r=0.0, q=0; q<DIM; q++)
            u_dot_r += rij[DIM*l+q]*uij[q];
        for (q=0; q<DIM; q++){
            double val = second[l]*u_dot_r*rij[DIM*l+q] + first[l]*uij[q];
            res[DIM*i + q] -= val;
            res[DIM*j + q] += val;
        }
    }
}

/*
 * hessianOnTwoVectors:  returns  u . H . v  (scalar bilinear form)
 *
 * Evaluates <u|H|v> = sum_{l} ([2nd]*(r.uij)*(r.vij) + [1st]*(uij.vij))
 * Used e.g. to compute mode frequencies: omega^2 = <Psi|H|Psi>.
 *
 * Note: local parameter 'u' (double *) shadows the global 'u' (double).
 * This is intentional and harmless — the global is used only outside
 * this function's scope.
 */
double hessianOnTwoVectors(double *u, double *v){
    int l, i, j, q;
    double u_dot_r, v_dot_r, u_dot_v, res=0.0;
    double uij[DIM], vij[DIM];

    for (l=0; l<contacts; l++){
        i = lookupBond[2*l + I_INDEX];
        j = lookupBond[2*l + J_INDEX];
        for (q=0; q<DIM; q++){
            uij[q] = u[DIM*j+q] - u[DIM*i+q];
            vij[q] = v[DIM*j+q] - v[DIM*i+q];
        }
        for (v_dot_r=0.0, u_dot_r=0.0, u_dot_v=0.0, q=0; q<DIM; q++){
            u_dot_r += rij[DIM*l+q]*uij[q];
            v_dot_r += rij[DIM*l+q]*vij[q];
            u_dot_v += uij[q]*vij[q];
        }
        res += second[l]*u_dot_r*v_dot_r + first[l]*u_dot_v;
    }
    return res;
}


/* ======================================================================= */
/*  Tessian  (d^3U/dx^3)                                                   */
/* ======================================================================= */

/*
 * tessianOnTwoVectors:  res = T : u : v   (3N-vector result)
 *
 * Contracts the rank-3 tessian T_{(k,a),(l,b),(m,c)} = d^3U/dx^3 with
 * vectors u and v over its last two indices, yielding a 3N-vector:
 *
 *   res_{(k,a)} = sum_{l,b,m,c} T_{(k,a),(l,b),(m,c)} * u_{(l,b)} * v_{(m,c)}
 *
 * From PDF eq. (12), the bond-(i,j) contribution to res is:
 *
 *   (delta_{jk}-delta_{ik}) * {
 *       [3rd] * rij_a * (rij.uij) * (rij.vij)
 *     + [2nd] * [ uij_a*(rij.vij) + vij_a*(rij.uij) + (uij.vij)*rij_a ]
 *   }
 *
 * In the anharmonic dP/dT formula this is used as:
 *   (d^2U/deta dx) . H^{-1} . (T : H^{-1})
 * via sum_omega (1/omega^2) cdot( xi_eta, cgSolve(tessianOnTwoVectors(Psi,Psi)) )
 */
void tessianOnTwoVectors(double *u, double *v, double *res){
    int l, i, j, q, a;
    double uij[DIM], vij[DIM];
    double u_dot_r, v_dot_r, u_dot_v, val;

    for (a=0; a<DIM*N; a++) res[a] = 0.0;

    for (l=0; l<contacts; l++){
        i = lookupBond[2*l + I_INDEX];
        j = lookupBond[2*l + J_INDEX];
        for (q=0; q<DIM; q++){
            uij[q] = u[DIM*j+q] - u[DIM*i+q];
            vij[q] = v[DIM*j+q] - v[DIM*i+q];
        }
        for (u_dot_r=0.0, v_dot_r=0.0, u_dot_v=0.0, q=0; q<DIM; q++){
            u_dot_r += rij[DIM*l+q]*uij[q];
            v_dot_r += rij[DIM*l+q]*vij[q];
            u_dot_v += uij[q]*vij[q];
        }
        for (q=0; q<DIM; q++){
            val = third[l]  * rij[DIM*l+q] * u_dot_r * v_dot_r
                + second[l] * ( uij[q]*v_dot_r
                              + vij[q]*u_dot_r
                              + u_dot_v*rij[DIM*l+q] );
            res[DIM*i + q] -= val;
            res[DIM*j + q] += val;
        }
    }
}

/*
 * tessianOnThreeVectors:  returns  T : u : v : w  (fully contracted scalar)
 *
 *   result = sum_{k,a,l,b,m,c} T_{(k,a),(l,b),(m,c)} * u_{k,a} * v_{l,b} * w_{m,c}
 *
 * Bond (i,j) contribution:
 *   [3rd]*(rij.uij)*(rij.vij)*(rij.wij)
 * + [2nd]*[ (rij.uij)*(vij.wij) + (rij.vij)*(uij.wij) + (rij.wij)*(uij.vij) ]
 *
 * Used to build U''' : H^{-1} mode by mode:
 *   sum_omega (1/omega^2) * tessianOnThreeVectors(Psi_omega, Psi_omega, xi)
 */
double tessianOnThreeVectors(double *u, double *v, double *w){
    int l, i, j, q;
    double uij[DIM], vij[DIM], wij[DIM];
    double u_dot_r, v_dot_r, w_dot_r;
    double u_dot_v, u_dot_w, v_dot_w;
    double res = 0.0;

    for (l=0; l<contacts; l++){
        i = lookupBond[2*l + I_INDEX];
        j = lookupBond[2*l + J_INDEX];
        for (q=0; q<DIM; q++){
            uij[q] = u[DIM*j+q] - u[DIM*i+q];
            vij[q] = v[DIM*j+q] - v[DIM*i+q];
            wij[q] = w[DIM*j+q] - w[DIM*i+q];
        }
        for (u_dot_r=0.0, v_dot_r=0.0, w_dot_r=0.0,
             u_dot_v=0.0, u_dot_w=0.0, v_dot_w=0.0, q=0; q<DIM; q++){
            u_dot_r += rij[DIM*l+q]*uij[q];
            v_dot_r += rij[DIM*l+q]*vij[q];
            w_dot_r += rij[DIM*l+q]*wij[q];
            u_dot_v += uij[q]*vij[q];
            u_dot_w += uij[q]*wij[q];
            v_dot_w += vij[q]*wij[q];
        }
        res += third[l]  * u_dot_r * v_dot_r * w_dot_r
             + second[l] * ( u_dot_r*v_dot_w
                           + v_dot_r*u_dot_w
                           + w_dot_r*u_dot_v );
    }
    return res;
}


/* ======================================================================= */
/*  Mixed eta-x derivatives: pressure gradient / pressure curvature         */
/* ======================================================================= */

/*
 * xiEtaOnVector:  res = d^2U / (deta dx)     ("xi" vector, PDF eq. 14)
 *
 * This 3N-vector couples volumetric strain eta to atomic displacements:
 *
 *   res_{(k,a)} = sum_{i<j} (delta_{jk}-delta_{ik})
 *                   * ([2nd]*r^2_ij + [1st]) * rij_a
 *
 * At the spring inherent state ([1st]=0, [2nd]*r^2=1):
 *   res_{(k,a)} = sum_{bonds on k}  (+-) rij_a
 *
 * Physical meaning: xi is the "affine force" that drives nonaffine
 * displacements when the box is strained.  The nonaffine velocity field
 * (needed for the anharmonic dP/dT term) is obtained via:
 *
 *   v_na = -H^{-1} . xi      (solve with cgSolver)
 *
 * This is also -(V*d) * dp/dx  (the pressure gradient, PDF eq. 19).
 */
void xiEtaOnVector(double *res){
    int l, i, j, q, a;
    double r2_bond, coeff;

    for (a=0; a<DIM*N; a++) res[a] = 0.0;

    for (l=0; l<contacts; l++){
        i = lookupBond[2*l + I_INDEX];
        j = lookupBond[2*l + J_INDEX];
        r2_bond = 0.0;
        for (q=0; q<DIM; q++) r2_bond += rij[DIM*l+q]*rij[DIM*l+q];

        /* [2nd]*r^2 + [1st]  (= 1 for springs at equilibrium) */
        coeff = second[l]*r2_bond + first[l];

        for (q=0; q<DIM; q++){
            res[DIM*j + q] += coeff * rij[DIM*l+q];
            res[DIM*i + q] -= coeff * rij[DIM*l+q];
        }
    }
}

/*
 * pressureCurvatureOnTwoVectors:  returns  (d^3U/deta dxx) : u : v
 *                                          (PDF eq. 15, bilinear scalar)
 *
 * Needed for the harmonic contribution to dP/dT:
 *   -(1/(2KVd)) * sum_omega (1/omega^2)
 *                * pressureCurvatureOnTwoVectors(Psi_omega, Psi_omega)
 *
 * Bond (i,j) contribution:
 *   ([3rd]*r^2 + 2[2nd]) * (rij.uij) * (rij.vij)
 * + [2nd]*r^2 * (uij.vij)
 *
 * At the spring inherent state:
 *   [3rd]*r^2 + 2[2nd] = -3[2nd] + 2[2nd] = -[2nd]
 *   [2nd]*r^2          = 1
 * so the bond contribution simplifies to:
 *   (uij.vij) - [2nd]*(rij.uij)*(rij.vij)   <- "transverse" dot product
 */
double pressureCurvatureOnTwoVectors(double *u, double *v){
    int l, i, j, q;
    double uij[DIM], vij[DIM];
    double u_dot_r, v_dot_r, u_dot_v, r2_bond, coeff1, coeff2, res=0.0;

    for (l=0; l<contacts; l++){
        i = lookupBond[2*l + I_INDEX];
        j = lookupBond[2*l + J_INDEX];
        for (q=0; q<DIM; q++){
            uij[q] = u[DIM*j+q] - u[DIM*i+q];
            vij[q] = v[DIM*j+q] - v[DIM*i+q];
        }
        u_dot_r = v_dot_r = u_dot_v = r2_bond = 0.0;
        for (q=0; q<DIM; q++){
            u_dot_r += rij[DIM*l+q]*uij[q];
            v_dot_r += rij[DIM*l+q]*vij[q];
            u_dot_v += uij[q]*vij[q];
            r2_bond += rij[DIM*l+q]*rij[DIM*l+q];
        }
        coeff1 = third[l]*r2_bond + 2.0*second[l]; /* [3rd]*r^2 + 2[2nd] */
        coeff2 = second[l]*r2_bond;                  /* [2nd]*r^2           */
        res += coeff1*u_dot_r*v_dot_r + coeff2*u_dot_v;
    }
    return res;
}


/* ======================================================================= */
/*  CG Solver:  linearOperator(x) = b                                      */
/* ======================================================================= */

/*
 * cgSolver:  Conjugate-gradient iterative solver.
 *
 * Solves H.x = b (or whatever linearOperator is set to) to tolerance
 * CGSOLVER_TOL on the RMS residual relative to the RMS of b.
 *
 * Primary use: obtain the nonaffine velocity field
 *
 *   v_na = H^{-1} . xi_eta
 *
 * by calling  cgSolver(v_na, xi_eta)  after setting
 *   linearOperator = &hessianOnVector.
 *
 * Translations are projected out from the solution at convergence.
 * Returns 1 on success, 0 if the iteration limit was reached.
 */
int cgSolver(double *x, double *b){
    double *r, *w, *z;
    double a, B, sum, scale;
    double maxNorm, typicalNorm;
    int i, k, q;

    scale = sqrt(cdot(b,b)/((double)DIM*DOUBLE_N));

    r = (double *)malloc(sizeof(double)*DIM*N);
    w = (double *)malloc(sizeof(double)*DIM*N);
    z = (double *)malloc(sizeof(double)*DIM*N);

    linearOperator(x,r);
    for (i=0; i<DIM*N; i++){
        r[i] = b[i] - r[i];
        w[i] = -r[i];
    }
    linearOperator(w,z);
    a = cdot(r,w)/cdot(w,z);
    for (i=0; i<DIM*N; i++) x[i] += a*w[i];
    B = 0.0;

    for (k=0; k<HUGE_INT; k++){
        maxNorm = typicalNorm = 0.0;
        for (i=0; i<DIM*N; i++){
            r[i] -= a*z[i];
            if (r[i]*r[i] > maxNorm) maxNorm = r[i]*r[i];
            typicalNorm += r[i]*r[i];
        }
        maxNorm     = sqrt(maxNorm);
        typicalNorm = sqrt(typicalNorm/((double)DIM*DOUBLE_N));
        if (typicalNorm/scale < CGSOLVER_TOL) break;

        B = cdot(r,z)/cdot(w,z);
        for (i=0; i<DIM*N; i++) w[i] = -r[i] + B*w[i];
        linearOperator(w,z);
        a = cdot(r,w)/cdot(w,z);
        for (i=0; i<DIM*N; i++) x[i] += a*w[i];
    }

    /* project out translations from the solution */
    for (q=0; q<DIM; q++){
        sum = 0.0;
        for (i=0; i<N; i++) sum += x[DIM*i+q];
        sum /= DOUBLE_N;
        for (i=0; i<N; i++) x[DIM*i+q] -= sum;
    }
    free(r); free(w); free(z);
    return (k < HUGE_INT) ? 1 : 0;
}


/* ======================================================================= */
/*  Harmonic (normal-mode) dynamics                                         */
/* ======================================================================= */

/*
 * run_harmonic_dynamics
 *
 * Integrates  z'' = -H.z  (velocity Verlet) starting from
 *   z(0) = 0,   z'(0) = initial_velocity
 * for 'duration' time units.
 *
 * Records the velocity-velocity autocorrelation:
 *   c[q] = <v(t).v(0)> / (|v(t)||v(0)|)
 *
 * The Fourier transform of c(t) gives the projected density of states
 * weighted by the initial-velocity mode occupancy — useful for checking
 * the vibrational spectrum of the network.
 *
 * The array c must be pre-allocated with at least (int)(duration/TIME_STEP)
 * elements.
 *
 * Note: this evolves in the harmonic well around the inherent state using
 * the stored rij/second/first, NOT the full anharmonic calculateForces().
 */
void run_harmonic_dynamics(double *initial_velocity, double duration, double *c){
    double *z, *p, *hesOnZ;
    int a, q, total_run_steps;

    total_run_steps = (int)(duration/TIME_STEP);
    z      = (double *)malloc(sizeof(double)*DIM*N);
    p      = (double *)malloc(sizeof(double)*DIM*N);
    hesOnZ = (double *)malloc(sizeof(double)*DIM*N);

    for (a=0; a<DIM*N; a++){
        p[a] = initial_velocity[a];
        z[a] = 0.0;
    }
    hessianOnVector(z, hesOnZ);

    c[0] = 1.0;
    for (q=1; q<total_run_steps; q++){
        for (a=0; a<DIM*N; a++){
            p[a] -= HALF_TIME_STEP * hesOnZ[a];
            z[a] += TIME_STEP * p[a];
        }
        hessianOnVector(z, hesOnZ);
        for (a=0; a<DIM*N; a++)
            p[a] -= HALF_TIME_STEP * hesOnZ[a];

        /* normalised velocity-velocity autocorrelation */
        c[q] = cdot(p, initial_velocity) /
               (sqrt(cdot(initial_velocity, initial_velocity)) * sqrt(cdot(p,p)));
    }
    free(z); free(p); free(hesOnZ);
}


/* ======================================================================= */
/*  Force calculation and NVT integrator                                    */
/* ======================================================================= */

/*
 * calculateForces
 *
 * Computes forces, potential energy u, and virial from CURRENT positions.
 * Uses the full spring potential phi(r) = 0.5*(r-r0)^2 with the stored
 * restlength[]; does NOT use first/second/third (those are inherent-state).
 *
 * Force on particle j from bond l (i,j):
 *   F_{j,a} = (r0-r)/r * rij_a   (harmonic restoring, in real units)
 *
 * Virial:
 *   pressure += g*r^2 = (r0-r)*r     (un-normalised; physical P = pressure/(DIM*V))
 *
 * IMPORTANT: u and pressure are reset to zero at the start of each call.
 * Earlier versions accumulated across calls, which was a bug — fixed here.
 */
void calculateForces(){
    int l, i, j;
    double dx, dy, dz, r, g, temp;

    /* reset accumulators each call */
    u        = 0.0;
    pressure = 0.0;

    for (i=0; i<N; i++){ fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0; }

    for (l=0; l<contacts; l++){
        i = lookupBond[2*l + I_INDEX];
        j = lookupBond[2*l + J_INDEX];

        /* minimum-image bond vector from current (box-unit) positions */
        dx = rx[j] - rx[i];
        dy = ry[j] - ry[i];
        dz = rz[j] - rz[i];
        if      (dz >=  HALF_LENGTH) dz -= LENGTH;
        else if (dz < -HALF_LENGTH ) dz += LENGTH;
        if      (dy >=  HALF_LENGTH) dy -= LENGTH;
        else if (dy < -HALF_LENGTH ) dy += LENGTH;
        if      (dx >=  HALF_LENGTH) dx -= LENGTH;
        else if (dx < -HALF_LENGTH ) dx += LENGTH;

        r = L * sqrt(dx*dx + dy*dy + dz*dz);   /* actual bond length (real) */
        g = (restlength[l] - r) / r;            /* = -[1st] at current r    */

        temp = L*dx*g; fx[j] += temp; fx[i] -= temp;
        temp = L*dy*g; fy[j] += temp; fy[i] -= temp;
        temp = L*dz*g; fz[j] += temp; fz[i] -= temp;

        u        += 0.5*(restlength[l]-r)*(restlength[l]-r); /* phi(r)     */
        pressure += g*r*r;    /* virial contribution: F*r = (r0-r)*r        */
    }
}

/*
 * advanceTimeBerendsen
 *
 * One velocity-Verlet step with Berendsen NVT rescaling.
 *
 * Algorithm:
 *   1. Half-kick:  p += (dt/2)*F
 *   2. Drift:      x += dt*p   (with PBC, positions in box units)
 *   3. Recompute forces at new positions.
 *   4. Half-kick:  p += (dt/2)*F
 *   5. Berendsen rescale:  p *= xi_T = sqrt(1 + dt/tau*(T/T_inst - 1))
 *
 * T_inst = 2*KE/(d*N).  The rescaling weakly couples to a heat bath with
 * coupling time TAU_B.  Not canonical-ensemble exact, but ideal for fast
 * equilibration.
 */
void advanceTimeBerendsen(){
    int i;
    double temp, dx, dy, dz, inst_T, xi_T;

    /* half-kick + drift */
    for (i=0; i<N; i++){
        px[i] += HALF_TIME_STEP*fx[i];
        py[i] += HALF_TIME_STEP*fy[i];
        pz[i] += HALF_TIME_STEP*fz[i];

        /* x displacement (real) -> box-unit shift */
        dx   = TIME_STEP*px[i];
        temp = rx[i] + dx*invL;
        if      (temp >= LENGTH) rx[i] = temp - LENGTH;
        else if (temp  < 0.0   ) rx[i] = temp + LENGTH;
        else                     rx[i] = temp;

        dy   = TIME_STEP*py[i];
        temp = ry[i] + dy*invL;
        if      (temp >= LENGTH) ry[i] = temp - LENGTH;
        else if (temp  < 0.0   ) ry[i] = temp + LENGTH;
        else                     ry[i] = temp;

        dz   = TIME_STEP*pz[i];
        temp = rz[i] + dz*invL;
        if      (temp >= LENGTH) rz[i] = temp - LENGTH;
        else if (temp  < 0.0   ) rz[i] = temp + LENGTH;
        else                     rz[i] = temp;
    }

    calculateForces();

    /* second half-kick + accumulate kinetic energy */
    kinetic = 0.0;
    for (i=0; i<N; i++){
        px[i] += HALF_TIME_STEP*fx[i];
        py[i] += HALF_TIME_STEP*fy[i];
        pz[i] += HALF_TIME_STEP*fz[i];
        kinetic += px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
    }
    kinetic *= 0.5;

    /* Berendsen thermostat rescaling */
    inst_T = 2.0*kinetic/(DIM*DOUBLE_N);
    xi_T   = sqrt(1.0 + TIME_STEP*(T/inst_T - 1.0)/TAU_B);
    for (i=0; i<N; i++){ px[i] *= xi_T; py[i] *= xi_T; pz[i] *= xi_T; }
}

/*
 * run_berendsen
 *
 * Drives the system at target temperature T for 'duration' time units.
 * Prints  (t  KE/N  PE/N)  every PRINT_INTERVAL steps to stdout.
 * Removes center-of-mass drift every 16384 steps.
 */
void run_berendsen(double duration){
    int k, steps;
    steps = (int)(duration/TIME_STEP);
    for (k=0; k<steps; k++){
        advanceTimeBerendsen();
        if (!(k % PRINT_INTERVAL))
            printf("%.8g\t%.8g\t%.8g\n",
                   (double)k*TIME_STEP,
                   kinetic/DOUBLE_N,
                   u/DOUBLE_N);
        if (!(k % 16384))
            fixDrift();
    }
}


/* ======================================================================= */
/*  Memory management                                                       */
/* ======================================================================= */

void free_everything(){
    free(rx); free(ry); free(rz);
    free(px); free(py); free(pz);
    free(fx); free(fy); free(fz);
    free(lookupBond);
    free(first); free(second); free(third); /* third is new */
    free(rij); free(inv_r); free(restlength);
}


/* ======================================================================= */
/*  Main                                                                    */
/* ======================================================================= */

/*
 * Usage:  ./prog <serial> <temperature>
 *
 *   serial      : network file index (int), builds filename
 *                 network_N<N>_z<coord>_s<serial>.dat
 *   temperature : target temperature T  (k_B=1, mass=1, spring constant k=1)
 *                 Typical range for glass simulations: 0.001 to 0.1
 *
 * Workflow:
 *   1. Parse arguments.
 *   2. readNetwork()          — load geometry; sets first/second/third.
 *   3. initializeVelocities() — Maxwell-Boltzmann at target T.
 *   4. calculateForces()      — initial forces at inherent-state config.
 *   5. run_berendsen()        — NVT equilibration.
 *
 * After equilibration, the arrays rij/first/second/third remain valid for
 * inherent-state Hessian and tessian computations (dP/dT analysis).
 */
int main(int n, char **inputStrings){
    char networkFileName[1024];
    double coordination;

    srand(time(0));

    if (n < 3){
        printf("Usage: %s <serial> <temperature>\n", inputStrings[0]);
        printf("  serial      : network serial number (integer)\n");
        printf("  temperature : target temperature T (k_B=mass=k_spring=1 units)\n");
        return 1;
    }

    sscanf(inputStrings[1], "%d",  &serial);
    sscanf(inputStrings[2], "%lf", &T);

    linearOperator = &hessianOnVector;
    coordination   = 16.24;

    sprintf(networkFileName, "network_N%d_z%.2f_s%.4d.dat", N, coordination, serial);

    if (readNetwork(networkFileName)){
        printf("# Network: %s\n", networkFileName);
        printf("# N=%d  contacts=%d  z=%.4f  L=%.6g  V=%.6g\n",
               N, contacts, 2.0*(double)contacts/DOUBLE_N, L, V);
        printf("# Target T=%.4g  TIME_STEP=%.3g  TAU_B=%.3g\n",
               T, (double)TIME_STEP, (double)TAU_B);

        /* --- Step 1: Maxwell-Boltzmann velocity initialisation --- */
        initializeVelocities(T);
        printf("# Velocities initialized from Maxwell-Boltzmann at T=%.4g\n", T);

        /* --- Step 2: Initial forces at (near-)inherent-state geometry --- */
        calculateForces();
        printf("# Initial state: u/N=%.8g  P=%.8g  (P = pressure/(DIM*V))\n",
               u/DOUBLE_N, pressure/(DIM*V));

        /* --- Step 3: Berendsen NVT equilibration --- */
        printf("# Berendsen equilibration for %.1f time units\n",
               (double)EQUIL_DURATION);
        printf("# t\t\tKE/N\t\tPE/N\n");
        run_berendsen(EQUIL_DURATION);

        printf("# Equilibration complete. Final: KE/N=%.8g  PE/N=%.8g  P=%.8g\n",
               kinetic/DOUBLE_N, u/DOUBLE_N, pressure/(DIM*V));

        /*
         * At this point:
         *   - rij[], first[], second[], third[] = inherent-state geometry
         *     -> ready for hessianOnVector, tessianOnTwoVectors, xiEtaOnVector,
         *        pressureCurvatureOnTwoVectors (for dP/dT computation)
         *   - rx/ry/rz, px/py/pz = equilibrated at temperature T
         *
         * Add production-run / analysis code here as needed.
         */

        free_everything();
    }
    else {
        printf("%s not found\n", networkFileName);
        return 1;
    }

    return 0;
}