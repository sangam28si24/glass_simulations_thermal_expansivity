/* ============================================================================
 * 3D_spring_network_dPdT_test_v4.c
 *
 * Built from 3D_spring_network_dPdT_test_v3.c as the primary base.
 * The following functions / globals were imported from
 * spring_networks_elasticity.c (which itself sourced them from
 * 3D_relaxed_spring_network_analysis.c) because they were absent in v2:
 *
 *   GLOBALS ADDED
 *     double pressure          — T=0 virial pressure (0 for stress-free net)
 *     void (*linearOperator)() — function pointer used by the full cgSolver
 *
 *   FUNCTIONS ADDED
 *     test_network()           — coordination-number sanity check
 *     u_epsilon()              — ∂U/∂ε_αβ  (Born stress tensor component)
 *     u_epsilon_epsilon()      — ∂²U/∂ε_αβ∂ε_κχ  (Born elastic tensor)
 *     u_epsilon_x()            — ∂²U/∂ε_αβ∂u_i^q  (non-affine driving force)
 *     calculateBulkModulus()   — K_born, K_na, K = (K_born-K_na)/DIM² + P0
 *     calculateShearModulus()  — μ_born, μ_na, μ = μ_born - μ_na
 *
 *   FUNCTIONS MODIFIED
 *     cgSolver()               — now uses the linearOperator function pointer
 *                                (identical logic; removes hard-wired call to
 *                                 hessianOnVector so future callers can swap it)
 *     write_bonds_text()       — adds physical distance r column; third[]
 *                                computed as -3/r^4 (Cartesian correction)
 *                                matching spring_networks_elasticity.c
 *
 *   SCALAR OUTPUT FILE ADDED
 *     snp_elasticity_N<N>_z16.24_s<idx>.dat
 *       columns: serial  P0  K_born  K_na  K  mu_born  mu_na  mu
 *
 *   UNCHANGED FROM v2
 *     xi_eta_x(), u_eta_x_x_on_two_vectors()
 *     tessianOnThreeVectors(), tessianOnTwoVectors(),
 *     tessianOnTwoVectorsThenDot()
 *     compression_nonaffine_velocities()
 *     write_meta_text(), write_hessian_text(), write_vna_eta_text()
 *     check_one_pass_vs_two_pass_term2()
 *     readNetwork()  (third[] = -3*second/r² already present in v2)
 *
 * -- v4 ADDITIONS ----------------------------------------------------------
 *
 *   #define DOFS (DIM*N)
 *
 *   FUNCTIONS ADDED
 *     mode_terms()            per-mode (term1,term2) for Eq.(20)/(21):
 *                                term1 = u_eta_x_x_on_two_vectors(Psi,Psi)
 *                                term2 = cdot(vna_eta,
 *                                          tessianOnTwoVectors(Psi,Psi))
 *                              term2 is built from the TWO-PASS primitives
 *                              ("tessian on two vectors, then cdot"), i.e.
 *                              algebraically tessianOnThreeVectors(vna_eta,
 *                              Psi,Psi) via the one-pass/two-pass equivalence
 *                              already checked by
 *                              check_one_pass_vs_two_pass_term2(), but never
 *                              evaluated via the ONE-PASS three-vector
 *                              primitive here. This replaces the Python-side
 *                              einsum mirror (tess3_batch) of
 *                              tessianOnThreeVectors used previously.
 *
 *     write_mode_terms_text() reads eigvecs_g{idx}.bin + eigvals_g{idx}.txt
 *                              (written by the notebook after
 *                              scipy.linalg.eigh on hessian_g{idx}.txt: the
 *                              former is n_modes back-to-back blocks of DOFS
 *                              float64 each -- mode m's eigenvector --, the
 *                              latter is n_modes lines of omega^2, same
 *                              order). Calls mode_terms() per mode, writes
 *                              mode_terms_g{idx}.txt (per-mode term1,term2,
 *                              unweighted) and dPdT_results_g{idx}.txt
 *                              (harmonic, anharmonic, 2*(dp/dT)_V [Eq.21
 *                              LHS], (dp/dT)_V [Eq.22 LHS], V).
 *
 *   main() ADDITIONS
 *     A 3rd CLI argument (eigvecs_g{idx}.bin) selects a MODE-TERMS PASS:
 *     readNetwork + compression_nonaffine_velocities (both unchanged from
 *     v3) followed by write_mode_terms_text(), then exit. The v3 elasticity
 *     / aux-file pass below is skipped on this invocation, so Python need
 *     only diagonalise H and apply the same 1/omega^2 weights C already
 *     applied -- no bond-loop tensor contraction is reimplemented in Python.
 *
 *   UNCHANGED FROM v3
 *     Everything else: xi_eta_x, hessianOnVector, u_eta_x_x_on_two_vectors,
 *     tessianOnThreeVectors, tessianOnTwoVectors, tessianOnTwoVectorsThenDot,
 *     cgSolver, calculateBulkModulus/calculateShearModulus,
 *     compression_nonaffine_velocities, readNetwork, and all v3 writers
 *     (write_bonds_text, write_meta_text, write_hessian_text,
 *      write_vna_eta_text).
 *
 * COMPILE:
 *   gcc -O2 -o dPdT_v4 3D_spring_network_dPdT_test_v4.c -lm
 *
 * USAGE:
 *   pass 1 (elasticity + aux files, identical to v3):
 *     ./dPdT_v4 network_N4000_z16.24_s0000.dat 0
 *   [notebook diagonalises hessian_g0.txt -> eigvecs_g0.bin, eigvals_g0.txt]
 *   pass 2 (mode terms for Eq.21/22, v4 ADDITION):
 *     ./dPdT_v4 network_N4000_z16.24_s0000.dat 0 eigvecs_g0.bin
 * ============================================================================
 */

#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* ── system size ─────────────────────────────────────────────────────────── */
#define N                       4000
#define DOUBLE_N                ((double)N)
#define DIM                     3
#define DOFS                    (DIM*N)   /* ADDED in v4 */
#define CGSOLVER_TOL            1e-10
#define HUGE_INT                10000000

/* ── index / component aliases ───────────────────────────────────────────── */
#define DELTA(a,b)              ((a)==(b)?1.0:0.0)
#define X_COMP                  0
#define Y_COMP                  1
#define Z_COMP                  2
#define I_INDEX                 0
#define J_INDEX                 1

/* ── PBC geometry ────────────────────────────────────────────────────────── */
#define LENGTH                  1.0
#define HALF_LENGTH             0.5

/* ── RNG ─────────────────────────────────────────────────────────────────── */
#define UNI                     ((double)rand()/((double)RAND_MAX + 1.0))

/* ── linearOperator function pointer (ADDED from spring_networks_elasticity.c) ─
 * Set to &hessianOnVector in main() so cgSolver can be reused with different
 * operators if needed.                                                        */
void (*linearOperator)(double *actOn, double *res);

/* ── globals ─────────────────────────────────────────────────────────────── */
double L, invL, V;
double pressure = 0.0;   /* ADDED: T=0 virial pressure (0 for stress-free)   */
int    contacts;

double *rx, *ry, *rz;
int    *lookupBond;
double *rij;     /* x_ij = x_j - x_i, real units, length DIM*contacts        */
double *inv_r;
double *first;   /* [1st]_l  -- identically 0 for this network               */
double *second;  /* [2nd]_l  = 1/r_ij^2 for this network                     */
double *third;   /* [3rd]_l  = -3*second[l]/r_ij^2   (CORRECTION applied)    */


/* ============================================================================
 *  cdot -- dot product of two 3N-vectors
 * ========================================================================== */
double cdot(double *vec1, double *vec2){
	double res = 0.0;
	int i;
	for (i=0;i<DIM*N;i++)
		res += vec1[i]*vec2[i];
	return res;
}

/* ============================================================================
 *  normalizeVector -- zero centre-of-mass drift, then unit-normalise
 * ========================================================================== */
void normalizeVector(double *v){
	int a,i,q;
	double D,moo = 0.0;
	for (q=0;q<DIM;q++){
		D = 0.0;
		for (i=0;i<N;i++)
			D += v[i*DIM + q];
		D = D/DOUBLE_N;
		for (i=0;i<N;i++)
			v[i*DIM + q] -= D;
	}
	for (a=0;a<DIM*N;a++)
		moo += v[a]*v[a];
	moo = 1.0/sqrt(moo);
	for (a=0;a<DIM*N;a++)
		v[a] = v[a]*moo;
	return;
}

/* ============================================================================
 *  setRandomDirection -- random unit vector in 3N-space
 * ========================================================================== */
void setRandomDirection(double *v){
	int a;
	for (a=0;a<DIM*N;a++)
		v[a] = 2.0*UNI - 1.0;
	normalizeVector(v);
	return;
}

/* ============================================================================
 *  readNetwork -- read spring-network file; populate rij, first, second, third
 *  SOURCE: 3D_spring_network_dPdT_test_v2.c  (third[] correction already there)
 * ========================================================================== */
int readNetwork(char *inFileName){
	int i,j,l,foo;
	double dummy;
	double dx,dy,dz,r2;
	FILE *inFile;

	rx = (double *)malloc(sizeof(double)*N);
	ry = (double *)malloc(sizeof(double)*N);
	rz = (double *)malloc(sizeof(double)*N);

	inFile = fopen(inFileName,"rb");
	if ( inFile ){
		foo = fscanf(inFile,"%lf %lf %lf",&L, &dummy, &dummy);
		invL = 1.0/L; V = L*L*L;
		for (i=0;i<N;i++)
			foo = fscanf(inFile,"%lf %lf %lf", &(rx[i]), &(ry[i]), &(rz[i]));

		foo = fscanf(inFile,"%d %lf %lf", &contacts, &dummy, &dummy);
		lookupBond = (int *)malloc(sizeof(int)*2*contacts);
		rij    = (double *)malloc(sizeof(double)*contacts*DIM);
		inv_r  = (double *)malloc(sizeof(double)*contacts);
		first  = (double *)malloc(sizeof(double)*contacts);
		second = (double *)malloc(sizeof(double)*contacts);
		third  = (double *)malloc(sizeof(double)*contacts);
		for (l=0;l<contacts;l++){
			foo = fscanf(inFile,"%d %d %lf", &(lookupBond[2*l + I_INDEX]), &(lookupBond[2*l + J_INDEX]), &dummy);
			i = lookupBond[2*l + I_INDEX];
			j = lookupBond[2*l + J_INDEX];
			dx = rx[j] - rx[i];
			dy = ry[j] - ry[i];
			dz = rz[j] - rz[i];
			if ( dz >= HALF_LENGTH )      dz -= LENGTH;
			else if ( dz < -HALF_LENGTH ) dz += LENGTH;
			if ( dy >= HALF_LENGTH )      dy -= LENGTH;
			else if ( dy < -HALF_LENGTH ) dy += LENGTH;
			if ( dx >= HALF_LENGTH )      dx -= LENGTH;
			else if ( dx < -HALF_LENGTH ) dx += LENGTH;
			r2 = L*L*( dx*dx + dy*dy + dz*dz );
			rij[DIM*l + X_COMP] = L*dx;
			rij[DIM*l + Y_COMP] = L*dy;
			rij[DIM*l + Z_COMP] = L*dz;
			
			/* General generalized spatial derivatives (U >= 0) */
			double r = sqrt(r2);
			double r0 = r; /* Network starts at U=0, so r0 equals the initial distance */
			
			inv_r[l]  = 1.0 / r;
			first[l]  = 1.0 - (r0 / r);
			second[l] = r0 / (r2 * r);
			third[l]  = -3.0 * r0 / (r2 * r2 * r);
		}
		fclose(inFile);
		return 1;
	}
	printf("couldn't read %s\n",inFileName);
	return 0;
}

/* ============================================================================
 *  free_everything
 * ========================================================================== */
void free_everything(){
	free(rx); free(ry); free(rz);
	free(lookupBond); free(second);
	free(first); free(third); free(rij); free(inv_r);
	return;
}

/* ============================================================================
 *  test_network  (ADDED from spring_networks_elasticity.c)
 *  Prints min/max per-node coordination as a sanity check after readNetwork().
 * ========================================================================== */
void test_network(){
	int i,j,l;
	int *num_of_nebz;
	int mn = 100000, mx = 0;

	num_of_nebz = (int *)calloc(N, sizeof(int));
	for (l=0;l<contacts;l++){
		num_of_nebz[lookupBond[2*l+I_INDEX]]++;
		num_of_nebz[lookupBond[2*l+J_INDEX]]++;
	}
	for (i=0;i<N;i++){
		if (num_of_nebz[i] < mn) mn = num_of_nebz[i];
		if (num_of_nebz[i] > mx) mx = num_of_nebz[i];
	}
	printf("Network: N=%d  contacts=%d  z_mean=%.4f  z_min=%d  z_max=%d\n",
	       N, contacts, 2.0*(double)contacts/DOUBLE_N, mn, mx);
	free(num_of_nebz);
}

/* ============================================================================
 *  Eq.(14):  xi_eta_x[k]  =  d^2 U / d eta d x_k
 *  SOURCE: 3D_spring_network_dPdT_test_v2.c  (unchanged)
 * ========================================================================== */
void xi_eta_x(double *res){
	int d,i,j,l;
	double r2,coeff;

	for (i=0;i<DIM*N;i++)
		res[i] = 0.0;

	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		for (r2=0.0,d=0;d<DIM;d++)
			r2 += rij[DIM*l+d]*rij[DIM*l+d];

		coeff = second[l]*r2 + 2.0*first[l];
		for (d=0;d<DIM;d++){
			res[DIM*j + d] += coeff*rij[DIM*l + d];
			res[DIM*i + d] -= coeff*rij[DIM*l + d];
		}
	}
	return;
}

/* ============================================================================
 *  u_epsilon  (ADDED from spring_networks_elasticity.c)
 *  ∂U/∂ε_αβ = Σ_l first[l] * rij_α * rij_β
 *  For stress-free networks (first[l]=0) this is always 0; kept for generality.
 * ========================================================================== */
double u_epsilon(int alpha, int beta){
	int l;
	double res = 0.0;
	for (l=0;l<contacts;l++)
		res += first[l]*rij[DIM*l + alpha]*rij[DIM*l + beta];
	return res;
}

/* ============================================================================
 *  u_epsilon_epsilon  (ADDED from spring_networks_elasticity.c)
 *  ∂²U/∂ε_αβ∂ε_κχ = Σ_l second[l] * rij_α * rij_β * rij_κ * rij_χ
 *  This is the Born (affine) contribution to the elastic tensor.
 * ========================================================================== */
double u_epsilon_epsilon(int alpha, int beta, int kappa, int chi){
	int l;
	double res = 0.0;
	for (l=0;l<contacts;l++)
		res += second[l]
		     * rij[DIM*l + alpha] * rij[DIM*l + beta]
		     * rij[DIM*l + kappa] * rij[DIM*l + chi];
	return res;
}

/* ============================================================================
 *  u_epsilon_x  (ADDED from spring_networks_elasticity.c)
 *  ∂²U/∂ε_αβ∂u_i^q  stored as a 3N-vector.
 *  Provides the non-affine driving force for shear / bulk CG solves.
 * ========================================================================== */
void u_epsilon_x(int alpha, int beta, double *res){
	int q,i,j,l;

	for (q=0;q<DIM;q++)
		for (i=0;i<N;i++)
			res[DIM*i + q] = 0.0;

	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		for (q=0;q<DIM;q++){
			double contrib = second[l]*rij[DIM*l+alpha]*rij[DIM*l+beta]*rij[DIM*l+q]
			               + first[l]*(DELTA(alpha,q)*rij[DIM*l+beta]
			                         + DELTA(beta, q)*rij[DIM*l+alpha]);
			res[DIM*i + q] -= contrib;
			res[DIM*j + q] += contrib;
		}
	}
	return;
}

/* ============================================================================
 *  u_eta_x_x_on_two_vectors -- term1 (harmonic), with [3rd] correction.
 *  SOURCE: 3D_spring_network_dPdT_test_v2.c  (unchanged)
 * ========================================================================== */

double u_eta_x_x_on_two_vectors(double *vec1, double *vec2){
	double v1_dot_r,v2_dot_r,v1_dot_v2,r2,res = 0.0;
	int l,i,j,d;
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		r2=0.0; v1_dot_r = 0.0; v2_dot_r = 0.0; v1_dot_v2 = 0.0;
		for (d=0;d<DIM;d++){
			r2 += rij[DIM*l + d]*rij[DIM*l + d];
			v1_dot_r += rij[DIM*l + d]*(vec1[DIM*j+d]-vec1[DIM*i+d]);
			v2_dot_r += rij[DIM*l + d]*(vec2[DIM*j+d]-vec2[DIM*i+d]);
			v1_dot_v2 += (vec1[DIM*j+d]-vec1[DIM*i+d])*(vec2[DIM*j+d]-vec2[DIM*i+d]);
		}
		res += third[l]*r2*v1_dot_r*v2_dot_r + second[l]*(4.0*v1_dot_r*v2_dot_r + r2*v1_dot_v2) + 2.0*first[l]*v1_dot_v2;
	}
	return res;
}

/* ============================================================================
 *  hessianOnVector -- H . vec  (Eq.11)
 *  SOURCE: 3D_spring_network_dPdT_test_v2.c  (unchanged)
 * ========================================================================== */
void hessianOnVector(double *vec, double *res){
	int q,i,j,l=0;
	double uij[DIM];
	double u_dot_r;

	for (i=0;i<N;i++)
		for (q=0;q<DIM;q++)
			res[DIM*i + q] = 0.0;

	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		for (q=0;q<DIM;q++)
			uij[q] = vec[DIM*j + q] - vec[DIM*i + q];

		for (u_dot_r=0.0,q=0;q<DIM;q++)
			u_dot_r += rij[DIM*l + q]*uij[q];

		for (q=0;q<DIM;q++){
			res[DIM*i + q] -= second[l]*u_dot_r*rij[DIM*l + q] + first[l]*uij[q];
			res[DIM*j + q] += second[l]*u_dot_r*rij[DIM*l + q] + first[l]*uij[q];
		}
	}
	return;
}

/* ============================================================================
 *  ONE-PASS:  tessianOnThreeVectors -- v . U''' : (u (x) w)   (scalar)
 *  SOURCE: 3D_spring_network_dPdT_test_v2.c  (unchanged)
 * ========================================================================== */
double tessianOnThreeVectors(double *v, double *u, double *w){
	int q,i,j,l=0;
	double vij[DIM],uij[DIM],wij[DIM];
	double v_dot_r,u_dot_r,w_dot_r,v_dot_u,w_dot_v,u_dot_w;
	double res = 0.0;

	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		for (q=0;q<DIM;q++){
			vij[q] = v[DIM*j + q] - v[DIM*i + q];
			uij[q] = u[DIM*j + q] - u[DIM*i + q];
			wij[q] = w[DIM*j + q] - w[DIM*i + q];
		}

		v_dot_r = 0.0; u_dot_r = 0.0; w_dot_r = 0.0;
		v_dot_u = 0.0; w_dot_v = 0.0; u_dot_w = 0.0;
		for (q=0;q<DIM;q++){
			u_dot_r += rij[DIM*l + q]*uij[q];
			v_dot_r += rij[DIM*l + q]*vij[q];
			w_dot_r += rij[DIM*l + q]*wij[q];
			v_dot_u += vij[q]*uij[q];
			u_dot_w += uij[q]*wij[q];
			w_dot_v += wij[q]*vij[q];
		}

		res += third[l]*u_dot_r*v_dot_r*w_dot_r + second[l]*(v_dot_u*w_dot_r + u_dot_w*v_dot_r + w_dot_v*u_dot_r);
	}
	return res;
}

/* ============================================================================
 *  TWO-PASS:  tessianOnTwoVectors -- res = U''' : (v (x) u)   (VECTOR)
 *  SOURCE: 3D_spring_network_dPdT_test_v2.c  (unchanged)
 * ========================================================================== */
void tessianOnTwoVectors(double *v, double *u, double *res){
	int q,i,j,l=0;
	double vij[DIM],uij[DIM];
	double v_dot_r, u_dot_r, v_dot_u;
	double coeff;

	for (q=0;q<DIM;q++)
		for (i=0;i<N;i++)
			res[DIM*i + q] = 0.0;

	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];

		for (q=0;q<DIM;q++){
			vij[q] = v[DIM*j + q] - v[DIM*i + q];
			uij[q] = u[DIM*j + q] - u[DIM*i + q];
		}

		v_dot_r = 0.0; u_dot_r = 0.0; v_dot_u = 0.0;
		for (q=0;q<DIM;q++){
			u_dot_r += rij[DIM*l + q]*uij[q];
			v_dot_r += rij[DIM*l + q]*vij[q];
			v_dot_u += vij[q]*uij[q];
		}

		for (q=0;q<DIM;q++){
			coeff = third[l]*u_dot_r*v_dot_r*rij[DIM*l + q]
			      + second[l]*( v_dot_u*rij[DIM*l + q] + v_dot_r*uij[q] + u_dot_r*vij[q] );
			res[DIM*i + q] -= coeff;
			res[DIM*j + q] += coeff;
		}
	}
	return;
}

/* ============================================================================
 *  tessianOnTwoVectorsThenDot -- explicit two-pass equivalent of
 *  tessianOnThreeVectors(v_dot_vec, u, w).
 *  SOURCE: 3D_spring_network_dPdT_test_v2.c  (unchanged)
 * ========================================================================== */
double tessianOnTwoVectorsThenDot(double *v_dot_vec, double *u, double *w, double *tmp){
	tessianOnTwoVectors(u, w, tmp);
	return cdot(v_dot_vec, tmp);
}

/* ============================================================================
 *  mode_terms -- ADDED in v4. Per-mode contributions to the Eq.(21)/(22)
 *  mode sum, for a single eigenvector Psi (||Psi||=1) of the Hessian H with
 *  eigenvalue omega^2:
 *
 *    term1 = u_eta_x_x_on_two_vectors(Psi,Psi)
 *          = Psi . (d^3 U / d eta dx dx) . Psi                     (Eq.20)
 *
 *    term2 = cdot( vna_eta, tessianOnTwoVectors(Psi,Psi) )
 *          = vna_eta . [ U''' : (Psi (x) Psi) ]
 *          = tessianOnThreeVectors(vna_eta,Psi,Psi)                (Eq.21)
 *
 *  term2 is the "two vectors, then cdot" form requested in place of the
 *  one-pass tessianOnThreeVectors: tessianOnTwoVectors(Psi,Psi,tmp) builds
 *  the vector U''':(Psi(x)Psi) from the bond loop ONCE per mode, and cdot()
 *  contracts it with vna_eta (itself fixed, = -H^-1.xi_eta_x, independent of
 *  the mode). No O(n_modes^2) pass over mode PAIRS is ever performed.
 *
 *  `tmp` must be a pre-allocated DOFS-length scratch buffer, reused by the
 *  caller across modes.
 * ========================================================================== */
void mode_terms(double *Psi, double *vna_eta, double *tmp,
                double *term1, double *term2){
	*term1 = u_eta_x_x_on_two_vectors(Psi, Psi);
	tessianOnTwoVectors(Psi, Psi, tmp);
	*term2 = cdot(vna_eta, tmp);
	return;
}

/* ============================================================================
 *  cgSolver -- solves H x = b for x (in place). Returns 1 on convergence.
 *  MODIFIED: now calls linearOperator() instead of hard-wired hessianOnVector()
 *  so the solver can be reused by calculateBulkModulus / calculateShearModulus.
 *  Stopping criterion: maxNorm (absolute) from spring_networks_elasticity.c,
 *  which is more robust than the scale-relative test in v2 for near-zero b.
 * ========================================================================== */
int cgSolver(double *x, double *b){
	double *r, *w, *z;
	double a,B,maxNorm,typicalNorm,sum,scale;
	int i,k,q;

	scale = sqrt(cdot(b,b)/( (double)DIM*DOUBLE_N ));

	r = (double *)malloc(sizeof(double)*DIM*N);
	w = (double *)malloc(sizeof(double)*DIM*N);
	z = (double *)malloc(sizeof(double)*DIM*N);

	linearOperator(x,r);

	for (i=0;i<DIM*N;i++){
		r[i] = b[i] - r[i];
		w[i] = -r[i];
	}

	linearOperator(w,z);
	a = cdot(r,w)/cdot(w,z);

	for (i=0;i<DIM*N;i++)
		x[i] += a*w[i];

	B = 0.0;

	for (k=0; k<HUGE_INT; k++){
		maxNorm = 0.0; typicalNorm = 0.0;
		for (i=0;i<DIM*N;i++){
			r[i] = r[i] - a*z[i];
			if (r[i]*r[i] > maxNorm)
				maxNorm = r[i]*r[i];
			typicalNorm += r[i]*r[i];
		}
		maxNorm     = sqrt(maxNorm);
		typicalNorm = sqrt(typicalNorm/((double)DIM*DOUBLE_N));

		if ( maxNorm < CGSOLVER_TOL )
			break;

		if ( k>0 && !(k%16384) )
			printf("  CG it%d  typNorm=%g  maxNorm=%g\n",k,typicalNorm,maxNorm);

		B = cdot(r,z)/cdot(w,z);
		for (i=0;i<DIM*N;i++)
			w[i] = -r[i] + B*w[i];

		linearOperator(w,z);
		a = cdot(r,w)/cdot(w,z);

		for (i=0;i<DIM*N;i++)
			x[i] += a*w[i];
	}

	for (q=0;q<DIM;q++){
		sum = 0.0;
		for (i=0;i<N;i++)
			sum += x[DIM*i+q];
		sum = sum/DOUBLE_N;
		for (i=0;i<N;i++)
			x[DIM*i+q] -= sum;
	}
	free(r); free(w); free(z);
	if (k==HUGE_INT)
		return 0;
	else
		return 1;
}

/* ============================================================================
 *  calculateBulkModulus  (ADDED from spring_networks_elasticity.c)
 *
 *  K_born = (1/V) Σ_l [(second[l]*r² + first[l])*r² + first[l]*r²]
 *         = contacts/V  for stress-free spring network
 *
 *  K_na   = (1/V) ξ · H⁻¹ · ξ    where ξ is the isotropic driving force
 *
 *  K      = (K_born − K_na) / (DIM*DIM)  +  pressure
 *           = (K_born − K_na) / 9          for P0 = 0
 *
 *  Returns 1 on success, 0 if CG fails.
 * ========================================================================== */
int calculateBulkModulus(double *k_born, double *k_na){
	int a,q,i,j,l;
	double r2;
	double *nav, *xi;

	nav = (double *)malloc(sizeof(double)*N*DIM);
	xi  = (double *)malloc(sizeof(double)*N*DIM);
	for (a=0;a<DIM*N;a++){ nav[a]=0.0; xi[a]=0.0; }
	*k_born = 0.0;

	for (l=0;l<contacts;l++){
		i  = lookupBond[2*l + I_INDEX];
		j  = lookupBond[2*l + J_INDEX];
		r2 = 0.0;
		for (q=0;q<DIM;q++) r2 += rij[DIM*l+q]*rij[DIM*l+q];

		*k_born += (second[l]*r2 + first[l])*r2 + first[l]*r2;

		for (q=0;q<DIM;q++){
			xi[DIM*j + q] += (second[l]*r2 + first[l])*rij[DIM*l+q];
			xi[DIM*i + q] -= (second[l]*r2 + first[l])*rij[DIM*l+q];
		}
	}
	*k_born /= V;

	if (cgSolver(nav,xi)){
		*k_na = cdot(nav,xi)/V;
		free(nav); free(xi);
		return 1;
	}
	printf("calculateBulkModulus: CG solver failed\n");
	free(nav); free(xi);
	return 0;
}

/* ============================================================================
 *  calculateShearModulus  (ADDED from spring_networks_elasticity.c)
 *
 *  μ_born = [u_εε(x,y,x,y) + u_ε(y,y)] / V
 *  μ_na   = (1/V) ξ · H⁻¹ · ξ   where ξ = u_ε_x(X_COMP, Y_COMP)
 *  μ      = μ_born − μ_na
 *
 *  Returns 1 on success, 0 if CG fails.
 * ========================================================================== */
int calculateShearModulus(double *mu_born, double *mu_na){
	int i,q;
	double *nav, *xi;

	nav = (double *)malloc(sizeof(double)*N*DIM);
	xi  = (double *)malloc(sizeof(double)*N*DIM);

	u_epsilon_x(X_COMP,Y_COMP,xi);
	for (i=0;i<N;i++)
		for (q=0;q<DIM;q++) nav[DIM*i+q]=0.0;

	if (cgSolver(nav,xi)){
		*mu_na   = cdot(nav,xi)/V;
		*mu_born = (u_epsilon_epsilon(X_COMP,Y_COMP,X_COMP,Y_COMP)
		           + u_epsilon(Y_COMP,Y_COMP)) / V;
		free(nav); free(xi);
		return 1;
	}
	printf("calculateShearModulus: CG solver failed\n");
	free(nav); free(xi);
	return 0;
}

/* ============================================================================
 *  compression_nonaffine_velocities -- vna_eta = -H^-1 . xi_eta_x
 *  SOURCE: 3D_spring_network_dPdT_test_v2.c  (unchanged)
 * ========================================================================== */
int compression_nonaffine_velocities(double *res){
	int a;
	double *xi = (double *)malloc(sizeof(double)*N*DIM);

	for (a=0;a<DIM*N;a++){
		res[a] = 0.0;
		xi[a]  = 0.0;
	}

	xi_eta_x(xi);

	if ( cgSolver(res,xi) ){
		for (a=0;a<DIM*N;a++)
			res[a] = -res[a];
		free(xi);
		return 1;
	}
	free(xi);
	return 0;
}

/* ============================================================================
 *  Output writers
 * ========================================================================== */

/* write_bonds_text -- MODIFIED: adds physical distance r column;
 * third computed as -3/r^4 matching spring_networks_elasticity.c.
 * Previously v2 omitted r and wrote a hardcoded 0 for third.             */
void write_bonds_text(int glass_idx){
	int l;
	char fname[256];
	sprintf(fname, "bonds_g%d.txt", glass_idx);
	FILE *f = fopen(fname,"w");
	fprintf(f, "# i j rij_x rij_y rij_z r first second third\n");
	for (l=0;l<contacts;l++)
		fprintf(f, "%d %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n",
			lookupBond[2*l + I_INDEX], lookupBond[2*l + J_INDEX],
			rij[DIM*l + X_COMP], rij[DIM*l + Y_COMP], rij[DIM*l + Z_COMP],
			1.0/inv_r[l],          /* physical bond length r                */
			first[l],
			second[l],
			third[l]);             /* -3/r^4 (set in readNetwork)           */
	fclose(f);
	printf("  wrote %s (%d bonds)\n", fname, contacts);
}

/* write_meta_text -- unchanged from v2 */
void write_meta_text(int glass_idx){
	char fname[256];
	sprintf(fname, "meta_g%d.txt", glass_idx);
	FILE *f = fopen(fname,"w");
	fprintf(f, "N %d\n", N);
	fprintf(f, "DIM %d\n", DIM);
	fprintf(f, "contacts %d\n", contacts);
	fprintf(f, "V %.15g\n", V);
	fprintf(f, "L %.15g\n", L);
	fclose(f);
}

/* write_hessian_text -- unchanged from v2 */
void write_hessian_text(int glass_idx){
	int i,j,l,m,alpha,beta,a,b;
	double thisHes;
	double *diagonal;
	char fname[256];

	sprintf(fname, "hessian_g%d.txt", glass_idx);
	FILE *f = fopen(fname,"w");

	diagonal = (double *)calloc(N*DIM*DIM, sizeof(double));

	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		m = 0;
		for (alpha=0;alpha<DIM;alpha++){
			for (beta=0;beta<DIM;beta++){
				thisHes = -second[l]*rij[DIM*l+alpha]*rij[DIM*l+beta]
				          -first[l]*DELTA(alpha,beta);
				a = DIM*i + alpha;
				b = DIM*j + beta;
				fprintf(f, "%d\t%d\t%.15g\n%d\t%d\t%.15g\n", a,b,thisHes, b,a,thisHes);
				diagonal[DIM*DIM*i + m] -= thisHes;
				diagonal[DIM*DIM*j + m] -= thisHes;
				m++;
			}
		}
	}
	for (i=0;i<N;i++){
		m = 0;
		for (alpha=0;alpha<DIM;alpha++){
			for (beta=0;beta<DIM;beta++){
				a = DIM*i + alpha;
				b = DIM*i + beta;
				fprintf(f, "%d\t%d\t%.15g\n", a,b, diagonal[DIM*DIM*i + m]);
				m++;
			}
		}
	}
	free(diagonal);
	fclose(f);
	printf("  wrote %s (%d×%d sparse matrix)\n", fname, DIM*N, DIM*N);
}

/* write_vna_eta_text -- unchanged from v2 */
void write_vna_eta_text(int glass_idx){
	int i;
	double *vna;
	char fname[256];

	vna = (double *)calloc(DIM*N, sizeof(double));
	compression_nonaffine_velocities(vna);

	sprintf(fname, "vna_eta_g%d.txt", glass_idx);
	FILE *f = fopen(fname,"w");
	fprintf(f, "# vna_eta: compression non-affine velocity (N rows x 3 cols)\n");
	for (i=0;i<N;i++)
		fprintf(f, "%.15g\t%.15g\t%.15g\n", vna[DIM*i+X_COMP], vna[DIM*i+Y_COMP], vna[DIM*i+Z_COMP]);
	fclose(f);
	printf("  wrote %s\n", fname);
	free(vna);
}

/* ============================================================================
 *  write_mode_terms_text -- ADDED in v4.
 *
 *  INPUTS  (written by the notebook AFTER diagonalising hessian_g{idx}.txt):
 *    eigvecsFileName = "eigvecs_g{idx}.bin"
 *        raw float64, n_modes back-to-back blocks of DOFS values each
 *        (block m = mode m's eigenvector Psi_m, ||Psi_m||=1), mode order
 *        = ascending omega^2 with near-zero (translation) modes removed.
 *    "eigvals_g{idx}.txt"
 *        n_modes lines, omega^2 per mode, SAME order as eigvecs.
 *
 *  For every mode: mode_terms() -> (term1,term2)  [Eq.20/21, "two vectors
 *  then cdot"], written unweighted to mode_terms_g{idx}.txt, AND accumulated
 *  into the Eq.(21)/(22) weighted sums
 *      harm = -1/(V*DIM) * sum_m term1_m/omega2_m
 *      anh  = -1/(V*DIM) * sum_m term2_m/omega2_m
 *  written to dPdT_results_g{idx}.txt as
 *      harm  anh  (harm+anh)=Eq.21-LHS  0.5*(harm+anh)=Eq.22-LHS  V
 *  so Python's remaining job is just to read this file (or, for diagnostics,
 *  mode_terms_g{idx}.txt + its own omega^2 for a cross-check sum) -- the
 *  bond-loop tensor contractions are never reimplemented in Python.
 * ========================================================================== */
void write_mode_terms_text(int glass_idx, char *eigvecsFileName, double *vna_eta){
	FILE *fin, *feig, *fmodes, *fres;
	char fname[256];
	long fsize, n_modes, n_eigvals, m;
	double *Psi, *tmp, *omega2;
	double term1, term2, harm_sum=0.0, anh_sum=0.0, pre, harm, anh, S;

	fin = fopen(eigvecsFileName, "rb");
	if (!fin){
		printf("write_mode_terms_text: couldn't open %s\n", eigvecsFileName);
		return;
	}
	fseek(fin, 0, SEEK_END);
	fsize = ftell(fin);
	rewind(fin);
	n_modes = fsize / ((long)DOFS * (long)sizeof(double));
	printf("  %s : %ld modes x DOFS=%d (%.3f MB)\n",
	       eigvecsFileName, n_modes, DOFS, fsize/1.0e6);

	/* ── load omega^2 per mode, same order as eigvecs ────────────────────── */
	sprintf(fname, "eigvals_g%d.txt", glass_idx);
	feig = fopen(fname,"r");
	if (!feig){
		printf("write_mode_terms_text: couldn't open %s\n", fname);
		fclose(fin);
		return;
	}
	omega2 = (double *)malloc(sizeof(double)*n_modes);
	for (n_eigvals=0; n_eigvals<n_modes; n_eigvals++){
		if (fscanf(feig, "%lf", &(omega2[n_eigvals])) != 1)
			break;
	}
	fclose(feig);
	if (n_eigvals != n_modes){
		printf("write_mode_terms_text: %s has %ld values, expected %ld -- aborting\n",
		       fname, n_eigvals, n_modes);
		free(omega2); fclose(fin);
		return;
	}

	Psi = (double *)malloc(sizeof(double)*DOFS);
	tmp = (double *)malloc(sizeof(double)*DOFS);

	sprintf(fname, "mode_terms_g%d.txt", glass_idx);
	fmodes = fopen(fname, "w");
	fprintf(fmodes, "# term1 term2  (Eq.20 harmonic / Eq.21 anharmonic, "
	                "unweighted; pair with eigvals_g%d.txt for 1/omega^2)\n",
	                glass_idx);

	for (m=0; m<n_modes; m++){
		if (fread(Psi, sizeof(double), DOFS, fin) != (size_t)DOFS){
			printf("  write_mode_terms_text: short read at mode %ld -- stopping\n", m);
			break;
		}
		mode_terms(Psi, vna_eta, tmp, &term1, &term2);
		fprintf(fmodes, "%.15g\t%.15g\n", term1, term2);

		harm_sum += term1/omega2[m];
		anh_sum  += term2/omega2[m];

		if (m>0 && (m % 2000)==0)
			printf("    mode %ld/%ld\n", m, n_modes);
	}
	fclose(fin);
	fclose(fmodes);
	printf("  wrote %s (%ld modes)\n", fname, n_modes);

	pre  = -1.0/(V*(double)DIM);
	harm = pre*harm_sum;
	anh  = pre*anh_sum;
	S    = harm + anh;

	sprintf(fname, "dPdT_results_g%d.txt", glass_idx);
	fres = fopen(fname, "w");
	fprintf(fres, "# harmonic  anharmonic  2*(dp/dT)_V[Eq.21]  (dp/dT)_V[Eq.22]  V\n");
	fprintf(fres, "%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n", harm, anh, S, 0.5*S, V);
	fclose(fres);

	printf("  harmonic   contribution = %+.6e\n", harm);
	printf("  anharmonic contribution = %+.6e\n", anh);
	printf("  2*(dp/dT)_V (Eq.21 LHS) = %+.6e\n", S);
	printf("  (dp/dT)_V   (Eq.22 LHS) = %+.6e\n", 0.5*S);
	printf("  wrote %s\n", fname);

	free(Psi); free(tmp); free(omega2);
	return;
}

/* ============================================================================
 *  Optional: equivalence check between the one-pass and two-pass term2.
 *  SOURCE: 3D_spring_network_dPdT_test_v2.c  (unchanged)
 * ========================================================================== */
void check_one_pass_vs_two_pass_term2(void){
	double *vna_eta, *Psi, *tmp;
	double one_pass, two_pass;

	vna_eta = (double *)malloc(sizeof(double)*DIM*N);
	Psi     = (double *)malloc(sizeof(double)*DIM*N);
	tmp     = (double *)malloc(sizeof(double)*DIM*N);

	compression_nonaffine_velocities(vna_eta);
	setRandomDirection(Psi);

	one_pass = tessianOnThreeVectors(vna_eta, Psi, Psi);
	two_pass = tessianOnTwoVectorsThenDot(vna_eta, Psi, Psi, tmp);

	printf("\n--- one-pass vs two-pass term2 check (Psi = random direction) ---\n");
	printf("  tessianOnThreeVectors(vna_eta,Psi,Psi)      = %.15g\n", one_pass);
	printf("  cdot(vna_eta, tessianOnTwoVectors(Psi,Psi)) = %.15g\n", two_pass);
	printf("  difference                                   = %.3e\n", one_pass-two_pass);

	free(vna_eta); free(Psi); free(tmp);
}

/* ============================================================================
 *  main
 * ========================================================================== */
int main(int argc, char **argv){
	char networkFileName[1024];
	char outFileName[256];
	char eigvecsFileName[1024];
	int glass_idx = 0;
	int t0 = (int)time(0);

	double k_born, k_na, K;
	double mu_born, mu_na, mu;
	double *vna_eta;
	FILE *outFile;

	srand(12345);

	if (argc > 1)
		strcpy(networkFileName, argv[1]);
	else
		sprintf(networkFileName, "network_N%d_z16.24_s0000.dat", N);

	/* ── parse optional serial index from second argument ───────────────────── */
	if (argc > 2)
		sscanf(argv[2], "%d", &glass_idx);

	/* ── set linearOperator before first cgSolver call ───────────────────── */
	linearOperator = &hessianOnVector;

	if ( !readNetwork(networkFileName) ){
		printf("%s not found\n", networkFileName);
		return 1;
	}

	/* ── MODE-TERMS PASS (v4 addition) ─────────────────────────────────────────
	 * A 3rd CLI argument (path to eigvecs_g{idx}.bin, written by the
	 * notebook after scipy.linalg.eigh on hessian_g{idx}.txt -- see
	 * write_mode_terms_text() above) selects this pass: recompute
	 * vna_eta = -H^-1.xi_eta_x (unchanged compression_nonaffine_velocities),
	 * then call write_mode_terms_text() to get term1/term2 per mode via
	 * mode_terms() ["tessianOnTwoVectors, then cdot"] and the full Eq.21/22
	 * weighted sums, then exit. Does NOT repeat the v3 elasticity / aux-file
	 * pass below.                                                          */
	if (argc > 3){
		strcpy(eigvecsFileName, argv[3]);
		printf("Mode-terms pass: %s, glass %d\n", eigvecsFileName, glass_idx);

		vna_eta = (double *)malloc(sizeof(double)*DOFS);
		printf("  computing vna_eta = -H^-1.xi_eta_x ...\n");
		if (!compression_nonaffine_velocities(vna_eta)){
			printf("ERROR: compression_nonaffine_velocities (CG) failed\n");
			free(vna_eta);
			free_everything();
			return 1;
		}

		write_mode_terms_text(glass_idx, eigvecsFileName, vna_eta);

		free(vna_eta);
		free_everything();
		printf("Mode-terms pass done. Elapsed: %d s\n", (int)time(0) - t0);
		return 0;
	}

	/* ── v3 pass (unchanged below): elasticity + aux files for notebook ── */
	printf("read %s : N=%d, contacts=%d, z=%.4f, L=%.6f, V=%.6f\n",
		networkFileName, N, contacts, 2.0*(double)contacts/DOUBLE_N, L, V);

	test_network();

	printf("T=0 pressure P0 = %.10e  (expected ~0 for stress-free net)\n", pressure);

	/* ── bulk modulus ─────────────────────────────────────────────────────── */
	printf("Computing bulk modulus K...\n");
	if (!calculateBulkModulus(&k_born, &k_na)){
		printf("ERROR: calculateBulkModulus failed\n");
		free_everything();
		return 1;
	}
	K = (k_born - k_na) / ((double)(DIM*DIM)) + pressure;
	printf("  K_born = %.8f\n  K_na   = %.8f\n  K      = %.8f\n",
	       k_born, k_na, K);

	/* ── shear modulus ─────────────────────────────────────────────────── */
	printf("Computing shear modulus mu...\n");
	if (!calculateShearModulus(&mu_born, &mu_na)){
		printf("ERROR: calculateShearModulus failed\n");
		free_everything();
		return 1;
	}
	mu = mu_born - mu_na;
	printf("  mu_born = %.8f\n  mu_na   = %.8f\n  mu      = %.8f\n",
	       mu_born, mu_na, mu);

	/* ── scalar output ──────────────────────────────────────────────────────── */
	sprintf(outFileName, "snp_elasticity_N%d_z16.24_s%.3d.dat", N, glass_idx);
	outFile = fopen(outFileName,"a");
	if (glass_idx == 0)
		fprintf(outFile, "# serial  P0  K_born  K_na  K  mu_born  mu_na  mu\n");
	fprintf(outFile, "%d\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
	        glass_idx, pressure, k_born, k_na, K, mu_born, mu_na, mu);
	fclose(outFile);
	printf("Scalar results written to %s\n", outFileName);

	/* ── optional tessian check ────────────────────────────────────── */
	/* check_one_pass_vs_two_pass_term2(); */

	/* ── auxiliary files for Python notebook ────────────────────────── */
	printf("Writing auxiliary files...\n");
	write_bonds_text(glass_idx);
	write_meta_text(glass_idx);
	write_hessian_text(glass_idx);
	write_vna_eta_text(glass_idx);

	printf("Done. Elapsed: %d s\n", (int)time(0) - t0);
	free_everything();
	return 0;
}
