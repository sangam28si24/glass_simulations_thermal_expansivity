#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define	N							64000
#define	DOUBLE_N						((double)N)
#define	DIM							3
#define	CGSOLVER_TOL				1e-10

/* stuff that we don't touch usually */
#define	DELTA(a,b)					((a)==(b)?1.0:0.0)
#define	X_COMP					0
#define 	Y_COMP 					1
#define 	Z_COMP 					2
#define	I_INDEX					0
#define	J_INDEX					1


/* constant constants... :P*******************************************************/
#define	LENGTH				1.0
#define	HALF_LENGTH			0.5
#define	TWO_PI					6.28318530717958
#define	UNI		((double)rand()/((double)RAND_MAX + 1.0))
#define	HUGE_INT				10000000


void initializeSystem();
void updateNebzLists();
void calculateEverything();
void (*linearOperator)(double *actOn, double *res);	/* make this global, set it and then everyone knows it	*/

int serial;

/*globals*/
double L; 			/*	length of <square> box 	*/
double invL;			/* inverse length of box	*/
double V;			/*	volume			*/
double u;			/*	potential energy 		*/
double pressure = 0.0;

int contacts;

double *rx, *ry, *rz;

/************************ new bondnetwork nebz ******************************************************/
int *lookupBond;
double *second;
double *first;
double *rij;
double *inv_r;
/*******************************************************************************************************/

double getParticipationRatio(double *z){
	int i,q;
	double sum,temp;
	
	sum = 0.0;
	for (i=0;i<N;i++){
		temp = 0.0;
		for (q=0;q<DIM;q++)
			temp += z[DIM*i+q]*z[DIM*i+q];
		sum += temp*temp;
	}
	return 1.0/( DOUBLE_N * sum );
}

//also removes translations
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

void setRandomDirection(double *v){
	int a;
	for (a=0;a<DIM*N;a++)
		v[a] = 2.0*UNI - 1.0;
	normalizeVector(v);
	return;	
}

void writeVector(double *vec, char *outFileName){
	FILE *outFile;
	int i;
	outFile = fopen(outFileName,"wb");
	for (i=0;i<N;i++)
		fprintf(outFile,"%.8g\t%.8g\t%.8g\t%.15g\t%.15g\t%.15g\n",rx[i],ry[i],rz[i],
												vec[DIM*i+X_COMP],
												vec[DIM*i+Y_COMP],
												vec[DIM*i+Z_COMP] );
	fclose(outFile);
	return;
}

int readNetwork(char *inFileName){
	int i,j,l,foo;
	double dummy;
	double dx,dy,dz,r2;
	FILE *inFile;
	
	//allocate coordinates and node-forces
	rx = (double *)malloc(sizeof(double)*N);
	ry = (double *)malloc(sizeof(double)*N);
	rz = (double *)malloc(sizeof(double)*N);
	
	inFile = fopen(inFileName,"rb");
	if ( inFile ){
		foo = fscanf(inFile,"%lf %lf %lf",&L, &dummy, &dummy);
		invL = 1.0/L; V = L*L*L;
		for (i=0;i<N;i++)
			foo = fscanf(inFile,"%lf %lf %lf", &(rx[i]), &(ry[i]), &(rz[i])); //read nodes positions
				
		foo = fscanf(inFile,"%d %lf %lf", &contacts, &dummy, &dummy); 	// read number of bonds in the network
		lookupBond = (int *)malloc(sizeof(int)*2*contacts);				// allocate lookup-bond table memory
		rij = (double *)malloc(sizeof(double)*contacts*DIM);
		inv_r = (double *)malloc(sizeof(double)*contacts);
		first = (double *)malloc(sizeof(double)*contacts);
		second = (double *)malloc(sizeof(double)*contacts);
		for (l=0;l<contacts;l++){
			foo = fscanf(inFile,"%d %d %lf", &(lookupBond[2*l + I_INDEX]), &(lookupBond[2*l + J_INDEX]), &dummy ); //read bonds
			// I assume here that i < j
			i = lookupBond[2*l + I_INDEX];
			j = lookupBond[2*l + J_INDEX];
			dx = rx[j] - rx[i];
			dy = ry[j] - ry[i];
			dz = rz[j] - rz[i];
			// mess due to periodic boundary conditions 
			if ( dz >= HALF_LENGTH )
				dz -= LENGTH;
			else if ( dz < -HALF_LENGTH )
				dz += LENGTH;
			if ( dy >= HALF_LENGTH )
				dy -= LENGTH;
			else if ( dy < -HALF_LENGTH )
				dy += LENGTH;
			if ( dx >= HALF_LENGTH )
				dx -= LENGTH;
			else if ( dx < -HALF_LENGTH )
				dx += LENGTH;
			// end of mess 
			r2 = L*L*( dx*dx + dy*dy + dz*dz );
			rij[DIM*l + X_COMP] = L*dx;
			rij[DIM*l + Y_COMP] = L*dy;
			rij[DIM*l + Z_COMP] = L*dz;
			inv_r[l] = 1.0/sqrt(r2);
			first[l] = 0.0;
			second[l] = 1.0/r2;
		}
		fclose(inFile);
		return 1;
	}
	printf("couldn't read %s\n",inFileName);
	return 0;
}

//this is only useful is a snapshot has already been loaded...
int readDirection(double *z, char *fileName){
	int i,goo;
	double dummy;
	FILE *inFile;
	inFile = fopen(fileName, "rb");
	if ( inFile ){
		for (i=0;i<N;i++)
			goo = fscanf(inFile, "%lf %lf %lf %lf %lf %lf",&dummy,&dummy,&dummy,&(z[DIM*i + X_COMP]),&(z[DIM*i + Y_COMP]),&(z[DIM*i + Z_COMP]));
		fclose(inFile);
		return 1;
	}
	return 0;
}


int readMode(double *z, char *fileName){
	int a,goo;
	FILE *inFile;
	inFile = fopen(fileName, "rb");
	if ( inFile ){
		for (a=0;a<DIM*N;a++)
			goo = fscanf(inFile, "%lf",&(z[a]));
		fclose(inFile);
		return 1;
	}
	return 0;
}



/***********************************************************************************************************************/

double cdot( double *vec1, double *vec2 ){
	double res = 0.0;
	int i;
	for (i=0; i<DIM*N; i++)
		res += vec1[i]*vec2[i];
	return res;
}

double u_epsilon(int alpha, int beta){
	int i,k,m,l;
	double res = 0.0;
	
	for (l=0;l<contacts;l++)
		res += first[l]*rij[DIM*l + alpha]*rij[DIM*l + beta];
	
	return res;
}

double u_epsilon_epsilon(int alpha, int beta, int kappa, int chi){
	int i,j,k,m,l;
	double res = 0.0;
	
	for (l=0;l<contacts;l++)
		res += second[l]*rij[DIM*l + alpha]*rij[DIM*l + beta]*rij[DIM*l + kappa]*rij[DIM*l + chi];
	
	return res;
}

void u_epsilon_x(int alpha, int beta, double *res){
	int q,i,j,k,m,l;
	
	for (q=0;q<DIM;q++)
		for (i=0;i<N;i++)
			res[DIM*i + q] = 0.0;
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		for (q=0;q<DIM;q++){
			res[DIM*i + q] -= second[l]*rij[DIM*l + alpha]*rij[DIM*l + beta]*rij[DIM*l + q] + first[l]*( DELTA(alpha,q)*rij[DIM*l + beta] + DELTA(beta,q)*rij[DIM*l + alpha] );
			res[DIM*j + q] += second[l]*rij[DIM*l + alpha]*rij[DIM*l + beta]*rij[DIM*l + q] + first[l]*( DELTA(alpha,q)*rij[DIM*l + beta] + DELTA(beta,q)*rij[DIM*l + alpha] );
		}
	}
	return;
}

void hessianOnVector(double *vec, double *res){
	int q,i,j,k,m,l=0;
	double uij[DIM];
	double u_dot_r;
	
	for (q=0;q<DIM;q++)
		for (i=0;i<N;i++)
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

double hessianOnTwoVectors(double *u, double *v){
	int l,i,j,q;
	double u_dot_r,v_dot_r,u_dot_v,res=0.0;
	double uij[DIM];
	double vij[DIM];
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		for (q=0;q<DIM;q++){
			uij[q] = u[DIM*j + q] - u[DIM*i + q];
			vij[q] = v[DIM*j + q] - v[DIM*i + q];
		}
		
		for (v_dot_r=0.0,u_dot_r=0.0,u_dot_v=0.0,q=0;q<DIM;q++){
			u_dot_r += rij[DIM*l + q]*uij[q];
			v_dot_r += rij[DIM*l + q]*vij[q];
			u_dot_v += uij[q]*vij[q];
		}
		
		res += second[l]*u_dot_r*v_dot_r + first[l]*u_dot_v;
	}
	return res;
}



/***********************************************************************************************************************/

int cgSolver(double *x, double *b){
	double *r, *w, *z;
	double a,B,sum,scale;
	double maxNorm,typicalNorm;
	int i,j,k,t_0,q;
	
	t_0 = (int)time(0);
	
	scale = sqrt(cdot(b,b)/( (double)DIM*DOUBLE_N ) );
	
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
		maxNorm = 0.0;
		typicalNorm = 0.0;
		for (i=0;i<DIM*N;i++){
			r[i] = r[i] - a*z[i];
			if (r[i]*r[i] > maxNorm)
				maxNorm = r[i]*r[i];
			typicalNorm += r[i]*r[i];
		}
		maxNorm = sqrt(maxNorm);
		typicalNorm = sqrt(typicalNorm/((double)DIM*DOUBLE_N));
		
		//here the typical scale is always 1 = spring stiffness. 
		if ( maxNorm < CGSOLVER_TOL ) // || ( (typicalNorm < CGSOLVER_TOL) && (maxNorm < 10.0*CGSOLVER_TOL) )  )
			break;
		
		if ( k>0 && !(k%16384) )
			printf("it%d, typicalNorm %g, maxNorm %g\n",k,typicalNorm,maxNorm);
		
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

void writeHessian(char *outFileName){
	int i,m,l,j,alpha,beta,a,b;
	double *diagonal;
	FILE *outFile = fopen(outFileName,"wb");
	double thisHes;

	diagonal = (double *)malloc(sizeof(double)*N*DIM*DIM);

	for (i = 0; i < N; i++) {
		for (m = 0; m < DIM*DIM; m++) {
			diagonal[DIM * DIM * i + m] = 0.0;
		}
	}

	for (l=0; l<contacts; l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		m = 0;

		for (alpha=0;alpha<DIM;alpha++){
			for (beta=0;beta<DIM;beta++){
				thisHes = -second[l]*rij[DIM*l + alpha]*rij[DIM*l + beta] -first[l]*DELTA(alpha,beta);
				a = DIM*i + alpha;  // Flat index (i, \alpha)
				b = DIM*j + beta;   // Flat index (j, \beta)
				fprintf(outFile,"%d\t%d\t%.15g\n%d\t%d\t%.15g\n",a,b,thisHes,b,a,thisHes);
				diagonal[DIM*DIM*i + m] += second[l]*rij[DIM*l + alpha]*rij[DIM*l + beta] + first[l]*DELTA(alpha,beta);
				diagonal[DIM*DIM*j + m] += second[l]*rij[DIM*l + alpha]*rij[DIM*l + beta] + first[l]*DELTA(alpha,beta);
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
				fprintf(outFile,"%d\t%d\t%.15g\n", a, b, diagonal[DIM*DIM*i + m]);
				m++;
			}
		}
	}
	fclose(outFile);
	free(diagonal);
	return;
}

void free_everything(){
	free(rx); free(ry); free(rz); 
	free(lookupBond); free(second);
	free(first); free(rij); free(inv_r);
	return;
}



int calculateShearModulus(double *mu_born, double *mu_na){
        int i,l,q;
        double *nav, *xi;

        nav = (double *)malloc(sizeof(double)*N*DIM);
        xi = (double *)malloc(sizeof(double)*N*DIM);

        u_epsilon_x(X_COMP, Y_COMP, xi);
	
	for (i=0;i<N;i++)
                for (q=0;q<DIM;q++)
                        nav[DIM*i+q] = 0.0;

        if ( cgSolver(nav,xi) ){
                *mu_na = cdot(nav,xi)/V;
                *mu_born = ( u_epsilon_epsilon(X_COMP,Y_COMP,X_COMP,Y_COMP) + u_epsilon(Y_COMP,Y_COMP) )/V;
                free(nav); free(xi);
                return 1;
        }
        printf("failed to calculate mu\n");
	free(nav); free(xi);
        return 0;
}

int calculateBulkModulus(double *k_born, double *k_na){
        int a,q,i,j,l=0;
        double r2;
        double *nav, *xi;

        nav = (double *)malloc(sizeof(double)*N*DIM);
        xi = (double *)malloc(sizeof(double)*N*DIM);

        for (a=0;a<DIM*N;a++){
                nav[a] = 0.0;
                xi[a] = 0.0;
        }

        *k_born = 0.0;

        for (l=0;l<contacts;l++){
                i = lookupBond[2*l + I_INDEX];
                j = lookupBond[2*l + J_INDEX];

                r2 = 0.0;
                for (q=0;q<DIM;q++)
                        r2 += rij[DIM*l+q]*rij[DIM*l+q];

                *k_born += (second[l]*r2 + first[l])*r2 + first[l]*r2;

                for (q=0;q<DIM;q++){
                        xi[DIM*j + q] += (second[l]*r2 + first[l])*rij[DIM*l + q];
                        xi[DIM*i + q] -= (second[l]*r2 + first[l])*rij[DIM*l + q];
                }
	}

        *k_born = *k_born/V;

        if ( cgSolver(nav,xi) ){
                *k_na = cdot(nav,xi)/V;
                free(nav); free(xi);
                return 1;
        }
        free(nav); free(xi);
        return 0;
}


void test_network(){
	int i,j,l;
	int *num_of_nebz;
	int min_connectivity,max_connectivity;
	
	num_of_nebz = (int *)malloc(sizeof(int)*N);
	for (i=0;i<N;i++)
		num_of_nebz[i] = 0;
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		num_of_nebz[i]++;
		num_of_nebz[j]++;
	}
	
	min_connectivity = 100;
	max_connectivity = 0;
	for (i=0;i<N;i++){
		if ( num_of_nebz[i] > max_connectivity )
			max_connectivity = num_of_nebz[i];
		if ( num_of_nebz[i] < min_connectivity )
			min_connectivity = num_of_nebz[i];
	}
	printf("min = %d\n",min_connectivity);
	printf("max = %d\n",max_connectivity);
	
	free(num_of_nebz);
	return;
}

double get_median(double *numbers, int how_many){
	int q=0;
	int i;
	int largest_index;
	double largest;
	
	//~ for (q=0;q<how_many;q++)
		//~ printf("%g\n",numbers[q]);
	//~ printf("*********************\n");
	//~ q=0;
	while (q<how_many/2){
		//get next largest
		largest = numbers[0];
		largest_index = 0;
		for (i=1;i<how_many;i++){
			if ( numbers[i] > largest ){
				largest = numbers[i];
				largest_index = i;
			}
		}
		numbers[largest_index] = 0.0;
		q++;
	}
	//~ for (q=0;q<how_many;q++)
		//~ printf("%g\n",numbers[q]);
	//~ printf("\n\nmedian:%g\nhit enter\n",largest);
	//~ getchar();
	return largest;
}

double get_mean(double *numbers, int how_many){
	int q=0;
	double sum=0.0;
	
	for (q=0;q<how_many;q++)
		sum += numbers[q];
	
	return sum/(double)how_many;
}


//result goes into x, only work in rescaled coordinates.
void get_middle(int i, int j, double *x){
	double dx,dy,dz;
	double temp;
	
	dx = rx[j] - rx[i];
	dy = ry[j] - ry[i];
	dz = rz[j] - rz[i];
	
	// mess due to periodic boundary conditions 
	if ( dz >= HALF_LENGTH )
		dz -= LENGTH;
	else if ( dz < -HALF_LENGTH )
		dz += LENGTH;
	if ( dx >= HALF_LENGTH )
		dx -= LENGTH;
	else if ( dx < -HALF_LENGTH )
		dx += LENGTH;
	if ( dy >= HALF_LENGTH )
		dy -= LENGTH;
	else if ( dy < -HALF_LENGTH )
		dy += LENGTH;
	// end of mess
	//r = sqrt(dx*dx + dy*dy + dz*dz);
	
	temp = rx[i] + 0.5*dx;
	if ( temp > LENGTH )
		x[X_COMP] = temp - LENGTH;
	else if ( temp < 0.0 )
		x[X_COMP] = temp + LENGTH;
	else
		x[X_COMP] = temp;
	
	temp = ry[i] + 0.5*dy;
	if ( temp > LENGTH )
		x[Y_COMP] = temp - LENGTH;
	else if ( temp < 0.0 )
		x[Y_COMP] = temp + LENGTH;
	else
		x[Y_COMP] = temp;
	
	temp = rz[i] + 0.5*dz;
	if ( temp > LENGTH )
		x[Z_COMP] = temp - LENGTH;
	else if ( temp < 0.0 )
		x[Z_COMP] = temp + LENGTH;
	else
		x[Z_COMP] = temp;
	
	//printf("middle of %g,%g is %g\n",rx[i],rx[j],x[X_COMP]);
	//getchar();
	
	return;
}

void get_decay(int contact_index, double *z, double *c){
	int q,bin_index,d,w;
	int l; //iterate over contacts
	int i,j; //the pushed dipole
	double mid_pushed_r[DIM];
	double mid_projected_r[DIM];
	double r,rmn;
	int m,n; //the projected pair
	double dx,dy,dz;
	double binSize = 1.0;
	int numOfBins = (int)(L/binSize) + 1;
	double **projections;
	int *counters;
	double projection,zij;
	//FILE *dumpfile = fopen("dumpfile.dat","wb");
	
	i = lookupBond[2*contact_index + I_INDEX];
        j = lookupBond[2*contact_index + J_INDEX];
	
	get_middle(i,j,mid_pushed_r);
		
	projections = (double **)malloc(sizeof(double *)*numOfBins);
	counters = (int *)malloc(sizeof(int)*numOfBins);
	for (q=0;q<numOfBins;q++){
		projections[q] = (double *)malloc(sizeof(double)*contacts); //this is perfectly safe, but could be wasteful for large systems...
		counters[q] = 0;
	}
	
	for (l=0;l<contacts;l++){
		if ( l != contact_index ){
			m = lookupBond[2*l + I_INDEX];
			n = lookupBond[2*l + J_INDEX];
			get_middle(m,n,mid_projected_r);
			
			dx = mid_projected_r[X_COMP] - mid_pushed_r[X_COMP];
			dy = mid_projected_r[Y_COMP] - mid_pushed_r[Y_COMP];
			dz = mid_projected_r[Z_COMP] - mid_pushed_r[Z_COMP];
			// mess due to periodic boundary conditions 
			if ( dz >= HALF_LENGTH )
				dz -= LENGTH;
			else if ( dz < -HALF_LENGTH )
				dz += LENGTH;
			if ( dx >= HALF_LENGTH )
				dx -= LENGTH;
			else if ( dx < -HALF_LENGTH )
				dx += LENGTH;
			if ( dy >= HALF_LENGTH )
				dy -= LENGTH;
			else if ( dy < -HALF_LENGTH )
				dy += LENGTH;
			// end of mess 
			r = L*sqrt(dx*dx + dy*dy + dz*dz);
			bin_index = (int)(r/binSize);
			
			//calculate projection. here r can be abused
			rmn = sqrt(rij[DIM*l + X_COMP]*rij[DIM*l + X_COMP] + rij[DIM*l + Y_COMP]*rij[DIM*l + Y_COMP] + rij[DIM*l + Z_COMP]*rij[DIM*l + Z_COMP]);
			//get projection
			projection = 0.0;
			for (d=0;d<DIM;d++){
				zij = z[DIM*n + d] - z[DIM*m + d];
				projection += rij[DIM*l + d]*zij/rmn;
			}
			projections[bin_index][counters[bin_index]] = projection*projection;
			//fprintf(dumpfile,"%g %g\n",r,projections[bin_index][counters[bin_index]]);
			counters[bin_index]++;
		}
	}
	//get medians
	for (q=0;q<numOfBins;q++)
		c[q] += get_median(projections[q],counters[q]);
	
		
	for (q=0;q<numOfBins;q++)
		free(projections[q]);
	free(counters);
	free(projections);	
	//fclose(dumpfile);
	return;
}

int get_dipole_response(int contact_index, double *res){
        int a,i,j,q;
        double r;
        double *dipole_force = (double *)malloc(sizeof(double)*N*DIM);

        for (a=0;a<DIM*N;a++){
                dipole_force[a] = 0.0;
                res[a] = 0.0;
        }

        i = lookupBond[2*contact_index + I_INDEX];
        j = lookupBond[2*contact_index + J_INDEX];

        for (r=0.0,q=0;q<DIM;q++)
                r += rij[DIM*contact_index + q]*rij[DIM*contact_index + q];
        r = sqrt(r);

        for (q=0;q<DIM;q++){
                dipole_force[DIM*j + q] = rij[DIM*contact_index + q]/r;
                dipole_force[DIM*i + q] = -rij[DIM*contact_index + q]/r;
        }
        if ( cgSolver(res,dipole_force) ){
                normalizeVector(res);
                free(dipole_force);
                return 1;
	}
	free(dipole_force);
        return 0;
}


int main(int n,char **inputStrings){
	int q,l,numOfNetworks;
	int nx,ny,nz;
	char networkFileName[1024];
	char outFileName[128];
	double bulk_modulus, shear_modulus;
        double k_born, k_na, g_born, g_na;
	double coordination;
	FILE *outFile;
	char command[1024];
	//srand(time(0));
	
	//~ sscanf(inputStrings[1],"%d",&serial);
	//~ sscanf(inputStrings[2],"%d",&numOfNetworks);
	//~ sscanf(inputStrings[3],"%lf",&coordination);
	
	//srand( 1000 + serial );
	
	linearOperator = &hessianOnVector;
	//coordination = 6.64;
	//coordination = 7.28;
	//coordination = 8.56;
	coordination = 11.12;
	
	sprintf(outFileName,"G_elasticNetworks_N%d_z%.2f.dat", N, coordination );
	outFile = fopen(outFileName,"wb");
	
	for (q=serial*numOfNetworks; q < (serial+1)*numOfNetworks; q++){
		sprintf(networkFileName,"network_N%d_z%.2f_s%.4d.dat", N, coordination, q);
		if ( readNetwork(networkFileName) ){
			printf("read %s, coordination %g\n",networkFileName, (double)(2*contacts)/DOUBLE_N);
			
			calculateShearModulus(&g_born, &g_na);
			shear_modulus = g_born-g_na;
			
			fprintf(outFile,"%.11g\n", shear_modulus);
			
			/*calculateBulkModulus(&k_born, &k_na);
			bulk_modulus = (k_born-k_na)/( (double)(DIM*DIM) ) + pressure; */

		}
		else
			printf("%s not about\n",networkFileName);
	}
	fclose(outFile);
	//free stuff
	free_everything();
	return 0;
}
