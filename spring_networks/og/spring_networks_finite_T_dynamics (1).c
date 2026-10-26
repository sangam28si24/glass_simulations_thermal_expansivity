#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define	N							4000
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
#define	TIME_STEP				0.1
#define	HALF_TIME_STEP			(0.5*TIME_STEP)
#define	TAU_B					1.0


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
double kinetic;
double T;


int contacts;

double *rx, *ry, *rz;
double *px, *py, *pz;
double *fx, *fy, *fz;

/************************ new bondnetwork nebz ******************************************************/
int *lookupBond;
double *second;
double *first;
double *rij;
double *inv_r;
double *restlength;
/*******************************************************************************************************/

void fixDrift(){
	int i;
	double Px=0.0,Py=0.0,Pz=0.0;
	for (i=0;i<N;i++){ Px += px[i]; Py += py[i]; Pz += pz[i]; }
	Px = Px/DOUBLE_N; Py = Py/DOUBLE_N; Pz = Pz/DOUBLE_N; 
	for (i=0;i<N;i++){ px[i] -= Px; py[i] -= Py; pz[i] -= Pz; }
	return;
}


void removeTranslation(double *v){
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
	return;
}

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
	px = (double *)malloc(sizeof(double)*N);
	py = (double *)malloc(sizeof(double)*N);
	pz = (double *)malloc(sizeof(double)*N);
	fx = (double *)malloc(sizeof(double)*N);
	fy = (double *)malloc(sizeof(double)*N);
	fz = (double *)malloc(sizeof(double)*N);
	
	
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
		restlength = (double *)malloc(sizeof(double)*contacts);
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
			restlength[l] = sqrt(r2);
		}
		fclose(inFile);
		return 1;
	}
	printf("couldn't read %s\n",inFileName);
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
		//~ if ( maxNorm < CGSOLVER_TOL ) // || ( (typicalNorm < CGSOLVER_TOL) && (maxNorm < 10.0*CGSOLVER_TOL) )  )
			//~ break;
		
		if ( typicalNorm/scale < CGSOLVER_TOL )
			break;
		
		//~ if ( k>0 && !(k%16384) )
			//~ printf("it%d, typicalNorm %g, maxNorm %g\n",k,typicalNorm,maxNorm);
		
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



void run_harmonic_dynamics(double *initial_velocity, double duration, double *c){
	double *z, *p, *hesOnZ;
	double dz,this_c;
	int a,q,total_run_steps;
	
	total_run_steps = (int)(duration/TIME_STEP);
	z = (double *)malloc(sizeof(double)*DIM*N);
	p = (double *)malloc(sizeof(double)*DIM*N);
	hesOnZ = (double *)malloc(sizeof(double)*DIM*N);
	for (a=0;a<DIM*N;a++){
		p[a] = initial_velocity[a];
		z[a] = 0.0;
	}
	hessianOnVector(z,hesOnZ);
	
	c[0] = 1.0;
	for (q=1;q<total_run_steps;q++){
		for (a=0;a<DIM*N;a++){
			p[a] -= HALF_TIME_STEP*hesOnZ[a];
			z[a] += TIME_STEP*p[a];		
		}
		hessianOnVector(z,hesOnZ);
		for (a=0;a<DIM*N;a++)
			p[a] -= HALF_TIME_STEP*hesOnZ[a];
		
		//calculate correlations
		c[q] = cdot(p,initial_velocity)/( sqrt(cdot(initial_velocity,initial_velocity))*sqrt(cdot(p,p)) );
	}
	free(z); free(p); free(hesOnZ);
	return;
}

void calculateForces(){
	int l,i,j;
	double dx,dy,dz,r,g,temp;
	
	for (i=0;i<N;i++){ fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;}
	
	for (l=0;l<contacts;l++){
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
		r = L*sqrt( dx*dx + dy*dy + dz*dz );
		g = (restlength[l]-r)/r;
		temp = L*dx*g;
		fx[j] += temp;
		fx[i] -= temp;
		temp = L*dy*g;
		fy[j] += temp;
		fy[i] -= temp;
		temp = L*dz*g;
		fz[j] += temp;
		fz[i] -= temp;
		u += 0.5*(restlength[l] - r)*(restlength[l] - r);
		pressure += g*r*r; // notice that that g = f/r, so this is f x r
	}
	return;
}



void advanceTimeBerendsen(){
	int i;
	double temp,dx,dy,dz,inst_T,xi_T;
	
	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		pz[i] += HALF_TIME_STEP*fz[i];
		
		dx = TIME_STEP*px[i]; //real space displacement
		temp = rx[i] + dx*invL;
		if (temp >= LENGTH)
			rx[i] = temp - LENGTH;
		else if (temp < 0.0)
			rx[i] = temp + LENGTH;
		else
			rx[i] = temp;
		
		dy = TIME_STEP*py[i]; //real space displacement
		temp = ry[i] + dy*invL;
		if (temp >= LENGTH)
			ry[i] = temp - LENGTH;
		else if (temp < 0.0)
			ry[i] = temp + LENGTH;
		else
			ry[i] = temp;
		
		dz = TIME_STEP*pz[i]; //real space displacement
		temp = rz[i]  + dz*invL;
		if (temp >= LENGTH)
			rz[i] = temp - LENGTH;
		else if (temp < 0.0)
			rz[i] = temp + LENGTH;
		else
			rz[i] = temp;
	}
	
	calculateForces();
	kinetic = 0.0;
	for (i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		pz[i] += HALF_TIME_STEP*fz[i];
		kinetic += px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
	}
	kinetic = 0.5*kinetic;
	
	/************************** thermostating part *************************************/
	inst_T = 2.0*kinetic/( ((double)DIM)*DOUBLE_N );
	xi_T = sqrt(1.0 + TIME_STEP*(T/inst_T - 1.0)/TAU_B);
	for (i=0;i<N;i++){ px[i] *= xi_T; py[i] *= xi_T; pz[i] *= xi_T; }
	/********************************************************************************/
	return;
}


void run_berendsen(double duration){
        int k,steps,between=128;
        
        steps = (int)(duration/TIME_STEP);
        
        for (k=0; k<steps; k++){
                advanceTimeBerendsen();
		if ( !(k%between) )
                	printf("%.8g\t%.8g\t%.8g\n",(double)k*TIME_STEP, kinetic/DOUBLE_N, u/DOUBLE_N );
                
		if ( !(k%16384) )
                        fixDrift();
        }
        return;
}

void free_everything(){
	free(rx); free(ry); free(rz); 
	free(px); free(py); free(pz); 
	free(fx); free(fy); free(fz);
	free(lookupBond); free(second);
	free(first); free(rij); free(inv_r); free(restlength);
	return;
}



int main(int n,char **inputStrings){
	char networkFileName[1024];
	char outFileName[1024];
	double coordination;
	
	srand(time(0));
	
	//~ sscanf(inputStrings[1],"%d",&serial);
	//~ sscanf(inputStrings[3],"%lf",&coordination);
	
	linearOperator = &hessianOnVector;
	
	coordination = 16.24;
	
	sprintf(networkFileName,"network_N%d_z%.2f_s%.4d.dat",N,coordination,serial);
	if ( readNetwork(networkFileName) ){
		
		
		
		free_everything();
	}
	else
		printf("%s not about\n",networkFileName);
	
	return 0;
}
