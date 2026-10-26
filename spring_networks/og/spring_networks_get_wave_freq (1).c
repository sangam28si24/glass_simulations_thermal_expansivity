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
#define	TIME_STEP				0.1
#define	HALF_TIME_STEP			(0.5*TIME_STEP)


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

/********************* dont touch these!!! used as globals ******************************************/
int waves_checked[4096][3];
int num_of_waves_checked = 0;
/*************************************************************************************************************/

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

void free_everything(){
	free(rx); free(ry); free(rz); 
	free(lookupBond); free(second);
	free(first); free(rij); free(inv_r);
	return;
}



void build_transverse_wave(int nx, int ny, int nz, int mode_type, double *force){
	int i;
	double kx,ky,kz,mag_k;
	double perp1_x, perp1_y, perp1_z;
	double perp2_x, perp2_y, perp2_z;
	double perp_x, perp_y, perp_z;
	double ux,uy,uz,mag_u;
	double cos_theta,sin_theta,mag;
	
	
	kx = TWO_PI*( (double)nx )/L;
	ky = TWO_PI*( (double)ny )/L; 
	kz = TWO_PI*( (double)nz )/L;
	mag_k = sqrt(kx*kx + ky*ky + kz*kz);

	if ( ( nx == 0 ) && ( ny == 0 ) ){
		perp1_x = 1.0;
		perp1_y = 0.0;
		perp1_z = 0.0;
		perp2_x = 0.0;
		perp2_y = 1.0;
		perp2_z = 0.0;
	}
	else{
		//u is the cross-product between k and z
		ux = -ky;
		uy = kx;
		uz = 0.0;
		//normalize it
		mag_u = sqrt(ux*ux + uy*uy + uz*uz);
		ux = ux/mag_u;
		uy = uy/mag_u;
		uz = uz/mag_u;
		
		cos_theta = kz/mag_k;
		sin_theta = sqrt(1.0 - cos_theta*cos_theta);
		
		perp1_x = cos_theta + ux*ux*(1.0 - cos_theta);
		perp1_y = uy*ux*(1.0 - cos_theta) + uz*sin_theta;
		perp1_z = uz*ux*(1.0 - cos_theta) - uy*sin_theta;
		
		perp2_x = ux*uy*(1.0 - cos_theta) - uz*sin_theta;
		perp2_y = cos_theta + uy*uy*(1.0 - cos_theta);
		perp2_z = uz*uy*(1.0 - cos_theta) + ux*sin_theta;
	}
	//chose which transverse wave we want
	if ( mode_type ){
		perp_x = perp1_x;
		perp_y = perp1_y;
		perp_z = perp1_z;
	}
	else{
		perp_x = perp2_x;
		perp_y = perp2_y;
		perp_z = perp2_z;
	}
	
	
	for (i=0;i<N;i++){
		mag = sin(kx*(rx[i]*L) + ky*(ry[i]*L) + kz*(rz[i]*L));
		force[DIM*i + X_COMP] = perp_x*mag;
		force[DIM*i + Y_COMP] = perp_y*mag;
		force[DIM*i + Z_COMP] = perp_z*mag;
	}
	removeTranslation(force);
	return;
}




int add_wave(int nx, int ny, int nz){
	int w,found = 0;
	int minus_nx, minus_ny, minus_nz;
	
	minus_nx = -nx;
	minus_ny = -ny;
	minus_nz = -nz;
	
	for ( w=0; w<num_of_waves_checked; w++){
		if ( ( waves_checked[w][X_COMP] == minus_nx ) && ( waves_checked[w][Y_COMP] == minus_ny ) && ( waves_checked[w][Z_COMP] == minus_nz ) ){
			found = 1;
			break;
		}
	}
	if ( !found ){
		waves_checked[num_of_waves_checked][X_COMP] = nx;
		waves_checked[num_of_waves_checked][Y_COMP] = ny;
		waves_checked[num_of_waves_checked][Z_COMP] = nz;
		num_of_waves_checked++;
		return 1;
	}
	return 0;
}

void run_dynamics(double *initial_velocity, double duration, double *c){
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
	



int main(int n,char **inputStrings){
	int kyooz[32] = {48};
	int numOfKyooz = 1;
	int q,numOfNetworks,int_sum_of_sqrz,m,w,i,counter;
	int nx,ny,nz;
	int total_run_steps;
	double coordination,duration,frequency;
	double *z, *force, *p;
	double *c, *mean_c;
	double omega, affine_omega, overlap,k;
	char networkFileName[1024];
	char outFileName[1024];
	FILE *outFile;
	srand(time(0));
	
	sscanf(inputStrings[1],"%d",&serial);
	//~ sscanf(inputStrings[3],"%lf",&coordination);
	
	//srand( 1000 + serial );
	
	force = (double *)malloc(sizeof(double)*DIM*N);
	
	linearOperator = &hessianOnVector;
	//coordination = 6.64;
	//coordination = 7.28;
	//coordination = 8.56;
	//coordination = 11.12;
	coordination = 16.24;
	frequency = 0.75;
	duration = 5.0*TWO_PI/frequency;
	total_run_steps = (int)(duration/TIME_STEP);
	c = (double *)malloc(sizeof(double)*total_run_steps);
	mean_c = (double *)malloc(sizeof(double)*total_run_steps);
	for (w=0;w<total_run_steps;w++)
		mean_c[w] = 0.0;
	
	
	//sprintf(networkFileName,"/work/edanl/networks/N%d/network_N%d_z%.2f_s%.4d.dat",N,N,coordination,serial);
	sprintf(networkFileName,"network_N%d_z%.2f_s%.4d.dat",N,coordination,serial);
	if ( readNetwork(networkFileName) ){
		for (i=0;i<numOfKyooz;i++){
			q = kyooz[i];
			counter = 0;
			sprintf(outFileName,"waveFreqs_N%d_z%.2f_q%.4d_s%.4d.dat",N,coordination,q,serial);
			m = (int)(sqrt((double)q)+1.0);
			for (nx=-m;nx<=m;nx++){
				for (ny=-m;ny<=m;ny++){
					for (nz=-m;nz<=m;nz++){
						int_sum_of_sqrz = nx*nx + ny*ny + nz*nz;
						if ( int_sum_of_sqrz == q ){
							if ( add_wave( nx, ny, nz ) ){
								
								k = TWO_PI*sqrt((double)q)/L;
								
								build_transverse_wave( nx, ny, nz, 0, force);
								run_dynamics(force,duration,c); //put the correlation function in c
								counter++;
								for (w=0;w<total_run_steps;w++)
									mean_c[w] += c[w];
								
								printf("finished %d,%d,%d\n",nx,ny,nz);
							}
						}
					}
				}
			}
			for (w=0;w<total_run_steps;w++)
				mean_c[w] = mean_c[w]/(double)counter;
			//write
			outFile = fopen(outFileName,"wb");
			for (w=0;w<total_run_steps;w++)
				fprintf(outFile,"%.8g\t%.8g\n", TIME_STEP*(double)w, mean_c[w]);
			fclose(outFile);
		}
		free_everything();
	}
	else
		printf("%s not about\n",networkFileName);
	
	free(force);
	return 0;
}
