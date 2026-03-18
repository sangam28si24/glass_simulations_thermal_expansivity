#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


//all these are for initializing an FCC
//notice that here N = 4n^3 with n some integer. here n^3 = 1000 = 10^3 giving N=4000
//notice that here N = 4n^3 with n some integer. here n^3 = 343 = 7^3 giving N=1372
//notice that here N = 4n^3 with n some integer. here n^3 = 1728 = 12^3 giving N=6912
//end of FCC numbers


#define	N							4000
#define	DOUBLE_N						((double)N)
#define	DIM							3
#define	CGSOLVER_TOL				1e-10


/* stuff that we don't touch usually */
#define	TIME_STEP                                    	0.001
#define	HALF_TIME_STEP                          (0.5*TIME_STEP)
#define	DELTA(a,b)					((a)==(b)?1.0:0.0)
#define	MAX_NEBZ						256
#define	CONTACTS					32
#define	MAX_CELLS					(N)
#define	X_COMP					0
#define 	Y_COMP 					1
#define 	Z_COMP 					2
#define	I_INDEX					0
#define	J_INDEX					1
#define	CUTOFF					(2.5)
#define	CUTOFF_SQRD					(CUTOFF * CUTOFF)
#define	C_0							(0.075421526038318)
#define	C_2							(-0.035416041040039)
#define	C_4							(0.005708815153688)
#define	C_6							(-0.000318029410118)
#define	N_REP						12.0
#define	N_ATT						6.0
#define	DR							0.01


/* constant constants... :P*******************************************************/
#define	LENGTH				1.0
#define	HALF_LENGTH			0.5

#define	TWO_PI					6.28318530717958
#define	UNI		((double)rand()/((double)RAND_MAX + 1.0))
#define	HUGE_INT				10000000

void initializeSystem();
void updateNebzLists();
void calculateEverything();


int serial;

/*globals*/
double L; 			/*	length of <square> box 	*/
double invL;			/* inverse length of box	*/
double V;			/*	volume			*/
double T;			/*	temperature			*/
double t;			/*	time				*/
double u;			/*	potential energy 		*/
double kinetic;		/*	kinetic energy 		*/
double stress;
double pressure;
double typicalGrad;
double typicalInteractionStrength;

int contacts;

double *rx, *ry, *rz;	
double *fx, *fy, *fz;
int *type; 			/*	for binary system	*/

/******************************* for nebz list & cell subdivision *******************************/
int **nebz;
int *numOfNebz;
int *numInCell;
int **cells;
double cellCutOff,maxD,listCutOffSqrd;
int nebListCounter;
const int shift[13][3] = {{-1,-1,1},{-1,0,1},{-1,1,1},{0,-1,1},{0,0,1},{0,1,1},{1,-1,1},{1,0,1},{1,1,1},{-1,-1,0},{-1,0,0},{-1,1,0},{0,1,0}};
/************************************************************************************/

/************************ new bondnetwork nebz ******************************************************/
int *lookupBond;
double *third;
double *second;
double *first;
double *rij;
/*******************************************************************************************************/

const double invSizeSqrd[2][2] = {{1.00,0.718184429761563},{0.718184429761563,0.510204081632653}}; // [0][0] = small-small; [0][1] = [1][0] = large-small; [1][1] = large-large;

double cdot( double *vec1, double *vec2 ){
	double res = 0.0;
	int i;
	for (i=0; i<DIM*N; i++)
		res += vec1[i]*vec2[i];
	return res;
}


void initializeSystem(){
	cellCutOff = 1.4*sqrt(CUTOFF_SQRD) + DR;
	listCutOffSqrd = cellCutOff*cellCutOff;
	updateNebzLists();
	calculateEverything();
	return;
}


int readSnapShot(char *fileName){
	int i,k,q;
	double dummy;
	FILE *snapFile;
	snapFile = fopen(fileName, "rb");
	if ( snapFile ){
		q = fscanf(snapFile, "%lf %lf %lf %lf",&L,&dummy,&dummy,&dummy);
		V = L*L*L; invL = 1.0/L;
		for (i=0;i<N;i++)
			q = fscanf(snapFile, "%lf %lf %lf %d",&(rx[i]),&(ry[i]),&(rz[i]),&(type[i]));
		fclose(snapFile);
		initializeSystem();
		return 1;
	}
	return 0;
}


//this uses the cell subdivision.
void updateNebzLists(){
	int i,j,x,y,z,current,m,m2,m3;
	int a,b,c,k,l,numHere,numThere,target,w,xIndex,yIndex;
	double dx,dy,dz,r2,invCellSize;
	
	nebListCounter = 0;
	maxD = 0.0;
	
	m =(int)(L/cellCutOff);
	m3 = m*m*m; 		//the length of cells and numInCell...
	m2 = m*m;
	invCellSize = (double)m; 	// for reduced coordinates...
	
	//set number in each cell to zero
	for (i=0; i<m3; i++)
		numInCell[i] = 0;
	
	//********cells sorting********
	for (i=0;i<N;i++){
		current = m2*(int)(rz[i]*invCellSize) + m*(int)(ry[i]*invCellSize) + (int)(rx[i]*invCellSize);
		cells[current][numInCell[current]] = i;
		numInCell[current]++;
	}
	
	//set all neb lists to zero
	for (i=0; i<N; i++)
		numOfNebz[i] = 0;
	
	//z is layer... 
	for (z=0;z<m;z++){
		for (y=0;y<m;y++){
			for (x=0;x<m;x++){
				current = m2*z + m*y + x; //this cell index.
				numHere = numInCell[current];
				//first check interactions within a cell
				for (k=0;k<numHere-1;k++){
					for (l=k+1;l<numHere;l++){
						i = cells[current][k];
						j = cells[current][l];
						dx = rx[j] - rx[i];
						dy = ry[j] - ry[i];
						dz = rz[j] - rz[i];
						r2 = L*L*( dx*dx + dy*dy + dz*dz );
						if (r2 < listCutOffSqrd){
							nebz[i][numOfNebz[i]] = j;
							nebz[j][numOfNebz[j]] = i;
							numOfNebz[i]++;
							numOfNebz[j]++;
						}
					}
				}
			
				//shift[13][3] = {{-1,-1,1},{-1,0,1},{-1,1,1},{0,-1,1},{0,0,1},{0,1,1},{1,-1,1},{1,0,1},{1,1,1},{-1,-1,0},{-1,0,0},{-1,1,0},{0,1,0}}
				//now there are 13 nebz that need to be checked.
				for (w=0;w<13;w++){ 
					a = x + shift[w][0];
					if (a<0) a = m-1;
					else if (a==m) a = 0;
					b = y + shift[w][1];
					if (b<0) b = m-1;
					else if (b==m) b = 0;
					c = z + shift[w][2];
					if (c<0) c = m-1;
					else if (c==m) c = 0;
					target = a + m*b + m2*c;
					numThere = numInCell[target];
					for (k=0;k<numHere;k++){
						for (l=0;l<numThere;l++){
							i = cells[current][k];
							j = cells[target][l];
							dx = rx[j] - rx[i];
							dy = ry[j] - ry[i];
							dz = rz[j] - rz[i];
							// mess due to periodic boundary conditions 
							if ( dx >= HALF_LENGTH )
								dx -= LENGTH;
							else if ( dx < -HALF_LENGTH )
								dx += LENGTH;
							if ( dy >= HALF_LENGTH )
								dy -= LENGTH;
							else if ( dy < -HALF_LENGTH )
								dy += LENGTH;
							if ( dz >= HALF_LENGTH )
								dz -= LENGTH;
							else if ( dz < -HALF_LENGTH )
								dz += LENGTH;
							// end of mess 
							r2 = L*L*( dx*dx + dy*dy + dz*dz );
							if (r2 < listCutOffSqrd){
								nebz[i][numOfNebz[i]] = j;
								nebz[j][numOfNebz[j]] = i;
								numOfNebz[i]++;
								numOfNebz[j]++;
							}
						}
					}
				}
			}
		}
	}
	return;
}

void calculateEverything(){
	int i,j,m,k;
	double goo,foo,dx,dy,dz,invSigmaSqrd,r2OverSigma2;
	double sigma2OverR2,sigma6OverR6,sigma8OverR8,temp;
		
	//set stuff to zero
	u = 0.0; stress = 0.0; pressure = 0.0; contacts = 0; typicalInteractionStrength = 0.0;
	for (i=0;i<N;i++){ fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0; }
	
	for (i=0;i<N;i++){
		m = numOfNebz[i];
		for (k=0;k<m;k++){
			j = nebz[i][k];
			if (j > i){
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
				invSigmaSqrd = invSizeSqrd[type[i]][type[j]];
				r2OverSigma2 = L*L*( dx*dx + dy*dy + dz*dz )*invSigmaSqrd;
				if ( r2OverSigma2 < CUTOFF_SQRD){
					lookupBond[2*contacts + I_INDEX] = i;
					lookupBond[2*contacts + J_INDEX] = j;
					rij[DIM*contacts + X_COMP] = L*dx;
					rij[DIM*contacts + Y_COMP] = L*dy;
					rij[DIM*contacts + Z_COMP] = L*dz;
					
					sigma2OverR2 = 1.0/r2OverSigma2;
					sigma6OverR6 = sigma2OverR2*sigma2OverR2*sigma2OverR2;
					sigma8OverR8 = sigma6OverR6*sigma2OverR2; 
					goo = 4.0*(sigma8OverR8*(N_REP*sigma6OverR6 - N_ATT) - 2.0*C_2 - r2OverSigma2*(4.0*C_4 + 6.0*r2OverSigma2*C_6) )*invSigmaSqrd;
					
					third[contacts] = 4.0*( -2688.0*sigma8OverR8*sigma8OverR8*sigma2OverR2 + 480.0*sigma6OverR6*sigma6OverR6 + 48.0*C_6)*(invSigmaSqrd*invSigmaSqrd*invSigmaSqrd);
					second[contacts] = 4.0*( 168.0*sigma8OverR8*sigma8OverR8 - 48.0*sigma8OverR8*sigma2OverR2 + 24.0*C_6*r2OverSigma2 + 8.0*C_4)*(invSigmaSqrd*invSigmaSqrd);
					first[contacts] = -goo;
					contacts++;
					
					pressure += goo*L*L*( dx*dx + dy*dy + dz*dz ); // notice that that g = f/r, so this is f x r
					typicalInteractionStrength += goo*goo*L*L*( dx*dx + dy*dy + dz*dz ); // notice that g = f/r, so this is f^2
					
					temp = L*dx*goo;
					fx[j] += temp;
					fx[i] -= temp;
					temp = L*dy*goo;
					fy[j] += temp;
					fy[i] -= temp;
					stress += -temp*L*dx;
					temp = L*dz*goo;
					fz[j] += temp;
					fz[i] -= temp;
					
					u += 4.0*(sigma6OverR6*(sigma6OverR6 - 1.0) + r2OverSigma2*(C_2 + r2OverSigma2*(C_4 + r2OverSigma2*C_6) ) + C_0);
				}
			}
		}
	}
	stress = stress/V;
	pressure = pressure/( V*(double)DIM);
	typicalInteractionStrength = sqrt(typicalInteractionStrength/DOUBLE_N);
	typicalGrad = 0.0;
	for (i=0;i<N;i++)
		typicalGrad += fx[i]*fx[i] + fy[i]*fy[i] + fz[i]*fz[i];
	typicalGrad = sqrt(typicalGrad/DOUBLE_N);
	return;
}



void allocate_memory(){
	int i;
	
	rx = (double *)malloc(sizeof(double)*N);
	ry = (double *)malloc(sizeof(double)*N);
	rz = (double *)malloc(sizeof(double)*N);
	fx = (double *)malloc(sizeof(double)*N);
	fy = (double *)malloc(sizeof(double)*N);
	fz = (double *)malloc(sizeof(double)*N);
	type = (int *)malloc(sizeof(int)*N);
		
	numOfNebz = (int *)malloc(sizeof(int)*N);
	nebz = (int **)malloc(sizeof(int *)*N);
	for (i=0;i<N;i++)
		nebz[i] = (int *)malloc(sizeof(int)*MAX_NEBZ);
	
	numInCell = (int *)malloc(sizeof(int)*MAX_CELLS);
	cells = (int **)malloc(sizeof(int *)*MAX_CELLS);
	for (i=0;i<MAX_CELLS;i++)
		cells[i] = (int *)malloc(sizeof(int)*MAX_NEBZ);
		
	lookupBond = (int *)malloc(sizeof(int)*N*CONTACTS*2);
	third = (double *)malloc(sizeof(double)*N*CONTACTS);
	second = (double *)malloc(sizeof(double)*N*CONTACTS);
	first = (double *)malloc(sizeof(double)*N*CONTACTS);
	rij = (double *)malloc(sizeof(double)*N*CONTACTS*DIM);
	
	return;
}

void free_everything(){
	int i;
	free(rx); free(ry); free(rz); 
	free(fx); free(fy); free(fz);
	free(type);
	
	for (i=0;i<N;i++)
		free(nebz[i]);
	free(nebz);
	free(numOfNebz);
	
	for (i=0;i<MAX_CELLS;i++)
		free(cells[i]);
	free(cells);
	free(numInCell);
	free(lookupBond); free(third); free(second); free(first); free(rij); 
	return;
}

/******************************all energy derivatives************************************************/

double u_gamma(){
	int l;
	double res = 0.0;
	
	for (l=0;l<contacts;l++)
		res += first[l]*rij[DIM*l+X_COMP]*rij[DIM*l+Y_COMP];
	
	return res;
}

double u_gamma_gamma(){
	int l;
	double res = 0.0;
	
	for (l=0;l<contacts;l++){
		res += second[l]*rij[DIM*l+X_COMP]*rij[DIM*l+X_COMP]*rij[DIM*l+Y_COMP]*rij[DIM*l+Y_COMP];
		res += first[l]*rij[DIM*l+Y_COMP]*rij[DIM*l+Y_COMP];
	}
	
	return res;
}

void u_gamma_x(double *res){
	int q,i,j,k,m,l=0;
	double rx,ry;
	
	for (i=0;i<N;i++)
		for (q=0;q<DIM;q++)
			res[DIM*i + q] = 0.0;
		
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		rx = rij[DIM*l + X_COMP];
		ry = rij[DIM*l + Y_COMP];
		for (q=0;q<DIM;q++){
			res[DIM*i + q] -= second[l]*rx*ry*rij[DIM*l + q];
			res[DIM*j + q] += second[l]*rx*ry*rij[DIM*l + q];
			
			res[DIM*i + q] -= first[l]*ry*DELTA(q,X_COMP);
			res[DIM*j + q] += first[l]*ry*DELTA(q,X_COMP);
		}
	}
	return;
}


double u_gamma_x_on_vector(double *vec){
	int q,i,j,l=0;
	double rx,ry;
	double res = 0.0;
	double vij[DIM];
	double v_dot_r;
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		rx = rij[DIM*l + X_COMP];
		ry = rij[DIM*l + Y_COMP];
		
		for (v_dot_r=0.0,q=0;q<DIM;q++){
			vij[q] = vec[DIM*j + q] - vec[DIM*i + q];
			v_dot_r += vij[q]*rij[DIM*l + q];
		}
		
		res += second[l]*rx*ry*v_dot_r + first[l]*ry*vij[X_COMP];
	}
	return res;
}

double u_gamma_eta(){
	int l,d;
	double xij,yij,r2,res = 0.0;
	
	for (l=0;l<contacts;l++){
		for (r2=0.0,d=0;d<DIM;d++)
			r2 += rij[DIM*l + d]*rij[DIM*l + d];
		xij = rij[DIM*l + X_COMP];
		yij = rij[DIM*l + Y_COMP];
		
		res += second[l]*r2*xij*yij + first[l]*2.0*xij*yij;
	}
	
	return res;
}



double u_gamma_gamma_eta(){
	int l,d;
	double xij,yij,r2,res = 0.0;
	
	for (l=0;l<contacts;l++){
		for (r2=0.0,d=0;d<DIM;d++)
			r2 += rij[DIM*l + d]*rij[DIM*l + d];
		xij = rij[DIM*l + X_COMP];
		yij = rij[DIM*l + Y_COMP];
		
		res += third[l]*r2*xij*xij*yij*yij;
		
		res += second[l]*(4.0*xij*xij*yij*yij + r2*yij*yij);
		
		res += first[l]*2.0*yij*yij;
	}
		
	return res;
}

double u_gamma_eta_x_on_vector(double *vec){
	double vxij,v_dot_r,xij,yij,r2,res = 0.0;
	int l,i,j,d;
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		r2=0.0; v_dot_r = 0.0;
		for (d=0;d<DIM;d++){
			r2 += rij[DIM*l + d]*rij[DIM*l + d];
			v_dot_r += rij[DIM*l + d]*(vec[DIM*j+d]-vec[DIM*i+d]);
		}
		xij = rij[DIM*l + X_COMP];
		yij = rij[DIM*l + Y_COMP];
		vxij = vec[DIM*j+X_COMP]-vec[DIM*i+X_COMP];
		
		res += third[l]*r2*xij*yij*v_dot_r + second[l]*(3.0*xij*yij*v_dot_r + r2*yij*vxij) + first[l]*yij*vxij;
	}
	return res;
}

double u_gamma_gamma_x_on_vector(double *vec){
	double vxij,v_dot_r,xij,yij,r2,res = 0.0;
	int l,i,j,d;
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		r2=0.0; v_dot_r = 0.0;
		for (d=0;d<DIM;d++){
			r2 += rij[DIM*l + d]*rij[DIM*l + d];
			v_dot_r += rij[DIM*l + d]*(vec[DIM*j+d]-vec[DIM*i+d]);
		}
		xij = rij[DIM*l + X_COMP];
		yij = rij[DIM*l + Y_COMP];
		vxij = vec[DIM*j+X_COMP]-vec[DIM*i+X_COMP];
		
		res += third[l]*xij*xij*yij*yij*v_dot_r + second[l]*(yij*yij*v_dot_r + 2.0*xij*yij*yij*vxij);
	}
	return res;
}


double u_gamma_x_x_on_two_vectors(double *vec1, double *vec2){
	double v1xij,v2xij,v1_dot_r,v2_dot_r,v1_dot_v2,xij,yij,r2,res = 0.0;
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
		xij = rij[DIM*l + X_COMP];
		yij = rij[DIM*l + Y_COMP];
		v1xij = vec1[DIM*j+X_COMP]-vec1[DIM*i+X_COMP];
		v2xij = vec2[DIM*j+X_COMP]-vec2[DIM*i+X_COMP];
		
		res += third[l]*xij*yij*v1_dot_r*v2_dot_r + second[l]*(yij*(v1xij*v2_dot_r + v2xij*v1_dot_r) + xij*yij*v1_dot_v2);
	}
	return res;
}


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
		res += third[l]*r2*v1_dot_r*v2_dot_r + second[l]*(2.0*v1_dot_r*v2_dot_r + r2*v1_dot_v2);
	}
	return res;
}




	
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

void tessianOnTwoVectors(double *v, double *u, double *res){
	int q,i,j,l=0;
	double vij[DIM],uij[DIM];
	double v_dot_r, u_dot_r, v_dot_u;
	
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
		
		
		for (u_dot_r=0.0,q=0;q<DIM;q++)
			u_dot_r += rij[DIM*l + q]*uij[q];
		
		for (q=0;q<DIM;q++){
			res[DIM*i + q] -= third[l]*u_dot_r*v_dot_r*rij[DIM*l + q] + second[l]*( v_dot_u*rij[DIM*l + q] + v_dot_r*uij[q] + u_dot_r*vij[q] );
			res[DIM*j + q] += third[l]*u_dot_r*v_dot_r*rij[DIM*l + q] + second[l]*( v_dot_u*rij[DIM*l + q] + v_dot_r*uij[q] + u_dot_r*vij[q] );
		}
	}
	return;
}

double tessianOnThreeVectors(double *v, double *u, double *w){
	int q,i,j,k,m,l=0;
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


/******************************end of all energy derivatives************************************************/

int cgSolver(double *x, double *b){
	double *r, *w, *z;
	double a,B,maxNorm,sum,scale;
	int i,j,k,t_0,q;
	
	t_0 = (int)time(0);
	
	scale = sqrt(cdot(b,b)/( (double)DIM*DOUBLE_N ) );
	
	r = (double *)malloc(sizeof(double)*DIM*N);
	w = (double *)malloc(sizeof(double)*DIM*N);
	z = (double *)malloc(sizeof(double)*DIM*N);
	
	hessianOnVector(x,r);
	
	for (i=0;i<DIM*N;i++){
		r[i] = b[i] - r[i];
		w[i] = -r[i];
	}
	
	hessianOnVector(w,z);
	a = cdot(r,w)/cdot(w,z);
	
	for (i=0;i<DIM*N;i++)
		x[i] += a*w[i];
	
	B = 0.0;
	
	for (k=0; k<HUGE_INT; k++){
		maxNorm = 0.0;
		for (i=0;i<DIM*N;i++){
			r[i] = r[i] - a*z[i];
			if (r[i]*r[i] > maxNorm)
				maxNorm = r[i]*r[i];
		}
		maxNorm = sqrt(maxNorm);
		
		//~ if ( !(k%32) )
			//~ printf("maxNorm:%g, scale:%g, ratio:%g\n",maxNorm,scale, maxNorm/scale);
			
		if ( maxNorm/scale < CGSOLVER_TOL ) // make sure this is the correct stopping condition!
			break;
		
		B = cdot(r,z)/cdot(w,z);
		for (i=0;i<DIM*N;i++)
			w[i] = -r[i] + B*w[i];
		
		hessianOnVector(w,z);
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

	calculateEverything();
	
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
				a = DIM*i + alpha;
				b = DIM*j + beta;
				fprintf(outFile,"%d\t%d\t%.15g\n%d\t%d\t%.15g\n",a,b,thisHes,b,a,thisHes);
				diagonal[DIM*DIM*i + m] += second[l]*rij[DIM*l + alpha]*rij[DIM*l + beta] + first[l]*DELTA(alpha,beta);
				diagonal[DIM*DIM*j + m] += second[l]*rij[DIM*l + alpha]*rij[DIM*l + beta] + first[l]*DELTA(alpha,beta);
				m++;
			}
		}
	}
	//diagonal blocks
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
	free(diagonal);
	fclose(outFile);
	return;
}


int compression_nonaffine_velocities(double *res){
	int a,i,j,l,d;
	double r2;
	double *xi = (double *)malloc(sizeof(double)*N*DIM);
	
	for (a=0;a<DIM*N;a++){
		res[a] = 0.0;
		xi[a] = 0.0;
	}
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		
		for (r2 = 0.0,d=0;d<DIM;d++)
			r2 += rij[DIM*l+d]*rij[DIM*l+d];
		
		for (d=0;d<DIM;d++){
			xi[DIM*j + d] += (second[l]*r2 + first[l])*rij[DIM*l + d];
			xi[DIM*i + d] -= (second[l]*r2 + first[l])*rij[DIM*l + d];
		}
	}
	
	if ( cgSolver(res,xi) ){
		//invert sign
		for (a=0;a<DIM*N;a++)
			res[a] = -res[a];
		free(xi);
		return 1;
	}
	free(xi);
	return 0;
}


int shear_nonaffine_velocities(double *res){
	int a;
	double *xi = (double *)malloc(sizeof(double)*N*DIM);
	
	for (a=0;a<DIM*N;a++){
		res[a] = 0.0;
		xi[a] = 0.0;
	}
	u_gamma_x(xi);
	
	if ( cgSolver(res,xi) ){
		for (a=0;a<DIM*N;a++)
			res[a] = -res[a];
		free(xi);
		return 1;
	}
	free(xi);
	return 0;
}


//notice that g_na calculated and returned is POSITIVE
int get_shear_modulus(double *g_born, double *g_na){
	int i,q;
        double *nav, *xi;
	
	nav = (double *)malloc(sizeof(double)*N*DIM);
        xi = (double *)malloc(sizeof(double)*N*DIM);
	
	for (i=0;i<N;i++)
                for (q=0;q<DIM;q++)
                        nav[DIM*i+q] = 0.0;
	
	u_gamma_x(xi);
	if ( cgSolver(nav,xi) ){
                *g_na = cdot(nav,xi)/V;
                *g_born = u_gamma_gamma()/V;
                free(nav); free(xi);
                return 1;
        }
        printf("failed to calculate mu\n");
	free(nav); free(xi);
        return 0;
}

//notice that k_na calculated and returned is POSITIVE
int get_bulk_modulus(double *k_born, double *k_na){
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

                *k_born += (second[l]*r2 + first[l])*r2 + first[l]*r2; //check notes, this looks weird but its correct

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

double get_second_order_dilatancy(){
	double res = 0.0;
	double *vna_eta, *vna_gamma;
	
	vna_eta = (double *)malloc(sizeof(double)*DIM*N);
	vna_gamma = (double *)malloc(sizeof(double)*DIM*N);
	
	shear_nonaffine_velocities(vna_gamma);
	compression_nonaffine_velocities(vna_eta);
	
	res += u_gamma_gamma_eta();
	//printf("1: %g\n",u_gamma_gamma_eta()/DOUBLE_N);
	
	res += 2.0*u_gamma_eta_x_on_vector(vna_gamma);
	//printf("2: %g\n",2.0*u_gamma_eta_x_on_vector(vna_gamma)/DOUBLE_N);
	
	res += u_gamma_gamma_x_on_vector(vna_eta);
	//printf("3: %g\n",u_gamma_gamma_x_on_vector(vna_eta)/DOUBLE_N);
	
	res += u_eta_x_x_on_two_vectors(vna_gamma, vna_gamma);
	//printf("4: %g\n", u_eta_x_x_on_two_vectors(vna_gamma, vna_gamma)/DOUBLE_N);
	
	res += 2.0*u_gamma_x_x_on_two_vectors(vna_eta,vna_gamma);
	//printf("5: %g\n", 2.0*u_gamma_x_x_on_two_vectors(vna_eta,vna_gamma)/DOUBLE_N);
	
	res += tessianOnThreeVectors(vna_eta,vna_gamma,vna_gamma);
	//printf("6: %g\n", tessianOnThreeVectors(vna_eta,vna_gamma,vna_gamma)/DOUBLE_N);
	
	free(vna_eta); free(vna_gamma);
	return -res/(V*(double)DIM);
}

/* Write H (COO sparse) and one CG solve to text files for Python cross-check.
 * Writes: hessian_check_H.txt   (row col val, one triplet per line)
 *         hessian_check_b.txt   (3N values, one per line)
 *         hessian_check_Hinvb.txt  (3N values: w = H^-1 * b) */
void write_h_check(int g){
	int i,l,m,alpha,beta,a,b;
	double *diagonal, *bvec, *wvec;
	double thisHes;
	char nameH[64], nameb[64], namew[64];
	FILE *fH, *fb, *fw;

	sprintf(nameH, "hessian_check_H_g%d.txt",    g);
	sprintf(nameb, "hessian_check_b_g%d.txt",    g);
	sprintf(namew, "hessian_check_Hinvb_g%d.txt", g);

	/* ── write H as sparse COO (same logic as writeHessian) ── */
	fH = fopen(nameH, "w");
	fb = fopen(nameb, "w");
	fw = fopen(namew, "w");
	diagonal = (double *)malloc(sizeof(double)*N*DIM*DIM);
	for (i=0;i<N*DIM*DIM;i++) diagonal[i] = 0.0;

	for (l=0;l<contacts;l++){
		int ii = lookupBond[2*l + I_INDEX];
		int jj = lookupBond[2*l + J_INDEX];
		m = 0;
		for (alpha=0;alpha<DIM;alpha++){
			for (beta=0;beta<DIM;beta++){
				thisHes = -second[l]*rij[DIM*l+alpha]*rij[DIM*l+beta]
				          -first[l]*DELTA(alpha,beta);
				a = DIM*ii + alpha;
				b = DIM*jj + beta;
				fprintf(fH,"%d\t%d\t%.15g\n", a, b, thisHes);
				fprintf(fH,"%d\t%d\t%.15g\n", b, a, thisHes);
				diagonal[DIM*DIM*ii + m] +=  second[l]*rij[DIM*l+alpha]*rij[DIM*l+beta] + first[l]*DELTA(alpha,beta);
				diagonal[DIM*DIM*jj + m] +=  second[l]*rij[DIM*l+alpha]*rij[DIM*l+beta] + first[l]*DELTA(alpha,beta);
				m++;
			}
		}
	}
	for (i=0;i<N;i++){
		m=0;
		for (alpha=0;alpha<DIM;alpha++){
			for (beta=0;beta<DIM;beta++){
				a = DIM*i+alpha; b = DIM*i+beta;
				fprintf(fH,"%d\t%d\t%.15g\n", a, b, diagonal[DIM*DIM*i+m]);
				m++;
			}
		}
	}
	free(diagonal);
	fclose(fH);

	/* ── fixed random b, solve H*w = b, write both ── */
	bvec = (double *)malloc(sizeof(double)*DIM*N);
	wvec = (double *)malloc(sizeof(double)*DIM*N);
	srand(12345);
	for (i=0;i<DIM*N;i++) bvec[i] = (double)rand()/RAND_MAX - 0.5;
	/* project out translations so cgSolver converges */
	{
		double sx=0,sy=0,sz=0;
		for (i=0;i<N;i++){sx+=bvec[DIM*i];sy+=bvec[DIM*i+1];sz+=bvec[DIM*i+2];}
		sx/=N; sy/=N; sz/=N;
		for (i=0;i<N;i++){bvec[DIM*i]-=sx;bvec[DIM*i+1]-=sy;bvec[DIM*i+2]-=sz;}
	}
	for (i=0;i<DIM*N;i++) wvec[i]=0.0;
	cgSolver(wvec, bvec);

	for (i=0;i<DIM*N;i++) fprintf(fb, "%.15g\n", bvec[i]);
	for (i=0;i<DIM*N;i++) fprintf(fw, "%.15g\n", wvec[i]);
	fclose(fb); fclose(fw);
	free(bvec); free(wvec);
	printf("  wrote hessian_check_H.txt, hessian_check_b.txt, hessian_check_Hinvb.txt\n");
}

/* Write H (sparse COO) and vna_eta for glass g so Python can diagonalise H,
 * compute the mode sums, and evaluate alpha from eq. 24. */
void write_for_python(int g){
	int i, l, m, alpha, beta, a, b;
	double *diagonal, *vna_eta;
	double thisHes;
	char nameH[64], nameV[64];
	FILE *fH, *fV;

	sprintf(nameH, "hessian_g%d.txt", g);
	sprintf(nameV, "vna_eta_g%d.txt", g);
	fH = fopen(nameH, "w");
	fV = fopen(nameV, "w");

	/* H sparse COO — same logic as writeHessian */
	diagonal = (double *)calloc(N*DIM*DIM, sizeof(double));
	for (l = 0; l < contacts; l++){
		int ii = lookupBond[2*l + I_INDEX];
		int jj = lookupBond[2*l + J_INDEX];
		m = 0;
		for (alpha = 0; alpha < DIM; alpha++){
			for (beta = 0; beta < DIM; beta++){
				thisHes = -second[l]*rij[DIM*l+alpha]*rij[DIM*l+beta]
				          -first[l]*DELTA(alpha,beta);
				a = DIM*ii + alpha;
				b = DIM*jj + beta;
				fprintf(fH, "%d\t%d\t%.15g\n", a, b, thisHes);
				fprintf(fH, "%d\t%d\t%.15g\n", b, a, thisHes);
				diagonal[DIM*DIM*ii + m] += second[l]*rij[DIM*l+alpha]*rij[DIM*l+beta] + first[l]*DELTA(alpha,beta);
				diagonal[DIM*DIM*jj + m] += second[l]*rij[DIM*l+alpha]*rij[DIM*l+beta] + first[l]*DELTA(alpha,beta);
				m++;
			}
		}
	}
	for (i = 0; i < N; i++){
		m = 0;
		for (alpha = 0; alpha < DIM; alpha++){
			for (beta = 0; beta < DIM; beta++){
				a = DIM*i + alpha; b = DIM*i + beta;
				fprintf(fH, "%d\t%d\t%.15g\n", a, b, diagonal[DIM*DIM*i + m]);
				m++;
			}
		}
	}
	free(diagonal);
	fclose(fH);

	/* vna_eta = H^-1 * xi_eta via existing cgSolver */
	vna_eta = (double *)malloc(sizeof(double)*DIM*N);
	compression_nonaffine_velocities(vna_eta);
	for (i = 0; i < DIM*N; i++) fprintf(fV, "%.15g\n", vna_eta[i]);
	free(vna_eta);
	fclose(fV);

	/* random b and w = H^-1 b — for Python identity check */
	{
		char nameb[64], namew[64];
		double *bvec = (double *)malloc(sizeof(double)*DIM*N);
		double *wvec = (double *)calloc(DIM*N, sizeof(double));
		FILE *fb, *fw;
		double sx=0, sy=0, sz=0;
		sprintf(nameb, "check_b_g%d.txt", g);
		sprintf(namew, "check_Hinvb_g%d.txt", g);
		srand(12345);
		for (i=0; i<DIM*N; i++) bvec[i] = (double)rand()/RAND_MAX - 0.5;
		for (i=0; i<N; i++){ sx+=bvec[DIM*i]; sy+=bvec[DIM*i+1]; sz+=bvec[DIM*i+2]; }
		sx/=N; sy/=N; sz/=N;
		for (i=0; i<N; i++){ bvec[DIM*i]-=sx; bvec[DIM*i+1]-=sy; bvec[DIM*i+2]-=sz; }
		cgSolver(wvec, bvec);
		fb = fopen(nameb, "w"); for (i=0;i<DIM*N;i++) fprintf(fb,"%.15g\n",bvec[i]); fclose(fb);
		fw = fopen(namew, "w"); for (i=0;i<DIM*N;i++) fprintf(fw,"%.15g\n",wvec[i]); fclose(fw);
		free(bvec); free(wvec);
	}
	printf("  wrote hessian_g%d.txt  vna_eta_g%d.txt  check_b_g%d.txt  check_Hinvb_g%d.txt\n",g,g,g,g);
}

/* alpha computed in Python — this stub kept so main compiles unchanged */
double get_alpha_theory(double K){
	return 0.0;
}

double get_first_order_dilatancy(){
	double res = 0.0;
	double *xi_gamma, *vna_eta;
	
	vna_eta = (double *)malloc(sizeof(double)*DIM*N);
	xi_gamma = (double *)malloc(sizeof(double)*DIM*N);
	
	u_gamma_x(xi_gamma);
	compression_nonaffine_velocities(vna_eta);
	
	res += u_gamma_eta();
	
	res += cdot(vna_eta,xi_gamma);
	
	free(vna_eta); free(xi_gamma);
	return -res/(V*(double)DIM);
}




/* Main: loop over all 10 LJ glass inherent states (s00000 to s00009).
   Reads  3dlj_N4000_sXXXXX.dat  from the current directory.
   Writes elasticity_results.txt with columns:
      serial  pressure  u_per_N  G  K  first_order_dilatancy  second_order_dilatancy
   No command-line arguments required.  Run as:  ./elasticity
*/
int main(int n, char **inputStrings){
	int t0, g;
	double bulk_modulus, shear_modulus;
	double k_born, k_na, g_born, g_na;
	double first_order, second_order;
	char snapFileName[1024];
	FILE *outFile;

	allocate_memory();
	srand(time(0));
	t0 = (int)time(0);

	outFile = fopen("elasticity_results.txt", "w");
	fprintf(outFile, "# serial   pressure    u_per_N      G           K           first_order  second_order\n");

	for (g = 0; g < 10; g++){
		sprintf(snapFileName, "3dlj_N%d_s%.5d.dat", N, g);
		if ( readSnapShot(snapFileName) ){
			printf("read %s,  |grad|/|f| = %g\n", snapFileName,
			       typicalGrad / typicalInteractionStrength);

			printf("  writing H and vna_eta for Python...\n");
			write_for_python(g);

			get_shear_modulus(&g_born, &g_na);
			shear_modulus = g_born - g_na;

			get_bulk_modulus(&k_born, &k_na);
			bulk_modulus = (k_born - k_na) / ((double)(DIM * DIM)) + pressure;

			first_order  = get_first_order_dilatancy();
			second_order = get_second_order_dilatancy();

			/* columns: serial  pressure  u/N  G  K  dil1  dil2 */
			fprintf(outFile, "%d\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\n",
			        g, pressure, u/DOUBLE_N,
			        shear_modulus, bulk_modulus,
			        first_order, second_order);
			fflush(outFile);
		}
		else
			printf("WARNING: %s not found, skipping.\n", snapFileName);
	}

	fclose(outFile);
	printf("Done. Took %d seconds.\n", (int)time(0) - t0);
	free_everything();
	return 0;
}
