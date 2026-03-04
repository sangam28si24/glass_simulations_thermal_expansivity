#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define	N							16000
#define	DOUBLE_N						((double)N)
#define	DIM							3
#define	CGSOLVER_TOL				1e-11
#define	DR							0.1

/* stuff that we don't touch usually */
#define	DELTA(a,b)					((a)==(b)?1.0:0.0)
#define	MAX_NEBZ						128
#define	CONTACTS					32
#define	MAX_CELLS					(N/3)
#define	X_COMP					0
#define 	Y_COMP 					1
#define 	Z_COMP 					2
#define	I_INDEX					0
#define	J_INDEX					1
#define	CUTOFF_SQRD					(1.48*1.48)
#define	C_0							(-1.1106337662511798)
#define	C_2							1.2676152372297065
#define	C_4							(-0.4960406072849212)
#define	C_6							0.0660511826415732
#define	N_0							10.0


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
double u;			/*	potential energy 		*/
double stress;
double pressure;
double typicalGrad;
double typicalContactForce;

int contacts;

double *rx, *ry, *rz, *fx, *fy, *fz;
int *type; 		/*	polydispersity */

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
double *fourth;
double *third;
double *second;
double *first;
double *rij;
/*******************************************************************************************************/
	
const double invSizeSqrd[2][2] = {{1.00,0.718184429761563},{0.718184429761563,0.510204081632653}}; // [0][0] = small-small; [0][1] = [1][0] = large-small; [1][1] = large-large;

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

void initializeSystem(){
	cellCutOff = 1.4*sqrt(CUTOFF_SQRD) + DR;
	listCutOffSqrd = cellCutOff*cellCutOff;
	updateNebzLists();
	calculateEverything();
	return;
}

void saveSnapShot(char *snapFileName){
	int i;
	FILE *snapFile;
	snapFile = fopen(snapFileName,"wb");
	fprintf(snapFile,"%.15g\t%.15g\t%.10g\t%.10g\n",L,0.0,stress,u/DOUBLE_N);
	for (i=0; i<N; i++)
		fprintf(snapFile,"%.15g\t%.15g\t%.15g\t%d\n",rx[i],ry[i],rz[i],type[i]);
	fclose(snapFile);
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


void writeVector(double *vec, char *outFileName){
	FILE *outFile;
	int i;
	outFile = fopen(outFileName,"wb");
	for (i=0;i<N;i++)
		fprintf(outFile,"%.8g\t%.8g\t%.8g\t%.15g\t%.15g\t%.15g\t%d\n",rx[i],ry[i],rz[i],
												vec[DIM*i+X_COMP],
												vec[DIM*i+Y_COMP],
												vec[DIM*i+Z_COMP],type[i]);
	fclose(outFile);
	return;
}


//this is only useful is a snapshot has already been loaded...
int readDirection(double *z, char *fileName){
	int i, int_dummy,goo;
	double dummy;
	FILE *inFile;
	inFile = fopen(fileName, "rb");
	if ( inFile ){
		for (i=0;i<N;i++)
			goo = fscanf(inFile, "%lf %lf %lf %lf %lf %lf %lf",&dummy,&dummy,&dummy,&(z[DIM*i + X_COMP]),&(z[DIM*i + Y_COMP]),&(z[DIM*i + Z_COMP]),&dummy);
		fclose(inFile);
		return 1;
	}
	return 0;
}

void updateNebzLists(){
	int i,j,x,y,z,current,m,m2,m3;
	int a,b,c,k,l,numHere,numThere,target,w,xIndex,yIndex;
	double dx,dy,dz,r2,invCellSize,listCutOff;
	
	nebListCounter = 0;
	maxD = 0.0;
	
	m =(int)(L/cellCutOff);
	
	//printf("\n\nm is %d\n\n",m); getchar();
	
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
	double sigmaOverR12,sigmaOverR14,sigmaOverR16,sigmaOverR18;
	double g,dx,dy,dz,r2OverSigma2,invSigmaSqrd;
	double temp,phi_r_over_r;
	
	u = 0.0; stress = 0.0; pressure = 0.0; contacts = 0; typicalContactForce = 0.0;
	for (i=0;i<N;i++){
		fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0; 
	}
	
	for (i=0;i<N-1;i++){
		m = numOfNebz[i];
		for (k=0;k<m;k++){
			j = nebz[i][k];
			if (j > i){
				dx = rx[j] - rx[i];
				dy = ry[j] - ry[i];
				dz = rz[j] - rz[i];
				// mess due to periodic boundary conditions 
				if ( dy >= HALF_LENGTH )
					dy -= LENGTH;
				else if ( dy < -HALF_LENGTH )
					dy += LENGTH;
				if ( dx >= HALF_LENGTH )
					dx -= LENGTH;
				else if ( dx < -HALF_LENGTH )
					dx += LENGTH;
				if ( dz >= HALF_LENGTH )
					dz -= LENGTH;
				else if ( dz < -HALF_LENGTH )
					dz += LENGTH;
				// end of mess 
				invSigmaSqrd = invSizeSqrd[type[i]][type[j]];
				r2OverSigma2 = L*L*( dx*dx + dy*dy + dz*dz )*invSigmaSqrd;
				if (r2OverSigma2 < CUTOFF_SQRD){
					lookupBond[2*contacts + I_INDEX] = i;
					lookupBond[2*contacts + J_INDEX] = j;
					rij[DIM*contacts + X_COMP] = L*dx;
					rij[DIM*contacts + Y_COMP] = L*dy;
					rij[DIM*contacts + Z_COMP] = L*dz;
					
					sigmaOverR18 = 1.0/(r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2*r2OverSigma2);
					sigmaOverR16 = sigmaOverR18*r2OverSigma2;
					sigmaOverR14 = sigmaOverR16*r2OverSigma2;
					sigmaOverR12 = sigmaOverR14*r2OverSigma2;
					
					phi_r_over_r = ( 2.0*C_2 + r2OverSigma2*(4.0*C_4 + r2OverSigma2*6.0*C_6) - N_0*sigmaOverR12 )*invSigmaSqrd; //this is -g in my previous language
					
					fourth[contacts] = N_0*(N_0+2.0)*(N_0+4.0)*(N_0+6.0)*sigmaOverR18*(invSigmaSqrd*invSigmaSqrd*invSigmaSqrd*invSigmaSqrd);
					third[contacts] = (48.0*C_6 - N_0*(N_0+2.0)*(N_0+4.0)*sigmaOverR16 )*(invSigmaSqrd*invSigmaSqrd*invSigmaSqrd);
					second[contacts] = ( 8.0*( C_4 + r2OverSigma2*3.0*C_6) + N_0*(N_0+2.0)*sigmaOverR14 )*(invSigmaSqrd*invSigmaSqrd);
					first[contacts] = phi_r_over_r;
					
					g = -phi_r_over_r;
					temp = L*dx*g;
					fx[j] += temp;
					fx[i] -= temp;
					temp = L*dy*g;
					fy[j] += temp;
					fy[i] -= temp;
					u += ( r2OverSigma2*(sigmaOverR12 + C_2 + r2OverSigma2*(C_4 + r2OverSigma2*C_6)) + C_0 );
					stress += -temp*L*dx;
					temp = L*dz*g;
					fz[j] += temp;
					fz[i] -= temp;
					temp = g*L*L*( dx*dx + dy*dy + dz*dz );
					pressure += temp;
					typicalContactForce += g*temp; //this is f^2
					contacts++;
				}
			}
		}
	}
	for (typicalGrad = 0.0,i=0;i<N;i++)
		typicalGrad += fx[i]*fx[i] + fy[i]*fy[i] + fz[i]*fz[i];
	typicalGrad = sqrt(typicalGrad/DOUBLE_N);
	typicalContactForce = sqrt(typicalContactForce/DOUBLE_N);
	stress = stress/V;
	pressure = pressure/( (double)DIM*V );
	return;
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
	
	for (i=0;i<N;i++)
		for (q=0;q<DIM;q++)
			res[DIM*i + q] = 0.0;
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		for (q=0;q<DIM;q++){
			res[DIM*i + q] -= second[l]*rij[DIM*l + alpha]*rij[DIM*l + beta]*rij[DIM*l + q];
			res[DIM*j + q] += second[l]*rij[DIM*l + alpha]*rij[DIM*l + beta]*rij[DIM*l + q];
		}
	}
	return;
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

void messianOnThreeVectors(double *v, double *u, double *w, double *res){
	int q,i,j,l;
	double vij[DIM],uij[DIM],wij[DIM];
	double v_dot_u, v_dot_w, u_dot_w, v_dot_r, u_dot_r, w_dot_r;
	
	for (i=0;i<N;i++)
		for (q=0;q<DIM;q++)
			res[DIM*i+q] = 0.0;
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		for (q=0;q<DIM;q++){
			vij[q] = v[DIM*j + q] - v[DIM*i + q];
			uij[q] = u[DIM*j + q] - u[DIM*i + q];
			wij[q] = w[DIM*j + q] - w[DIM*i + q];
		}
		
		v_dot_r = 0.0; u_dot_r = 0.0; w_dot_r = 0.0;
		v_dot_u = 0.0; v_dot_w = 0.0; u_dot_w = 0.0;
		for (q=0;q<DIM;q++){
			u_dot_r += rij[DIM*l + q]*uij[q];
			v_dot_r += rij[DIM*l + q]*vij[q];
			w_dot_r += rij[DIM*l + q]*wij[q];
			
			v_dot_u += vij[q]*uij[q];
			v_dot_w += vij[q]*wij[q];
			u_dot_w += uij[q]*wij[q];
		}
		
		for (q=0;q<DIM;q++){
			res[DIM*i + q] -= fourth[l]*v_dot_r*u_dot_r*w_dot_r*rij[DIM*l+q];
			res[DIM*j + q] += fourth[l]*v_dot_r*u_dot_r*w_dot_r*rij[DIM*l+q];
			
			res[DIM*i + q] -= third[l]*(v_dot_r*u_dot_r*wij[q] + v_dot_r*w_dot_r*uij[q] + v_dot_r*u_dot_w*rij[DIM*l+q] + u_dot_r*w_dot_r*vij[q] + u_dot_r*v_dot_w*rij[DIM*l+q] + w_dot_r*v_dot_u*rij[DIM*l+q]);
			res[DIM*j + q] += third[l]*(v_dot_r*u_dot_r*wij[q] + v_dot_r*w_dot_r*uij[q] + v_dot_r*u_dot_w*rij[DIM*l+q] + u_dot_r*w_dot_r*vij[q] + u_dot_r*v_dot_w*rij[DIM*l+q] + w_dot_r*v_dot_u*rij[DIM*l+q]);
			
			res[DIM*i + q] -= second[l]*(v_dot_u*wij[q] + v_dot_w*uij[q] + u_dot_w*vij[q]);
			res[DIM*j + q] += second[l]*(v_dot_u*wij[q] + v_dot_w*uij[q] + u_dot_w*vij[q]);
		}
	}
	return;
}

double messianOnFourVectors(double *v, double *u, double *w, double *y){
	int q,i,j,l=0;
	double vij[DIM],uij[DIM],wij[DIM],yij[DIM];
	double v_dot_u, v_dot_w, v_dot_y, u_dot_w, u_dot_y, w_dot_y;
	double v_dot_r,u_dot_r,w_dot_r,y_dot_r;
	double res = 0.0;
	
	for (l=0;l<contacts;l++){
		i = lookupBond[2*l + I_INDEX];
		j = lookupBond[2*l + J_INDEX];
		for (q=0;q<DIM;q++){
			vij[q] = v[DIM*j + q] - v[DIM*i + q];
			uij[q] = u[DIM*j + q] - u[DIM*i + q];
			wij[q] = w[DIM*j + q] - w[DIM*i + q];
			yij[q] = y[DIM*j + q] - y[DIM*i + q];
		}
		
		v_dot_r = 0.0; u_dot_r = 0.0; w_dot_r = 0.0; y_dot_r = 0.0;
		v_dot_u = 0.0; v_dot_w = 0.0; v_dot_y = 0.0; u_dot_w = 0.0; u_dot_y = 0.0; w_dot_y = 0.0;
		for (q=0;q<DIM;q++){
			u_dot_r += rij[DIM*l + q]*uij[q];
			v_dot_r += rij[DIM*l + q]*vij[q];
			w_dot_r += rij[DIM*l + q]*wij[q];
			y_dot_r += rij[DIM*l + q]*yij[q];
			
			v_dot_u += vij[q]*uij[q];
			v_dot_w += vij[q]*wij[q];
			v_dot_y += vij[q]*yij[q];
			u_dot_w += uij[q]*wij[q];
			u_dot_y += uij[q]*yij[q];
			w_dot_y += wij[q]*yij[q];
		}
		res += fourth[l]*v_dot_r*u_dot_r*w_dot_r*y_dot_r;
		res += third[l]*(v_dot_r*u_dot_r*w_dot_y + v_dot_r*w_dot_r*u_dot_y + v_dot_r*y_dot_r*u_dot_w + u_dot_r*w_dot_r*v_dot_y + u_dot_r*y_dot_r*v_dot_w + w_dot_r*y_dot_r*v_dot_u);
		res += second[l]*(v_dot_u*w_dot_y + v_dot_w*u_dot_y + v_dot_y*u_dot_w);
	}
	return res;
}


/***********************************************************************************************************************/
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
		/*hessianOnVector(x,pi);
		for (i=0;i<DIM*N;i++)
			printf("%g,%g\n",r[i],pi[i] - b[i]);
		getchar();
		scale1 = sqrt(scale1/( (double)(DIM*N) ) );
		scale2 = sqrt(scale2/( (double)(DIM*N) ) );*/
		
		
		
		//if ( !(k%1024) )
		//	printf("maxNorm:%g, scale:%g, ratio:%g\n",maxNorm,scale, maxNorm/scale);
		//	getchar();
		//}
		//~ fprintf(stderr, "\r",maxNorm);
		//~ fprintf(stderr, "%.10g\n",maxNorm);
		
		//if ( maxNorm < CGSOLVER_TOL )
		if ( maxNorm/scale < CGSOLVER_TOL )
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
	//printf("maxNorm:%g, scale:%g, ratio:%g\n",maxNorm,scale, maxNorm/scale);
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


void allocate_memory(){
	int i;
	
	rx = (double *)malloc(sizeof(double)*N);
	ry = (double *)malloc(sizeof(double)*N);
	rz = (double *)malloc(sizeof(double)*N);
	fx = (double *)malloc(sizeof(double)*N);
	fy = (double *)malloc(sizeof(double)*N);
	fz = (double *)malloc(sizeof(double)*N);
	type = (int *)malloc(sizeof(double)*N);
	
	numOfNebz = (int *)malloc(sizeof(int)*N);
	nebz = (int **)malloc(sizeof(int *)*N);
	for (i=0;i<N;i++)
		nebz[i] = (int *)malloc(sizeof(int)*MAX_NEBZ);
	
	numInCell = (int *)malloc(sizeof(int)*MAX_CELLS);
	cells = (int **)malloc(sizeof(int *)*MAX_CELLS);
	for (i=0;i<MAX_CELLS;i++)
		cells[i] = (int *)malloc(sizeof(int *)*MAX_NEBZ);
	
	lookupBond = (int *)malloc(sizeof(int)*N*CONTACTS*2);
	fourth = (double *)malloc(sizeof(double)*N*CONTACTS);
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
		
	free(lookupBond); free(fourth); free(third); 
	free(second); free(first); free(rij); 
	
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
		//~ printf("got mu_born = %.10g,  mu_na = %.10g, ratio %.10g and mu = %.10g\n",*mu_born,*mu_na,-mu_na/mu_born,mu_born+mu_na);
		//printf("got mu = %.10g\n",*mu_born-*mu_na);
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
	u_epsilon_x(X_COMP, Y_COMP, xi);
	
	if ( cgSolver(res,xi) ){
		for (a=0;a<DIM*N;a++)
			res[a] = -res[a];
		free(xi);
		return 1;
	}
	free(xi);
	return 0;
}

int main(int n,char **inputStrings){
	int q;
	double bulk_modulus, shear_modulus;
	double k_born, k_na, g_born, g_na;
	char snapFileName[1024];
	FILE *outFile;
	srand(time(0));
	
	allocate_memory();
	
	/* open output file */
	outFile = fopen("elasticity_results.txt","w");
	fprintf(outFile,"# serial\tpressure\tu_per_N\tG\tK\n");
	printf(         "# serial\tpressure\tu_per_N\tG\tK\n");
	
	/* loop over all 10 glass configurations */
	for(q=0; q<=9; q++){
		sprintf(snapFileName,"3dipl_glass_N%d_0.56_%.5d.dat", N, q);
		if ( readSnapShot(snapFileName) ){
			calculateShearModulus(&g_born, &g_na);
			shear_modulus = g_born - g_na;
			
			calculateBulkModulus(&k_born, &k_na);
			bulk_modulus = (k_born - k_na)/( (double)(DIM*DIM) ) + pressure;
			
			/* print to screen and to file: serial  pressure  u/N  G  K */
			printf(  "%d\t%.10g\t%.10g\t%.10g\t%.10g\n", q, pressure, u/DOUBLE_N, shear_modulus, bulk_modulus);
			fprintf(outFile, "%d\t%.10g\t%.10g\t%.10g\t%.10g\n", q, pressure, u/DOUBLE_N, shear_modulus, bulk_modulus);
		}
		else
			printf("could not read %s\n", snapFileName);
	}
	
	fclose(outFile);
	free_everything();
	return 0;
}

