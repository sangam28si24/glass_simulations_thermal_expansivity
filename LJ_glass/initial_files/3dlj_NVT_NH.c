/* 3d LJ (Kob-Andersen-style binary mixture) NVT code with Nose-Hoover thermostat.
   Potential: phi(r) = 4[ (sigma/r)^12 - (sigma/r)^6 + C0 + C2*(r/sigma)^2
                                                       + C4*(r/sigma)^4 + C6*(r/sigma)^6 ]
   Cutoff: r/sigma = 2.5  (standard LJ cutoff with smooth polynomial correction).
   N=4000, rho=0.541.  Glass files: 3dlj_N4000_sXXXXX.dat */
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


//notice that here N = 4n^3 with n some integer. here n^3 = 1000 = 10^3 giving N=4000
#define	N							4000
#define	DOUBLE_N					((double)N)
#define	DIM							3
#define	RHO							0.541
#define	TAU_NH						1.0		//for NH thermostat
#define	TAU_B						3.0		//for berendsen thermostat
#define	NUM_OF_DOF					((double)(DIM*N - DIM - 1)) //momentum & energy conservation

/* LJ potential parameters *****************************************************/
#define	CUTOFF					2.5
#define	CUTOFF_SQRD				(CUTOFF * CUTOFF)
#define	C_0						(0.075421526038318)
#define	C_2						(-0.035416041040039)
#define	C_4						(0.005708815153688)
#define	C_6						(-0.000318029410118)
#define	N_REP					12.0		/* repulsive exponent   */
#define	N_ATT					6.0		/* attractive exponent  */
/*******************************************************************************/

/* constant constants... :P*******************************************************/
#define	LENGTH				1.0
#define	HALF_LENGTH			0.5

#define	TWO_PI					6.28318530717958
#define	UNI		((double)rand()/((double)RAND_MAX + 1.0))

/* time-step and structural defines *******************************************/
#define	TIME_STEP				0.001
#define	HALF_TIME_STEP				(0.5*TIME_STEP)
#define	DELTA(a,b)				((a)==(b)?1.0:0.0)
#define	MAX_NEBZ				128
#define	MAX_CELLS				(N/3)
#define	X_COMP					0
#define	Y_COMP					1
#define	Z_COMP					2
#define	DR					0.10


void initializeSystem();
void updateNebzLists();
void calculateForces();

int serial;

/*globals*/
double L; 			/*	length of <square> box 	*/
double invL;			/* inverse length of box	*/
double V;			/*	volume			*/
double T;			/*	temperature			*/
double t;				/*	time				*/
double u;			/*	potential energy 		*/
double kinetic;		/*	kinetic energy 		*/
double stress;			/*	x-y plane stress		*/
double pressure;
double typicalInteractionStrength;
double typicalGrad;

double *rx, *ry, *rz;		/*	(rescaled) particle coordinates rx/L, ry/L, rz/L */
double *px, *py, *pz;	/*	momentum (here = velocities because all m=1)	*/
double *fx, *fy, *fz;	/*	net forces */
int *type; 			/*	for binary system	*/

/******************************* for nebz list & cell subdivision *******************************/
int **nebz;
int *numOfNebz;
int *numInCell;
int **cells;
double cellCutOff,maxD,listCutOffSqrd; //for updating nebz list.
int nebListCounter; //for counting how many steps have took place for nebz list updating.
const int shift[13][3] = {{-1,-1,1},{-1,0,1},{-1,1,1},{0,-1,1},{0,0,1},{0,1,1},{1,-1,1},{1,0,1},{1,1,1},{-1,-1,0},{-1,0,0},{-1,1,0},{0,1,0}};
/************************************************************************************/

/*******************************  Nose-Hoover variables  *******************************/
double eta = 0.0;   /* Nose-Hoover friction parameter */
/***************************************************************************************************/


const double invSizeSqrd[2][2] = {{1.00,0.718184429761563},{0.718184429761563,0.510204081632653}}; // [0][0] = small-small; [0][1] = [1][0] = large-small; [1][1] = large-large;


/************************aux functions*****************************************************/

double normal(double mean, double variance){
	double v1,v2,rsq;
	do{
		v1=2.0*UNI-1.0;
		v2=2.0*UNI-1.0;
		rsq =v1*v1+v2*v2;
	} while (rsq >= 1.0 || rsq == 0.0);
	return sqrt(variance)*( v1*sqrt(-2.0*log(rsq)/rsq) ) + mean;
}

int getUpperThrdRoot(int number){
	int q=1;
	while (q*q*q<number)
		q++;
	return q;
}

/************************end of aux functions*****************************************************/


/************************input/output functions*****************************************************/

void saveSnapShot(char *snapFileName){
	int i;
	FILE *snapFile;
	snapFile = fopen(snapFileName,"wb");
	fprintf(snapFile,"%.15g\t%.15g\t%.10g\t%.10g\n",L,0.0,u/DOUBLE_N,pressure);
	for (i=0; i<N; i++)
		fprintf(snapFile,"%.15g\t%.15g\t%.15g\t%d\n",rx[i],ry[i],rz[i],type[i]);
	fclose(snapFile);
	return;
}

// returns 1 upon successful reading of fileName
int readSnapShot(char *snapFileName){
	int i,q;
	double dummy;
	FILE *snapFile;
	snapFile = fopen(snapFileName, "rb");
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

/************************end of input/output functions*****************************************************/

/************************ memory management *****************************************************/

void allocate_memory(){
	int i;
	
	rx = (double *)malloc(sizeof(double)*N);
	ry = (double *)malloc(sizeof(double)*N);
	rz = (double *)malloc(sizeof(double)*N);
	px = (double *)malloc(sizeof(double)*N);
	py = (double *)malloc(sizeof(double)*N);
	pz = (double *)malloc(sizeof(double)*N);
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
		
	return;
}

void free_everything(){
	int i;
	free(rx); free(ry); free(rz); 
	free(px); free(py); free(pz); 
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
	return;
}
/************************ end of memory management *****************************************************/



void initialize_fcc_grid(){
        int i, j, M, thrdrtM;
        double space, x0, y0, z0;
	
	//allocate particle types
        for (i=0; i<N; i++)
                type[i] = 0;

        for (i=0; i<N/2; i++){
                j = (int)(UNI*DOUBLE_N);
                while ( type[j] )
                        j = (int)(UNI*DOUBLE_N);
                type[j] = 1;
                px[i] = 0.0; py[i] = 0.0; pz[i] = 0.0;
        }
        
        V = DOUBLE_N/RHO;
        M = N/4;
        thrdrtM = getUpperThrdRoot(M);
	L = pow(V,1.0/3.0); invL = 1.0/L;
        space = 1.0/(double)thrdrtM;
        for (i=0; i<M; i++){
                x0 = space*(double)( (i%(thrdrtM*thrdrtM)) % thrdrtM );
                y0 = space*(double)((i%(thrdrtM*thrdrtM)) / thrdrtM );
                z0 = space*(double)(i/(thrdrtM*thrdrtM));

		rx[4*i] = x0; ry[4*i] = y0; rz[4*i] = z0;

                //3 more:
                rx[4*i+1] = x0 + 0.5*space;
                ry[4*i+1] = y0;
                rz[4*i+1] = z0 + 0.5*space;

                rx[4*i+2] = x0 + 0.5*space;
                ry[4*i+2] = y0 + 0.5*space;
                rz[4*i+2] = z0;

                rx[4*i+3] = x0;
                ry[4*i+3] = y0 + 0.5*space;
                rz[4*i+3] = z0 + 0.5*space;
        }
        initializeSystem();
        return;
}

void initializeSystem(){
	cellCutOff = 1.4*sqrt(CUTOFF_SQRD) + DR;
	listCutOffSqrd = cellCutOff*cellCutOff;
	updateNebzLists();
	calculateForces();
	return;
}

void fixDrift(){
	int i;
	double Px=0.0,Py=0.0,Pz=0.0;
	for (i=0;i<N;i++){ Px += px[i]; Py += py[i]; Pz += pz[i]; }
	Px = Px/DOUBLE_N; Py = Py/DOUBLE_N; Pz = Pz/DOUBLE_N; 
	for (i=0;i<N;i++){ px[i] -= Px; py[i] -= Py; pz[i] -= Pz; }
	return;
}

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

void calculateForces(){
	int i,j,m,k;
	double g,dx,dy,dz,r2OverSigma2,invSigmaSqrd;
	double sigma2OverR2,sigma6OverR6,sigma8OverR8,temp;

	//set stuff to zero
	u = 0.0; stress = 0.0; pressure = 0.0; typicalInteractionStrength=0.0;
	for (i=0;i<N;i++){ fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0; }

	for (i=0;i<N-1;i++){
		m = numOfNebz[i];
		for (k=0;k<m;k++){
			j = nebz[i][k];
			if (j > i){
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
				invSigmaSqrd = invSizeSqrd[type[i]][type[j]];
				r2OverSigma2 = L*L*( dx*dx + dy*dy + dz*dz )*invSigmaSqrd;
				if ( r2OverSigma2 < CUTOFF_SQRD ){
					/* LJ force: g = -(1/r)(dphi/dr), so F_alpha = g * r_alpha */
					sigma2OverR2 = 1.0/r2OverSigma2;
					sigma6OverR6 = sigma2OverR2*sigma2OverR2*sigma2OverR2;
					sigma8OverR8 = sigma6OverR6*sigma2OverR2;
					g = 4.0*( sigma8OverR8*(N_REP*sigma6OverR6 - N_ATT)
					          - 2.0*C_2
					          - r2OverSigma2*(4.0*C_4 + 6.0*r2OverSigma2*C_6)
					        )*invSigmaSqrd;

					temp = L*dx*g;
					fx[j] += temp;
					fx[i] -= temp;
					temp = L*dy*g;
					fy[j] += temp;
					fy[i] -= temp;
					stress += -temp*L*dx;
					temp = L*dz*g;
					fz[j] += temp;
					fz[i] -= temp;

					/* LJ energy */
					u += 4.0*( sigma6OverR6*(sigma6OverR6 - 1.0)
					           + r2OverSigma2*(C_2 + r2OverSigma2*(C_4 + r2OverSigma2*C_6))
					           + C_0 );

					pressure += g*L*L*( dx*dx + dy*dy + dz*dz );
					typicalInteractionStrength += g*g*L*L*( dx*dx + dy*dy + dz*dz );
				}
			}
		}
	}
	stress = stress/V;
	pressure = pressure/( (double)DIM*V );
	for (typicalGrad = 0.0,i=0;i<N;i++)
		typicalGrad += fx[i]*fx[i] + fy[i]*fy[i] + fz[i]*fz[i];
	typicalGrad = sqrt(typicalGrad/DOUBLE_N);
	typicalInteractionStrength = sqrt(typicalInteractionStrength/DOUBLE_N);
	return;
}

/***********************************************************************************************************************/

void advanceTimeNVE(){
	int i;
	double temp,dx,dy,dz,inst_T;
	
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
		
		temp = dx*dx + dy*dy + dz*dz;
		if (temp > maxD)
			maxD = temp;
	}
	
	nebListCounter++;
	if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR )
		updateNebzLists();
	
	calculateForces();
	for (kinetic = 0.0,i=0; i<N; i++){
		px[i] += HALF_TIME_STEP*fx[i];
		py[i] += HALF_TIME_STEP*fy[i];
		pz[i] += HALF_TIME_STEP*fz[i];
		kinetic += px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
	}
	kinetic = 0.5*kinetic;
	
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
		
		temp = dx*dx + dy*dy + dz*dz;
		if (temp > maxD)
			maxD = temp;
	}
	
	nebListCounter++;
	if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR )
		updateNebzLists();
	
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


void advanceTimeNoseHoover() {
	//following Toxvaerd 1991 Molecular Physics Vol. 72.
	int i;
	double px_midpoint, py_midpoint, pz_midpoint;
	double factor, temp, dx, dy, dz, sumOfSquaredMidpoints = 0.0;
	
	for (i=0; i<N; i++){
		// calculate v(t+0.5dt)
		px_midpoint = (1.0 - eta*HALF_TIME_STEP)*px[i] + HALF_TIME_STEP*fx[i];
		py_midpoint = (1.0 - eta*HALF_TIME_STEP)*py[i] + HALF_TIME_STEP*fy[i];
		pz_midpoint = (1.0 - eta*HALF_TIME_STEP)*pz[i] + HALF_TIME_STEP*fz[i];
		
		sumOfSquaredMidpoints += px_midpoint*px_midpoint + py_midpoint*py_midpoint + pz_midpoint*pz_midpoint;
		
		// calculate r(t+dt)
		dx = TIME_STEP*px_midpoint; //real space displacement
		temp = rx[i] + dx*invL;
		if (temp >= LENGTH)
			rx[i] = temp - LENGTH;
		else if (temp < 0.0)
			rx[i] = temp + LENGTH;
		else
			rx[i] = temp;

		dy = TIME_STEP*py_midpoint; //real space displacement
		temp = ry[i] + dy*invL;
		if (temp >= LENGTH)
			ry[i] = temp - LENGTH;
		else if (temp < 0.0)
			ry[i] = temp + LENGTH;
		else
			ry[i] = temp;
		
		dz = TIME_STEP*pz_midpoint; //real space displacement
		temp = rz[i] + dz*invL;
		if (temp >= LENGTH)
			rz[i] = temp - LENGTH;
		else if (temp < 0.0)
			rz[i] = temp + LENGTH;
		else
			rz[i] = temp;

		temp = dx*dx + dy*dy + dz*dz;
		if (temp > maxD)
			maxD = temp;
		
		//don't need old velocities anymore, will need midpoints
		px[i] = px_midpoint;
		py[i] = py_midpoint;
		pz[i] = pz_midpoint;
	}
	
	//evolve eta
	eta += TIME_STEP*( sumOfSquaredMidpoints - NUM_OF_DOF*T ) / ( NUM_OF_DOF*T*TAU_NH*TAU_NH );
	
	nebListCounter++;
	if ( 2.0*((double)nebListCounter)*sqrt(maxD) > DR ) {
		updateNebzLists();
	}
	calculateForces();
	
	//now advance velocities to t + \delta t
	kinetic = 0.0; factor = 1.0/( 1.0 + HALF_TIME_STEP*eta );
	for (i=0; i<N; i++){
		px[i] = factor*(px[i] + HALF_TIME_STEP*fx[i]);
		py[i] = factor*(py[i] + HALF_TIME_STEP*fy[i]);
		pz[i] = factor*(pz[i] + HALF_TIME_STEP*fz[i]);
		kinetic += px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
	}
	kinetic = 0.5*kinetic;
	return;
}

void run_NVE(double duration){
        int k,steps,between=256;
        //FILE *outFile;
	
        //outFile = fopen("conservation_test.dat","wb");
        steps = (int)(duration/TIME_STEP);
        
        for (k=0; k<steps; k++){
                advanceTimeNVE();
                if ( !(k%between) ){
                        //fprintf(outFile,"%.8g\t%.8g\t%.8g\t%.12g\n",(double)k*TIME_STEP, kinetic/DOUBLE_N, u/DOUBLE_N, (kinetic + u)/DOUBLE_N );
			printf("%.8g\t%.8g\t%.8g\t%.12g\n",(double)k*TIME_STEP, kinetic/DOUBLE_N, u/DOUBLE_N, (kinetic + u)/DOUBLE_N );
                }
                if ( !(k%16384) )
                        fixDrift();
        }
        //fclose(outFile);
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


void run_NH(double duration, char *outFileName){
        int k,steps,between=128;
        FILE *outFile;
	
        outFile = fopen(outFileName,"w");
        steps = (int)(duration/TIME_STEP);
        
        for (k=0; k<steps; k++){
                advanceTimeNoseHoover();
                if ( !(k%between) ){
			/* columns: time   KE/N   PE/N   pressure(virial, no kinetic term) */
                        fprintf(outFile,"%.8g\t%.8g\t%.8g\t%.8g\n",(double)k*TIME_STEP, kinetic/DOUBLE_N, u/DOUBLE_N, pressure);
                }
                if ( !(k%16384) )
                        fixDrift();
        }
        fclose(outFile);
        return;
}



int main(int argc, char **argv){
	int i, glass_serial;
	char snapFileName[256];
	char outFileName[256];

	/* ---- read inputs from command line ---- */
	/* Usage:  ./nvt  <temperature>  <glass_serial>       */
	/* Example: ./nvt 0.01 3                               */
	if(argc < 3){
		printf("Usage: ./nvt <temperature> <glass_serial_0_to_9>\n");
		return 1;
	}
	sscanf(argv[1], "%lf", &T);
	sscanf(argv[2], "%d",  &glass_serial);

	allocate_memory();
	srand(time(0) + glass_serial);   /* different random seed per glass */

	/* load the LJ glass inherent state */
	sprintf(snapFileName, "3dlj_N%d_s%.5d.dat", N, glass_serial);
	if( !readSnapShot(snapFileName) ){
		printf("ERROR: could not read %s\n", snapFileName);
		return 1;
	}
	printf("Loaded %s,  running at T=%.4f\n", snapFileName, T);

	/* initialise velocities from Maxwell-Boltzmann at temperature T */
	for(i=0; i<N; i++){
		px[i] = normal(0.0, T);
		py[i] = normal(0.0, T);
		pz[i] = normal(0.0, T);
	}
	fixDrift();

	/* name the output file */
	sprintf(outFileName, "nvt_lj_T%.4f_g%d.dat", T, glass_serial);

	/* run */
	run_berendsen(20.0);         /* fast thermalisation          */
	run_NH(500.0, outFileName);  /* production run, saves output */

	printf("Done. Output written to %s\n", outFileName);
	free_everything();
	return 0;
}


