/*
 * spring_networks_finite_T_dynamics.c
 *
 * Runs Berendsen NVT on the spring network at a set of temperatures.
 * Writes one output file per temperature, exactly matching the
 * 3dlj_NVT_NH.c column format so the notebook can reuse read_nvt().
 *
 * OUTPUT per temperature:
 *   spring_nvt_T{T:.6f}_s{serial:04d}.dat
 *   Columns: time   KE/N   PE/N   P_vir
 *   P_vir = pressure / (3*V)   [virial only, no kinetic term]
 *   — same as column 3 in LJ NVT output
 *
 * ═══════════════════════════════════════════════════════════
 * SOURCING
 * ═══════════════════════════════════════════════════════════
 * normal(mean, variance)
 *   SOURCED verbatim: 3dlj_NVT_NH.c Marsaglia polar RNG
 *
 * initVelocities(T)
 *   SOURCED: draw from normal(0,T) + fixDrift() pattern
 *   from 3dlj_NVT_NH.c main() lines 728-733
 *
 * u=0; pressure=0 in calculateForces()
 *   SOURCED: 3dlj_elasticity_analysis.c calculateEverything()
 *   change-log item 13 in 3dlj_NVT_NH.c
 *
 * run_berendsen_write() loop structure
 *   SOURCED: run_berendsen() in spring_networks_finite_T_dynamics.c
 *   + run_NH() output pattern in 3dlj_NVT_NH.c
 * ═══════════════════════════════════════════════════════════
 */

#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N               4000
#define DOUBLE_N        ((double)N)
#define DIM             3
#define X_COMP          0
#define Y_COMP          1
#define Z_COMP          2
#define I_INDEX         0
#define J_INDEX         1
#define TIME_STEP       0.1
#define HALF_TIME_STEP  (0.5*TIME_STEP)
#define TAU_B           1.0
#define LENGTH          1.0
#define HALF_LENGTH     0.5
#define TWO_PI          6.28318530717958
#define UNI             ((double)rand()/((double)RAND_MAX + 1.0))

int serial;
double L, invL, V, u, pressure=0.0, kinetic, T;
int contacts;

double *rx, *ry, *rz;
double *px, *py, *pz;
double *fx, *fy, *fz;

int    *lookupBond;
double *second, *first, *rij, *inv_r, *restlength;


/* ─── ORIGINAL functions — unchanged ──────────────────────────── */

void fixDrift(){
    int i; double Px=0,Py=0,Pz=0;
    for(i=0;i<N;i++){Px+=px[i];Py+=py[i];Pz+=pz[i];}
    Px/=DOUBLE_N; Py/=DOUBLE_N; Pz/=DOUBLE_N;
    for(i=0;i<N;i++){px[i]-=Px; py[i]-=Py; pz[i]-=Pz;}
}

int readNetwork(char *inFileName){
    int i,j,l,foo; double dummy,dx,dy,dz,r2; FILE *inFile;
    rx=(double*)malloc(sizeof(double)*N); ry=(double*)malloc(sizeof(double)*N);
    rz=(double*)malloc(sizeof(double)*N); px=(double*)malloc(sizeof(double)*N);
    py=(double*)malloc(sizeof(double)*N); pz=(double*)malloc(sizeof(double)*N);
    fx=(double*)malloc(sizeof(double)*N); fy=(double*)malloc(sizeof(double)*N);
    fz=(double*)malloc(sizeof(double)*N);
    inFile=fopen(inFileName,"rb");
    if(!inFile){printf("couldn't read %s\n",inFileName);return 0;}
    foo=fscanf(inFile,"%lf %lf %lf",&L,&dummy,&dummy);
    invL=1.0/L; V=L*L*L;
    for(i=0;i<N;i++) foo=fscanf(inFile,"%lf %lf %lf",&rx[i],&ry[i],&rz[i]);
    foo=fscanf(inFile,"%d %lf %lf",&contacts,&dummy,&dummy);
    lookupBond=(int*)   malloc(sizeof(int)   *2*contacts);
    rij       =(double*)malloc(sizeof(double)*contacts*DIM);
    inv_r     =(double*)malloc(sizeof(double)*contacts);
    first     =(double*)malloc(sizeof(double)*contacts);
    second    =(double*)malloc(sizeof(double)*contacts);
    restlength=(double*)malloc(sizeof(double)*contacts);
    for(l=0;l<contacts;l++){
        foo=fscanf(inFile,"%d %d %lf",&lookupBond[2*l+I_INDEX],&lookupBond[2*l+J_INDEX],&dummy);
        i=lookupBond[2*l+I_INDEX]; j=lookupBond[2*l+J_INDEX];
        dx=rx[j]-rx[i]; dy=ry[j]-ry[i]; dz=rz[j]-rz[i];
        if(dz>=.5)dz-=1;else if(dz<-.5)dz+=1;
        if(dy>=.5)dy-=1;else if(dy<-.5)dy+=1;
        if(dx>=.5)dx-=1;else if(dx<-.5)dx+=1;
        r2=L*L*(dx*dx+dy*dy+dz*dz);
        rij[DIM*l+X_COMP]=L*dx; rij[DIM*l+Y_COMP]=L*dy; rij[DIM*l+Z_COMP]=L*dz;
        inv_r[l]=1.0/sqrt(r2); first[l]=0.0;
        second[l]=1.0/r2; restlength[l]=sqrt(r2);
    }
    fclose(inFile); (void)foo; return 1;
}

/* ─── FIX — sourced from 3dlj_elasticity_analysis.c ───────────
 * u and pressure reset to zero at top of each call.          */
void calculateForces(){
    int l,i,j; double dx,dy,dz,r,g,tmp;
    u=0.0; pressure=0.0;
    for(i=0;i<N;i++){fx[i]=0.0; fy[i]=0.0; fz[i]=0.0;}
    for(l=0;l<contacts;l++){
        i=lookupBond[2*l+I_INDEX]; j=lookupBond[2*l+J_INDEX];
        dx=rx[j]-rx[i]; dy=ry[j]-ry[i]; dz=rz[j]-rz[i];
        if(dz>=.5)dz-=1;else if(dz<-.5)dz+=1;
        if(dy>=.5)dy-=1;else if(dy<-.5)dy+=1;
        if(dx>=.5)dx-=1;else if(dx<-.5)dx+=1;
        r=L*sqrt(dx*dx+dy*dy+dz*dz); g=(restlength[l]-r)/r;
        tmp=L*dx*g; fx[j]+=tmp; fx[i]-=tmp;
        tmp=L*dy*g; fy[j]+=tmp; fy[i]-=tmp;
        tmp=L*dz*g; fz[j]+=tmp; fz[i]-=tmp;
        u        += 0.5*(restlength[l]-r)*(restlength[l]-r);
        pressure += g*r*r;  /* W = sum of f_l * r_l over bonds */
    }
}

void advanceTimeBerendsen(){
    int i; double tmp,dx,dy,dz,inst_T,xi_T;
    for(i=0;i<N;i++){
        px[i]+=HALF_TIME_STEP*fx[i]; py[i]+=HALF_TIME_STEP*fy[i]; pz[i]+=HALF_TIME_STEP*fz[i];
        dx=TIME_STEP*px[i]; tmp=rx[i]+dx*invL;
        if(tmp>=1)rx[i]=tmp-1;else if(tmp<0)rx[i]=tmp+1;else rx[i]=tmp;
        dy=TIME_STEP*py[i]; tmp=ry[i]+dy*invL;
        if(tmp>=1)ry[i]=tmp-1;else if(tmp<0)ry[i]=tmp+1;else ry[i]=tmp;
        dz=TIME_STEP*pz[i]; tmp=rz[i]+dz*invL;
        if(tmp>=1)rz[i]=tmp-1;else if(tmp<0)rz[i]=tmp+1;else rz[i]=tmp;
    }
    calculateForces(); kinetic=0.0;
    for(i=0;i<N;i++){
        px[i]+=HALF_TIME_STEP*fx[i]; py[i]+=HALF_TIME_STEP*fy[i]; pz[i]+=HALF_TIME_STEP*fz[i];
        kinetic+=px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i];
    }
    kinetic*=0.5;
    inst_T=2.0*kinetic/((double)DIM*DOUBLE_N);
    xi_T=sqrt(1.0+TIME_STEP*(T/inst_T-1.0)/TAU_B);
    for(i=0;i<N;i++){px[i]*=xi_T; py[i]*=xi_T; pz[i]*=xi_T;}
}

void free_everything(){
    free(rx);free(ry);free(rz);free(px);free(py);free(pz);
    free(fx);free(fy);free(fz);free(lookupBond);free(second);
    free(first);free(rij);free(inv_r);free(restlength);
}


/* ─── ADD — sourced from 3dlj_NVT_NH.c ──────────────────────── */

/* normal(): Marsaglia polar RNG — copied verbatim from 3dlj_NVT_NH.c */
double normal(double mean, double variance){
    double v1,v2,rsq;
    do{v1=2.0*UNI-1.0; v2=2.0*UNI-1.0; rsq=v1*v1+v2*v2;}
    while(rsq>=1.0||rsq==0.0);
    return sqrt(variance)*(v1*sqrt(-2.0*log(rsq)/rsq))+mean;
}

/* initVelocities(): pattern from 3dlj_NVT_NH.c main() lines 728-733
 * Draw each component from N(0,T) then remove COM drift. */
void initVelocities(double T_target){
    int i;
    for(i=0;i<N;i++){
        px[i]=normal(0.0,T_target);
        py[i]=normal(0.0,T_target);
        pz[i]=normal(0.0,T_target);
    }
    fixDrift();
}

/* run_berendsen_write(): writes time series to file one row per PRINT_EVERY steps.
 *
 * Loop structure sourced from run_berendsen() in original spring file
 * and from run_NH() in 3dlj_NVT_NH.c.
 *
 * Column format: time  KE/N  PE/N  P_vir
 * P_vir = pressure/(3V) = W/(3V)  [virial only, no kinetic term]
 * Matches column 3 of 3dlj_NVT_NH.c output so notebook read_nvt() works directly.
 *
 * equil_steps: discarded (no file output) — equilibration
 * meas_steps:  written to file — production run
 */
void run_berendsen_write(double equil_duration, double meas_duration,
                         int print_every, FILE *outFile){
    int k;
    int equil_steps = (int)(equil_duration / TIME_STEP);
    int meas_steps  = (int)(meas_duration  / TIME_STEP);

    calculateForces();  /* initialise fx/fy/fz before first Verlet step */

    for(k=0;k<equil_steps;k++){
        advanceTimeBerendsen();
        if(!(k%16384)) fixDrift();
    }

    for(k=0;k<meas_steps;k++){
        advanceTimeBerendsen();
        if(!(k%16384)) fixDrift();

        if(!(k%print_every)){
            double t     = (double)k * TIME_STEP;
            double ke_n  = kinetic / DOUBLE_N;
            double pe_n  = u       / DOUBLE_N;
            double p_vir = pressure / (3.0*V);  /* W/(3V): virial pressure */
            fprintf(outFile, "%.8g\t%.8g\t%.8g\t%.8g\n", t, ke_n, pe_n, p_vir);
        }
    }
}


int main(int n, char **inputStrings){
    char networkFileName[1024];
    char outFileName[1024];
    double coordination;
    FILE *outFile;
    int t;

    /* Temperature array — same low-T range as LJ NVT notebook
     * to stay in the harmonic regime for meaningful theory comparison */
    double T_array[] = {0.000001, 0.000005, 0.000010, 0.000015, 0.000020, 0.000025, 0.000030, 0.000035, 0.000040};
    int    num_T     = 9;

    /* Print every 128 steps during production.
     * At TIME_STEP=0.1 and meas=2000 time units: 20000 steps / 128 = 156 rows per file. */
    int print_every = 128;

    srand(time(0));
    sscanf(inputStrings[1],"%d",&serial);

    coordination = 16.24;
    sprintf(networkFileName,"network_N%d_z%.2f_s%.4d.dat",N,coordination,serial);
    printf("network: %s\n",networkFileName);

    if(!readNetwork(networkFileName)){
        printf("not found: %s\n",networkFileName);
        return 1;
    }

    for(t=0;t<num_T;t++){
        T = T_array[t];

        sprintf(outFileName,"spring_nvt_T%.6f_s%.4d.dat",T,serial);

        printf("T=%.4g  ->  %s\n",T,outFileName); fflush(stdout);

        initVelocities(T);

        outFile = fopen(outFileName,"w");
        fprintf(outFile,"# spring network NVT  N=%d  z=%.2f  serial=%d  T=%.6f\n",
                N,coordination,serial,T);
        fprintf(outFile,"# Columns: time   KE/N   PE/N   P_vir\n");
        fprintf(outFile,"# P_vir = W/(3V): virial pressure (no kinetic term)\n");

        run_berendsen_write(500.0, 2000.0, print_every, outFile);

        fclose(outFile);
    }

    free_everything();
    return 0;
}
