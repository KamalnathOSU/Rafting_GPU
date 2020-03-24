#include<stdio.h>
#include<stdlib.h>
#include <math.h>
//#include <fftw3.h>
#include "constants.h"
#include "shortcut_wrappers.h"
#include "pfr_cufft1.h"
#include "curand_cuda.h"

//Concentration field
extern double ec_r[L][M][N];
extern double *ecr;

//Order-parameter field
extern double phi_r[V][L][M][N];
extern double *phir;

#ifdef DISLOCATION_ACTIVITY_ON
//Dislocation density field
extern double eta_r[P][L][M][N];
extern double *etar;
#endif

//Field to store the value of \vec{g} and \mod{g}
extern double g_v[L][M][N/2+1][3];
extern double g_mod2[L][M][N/2+1];

#ifdef DISLOCATION_ACTIVITY_ON
//Bpq matrix
extern double B_pq[L][M][N/2+1][P+1][P+1];
extern double *Bpq;
#else
extern double Bhh[L][M][N/2+1];
#endif

//h fields
extern double hphi_r[L][M][N];
extern double *hphir;

//To store the transformation matrix A
extern double A[3][3];

#ifdef DISLOCATION_ACTIVITY_ON
//Slip system parameters
extern double n[8][3];// later *1/sqrt(3)
extern double b[8][3];// later *1/sqrt(2)
extern double te_p[P+1][3][3]; // Transformation strain
extern double ts_p[P+1][3][3];// Transformation stress
#else
extern double te_p[3][3];
extern double ts_p[3][3];
#endif
extern double s_app[3][3];
extern int LMN;
extern int LMNH;
extern double temp1_r[L][M][N];
extern double temp2_r[L][M][N];
  extern double dxyz;
  extern double dxyz_V;
// P dislocation + 1 misfit

//For performing fftw
extern cufftDoubleComplex ec_k[L][M][N/2+1];
extern cufftDoubleComplex *eck;
extern cufftDoubleComplex phi_k[V][L][M][N/2+1];
extern cufftDoubleComplex *phik;

#ifdef DISLOCATION_ACTIVITY_ON
extern cufftDoubleComplex eta_k[P][L][M][N/2+1];
extern cufftDoubleComplex *etak;
#endif
extern cufftDoubleComplex hphi_k[L][M][N/2+1];
extern cufftDoubleComplex *hphik;
extern cufftDoubleComplex temp1_k[L][M][N/2+1];
extern cufftDoubleComplex temp2_k[L][M][N/2+1];



