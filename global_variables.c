#include "global_variables.h"
//Concentration field
double ec_r[L][M][N];
double *ecr=&ec_r[0][0][0];


//Order-parameter field
double phi_r[V][L][M][N];
double *phir=&phi_r[0][0][0][0];

#ifdef DISLOCATION_ACTIVITY_ON
//Dislocation density field
double eta_r[P][L][M][N];
double *etar=&eta_r[0][0][0][0];
#endif

//Field to store the value of \vec{g} and \mod{g}
double g_v[L][M][N/2+1][3];
double g_mod2[L][M][N/2+1];

#ifdef DISLOCATION_ACTIVITY_ON
//Bpq matrix
double B_pq[L][M][N/2+1][P+1][P+1];
double *Bpq=&B_pq[0][0][0][0][0];
#else
double Bhh[L][M][N/2+1];
#endif

//h fields
double hphi_r[L][M][N];
double *hphir=&hphi_r[0][0][0];

//To store the transformation matrix A
double A[3][3];

#ifdef DISLOCATION_ACTIVITY_ON
//Slip system parameters
double n[][3]={{1,-1,1},{-1,-1,1},{-1,1,1},{1,1,1}
		,{1,1,1},{-1,1,1},{-1,-1,1},{1,-1,1}};// later *1/sqrt(3)
double b[][3]={{0,1,1},{1,0,1},{0,-1,1},{-1,0,1}
		,{0,-1,1},{1,0,1},{0,1,1},{-1,0,1}};// later *1/sqrt(2)
double te_p[P+1][3][3]; // Transformation strain
double ts_p[P+1][3][3];// Transformation stress
#else
double te_p[3][3];
double ts_p[3][3];
#endif
double s_app[3][3];
int LMN=L*M*N;
int LMNH=L*M*(N/2+1);
double temp1_r[L][M][N];
double temp2_r[L][M][N];
double temp3_r[L*M*N];//used in the fftw functions
  double dxyz=dx*dy*dz;
  double dxyz_V=dx*dy*dz*L*M*N;
// P dislocation + 1 misfit

//For performing fftw
cufftDoubleComplex ec_k[L][M][N/2+1];
cufftDoubleComplex *eck=&ec_k[0][0][0];
cufftDoubleComplex phi_k[V][L][M][N/2+1];
cufftDoubleComplex *phik=&phi_k[0][0][0][0];

#ifdef DISLOCATION_ACTIVITY_ON
cufftDoubleComplex eta_k[P][L][M][N/2+1];
cufftDoubleComplex *etak=&eta_k[0][0][0][0];
#endif
cufftDoubleComplex hphi_k[L][M][N/2+1];
cufftDoubleComplex *hphik=&hphi_k[0][0][0];
cufftDoubleComplex temp1_k[L][M][N/2+1];
cufftDoubleComplex temp2_k[L][M][N/2+1];

