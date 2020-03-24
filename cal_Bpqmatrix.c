#include "global_variables.h"
//Using non-symmetric Bpq
void compute_inv(double m[3][3], double minv[3][3]);
void StrainBpq()
{
  int g_i,g_j,g_k,gijk;//Used in the shortcut wrapper
  int i,j,k,l;
  int ii_i,ii_j,ii_k,ii_l;
  double omega_inv[3][3],omega[3][3],Bpqsum=0;
#ifdef DISLOCATION_ACTIVITY_ON
  double Bpq_const[P+1][P+1];
#else
	//double Bhh_const;
	//No need to calculate Bhh_const as anyway we are subtracting the average Bpq value in the end
#endif
  double C11,C12,C44;

#ifdef DISLOCATION_ACTIVITY_ON
  double *tep=&te_p[0][0][0];
  double *tsp=&ts_p[0][0][0];
#else
	double *tep=&te_p[0][0];
	double *tsp=&ts_p[0][0];
#endif
  double *g=&g_v[0][0][0][0];
  double *g2=&g_mod2[0][0][0];

#ifdef ISOTROPIC_ELASTICITY
  // Initialize Cijkl
C11=2*G*(1-v_poisson)/(1-2*v_poisson);
C12=2*G*v_poisson/(1-2*v_poisson);
C44=G;
#else
C11=Cijkla;
C12=Cijklb;
C44=Cijklc;
#endif

  //calcuating "A" marix to transform crystal co-ordinates to global co-ordinates
  //By default "A" is identity
  A[0][0]=A[1][1]=A[2][2]=1;

  //Calculate the applied stress in global coordinate
  s_app[2][2]=s_app_mag;

#ifdef DISLOCATION_ACTIVITY_ON
  //Including the factor sqrt(3),sqrt(2)
  //Normalizing the slip system parameters
for(i=0;i<P;i++)
	for(j=0;j<3;j++)
	{
		n[i][j]=n[i][j]/sqrt(3);
		b[i][j]=b[i][j]/sqrt(2);
	}
  //Compute the transformation strain for dislocation
  for(k=0;k<P;k++)
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	  te_p[k][i][j]=(n[k][i]*b[k][j]+n[k][j]*b[k][i])/2;
   //Compute the transformation strain for misfit
  te_p[P][0][0]=MISFIT_STRAIN;
  te_p[P][1][1]=MISFIT_STRAIN;
  te_p[P][2][2]=MISFIT_STRAIN;
 //Compute transformation stress for all the defects (P dislocations + 1 misfit)
  for(k=0;k<=P;k++)
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{if(i==j)
			ts_p[k][i][j]=C11*te_p[k][i][j]+C12*(te_p[k][(i+1)%3][(j+1)%3]+te_p[k][(i+2)%3][(j+2)%3]);
		 else
		 	ts_p[k][i][j]=2*C44*te_p[k][i][j];
		}
#else
 //Compute the transformation strain for misfit
  te_p[0][0]=MISFIT_STRAIN;
  te_p[1][1]=MISFIT_STRAIN;
  te_p[2][2]=MISFIT_STRAIN;
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{if(i==j)
			ts_p[i][j]=C11*te_p[i][j]+C12*(te_p[(i+1)%3][(j+1)%3]+te_p[(i+2)%3][(j+2)%3]);
		 else
		 	ts_p[i][j]=2*C44*te_p[i][j];
		}

#endif

#ifdef DISLOCATION_ACTIVITY_ON
//Calculate the constant part of Bpq
for(i=0;i<=P;i++)
   for(j=0;j<=P;j++)
   {Bpq_const[i][j]=0;
     for(l=0;l<9;l++)
       Bpq_const[i][j]+=tsp[i*9+l]*tep[j*9+l];
   }
#endif
 //Bpq calculation
#pragma omp parallel for private(g_i,g_j,g_k,gijk,i,j,ii_i,ii_j,ii_k,ii_l,omega_inv,omega)
 for_gijk
	 if(gijk==0)	continue;//skipping the origin of reciprocal space

    for(i=0;i<3;i++)
	for(j=0;j<3;j++)
	{if(i==j) omega_inv[i][j] = C44 + (C11-C44)*g[gijk*3+i]*g[gijk*3+i]/g2[gijk];
	else	  omega_inv[i][j] = (C12+C44)*g[gijk*3+i]*g[gijk*3+j]/g2[gijk];
	}

//Computing the inverse
    compute_inv(omega_inv,omega);
#ifdef DISLOCATION_ACTIVITY_ON
 //Calculating the Bpq values
	for(i=0;i<=P;i++)  //looping over the defect
	  for(j=0;j<=P;j++)//looping over the defect
	  {
	    B_pq[g_i][g_j][g_k][i][j]=Bpq_const[i][j];
			  for(ii_i=0;ii_i<3;ii_i++)
			    for(ii_j=0;ii_j<3;ii_j++)
			      for(ii_k=0;ii_k<3;ii_k++)
				for(ii_l=0;ii_l<3;ii_l++)
						B_pq[g_i][g_j][g_k][i][j]-=g[gijk*3+ii_i]*ts_p[i][ii_i][ii_j]* \
											omega[ii_j][ii_k]*                    \
											ts_p[j][ii_k][ii_l]*g[gijk*3+ii_l] \
											/g2[gijk];
	  }
#else
	Bhh[g_i][g_j][g_k]=3*(C11+2*C12)*MISFIT_STRAIN*MISFIT_STRAIN;
	for(ii_i=0;ii_i<3;ii_i++)
	for(ii_j=0;ii_j<3;ii_j++)
	for(ii_k=0;ii_k<3;ii_k++)
	for(ii_l=0;ii_l<3;ii_l++)
	Bhh[g_i][g_j][g_k]-=g[gijk*3+ii_i]*ts_p[ii_i][ii_j]*omega[ii_j][ii_k]* \
				ts_p[ii_k][ii_l]*g[gijk*3+ii_l] \
				/g2[gijk];
#endif
 efor_gijk
	
#ifdef DISLOCATION_ACTIVITY_ON
 //Subtracting the average Bpq for misfit defects
 for(i=0;i<=P;i++)
 {
   Bpqsum=0;
#pragma omp parallel for private(g_i,g_j,g_k,gijk) reduction(+: Bpqsum)
   for_gijk
// tuned to match Ning's code
     Bpqsum+=B_pq[g_i][g_j][g_k][i][P]*( 1.0 + (g_k!=0)*1.0 );
   efor_gijk
   Bpqsum/=(L*M*N);
#pragma omp parallel for private(g_i,g_j,g_k,gijk)
   for_gijk
	   if(gijk==0) continue;//skipping the origin
     B_pq[g_i][g_j][g_k][P][i]-=Bpqsum;
   efor_gijk
 }
#else
 Bpqsum=0;//slight modification is done to avoid overflow problem
#pragma omp parallel for private(g_i,g_j,g_k,gijk) reduction(+: Bpqsum)
   for_gijk
// tuned to match Ning's code
     Bpqsum+=Bhh[g_i][g_j][g_k]*( 1.0 + ( (g_k!=0)&&(g_k!=N/2) )*1.0 );
   efor_gijk
   Bpqsum/=(L*M*N);
#pragma omp parallel for private(g_i,g_j,g_k,gijk)
   for_gijk
	   if(gijk==0) continue;//skipping the origin
     Bhh[g_i][g_j][g_k]-=Bpqsum;
   efor_gijk
#endif
 //  Use the applied stress routine to calculate the stress in global coordinates

}//end of StrainBpq

void compute_inv(double m[3][3], double minv[3][3])
{
//This function computes inverse of 3*3 matrix
double det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) - \
             m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + \
	     m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
double invdet = 1 / det;
if(det>1e100)
	fprintf(stderr,"#Warning:Determinant is too huge. The matrix inversion may not be correct.\n");
  // inverse of matrix m
  minv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
  minv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
  minv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
  minv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
  minv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
  minv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
  minv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
  minv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
  minv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
}
