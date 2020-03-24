#include "global_variables.h"

double temp3_r[L*M*N];//used in the fftw functions
fftw_complex temp3_k[L*M*(N/2+1)];//used in the fftw functions
fftw_plan plan_forward;
fftw_plan plan_backward;
char fftw_thread_on=0;
void fftrc3(double *tpr, fftw_complex *tpk)
{
  int i;
  //Copying the input array to a seperate place
#pragma omp parallel for private(i)
  for(i=0;i<L*M*N;i++) 
  	temp3_r[i]=tpr[i];
  fftw_execute(plan_forward);
  //Copying the output array to the desired location
#pragma omp parallel for private(i)
  for(i=0;i<L*M*(N/2+1);i++) 
  {tpk[i][0]=temp3_k[i][0]*dxyz;tpk[i][1]=temp3_k[i][1]*dxyz;}
}//end of fftrc3

void fftcr3(fftw_complex *tpk, double *tpr)
{
  int i;
  //Copying the input array to a seperate place
#pragma omp parallel for private(i)
  for(i=0;i<L*M*(N/2+1);i++) 
  {temp3_k[i][0]=tpk[i][0];temp3_k[i][1]=tpk[i][1];}
  fftw_execute ( plan_backward );
#pragma omp parallel for private(i)
  for(i=0;i<LMN;i++)
  	tpr[i]=temp3_r[i]/(dxyz_V);
}//end of fftcr3
void fftw_finish()
{
fftw_destroy_plan(plan_forward);
fftw_destroy_plan(plan_backward);
}
void fftw_start(int nthreads, char option)
{
	if(fftw_init_threads())
	{
	fftw_thread_on=1;
	fftw_plan_with_nthreads(nthreads);
switch(option)
{
	case 3:
	plan_backward = fftw_plan_dft_c2r_3d( L, M, N, temp3_k, temp3_r, FFTW_PATIENT);
  	plan_forward  = fftw_plan_dft_r2c_3d( L, M, N, temp3_r, temp3_k, FFTW_PATIENT);
	break;
	case 2:
	plan_backward = fftw_plan_dft_c2r_3d( L, M, N, temp3_k, temp3_r, FFTW_MEASURE);
  	plan_forward  = fftw_plan_dft_r2c_3d( L, M, N, temp3_r, temp3_k, FFTW_MEASURE);
	break;
	case 1:
	plan_backward = fftw_plan_dft_c2r_3d( L, M, N, temp3_k, temp3_r, FFTW_ESTIMATE);
  	plan_forward  = fftw_plan_dft_r2c_3d( L, M, N, temp3_r, temp3_k, FFTW_ESTIMATE);
	break;
	default:
	fprintf(stderr,"#Error:Not an valid option for fftw\n");
}

  	if(!plan_backward) 
		{fprintf(stderr,"#Error:fftw plan backward could not be created\n");
		exit(0);}
  	if(!plan_forward) 
		{fprintf(stderr,"#Error:fftw plan forward could not be created\n");
		exit(0);}
  	}
	else
	{fprintf(stderr,"#Error:fftw can't be used with threads\n");
		exit(0);
	}
}

