#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <omp.h>
#include "global_variables.h"
#include "constants.h"
#include "shortcut_wrappers.h"

void grad(); // To calculate the k^2 table
void StrainBpq();
void OutputMed(int timestep);
void GetEintC( int timestep);
void init_field();

int main(int argc, char *argv[])
{
int total_step=2000;
int timestep;
time_t t1,t2;
printf("Starting the program ...\n");
t1=time(NULL);
cufft_start();
curand_init();

//To store output data
if(!opendir("output_files"))  mkdir("output_files",0777);

#ifndef DISLOCATION_ACTIVITY_ON
printf("Dislocation activity is OFF\n");
#else
printf("Dislocation activity is ON\n");
#endif
total_step=TOTAL_TIMESTEP;
if(argc>1)
	total_step=atoi(argv[1]);
printf("total_step = %d\n",total_step);
printf("Initializing K^2 table ...\n");
grad();//Initializing k^2 table
printf("K^2 table initialized \n");
printf("Initializing Bpq matrix ...\n");
StrainBpq();//Initializing Bpq
printf("Bpq matrix initialized \n");
printf("Initializing the field ...\n");
init_field();
printf("All fields are initialized \n");
OutputMed(0);
t2=time(NULL);
printf("Cpu time for initialization = %lf\n",difftime(t2,t1));
printf("Starting the time evolution...\n");
	printf("time = %9d\n",0);
t1=time(NULL);

#pragma acc data copy(ec_r,hphi_r,phi_r,g_v,g_mod2)
#pragma acc data copyin(ec_k,hphi_k,phi_k,dxyz,dxyz_V)
#pragma acc data create(temp1_r,temp1_k,temp2_r,temp2_k)
#ifdef DISLOCATION_ACTIVITY_ON
	#pragma acc data copyin(B_pq,eta_r,eta_k,n)
#else
	#pragma acc data copyin(Bhh)
#endif
{// Data region for openacc
for(timestep=1;timestep<=total_step;timestep++)
{
	GetEintC(timestep);
	if( (timestep % SAVE_CONFIG )==0 )
	{
	#pragma acc update host(ec_r,hphi_r,phi_r)
#ifdef DISLOCATION_ACTIVITY_ON
	#pragma acc update host(eta_r)
#endif
			OutputMed(timestep);
	}
}
}//end of Data region
printf("\n");
t2=time(NULL);
OutputMed(total_step);
cufft_finish();
curand_finish();
printf("Total cpu time utilized for time evolution = %lf sec\n",difftime(t2,t1));

printf("#All is well \n");
return 0;
}//end of main()

//a function to print the relavant fields at the given time
void usage()
{

}







	
