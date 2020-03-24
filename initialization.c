#include "global_variables.h"
void init_field();
void init_from_file(double *tpr, char *filename);// this function is in io.c
void initial_condtion11(double);

void init_field()
{
int i;
int x_i,x_j,x_k,xijk;
double p;
#ifdef FROM_INPUT_FILE
char filename[500];
//Reading for phi
for(i=0;i<V;i++)
{
	sprintf(filename,"input_files/phi_%d_t18800.vtk",i+1);
	init_from_file(&phi_r[i][0][0][0],filename);
}
/*
	#ifdef DISLOCATION_ACTIVITY_ON
	//Reading for eta
	for(i=0;i<P;i++)
	{
		sprintf(filename,"input_files/eta_%d.vtk",i+1);
		init_from_file(&eta_r[i][0][0][0],filename);
	}
	#endif
*/
//Reading for ecr
	sprintf(filename,"input_files/ecr_t18800.vtk");
	init_from_file(&ec_r[0][0][0],filename);
//Initializing eta
	#ifdef DISLOCATION_ACTIVITY_ON
	for_xijk
	if(ecr[xijk]<0.2)
	{	for(i=0;i<P;i++)
			eta_r[i][x_i][x_j][x_k]=0.001*drand48();
	}
	else
	{	for(i=0;i<P;i++)
			eta_r[i][x_i][x_j][x_k]=0;
	}
	efor_xijk
	#endif

#else
//  initial_condition11(0.1945);
	initial_condition1();
#endif
//Initialzing h(phi)
for_xijk
	hphir[xijk]=0;
	for(i=0;i<V;i++)
		{p=phi_r[i][x_i][x_j][x_k];
		hphir[xijk]+=p*p*p*(10-15*p+6*p*p);
		}
efor_xijk

//Initializing stuffs in fourier space
cufftrc3(&ec_r[0][0][0],&ec_k[0][0][0]);
for(i=0;i<V;i++)
	cufftrc3(&phi_r[i][0][0][0],&phi_k[i][0][0][0]);
#ifdef DISLOCATION_ACTIVITY_ON
for(i=0;i<P;i++)
	cufftrc3(&eta_r[i][0][0][0],&eta_k[i][0][0][0]);
#endif
cufftrc3(&hphi_r[0][0][0],&hphi_k[0][0][0]);

}//end of init_field()



