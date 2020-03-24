//function declaration
#include "global_variables.h"
void print_rfield( double *tpr, char *filename);
void print_kfield( cufftDoubleComplex *tpk, char *filename);
void cal_vfrac();

void cal_vfrac()
{
double vfrac,cavg;
int x_i,x_j,x_k,xijk;
//Calculating volume fraction
vfrac=0;cavg=0;
for_xijk
	vfrac+=hphir[xijk];
	cavg+=ecr[xijk];
efor_xijk
vfrac/=L*M*N;
cavg/=L*M*N;
printf("Final volume fraction       = %lf\n Avg conc = %lf\n",vfrac,cavg);
}

