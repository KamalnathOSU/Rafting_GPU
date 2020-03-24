#include "global_variables.h"
void initial_condition1()
{// Gamma prime particle in gamma
int x_i,x_j,x_k,xijk;
double ah=N/1.26/2.0;
double C=(N-1)/2.0;
	for_xijk
	if( fabs(x_i-C)<ah && fabs(x_j-C)<ah && fabs(x_k-C)<ah )
	{
	#ifdef DISLOCATION_ACTIVITU_ON
	etar[xijk]=0;
	#endif
	ecr[xijk]=cgp_e;
	phir[xijk]=1;
	}
	else
	{
	#ifdef DISLOCATION_ACTIVITY_ON
	etar[xijk]=0.01;
	#endif
	ecr[xijk]=cg_e;
	phir[xijk]=0;
	}
	efor_xijk
}//end of initial condition1

void initial_condition2()
{// Gamma prime , gamma interface
int x_i,x_j,x_k,xijk;
double C=(N-1)/2.0;
	for_xijk
	if( x_i<C )
	{
	#ifdef DISLOCATION_ACTIVITY_ON
	etar[xijk]=0;
	#endif
	ecr[xijk]=cgp_e;
	phir[xijk]=1;
	}
	else
	{
	#ifdef DISLOCATION_ACTIVITY_ON
	etar[xijk]=0.0;
	#endif
	ecr[xijk]=cg_e;
	phir[xijk]=0;
	}
	efor_xijk
}//end of initial_condition2

void initial_condition10()
{// 51x51x51 particle in 64x64x64 matrix
int x_i,x_j,x_k,xijk;
	for_xijk
	if( (x_i<=50) && (x_j<=50) && (x_k<=50) )
	{
	#ifdef DISLOCATION_ACTIVITY_ON
	etar[xijk]=0;
	#endif
	ecr[xijk]=cgp_e;
	phir[xijk]=1;
	}
	else
	{
	#ifdef DISLOCATION_ACTIVITY_ON
	etar[xijk]=0.0;
	#endif
	ecr[xijk]=cg_e;
	phir[xijk]=0;
	}
	efor_xijk
}//end of initial condition10

void initial_condition11(double chomo)
{//homogenous supersaturated gamma phase 
int x_i,x_j,x_k,xijk;
int i;
	for_xijk
	#ifdef DISLOCATION_ACTIVITY_ON
	etar[xijk]=0.00001;
	#endif
	ecr[xijk]=chomo;
	for(i=0;i<V;i++)
	phir[i*L*M*N+xijk]=0;
	efor_xijk
}//end of initial condition11

