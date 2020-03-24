#include "global_variables.h"
void grad()
{
  int g_i,g_j,g_k,gijk;
#pragma omp parallel for private(g_i,g_j,g_k,gijk)
for_gijk
	if(g_i<L/2)	g_v[g_i][g_j][g_k][0]=2*M_PI/(L*dx) *g_i;
	else			g_v[g_i][g_j][g_k][0]=2*M_PI/(L*dx) *(g_i-L);
	if(g_j<M/2)	g_v[g_i][g_j][g_k][1]=2*M_PI/(M*dy) *g_j;
	else			g_v[g_i][g_j][g_k][1]=2*M_PI/(M*dy) *(g_j-M);
	g_v[g_i][g_j][g_k][2]=2*M_PI/(N*dz) *g_k;
	g_mod2[g_i][g_j][g_k]=g_v[g_i][g_j][g_k][0]* \
			      g_v[g_i][g_j][g_k][0]  \
			     +g_v[g_i][g_j][g_k][1]* \
			      g_v[g_i][g_j][g_k][1] \
			     +g_v[g_i][g_j][g_k][2]* \
			      g_v[g_i][g_j][g_k][2];
efor_gijk
}//end of grad()

