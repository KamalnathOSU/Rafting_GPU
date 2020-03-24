#include <stdio.h>
#include <stdlib.h>
#include "shortcut_wrappers.h"
#include "pfr_cufft1.h"
double fr[L][M][N];
double gr[L][M][N];
cufftDoubleComplex fk[L][M][N/2+1];
int main()
{
int x_i,x_j,x_k;
double sum,tol=1e-6;
cufft_start();

fxijk
fr[x_i][x_j][x_k]=drand48();
efxijk
fr[0][0][0]=1;		

cufftrc3(&fr[0][0][0],&fk[0][0][0]);
cufftcr3(&fk[0][0][0],&gr[0][0][0]);

sum=0;
fxijk
sum+=(fr[x_i][x_j][x_k]-gr[x_i][x_j][x_k])*(fr[x_i][x_j][x_k]-gr[x_i][x_j][x_k]);
efxijk

printf("sum=%lf,tol=%lf\n",sum,tol);
if(sum<tol)
	printf("TEST PASSED\n");
else
	printf("TEST FAILED\n");
cufft_finish();

return 0;
}
