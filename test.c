#include <stdio.h>
#include <stdlib.h>
#define N (100)
double fr[N][N];
double *f=&fr[0][0];
void mulscale();
int main()
{
 int i;
#pragma acc parallel loop copy(f[0:N*N-1])
 for(i=0;i<N*N;i++)
		 f[i]=2.0;

#pragma acc data copy(f[0:N*N-1])
 {
mulscale();
 }
 for(i=0;i<5;i++)
		 printf("%lf \n",f[i]);
 return 0;
}
void mulscale()
{
int i;
double mu_factor[N*N],mu;
for(i=0;i<N*N-1;i++)
		mu_factor[i]=drand48();
#pragma acc parallel loop private (i,mu)
for(i=0;i<N*N-1;i++)
{
		mu=mu_factor[i];
		f[i]*=mu;
}

#pragma acc parallel loop private (i,mu)
for(i=0;i<N*N-1;i++)
{
		mu=mu_factor[i];
		f[i]/=mu;
}

}
