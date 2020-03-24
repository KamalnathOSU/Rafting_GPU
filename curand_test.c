#include <stdio.h>
#include <openacc.h>
#include "curand_cuda.h"
#define N (1000)
float randnum[N];
int main()
{
		int i;
		curand_init();
#pragma acc data copyout(randnum)
		{
		curand_generate((float*)acc_deviceptr(randnum),N);
		}
		curand_finish();
	for(i=0;i<10;i++)
			printf("%f\n",randnum[i]);
return 0;
}
