#ifndef CURAND_CUDA
#define CURAND_CUDA
#include <curand.h>

#define CUDA_CALL(x) if((x)!=cudaSuccess) { printf("Error at %s:%d\n",__FILE__,__LINE__);} 

#define CURAND_CALL(x) if((x)!=CURAND_STATUS_SUCCESS) {  printf("Error at %s:%d\n",__FILE__,__LINE__);	}

void curand_init();
void curand_generate(float *devData,int n);
void curand_finish();

#endif
