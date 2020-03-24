#include "curand_cuda.h"

curandGenerator_t gen;			

void curand_init()
{
/* Create pseudo-random number generator */ 
CURAND_CALL(curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT));
/* Set seed */
    CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen,1234ULL));
}
void curand_generate(float *devData,int n)
{
/* Generate n floats on device */ 
		CURAND_CALL(curandGenerateUniform(gen, devData, n));
}
void curand_finish()
{
CURAND_CALL(curandDestroyGenerator(gen));
}
