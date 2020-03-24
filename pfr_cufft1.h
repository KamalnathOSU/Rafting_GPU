#ifndef PFR_CUFFT
#define PFR_CUFFT
#include <cufft.h>
#include <complex.h>
#include <openacc.h>
#include "constants.h"

void cufft_error(cufftResult ierr);
void cufft_start();
void cufftrc3(double *restrict tpr, cufftDoubleComplex *restrict tpk);
void cufftcr3(cufftDoubleComplex *restrict tpk,double *restrict tpr );
void cufftrc3k(double *restrict tpr, cufftDoubleComplex *restrict tpk);
void cufftcr3k(cufftDoubleComplex *restrict tpk,double *restrict tpr );
void cufft_finish();

#endif
