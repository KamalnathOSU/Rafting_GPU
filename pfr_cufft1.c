#include "pfr_cufft1.h"

//double temp3_r[L*M*N];//used in the cufft functions
//cufftDoubleComplex temp3_k[L*M*(N/2+1)];
double *temp3_r=NULL;
cufftDoubleComplex *temp3_k=NULL;
cufftHandle plan2, plan3;
char fftw_thread_on=0;
int ierr;

void cufft_start()
{
cufft_error(cufftPlan3d(&plan2, L,M,N , CUFFT_D2Z));
cufft_error(cufftPlan3d(&plan3, L,M,N , CUFFT_Z2D));
cufft_error(cufftSetStream(plan2,(cudaStream_t) acc_get_cuda_stream(acc_async_sync)));
cufft_error(cufftSetStream(plan3,(cudaStream_t) acc_get_cuda_stream(acc_async_sync)));
temp3_r=(double*)acc_malloc(sizeof(double)*L*M*N);
temp3_k=(cufftDoubleComplex*)acc_malloc(sizeof(cufftDoubleComplex)*L*M*(N/2+1));
}//end of cufft_start()

void cufftrc3(double *restrict tpr, cufftDoubleComplex *restrict tpk)
{
int i;
#pragma acc data copyin(tpr[0:L*M*N]) copyout(tpk[0:L*M*(N/2+1)])
{
#pragma acc host_data use_device(tpr,tpk)
cufft_error(cufftExecD2Z(plan2, tpr, tpk ));
}
}//end of cufftrc3()

void cufftcr3(cufftDoubleComplex *restrict tpk,double *restrict tpr )
{
int i;
#pragma acc data copyout(tpr[0:L*M*N]) copyin(tpk[0:L*M*(N/2+1)])
{
#pragma acc host_data use_device(tpr,tpk)
cufft_error(cufftExecZ2D(plan3, tpk, tpr ));
#pragma acc parallel loop
for(i=0;i<L*M*N;i++)
		tpr[i]=tpr[i]/(L*M*N);
}
}//end of cufftcr3()

void cufftrc3k(double *restrict tpr, cufftDoubleComplex *restrict tpk)
{
int i;
//#pragma acc data  create(temp3_r,temp3_k)
#pragma acc data  deviceptr(temp3_r,temp3_k)
{
#pragma acc data present(tpr)
#pragma acc parallel loop 
for(i=0;i<L*M*N;i++)
		temp3_r[i]=tpr[i];
//#pragma acc host_data use_device(temp3_r,temp3_k)
cufft_error(cufftExecD2Z(plan2, temp3_r, temp3_k ));
#pragma acc data present(tpk)
#pragma acc parallel loop 
for(i=0;i<L*M*(N/2+1);i++)
		{tpk[i].x=temp3_k[i].x;
		 tpk[i].y=temp3_k[i].y;
		}
}
}//end of cufftrc3k()

void cufftcr3k(cufftDoubleComplex *restrict tpk,double *restrict tpr )
{
int i;
//#pragma acc data create(temp3_r,temp3_k)
#pragma acc data  deviceptr(temp3_r,temp3_k)
{
#pragma acc data present(tpk)
#pragma acc parallel loop
for(i=0;i<L*M*(N/2+1);i++)
	{temp3_k[i].x=tpk[i].x;
	 temp3_k[i].y=tpk[i].y;
	}
//#pragma acc host_data use_device(temp3_r,temp3_k)
cufft_error(cufftExecZ2D(plan3, temp3_k, temp3_r ));
#pragma acc data present(tpr)
#pragma acc parallel loop
for(i=0;i<L*M*N;i++)
		tpr[i]=temp3_r[i]/(L*M*N);
}
}//end of cufftcr3k()

void cufft_finish()
{
cufft_error(cufftDestroy(plan2));
cufft_error(cufftDestroy(plan3));
acc_free(temp3_r);
acc_free(temp3_k);
}


void cufft_error(cufftResult ierr)
{
 if(ierr==CUFFT_SUCCESS)   ;//The CUFFT operation was successful
 else
 {       fprintf(2,"#Error occured at %s in line %d\n",__FILE__,__LINE__);
	 if(ierr==CUFFT_INVALID_PLAN                ) {printf("#Error:CUFFT was passed an invalid plan hand \n");exit(0);} 
	 else if(ierr==CUFFT_ALLOC_FAILED                ) {printf("#Error:CUFFT failed to allocate GPU or CPU memory \n");exit(0);}
	 else if(ierr==CUFFT_INVALID_TYPE                ) {printf("#Error:No longer used \n");exit(0);}
	 else if(ierr==CUFFT_INVALID_VALUE               ) {printf("#Error:User specified an invalid pointer or parameter \n");exit(0);}
	 else if(ierr==CUFFT_INTERNAL_ERROR              ) {printf("#Error:Driver or internal CUFFT library error \n");exit(0);}
	 else if(ierr==CUFFT_EXEC_FAILED                 ) {printf("#Error:Failed to execute an FFT on the GPU \n");exit(0);}
	 else if(ierr==CUFFT_SETUP_FAILED                ) {printf("#Error:The CUFFT library failed to initialize \n");exit(0);}
	 else if(ierr==CUFFT_INVALID_SIZE                ) {printf("#Error:User specified an invalid transform size \n");exit(0);}
	 else if(ierr==CUFFT_UNALIGNED_DATA              ) {printf("#Error: \n");exit(0);}
	 else if(ierr==CUFFT_INCOMPLETE_PARAMETER_LIST   ) {printf("#Error: \n");exit(0);} 
	 else if(ierr==CUFFT_INVALID_DEVICE              ) {printf("#Error:Plan creation and execution are on different device \n");exit(0);}
	 else if(ierr==CUFFT_PARSE_ERROR                 ) {printf("#Error: \n");exit(0);}
	 else if(ierr==CUFFT_NO_WORKSPACE                ) {printf("#Error: \n");exit(0);}
	 else if(ierr==CUFFT_NOT_IMPLEMENTED             ) {printf("#Error: \n");exit(0);}
	 else if(ierr==CUFFT_LICENSE_ERROR               ) {printf("#Error: \n");exit(0);}
	 else
	 {fprintf(2,"#Error:Unknown error code. please update the check_error_function\n");
	   exit(0);}
 }//end of else
}


