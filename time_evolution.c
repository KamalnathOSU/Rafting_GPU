#include "global_variables.h"
void OutputMed(int timestep);
  int i,j,l;
  int g_i,g_j,g_k;
  int x_i,x_j,x_k;
  double temp1,temp2;
  double p=0;
  double g2,Wo,hr;
  cufftDoubleComplex nk,ck,hk;
  float randnum[L][M][N];

#ifdef DISLOCATION_ACTIVITY_ON
  double etasum[P]={0};
  double temp_etasum=0;
  double W=s_app_mag/SQRT_SIX;//This is the plastic work done on a single system
  char stress_flag=1;
//Function to compute plastic strain
inline double cal_eps_plastic()
{
double sum=0;
for(int i_var=0;i_var<P;i_var++)
		sum+=etasum[i_var]*te_p[i_var][2][2];
return sum/(L*M*N);
}
#endif



void GetEintC(int timestep)
{

#ifdef DISLOCATION_ACTIVITY_ON
	  Wo=Leta*dt*W*LMN*dxyz;
  //**********Evolving Eta **************
  for(i=0;i<P;i++)//loop over eta's
  {
#pragma omp parallel for private(g_i,g_j,g_k,gijk,temp1,temp2,l,j)
#pragma acc data present(g_mod2,hphi_k,g_v,n,B_pq,eta_k)
#pragma acc parallel loop collapse(3) private(g_i,g_j,g_k,temp1,temp2,g2,hk)
	  fgijk
	  //calculating the denominator in the evolution equation
	  g2=g_mod2[g_i][g_j][g_k];
	  hk=hphi_k[g_i][g_j][g_k];
	  temp1=g_v[g_i][g_j][g_k][0]*n[i][0]+g_v[g_i][g_j][g_k][1]*n[i][1]+g_v[g_i][g_j][g_k][2]*n[i][2];
	  temp1*=temp1;
	  temp1=1+Leta*dt*( B_pq[g_i][g_j][g_k][i][i]+keta*(g2-temp1) );
          {
		  //Real
		  temp2=0;
		  for(j=0;j<P;j++)
			  temp2+=B_pq[g_i][g_j][g_k][i][j]*eta_k[j][g_i][g_j][g_k].x;
		  temp2= temp2 -B_pq[g_i][g_j][g_k][i][i]*eta_k[i][g_i][g_j][g_k].x \
			 + B_pq[g_i][g_j][g_k][i][P]*hk.x;
		  eta_k[i][g_i][g_j][g_k].x= (eta_k[i][g_i][g_j][g_k].x - Leta*dt*temp2 +  ( (g_i==0)&&(g_j==0)&&(g_k==0) )*Wo  )/temp1;
		  //Imaginary
		  temp2=0;
		  for(j=0;j<P;j++)
			  temp2+=B_pq[g_i][g_j][g_k][i][j]*eta_k[j][g_i][g_j][g_k].y;
		  temp2= temp2 -B_pq[g_i][g_j][g_k][i][i]*eta_k[i][g_i][g_j][g_k].y \
			 + B_pq[g_i][g_j][g_k][i][P]*hk.y;
		  eta_k[i][g_i][g_j][g_k].y= (eta_k[i][g_i][g_j][g_k].y - Leta*dt*temp2)/temp1;
	  	  }
	  efgijk
	  //Converting eta back to real space
	  cufftcr3k(&eta_k[i][0][0][0],&eta_r[i][0][0][0]);

	  //Cutting off eta  
	  temp_etasum=0;
#pragma omp parallel for private(x_i,x_j,x_k,xijk)
#pragma acc data present(eta_r,ec_r)
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k) reduction(+:temp_etasum)
      fxijk
	  if( (eta_r[i][x_i][x_j][x_k]<0) | (ec_r[x_i][x_j][x_k]>(cg_e+cgp_e)/2.0) )	
	  		eta_r[i][x_i][x_j][x_k]=0;
	  temp_etasum+=eta_r[i][x_i][x_j][x_k];
	  efxijk
	  etasum[i]=temp_etasum;

          //Converting to k space
          cufftrc3k(&eta_r[i][0][0][0],&eta_k[i][0][0][0]);
	  
  }//end of loop over eta's
#endif
  //*******Evolving phi ***************
  // phi loop
  for(i=0;i<V;i++){
#pragma omp parallel for private(g_i,g_j,g_k,gijk,l,j,hk)
#pragma acc data present(hphi_k,temp1_k)
#ifdef DISLOCATION_ACTIVITY_ON
#pragma acc data present(eta_k,B_pq)
#else
#pragma acc data present(Bhh)
#endif
#pragma acc parallel loop collapse(3) private(g_i,g_j,g_k,j,hk)
	  fgijk
		{
	    hk=hphi_k[g_i][g_j][g_k];
		temp1_k[g_i][g_j][g_k].x=0;
#ifdef DISLOCATION_ACTIVITY_ON
		for(j=0;j<P;j++)
			temp1_k[g_i][g_j][g_k].x+=B_pq[g_i][g_j][g_k][P][j]*eta_k[j][g_i][g_j][g_k].x;
		temp1_k[g_i][g_j][g_k].x+=B_pq[g_i][g_j][g_k][P][P]*hk.x;
#else
		temp1_k[g_i][g_j][g_k].x+=Bhh[g_i][g_j][g_k]*hk.x;
#endif
		temp1_k[g_i][g_j][g_k].y=0;
#ifdef DISLOCATION_ACTIVITY_ON
		for(j=0;j<P;j++)
			temp1_k[g_i][g_j][g_k].y+=B_pq[g_i][g_j][g_k][P][j]*eta_k[j][g_i][g_j][g_k].y;
		temp1_k[g_i][g_j][g_k].y+=B_pq[g_i][g_j][g_k][P][P]*hk.y;
#else
		temp1_k[g_i][g_j][g_k].y+=Bhh[g_i][g_j][g_k]*hk.y;
#endif
		}
	  efgijk
	  cufftcr3k(&temp1_k[0][0][0],&temp1_r[0][0][0]);

#pragma omp parallel for private(x_i,x_j,x_k,xijk,p)
#pragma acc data present(phi_r,temp1_r,ec_r,hphi_r)
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k,p)
	  fxijk
	  p=phi_r[i][x_i][x_j][x_k];
	  temp1_r[x_i][x_j][x_k]= 30*p*p*(1-p)*(1-p)* \
		(temp1_r[x_i][x_j][x_k]+2*fo/Vm*(cg_e-cgp_e)*( ec_r[x_i][x_j][x_k] + hphi_r[x_i][x_j][x_k]*(cg_e-cgp_e) - cg_e )  );
#ifdef INCLUDE_MISFIT_WORK_TERM
	  temp1_r[x_i][x_j][x_k] += -30*phir[i*LMN+xijk]*phir[i*LMN+xijk]*(1-phir[i*LMN+xijk])*(1-phir[i*LMN+xijk])*s_app_mag*MISFIT_STRAIN; 
#endif
	  efxijk
	  cufftrc3k(&temp1_r[0][0][0],&temp1_k[0][0][0]);
	  
#pragma omp parallel for private(x_i,x_j,x_k,temp1,j,p)
#pragma acc data present(phi_r,temp2_r)
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k,p,temp1,j)
	  fxijk
	  p=phi_r[i][x_i][x_j][x_k];
		  temp1=0;
		  for(j=0;j<V;j++)	temp1+=phi_r[j][x_i][x_j][x_k]*phi_r[j][x_i][x_j][x_k];
		  temp1-=p*p;
	  temp2_r[x_i][x_j][x_k]=p*(2+p*(-6+4*p)) + 2*alpha_cross*p*temp1;
	  efxijk

	  cufftrc3k(&temp2_r[0][0][0],&temp2_k[0][0][0]);

#pragma omp parallel for private(g_i,g_j,g_k,l)
#pragma acc data present(g_mod2,temp1_k,phi_k,temp2_k)
#pragma acc parallel loop collapse(3) private(g_i,g_j,g_k,g2)
	  fgijk
		{
	  g2=g_mod2[g_i][g_j][g_k];
	  temp1_k[g_i][g_j][g_k].x=temp2_k[g_i][g_j][g_k].x*interface_omega+temp1_k[g_i][g_j][g_k].x;
	  phi_k[i][g_i][g_j][g_k].x = (phi_k[i][g_i][g_j][g_k].x- Lphi*dt*temp1_k[g_i][g_j][g_k].x)/(1+Lphi*kphi*g2*dt);
	  temp1_k[g_i][g_j][g_k].y=temp2_k[g_i][g_j][g_k].y*interface_omega+temp1_k[g_i][g_j][g_k].y;
	  phi_k[i][g_i][g_j][g_k].y = (phi_k[i][g_i][g_j][g_k].y- Lphi*dt*temp1_k[g_i][g_j][g_k].y)/(1+Lphi*kphi*g2*dt);
		  }
	  efgijk
         
	//Converting back to real space
	cufftcr3k(&phi_k[i][0][0][0],&phi_r[i][0][0][0]);
#pragma acc data present(phi_r)
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
	fxijk
		if(phi_r[i][x_i][x_j][x_k]>1.0)
			phi_r[i][x_i][x_j][x_k]=1.0;
		if(phi_r[i][x_i][x_j][x_k]<0.0)
			phi_r[i][x_i][x_j][x_k]=0;
	efxijk
        }// end of phi loop 

  //Compute h(r)
#pragma omp parallel for private(x_i,x_j,x_k,xijk,j)
#pragma acc data present(hphi_r,phi_r)
#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k,j)
  fxijk
  hphi_r[x_i][x_j][x_k]=0;
  for(j=0;j<V;j++)
       	  hphi_r[x_i][x_j][x_k]+=phi_r[j][x_i][x_j][x_k]*phi_r[j][x_i][x_j][x_k]*phi_r[j][x_i][x_j][x_k] \
       		    *(10-15*phi_r[j][x_i][x_j][x_k]+6*phi_r[j][x_i][x_j][x_k]*phi_r[j][x_i][x_j][x_k]); 
  efxijk
  //Compute \tilde{h}
	  cufftrc3k(&hphi_r[0][0][0],&hphi_k[0][0][0]);
  

// Evolution of concentration
temp1=2*fo*Vm*Mc*dt;
#pragma omp parallel for private(g_i,g_j,g_k,gijk,l,temp2)
#pragma acc data present(g_mod2,hphi_k,ec_k)
#pragma acc parallel loop collapse(3) private(g_i,g_j,g_k,g2,hk,temp1)
fgijk
{
	  g2=g_mod2[g_i][g_j][g_k];  
	  hk=hphi_k[g_i][g_j][g_k];
	ec_k[g_i][g_j][g_k].x=(ec_k[g_i][g_j][g_k].x+temp1*(cgp_e-cg_e)*g2*hk.x)/(1+temp1*g2);
	ec_k[g_i][g_j][g_k].y=(ec_k[g_i][g_j][g_k].y+temp1*(cgp_e-cg_e)*g2*hk.y)/(1+temp1*g2);
}
efgijk
cufftcr3k(&ec_k[0][0][0],&ec_r[0][0][0]);

#ifdef NOISE_ON
if(timestep < RAND_FLAG)
{
	#pragma acc data create(randnum) present(phi_r)
	for(i=0;i<V;i++)
	{
	    curand_generate((float*)acc_deviceptr(&randnum[0][0][0]),L*M*N);
		#pragma acc parallel loop collapse(3) private(x_i,x_j,x_k)
		fxijk
		phi_r[i][x_i][x_j][x_k]+=((float)drand48()-0.5)/20 ;
		if(phi_r[i][x_i][x_j][x_k]>1.0)
			phi_r[i][x_i][x_j][x_k]=1.0;
		if(phi_r[i][x_i][x_j][x_k]<0.0)
			phi_r[i][x_i][x_j][x_k]=0;
		//add noise also for composition field
		efxijk
		//Convert to fourier space
		cufftrc3k(&phi_r[i][0][0][0],&phi_k[i][0][0][0]);
	}//end of phi loop
}
#endif

#ifdef DISLOCATION_ACTIVITY_ON
if(stress_flag)
{
		if(cal_eps_plastic()>0.0015)
		{
		stress_flag=0;
		W=0.0;
		printf("Stress is turned off at timestep=%d; plastic_strain=%lf\n",timestep,cal_eps_plastic());
		
#pragma acc update host(ec_r,hphi_r,phi_r)
#ifdef DISLOCATION_ACTIVITY_ON
	#pragma acc update host(eta_r)
#endif
			OutputMed(timestep);
		}
}
#endif


}//end of GetEintC

