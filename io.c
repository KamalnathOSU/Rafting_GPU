#include "global_variables.h"
#include <string.h>
void write_vtk(double *rfield, char filename[], char fieldname[] );
void cal_vfrac();
void OutputMed(int timestep)
{
char filename[500];
int i;
#ifdef BIN_FILE
FILE *fp=NULL;
sprintf(filename,"output_files/time_%d.bin",timestep);
fp=fopen(filename,"wb");
if(fp)
{	fwrite(ecr,sizeof(double),LMN,fp);
	fwrite(phir,sizeof(double),V*LMN,fp);
	#ifdef DISLOCATION_ACTIVITY_ON
		fwrite(etar,sizeof(double),P*LMN,fp);
	#endif
	fwrite(hphir,sizeof(double),LMN,fp);
fclose(fp);
}
else
{fprintf(stderr,"#Error:Could open the file %s to write\n",filename);exit(0);}

#else 
//For producing vtk files
sprintf(filename,"output_files/ecr_t%d.vtk",timestep);
write_vtk(ecr,filename,"ecr");
#pragma omp parallel for private(filename)
for(i=0;i<V;i++)
	{sprintf(filename,"output_files/phi_%d_t%d.vtk",i+1,timestep);
	write_vtk(&phi_r[i][0][0][0],filename,"phi");
	}
	#ifdef DISLOCATION_ACTIVITY_ON
#pragma omp parallel for private(filename)
for(i=0;i<P;i++)
	{sprintf(filename,"output_files/eta_%d_t%d.vtk",i+1,timestep);
	write_vtk(&eta_r[i][0][0][0],filename,"eta");
	}
	#endif
sprintf(filename,"output_files/hphi_t%d.vtk",timestep);
write_vtk(hphir,filename,"hphi");
printf("timestep=%d completed \n",timestep);
fflush(stdout);
cal_vfrac();
#endif
}

void write_vtk(double *rfield, char filename[], char fieldname[] )
{
FILE *fp=NULL;
int i,j,k;
fp=fopen(filename,"w");
if(!fp)
	{fprintf(stderr,"#Error: %s could not be opened for writing\n",filename);exit(0);}

//writing vtk header
fprintf(fp,"# vtk DataFile Version 2.0\n");
fprintf(fp,"VTK from C-program\n");
fprintf(fp,"ASCII\n");
fprintf(fp,"DATASET STRUCTURED_POINTS\n");
fprintf(fp,"DIMENSIONS %d %d %d\n",L,M,N);
fprintf(fp,"SPACING 1 1 1\n");
fprintf(fp,"ORIGIN 0 0 0\n");
fprintf(fp,"POINT_DATA %d\n",L*M*N);
fprintf(fp,"SCALARS %s float 1\n",fieldname);
fprintf(fp,"LOOKUP_TABLE default\n");
for(k=0; k<N; k++)
	for(j=0; j<M; j++)
		for(i=0; i<L; i++)
			fprintf(fp,"%.6lf \n",rfield[i*M*N+j*N+k]);
fclose(fp);
}

void init_from_file(double *tpr, char *filename)
{
	int x_i,x_j,x_k,xijk,i;
	int l;
	FILE *fp=NULL;
	char ignore[40];
fp=fopen(filename,"r");
	if(fp)
	{
		l=strlen(filename);
		if(strcmp(&filename[l-4],".dat")==0)
		{	for_xijk
			fscanf(fp,"%lf",&tpr[xijk]);
			efor_xijk
		}
		else if(strcmp(&filename[l-4],".vtk")==0)
		{	for(i=0;i<10;i++)
				fgets(ignore,40,fp);//ignoring first ten lines	
			for_xijk
			fscanf(fp,"%lf",&tpr[x_i + x_j*L + x_k*L*M]);
			efor_xijk
		}
		else
		{fprintf(stderr,"#Error:Unknown input file format,  %s",filename);exit(0);}
		fclose(fp);
	}
	else
	{fprintf(stderr,"#Error: the file %s could not be opened for reading\n",filename);exit(0);
	}
}//end of init_from_file

