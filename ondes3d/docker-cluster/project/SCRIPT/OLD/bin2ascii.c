//	F.D
//
//	Passage de obs1.bin a obs1.dat
//	Passage de PGV.bin a PGV.vtk
//
//
//	Juillet 2017
/****************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include "assert.h"
#define NR_END 0
#define DELTA 10
#define XMIN -97
#define XMAX 100
#define YMIN -100
#define YMAX 100


static const int timestep_number=1000;
static const int freq_io=1;


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}





float **fmatrix(long nrl, long nrh, long ncl, long nch)
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
 
	/* allocate pointers to rows */
	m=malloc((size_t)((nrow+NR_END)*sizeof(float*)));

	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
 	m[nrl]= malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}





int main()
{
 
	FILE *  fp1;
	double * buff1,*buff2,*buff3;
	double * sismo;
	float ** PGVglobal;
	int l1,l2,nb_io_ops;
	int i,j;

//	Allocation
	buff1 = (double *)malloc(freq_io*sizeof(double));
	buff2 = (double *)malloc(freq_io*sizeof(double));
	buff3 = (double *)malloc(freq_io*sizeof(double));
	sismo = (double *)malloc(timestep_number*sizeof(double));

//	Lecture sismo complet
	fp1 = fopen("obs1.bin", "r");


	nb_io_ops = timestep_number / freq_io;
	
	for ( l1 = 0; l1 < nb_io_ops; l1++)
	{ 
	fread( buff1 , sizeof(double) , freq_io , fp1);	
	fread( buff2 , sizeof(double) , freq_io , fp1);
	fread( buff3 , sizeof(double) , freq_io , fp1);
//	printf ("%d %15e \n",l1,buff1[0]);
//	printf ("%d %15e \n",l1,buff2[0]);
//	printf ("%d %15e \n",l1,buff3[0]);

	      	for ( l2 = 0; l2 < freq_io; l2++)
	        {
		sismo[l2+l1*freq_io] = buff1[l2];
		} 

	}
	fclose(fp1);

//	Ecriture fichier ascii
	fp1 = fopen("obs1bin.dat", "w");
	fprintf(fp1, "%d  \n", nb_io_ops);
	fprintf(fp1, "%d \n", nb_io_ops);

	for ( l1 = 0; l1 < timestep_number; l1++)
	{
		fprintf(fp1, "%15e \n", sismo[l1]);
	}
	fclose(fp1);
	free (buff1);
	free(buff2);
	free(buff3);

////////////////////////////////////////////////////////////////////////::


//	Lecture PGV
	fp1 = fopen("PGV.bin", "r");


	 PGVglobal = fmatrix(XMIN-DELTA, XMAX+DELTA+2, YMIN-DELTA, YMAX+DELTA+2);

	for ( i = XMIN-DELTA+1; i <= XMAX+DELTA+1; i++){
        for ( j = YMIN-DELTA+1; j <= YMAX+DELTA+1; j++){
	fread(&PGVglobal[i][j] , sizeof(float) , 1 , fp1);	
	printf ("%d %d %15e \n",i,j,PGVglobal[i][j]);
	}
	}

	fclose(fp1);
	return 0;
}
