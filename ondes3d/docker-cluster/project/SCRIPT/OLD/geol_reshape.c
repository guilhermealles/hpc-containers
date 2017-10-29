//	F.D
//	Reshape du Modele Geol initial (Nice)
//	- on ne sur-echantillone pas le modele geologique
//	- on etend pas les limites du modeles geologiques
//	- le ratio ds_ini ds_new est un entier 
//
//	Generation du fichier np_max.in
//	Valeur maximale du nombre de coeurs admissible en fonction de la taille du domaine
//	La deuxieme valeur du fichier est l estimation de l elapsed time min en fonction du nombre max de procs 
//	A relire par OARSUB afin de fixer au mieux les ressources.
//
//	Fevrier 2017
/****************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include "assert.h"
static const int STRL=256;
int px,py,np;
int flag_np;
int xmin_new,xmax_new,ymin_new,ymax_new,zmin_new,zmax_new;



int compute_px_py()

{

	if (sqrt(np) == floor(sqrt(np)) ) 
	{
	px = sqrt (np);
	py = sqrt (np);
	} else {
	px = np;
	py = 1 ;

	 while ( px > py ) {

                        if (px % 2 == 0 ) {
                        py = 2*py;
                        px = px/2;
                        } else { break;}
	
			}


	 while ( px > py ) {

                        if (px % 3 == 0 ) {
                        py = 3*py;
                        px = px/3;
                        } else { break;}
	
			}



	while ( px > py ) {
		
			if (px % 5 == 0 ) {
			py = 5*py;
			px = px/5;
			} else  { break;}
			
			}
	}

return (0);
}


int check_np()


{
int XMIN,XMAX,YMIN,YMAX;
int DELTA=10;
int PX,PY,mpmx_tmp1,mpmy_tmp1,mpmx,mpmy;

flag_np = 1;

	XMIN = xmin_new;
	XMAX = xmax_new;
	YMIN = ymin_new;
	YMAX = ymax_new;


  mpmx = (double)( XMAX - XMIN + 2*DELTA + 3 )/(double)px;
  mpmx = floor(mpmx);
	
  if ( mpmx <= 10 )
  {
	flag_np = 0;
	//printf (" Reduce the number of processes used \n");
  }

  mpmy = (double)( YMAX - YMIN + 2*DELTA + 3 )/(double)py;
  mpmy = (int)floor(mpmy);

  if ( mpmy <= 10 )
  {
	flag_np = 0;
	//printf (" Reduce the number of processes used \n");
  }



}



char* File2str(const char* flname)
{
  FILE *fp=NULL;
  char* out=NULL;
  int size=0;

  char cur[STRL];

  fp=fopen(flname,"r");
  if (fp == NULL) {
    fprintf(stderr,"error when opening file: %s", flname);
    exit (EXIT_FAILURE) ;
  }

  int i=1;
  size = 1;
  out = malloc ( i* STRL * sizeof(char) );
  while ( fscanf(fp,"%s",&cur) != EOF ){
    size += strlen(cur);
    if ( size >= STRL-2 ){
      i = i+1;
      out=realloc(out, i* STRL * sizeof(char) );
      if ( out == NULL ){ perror("Error while converting file to string");}
      size=strlen(out)%STRL + 1;
    }
    strncat(out,cur,strlen(cur));
  }
  fclose(fp);

  return out;
}


/* OTHER FUNCTIONS (TO read PRMFILE ) */
/* =============== */
char* FindField(const char *cs, const char *ct)
{
  int indexB = -1;
  int indexE = -1;
  size_t sizect=strlen(ct);
  char *strB=NULL;
  char *strE=NULL;
  char *out=NULL;
  char *ctBegin=NULL;
  char *ctEnd=NULL;

  ctBegin=malloc(strlen(ct)+3);
  strcpy(ctBegin,"<");
  strcat(ctBegin,ct);
  strcat(ctBegin,">");

  ctEnd=malloc(strlen(ct)+4);
  strcpy(ctEnd,"</");
  strcat(ctEnd,ct);
  strcat(ctEnd,">");


	
  if (cs != NULL && ct != NULL)
    {
      strB = strstr(cs, ctBegin);

      if (strB != NULL) indexB = strB - cs + strlen(ctBegin);
      else{
        fprintf(stderr,"cannot find : %s in %s -- %s -- %s", ctBegin, cs,ct,ctEnd);
        exit (EXIT_FAILURE);
      }

      strE = strstr (cs, ctEnd);
      if (strE != NULL) indexE = strE - cs;
      else{
        fprintf(stderr,"cannot find : %s in %s", ctEnd, cs);
        exit (EXIT_FAILURE);
      }

    }
  else
    {
      assert(EXIT_SUCCESS==EXIT_FAILURE);
      return NULL;
    }
  out=malloc((indexE-indexB+1)*sizeof(char));

  strncpy(out,&cs[indexB],(indexE-indexB+1)*sizeof(char));
  out[indexE-indexB]='\0';

  free(ctBegin);
  free(ctEnd);

  return out;
}


double FindDble(const char* strIn, const char* srch ) /* work */
{
  double out=0.;
  char *cur=NULL;
  cur=FindField(strIn, srch);
  if ( cur != NULL ){ out=atof(cur); }
  else{
    fprintf(stderr,"Error in '%s' in file '%s', line %i", __func__, __FILE__, __LINE__);
    assert(EXIT_SUCCESS==EXIT_FAILURE);
    out = -1.;
  }
  free(cur);
  return out;
}


int FindInt(const char* strIn, const char* srch ) /* work */
{
  int out=0;
  char * cur=NULL;
  cur=FindField(strIn, srch);
  if ( cur != NULL ){ out=atoi(cur); }
  else{
    fprintf(stderr,"Error in '%s' in file '%s', line %i", __func__, __FILE__, __LINE__);
    assert(EXIT_SUCCESS==EXIT_FAILURE);
    out = -1;
  }
  free(cur);
  return out;
}






int main(int argc,char *argv[])
{


	int ratio_geol;

	int ngridx0,ngridy0,ngridz0;
	int ngridx_new,ngridy_new,ngridz_new;


  	int i,j,k,kz,jy,ix;
	float dx,dy,dz,x0,y0,z0;
	float ratio;

	float itemax,zcste;

	float x0_ini,y0_ini,z0_ini;
	float ds0,ds_new;

	float dummy_f;

	int cpt,cpt2,cpt1;
 
  	char prop_name[256];
	FILE *  fp1;
	FILE *  fp2;
	FILE *  fp3;

	int np_max;

/////////////////////////////////////////////////////////////
//	Lecture dans le fichier nice_geol.prm	
/////////////////////////////////////////////////////////////


	fp1 = fopen("nice_geol.prm", "r");
	fscanf (fp1, "%s %d", prop_name,&ngridx0);
	fscanf (fp1, "%s %d", prop_name,&ngridy0);
	fscanf (fp1, "%s %d", prop_name,&ngridz0);
	fscanf (fp1, "%s %f", prop_name,&x0_ini);
	fscanf (fp1, "%s %f", prop_name,&y0_ini);
	fscanf (fp1, "%s %f", prop_name,&z0_ini);
	fscanf (fp1, "%s %f", prop_name,&ds0);
	fscanf (fp1, "%s %f", prop_name,&ds0);
	fscanf (fp1, "%s %f", prop_name,&ds0);



//	Skip no_data value
	fscanf (fp1, "%s", prop_name);
	fscanf (fp1, "%s", prop_name);


  printf("Old BBOX grid points - X = %d - Y = %d - Z = %d \n",ngridx0,ngridy0,ngridz0);
  printf("Old BBOX origin - X = %f - Y = %f - Z = %f \n",x0_ini,y0_ini,z0_ini);

//////////////////////////////////////////////////////////////
//	Lecture dans le fichier nice2.prm
/////////////////////////////////////////////////////////////////
  char *filestr;

  filestr = File2str("nice2.prm");


  dummy_f=FindDble(filestr, "dt");
  ds_new=FindDble(filestr, "ds");



  itemax=FindInt(filestr, "tmax");
  xmin_new=FindInt(filestr, "xMin");
  xmax_new=FindInt(filestr, "xMax");

  ymin_new=FindInt(filestr, "yMin");
  ymax_new=FindInt(filestr, "yMax");

  zmin_new=FindInt(filestr, "zMin");
  zmax_new=FindInt(filestr, "zMax");

  x0=FindDble(filestr, "x0");
  y0=FindDble(filestr, "y0");
  z0=FindDble(filestr, "z0");




  printf("New BBOX grid points -- X = %d - Y = %d - Z = %d  \n",xmax_new,ymax_new,zmax_new);
  printf("New BBOX origin -- X = %f - Y = %f - Z = %f  \n",x0,y0,z0);


/////////////////////////////////////////////////////////////////////////////////::
//	Determination np_max

	np_max = 1;
	for (np=1;np<=256;np=np*2)
	{
		compute_px_py();
		check_np();
		if (flag_np == 1 ) np_max=np;
	}
	printf (" np max == %d ",np_max);

//	fp3 = fopen("np_max.in", "w");
//	fprintf(fp3,"%d \n",np_max);
//	fclose(fp3);	

///////////////////////////////////////////////////////////////////////////////////////////////////////////::


//	Generation du nouveau modele geologique
//	On pose une restriction forte en imposant
//	IL serait beaucoup plus simple de passer par les interfaces scud pour la gesion du modele Geol....
//	Probalement a faire evoluer ...


	ngridx_new=xmax_new-xmin_new+1;
	ngridy_new=ymax_new-ymin_new+1;
	ngridz_new=zmax_new-zmin_new+1;

	
	fp2 = fopen("nice_geol_2.prm","w");


	fprintf(fp2,"nx %d\n",ngridx_new);
	fprintf(fp2,"ny %d\n",ngridy_new);
	fprintf(fp2,"nz %d\n",ngridz_new);



	fprintf(fp2,"x0 %f\n",x0);
	fprintf(fp2,"y0 %f\n",y0);
	fprintf(fp2,"z0 %f\n",z0);

	fprintf(fp2,"dx %f\n",ds_new);
	fprintf(fp2,"dy %f\n",ds_new);
	fprintf(fp2,"dz %f\n",ds_new);

	fprintf(fp2,"nodata_value out \n");


//////////////////////////////////////////////////////////////////////////////////////////


//	Supposition pour les parametre du modele geol initial
//	xmin = 1 // ymin = 1 // zmin = 1 


//	Calcul en termes de COORD et pas de grid point

//	Parcours ancienne boite
//	Verification que le point courant appartient a la nouvelle emprise
//	Ecriture de la propriete (si point dans nouvelle boite) afin d extraire le nouveau modele


	ratio_geol = ds_new/ds0;
//	printf ("grid space -  old BBOX = %f // new BBOX = %f \n",ds0,ds_new);
//	printf ("grid space ratio = %d \n",ratio_geol); 

	cpt = 0;
	for ( k = 0; k < ngridz0; k+=ratio_geol)
	{
//	printf("Debug -Z - %d - %f -  %f  - %f \n",k,(z0 + (zmin_new-1)*ds_new),(z0_ini+k*ds0),(z0 + (zmax_new-1)*ds_new));

    	for ( j = 0; j < ngridy0; j+=ratio_geol)
	{
    	for ( i = 0; i < ngridx0; i+=ratio_geol)
	{

	
		
		fscanf (fp1, "%s", prop_name);

//		Controle que le point courant (old BBOX) est bien compris dans la nouvelle BBOX
		if ( ((z0_ini+k*ds0) <= (z0 + (zmax_new-1)*ds_new)) && ((z0_ini+k*ds0) >= (z0 + (zmin_new-1)*ds_new)) )
		{
		if ( ((y0_ini+j*ds0) < (y0 + ymax_new*ds_new)) && ((y0_ini+j*ds0) >= (y0 + (ymin_new-1)*ds_new)) )
		{
		if ( ((x0_ini+i*ds0) < (x0 + xmax_new*ds_new)) && ((x0_ini+i*ds0) >= (x0 + (xmin_new-1)*ds_new)) )
		{


		fprintf(fp2, "%s",prop_name);
		fprintf(fp2,"\n");
		cpt++;
		
		
		}	
		}
		}


	}
	}
	}


//	Controle nbre de pts dans fichier geol = nbre de pts fichier prm

	cpt2 = (zmax_new-zmin_new+1)*(ymax_new-ymin_new+1)*(xmax_new - xmin_new +1);
	if (cpt != cpt2 )
	{
	printf("ERROR -  Total points extraits de l'ancienne emprise  =%d /// Total theorique points nouvelle emprise = %d \n",cpt, cpt2);
	exit(0);
	}

	
	ratio = (ngridx0*ngridy0*ngridz0)/cpt2;
	printf("//// RATIO espace //// %f  \n",ratio);
		

        fp3 = fopen("np_max.in", "w");
	zcste=3e-10;
	printf (" np max == %d -- elapsed_time min= %e ",np_max,zcste*np_max*itemax*cpt2);
        fprintf(fp3,"%d %e  \n",np_max,zcste*np_max*itemax*cpt2);
        fclose(fp3);



	fclose(fp1);
	fclose(fp2);

}

