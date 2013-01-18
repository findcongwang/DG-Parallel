#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "FunctionSpace.h"
#include "Integrator.h"
#include "mEntity.h"
#include "Mapping.h"
#include <vector>
using namespace std;

extern  double ** allocateMatrix(int n);
extern  void freeMatrix( double **v);

void ScalarProducts(mEntity *e, 
		    Mapping *Map,
		    FunctionSpace *fs, 
		    double **matr)
{
  int i,j,k;
  int fSize = fs->size();
  for(i=0;i<fSize;i++)
    for(j=0;j<fSize;j++)(matr)[i][j] = 0.0;
  GaussIntegrator gauss;  
  double u,v,w,weight;
  int order = 2*fs->order();
  double *fct;
  fct = new double [fSize];
  double volume = 0.0;

  for(i=0;i<gauss.nbIntegrationPoints(e,order);i++)
    {
      gauss.iPoint(e,i,order,u,v,w,weight); 
      double detJac = Map->detJac(u,v,w);
      fs->fcts(u,v,w,fct);
      volume += detJac * weight;
      for(j=0;j<fSize;j++)
	{
	  for(k=0;k<fSize;k++)
	    {
	      (matr) [j][k] += detJac * fct[k] * fct[j] * weight;
	    }
	}
    }
  delete [] fct;
  printf("volume = %12.5E\n",volume);
}

static void print (int n,double **mat,FILE *f,char name[256])
{
  fprintf(f,"const int maxsize%s = %d;\n",name,n);
  fprintf(f,"static double %s[maxsize%s][maxsize%s] = {\n",name,name,name);
  for(int i=0;i<n;i++)
    {
      fprintf(f," {");
      for(int j=0;j<n;j++)
	{
	  if(j)fprintf(f,",%22.15E",mat[i][j]);
	  else fprintf(f,"%22.15E",mat[i][j]);
	}
      if(i!=n-1)fprintf(f,"},\n");
      else fprintf(f,"}\n");
    }
  fprintf(f,"};\n");
}

double ScalarProduct (double ** mat , double **sp, int f1, int f2)
{
  double SP = 0.0;
  for(int i = 0; i<= f1; i++)
    {
      for(int j = 0; j<= f2; j++)
	{
	  SP += sp[i][j] * mat[f1][i] * mat[f2][j];
	}
    }
  return SP;
}


double ** Orthogonalize(mEntity *e, 
			Mapping *Map,
			FunctionSpace *fs)
{
  double **gs;
  double **sp;
  double **figj;

  gs = allocateMatrix(fs->size());
  sp = allocateMatrix(fs->size());
  figj = allocateMatrix(fs->size());

  int order = fs->order();
  int n = fs->size();

  // Init scalar products and gauss integrations

  ScalarProducts(e,Map,fs,sp);
  print (n,sp,stdout,"toto");
  int i,k,j,l;
  for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
	{
	  gs[i][j] = 0.0;
	}
    }
  
  double normgj2[10000];
  
 
  for(i=0;i<n;i++)
    {
      gs[i][i] = 1.0;
      for(j=0;j<i;j++)
	{
	  for(k=j;k<i;k++)
	    {
	      gs[i][j] -= (gs[k][j]*figj[k][i]/normgj2[k]);
	    }
	}
      normgj2[i] = 0.0;

      // compute the norm of the
      // new vector using scalar products

      for(k=0;k<=i;k++) 
	for(l=0;l<=i;l++)
	  {
	    normgj2[i] +=  (gs[i][k] * gs[i][l] * sp[k][l]);
	  }	  

      // let us now orthonorm the line
      for(k=0;k<=i;k++) gs[i][k] /= sqrt(normgj2[i]);
      //      printf("norm[%d] = %22.15e\n",i,normgj2[i]);
      normgj2[i] = 1.0;


      // we have computed gi
      // we compute now <gi,fj> = Sum_k gs[i][k] * <fk,fj>

      for(j=0;j<n;j++)
	{
	  figj[i][j] = 0.0;
	  for(k=0;k<=i;k++)
	    {
	      figj[i][j] += sp[k][j] * gs[i][k];
	    }
	}
    }    
  freeMatrix(figj);
  freeMatrix(sp);

  FILE *f = fopen("TetOrtho.h","w");

  print(n,gs,f,"myTetOrtho");

  fclose(f);

  return gs;
}
