
#include "IntPt.h"
#include "GaussLegendreSimplex.h"
#include <stdio.h>

IntPt3d *getGQHPts(int order);
int getNGQHPts(int order);

IntPt3d * GQH[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int GaussLegendre1D(int, double **,double **);

IntPt3d *getGQHPts(int order)
{ 
  order = (order==0)?1:order;
  int n = (order+3)/2;
  //  printf("n = %d\n",n);
  int index = n-2;
  if(!GQH[index])
    {
      double *pt,*wt;
      GaussLegendre1D(n,&pt,&wt);
      GQH[index] = new IntPt3d[n*n*n];
      int l = 0;
      for(int i=0; i < n; i++) {
	for(int j=0; j < n; j++) {
	  for(int k=0; k < n; k++) {
	    GQH[index][l].pt[0] = pt[i];
	    GQH[index][l].pt[1] = pt[j];
	    GQH[index][l].pt[2] = pt[k];
	    GQH[index][l++].weight = wt[i]*wt[j]*wt[k];
	    //	    printf ("%f %f %f %f\n",pt[i],pt[j],pt[k],wt[i]*wt[j]*wt[k]);
	  }
	}
      }
    }
  return GQH[index];
}

int getNGQHPts(int order)
{ 
  order = (order==0)?1:order;
  return (order+3)*(order+3)*(order+3)/8;
}

