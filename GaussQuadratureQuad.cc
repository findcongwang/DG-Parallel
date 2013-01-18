#include "IntPt.h"
#include "GaussLegendreSimplex.h"
#include <stdio.h>

IntPt2d *getGQQPts(int order);
int getNGQQPts(int order);

IntPt2d * GQQ[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int GaussLegendre1D(int, double **,double **);

IntPt2d *getGQQPts(int order)
{ 
  order = (order==0)?1:order;
  int n = (order+3)/2;
  //  printf("n = %d\n",n);
  int index = n-2;
  if(!GQQ[index])
    {
      double *pt,*wt;
      GaussLegendre1D(n,&pt,&wt);
      GQQ[index] = new IntPt2d[n*n];
      int k = 0;
      for(int i=0; i < n; i++) {
	for(int j=0; j < n; j++) {
	  GQQ[index][k].pt[0] = pt[i];
	  GQQ[index][k].pt[1] = pt[j];
	  GQQ[index][k++].weight = wt[i]*wt[j];
	  //	  printf ("%f %f %f\n",pt[i],pt[j],wt[i]*wt[j]);
	}
      }
    }
  return GQQ[index];
}

int getNGQQPts(int order)
{ 
  order = (order==0)?1:order;
  return ((order+3)/2)*((order+3)/2);
}

