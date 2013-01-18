#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TINY 1.0e-20;

void ludcmp(double **a, int n, int *indx,double *d);
void lubksb(double **a, int n, int *indx, double b[]);
void invmat(double **a, double **y, int n);
void nrerror(char *error);
double *vector(int nl,int nh);
void free_vector(double *v,int nl,int nh);

void invmat(double **a, double **y, int n)
{

  double d,*col;
  int i,j,*indx;

  col = (double *)malloc(n*sizeof(double));
  indx = (int *)malloc(n*sizeof(int));
  ludcmp(a,n,indx,&d);
  for(j=0;j<n;j++){
    for(i=0; i < n; i++) col[i] = 0.0;
    col[j]=1.0;
    lubksb(a,n,indx,col);
    for(i=0;i<n;i++) y[i][j]=col[i];
  }
  free(col);
  free(indx);
}

void ludcmp(double **a, int n, int *indx,double *d)
{
  int i, imax, j,k;
  double big,dum,sum,temp;
  double *vv;
  
  vv=vector(0,n-1);
  *d=1.0;
  for(i=0; i<n; i++){
    big=0.0;
    for(j=0;j<n; j++)
      if((temp=fabs(a[i][j])) > big) big=temp;
    if(big==0.0) nrerror("Singular matrix in routine LUDCMP");
    vv[i]=1.0/big;
  }
  for(j=0; j<n;j++){
    for(i=0;i<j;i++){
      sum=a[i][j];
      for(k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big=0.0;
    for(i=j; i<n;i++){
      sum=a[i][j];
      for(k=0;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if((dum=vv[i]*fabs(sum)) >= big){
	big=dum;
	imax=i;
      }
    }
    if(j != imax){
      for(k=0; k <n; k++){
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if( a[j][j] == 0.0) a[j][j] = TINY;
    if(j !=n){
      dum=1.0/(a[j][j]);
      for(i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
  free_vector(vv,0,n-1);
}



void nrerror(char *error)
{
  fprintf(stderr,"error\n");
  fprintf(stderr,"%s\n",error);
  exit(1);
}

double *vector(int nl,int nh)
{
  double *v;
  
  v=(double*)malloc((unsigned)(nh-nl+1)*sizeof(double));
  if (!v) nrerror("allocation failure in vector()");
  return (v-nl);
}

void free_vector(double *v,int nl,int )
{
  free((char*)(v+nl));
}

void lubksb(double **a, int n, int *indx, double b[])
{
  int i,ii=-1,ip,j;
  double sum;
  
  for(i=0; i< n; i++){
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if(ii != -1)
      for(j=ii; j <= i-1; j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for(i=n-1; i>=0; i--){
    sum=b[i];
    for(j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}
