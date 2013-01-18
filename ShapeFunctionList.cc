#include <math.h>
#include "mEntity.h"
#include "Edge.h"
#include "Face.h"
#include "mRegion.h"
#include "ShapeFunctionList.h"
#include <assert.h>

void PascalgetIndices(int iFct, int &n, int &i)
{
  int k = 0;
  int l = 0;
  while(k<=iFct)
    {
      l++;
      k +=l;
    };
  n = l - 1;
  i = l - k + iFct;
}

int SHF::getSize(int order, int dim)
{
	switch (dim) 
	{
	case 1 : return order + 1;
	case 2 : return (order+1)*(order+2)/2;
	case 3 : return (order+1)*(order+2)*(order+3)/6;
	default : throw "Error : unknown dimentions";
	}
}

int SHF::getQSize(int order, int dim)
{
	switch (dim) 
	{
	case 1 : return order + 1;
	case 2 : return (order+1)*(order+1);
	case 3 : return (order+1)*(order+1)*(order+1);
	default : throw "Error:dimentions";
	}
}



/*
  Shape functions for discontinuous elements
  -------------------------------------------------
  - Pascal : using simple Pascal simplex
    For triangles , fi are:

          |  i=0         i=1        i=2     ...
      ----+----------------------------------------
      n=0 |   1
      n=1 |   u           v
      n=2 |   u*u         u*v        v*v 
	  |   u*u*u       u*u*v      u*v*v    v*v*v
          |
	  |                n-i-1 i
      n=n |               u     v

    Note that some fcts are precalculated for sake of
    performance.

  We have also a set of orthonormal shape functions where
  the norm is evaluated in the reference element. This
  orthogonality is transfered to the real element but not
  orthonormality. I have divided rows of the gram-schmidt
  matrix by the norm in order to improve conditioning.


  Hilbert space L2(Triangle) with inner product :

          1  1-u
           - -
          |  |
  <a,b> = |  | a b du dv
          |  |
         -  -
	 0  0

  and the induced norm :

       2
  ||a||   = <a,a>


  orthogonalization procedure (Gram-Schmidt) :

  for i = 0 ... N if N is the size of the function space

             i-1
  gi' = fi - Sum <fi,gj> gj
             j=0

  gi = gi'/||gi'||

  This calculation are made in an exterior program. The result
  is an orthogonalization matrix which is given on a header file :
  
  gi = Sum mat[i][j] * fj

*/

SHF::SHF_T SHF::FunctionSpace = SHF_Pascal;

// GramSchmidt for triangle"
#include "OrthoSFTri.h"
// GramSchmidt for tet"
#include "TetOrtho.h"

double SHF::getTetFctOrthogonal (int iFct, double u, double v, double w)
{
  assert (iFct < maxsizemyTetOrtho);
  double r = 0.0;
  for(int j=0;j<=iFct;j++) r += SHF::getTetFct(j,u,v,w) * myTetOrtho[iFct][j];
  return r;
}

double SHF::PascalgetFctOrthogonal (int iFct, double u, double v)
{
  assert (iFct < maxsizemyGramSchmidtTriangle);
  double r = 0.0;
  for(int j=0;j<=iFct;j++) r += SHF::PascalgetFct(j,u,v) * myGramSchmidtTriangle[iFct][j];
  return r;
}

double SHF::getTetdFctduOrthogonal (int iFct, double u, double v, double w)
{
  assert (iFct < maxsizemyTetOrtho);
  double r = 0.0;
  for(int j=0;j<=iFct;j++) r += SHF::getTetdFctdu(j,u,v,w) * myTetOrtho[iFct][j];
  return r;
}

double SHF::getTetdFctdvOrthogonal (int iFct, double u, double v, double w)
{
  assert (iFct < maxsizemyTetOrtho);
  double r = 0.0;
  for(int j=0;j<=iFct;j++) r += SHF::getTetdFctdv(j,u,v,w) * myTetOrtho[iFct][j];
  return r;
}

double SHF::getTetdFctdwOrthogonal (int iFct, double u, double v, double w)
{
  assert (iFct < maxsizemyTetOrtho);
  double r = 0.0;
  for(int j=0;j<=iFct;j++) r += SHF::getTetdFctdw(j,u,v,w) * myTetOrtho[iFct][j];
  return r;
}

double SHF::PascalgetdFctduOrthogonal (int iFct, double u, double v)
{
  assert (iFct < maxsizemyGramSchmidtTriangle);
  double r = 0.0;
  for(int j=0;j<=iFct;j++) r += SHF::PascalgetdFctdu(j,u,v) * myGramSchmidtTriangle[iFct][j];
  return r;
}

double SHF::PascalgetdFctGeneralOrthogonal (int iFct, int ithDeriv, int jth, double u, double v)
{
  assert (iFct < maxsizemyGramSchmidtTriangle);
  double r = 0.0;
  for(int j=0;j<=iFct;j++) r += SHF::PascalgetdFctGeneral(j,ithDeriv,jth,u,v) * myGramSchmidtTriangle[iFct][j];
  return r;
}


double SHF::PascalgetdFctdvOrthogonal (int iFct, double u, double v)
{
  assert (iFct < maxsizemyGramSchmidtTriangle);
  double r = 0.0;
  for(int j=0;j<=iFct;j++) r += SHF::PascalgetdFctdv(j,u,v) * myGramSchmidtTriangle[iFct][j];
  return r;
}

void SHF::init (SHF::SHF_T fcts)
{
  FunctionSpace = fcts;
}

double SHF::getFct (int i, double u, Edge *e)
{
  switch (FunctionSpace)
    {
    case SHF_Pascal :
      return PascalgetFct(i,u);
    case SHF_Pascal_Orthogonal :
      return PascalgetFct(i,u);
	default : throw "Error: unknown space";
    }
}

double SHF::getFct (int i, double u, double v, Face *f)
{
  switch (FunctionSpace)
    {
    case SHF_Pascal :
      return PascalgetFct(i,u,v);
    case SHF_Pascal_Orthogonal :
      return PascalgetFctOrthogonal(i,u,v);
	default : throw "Error: unknown space";
    }
}


double SHF::getdFctGeneral (int i, int derivOrder, int j, double u, double v, Face *f)
{
  switch (FunctionSpace)
    {
    case SHF_Pascal :
      return PascalgetdFctGeneral(i,derivOrder,j,u,v);
    case SHF_Pascal_Orthogonal :
      return PascalgetdFctGeneralOrthogonal(i,derivOrder,j,u,v);
	default : throw "Error: unknown space ";
    }
}

double SHF::getdFctdu (int i, double u, double v, Face *f)
{
  switch (FunctionSpace)
    {
    case SHF_Pascal :
      return PascalgetdFctdu(i,u,v);
    case SHF_Pascal_Orthogonal :
      return PascalgetdFctduOrthogonal(i,u,v);
	default : throw "Error: unknown space type";
    }
}

double SHF::getdFctdv (int i, double u, double v, Face *f)
{
  switch (FunctionSpace)
    {
    case SHF_Pascal :
      return PascalgetdFctdv(i,u,v);
    case SHF_Pascal_Orthogonal :
      return PascalgetdFctdvOrthogonal(i,u,v);
	default : throw "Error: unknown space type";
    }
}

double SHF::getdFctdu (int i, double u, Edge *e)
{
  switch (FunctionSpace)
    {
    case SHF_Pascal_Orthogonal :
    case SHF_Pascal : return PascalgetdFctdu(i,u);
    default : throw "Error: unknown space type";
    }
}


double SHF::getdFctduu (int i, double u, double v, Face *f)
{
  throw "sh fct not coded yet";  
}

double SHF::getdFctduv (int i, double u, double v, Face *f)
{
      throw "sh fct not coded yet";
}

double SHF::getdFctdvv (int i, double u, double v, Face *f)
{
  throw "sh fct not coded yet"; 
}


double SHF::getdFctduu (int i, double u, Edge *e)
{
  switch (FunctionSpace)
    {
    case SHF_Pascal :
    case SHF_Pascal_Orthogonal :
    case SHF_Fourier : return PascalgetdFctduu(i,u);
    default : throw "Error: unknown space type";
    }
}


/********************************************
  P a s c a l   S h a pe   F u n c t i o n s					
*********************************************/

double SHF::PascalgetFct (int i, double u)
{
  switch (i)
    {
    case -2: return 0.0;
    case -1: return 0.0;
    case 0 : return 1.0;
    case 1 : return u;
    case 2 : return u*u;
    case 3 : return u*u*u;
    default: return u* PascalgetFct(i-1,u);
    }
}

double SHF::PascalgetdFctdu (int i, double u)
{
  return (double) i * PascalgetFct(i-1,u);
}

double SHF::PascalgetdFctduu (int i, double u)
{
  return (double) (i-1) * PascalgetdFctdu(i-1,u);
}


double SHF::PascalgetFct (int iFct, double u, double v)
{
  switch (iFct)
    {
    case 0 : return 1.0;
    case 1 : return u;
    case 2 : return v;
    case 3 : return u*u;
    case 4 : return u*v;
    case 5 : return v*v;
    case 6 : return u*u*u;
    case 7 : return u*u*v;
    case 8 : return v*v*u;
    case 9 : return v*v*v;
    case 10 : return u*u*u*u;
    case 11 : return u*u*u*v;
    case 12 : return v*v*u*u;
    case 13 : return v*v*v*u;
    case 14 : return v*v*v*v;
    case 15 : return u*u*u*u*u;
    case 16 : return u*u*u*u*v;
    case 17 : return v*v*u*u*u;
    case 18 : return v*v*v*u*u;
    case 19 : return v*v*v*v*u;
    case 20 : return v*v*v*v*v;
    case 21 : return u*u*u*u*u*u;
    case 22 : return u*u*u*u*u*v;
    case 23 : return v*v*u*u*u*u;
    case 24 : return v*v*v*u*u*u;
    case 25 : return v*v*v*v*u*u;
    case 26 : return v*v*v*v*v*u;
    case 27 : return v*v*v*v*v*v;
    default:
     // PascalgetIndices(iFct,n,i);
	  //returns garbage!!! n is underfined
     // return pow(u,n-i) * pow(v,i);
	  return 0.0;
    }
}

double SHF::PascalgetdFctGeneral (int iFct, int ithDeriv, int jthTerm, double u, double v)
{

  if(ithDeriv == 1)
    {
      if(jthTerm == 0)return PascalgetdFctdu(iFct,u,v);
      if(jthTerm == 1)return PascalgetdFctdv(iFct,u,v);
    }

  if(ithDeriv == 2)
    {
      if(jthTerm == 0)return PascalgetdFctduu(iFct,u,v);
      if(jthTerm == 1)return PascalgetdFctduv(iFct,u,v);
      if(jthTerm == 2)return PascalgetdFctdvv(iFct,u,v);
    }


  int deriv_u = ithDeriv - jthTerm;
  int deriv_v = jthTerm;
  int n,i,k;
  PascalgetIndices(iFct,n,i);
  // f = pow (u,n-i)pow(v,i)
  if(deriv_u > n-i)return 0.0;
  if(deriv_v > i)return 0.0;

  double coef1 = 1.0;
  for(k=n-i;k>deriv_u;k--)coef1*=k;
  double coef2 = 1.0;
  for(k=i;k>deriv_v;k--)coef2*=k;

  return pow(u,n-i-deriv_u) * pow(u,i-deriv_v) * coef1 * coef2;

}

double SHF::PascalgetdFctduu (int iFct, double u, double v)
{
  int n,i;
  switch (iFct)
    {
    case 0 : return 0.0;
    case 1 : return 0.0;
    case 2 : return 0.0;
    case 3 : return 2.0;
    case 4 : return 0.0;
    case 5 : return 0.0;
    case 6 : return 6*u;
    case 7 : return 2*v;
    case 8 : return 0.0;
    case 9 : return 0.0;
    default:
      PascalgetIndices(iFct,n,i);
      if(n-i-2 <=0)return 0.0;
      return (double)((n-i)*(n-i-1))*pow(u,n-i-2) * pow(v,i);
    }
}

double SHF::PascalgetdFctduv (int iFct, double u, double v)
{
  int n,i;
  switch (iFct)
    {
    case 0 : return 0.0;
    case 1 : return 0.0;
    case 2 : return 0.0;
    case 3 : return 0.0;
    case 4 : return 1.0;
    case 5 : return 0.0;
    case 6 : return 0.0;
    case 7 : return 2*u;
    case 8 : return 2.*v;
    case 9 : return 0.0;
    default:
      PascalgetIndices(iFct,n,i);
      if(n-i-1 <=0)return 0.0;
      if(i-1 <=0)return 0.0;
      return (double)((n-i)*(i-1))*pow(u,n-i-1) * pow(v,i-1);
    }
}

double SHF::PascalgetdFctdu (int iFct, double u, double v)
{
  int n,i;
  switch (iFct)
    {
    case 0 : return 0.0;
    case 1 : return 1.0;
    case 2 : return 0.0;
    case 3 : return 2.0*u;
    case 4 : return v;
    case 5 : return 0.0;
    case 6 : return 3*u*u;
    case 7 : return 2*u*v;
    case 8 : return v*v;
    case 9 : return 0.0;
    case 10 : return 4*u*u*u;
    case 11 : return 3*u*u*v;
    case 12 : return 2*v*v*u;
    case 13 : return v*v*v;
    case 14 : return 0.0;
    case 15 : return 5*u*u*u*u;
    case 16 : return 4*u*u*u*v;
    case 17 : return 3*v*v*u*u;
    case 18 : return 2*v*v*v*u;
    case 19 : return v*v*v*v;
    case 20 : return 0.0;
    case 21 : return 6*u*u*u*u*u;
    case 22 : return 5*u*u*u*u*v;
    case 23 : return 4*v*v*u*u*u;
    case 24 : return 3*v*v*v*u*u;
    case 25 : return 2.*v*v*v*v*u;
    case 26 : return v*v*v*v*v;
    case 27 : return 0.0;
    default:
      PascalgetIndices(iFct,n,i);
      if(n-i-1 <0)return 0.0;
      return (double)(n-i)*pow(u,n-i-1) * pow(v,i);
    }
}

double SHF::PascalgetdFctdv (int iFct, double u, double v)
{
  int n,i;
  switch (iFct)
    {
    case 0 : return 0.0;
    case 1 : return 0.0;
    case 2 : return 1.0;
    case 3 : return 0.0;
    case 4 : return u;
    case 5 : return 2.*v;
    case 6 : return 0.;
    case 7 : return u*u;
    case 8 : return 2*u*v;
    case 9 : return 3*v*v;
    case 10 : return 0.0;
    case 11 : return u*u*u;
    case 12 : return 2*v*u*u;
    case 13 : return 3*v*v*u;
    case 14 : return 4*v*v*v;
    case 15 : return 0.0;
    case 16 : return u*u*u*u;
    case 17 : return 2*v*u*u*u;
    case 18 : return 3*v*v*u*u;
    case 19 : return 4*v*v*v*u;
    case 20 : return 5*v*v*v*v;
    case 21 : return 0.0;
    case 22 : return u*u*u*u*u;
    case 23 : return 2*v*u*u*u*u;
    case 24 : return 3*v*v*u*u*u;
    case 25 : return 4*v*v*v*u*u;
    case 26 : return 5*v*v*v*v*u;
    case 27 : return 6*v*v*v*v*v;
    default:
      PascalgetIndices(iFct,n,i);
      if(i-1 <0)return 0.0;
      return (double)(i)*pow(u,n-i) * pow(v,i-1);
    }
}

double SHF::PascalgetdFctdvv (int iFct, double u, double v)
{
  int n,i;
  switch (iFct)
    {
    case 0 : return 0.0;
    case 1 : return 0.0;
    case 2 : return 0.0;
    case 3 : return 0.0;
    case 4 : return 0.0;
    case 5 : return 2.0;
    case 6 : return 0.0;
    case 7 : return 0.0;
    case 8 : return 2*u;
    case 9 : return 6*v;
    default:
      PascalgetIndices(iFct,n,i);
      if(i-2 <=0)return 0.0;
      return (double)(i*(i-1))*pow(u,n-i) * pow(v,i-2);
    }
}


/*********** QUADS ****************/

double LegendrePolynomial (int i, double u)
{
  switch(i)
    {
    case 0: return 1.0;
    case 1: return u;
    case 2: return (1.5*u*u  - 0.5);
    case 3: return 2.5*u*u*u - 1.5 * u;
    case 4: return (35.*u*u*u*u - 30.*u*u + 3.)*.125;
    case 5: return (65.*u*u*u*u*u - 70.*u*u*u + 15.*u)*.125;
    case 6: return (231.*u*u*u*u*u*u - 315.*u*u*u*u + 105.*u*u-5.)*.0625;
	default: throw " Orders higher than 6 are not coded yet ";
    }
}

double dLegendrePolynomial (int i, double u)
{
  switch(i)
    {
    case 0: return 0.0;
    case 1: return 1.0;
    case 2: return 3.*u;
    case 3: return 7.5*u*u - 1.5;
    case 4: return (35.*u*u*u - 15.*u)*0.5;
    case 5: return 5.*(97./42.)*u*u*u*u - (75./7.) *u*u + (53./42.);
	default: throw " Orders higher than 5 are not coded yet"; 
    }
}

double SHF::getQFct (int iFct, double u, double v, Face *f)
{
  switch(iFct)
    {
    case 0:return LegendrePolynomial(0,u)*.5;
    case 1:return LegendrePolynomial(1,u)*.5*sqrt(3.);
    case 2:return LegendrePolynomial(1,v)*.5*sqrt(3.);
    case 3:return LegendrePolynomial(1,u) * LegendrePolynomial(1,v)*1.5;
    case 4:return LegendrePolynomial(2,u)*0.5*sqrt(5.);
    case 5:return LegendrePolynomial(2,v)*0.5*sqrt(5.);
    case 6:return LegendrePolynomial(2,u) * LegendrePolynomial(1,v)*sqrt(15.)*0.5;
    case 7:return LegendrePolynomial(1,u) * LegendrePolynomial(2,v)*sqrt(15.)*0.5;
    case 8:return LegendrePolynomial(2,u) * LegendrePolynomial(2,v)*2.5;
    case 9:return LegendrePolynomial(3,u)*.5*sqrt(7.);
    case 10:return LegendrePolynomial(3,v)*.5*sqrt(7.);
    case 11:return LegendrePolynomial(3,u) * LegendrePolynomial(1,v)*.5*sqrt(21.);
    case 12:return LegendrePolynomial(1,u) * LegendrePolynomial(3,v)*.5*sqrt(21.);
    case 13:return LegendrePolynomial(3,u) * LegendrePolynomial(2,v)*.5*sqrt(35.);
    case 14:return LegendrePolynomial(2,u) * LegendrePolynomial(3,v)*.5*sqrt(35.);
    case 15:return LegendrePolynomial(3,u) * LegendrePolynomial(3,v)*3.5;
    case 16:return LegendrePolynomial(4,u)*1.5;
    case 17:return LegendrePolynomial(4,v)*1.5;
    case 18:return LegendrePolynomial(4,u) * LegendrePolynomial(1,v)*1.5*sqrt(3.);
    case 19:return LegendrePolynomial(1,u) * LegendrePolynomial(4,v)*1.5*sqrt(3.);
    case 20:return LegendrePolynomial(4,u) * LegendrePolynomial(2,v)*1.5*sqrt(5.);
    case 21:return LegendrePolynomial(2,u) * LegendrePolynomial(4,v)*1.5*sqrt(5.);
    case 22:return LegendrePolynomial(4,u) * LegendrePolynomial(3,v)*1.5*sqrt(7.);
    case 23:return LegendrePolynomial(3,u) * LegendrePolynomial(4,v)*1.5*sqrt(7.);
    case 24:return LegendrePolynomial(4,u) * LegendrePolynomial(4,v)*4.5;
    default: throw " Orders higher than 4 are not coded yet "; 
    }
}

double SHF::getQdFctdu (int iFct, double u, double v, Face *f)
{
  switch(iFct)
    {
    case 0:return 0.0;
    case 1:return dLegendrePolynomial(1,u)*.5*sqrt(3.);
    case 2:return 0.0;
    case 3:return dLegendrePolynomial(1,u) * LegendrePolynomial(1,v)*1.5;
    case 4:return dLegendrePolynomial(2,u) * 0.5*sqrt(5.);
    case 5:return 0.0;
    case 6:return dLegendrePolynomial(2,u) * LegendrePolynomial(1,v)*sqrt(15.)*0.5;
    case 7:return dLegendrePolynomial(1,u) * LegendrePolynomial(2,v)*sqrt(15.)*0.5;
    case 8:return dLegendrePolynomial(2,u) * LegendrePolynomial(2,v)*2.5;
    case 9:return  dLegendrePolynomial(3,u)*.5*sqrt(7.);
    case 10:return 0.0;
    case 11:return dLegendrePolynomial(3,u) * LegendrePolynomial(1,v)*.5*sqrt(21.);
    case 12:return dLegendrePolynomial(1,u) * LegendrePolynomial(3,v)*.5*sqrt(21.);
    case 13:return dLegendrePolynomial(3,u) * LegendrePolynomial(2,v)*.5*sqrt(35.);
    case 14:return dLegendrePolynomial(2,u) * LegendrePolynomial(3,v)*.5*sqrt(35.);
    case 15:return dLegendrePolynomial(3,u) * LegendrePolynomial(3,v)*3.5;
	case 16:return dLegendrePolynomial(4,u)*1.5;
    case 17:return 0.0;
    case 18:return dLegendrePolynomial(4,u) * LegendrePolynomial(1,v)*1.5*sqrt(3.);
    case 19:return dLegendrePolynomial(1,u) * LegendrePolynomial(4,v)*1.5*sqrt(3.);
    case 20:return dLegendrePolynomial(4,u) * LegendrePolynomial(2,v)*1.5*sqrt(5.);
    case 21:return dLegendrePolynomial(2,u) * LegendrePolynomial(4,v)*1.5*sqrt(5.);
    case 22:return dLegendrePolynomial(4,u) * LegendrePolynomial(3,v)*1.5*sqrt(7.);
    case 23:return dLegendrePolynomial(3,u) * LegendrePolynomial(4,v)*1.5*sqrt(7.);
    case 24:return dLegendrePolynomial(4,u) * LegendrePolynomial(4,v)*4.5;
	default: throw " Orders higher than 4 are not coded yet "; 
    }
}

double SHF::getQdFctdv (int iFct, double u, double v, Face *f)
{
  switch(iFct)
    {
    case 0:return 0.0;
    case 1:return 0.0;
    case 2:return dLegendrePolynomial(1,v)*.5*sqrt(3.);
    case 3:return LegendrePolynomial(1,u) * dLegendrePolynomial(1,v)*1.5;
    case 4:return 0.0;
    case 5:return dLegendrePolynomial(2,v)* 0.5*sqrt(5.);
    case 6:return LegendrePolynomial(2,u) * dLegendrePolynomial(1,v)*sqrt(15.)*0.5;
    case 7:return LegendrePolynomial(1,u) * dLegendrePolynomial(2,v)*sqrt(15.)*0.5;
    case 8:return LegendrePolynomial(2,u) * dLegendrePolynomial(2,v)*2.5;
    case 9:return 0.0;
    case 10:return dLegendrePolynomial(3,v)*.5*sqrt(7.);
    case 11:return LegendrePolynomial(3,u) * dLegendrePolynomial(1,v)*.5*sqrt(21.);
    case 12:return LegendrePolynomial(1,u) * dLegendrePolynomial(3,v)*.5*sqrt(21.);
    case 13:return LegendrePolynomial(3,u) * dLegendrePolynomial(2,v)*.5*sqrt(35.);
    case 14:return LegendrePolynomial(2,u) * dLegendrePolynomial(3,v)*.5*sqrt(35.);
    case 15:return LegendrePolynomial(3,u) * dLegendrePolynomial(3,v)*3.5;
	case 16:return 0.0;
    case 17:return dLegendrePolynomial(4,v)*1.5;
    case 18:return LegendrePolynomial(4,u) * dLegendrePolynomial(1,v)*1.5*sqrt(3.);
    case 19:return LegendrePolynomial(1,u) * dLegendrePolynomial(4,v)*1.5*sqrt(3.);
    case 20:return LegendrePolynomial(4,u) * dLegendrePolynomial(2,v)*1.5*sqrt(5.);
    case 21:return LegendrePolynomial(2,u) * dLegendrePolynomial(4,v)*1.5*sqrt(5.);
    case 22:return LegendrePolynomial(4,u) * dLegendrePolynomial(3,v)*1.5*sqrt(7.);
    case 23:return LegendrePolynomial(3,u) * dLegendrePolynomial(4,v)*1.5*sqrt(7.);
    case 24:return LegendrePolynomial(4,u) * dLegendrePolynomial(4,v)*4.5;
	default: throw " Orders higher than 4 are not coded yet ";
    }
}


// T E T S

double SHF::getTetFct (int iFct, double u, double v, double w)
{
  switch (iFct)
    {
    case 0 : return 1.0;
    case 1 : return u;
    case 2 : return v;
    case 3 : return w;
    case 4 : return u*u;
    case 5 : return v*v;
    case 6 : return w*w;
    case 7 : return u*v;
    case 8 : return u*w;
    case 9 : return v*w;
    case 10 : return u*u*u;
    case 11 : return v*v*v;
    case 12 : return w*w*w;
    case 13 : return u*u*v;
    case 14 : return u*u*w;
    case 15 : return v*v*u;
    case 16 : return v*v*w;
    case 17 : return w*w*u;
    case 18 : return w*w*v;
    case 19 : return u*v*w;
    case 20 : return u*u*u*u;
    case 21 : return v*v*v*v;
    case 22 : return w*w*w*w;
    case 23 : return u*u*u*v;
    case 24 : return u*u*u*w;
    case 25 : return v*v*v*u;
    case 26 : return v*v*v*w;
    case 27 : return w*w*w*u;
    case 28 : return w*w*w*v;
    case 29 : return u*u*v*v;
    case 30 : return u*u*w*w;
    case 31 : return v*v*w*w;
    case 32 : return u*u*v*w;
    case 33 : return u*v*v*w;
    case 34 : return u*v*w*w;
	default : throw "The shape functions have not bee coded for this order";
    }
}

double SHF::getTetdFctdu (int iFct, double u, double v, double w)
{
  switch (iFct)
    {
    case 0 : return 0.0;
    case 1 : return 1.0;
    case 2 : return 0.0;
    case 3 : return 0.0;
    case 4 : return 2*u;
    case 5 : return 0.0;
    case 6 : return 0.0;
    case 7 : return v;
    case 8 : return w;
    case 9 : return 0.0;
    case 10 : return 3.*u*u;
    case 11 : return 0.0;
    case 12 : return 0.0;
    case 13 : return 2.*u*v;
    case 14 : return 2.*u*w;
    case 15 : return v*v;
    case 16 : return 0.0;
    case 17 : return w*w;
    case 18 : return 0.0;
    case 19 : return v*w;
    case 20 : return 4*u*u*u;
    case 21 : return 0.0;
    case 22 : return 0.0;
    case 23 : return 3*u*u*v;
    case 24 : return 3*u*u*w;
    case 25 : return v*v*v;
    case 26 : return 0.0;
    case 27 : return w*w*w;
    case 28 : return 0.0;
    case 29 : return 2*u*v*v;
    case 30 : return 2*u*w*w;
    case 31 : return 0.0;
    case 32 : return 2*u*v*w;
    case 33 : return v*v*w;
    case 34 : return v*w*w;
	default : throw "The shape functions have not bee coded for this order";
    }
}

double SHF::getTetdFctdv (int iFct, double u, double v, double w)
{
  switch (iFct)
    {
    case 0 : return 0.0;
    case 1 : return 0.0;
    case 2 : return 1.0;
    case 3 : return 0.0;
    case 4 : return 0;
    case 5 : return 2.*v;
    case 6 : return 0.;
    case 7 : return u;
    case 8 : return 0;
    case 9 : return w;
    case 10 : return 0.0;
    case 11 : return 3.*v*v;
    case 12 : return 0.0;
    case 13 : return u*u;
    case 14 : return 0.0;
    case 15 : return 2.*v*u;
    case 16 : return 2.*v*w;
    case 17 : return 0.0;
    case 18 : return w*w;
    case 19 : return u*w;
    case 20 : return 0.0;
    case 21 : return 4*v*v*v;
    case 22 : return 0.0;
    case 23 : return u*u*u;
    case 24 : return 0.0;
    case 25 : return 3*v*v*u;
    case 26 : return 3*v*v*w;
    case 27 : return 0.0;
    case 28 : return w*w*w;
    case 29 : return u*u*2*v;
    case 30 : return 0.0;
    case 31 : return 2*v*w*w;
    case 32 : return u*u*w;
    case 33 : return u*2*v*w;
    case 34 : return u*w*w;
	default : throw "The shape functions have not bee coded for this order";
    }
}

double SHF::getTetdFctdw (int iFct, double u, double v, double w)
{
  switch (iFct)
    {
    case 0 : return 0.0;
    case 1 : return 0.0;
    case 2 : return 0.0;
    case 3 : return 1.0;
    case 4 : return 0;
    case 5 : return 0.0;
    case 6 : return 2.*w;
    case 7 : return 0;
    case 8 : return u;
    case 9 : return v;
    case 10 : return 0.0;
    case 11 : return 0.0;
    case 12 : return 3.*w*w;
    case 13 : return 0.0;
    case 14 : return u*u;
    case 15 : return 0.0;
    case 16 : return v*v;
    case 17 : return 2.*w*u;
    case 18 : return 2.*w*v;
    case 19 : return u*v;
    case 20 : return 0.0;
    case 21 : return 0.0;
    case 22 : return 4*w*w*w;
    case 23 : return 0.0;
    case 24 : return u*u*u;
    case 25 : return 0.0;
    case 26 : return v*v*v;
    case 27 : return 3*w*w*u;
    case 28 : return 3*w*w*v;
    case 29 : return 0.0;
    case 30 : return u*u*2*w;
    case 31 : return v*v*2*w;
    case 32 : return u*u*v;
    case 33 : return u*v*v;
    case 34 : return u*v*2*w;
	default : throw "The shape functions have not bee coded for this order";
    }
}

// Hexes

double SHF::getHexFct (int iFct, double u, double v, double w)
{
  switch(iFct)
    {
    case 0:return 1.0;
    case 1:return LegendrePolynomial(1,u);
    case 2:return LegendrePolynomial(1,v);
    case 3:return LegendrePolynomial(1,w);
    case 4:return LegendrePolynomial(1,u)*LegendrePolynomial(1,v);
    case 5:return LegendrePolynomial(1,u)*LegendrePolynomial(1,w);
    case 6:return LegendrePolynomial(1,v)*LegendrePolynomial(1,w);
    case 7:return LegendrePolynomial(1,u)*LegendrePolynomial(1,v)*LegendrePolynomial(1,w);
	default : throw "The shape functions have not bee coded for this order";
    }
}

double SHF::getHexdFctdu (int iFct, double u, double v, double w)
{
  switch(iFct)
    {
    case 0:return 0.0;
    case 1:return dLegendrePolynomial(1,u);
    case 2:return 0.0;
    case 3:return 0.0;
    case 4:return dLegendrePolynomial(1,u)*LegendrePolynomial(1,v);
    case 5:return dLegendrePolynomial(1,u)*LegendrePolynomial(1,w);
    case 6:return 0.0;
    case 7:return dLegendrePolynomial(1,u)*LegendrePolynomial(1,v)*LegendrePolynomial(1,w);
	default : throw "The shape functions have not bee coded for this order";
    }
}

double SHF::getHexdFctdv (int iFct, double u, double v, double w)
{
  switch(iFct)
    {
    case 0:return 0.0;
    case 1:return 0.0;
    case 2:return dLegendrePolynomial(1,v);
    case 3:return 0.0;
    case 4:return LegendrePolynomial(1,u)*dLegendrePolynomial(1,v);
    case 5:return 0.0;
    case 6:return dLegendrePolynomial(1,v)*LegendrePolynomial(1,w);
    case 7:return LegendrePolynomial(1,u)*dLegendrePolynomial(1,v)*LegendrePolynomial(1,w);
	default : throw "The shape functions have not bee coded for this order";
    }
}

double SHF::getHexdFctdw (int iFct, double u, double v, double w)
{
  switch(iFct)
    {
    case 0:return 0.0;
    case 1:return 0.0;
    case 2:return 0.0;
    case 3:return dLegendrePolynomial(1,w);
    case 4:return 0.0;
    case 5:return LegendrePolynomial(1,u)*dLegendrePolynomial(1,w);
    case 6:return LegendrePolynomial(1,v)*dLegendrePolynomial(1,w);
    case 7:return LegendrePolynomial(1,u)*LegendrePolynomial(1,v)*dLegendrePolynomial(1,w);
	default : throw "The shape functions have not bee coded for this order";
    }
}







