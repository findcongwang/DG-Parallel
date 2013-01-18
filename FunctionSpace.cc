#include "Mapping.h"
#include "FunctionSpace.h"
#include "mTensor2.h"
#include <math.h>
#include <stdio.h>

double FunctionSpace::grads(double u, double v, double w, 
			  Mapping *m, vector<mVector> &gr)const 
{
  mTensor2 jInv;
  gradiends(u,v,w,gr);
  double detJac = m->jacInverse(u,v,w,jInv);
  for(int i=0;i<Size;i++)
      gr[i] *= jInv;
  return detJac;
}

void OrthogonalTriangleFunctionSpace::gradiends (double u, double v, double w, 
					     vector<mVector> &gr)const
{
  for(int i=0;i<Size;i++)
    {
      gr[i](0) = SHF::PascalgetdFctduOrthogonal(i,u,v);
      gr[i](1) = SHF::PascalgetdFctdvOrthogonal(i,u,v);
      gr[i](2) = 0.0;
    }
}

void QuadFunctionSpace::gradiends (double u, double v, double w, 
			       vector<mVector> &gr)const
{
  for(int i=0;i<Size;i++)
    {
      gr[i](0) = SHF::getQdFctdu(i,u,v,0);
      gr[i](1) = SHF::getQdFctdv(i,u,v,0);
      gr[i](2) = 0.0;
    }
}

void HexFunctionSpace::gradiends (double u, double v, double w, 
			      vector<mVector> &gr)const
{
  for(int i=0;i<Size;i++)
    {
      gr[i](0) = SHF::getHexdFctdu(i,u,v,w);
      gr[i](1) = SHF::getHexdFctdv(i,u,v,w);
      gr[i](2) = SHF::getHexdFctdw(i,u,v,w);
    }
}

void OrthogonalTetFunctionSpace::gradiends (double u, double v, double w, 
			       vector<mVector> &gr)const
{
  for(int i=0;i<Size;i++)
    {
      gr[i](0) = SHF::getTetdFctduOrthogonal(i,u,v,w);
      gr[i](1) = SHF::getTetdFctdvOrthogonal(i,u,v,w);
      gr[i](2) = SHF::getTetdFctdwOrthogonal(i,u,v,w);
    }
}

double* OrthogonalTetFunctionSpace::fcts(int order, int ithface, int ithpoint) const 
{
	switch(order)
	{
		case 0: return &TetFaceP0[ithface*1+ithpoint][0]; 
		break;
		case 1: return &TetFaceP1[ithface*4+ithpoint][0]; 
		break;
		case 2: return &TetFaceP2[ithface*7+ithpoint][0]; 
		break;
		case 3: return &(TetFaceP3[ithface*13+ithpoint][0]); 
		break;
		default: return 0;
	}
}