#include "BasisFunctions.h"
#include "mEntity.h"
#include "FunctionSpace.h"
#include "Integrator.h"

BasisFunctions:: BasisFunctions(mEntity* theMeshEntity, FunctionSpace *theFunctionSpace, int integOrder)
{
  GaussIntegrator gauss;
  int  nbPtGauss = gauss.nbIntegrationPoints(theMeshEntity,integOrder);
  int fSize = theFunctionSpace->size();
  int i;
  
  BF = new  double* [nbPtGauss] ;
  help = new double[nbPtGauss*fSize];
  
  for (i=0; i!= nbPtGauss; ++i)
    BF[i] = &help[i*fSize];
  double u,v,w,weight;
  for(i=0;i<nbPtGauss;++i)
    {
      gauss.iPoint(theMeshEntity,i,integOrder,u,v,w,weight);
      theFunctionSpace->fcts(u,v,w,BF[i]);
    }
}

double* BasisFunctions::ifct(int i) const
{return BF[i];}
