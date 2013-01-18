#ifndef _BASISFUNCTIONS_H_
#define _BASISFUNCTIONS_H_

class FunctionSpace;
class mEntity;

class BasisFunctions
{ 
 protected:
  double **BF;
  double **dBF;
  double* help;
 public:
  BasisFunctions(mEntity* theMeshEntity, FunctionSpace *theFunctionSpace,int);
  BasisFunctions(){};
  double* ifct(int) const;
 ~BasisFunctions(){};
};

#endif
