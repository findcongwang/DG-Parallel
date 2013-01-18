#ifndef _DGINTERPOLATION_H_
#define _DGINTERPOLATION_H_

#include <vector>
#include "mVector.h"
#include "ShapeFunctionList.h"
#include "BFunctions.h"
using namespace std;
class Mapping;

class FunctionSpace
{
 protected:
  short p;
  short Size;
 public:
  FunctionSpace(short order=1) : p(order){}
  short size () const {return Size;}
  virtual short size (short) const = 0;
  //returns Size(p) - Size(p-1) - true only for triangles!!!!;
  short pSize () const {return p+1;}
  short order() const {return p;}
  void setOrder(short order) { p=order; resetSize(); }
  virtual void resetSize() = 0;
  virtual void fcts  (double u, double v, double w, double *f) const = 0;
   virtual double* fcts  (int order, int q, int ith) const = 0;
  double grads (double u, double v, double w, 
		Mapping *, vector<mVector> &gr) const;
  virtual void  gradiends (double u, double v, double w, 
	        vector<mVector> &gr) const =0;
  virtual bool isOrthogonal() const = 0;
};

class OrthogonalTriangleFunctionSpace : public FunctionSpace
{
  public :
    OrthogonalTriangleFunctionSpace(short order) : FunctionSpace(order)
    {  Size = (short)(p+1)*(p+2)/2;  }
  virtual void resetSize(){Size = (short)(p+1)*(p+2)/2;}
  virtual short size (short q) const {return (short)(q+1)*(q+2)/2;}
  virtual void fcts  (double u, double v, double w, double *f) const
    { for(short i=0;i<Size;++i) f[i] = SHF::PascalgetFctOrthogonal(i,u,v); }
   virtual double* fcts  (int order, int q, int ith) const {return 0;}
  virtual void gradiends (double u, double v, double w, vector<mVector> &gr) const;
  virtual bool isOrthogonal() const {return 1;}  
};

class TriangleFunctionSpace : public OrthogonalTriangleFunctionSpace
{
  public :
	TriangleFunctionSpace(short order) : OrthogonalTriangleFunctionSpace(order)
    {  Size = (short)(p+1)*(p+2)/2;  }
  virtual bool isOrthogonal() const {return 0;}  
};

class QuadFunctionSpace : public FunctionSpace
{
  public :
    QuadFunctionSpace(short order) : FunctionSpace(order)
    {  Size = (p+1)*(p+1);  }
  virtual void resetSize(){Size =(p+1)*(p+1); }
  virtual short size (short q) const {return (q+1)*(q+1);}
  virtual void fcts  (double u, double v, double w, double *f) const
    { for(int i=0;i<Size;i++)f[i] = SHF::getQFct(i,u,v,0); }
  virtual double* fcts  (int order, int q, int ith) const {return 0;}
  virtual void gradiends (double u, double v, double w, vector<mVector> &gr) const;
  virtual bool isOrthogonal() const {return 1;}
};

class HexFunctionSpace : public FunctionSpace
{
  public :
    HexFunctionSpace(short order) : FunctionSpace(order)
    { Size = (p+1)*(p+1)*(p+1); }
  virtual void resetSize(){Size =(p+1)*(p+1)*(p+1); }
  virtual short size (short q) const {return (q+1)*(q+1)*(q+1);}
  virtual void fcts  (double u, double v, double w, double *f) const
    { for(int i=0;i<Size;i++)f[i] = SHF::getHexFct(i,u,v,w);  }
  virtual double* fcts  (int order, int q, int ith) const {return 0;}
  virtual void gradiends (double u, double v, double w, vector<mVector> &gr) const;
  virtual bool isOrthogonal() const {return 0;} 
};

class OrthogonalTetFunctionSpace : public FunctionSpace
{
  public :
    OrthogonalTetFunctionSpace(short order) : FunctionSpace(order)
    { Size = (p+1)*(p+2)*(p+3)/6; }
  virtual void resetSize(){Size = (p+1)*(p+1)*(p+1);}
  virtual short size (short q) const {return(q+1)*(q+2)*(q+3)/6 ;}
  virtual void fcts  (double u, double v, double w, double *f) const
    {  for(int i=0;i<Size;i++)f[i] = SHF::getTetFctOrthogonal(i,u,v,w); }
  virtual double* fcts  (int order, int q, int ith) const;
  virtual void gradiends (double u, double v, double w, vector<mVector> &gr) const;
  virtual bool isOrthogonal() const {return 1;} 
};

class TetFunctionSpace : public OrthogonalTetFunctionSpace
{
  public:
	TetFunctionSpace(short order) : OrthogonalTetFunctionSpace(order)
	{ Size = (p+1)*(p+2)*(p+3)/6; }
	virtual bool isOrthogonal() const {return 0;}
};
#endif
