#ifndef H_MeshMapping
#define H_MeshMapping
class mEntity;
#include "mTensor2.h"
class mVector;
#include <vector>
#include "mCompiler.h"
#include "Integrator.h"  //for arclength integration

class Mapping
{
 public:
  virtual ~Mapping(){};
  virtual void eval(double u, double v, double w,
		    double &x, double &y, double &z) const = 0;
  virtual bool invert(double x, double y, double z,
		      double &u, double &v, double &w) const = 0;
  virtual double jacInverse(double x, double y, double z, 
			    mTensor2&) const = 0;
  virtual double detJac(double u, double v, double w) const = 0;
  virtual void normalVector(mEntity *, double, double, double, mVector &) const = 0;
  virtual int detJacOrder() const = 0;
  virtual int geomOrder() const = 0;
  virtual bool inReferenceElement(double u, double v, double w) const = 0;
  virtual void COG (double &u, double &v, double &w) const = 0;
};


class MeshMapping : public Mapping
{
 public:
  MeshMapping(mEntity * me);
  virtual ~MeshMapping(){};
  virtual void eval(double u, double v, double w,
		    double &x, double &y, double &z) const;
  virtual bool invert(double x, double y, double z,
		      double &u, double &v, double &w) const;
  virtual double jacInverse(double x, double y, double z, 
				 mTensor2&) const;
  virtual double detJac(double u, double v, double w) const;
  virtual void normalVector(mEntity *, double, double, double, mVector &) const;
  virtual bool inReferenceElement(double u, double v, double w) const;  
  virtual int detJacOrder() const;
  virtual int geomOrder() const;
  virtual void COG (double &u, double &v, double &w) const;
  double GeomShapeFunction         (int iNod, double u, double v, double w ) const;
  void   GradGeomShapeFunction     (int iNod, double u, double v ,double w , mVector&) const;
  double GeomShapeFunctionTri      (int iNod, double u, double v, double w ) const;
  void   GradGeomShapeFunctionTri  (int iNod, double u, double v ,double w , mVector&) const;
  double GeomShapeFunctionQuad     (int iNod, double u, double v, double w ) const;
  void   GradGeomShapeFunctionQuad (int iNod, double u, double v ,double w , mVector&) const;
  double GeomShapeFunctionLine     (int iNod, double u, double v, double w ) const;
  void   GradGeomShapeFunctionLine (int iNod, double u, double v ,double w , mVector&) const;
  double GeomShapeFunctionTet     (int iNod, double u, double v, double w ) const;
  void   GradGeomShapeFunctionTet (int iNod, double u, double v ,double w , mVector&) const;
  double GeomShapeFunctionHex     (int iNod, double u, double v, double w ) const;
  void   GradGeomShapeFunctionHex (int iNod, double u, double v ,double w , mVector&) const;

protected:
  mEntity *ent;

};

class ConstantMeshMapping : public MeshMapping
{
  double det;
  mTensor2 jac;
 public:
  ConstantMeshMapping(mEntity * me);
  virtual ~ConstantMeshMapping(){};
  virtual double jacInverse(double x, double y, double z, 
				 mTensor2&) const;
  virtual double detJac(double u, double v, double w) const;
};

class CylindricalCoordinatesMeshMapping : public MeshMapping
{
 public:
  virtual ~CylindricalCoordinatesMeshMapping(){}
  CylindricalCoordinatesMeshMapping(mEntity * me);
  virtual double jacInverse(double x, double y, double z, 
				 mTensor2&) const;
  virtual double detJac(double u, double v, double w) const;
  virtual int detJacOrder() const;
};

#endif

