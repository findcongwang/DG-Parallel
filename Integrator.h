#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_
#include <stdio.h>
#include <vector>
class mEntity;
class BasisFunctions;
class FunctionSpace;

using namespace std;


class FctandDer
{
 protected:
 double** functions;
  double** derivatives;
 
 public:
  FctandDer(double **f):functions(f){}
  double* getFct(int i) const { return functions[i];}
  double*  getDer(int i) const {return derivatives[i];}
};



class FctandDerEdge
{
 protected:
  int edgeNb;
  vector<FctandDer*> fcts_and_der;
 public:
  FctandDerEdge(int edgeNb, mEntity *edge, mEntity *cell); 
  double* getEdgeFcts(int i, int integOrder) {return fcts_and_der[integOrder/2]->getFct(i);}
};

class FctandDerAllEdges
{
 protected:
  vector<FctandDerEdge*> edges;
 public:
  FctandDerAllEdges(mEntity*, mEntity*);
  FctandDerEdge* operator()(int i){return edges[i];}
};

class FctandDerAllEdgesQuad : public FctandDerAllEdges
{
public:
	FctandDerAllEdgesQuad(mEntity *edge, mEntity *cell) : FctandDerAllEdges(edge, cell) {}
};

class FctandDerAllEdgesTri : public FctandDerAllEdges
{
public:
	FctandDerAllEdgesTri(mEntity *edge, mEntity *cell) : FctandDerAllEdges(edge, cell) {}
};

class Intergrator
{
  public :
    virtual int nbIntegrationPoints(mEntity *ent, int order) const = 0;
  virtual void iPoint(mEntity *ent,int i, int order , double &u, double &v, double &w, double &weight) const = 0;  
};

class GaussIntegrator
{
  BasisFunctions *Edge0,*Edge1,*Edge2,*Edge3,*Edge4,*Edge5,*Edge6,*Edge7,*Edge8,*Edge9,*Edge10,*Edge11,*Edge12,*Edge13,*Edge14,*Edge15,*Edge16,*Edge17,*Edge18,*Edge19,*Edge20,*Edge21,*Edge22,*Edge23,*Edge24,*Edge25;
  BasisFunctions *Tri0,*Tri1,*Tri2,*Tri3,*Tri4,*Tri5,*Tri6,*Tri7,*Tri8,*Tri9,*Tri10,*Tri11,*Tri12,*Tri13,*Tri14,*Tri15,*Tri16,*Tri17,*Tri18,*Tri19,*Tri20,*Tri21,*Tri22,*Tri23,*Tri24,*Tri25;
  BasisFunctions *Quad0,*Quad1,*Quad2,*Quad3,*Quad4,*Quad5,*Quad6,*Quad7,*Quad8;
  BasisFunctions *Hex0,*Hex1,*Hex2,*Hex3,*Hex4,*Hex5,*Hex6,*Hex7;
  BasisFunctions *Tet0,*Tet1,*Tet2,*Tet3,*Tet4,*Tet5,*Tet6,*Tet7,*Tet8,*Tet9,*Tet10,*Tet11,*Tet12,*Tet13,*Tet14;
  FctandDerAllEdgesQuad* fctandDerAllEdgesQuad;
  FctandDerAllEdgesTri* fctandDerAllEdgesTri;

 public:
  GaussIntegrator();
  virtual int nbIntegrationPoints(mEntity *,int order) const;
  virtual void iPoint(mEntity *,int i,int order,double &u, double &v, double &w, double &weight) const;  
  double* iFct(mEntity *,FunctionSpace*,int ,int );
  double* iBoundFct(int i, int order, mEntity *boundary,mEntity *cell);
  double* iEdgeFct(int i, int order, mEntity* edge, mEntity* cell);
};

#endif
