#include "Integrator.h"
#include "mEntity.h"
#include "IntPt.h"
#include "BasisFunctions.h"
#include "FunctionSpace.h"
#include "Vertex.h"
#include "Constants.h"

int getNGQ1DPts(int);
int getNGQTPts(int);
int getNGQQPts(int);
int getNGQTetPts(int);
int getNGQHPts(int);
IntPt2d *getGQTPts(int order);
IntPt2d *getGQQPts(int order);
IntPt3d *getGQTetPts(int order);
IntPt3d *getGQHPts(int order);
int GaussLegendre1D(int,double **, double **);

GaussIntegrator::GaussIntegrator()
{
  Edge0=0;Edge1=0;Edge2=0;Edge3=0;Edge4=0;Edge5=0;Edge6=0;Edge7=0;Edge8=0;Edge9=0;Edge10=0;Edge11=0;Edge12=0;Edge13=0;Edge14=0;Edge15=0;Edge16=0;Edge17=0;Edge18=0;Edge19=0;Edge20=0;Edge21=0;Edge22=0;Edge23=0;Edge24=0;Edge25=0;
  Tri0=0; Tri1=0; Tri2=0;Tri3=0;Tri4=0;Tri5=0;Tri6=0;Tri7=0;Tri8=0,Tri9=0;Tri10=0;Tri11=0;Tri12=0;Tri13=0;Tri14=0;Tri15=0;Tri16=0;Tri17=0;Tri18=0,Tri19=0;Tri20=0; Tri21=0; Tri22=0;Tri23=0;Tri24=0;Tri25=0;
  Quad0=0;Quad1=0;Quad2=0;Quad3=0;Quad4=0;Quad5=0;Quad6=0;Quad7=0;Quad8=0;
  Hex0=0;Hex1=0;Hex2=0;Hex3=0;Hex4=0;Hex5=0;Hex6=0;Hex7=0;
  Tet0=0;Tet1=0;Tet2=0;Tet3=0;Tet4=0;Tet5=0;Tet6=0;Tet7=0;Tet8=0;Tet9=0;Tet10=0;Tet11=0;Tet12=0;Tet13=0;Tet14=0;
  fctandDerAllEdgesQuad = 0;
  fctandDerAllEdgesTri = 0;
}

int GaussIntegrator::nbIntegrationPoints(mEntity *ent,int order) const
{
  switch(ent->getType())
    {
    case mEntity::VERTEX : return 1;
    case mEntity::EDGE   : return (order<1)?1:getNGQ1DPts(order); 
    case mEntity::TRI    : return (order<1)?1:getNGQTPts(order);
    case mEntity::TET    : return (order<1)?1:getNGQTetPts(order);
    case mEntity::QUAD   : return getNGQQPts(order);
    case mEntity::HEX    : return getNGQHPts(order);
	default: printf("Not coded  yet\n"); 
		return 0;
    }
}

void GaussIntegrator::iPoint(mEntity *ent,int i, int order , double &u, double &v, double &w, double &weight)const
{
  switch(ent->getType())
    {
    case mEntity::VERTEX :
      u = v = w = 0.0;
      weight = 1.0;
      break;
    case mEntity::EDGE :
      {
	double *pt,*wt;
	GaussLegendre1D(nbIntegrationPoints(ent,order),&pt,&wt);
	u = pt[i];
	v = w = 0.0;
	weight = wt[i];
      }
      break;
    case mEntity::TRI  :
      {
	IntPt2d *pts = getGQTPts(order);
	u = pts[i].pt[0];
	v = pts[i].pt[1];
	w = 0.0;
	weight = 0.5 * pts[i].weight;
      }
      break;
    case mEntity::QUAD :
      {
	IntPt2d *pts = getGQQPts(order);
	u = pts[i].pt[0];
	v = pts[i].pt[1];
	w = 0.0;
	weight = pts[i].weight;
      }
      break;
    case mEntity::HEX :
      {
	IntPt3d *pts = getGQHPts(order);
	u = pts[i].pt[0];
	v = pts[i].pt[1];
	w = pts[i].pt[2];
	weight = pts[i].weight;
      }
      break;
    case mEntity::TET :
      {
	const double sixth = 1.0/6.0;
	IntPt3d *pts = getGQTetPts(order);
	u = pts[i].pt[0];
	v = pts[i].pt[1];
	w = pts[i].pt[2];
	weight = sixth*pts[i].weight;
      }
      break;
    }
}

double*  GaussIntegrator::iFct(mEntity *ent,FunctionSpace *FS,int i,int order) 
{
  switch(ent->getType())
    {
    case mEntity::VERTEX :
      printf("Error! Integration on a point is underfined!\n");
      return 0;
      break;
	case mEntity::EDGE :
      {
	if (order == 0) if (Edge0) return Edge0->ifct(i);
	else {Edge0 = new BasisFunctions(ent,FS,order); return Edge0->ifct(i);}
	if (order == 1) if (Edge1) return Edge1->ifct(i);
	else {Edge1 = new BasisFunctions(ent,FS,order); return Edge1->ifct(i);}
	if (order == 2) if (Edge2) return Edge2->ifct(i);
	else {Edge2 = new BasisFunctions(ent,FS,order); return Edge2->ifct(i);}
	if (order == 3) if (Edge3) return Edge3->ifct(i);
	else {Edge3 = new BasisFunctions(ent,FS,order); return Edge3->ifct(i);}
	if (order == 4) if (Edge4) return Edge4->ifct(i);
	else {Edge4 = new BasisFunctions(ent,FS,order); return Edge4->ifct(i);}
	if (order == 5) if (Edge5) return Edge5->ifct(i);
	else {Edge5 = new BasisFunctions(ent,FS,order); return Edge5->ifct(i);}
	if (order == 6) if (Edge6) return Edge6->ifct(i);
	else {Edge6 = new BasisFunctions(ent,FS,order); return Edge6->ifct(i);}
	if (order == 7) if (Edge7) return Edge7->ifct(i);
	else {Edge7 = new BasisFunctions(ent,FS,order); return Edge7->ifct(i);}
	if (order == 8) if (Edge8) return Edge8->ifct(i);
	else {Edge8 = new BasisFunctions(ent,FS,order); return Edge8->ifct(i);}
	if (order == 9) if (Edge9) return Edge9->ifct(i);
	else {Edge9 = new BasisFunctions(ent,FS,order); return Edge9->ifct(i);}
	if (order == 10) if (Edge10) return Edge10->ifct(i);
	else {Edge10 = new BasisFunctions(ent,FS,order); return Edge10->ifct(i);}
	if (order == 11) if (Edge11) return Edge11->ifct(i);							//newly added
	else {Edge11 = new BasisFunctions(ent,FS,order); return Edge11->ifct(i);}
	if (order == 12) if (Edge12) return Edge12->ifct(i);
	else {Edge12 = new BasisFunctions(ent,FS,order); return Edge12->ifct(i);}
	if (order == 13) if (Edge13) return Edge13->ifct(i);
	else {Edge13 = new BasisFunctions(ent,FS,order); return Edge13->ifct(i);}
	if (order == 14) if (Edge14) return Edge14->ifct(i);
	else {Edge14 = new BasisFunctions(ent,FS,order); return Edge14->ifct(i);}
	if (order == 15) if (Edge15) return Edge15->ifct(i);
	else {Edge15 = new BasisFunctions(ent,FS,order); return Edge15->ifct(i);}
	if (order == 16) if (Edge16) return Edge16->ifct(i);
	else {Edge16 = new BasisFunctions(ent,FS,order); return Edge16->ifct(i);}
	if (order == 17) if (Edge17) return Edge17->ifct(i);
	else {Edge17 = new BasisFunctions(ent,FS,order); return Edge17->ifct(i);}
	if (order == 18) if (Edge18) return Edge18->ifct(i);
	else {Edge18 = new BasisFunctions(ent,FS,order); return Edge18->ifct(i);}
	if (order == 19) if (Edge19) return Edge19->ifct(i);
	else {Edge19 = new BasisFunctions(ent,FS,order); return Edge19->ifct(i);}
	if (order == 20) if (Edge20) return Edge20->ifct(i);
	else {Edge20 = new BasisFunctions(ent,FS,order); return Edge20->ifct(i);}
	if (order == 21) if (Edge21) return Edge21->ifct(i);
	else {Edge21 = new BasisFunctions(ent,FS,order); return Edge21->ifct(i);}
	if (order == 22) if (Edge22) return Edge22->ifct(i);
	else {Edge22 = new BasisFunctions(ent,FS,order); return Edge22->ifct(i);}
	if (order == 23) if (Edge23) return Edge23->ifct(i);
	else {Edge23 = new BasisFunctions(ent,FS,order); return Edge23->ifct(i);}
	if (order == 24) if (Edge24) return Edge24->ifct(i);
	else {Edge24 = new BasisFunctions(ent,FS,order); return Edge24->ifct(i);}
	if (order == 25) if (Edge25) return Edge25->ifct(i);
	else {Edge25 = new BasisFunctions(ent,FS,order); return Edge25->ifct(i);}
      }
    break;
    case mEntity::TRI  :
      {
	if (order == 0) if (Tri0) return Tri0->ifct(i);
	else {Tri0 = new BasisFunctions(ent,FS,order); return Tri0->ifct(i);}
	if (order == 1) if (Tri1) return Tri1->ifct(i);
	else {Tri1 = new BasisFunctions(ent,FS,order); return Tri1->ifct(i);}
	if (order == 2) if (Tri2) return Tri2->ifct(i);
	else {Tri2 = new BasisFunctions(ent,FS,order); return Tri2->ifct(i);}
	if (order == 3) if (Tri3) return Tri3->ifct(i);
	else {Tri3 = new BasisFunctions(ent,FS,order); return Tri3->ifct(i);}
	if (order == 4) if (Tri4) return Tri4->ifct(i);
	else {Tri4 = new BasisFunctions(ent,FS,order); return Tri4->ifct(i);}
	if (order == 5) if (Tri5) return Tri5->ifct(i);
	else {Tri5 = new BasisFunctions(ent,FS,order); return Tri5->ifct(i);}
	if (order == 6) if (Tri6) return Tri6->ifct(i);
	else {Tri6 = new BasisFunctions(ent,FS,order); return Tri6->ifct(i);}
	if (order == 7) if (Tri7) return Tri7->ifct(i);
	else {Tri7 = new BasisFunctions(ent,FS,order); return Tri7->ifct(i);}
	if (order == 8) if (Tri8) return Tri8->ifct(i);
	else {Tri8 = new BasisFunctions(ent,FS,order); return Tri8->ifct(i);}
	if (order == 9) if (Tri9) return Tri9->ifct(i);
	else {Tri9 = new BasisFunctions(ent,FS,order); return Tri9->ifct(i);}
	if (order == 10) if (Tri10) return Tri10->ifct(i);
	else {Tri10 = new BasisFunctions(ent,FS,order); return Tri10->ifct(i);}
	if (order == 11) if (Tri11) return Tri11->ifct(i);                              //newly added
	else {Tri11 = new BasisFunctions(ent,FS,order); return Tri11->ifct(i);}
	if (order == 12) if (Tri12) return Tri12->ifct(i);
	else {Tri12 = new BasisFunctions(ent,FS,order); return Tri12->ifct(i);}
	if (order == 13) if (Tri13) return Tri13->ifct(i);
	else {Tri13 = new BasisFunctions(ent,FS,order); return Tri13->ifct(i);}
	if (order == 14) if (Tri14) return Tri14->ifct(i);
	else {Tri14 = new BasisFunctions(ent,FS,order); return Tri14->ifct(i);}
	if (order == 15) if (Tri15) return Tri15->ifct(i);
	else {Tri15 = new BasisFunctions(ent,FS,order); return Tri15->ifct(i);}
	if (order == 16) if (Tri16) return Tri16->ifct(i);
	else {Tri16 = new BasisFunctions(ent,FS,order); return Tri16->ifct(i);}
	if (order == 17) if (Tri17) return Tri17->ifct(i);
	else {Tri17 = new BasisFunctions(ent,FS,order); return Tri17->ifct(i);}
	if (order == 18) if (Tri18) return Tri18->ifct(i);
	else {Tri18 = new BasisFunctions(ent,FS,order); return Tri18->ifct(i);}
	if (order == 19) if (Tri19) return Tri19->ifct(i);
	else {Tri19 = new BasisFunctions(ent,FS,order); return Tri19->ifct(i);}
	if (order == 20) if (Tri20) return Tri20->ifct(i);
	else {Tri20 = new BasisFunctions(ent,FS,order); return Tri20->ifct(i);}
	if (order == 21) if (Tri21) return Tri21->ifct(i);
	else {Tri21 = new BasisFunctions(ent,FS,order); return Tri21->ifct(i);}
	if (order == 22) if (Tri22) return Tri22->ifct(i);
	else {Tri22 = new BasisFunctions(ent,FS,order); return Tri22->ifct(i);}
	if (order == 23) if (Tri23) return Tri23->ifct(i);
	else {Tri23 = new BasisFunctions(ent,FS,order); return Tri23->ifct(i);}
	if (order == 24) if (Tri24) return Tri24->ifct(i);
	else {Tri24 = new BasisFunctions(ent,FS,order); return Tri24->ifct(i);}
	if (order == 25) if (Tri25) return Tri25->ifct(i);
	else {Tri25 = new BasisFunctions(ent,FS,order); return Tri25->ifct(i);}
      }
    break;
    case mEntity::QUAD :
      {
	if (order == 0) if (Quad0) return Quad0->ifct(i);
	else {Quad0 = new BasisFunctions(ent,FS,order); return Quad0->ifct(i);}
	if (order == 1) if (Quad1) return Quad1->ifct(i);
	else {Quad1 = new BasisFunctions(ent,FS,order); return Quad1->ifct(i);}
	if (order == 2) if (Quad2) return Quad2->ifct(i);
	else {Quad2 = new BasisFunctions(ent,FS,order); return Quad2->ifct(i);}
	if (order == 3) if (Quad3) return Quad3->ifct(i);
	else {Quad3 = new BasisFunctions(ent,FS,order); return Quad3->ifct(i);}
	if (order == 4) if (Quad4) return Quad4->ifct(i);
	else {Quad4 = new BasisFunctions(ent,FS,order); return Quad4->ifct(i);}
	if (order == 5) if (Quad5) return Quad5->ifct(i);
	else {Quad5 = new BasisFunctions(ent,FS,order); return Quad5->ifct(i);}
	if (order == 6) if (Quad6) return Quad6->ifct(i);
	else {Quad6 = new BasisFunctions(ent,FS,order); return Quad6->ifct(i);}
	if (order == 7) if (Quad7) return Quad7->ifct(i);
	else {Quad7 = new BasisFunctions(ent,FS,order); return Quad7->ifct(i);}
	if (order == 8) if (Quad8) return Quad8->ifct(i);
	else {Quad8 = new BasisFunctions(ent,FS,order); return Quad8->ifct(i);}
      }
      break;
    case mEntity::HEX :
      {
	printf("not coded yet\n");
	return 0;
      }
      break;
    case mEntity::TET :
      {

	if (order == 0) if (Tet0) return Tet0->ifct(i);
	else {Tet0 = new BasisFunctions(ent,FS,order); return Tet0->ifct(i);}
	if (order == 1) if (Tet1) return Tet1->ifct(i);
	else {Tet1 = new BasisFunctions(ent,FS,order); return Tet1->ifct(i);}
	if (order == 2) if (Tet2) return Tet2->ifct(i);
	else {Tet2 = new BasisFunctions(ent,FS,order); return Tet2->ifct(i);}
	if (order == 3) if (Tet3) return Tet3->ifct(i);
	else {Tet3 = new BasisFunctions(ent,FS,order); return Tet3->ifct(i);}
	if (order == 4) if (Tet4) return Tet4->ifct(i);
	else {Tet4 = new BasisFunctions(ent,FS,order); return Tet4->ifct(i);}
	if (order == 5) if (Tet5) return Tet5->ifct(i);
	else {Tet5 = new BasisFunctions(ent,FS,order); return Tet5->ifct(i);}
	if (order == 6) if (Tet6) return Tet6->ifct(i);
	else {Tet6 = new BasisFunctions(ent,FS,order); return Tet6->ifct(i);}
	if (order == 7) if (Tet7) return Tet7->ifct(i);
	else {Tet7 = new BasisFunctions(ent,FS,order); return Tet7->ifct(i);}
	if (order == 8) if (Tet8) return Tet8->ifct(i);
	else {Tet8 = new BasisFunctions(ent,FS,order); return Tet8->ifct(i);}
	if (order == 9) if (Tet9) return Tet9->ifct(i);
	else {Tet9 = new BasisFunctions(ent,FS,order); return Tet9->ifct(i);}
	if (order == 10) if (Tet10) return Tet10->ifct(i);
	else {Tet10 = new BasisFunctions(ent,FS,order); return Tet10->ifct(i);}
	if (order == 11) if (Tet11) return Tet11->ifct(i);
	else {Tet11 = new BasisFunctions(ent,FS,order); return Tet11->ifct(i);}
	if (order == 12) if (Tet12) return Tet12->ifct(i);
	else {Tet12 = new BasisFunctions(ent,FS,order); return Tet12->ifct(i);}
	if (order == 13) if (Tet13) return Tet13->ifct(i);
	else {Tet13 = new BasisFunctions(ent,FS,order); return Tet13->ifct(i);}
	if (order == 14) if (Tet14) return Tet14->ifct(i);
	else {Tet14 = new BasisFunctions(ent,FS,order); return Tet14->ifct(i);}
	return 0;
      }
      break;
	default: 
		printf("not coded yet\n");
		return 0;
    }
return 0;
}


double*  GaussIntegrator::iBoundFct(int i, int integOrder, mEntity *boundary,mEntity *cell) 
{
  int nbEdges;
  FctandDerAllEdges* fctandDerAllEdges=0;
  switch(boundary->getType())						//needs higher order code?  //GRANT
    {
    case mEntity::VERTEX :
      printf("Error! Integration on a point is underfined!\n");
      return 0;
      break;
    case mEntity::EDGE :
      Vertex* v1 = (Vertex* ) boundary->get(0,0);
      Vertex* v2 = (Vertex* ) boundary->get(0,1);
	  //v1->print(); v2->print();
      
	  	  if (cell->getType()==mEntity::TRI) 
	  {
		  nbEdges=3;
		  if (!fctandDerAllEdgesTri) fctandDerAllEdgesTri=new FctandDerAllEdgesTri(boundary,cell);
		  fctandDerAllEdges = fctandDerAllEdgesTri;
	  }

      if (cell->getType()==mEntity::QUAD) 
	  {
		  nbEdges=4;
          if (!fctandDerAllEdgesQuad) fctandDerAllEdgesQuad=new FctandDerAllEdgesQuad(boundary,cell);
		  fctandDerAllEdges = fctandDerAllEdgesQuad;
      }
     
	  int nbPoints = nbIntegrationPoints(boundary,integOrder);
	  for (int j=0; j<nbEdges;j++)
	{
	  Vertex* V1 = (Vertex* ) cell->get(0,j);
	  Vertex* V2 = (Vertex* ) cell->get(0,(j+1)%nbEdges);
	  if (v1==V1 && v2==V2)  //edge 1
	    return (*fctandDerAllEdges)(j)->getEdgeFcts(i,integOrder);
     if (v1==V2 && v2==V1)  //edge 1
	    return (*fctandDerAllEdges)(j)->getEdgeFcts(nbPoints-i-1,integOrder);
	}
    }
  return 0;
}

#define MAX_INTEG_ORDER 9

FctandDerEdge::FctandDerEdge(int Nb, mEntity *edge, mEntity *cell)	//higher order code? dont think so...
  :edgeNb(Nb)
{
  double u,v,w,weight;
  int maxNbPoints = MAX_INTEG_ORDER/2 +1;
  fcts_and_der.reserve(maxNbPoints);
  
  for (int i=1;i<=maxNbPoints;i++) {
  GaussIntegrator gauss;
  int integOrder = 2*i-1;
  int  nbPtGauss = gauss.nbIntegrationPoints(edge,integOrder);

  double** functions;
  functions = new  double* [nbPtGauss] ;
  double* help = new double[nbPtGauss*MaxNbFunctions];
  for (int j=0; j < nbPtGauss; ++j)
    functions[j] = &help[j*MaxNbFunctions];
  FunctionSpace *fs;
  
  if (cell->getType()==mEntity::QUAD) {
	  fs = new QuadFunctionSpace(MAX_APPROX_ORDER_QUAD_TENSOR);
  for(i=0;i<nbPtGauss;++i)
    {
      if (edgeNb >1) gauss.iPoint(edge,i,integOrder,u,v,w,weight);
	  else gauss.iPoint(edge,nbPtGauss-i-1,integOrder,u,v,w,weight);
      if (edgeNb ==0) v = -1.;
      if (edgeNb ==1) {v = u; u = 1.;}
      if (edgeNb ==2) v = 1.;
      if (edgeNb ==3) {v=u; u =-1.;}	  
      fs->fcts(u,v,w,functions[i]);
	 /* printf("edge %d order %d point %d\n",edgeNb,integOrder, i);
	  for (int k=0; k<MaxNbFunctions; k++)
		  printf("%f ", functions[i][k]);
	  printf("\n");*/
    }
  }

  if (cell->getType()==mEntity::TRI) {
	  if (cell->getgeomOrder() == 1)
		fs= new OrthogonalTriangleFunctionSpace(MAX_APPROX_ORDER_QUAD_TENSOR);
	  else
	    fs= new TriangleFunctionSpace(MAX_APPROX_ORDER_QUAD_TENSOR);
	  for(i=0;i<nbPtGauss;++i)
    {
      gauss.iPoint(edge,i,integOrder,u,v,w,weight);
	  if (edgeNb ==0) {u=0.5*(1.-u); v = 0.;}
      if (edgeNb ==1) {w = u;
	                   v = 0.5*(1.-w); u = 0.5*(1.+w);}
      if (edgeNb ==2) {v = 0.5*(1.+u); u =0.;}
      fs->fcts(u,v,0,functions[i]);
/*	  printf("edge %d order %d point %d\n",edgeNb,integOrder, i);
	  for (int k=0; k<MaxNbFunctions; k++)
		  printf("%f ", functions[i][k]);
	  printf("\n");*/
    }
  }

    FctandDer* fd = new FctandDer(functions);
    fcts_and_der.push_back(fd);
  }
}

FctandDerAllEdges::FctandDerAllEdges(mEntity *edge, mEntity *cell)
{
	int size = cell->size(0);
  edges.reserve(size);
  for (int i=0; i<size;i++)
    {
      FctandDerEdge* e = new FctandDerEdge(i,edge,cell);
      edges.push_back(e);
    }
};
