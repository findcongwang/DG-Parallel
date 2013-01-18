// Face.cpp: implementation of the Face class.
#include "Face.h"
#include "mEntity.h"
#include "Vertex.h"
#include "Edge.h"
#include <stdio.h>

Face::Face(mDownwardAdjacencyContainer *lis ,gEntity *classification)
  : mEntity(lis,classification,2)
{
  theType = TRI;
}


//Generalized (to geomorder,thetype) face constructor schematic.
Face::Face(vector<Vertex*>& v, int Type, int geomO, gEntity *classification)
{
	vector<Vertex*>::iterator it;
	iD=0;
	theAdjacencies[0] = new mDownwardAdjacencyContainer ((int) v.size()); //Note this cast is safe since
																	//v.size is limited by geomOrder
	theType = (mType) Type;											//which has a small upper bound.
	geomOrder = geomO;
	theClassification = classification;

	for (it = v.begin(); it != v.end(); it++)
	{
		theAdjacencies[0]->add(*it);
		iD+=(*it)->getRAND();
	}

}

Face::Face(Vertex *v1, Vertex *v2, Vertex *v3,gEntity *classification) //kept for splitting compatibility
{
  theType = TRI;
  geomOrder = 1;
  theAdjacencies[0] = new mDownwardAdjacencyContainer (3);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);	
  theAdjacencies[0]->add(v3);
  theClassification = classification;
  iD = v1->getRAND()+ v2->getRAND() + v3->getRAND();
}

Face::Face(Edge *e1, Edge *e2, Edge *e3,gEntity *classification)
{
  theType = TRI;
  geomOrder = e1->getgeomOrder();
  if (!(e1->getgeomOrder() == e2->getgeomOrder() == e3->getgeomOrder()))
	  printf("this face has been constructed out of edges with different geomorders. This functionality is not currently implemented");
  theAdjacencies[1] = new mDownwardAdjacencyContainer (3);
  theAdjacencies[1]->add(e1);
  theAdjacencies[1]->add(e2);	
  theAdjacencies[1]->add(e3);	
  theClassification = classification;
  /*this will give the double of the result*/
  iD = e1->getId() + e2->getId() + e3->getId(); 
  iD /=2;
}

//only used for splitting
Face::Face(Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4, gEntity *classification)
{
  theType = QUAD;
  geomOrder = 1;
  theAdjacencies[0] = new mDownwardAdjacencyContainer (4);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);	
  theAdjacencies[0]->add(v3);	
  theAdjacencies[0]->add(v4);	
  theClassification = classification;
  iD = v1->getRAND()+ v2->getRAND() + v3->getRAND() + v4->getRAND();
}

Face::Face(Edge *e1, Edge *e2, Edge *e3, Edge *e4, gEntity *classification)
{
  theType = QUAD;
  geomOrder = e1->getgeomOrder();
  if (!(e1->getgeomOrder() == e2->getgeomOrder() == e3->getgeomOrder()))
	  printf("this face has been constructed out of edges with different geomorders. This functionality is not currently implemented");
  theAdjacencies[1] = new mDownwardAdjacencyContainer (4);
  theAdjacencies[1]->add(e1);
  theAdjacencies[1]->add(e2);	
  theAdjacencies[1]->add(e3);	
  theAdjacencies[1]->add(e4);	
  theClassification = classification;
  /*this will give the double of the result*/
  iD = e1->getId() + e2->getId() + e3->getId() + e4->getId(); 
  iD /=2;
}

int Face:: getNbTemplates (int what) const	//Note: Appears that the variable "what" serves no purpose where this is called
{                                         //Note2: Is the hierarchy redundant here?
	if (what == 1) //#of edges			  //       There's different virtual functions for types
	{									  //       as well as different options within some of the functions.
		if (getType() == TRI)
			return 3;
		if (getType() == QUAD)
			return 4;
	}
	if (what == 0) //# of vertices
	{
		if (getType() == TRI)
			return ((geomOrder+1)*(geomOrder+2))/2;
		if (getType() == QUAD)
			return (geomOrder+1)*(geomOrder+1); //twice that of TRI with the middle edge subtracted (since its counted twice)
	}												// IE: 0.5*(p+1)*(p+2)*2  -  (p+1) = p^2+2p+1 = (p+1)^2
		
	cout<<"Error: ill-formed, or not yet fully implemented type."<<endl;
	return 0;
}
//generalized
mEntity *Face::getTemplate(int ith, int what, int with)const
{
  switch(what)
    {
    case 0:  //Note: This currently does not access the secondary vertices...
		if(theAdjacencies[0]) return (*theAdjacencies[0])[ith];
		if(theAdjacencies[1])
		{
			Edge *e1 = (Edge*)get(1,ith);
			Edge *e2 = (Edge*)get(1,(ith+1)%getNbTemplates(1));
			return e1->commonVertex(e2);
		}
		break;
    case 1: 
		int iSide;
		int primarysize; //number of primary points in a face
		vector<Vertex*> v;

		switch (theType)  //probably an elegant way to avoid this case structure without much loss in efficiency
		{
			case TRI : primarysize = 3; break;
			case QUAD: primarysize = 4; break;
		}
		v.push_back((Vertex*)theAdjacencies[0]->get(ith));				//primary edge points
		v.push_back((Vertex*)theAdjacencies[0]->get((ith+1)%primarysize));
		for (int k=0;k<(geomOrder-1);k++) //loops through and adds secondary points for any
		{									//geomOrder.
			iSide=primarysize+(geomOrder-1)*(ith)+k;
			v.push_back((Vertex*)theAdjacencies[0]->get(iSide));
		}
		return new Edge (v,geomOrder,theClassification);
		break;
    }
  return (mEntity*) 0;
}

Vertex * Face::commonVertex (Face *f1, Face *f2)
{
  set<Vertex*,EntityLessThanKey> thisVertices;
  set<Vertex*,EntityLessThanKey> f1Vertices;
  set<Vertex*,EntityLessThanKey> f2Vertices;
  getVertices(thisVertices);
  f1->getVertices(f1Vertices);
  f2->getVertices(f2Vertices);
  for( set<Vertex*,EntityLessThanKey>::const_iterator iter = thisVertices.begin();
       iter != thisVertices.end();
       ++iter)
    {
      if(f1Vertices.find(*iter) != f1Vertices.end() &&
	 f2Vertices.find(*iter) != f2Vertices.end())return *iter;
    }
  return 0;
}

Edge * Face::commonEdge (Face *f1)
{
  if(!isAdjacencyCreated(1))return 0;

  for(int i=0;i<size(1);i++)
    {
      mEntity *e1 = get(1,i);
      for(int j=0;j<f1->size(1);j++)
	{
	  if(e1 == f1->get(1,j))return (Edge*)e1;
	}      
    }
  return 0;
}

Face::~Face()
{
//  if(theClassification)
    for(int i=0;i<4;i++)
      if(i <= getLevel() && theAdjacencies[i])
	for(int j=0;j<theAdjacencies[i]->size();j++)
	  theAdjacencies[i]->get(j)->del(this);
}


void Face::print()const
{
  printf("face %d with vertices ",iD);
  for (int i = 0;i<(geomOrder+1);i++)
	printf("%d ",get(0,i)->getId());
  printf("\n");
}

