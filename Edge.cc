// Edge.cpp: implementation of the Edge class.
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include "Edge.h"
#include "mEntity.h"
#include "Vertex.h"

Edge::Edge(Vertex *v1, Vertex *v2, gEntity *classification) //older function.
{															//currently being kept to keep splitting//refinement compatible
  edgeType=EDGE;											//GRANT
  theAdjacencies[0] = new mDownwardAdjacencyContainer (2);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);
  theClassification = classification;
  unsigned long int i1 = v1->getId();
  unsigned long int i2 = v2->getId();
  iD = v1->getRAND() + v2->getRAND();
}
Edge::Edge(vector<Vertex*>& v,int geomO, gEntity *classification)
{
	vector<Vertex*>::iterator it;
	iD=0;
	theAdjacencies[0] = new mDownwardAdjacencyContainer (v.size());
	
	edgeType=EDGE;	
	geomOrder = geomO;
	theClassification = classification;
												//Note:In the older version, id's were assigned as 
	for (it = v.begin(); it != v.end(); it++)	//unsigned long ints but not used. what was this for?
	{											//GRANT
		theAdjacencies[0]->add(*it);
		iD+=(*it)->getRAND();
	}
}
Edge::~Edge()
{
//printf("Deleting edge with id = %d \n",iD);
int j, nb;
if (theAdjacencies[0])
{
nb = theAdjacencies[0]->size();
for (j=0;j<nb;++j)                   //deleting pointers to this element from its vertices
theAdjacencies[0]->get(j)->del(this);
}
if (theAdjacencies[1])
{
nb = theAdjacencies[1]->size();     //deleting pointer to this element from its parent(?) or child(?)
for (j=0;j<nb;++j)
theAdjacencies[1]->get(j)->del(this);
}
if (theAdjacencies[2])
{
printf("faces should be destructed first \n");
exit(0);
}
if (theAdjacencies[3])
{
printf("~Edge: don't know how to destruct in3D\n");
exit(0);
}
}

int Edge::getNbTemplates(int what)const  //returns the number of vertices
{
  if(what == 0)
	  return (theAdjacencies[0]->size());
  else return 0;
}

Vertex *Edge::vertex(int i) const
{
  return (Vertex*)theAdjacencies[0]->get(i);
}

Vertex *Edge::commonVertex(Edge *other) const
{ //Note that the secondary vertices are not compared
  //as this information is currently irrelevant.
  if(other->vertex(0) == vertex(0))return vertex(0);
  if(other->vertex(1) == vertex(0))return vertex(0);
  if(other->vertex(0) == vertex(1))return vertex(1);
  if(other->vertex(1) == vertex(1))return vertex(1);
  return 0;
}

void Edge::print()const
{
  printf("edge %p with Id %d and vertices ",this,iD);
  for (int i = 0;i<(geomOrder+1);i++)
	printf("%d ",get(0,i)->getId());
  printf("\n");
}
