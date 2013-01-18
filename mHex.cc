// mHex.cpp: implementation of the mHex class.
#include "mHex.h"
//#include "mEntityContainer.h"
#include "mEntity.h"
#include "mException.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include <stdio.h>


mHex::mHex(Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4,
	   Vertex *v5, Vertex *v6, Vertex *v7, Vertex *v8, 
	   gEntity *classification)
{
  geomOrder=1;    //Only linear implementation right now.
  for(int i=1;i<4;i++)theAdjacencies[i] = (mAdjacencyContainer*)0;
  theAdjacencies[0] = new mDownwardAdjacencyContainer (8);
  theAdjacencies[0]->add(v1);
  theAdjacencies[0]->add(v2);	
  theAdjacencies[0]->add(v3);	
  theAdjacencies[0]->add(v4);	
  theAdjacencies[0]->add(v5);	
  theAdjacencies[0]->add(v6);	
  theAdjacencies[0]->add(v7);	
  theAdjacencies[0]->add(v8);	
  theClassification = classification;
  iD = v1->getId() + v2->getId() + v3->getId() + v4->getId() +
    v5->getId() + v6->getId() + v7->getId() + v8->getId();
  theAttachable=0;
}

int mHex:: getNbTemplates (int what) const
{
  switch(what)
    {
    case 1: return 12;
    case 2: return 6;
    default : throw new mException (__LINE__,__FILE__,"");
    }
}

static int Tev[12][2] = {{0,1},{1,2},{2,3},{3,0},
			 {0,4},{1,5},{2,6},{3,7},
			 {4,5},{5,6},{6,7},{7,4}};
static int Tfv[6][4] =   {{0,3,2,1},
			 {0,1,5,4},
			 {1,2,6,5},
			 {2,3,7,6},
			 {3,0,4,7},
			 {4,5,6,7}};

mEntity *mHex::getTemplate(int ith, int what, int with)const
{
  switch(what)
    {
    case 1:
      return new Edge ((Vertex*)theAdjacencies[0]->get(Tev[ith][0]),
			(Vertex*)theAdjacencies[0]->get(Tev[ith][1]),
			theClassification);
      break;
    case 2:
      if(with == 0)
	return new Face ((Vertex*)theAdjacencies[0]->get(Tfv[ith][0]),
			  (Vertex*)theAdjacencies[0]->get(Tfv[ith][1]),
			  (Vertex*)theAdjacencies[0]->get(Tfv[ith][2]),
			  (Vertex*)theAdjacencies[0]->get(Tfv[ith][3]),
			  theClassification);
      else if (with == 1) throw new mException (__LINE__,__FILE__,"not done yet");
    }
  return 0;
}

void mHex::print()const
{
  printf("hex %d with vertices %d %d %d %d %d %d %d %d\n",iD,get(0,0)->getId(),
	 get(0,1)->getId(),get(0,2)->getId(),get(0,3)->getId(),
	 get(0,4)->getId(),get(0,5)->getId(),get(0,6)->getId(),get(0,7)->getId());
}
