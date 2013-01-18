// mTet.cpp: implementation of the mTet class.
#include "mTet.h"
#include "mEntity.h"
#include "mException.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"

//generalized for any order
mTet::mTet(vector<Vertex*>& v, int geomO, gEntity *classification)
{
	vector<Vertex*>::iterator it;
	iD=0;
	theAdjacencies[0] = new mDownwardAdjacencyContainer ((int) v.size()); //Note this cast is safe since
																	//v.size is limited by geomOrder
																	//which has a small upper bound.
	geomOrder = geomO;
	theClassification = classification;

	for (it = v.begin(); it != v.end(); it++)
	{
		theAdjacencies[0]->add(*it);
		iD+=(*it)->getRAND();
	}
	theAttachable=0;
}
mTet::mTet(Face *f1, Face *f2, Face *f3, Face *f4, gEntity *classification)
{
  geomOrder = f1->getgeomOrder();
  if (!(f1->getgeomOrder() == f2->getgeomOrder() == f3->getgeomOrder() == f4->getgeomOrder()))
	  printf("This tet has been constructed out of faces with different geomorders. This mixed functionality is not currently implemented");
  for(int i=1;i<4;i++)theAdjacencies[i] = (mAdjacencyContainer*)0;
  theAdjacencies[2] = new mDownwardAdjacencyContainer (4);
  theAdjacencies[2]->add(f1);
  theAdjacencies[2]->add(f2);	
  theAdjacencies[2]->add(f3);	
  theAdjacencies[2]->add(f4);	
  theClassification = classification;
  computeId();
  theAttachable=0;
}

int mTet:: getNbTemplates (int what) const
{
  switch(what)
    {
    case 0: return ((geomOrder+1)*(geomOrder+2)*(geomOrder+3))/6; break;
    case 1: return 6; break;
    case 2: return 4; break;
    default : throw new mException (__LINE__,__FILE__,"");
    }
}
//Numbers in variable names represent geometric order.
static int Tev1[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
static int Tev2[6][3] = {{0,1,4},{0,2,6},{0,3,7},{1,2,5},{1,3,9},{2,3,8}};
static int Tfv1[4][3] = {  {0,1,2},
						   {0,1,3},  //should these geomorders be grouped into 1 thing
						   {0,2,3},  //or is it worth writing out each case explicitly for efficiency?
						   {1,2,3}        };
static int Tfv2[4][6] = {  {0,1,2,4,5,6},
						   {0,1,3,4,9,7},
						   {0,2,3,6,8,7},
						   {1,2,3,5,8,9}  };
static int Tve[4][2] = {{0,1},{0,3},{1,3},{2,4}};
static int Tvf[4][2] = {{0,1},{0,3},{1,3},{2,4}};

mEntity *mTet::getTemplate(int ith, int what, int with)const 
{
  switch(what)
    {
    case 0:  //Note: This currently only accesses primary vertices.
      if(theAdjacencies[1])
	return (((Edge*)theAdjacencies[1]->get(Tve[ith][0]))->commonVertex((Edge*)theAdjacencies[1]
									    ->get(Tve[ith][1])));
      else if (theAdjacencies[2])
	{
	  return 0;
	}
      break;
    case 1:
	  {
		  vector<Vertex*> v;
		  for (int k=0; k<(geomOrder+1);k++)
		  {
			  switch(geomOrder)
			  {
				case 1: v.push_back((Vertex*)theAdjacencies[0]->get(Tev1[ith][k])); break;
				case 2: v.push_back((Vertex*)theAdjacencies[0]->get(Tev2[ith][k])); break;
				default: printf("Error: Either geomOrder is not properly defined, or case has not been coded yet");
					return 0;
			  }
		  }
		  return new Edge(v,geomOrder,theClassification);
	  }
      break;
    case 2:
      if(with == 0)
	  {
		  vector<Vertex*> v;
		  for (int k=0; k<((geomOrder+1)*(geomOrder+2))/2;k++)
		  {
			  switch(geomOrder)
			  {
				case 1: v.push_back((Vertex*)theAdjacencies[0]->get(Tfv1[ith][k])); break;
				case 2: v.push_back((Vertex*)theAdjacencies[0]->get(Tfv2[ith][k])); break;
				default: printf("Error: Either geomOrder is not properly defined, or case has not been coded yet");
					return 0;
			  }
		  }
		  return new Face(v,TRI,geomOrder,theClassification);
	  }
      else if (with == 1) throw new mException (__LINE__,__FILE__,"not done yet");
    default : throw new mException (__LINE__,__FILE__,"this tet template is not done");
    }
  return 0;
}


