#include "mMesh.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "mTet.h"
#include "mHex.h"
#include "mException.h"
#include <stdio.h>
#include <math.h>
#include <cstdlib>

mMesh::~mMesh()
{
  for(int i=3;i>=0;i--)
    {
      printf("destroying level %d %d\n",i,size(i));
      for(iter it = begin(i);it != end(i);)
	{
	  switch(i)
	    {
	    case 0:
	      {
		Vertex *e = dynamic_cast<Vertex *>(*it);
		++it;
		delete e;
	      }
	      break;
	    case 1:
	      {
		Edge *e = dynamic_cast<Edge *>(*it);
		++it;
		delete e;
	      }
	      break;
	    case 2:
	      {
		Face *e = dynamic_cast<Face *>(*it);
		++it;
		delete e;
	      }
	      break;
	    case 3:
	      {
		mRegion *e = dynamic_cast<mRegion *>(*it);
		++it;
		delete e;
	      }
	      break;
	    }
	}
    }
}



Vertex *mMesh::getVertex (int theId)
{
  //Given ID, returns the vertex with this iD
  Vertex v(theId,mPoint(0,0,0),0);
  Vertex* found = (Vertex*) find(&v);
  if (found) return found;
  else {
	  printf("Can't find a vertex !");
	  exit(0);
  }
}
//generalized edge finder
Edge *mMesh::getEdge (vector<Vertex*>& v, int geomO)
{
  Edge e(v,geomO,0);
  Edge *theE = (Edge*)checkExisting(&e);
  if(theE != &e)return theE;
    else return (Edge*)find(&e);
}
//generalized face finder
Face* mMesh::getFace (vector<Vertex*>& v, int Type, int geomO)
{
  Face f(v,Type,geomO,0);
  Face *found = (Face*)find(&f);
  if (found) return found;
  else printf("Can't find a face!");
  return (Face*) 0;
}
//generalized tet finder
mTet* mMesh::getTet (vector<Vertex*>& v, int geomO)
{
	mTet t(v,geomO,0);
	mTet *theT = (mTet*)checkExisting(&t);
	if(theT != &t)return theT;
	return (mTet*)find(&t);
}
Vertex *mMesh::createVertex (int theId, const double x, const double y, const double z, gEntity *classif)
{
  Vertex *theNewVertex;
  theNewVertex = new Vertex (theId,mPoint(x,y,z),classif);
  theIdGenerator.setMaxValue(theId);
  allEntities.add(theNewVertex);
  return theNewVertex;
}

Vertex *mMesh::createVertex (const double x, const double y, const double z, gEntity *classif)
{
  return createVertex(theIdGenerator.generateId(),x,y,z,classif);
}

//generalized
Edge *mMesh::createEdge (vector<int>& vertexiD, int geomOrder, gEntity *classif)
{
  Edge *theNewEdge;
  vector<Vertex*> v;
  vector<int>::iterator it;

  for (it = vertexiD.begin(); it != vertexiD.end(); it++)
	v.push_back(getVertex(*it));

  bool compoundUnit=0;
  vector<Vertex*>::iterator it2;
  for (it2 = v.begin(); it2 != v.end(); it2++)
	compoundUnit=compoundUnit || !(*it2);
  if(compoundUnit)
    {
      throw new mException (__LINE__,__FILE__,"unknown verex id's in new edge creation");
    }
#ifdef _DEBUG_
  if(theNewEdge = getEdge(v,geomOrder))
    {
      throw new mException (__LINE__,__FILE__,"trying to create a vertex with an used ID");
    }
#endif
  theNewEdge = new Edge(v,geomOrder,classif);
  return theNewEdge;
}

//Left in to keep splitting/refinement compatible in an easy fashion. //GRANT
Edge* mMesh::createEdge (Vertex *v1, Vertex *v2, gEntity *classif)
{
  Edge *theNewEdge = new Edge(v1,v2,classif);
  Edge *e = (Edge*)checkExisting(theNewEdge);
  if(e != theNewEdge)
    {
      delete theNewEdge;
      theNewEdge = e;
    }
  else
    {
      add(theNewEdge);
    }
  return theNewEdge;
}

Face *mMesh::createFaceWithVertices (vector<int>& vertexiD, int theType, int geomOrder, gEntity *classif)
{
  Face *theNewFace;
  vector<Vertex*> v;
  vector<int>::iterator it;

  for (it = vertexiD.begin(); it != vertexiD.end(); it++)
	v.push_back(getVertex(*it));

  if (theType == 3) //for quad
  {
	    Vertex* v[4];	
  v[0] = getVertex(vertexiD[0]);
  v[1] = getVertex(vertexiD[1]);
  v[2] = getVertex(vertexiD[2]);
  v[3] = getVertex(vertexiD[3]);
  Vertex *v1, *v2, *v3, *v4;
  
  mPoint p; //find the geometric center  
  p(0)=0.25*(v[0]->point()(0)+v[1]->point()(0)+v[2]->point()(0)+v[3]->point()(0));
  p(1)=0.25*(v[0]->point()(1)+v[1]->point()(1)+v[2]->point()(1)+v[3]->point()(1));

  /* Reorder vertices so that 
    (x_{n-1},y_{n-1}) -> (-1,-1)
    (x_{n},y_{n-1})   -> (1,-1)
    (x_{n},y_{n})     -> (1,1)
    (x_{n-1},y_{n})   -> (-1,1) */

  /*for (int i=0;i<4;i++) // Lilia's version - might be important for limiting
	{
		if (v[i]->point()(0)<p(0)&&v[i]->point()(1)<p(1)) v1=v[i];
		if (v[i]->point()(0)>p(0)&&v[i]->point()(1)<p(1)) v2=v[i];
		if (v[i]->point()(0)>p(0)&&v[i]->point()(1)>p(1)) v3=v[i];
		if (v[i]->point()(0)<p(0)&&v[i]->point()(1)>p(1)) v4=v[i];
	}*/

  int i=0;  //newer version by Ruibin. needs to be here for quads.
  int myind=1;
	{
		if (v[i]->point()(0)<p(0)&&v[i]->point()(1)<p(1)) {v1=v[i]; myind=1;}
		if (v[i]->point()(0)>p(0)&&v[i]->point()(1)<p(1)) {v2=v[i]; myind=2;}
		if (v[i]->point()(0)>p(0)&&v[i]->point()(1)>p(1)) {v3=v[i]; myind=3;}
		if (v[i]->point()(0)<p(0)&&v[i]->point()(1)>p(1)) {v4=v[i]; myind=4;}
	}
	if(myind==1){
         v2=v[1];
		 v3=v[2];
		 v4=v[3];
	}
	if(myind==2){
		v3=v[1];
		v4=v[2];
		v1=v[3];
	}
	if(myind==3){
		v4=v[1];
		v1=v[2];
		v2=v[3];
	}
	if(myind==4){
		v1=v[1];
		v2=v[2];
		v3=v[3];
	}
  if(!v1 || !v2 || !v3 || !v4)
    {
      char text[256];
      sprintf(text,"trying to create a quad with vertices %d %d %d %d, search gives %p %p %p %p\n",
	      vertexiD[0],vertexiD[1],vertexiD[2],vertexiD[3],v1,v2,v3,v4);
      throw new mException (__LINE__,__FILE__,text);
    }
  }



#ifdef _DEBUG_

  bool compoundUnit=1;
  vector<Vertex*>::iterator it2;
  for (it2 = v.begin(); it2 != v.end(); it2++)
	compoundUnit=compoundUnit || !(*it2);

  if(compoundUnit)
    {
      throw new mException (__LINE__,__FILE__,"unknown verex id's in new face creation");
    }
  if(theNewFace = getFace(v))
    {
      throw new mException (__LINE__,__FILE__,"trying to create a face already existing");
    }
#endif
  theNewFace = new Face(v,theType,geomOrder,classif);
  return theNewFace;
}
mTet *mMesh::createTetWithVertices (vector<int>& vertexiD, int geomOrder, gEntity *classif)
{
  mTet *theNewTet;
  vector<Vertex*> v;
  vector<int>::iterator it;

  for (it = vertexiD.begin(); it != vertexiD.end(); it++)
	v.push_back(getVertex(*it));

  #ifdef _DEBUG_
  if(theNewTet = getTet(v,geomOrder))
    {
      throw new mException (__LINE__,__FILE__,"trying to create a tet already existing");
    }
  #endif

  theNewTet = new mTet(v,geomOrder,classif);
  return theNewTet;
}

mTet *mMesh::createTetWithFaces (Face *f1, Face *f2, Face *f3, Face *f4, gEntity *classif)
{
  mTet *theNewTet;
  theNewTet = new mTet(f1,f2,f3,f4,classif);
  allEntities.add(theNewTet);
  return theNewTet;
}

Face *mMesh::createFaceWithEdges (Edge *e1, Edge *e2, Edge *e3, Edge *e4, gEntity *classif)
{
  Face *theNewFace;
  theNewFace = new Face(e1,e2,e3,e4,classif);
  allEntities.add(theNewFace);
  return theNewFace;
}

Face *mMesh::createFaceWithEdges (Edge *e1, Edge *e2, Edge *e3, gEntity *classif)
{
  Face *theNewFace;
  theNewFace = new Face(e1,e2,e3,classif);
  allEntities.add(theNewFace);
  return theNewFace;
}


mHex *mMesh::createHexWithVertices (int i1, int i2, int i3, int i4, 
				    int i5, int i6, int i7, int i8, 
				    gEntity *classif)
{
  Vertex *v1 = getVertex(i1);
  Vertex *v2 = getVertex(i2);
  Vertex *v3 = getVertex(i3);
  Vertex *v4 = getVertex(i4);
  Vertex *v5 = getVertex(i5);
  Vertex *v6 = getVertex(i6);
  Vertex *v7 = getVertex(i7);
  Vertex *v8 = getVertex(i8);
  if(!v1 || !v2 || !v3 || !v4 ||
     !v5 || !v6 || !v7 || !v8)
    {
      throw new mException (__LINE__,__FILE__,"unknown verex id's in new face creation");
    }
#ifdef _DEBUG_
#endif
  mHex *theNewHex = new mHex(v1,v2,v3,v4,v5,v6,v7,v8,classif);
  allEntities.add(theNewHex);
  return theNewHex;
}

void mMesh::resize(int what, int size)
{
  allEntities.resize(what,size);
}

void mMesh::createCellBoundaries (int i,int j, int with)
{
  //comments are for 2D.... Same relationships hold in 3D
  // mesh elements have been created using appropriate number of vertices
  // Now we are creating edges for these elements
  //i=2(cells),j=1(edges) 
  iter end_i=end(i);
  for(iter it=begin(i);it!=end_i;++it)   //loop over all cells
    {
      mEntity *e = (*it); //get a cell "e"
      int NbBoundaries = e->getNbTemplates(j); //number of edges
      for(int k=0;k<NbBoundaries;++k)
	    {
	      mEntity *t = e->getTemplate(k,j,with);// create k-th edge using kth and (k+1)st vertices
	      mEntity *q;
		  //now we are checking if this edge has been already created by the neiboring face
	      if(q = allEntities.find(t)) // find returns t, if it already exists and zero if it doesn't
		{
		  e->add(q);   // add a pointer to this edge to the face 
		  if(q != t) 
		    {
			  delete t;
		    }
		  else {printf("weird edge in Mesh.cc\n");exit(0);}
		}
	      else
		{
		  e->add(t);  // adds a pointer to just created edge(face) t in element e
		  allEntities.add(t); // add a pointer to t to the set of all elements.
		}
	    }
    }
}

void mMesh::createConnections(int from, int to)
{
	//adding a pointer to higher order elements in a lower element( to faces in an edge)
	if (from < to)  // j corresponds to higher order elements
	{
	  iter end_j=end(to);
      for(iter it=begin(to);it!=end_j;it++)
	{
	  mEntity *e = (*it);  
	  const int s =  e->size(from);
	  for(int k = 0;k<s;k++)
	      e->get(from,k)->add(e);
	}
	}
	else
	{
/*iter end_j=end(from);
for(iter it=begin(to);it!=end_j;it++)
mEntity* e=(*it);
*/
		return;
	}

}

gEntity *mMesh::getGEntity(int id, int level)
{
  gEntity gg(id,level);    //rewrite
  set<gEntity*,gEntityLessThanKey>::const_iterator it = allGEntities.find(&gg);
  if(it != allGEntities.end())
    return *it;
  gEntity *pg = new gEntity(id,level);
  allGEntities.insert(pg);
  return pg;
}

mEntity* mMesh::checkExisting(mEntity *e)
{
  return e;
}

mEntity* mMesh::find(mEntity *e)
{
  return allEntities.find(e);
}

void mMesh:: writeMEntity(ostream &o, mEntity *m)
{
	/*type of an element : triangle, hex, etc. */
  o << m->getType() << " " ;
  int NbTemplates = m->getNbTemplates(0);
  /*writing NbTemplates vertices */
  for(int i=0;i<NbTemplates;i++)
    {
      o << m->get(0,i)->getId() << " " ;
    }
  o << "\n";
}
mEntity * mMesh:: readMEntity(istream &is)
{
  int typ,geomOrder,temp;
  vector<Vertex*> iv;
  vector<Vertex*>::iterator it;


  is >> typ;
  is >> geomOrder; //Note: recoveryFile needs to be created such that geomOrder is in it.
  switch(typ)
    {
    case mEntity::EDGE :
      for (int i=0;i<(geomOrder+1);i++)
	  {
		  is >> temp;
		  iv.push_back(getVertex(temp));
	  }
      return getEdge(iv,geomOrder);
      break;
    case mEntity::TRI :    //faces
	case mEntity::QUAD:
      {
      for (int i=0;i<0.5*(geomOrder+1)*(geomOrder+2);i++)
	  {
		  is >> temp;
		  iv.push_back(getVertex(temp));
	  }
	  Face *found = getFace(iv,typ,geomOrder);
      return found;
      }
      break;
	case mEntity::TET :
       for (int i=0;i<((geomOrder+1)*(geomOrder+2)*(geomOrder+3))/6;i++)  
	  {																		
		  is >> temp;														
		  iv.push_back(getVertex(temp));
	  }
      return getTet(iv,geomOrder);
      break;
	default:
		printf("Don't know how to read this element\n");
		exit(0);
    }
}
void mMesh::setPeriodicBC(int dim) 
{
	iter it, it2;
	iter mesh_end_nm1 = end(dim-1);
	double SN = 1.0e-10;
    for(it = begin(dim-1);it !=mesh_end_nm1 ;++it)
    {
      mEntity *ent = (*it);
  
      if (ent->getClassification()->getId() == 510)
	for(it2 = begin(dim-1);it2 != mesh_end_nm1;++it2)
	  {
		mEntity *another =(*it2);
		if (another->getClassification()->getId() == 510 && ent != another)
	      {
		//	ent->print();
		//	another->print();
			
		Vertex* v1 = (Vertex*) ent->get(0,0);
		Vertex* v2 = (Vertex*) ent->get(0,1);
		Vertex* va1 = (Vertex*) another->get(0,0);
		Vertex* va2 = (Vertex*) another->get(0,1);
		if ((fabs(v1->point()(1)-va1->point()(1))<SN && fabs(v2->point()(1)-va2->point()(1))<SN)|| //horizontal periodicity
		    (fabs(v1->point()(1)-va2->point()(1))<SN && fabs(v2->point()(1)-va1->point()(1))<SN)) 
		  {
			mAttachableEntity *ae = new mAttachableEntity;
             ae->e = *it2;
 			 ent->attachData(TYPE_CELL,ae);
		    /*mUpwardAdjacencyContainer::iter it3;
		    mUpwardAdjacencyContainer::iter end_np1 = another->end(dim);
		    
		    for( it3= another->begin(dim); it3 !=  end_np1; ++it3)
		      {
			if (*ent->begin(2) != (*it3))
			  ent->add(*it3);
			(*it3)->get(1)->replace(another,ent);
			  }
			allEntities.del(another);*/
		    break;
		  }
	      }
	  }

      if (ent->getClassification()->getId() == 610)
	for(it2 = begin(dim-1);it2 != mesh_end_nm1;++it2)
	  {
	    mEntity *another = (*it2);
	    //int n = ent->getLevel();

	    if (another->getClassification()->getId() == 610 && ent != another) //vertical periodicity
	      {
		Vertex* v1 = (Vertex*) ent->get(0,0);
		Vertex* v2 = (Vertex*) ent->get(0,1);
		Vertex* va1 = (Vertex*) another->get(0,0);
		Vertex* va2 = (Vertex*) another->get(0,1);
		if ((fabs(v1->point()(0)-va1->point()(0))<SN && fabs(v2->point()(0)-va2->point()(0))<SN)|| //horizontal periodicity
		    (fabs(v1->point()(0)-va2->point()(0))<SN && fabs(v2->point()(0)-va1->point()(0))<SN)) 
	  {
			mAttachableEntity *ae = new mAttachableEntity;
             ae->e = *it2;
 			 ent->attachData(TYPE_CELL,ae);
		  
		    break;
		  }
	      }
    }
  }
}



