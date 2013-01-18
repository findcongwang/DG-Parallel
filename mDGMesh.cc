#include "mDGMesh.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include <map>
#include <assert.h>
#include <stdio.h>
#include "mEntity.h"
using namespace std;

template <class T>
void myswap (T &t1, T &t2)
{
  T temp;
  temp = t1;
  t1 = t2;
  t2 = temp;
}

int mDGMesh::getRefinementLevel(mEntity *e)
{
if(!e->parent())return 0;
  else return 1 + getRefinementLevel(e->parent());
	//  int *i = (int*)e->getData(TYPE_RLV);
 // if(!i)return 0;
  //else return *i;
}

void attachRefinementLevel( mEntity *e, int ref)
{
  int *i = new int;
  *i = ref;
 // e->attachData("rlv",i);
}

int mDGMesh::getRefinementDepth(mEntity *e)
{
  int dim = e->getLevel();
  if(e->isAdjacencyCreated(dim))
    {
      int theMax=0;
      for(int i=0;i<e->size(dim);i++)
	{
	  mEntity *sub = e->get(dim,i);
	  int x  =  getRefinementDepth(sub) + 1;
	  theMax = (theMax>x)?theMax:x;
	}
      
       return theMax;
    }
  else return 0;
}

Vertex* mDGMesh::splitEdge(Edge *e) //PRESENT
{
  if(e->isAdjacencyCreated(1))return ((Edge*)e->get(1,0))->commonVertex((Edge*)e->get(1,1));
  mPoint pp = (e->vertex(0)->point()+e->vertex(1)->point())*0.5;
  Vertex *newv = createVertex (pp(0),pp(1),pp(2),e->getClassification());

  Edge *e1 = createEdge(e->vertex(0),newv,e->getClassification());
  Edge *e2 = createEdge(newv,e->vertex(1),e->getClassification());

  int r =  getRefinementLevel(e) + 1;
  attachRefinementLevel( e1, r);
  attachRefinementLevel( e2, r);

  e->add(e1); e->add(e2);
  e1->setParent(e);
  e2->setParent(e);

  Face *f0 = (Face*)e->get(2,0);
  Face *f1 = (Face*) 0;
  if (e->get(2,1)) f1 = (Face *) e->get(2,1);
  newv->add(f0);
  if (f1) newv->add(f1);
  newv->add((mEntity *) e1);
  newv->add((mEntity *) e2);

  int dim = f0->getLevel();
  int s = f0->size(dim-1);
  e1->add(f0); e2->add(f0);
  if (f1) { e1->add(f1);e2->add(f1);}
  for(int i=0;i<s;i++)
	{
	  //printf("sub entity :");sub->print();
	  list<mEntity *>leaves;
	  getLeaves(e,leaves);
	  list<mEntity*>::const_iterator leaves_end = leaves.end();
	  for(list<mEntity*>::const_iterator it2 = leaves.begin();it2 != leaves_end;++it2)
	    {
	      mEntity *sup = *it2;
	      //printf("leave entity :");sup->print();
	      //sup->add(f0);
	    }
  }
if (f1) 
      for(int i=0;i<s;i++)
	{
	  //printf("sub entity :");sub->print();
	  list<mEntity *>leaves;
	  getLeaves(e,leaves);
	  list<mEntity*>::const_iterator leaves_end = leaves.end();
	  for(list<mEntity*>::const_iterator it2 = leaves.begin();it2 != leaves_end;++it2)
	    {
	      mEntity *sup = *it2;
	      //printf("leave entity :");sup->print();
	      //sup->add(f1);
	    }
	  }
  allEntities.del(e);
  allSplittedEntities.add(e);
  return newv;
}

void mDGMesh::splitQuad(Face *e)
{
  if(e->isAdjacencyCreated(2))return;   //check if this cell has alredy been splitted
  Edge *e1 = (Edge*)e->get(1,0);		//get 4 edges
  Edge *e2 = (Edge*)e->get(1,1);
  Edge *e3 = (Edge*)e->get(1,2);
  Edge *e4 = (Edge*)e->get(1,3);
  Vertex* v12 = splitEdge(e1);						//split original edges in halves
  Vertex* v23 = splitEdge(e2);
  Vertex* v34 = splitEdge(e3);
  Vertex* v41 = splitEdge(e4);
  
  Vertex *v2 = e1->commonVertex(e2);     //settting notation
  Vertex *v3 = e2->commonVertex(e3);
  Vertex *v4 = e3->commonVertex(e4);
  Vertex *v1 = e4->commonVertex(e1);
  
  Edge *e11 = (Edge*)e1->get(1,0);		// get 2 new edges = 2 halves of the first edge
  Edge *e12 = (Edge*)e1->get(1,1);

  if(!e11->commonVertex(e4))			//order them
    {
      myswap<Edge*> (e11,e12);
    }

  Edge *e21 = (Edge*)e2->get(1,0);
  Edge *e22 = (Edge*)e2->get(1,1);

  if(!e21->commonVertex(e1))
    {
      myswap<Edge*> (e21,e22);
    }

  Edge *e31 = (Edge*)e3->get(1,0);
  Edge *e32 = (Edge*)e3->get(1,1);

  if(!e31->commonVertex(e2))
    {
      myswap<Edge*> (e31,e32);
    }

  Edge *e41 = (Edge*)e4->get(1,0);
  Edge *e42 = (Edge*)e4->get(1,1);

  if(!e41->commonVertex(e3))
    {
      myswap<Edge*> (e41,e42);
    }

  /*Vertex *v12 = e11->commonVertex(e12);
  Vertex *v23 = e21->commonVertex(e22);
  Vertex *v34 = e31->commonVertex(e32);
  Vertex *v41 = e41->commonVertex(e42);
*/
  Vertex *vnew  = new Vertex (theIdGenerator, 
			       (v1->point()+
			        v2->point()+
			        v3->point()+
			        v4->point()
				)*0.25,e->getClassification());
  
  Edge *enew1 = new Edge (vnew,v12,e->getClassification());
  Edge *enew2 = new Edge (vnew,v23,e->getClassification());
  Edge *enew3 = new Edge (vnew,v34,e->getClassification());
  Edge *enew4 = new Edge (vnew,v41,e->getClassification());
  
  Face *f1 = new Face (v1,v12,vnew,v41,e->getClassification());
  Face *f2 = new Face (v2,v23,vnew,v12,e->getClassification());
  Face *f3 = new Face (v3,v34,vnew,v23,e->getClassification());
  Face *f4 = new Face (v4,v41,vnew,v34,e->getClassification());

  int r =  getRefinementLevel(e) + 1;
  attachRefinementLevel( f1, r);
  attachRefinementLevel( f2, r);
  attachRefinementLevel( f3, r);
  attachRefinementLevel( f4, r);

  f1->add(e11);  f1->add(enew1); f1->add(enew4);  f1->add(e42);
  f2->add(e21);  f2->add(enew2); f2->add(enew1);  f2->add(e12);  
  f3->add(e31);  f3->add(enew3); f3->add(enew2); f3->add(e22);   
  f4->add(e41);  f4->add(enew4); f4->add(enew3);  f4->add(e32);

  enew1->add(f1); enew1->add(f2);    //new edges: add pointers to new faces
  enew2->add(f2); enew2->add(f3);
  enew3->add(f3); enew3->add(f4);
  enew4->add(f4); enew4->add(f1);

  e11->del(e);e11->add(f1);			//splitted edges: delete pointers to old elements, add pointers to new
  e12->del(e);e12->add(f2);
  e21->del(e);e21->add(f2);
  e22->del(e);e22->add(f3);
  e31->del(e);e31->add(f3);
  e32->del(e);e32->add(f4);
  e41->del(e);e41->add(f4);
  e42->del(e);e42->add(f1);
  
  allEntities.add (vnew);
  allEntities.add (enew1);
  allEntities.add (enew2);
  allEntities.add (enew3);
  allEntities.add (enew4);
  allEntities.add (f1);
  allEntities.add (f2);
  allEntities.add (f3);
  allEntities.add (f4);
  
  f1->setParent(e);
  f2->setParent(e);
  f3->setParent(e);
  f4->setParent(e);
  
  e->add(f1);
  e->add(f2);
  e->add(f3);
  e->add(f4);
  allSplittedEntities.add(e);
  allEntities.del(e);
}

void mDGMesh::splitTriangle(Face *e)
{
  if(e->isAdjacencyCreated(2))return;
  Edge *e1 = (Edge*)e->get(1,0);
  Edge *e2 = (Edge*)e->get(1,1);
  Edge *e3 = (Edge*)e->get(1,2);

  Vertex *v12 = splitEdge(e1);
  Vertex *v13 = splitEdge(e2);
  Vertex *v23 = splitEdge(e3);

  Vertex *v1 = e1->commonVertex(e2);
  Vertex *v2 = e1->commonVertex(e3);
  Vertex *v3 = e2->commonVertex(e3);

  Edge *e11 = (Edge*)e1->get(1,0);
  Edge *e12 = (Edge*)e1->get(1,1);

  if(!e11->commonVertex(e3))
    {
      myswap<Edge*> (e11,e12);
    }

  Edge *e21 = (Edge*)e2->get(1,0);
  Edge *e22 = (Edge*)e2->get(1,1);

  if(!e21->commonVertex(e1))
    {
      myswap<Edge*> (e21,e22);
    }

  Edge *e31 = (Edge*)e3->get(1,0);
  Edge *e32 = (Edge*)e3->get(1,1);

  if(!e32->commonVertex(e1))
    {
      myswap<Edge*> (e31,e32);
    }
  
  Edge *enew1 = new Edge (v12,v23,e->getClassification());
  Edge *enew2 = new Edge (v23,v13,e->getClassification());
  Edge *enew3 = new Edge (v13,v12,e->getClassification());
  
  Face *f1 = new Face (v2,v12,v23,e->getClassification());
  Face *f2 = new Face (v23,v3,v13,e->getClassification());
  Face *f3 = new Face (v13,v1,v12,e->getClassification());
  Face *f4 = new Face (v23,v12,v13,e->getClassification());

  /* Updating dependencies */

  f1->add(e11);    f1->add(enew1);  f1->add(e32);     //volume to edges
  f2->add(e31);    f2->add(e22);    f2->add(enew2);   
  f3->add(e21);    f3->add(e12);    f3->add(enew3);
  f4->add(enew1);  f4->add(enew3);  f4->add(enew2);

  f1->setParent(e);
  f2->setParent(e);
  f3->setParent(e);
  f4->setParent(e);

  e->add(f1);
  e->add(f2);
  e->add(f3);
  e->add(f4);

  v1->del(e); v1->add(f3);
  v2->del(e); v2->add(f1);
  v3->del(e); v3->add(f2);
  v12->add(f1); v12->add(f4); v12->add(f3);
  v13->add(f2); v12->add(f4); v12->add(f3);
  v23->add(f1); v12->add(f4); v12->add(f2);
  	
  v1->del(e1);  v1->del(e2); v1->add(e12); v1->add(e21);
  v2->del(e1);  v2->del(e3); v2->add(e11); v2->add(e32);
  v3->del(e2);  v3->del(e3); v3->add(e22); v1->add(e31);

  enew1->add(f1); enew1->add(f4);
  enew2->add(f2); enew2->add(f4);
  enew3->add(f3); enew3->add(f4);
  
	e11->del(e); 
	e11->add(f1);
    e12->del(e); 
	e12->add(f3);
    e21->del(e); 
	e21->add(f3);
	e22->del(e);
	e22->add(f2);
	e31->del(e); 
	e31->add(f2);
	e32->del(e); 
	e32->add(f1);

  allEntities.add (enew1);    //Adding new boundaries to the mesh  
  allEntities.add (enew2);
  allEntities.add (enew3);
  allEntities.add (f1);       //Adding new elements to the mesh
  allEntities.add (f2);
  allEntities.add (f3);
  allEntities.add (f4);

  allSplittedEntities.add(e);    
  allEntities.del(e);        //removing the parent
}

void mDGMesh::unsplitTriangle(Face *e)
{
  DEL(e->get(2,3)->get(1,0));
  DEL(e->get(2,3)->get(1,1));
  DEL(e->get(2,3)->get(1,2));
  e->get(2,0)->deleteParent();
  e->get(2,1)->deleteParent();
  e->get(2,2)->deleteParent();
  DEL(e->get(2,0));
  DEL(e->get(2,1));
  DEL(e->get(2,2));
  DEL(e->get(2,3));
  e->deleteAdjacencies(2);
}

void mDGMesh::split(mEntity* e)
{
      switch(e->getLevel())
	{
	case 0:
		printf("Error:Can't split a vertex \n");
	  break;
	case 1:
	  {
	    Edge *ed = (Edge*)e;
	    splitEdge(ed);												
	  }
	  break;
	case 2:
	  {
	    Face *f = (Face*)e;
	    if(f->getType() == mEntity::TRI)splitTriangle(f);    
	    if(f->getType() == mEntity::QUAD)splitQuad(f);
	  }
	  break;      
	case 3:
		printf("refining in 3D isn't done yet \n");
	}
}

void mDGMesh::split(int dim , list<mEntity*> &toSplit) //Only works for linear order shapes
{
	list<mEntity*>::const_iterator toSplit_end=toSplit.end();
  for(list<mEntity*>::const_iterator it = toSplit.begin();
      it != toSplit_end;
      ++it)
    {
      mEntity *e = *it;
      switch(e->getLevel())
	{
	case 0:
		printf("Error:Can't split a vertex \n");
	  break;
	case 1:
	  {
	    Edge *ed = (Edge*)e;
	    splitEdge(ed);												
	  }
	  break;
	case 2:
	  {
	    Face *f = (Face*)e;
	    if(f->getType() == mEntity::TRI)splitTriangle(f);				
	    if(f->getType() == mEntity::QUAD)splitQuad(f);
	  }
	  break;      
	case 3:
		printf("refining in 3D isn't done yet \n");
	}
    }
//  reconnect(dim);
}

void mDGMesh::reconnect (int dim)
{
	iter it;
	iter end_m1 = end(dim-1);
  for(it = begin(dim-1);it != end_m1;++it)
    {
      mEntity *e = *it;
      if(e->isAdjacencyCreated(dim-1))printf("some non final edges in the mesh !!!\n");
      e->deleteAdjacencies(dim);
    }
  iter end_n = end(dim);
  for(it = begin(dim);it != end_n;++it)
    {
      mEntity *e = *it;
      //      printf("entity :");e->print();
	  int s = e->size(dim-1);
      for(int i=0;i<s;i++)
	{
	  mEntity *sub = e->get(dim-1,i);
	  //printf("sub entity :");sub->print();
	  list<mEntity *>leaves;
	  getLeaves(sub,leaves);
	  list<mEntity*>::const_iterator leaves_end = leaves.end();
	  for(list<mEntity*>::const_iterator it2 = leaves.begin();it2 != leaves_end;++it2)
	    {
	      mEntity *sup = *it2;
	      //printf("leave entity :");sup->print();
	      sup->add(e);
	    }
	}
    }
  printf("size(%d) = %d\n",dim-1,allEntities.size(dim-1));
  list<mEntity *> toDelete;
  for(it = begin(dim-1);it != end_m1;)
    {
      mEntity *e = *it;
      ++it;
      if(!e->isAdjacencyCreated(dim))
	{	  
	  //printf("suppressing an unused edge\n");
	  toDelete.push_back(e);
	}
    }
	list<mEntity*>::const_iterator toDelete_end = toDelete.end();
  for(list<mEntity*>::const_iterator it2 = toDelete.begin();it2 != toDelete_end;++it2)
    {
      allEntities.del(*it2);
    }
  printf("size(%d) = %d\n",dim-1,allEntities.size(dim-1));
}

void mDGMesh::unsplit(int dim , list<mEntity*> &toUnsplit)
{
  int s;
  printf("phase 1...\n");
  list<mEntity*>::const_iterator toUnsplitEnd = toUnsplit.end();

  for(list<mEntity*>::const_iterator it = toUnsplit.begin();
      it != toUnsplitEnd; ++it)
    {
      mEntity *e = *it;
	  s = e->size(dim);
      for(int i=0;i<s;i++)
		  allEntities.del(e->get(dim,i));
      e->deleteAdjacencies(dim);
      allSplittedEntities.del(e);
      allEntities.add(e);
    }
  reconnect(dim);

  printf("phase 2...\n");

  list<mEntity *> toDelete;

  for(iter itt = beginsplit(dim-1);itt != endsplit(dim-1);)
    {
      mEntity *ent = *itt;
      ++itt;
      list<mEntity*> leaves;
      getLeaves(ent,leaves);
      mMeshEntityContainer c;
      int nbsides = 0;

      if(getRefinementDepth(ent) > 1)
	{
	  printf(" ref depth = %d %d leaves ",getRefinementDepth(ent),leaves.size());ent->print();
	}


      for(list<mEntity*>::const_iterator it2 = leaves.begin();it2 != leaves.end();++it2)
	{
	  mEntity *sup = *it2;
	  if(!allEntities.find(sup))
	  {
		  printf("a leave in NOT in the mesh !!!\n");
		  sup->print();
	  }
		  if(!sup->isAdjacencyCreated(dim))printf("lonely edge !!!\n");

	  if(getRefinementDepth(ent) > 1)
	    {
	      sup->print();
	    }
	  nbsides = (nbsides > sup->size(dim))?nbsides: sup->size(dim);
	  for(mUpwardAdjacencyContainer::iter it3 = sup->begin(dim);it3 != sup->end(dim);++it3)
	    {
	      c.add(*it3);
	    }
	}
      if(c.size(dim) == nbsides)
	{
	  //printf("merging an edge with a total of %d neighs\n",nbsides);
	  for(list<mEntity*>::const_iterator it2 = leaves.begin();it2 != leaves.end();++it2)
	    {
	      mEntity *sup = *it2;
	      allEntities.del(sup);
	    }
	  toDelete.push_back(ent);
	}
    }
  for(list<mEntity*>::const_iterator it2 = toDelete.begin();it2 != toDelete.end();++it2)
    {
      (*it2)->deleteAdjacencies(dim-1);
      allSplittedEntities.del(*it2);
      allEntities.add(*it2);
    }
  reconnect(dim);
}


static void recur_get_leaves(mEntity * e, list<mEntity*> &leaves, int n)
{
	int s = e->size(n);
	for(int i=0;i<s;i++)
    {
      mEntity *sub = e->get(n,i);
      if(sub->isAdjacencyCreated(n))recur_get_leaves(sub,leaves,n);
      else leaves.push_back(sub);
    }
}


void mDGMesh::getLeaves(mEntity * e, list<mEntity*> &leaves)
{
  int n = e->getLevel();
  if (e->isAdjacencyCreated(n))recur_get_leaves(e,leaves,n);
  else leaves.push_back(e);
}



