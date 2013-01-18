// mEntity.cpp: implementation of the mEntity class.
#include "mEntity.h"
#include "mException.h"
#include "gEntity.h"
#include "Vertex.h"
#include "Edge.h"
#include <algorithm>

// create an entity by giving sub-entities
mEntity::mEntity(mDownwardAdjacencyContainer *lis, gEntity *classification, int level)
  : theClassification(classification), theAttachable (0) 
{
  int i;
  for(i=0;i<4;i++)theAdjacencies[i] = (mAdjacencyContainer*)0;
  theAdjacencies[lis->get(0)->getLevel()] = new mDownwardAdjacencyContainer (*lis);
  for(i=0;i<lis->size();i++)
    {
      if(lis->get(i)->isAdjacencyCreated(level))lis->get(i)->add(this);
    }
  computeId();
}
// create an empty entity
mEntity::mEntity()
  : theClassification((gEntity*)0), theAttachable (0) 
{
  iD = 0;
  cell = 0;
  for(int i=0;i<4;i++)theAdjacencies[i] = (mAdjacencyContainer*)0;
}

mEntity::~mEntity()
{
  // cout << "destroying " << iD <<"\n";
  for(int i=0;i<4;i++)
    {
      if(theAdjacencies[i])
	delete theAdjacencies[i];
    }
  if(theAttachable)
    {
      mAttachableDataContainer::iter it = theAttachable->begin();
      mAttachableDataContainer::iter itEnd = theAttachable->end();
      for(;it!=itEnd;++it)
	{
	  mAttachableData *a = (*it).second;
	  delete a;
	}
      delete theAttachable;
    }
}

void mEntity::printAllAttachable () const
{
  if(theAttachable)
    {
      mAttachableDataContainer::iter it = theAttachable->begin();
      mAttachableDataContainer::iter itEnd = theAttachable->end();
      for(;it!=itEnd;++it)
	{
	  cout << (*it).first << " ";
	}
      cout << "\n";
    }  
}

int mEntity::getNbTemplates(int what) const
{
  throw new mException (__LINE__,__FILE__,"no template in mEntity base class");
}

mEntity* mEntity::getTemplate(int,int,int)const 
{
  throw new mException (__LINE__,__FILE__,"no template in mEntity base class");
}

mEntity * mEntity::find(mEntity*me)const
{
	mEntity *ret;
	int level = me->getLevel();

  if(theAdjacencies[level] &&
	  (ret = theAdjacencies[level]->find(me)))
	  return ret;
  return (mEntity *)0;
}

void mEntity::add(mEntity*me)
{
  int thisOrder = getLevel();
  int meOrder = me->getLevel();
  
  if(!theAdjacencies[meOrder] && meOrder > thisOrder) 
    theAdjacencies[meOrder] = new mUpwardAdjacencyContainer;
  else if(!theAdjacencies[meOrder]) 
    theAdjacencies[meOrder] = new mDownwardAdjacencyContainer(getNbTemplates(meOrder));
  theAdjacencies[meOrder]->add(me);
}

void mEntity::del(mEntity*me)
{
  int meOrder = me->getLevel();
  if(theAdjacencies[meOrder])theAdjacencies[meOrder]->del(me);
}

void mEntity::classify(gEntity *g)
{
  theClassification = g;
  if(parent())parent()->classify(g);
}

void mEntity::deleteAdjacencies(int what)
{
  if(theAdjacencies[what])
    {
      delete theAdjacencies[what];
      theAdjacencies[what] = (mAdjacencyContainer*)0;
    }
}

// are two mesh entities equal ?
bool mEntity::equal (mEntity *other) const
{
  int i;
  int level = getLevel();
  if(level != other->getLevel())return false;
  if(level == 0)return iD == other->getId();
  
  if(level == 1)
    {
      unsigned long int v11 = get(0,0)->getId();
      unsigned long int v12 = get(0,1)->getId();
      unsigned long int v21 = other->get(0,0)->getId();
      unsigned long int v22 = other->get(0,1)->getId();
      if(v11 == v21 && v12 == v22)return true;
      if(v11 == v22 && v12 == v21)return true;
      return false;
    }
  
  if(isAdjacencyCreated(0) && other->isAdjacencyCreated(0) && size(0) < 100)
    {
      unsigned long int v1[100],v2[100];
      int s1 = size(0);
      int s2 = other->size(0);
      if(s1 != s2)return false;
      for(i=0;i<s1;i++)
	{
	  v1[i]=get(0,i)->getId();
	  v2[i]=other->get(0,i)->getId();
	}
      sort(v1,v1+s1);
      sort(v2,v2+s2);
      for(i=0 ; i<s1 ;i++)
		if(v1[i] != v2[i])return false;
      return true;
    }
  
  set<Vertex*,EntityLessThanKey> thisVertices;
  set<Vertex*,EntityLessThanKey> otherVertices;
  getVertices(thisVertices);
  other->getVertices(otherVertices);
  if(thisVertices.size() != otherVertices.size())return false;
  set<Vertex*,EntityLessThanKey>::const_iterator iter1 = thisVertices.begin();
  set<Vertex*,EntityLessThanKey>::const_iterator iter2 = otherVertices.begin();
  for( ;
       iter1 != thisVertices.end();
       ++iter1,++iter2)
    if((*iter1)->getId() != (*iter2)->getId())return false;

  return true;
}

// compare two entities
bool mEntity::lessthan (mEntity *other) const
{
	//if (this <other) return true; else return false; 
  /* first do some fast checks**/
  if(getLevel() != other->getLevel())return getLevel() < other->getLevel();
  if(getLevel() == 0)return getId() < other->getId();
  if(getId() < other->getId())return true;
  if(getId() > other->getId())return false;

  /*then do the whole comparison on vertices*/
  set<Vertex*,EntityLessThanKey> thisVertices;
  set<Vertex*,EntityLessThanKey> otherVertices;
  getVertices(thisVertices);
  other->getVertices(otherVertices);
  if(thisVertices.size() < otherVertices.size())return true;
  if(thisVertices.size() > otherVertices.size())return false;
  set<Vertex*,EntityLessThanKey>::const_iterator iter1 = thisVertices.begin();
  set<Vertex*,EntityLessThanKey>::const_iterator iter2 = otherVertices.begin();
  for( ;
       iter1 != thisVertices.end();
       ++iter1,++iter2)
    {
      if((*iter1)->getId() < (*iter2)->getId())return true;
      if((*iter1)->getId() > (*iter2)->getId())return false;
    }

  return false;
}
//With addition of higher order entities, order of vertices becomes important. Perhaps the use
//of a set should be removed.
void mEntity::getVertices(set<Vertex*, EntityLessThanKey> &vert) const
{
	//remove if
  if(theAdjacencies[0])
    {
	  int Size = size(0);
      for(int i=0;i<Size;i++)vert.insert((Vertex*)theAdjacencies[0]->get(i));
      return;
    }
  //Is this needed? 06/10/03
  int level = getLevel();
  for(int j=1;j<level;j++)
    {
      if(theAdjacencies[j])
	{
	  for(int i=0;i<size(j);i++)theAdjacencies[j]->get(i)->getVertices(vert);
	  return;
	}
    }
  throw new mException (__LINE__,__FILE__, 
			"Manipulation of an entity defined without any downward adjacencies set");
}
//Overloaded: Uses a vector instead of a set.
void mEntity::getVertices(vector<Vertex*> &vert) const
{
	//remove if
  if(theAdjacencies[0])
    {
	  int Size = size(0);
      for(int i=0;i<Size;i++)vert.push_back((Vertex*)theAdjacencies[0]->get(i));
      return;
    }
  //Is this needed? 06/10/03
  int level = getLevel();
  for(int j=1;j<level;j++)
    {
      if(theAdjacencies[j])
	{
	  for(int i=0;i<size(j);i++)theAdjacencies[j]->get(i)->getVertices(vert);
	  return;
	}
    }
  throw new mException (__LINE__,__FILE__, 
			"Manipulation of an entity defined without any downward adjacencies set");
}
void mEntity::computeId()
{
  if(getLevel() == 0) return;
  iD = 0;
  set<Vertex*,EntityLessThanKey> thisVertices;
  getVertices(thisVertices);
  set<Vertex*,EntityLessThanKey>::const_iterator end = thisVertices.end();
  for( set<Vertex*,EntityLessThanKey>::const_iterator iter = thisVertices.begin();
       iter != end; ++iter)
      iD += (*iter)->getRAND();
}

/************   ATTACHABLE DATA  ******************************/

void mEntity::attachData(int c ,  mAttachableData *v)
{
  if(!theAttachable)theAttachable = new mAttachableDataContainer;
  theAttachable->set(c,v);
}

void mEntity::deleteData(int c)
{
  if(!theAttachable)return;
  theAttachable->del(c);
}

void mEntity::setOwner (mMesh *m)
{
  mAttachableMesh *am = new mAttachableMesh;
  am->m = m;
  attachData(TYPE_OWNER,am);
}


mMesh * mEntity::getOwner ()
{
  mAttachableMesh *am = (mAttachableMesh *)getData(TYPE_OWNER);
  if(!am)return 0;
  return am->m;
}

void mEntity::attachInt (int c, int i)
{
  mAttachableInt *ai = (mAttachableInt *)getData(c);
  if(!ai)
    {
      ai = new mAttachableInt;
      attachData(c,ai);
    }
  ai->i = i;
}

void mEntity::setParent (mEntity* e)
{
  mAttachableEntity *ai = (mAttachableEntity *)getData(TYPE_P);
  if(!ai)
    {
      ai = new mAttachableEntity;
      attachData(TYPE_P,ai);
    }
  ai->e = e;
}

mEntity* mEntity::parent ()
{
  mAttachableEntity *ai = (mAttachableEntity *)getData(TYPE_P);
  if(!ai)return 0;
  return ai->e;
}

void mEntity::deleteParent ()
{
  mAttachableEntity *ai = (mAttachableEntity *)getData(TYPE_P);
  if(!ai)return;
  deleteData(TYPE_P);
}

void mEntity::attachDouble (int c, double d)
{
  mAttachableDouble *ad = (mAttachableDouble *)getData(c);
  if(!ad)
    {
      ad = new mAttachableDouble;
      attachData(c,ad);
    }
  ad->d = d;
}

double mEntity::getAttachedDouble ( int c)
{
  mAttachableDouble *ad = (mAttachableDouble *)getData(c);
  if(!ad)return 0.0;
  return ad->d;
}

/**********************  LEAVES  ******************************/

static void recur_get_leaves(mEntity * e, list<mEntity*> &leaves, int n)
{
	int Size = e->size(n);
  for(int i=0;i<Size;i++)
    {
      mEntity *sub = e->get(n,i);
      if(sub->isAdjacencyCreated(n))recur_get_leaves(sub,leaves,n);
      else leaves.push_back(sub);
    }
}

void mEntity::getLeaves(list<mEntity*> &leaves)
{
  int n = getLevel();
  if (isAdjacencyCreated(n))recur_get_leaves(this,leaves,n);
  else leaves.push_back(this);
}

void recur_get_all(mEntity * e, list<mEntity*> &leaves, int n)
{
	int Size = e->size(n);
  for(int i=0;i<Size;i++)
    {
      mEntity *sub = e->get(n,i);
      leaves.push_back(sub);
      if(sub->isAdjacencyCreated(n))recur_get_all(sub,leaves,n);
    }
}

void mEntity::getAllSubTree(list<mEntity*> &leaves)
{
  int n = getLevel();
  leaves.push_back(this);
  if (isAdjacencyCreated(n))recur_get_all(this,leaves,n);
}




