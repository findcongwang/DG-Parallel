//#include "mEntityContainer.h"
#include "mEntity.h"
#include "mException.h"


/************* mDOWNWARDAdjacencyCONTAINER *******************/
mDownwardAdjacencyContainer::mDownwardAdjacencyContainer(int i) 
{ 
  mContainer.reserve(i);
}

void mDownwardAdjacencyContainer::del(mEntity* e)
{ 
  throw new mException(__LINE__,__FILE__,"deletion operation NOT PERMITTED in this context");
}

mEntity*  mDownwardAdjacencyContainer::find(mEntity* ent) const 
{ 
#ifdef _DEBUG_
  throw new mException(__LINE__,__FILE__,"inefficient way to find in a vector<>");
#endif
 int size = (int) mContainer.size();
 for(int i=0;i<size;++i)
   if(get(i)->equal(ent)) return get(i);
 return 0;
}

mEntity *mDownwardAdjacencyContainer::replace(mEntity *toreplace, mEntity* with) 
{
int size = (int) mContainer.size();
 for(int i=0;i<size;++i)
   if(get(i)->equal(toreplace)) 
     {
       mContainer[i] = with;
     }
 return 0;
}

/************* mUPWARDAdjacencyCONTAINER *******************/


mEntity*  mUpwardAdjacencyContainer::get(int i) const
{
#ifdef _DEBUG_
  throw new mException(__LINE__,__FILE__,"no random iterator in a set<>");
#endif
  iter it=begin();
  iter End =end();
  for(int j=0;it!=End; ++it, ++j)
    if(i==j) return(*it); 
return (mEntity* ) 0;
}

mEntity* mUpwardAdjacencyContainer::find(mEntity* ent) const 
{
  iter it=mContainer.find(ent);
  if(it != end ())return *it;
  return((mEntity*)0);
}

mEntity *mUpwardAdjacencyContainer::replace(mEntity *toreplace, mEntity* with) 
{
  mEntity *other = find(toreplace);
  if(!other)return other;
  del(toreplace);
  add(with);
  return 0;
}

/************* mMESHEntityCONTAINER *******************/
mMeshEntityContainer::mMeshEntityContainer()
{
  //for(int i=0;i<4;i++)mEntities[i] = new hash_set<mEntity*,EntityHashKey,EntityEqualKey>;
  for(int i=0;i<4;i++)mEntities[i] = new set<mEntity*,EntityLessThanKey>;
}
mMeshEntityContainer::~mMeshEntityContainer()
{
  for(int i=0;i<4;i++)delete mEntities[i];
}

void mMeshEntityContainer::add(mEntity* e)
{ 
  mEntities[e->getLevel()]->insert(e);
}

void mMeshEntityContainer::del(mEntity* e)
{ 
  mEntities[e->getLevel()]->erase(e); 
}

mEntity* mMeshEntityContainer::find (int what, mDownwardAdjacencyContainer* with) const 
{
  /*  mEntity ent(with,0);
  iter it=mEntities[what]->find(&ent);
  if(it!=end(what))
    return(*it);
    else*/
    return(mEntity*)0;
}

mEntity* mMeshEntityContainer::find (mEntity *e) const 
{
  iter it=mEntities[e->getLevel()]->find(e);
  if(it!=end(e->getLevel()))
    return(*it);
  else
    return(mEntity*)0;
}

/* is not used in a newer version*/
void mMeshEntityContainer::resize (int what, int size)
{
//   mEntities[what]->resize(size);
}

