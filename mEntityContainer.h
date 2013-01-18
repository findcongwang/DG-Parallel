#ifndef _MENTITY_CONTAINER_H_
#define _MENTITY_CONTAINER_H_

#ifdef WIN32
#include <vector>
#include <set>
#include <hash_set>
#else
#include <vector>
#include <set>
#include <hash_set.h>
#endif 
using namespace std;
class mEntity;

class EntityHashKey 
{
 public:
	 int operator()(mEntity* ent1) const;
};
class EntityEqualKey 
{
public:
	 bool operator()(mEntity* ent1, mEntity* ent2) const;
};

class EntityLessThanKey 
{
 public:
	 bool operator()(mEntity* ent1, mEntity* ent2) const; 
};

class mAdjacencyContainer {
 public:
  virtual mEntity* operator [] (int i)=0;
  virtual void add(mEntity*) =0;
  virtual void del(mEntity*) =0;
  virtual mEntity* get(int) const  =0;
  virtual int size() const =0;
  virtual mEntity* find(mEntity*) const =0;
  virtual mEntity* replace(mEntity *,mEntity*)=0;
};

class mDownwardAdjacencyContainer : public mAdjacencyContainer 
{
 private:
  vector<mEntity*> mContainer;
  int mContainerSize;
 public:
  typedef vector<mEntity*>::const_iterator iter;
  mDownwardAdjacencyContainer(int i=1);
  mEntity* operator [] (int i) {return(mContainer[i]);}
  void add(mEntity* e){mContainer.push_back(e);mContainerSize++;}
  void del(mEntity* e);//undefined in this case
  mEntity* get(int i) const {return(mContainer[i]);}
  int size() const {return mContainerSize;}
  mEntity* find(mEntity*) const;
  virtual mEntity* replace(mEntity *,mEntity*);
};

class mUpwardAdjacencyContainer : public mAdjacencyContainer 
{
 private:
  set<mEntity*,EntityLessThanKey> mContainer;
  int mContainerSize;
 public:
	 mUpwardAdjacencyContainer() {mContainerSize = 0;}
  typedef set<mEntity*,EntityLessThanKey>::const_iterator iter;
  iter begin() const {return mContainer.begin();}
  iter end() const {return mContainer.end();}
  void add(mEntity* e) {mContainer.insert(e);mContainerSize++;}
  void del(mEntity* e) {mContainer.erase(e);mContainerSize--;}
  int  size() const {return mContainerSize;}
  mEntity*  find(mEntity*) const ;
  mEntity* get(int i) const; 
  mEntity* operator [] (int i) {return get(i);}
  virtual mEntity* replace(mEntity *,mEntity*);
};

class mMeshEntityContainer
{
	//can't inline "add" etc. 05/02/03
 private:
  hash_set<mEntity*,EntityHashKey,EntityEqualKey> *mEntities[4];
  // set<mEntity*,EntityEqualKey> *mEntities[4];
 public:
  typedef hash_set<mEntity*,EntityHashKey,EntityEqualKey>::const_iterator iter;
  //typedef set<mEntity*,EntityHashKey,EntityEqualKey>::const_iterator iter;
  //typedef set<mEntity*,EntityEqualKey>::const_iterator iter;
   mMeshEntityContainer();
   virtual ~mMeshEntityContainer();
   mEntity* find (int what, mDownwardAdjacencyContainer* with) const;
   iter begin(int what) const { return mEntities[what]->begin();}
   iter end(int what) const {return mEntities[what]->end();}
   void add(mEntity* e);
   void del(mEntity* e);
   int  size(int what) const { return mEntities[what]->size();}
   void resize(int what, int size);
   mEntity* find(mEntity*) const; 
};

#endif









