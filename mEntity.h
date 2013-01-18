// mEntity.h: interface for the mEntity class.
#ifndef _MENTITY_H_
#define _MENTITY_H_

#include "mAttachableDataContainer.h"
#include <list>
#include <vector>
#include <set>
#include <stdio.h>

class gEntity;
class Vertex;
class mMesh;
class mEntity;

using namespace std;

class EntityHashKey 
{
 public:
 inline int operator()(mEntity* ent1) const;
};
class EntityEqualKey 
{
public:
	 inline bool operator()(mEntity* ent1, mEntity* ent2) const;
};

class EntityLessThanKey 
{
 public:
	 inline bool operator()(mEntity* ent1, mEntity* ent2) const; 
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
 public:
  typedef vector<mEntity*>::const_iterator iter;
  mDownwardAdjacencyContainer(int i=1);
  mEntity* operator [] (int i) {return(mContainer[i]);}
  void add(mEntity* e){mContainer.push_back(e);}
  void del(mEntity* e);//undefined in this case
  mEntity* get(int i) const {return(mContainer[i]);}
  int size() const {return (int) mContainer.size();}
  mEntity* find(mEntity*) const;
  virtual mEntity* replace(mEntity *,mEntity*);
};

class mUpwardAdjacencyContainer : public mAdjacencyContainer 
{
 private:
  set<mEntity*,EntityLessThanKey> mContainer;
 public:
  mUpwardAdjacencyContainer() {}
  typedef set<mEntity*,EntityLessThanKey>::const_iterator iter;
  iter begin() const {return mContainer.begin();}
  iter end() const {return mContainer.end();}
  void add(mEntity* e) {mContainer.insert(e);}
  void del(mEntity* e) {mContainer.erase(e);}
  int  size() const {return (int) mContainer.size();}
  mEntity*  find(mEntity*) const ;
  mEntity* get(int i) const; 
  mEntity* operator [] (int i) {return get(i);}
  virtual mEntity* replace(mEntity *,mEntity*);
};

class mMeshEntityContainer
{
 private:
//  hash_set<mEntity*,EntityHashKey,EntityEqualKey> *mEntities[4];
  set<mEntity*,EntityLessThanKey> *mEntities[4];
 public:
//  typedef hash_set<mEntity*,EntityHashKey,EntityEqualKey>::const_iterator iter;
   typedef set<mEntity*,EntityLessThanKey>::const_iterator iter;
   mMeshEntityContainer();
   virtual ~mMeshEntityContainer();
   mEntity* find (int what, mDownwardAdjacencyContainer* with) const;
   iter begin(int what) const { return mEntities[what]->begin();}
   iter end(int what) const {return mEntities[what]->end();}
   void add(mEntity* e);
   void del(mEntity* e);
   int  size(int what) const { return (int) mEntities[what]->size();}
   void resize(int what, int size);
   mEntity* find(mEntity*) const; 
};


class mEntity  //Note: geomOrder could perhaps be of use for vertices by representing whether it is
{			   //a primary (corner), secondary (on edge), tertiary (contained within area), quaternary
	           //(contained within volume) vertex. This suggestion doesn't have any purpose as of now
	           //but maybe later..
	public:
  typedef enum mType {VERTEX,EDGE,TRI,QUAD,TET,HEX,WEDGE,PYRAMID};
 protected:
  //Represents the geometric order of the entity (presently not of use for vertices)
  int geomOrder; 
  // The iD is a function  of Vertices Id's
  unsigned long iD;
  mAttachableData* cell;
  mAdjacencyContainer *theAdjacencies[4];
  // Classification of the mesh entity to a model entity
  gEntity *theClassification; 
  //Attachable datas
  mAttachableDataContainer *theAttachable;
  mEntity();
  void computeId();
 
  public:
	  mAttachableData* getCell(){return cell;}
	  void setCell(mAttachableData* c){cell=c;}

  mEntity(mDownwardAdjacencyContainer *lis, gEntity *classification, int level);
  virtual ~mEntity();
  
  int getId() const {return iD;}
  // how many sub entities of level "what"
  virtual int getNbTemplates (int what) const;

  // create the ith sub entity of level "what" using "with"
  virtual mEntity* getTemplate (int ith , int what , int with) const; 

  // gives the level of the entity
  virtual int getLevel() const =0;
  //gives the type of the entity
  virtual mType getType() const =0;
  //gives the geometric order of the entity
  int getgeomOrder() const {return geomOrder;}

  // Take care of Upward Entities
  void add (mEntity* m);
  void del (mEntity* m);
  void replace (mEntity* from, mEntity *with);
  // we give the flexibility to find downward entities
  mEntity *find(mEntity* m)const;
  int size (int what)const
	  {
  if(!theAdjacencies[what]) return getNbTemplates(what);
  return theAdjacencies[what]->size();
		}
  mUpwardAdjacencyContainer::iter begin(int what)
  { 
	  mUpwardAdjacencyContainer *c = (mUpwardAdjacencyContainer*)theAdjacencies[what];
	  return c->begin();
  }
  mUpwardAdjacencyContainer::iter end(int what) 
  {
   mUpwardAdjacencyContainer *c = (mUpwardAdjacencyContainer*)theAdjacencies[what];
  return c->end();
	}

  // Get the ith downward entity of level what
  // if what is upward, it will give an error message
  // in this case, use ITERATORS
  mEntity* get(int what, int i)const
	{
  if(!theAdjacencies[what]) return getTemplate(i,what,0);
  return (*theAdjacencies[what])[i];
	}
  mAdjacencyContainer* get(int i) const
  {
	  return theAdjacencies[i];
  }
  void classify(gEntity *);
  void deleteAdjacencies(int what);
  gEntity* getClassification() const {return theClassification;}
  
  // Attachable Datas
  void attachData(int c , mAttachableData *v);
  void deleteData(int c);
  mAttachableData *getData(int c){
  if(!theAttachable)return 0;
  return theAttachable->get(c);
  }
  mMesh *getOwner();
  void setOwner (mMesh*);
  void attachInt(int, int);
  int getAttachedInt(int c) 
    {
      mAttachableInt *ai = (mAttachableInt *)getData(c);
      if(!ai)return 0;
      return ai->i;
    }
  mEntity* getAttachedEntity(int c) 
    {
      mAttachableEntity *ae = (mAttachableEntity *)getData(c);
      if(!ae)return 0;
      else return ae->e;
    }
  void attachDouble(int , double);
  double getAttachedDouble(int);

  bool equal    (mEntity *) const;
  bool lessthan (mEntity *) const;

  int createAdjacency(int what);
  int isAdjacencyCreated(int what) const
	{return (theAdjacencies[what])?1:0;}

  void getVertices (set<Vertex*, EntityLessThanKey> &vert) const;
  void getVertices(vector<Vertex*> &vert) const; //overloaded for vector
  // octree functions
  void getLeaves(list<mEntity*> &leaves);
  void getAllSubTree(list<mEntity*> &family);
  void setParent(mEntity*);
  void deleteParent();
  mEntity *parent();

  // debug functions
  virtual void print() const {printf("Id = %d\n",iD);}
  void printAllAttachable () const;
 };


 inline int EntityHashKey::operator()(mEntity* ent1) const 
{
  return(ent1->getId()); 
}

inline bool EntityEqualKey::operator()(mEntity* ent1,mEntity* ent2) const 
{
  return ent1->equal(ent2);
}

inline bool EntityLessThanKey::operator()(mEntity* ent1, mEntity* ent2) const 
{
  return ent1->lessthan(ent2);
}
#endif 



