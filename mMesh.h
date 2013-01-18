// mMesh.h: interface for the mMesh class.

#ifndef _MMESH_H_
#define _MMESH_H_

#include <iostream>
#include "mEntity.h"
#include "mIdGenerator.h"
#include "gEntity.h"
#include "mAttachableDataContainer.h"

class Face;
class Edge;
class Vertex;
class mRegion;
class mTet;
class mHex;

using namespace std;

class mMesh  
{
 protected:
  // the containers
  mMeshEntityContainer allEntities;
  mIdGenerator theIdGenerator;
  set<gEntity*,gEntityLessThanKey> allGEntities;
  mMeshEntityContainer bins[10];
  mMeshEntityContainer interfaceBins[10];
  mMeshEntityContainer largeOnInterface[10];
  mMeshEntityContainer smallOnInterface[10];
  int nbBins;
 public:
  typedef mMeshEntityContainer::iter iter;
  mMesh(){nbBins=1;}
  virtual ~mMesh();
  // adding mesh entities, they are created with any other entities.

  Vertex* getVertex	(int iD);
  Edge*   getEdge (vector<Vertex*>& v, int geomO); 
  Face*   getFace (vector<Vertex*>& v, int Type, int geomO); 
  mTet*    getTet (vector<Vertex*>& v, int geomO); 
  Vertex* createVertex	(int iD, const double x, const double y, const double z, gEntity *classif);
  Vertex* createVertex	(const double x, const double y, const double z , gEntity *classif);
  Edge*   createEdge      (Vertex *v1 , Vertex *v2, gEntity *classif); //kept for split compatibility
  Edge*   createEdge      (vector<int>& vertexiD, int geomOrder, gEntity *classif);
  Face*   createFaceWithVertices (vector<int>& vertexiD, int theType, int geomOrder, gEntity *classif); //generalized
  Face*   createFaceWithEdges (Edge *e1 , Edge *e2, Edge *e3, gEntity *classif);
  Face*   createFaceWithEdges (Edge *e1 , Edge *e2, Edge *e3, Edge *e4, gEntity *classif);
  mTet*   createTetWithVertices (vector<int>& vertexiD, int geomOrder, gEntity *classif); //generalized
  mTet*   createTetWithFaces (Face *f1 , Face *f2, Face *f3, Face *f4, gEntity *classif);
  mHex*   createHexWithVertices (int iV1 , int iV2, int iV3, int iV4, 
				  int iV5 , int iV6, int iV7, int iV8, 
				  gEntity *classif);
  void createConnections (int from, int to);
  void createCellBoundaries (int i, int j, int with =0);
  void setPeriodicBC(int dim);
  iter begin(int what) { return allEntities.begin(what);}
  iter end(int what)   { return allEntities.end(what);}
  int size(int what)   { return allEntities.size(what);}
  void resize(int what, int size);

  virtual mEntity *checkExisting (mEntity *);
  virtual mEntity *find (mEntity *);
  void add (mEntity *e) { allEntities.add(e);}
  void del (mEntity *e) { allEntities.del(e);} //delete from mesh
  void DEL (mEntity *e) { del(e); delete e;}  //physically delete the element
  
  /********************Adaptive Time Steppping**********************************/

  void setNbBins(int n) {nbBins=n;}
  int getNbBins() {return nbBins;}
  iter beginBin(int binNb,int what) { return bins[binNb].begin(what);}
  iter endBin(int binNb,int what)   { return bins[binNb].end(what);}
  iter beginLargeOnInterface(int binNb,int what) { return largeOnInterface[binNb].begin(what);}
  iter endLargeOnInterface(int binNb,int what)   { return largeOnInterface[binNb].end(what);}
  iter beginSmallOnInterface(int binNb,int what) { return smallOnInterface[binNb].begin(what);}
  iter endSmallOnInterface(int binNb,int what)   { return smallOnInterface[binNb].end(what);}
  iter beginInterfaceBin(int binNb,int what) { return interfaceBins[binNb].begin(what);}
  iter endInterfaceBin(int binNb,int what)   { return interfaceBins[binNb].end(what);}
 
  void addToBin   (int i, mEntity *e) { bins[i].add(e);}
  void delFromBin (int i, mEntity *e) { bins[i].del(e);}
  int sizeOfBin  (int binNb, int levelNb) {return bins[binNb].size(levelNb);}
  
  void addToInterfaceBin   (int i, mEntity *e) { interfaceBins[i].add(e);}
  void delFromInterfaceBin (int i, mEntity *e) { interfaceBins[i].del(e);}
  int sizeOfInterfaceBin  (int binNb, int levelNb) {return interfaceBins[binNb].size(levelNb);}

  void addToLargeOnInterface   (int i, mEntity *e) { largeOnInterface[i].add(e);}
  void delFromLargeOnInterface (int i, mEntity *e) { largeOnInterface[i].del(e);}
  int sizeOfLargeOnInterface  (int binNb, int levelNb) {return largeOnInterface[binNb].size(levelNb);}

  void addToSmallOnInterface   (int i, mEntity *e) { smallOnInterface[i].add(e);}
  void delToSmallOnInterface (int i, mEntity *e) { smallOnInterface[i].del(e);}
  int sizeOfSmallOnInterface  (int binNb, int levelNb) {return smallOnInterface[binNb].size(levelNb);}
  
  /********************Input Output ********************************************/
  mEntity *readMEntity(istream &);
  void     writeMEntity(ostream &, mEntity*);

  // iterator on all entities
  // virtual iterall beginall(int what) const;
  //virtual iterall endall(int what) const;
  //virtual int size(int what) const;
  gEntity *getGEntity (int,int);
};

#endif 
