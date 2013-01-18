#ifndef _DGMESH_H_
#define _DGMESH_H_

#include <list>
#include "mMesh.h"

class mEntity;
using namespace std;
 
class mDGMesh : public mMesh
{
  Vertex* splitEdge (Edge *);
  void splitTriangle (Face *);
  void unsplitTriangle(Face *e);
  void splitQuad (Face *);
  void reconnect (int dim);
  protected:
   mMeshEntityContainer allSplittedEntities;
  public :
    iter beginsplit(int what){return allSplittedEntities.begin(what);} 
    iter endsplit(int what){return allSplittedEntities.end(what);}
    void split  (int dim, list<mEntity *> &toSplit);
	void split (mEntity* );
    void unsplit(int dim, list<mEntity *> &toUnSplit);
	void unsplit (mEntity* );
    void getLeaves(mEntity *, list<mEntity*> &leaves);
    int getRefinementDepth(mEntity *);
    int getRefinementLevel(mEntity *);
};

#endif
