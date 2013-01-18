// mAttachableData.h: interface for the mAttachableData class.
//
//////////////////////////////////////////////////////////////////////
#ifndef _MATTACHABLEDATACONTAINER_H_
#define _MATTACHABLEDATACONTAINER_H_

#include <map>
#include <vector>
#include <iostream>
using namespace std;
class mMesh;
class mEntity;

#define TYPE_CELL 1
#define TYPE_OWNER 2
#define TYPE_RLV 3
#define TYPE_P 4
#define TYPE_NO 5

class mAttachableData
{};

class mAttachableIntVector : public mAttachableData
{
 public:
    vector<int> v;
};

class mAttachableInt : public mAttachableData
{
  public :
    int i;
};

class mAttachableDouble : public mAttachableData
{
  public :
    double d;
};

class mAttachableMesh : public mAttachableData
{
  public :
    mMesh *m;
};

class mAttachableEntity : public mAttachableData
{
  public :
    mEntity *e;
};

class mAttachableDataContainer  
{
  map<int, mAttachableData *> tab;
public:
  typedef map<int ,mAttachableData *>::const_iterator iter;
  void set(int, mAttachableData *);
  void del(int);
  mAttachableData * get(int c) {return tab[c];}
  iter begin() const {return tab.begin();};
  iter end() const {return tab.end();};	
};

#endif

