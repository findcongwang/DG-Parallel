#ifndef _MIDGENERATOR_H
#define _MIDGENERATOR_H

#include <stack>
#include <queue>
//#include "mCompiler.h"
using namespace std;

class mIdGenerator {
 private :
   unsigned long int maximumValue;
   queue<unsigned long int> Q;
 public:
  mIdGenerator();
  void addId(unsigned long int Id);
  unsigned long int generateId();
  unsigned long int getMaxValue();
  void setMaxValue(unsigned long int);
};

#endif


