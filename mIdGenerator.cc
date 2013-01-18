#include "mIdGenerator.h"
using namespace std;

mIdGenerator::mIdGenerator(){maximumValue=0;}
void mIdGenerator::addId(unsigned long int Id){ Q.push(Id);}
unsigned long int mIdGenerator::generateId(){
  if(Q.size()==0) 
    return ++maximumValue; 
  unsigned long int k=Q.front(); 
  Q.pop();  
  return k;
}

unsigned long int mIdGenerator::getMaxValue()
{
  return maximumValue;
}

void mIdGenerator::setMaxValue(unsigned long int v)
{
  maximumValue = (maximumValue >v)?maximumValue:v;
}
