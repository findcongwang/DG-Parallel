// Vertex.cpp: implementation of the Vertex class.
//
//////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include "Vertex.h"
#include "mException.h"
#include "mVector.h"
#include "mIdGenerator.h"
#include "mEntity.h"
#include <stdlib.h>
Vertex::Vertex(mIdGenerator &theIdGenerator, const mPoint &pt, gEntity *classif)
  :mEntity(),p(pt)
{
  iD = theIdGenerator.generateId();
  theClassification = classif;
}


Vertex::Vertex(int theId, const mPoint &pt, gEntity *classif)
  :mEntity(),p(pt)
{
  iD = theId;
  theClassification = classif;
}

unsigned long Vertex::getRAND()
{
  srand(iD);
  unsigned long RAND = rand();
  return RAND; 
}

void Vertex::deleteId(mIdGenerator &theIdGenerator)
{
  theIdGenerator.addId(iD);
  iD = 0;
}

void Vertex :: print() const
{
  printf("vertex Id %d = (%12.5E,%12.5E,%12.5E)\n",getId(),p(0),p(1),p(2));
}

