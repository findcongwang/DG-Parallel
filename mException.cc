// mException.cpp: implementation of the mException class.
//
//////////////////////////////////////////////////////////////////////
#include <string.h>
#include <stdio.h>
#include "mException.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

mException::mException(int i, char *f, char *text)
: theLine(i)
{
   strcpy(theFile,f);
   printf("an exception has been thrown at line %d in file %s\n",i,f);
   printf("TEXT : %s\n",text);

}


