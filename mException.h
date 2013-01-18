// mException.h: interface for the mException class.
//
//////////////////////////////////////////////////////////////////////
#ifndef _MEXCEPTION_H_
#define _MEXCEPTION_H_

class mException  
{
	int theLine;
	char theFile[256];
public:
	mException(int Line,char *file, char *text); 
};

class mExceptionFileNotFound : public mException
{
public:
	mExceptionFileNotFound(int Line,char *file, char *text)
		: mException(Line,file,text)
	{}
};
#endif 

