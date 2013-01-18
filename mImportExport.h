// mImportExport.h: interface for the mImportExport class.

#ifndef _MIMPORTEXPORT_H_
#define _MIMPORTEXPORT_H_
#include <iostream>
class mMesh;

using namespace std;

int decodeGeomOrder (int geomCode);

class mImportExport  
{
 protected:
 public:
#ifdef PARALLEL
  // import in parallel, same structure as rpm but much more simple
  void importGmshParallel (char *probName, // the .partition file 
			   mMesh *theMesh); // the local mesh for the proc
  void exportGmshParallel (char *probName, // the .partition file 
			   mMesh *theMesh, // the local mesh for the proc
			   char *directory);// directory where to put the files 
#endif
  //IMPORT
  void import (char *, mMesh *);
 // specific to gmsh format, uses streams
  void importGmshFile  (istream &, mMesh *);
  // specific to gmsh format, uses <stdio>
  void importGmshFile  (char *, mMesh *);
  // specific to sms with <stdio>
  void importSmsFile	 (char *, mMesh *);
    void importDGFile	 (istream &, mMesh *);
  //EXPORT
  void   exportSmsFile   (char *,mMesh *);
  void   exportGmshFile  (char *,mMesh *);
  void   exportGmshFile  (ostream&,mMesh *) const;
  void   exportDGFile    (ostream&,mMesh *);

};

#endif 
