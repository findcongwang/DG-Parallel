#ifndef _POSTPRO_GMSH_
#define _POSTPRO_GMSH_

#include <fstream>
#include <vector>
#include "mCompiler.h"

class PostProGmshShape
{
 protected:
  double X[4],Y[4],Z[4];
  int NbFields;
  int meshEntity;
  vector<double> Values;
  //  double *TimeValues;
 public:
  virtual void addValue (int ,double, double, double){}
  virtual void getValue (int ,double&, double&, double&){};
  virtual void setValue (int ,double, double, double){}
  virtual void setValue (int ,double, double, double, double){}
  virtual void setValue (int ,double, double, double, double
			   ,double, double, double, double,double){}
  virtual void setValue (int ,double, double, double, double
			   ,double, double, double, double,double,
			   double,double,double){}
  virtual void print (ofstream & ofs, int) = 0; 
  virtual void print (ofstream & ofs, int q, int ITER){} 

  void setXYZ(const double &x1,const double &y1,const double &z1,
	  const double &x2,const double &y2,const double &z2,
	  const double &x3,const double &y3,const double &z3)
  {
	  X[0]=x1;	Y[0]=y1;	Z[0]=z1;
	  X[1]=x2;	Y[1]=y2;	Z[1]=z2;
	  X[2]=x3;	Y[2]=y3;	Z[2]=z3;

  }

  void setXYZ(const double &x1,const double &y1,const double &z1,
	  const double &x2,const double &y2,const double &z2,
	  const double &x3,const double &y3,const double &z3,
	  const double &x4,const double &y4,const double &z4)
  {
	  X[0]=x1;	Y[0]=y1;	Z[0]=z1;
	  X[1]=x2;	Y[1]=y2;	Z[1]=z2;
	  X[2]=x3;	Y[2]=y3;	Z[2]=z3;
	  X[3]=x4;	Y[3]=y4;	Z[3]=z4;

  }
  void setMeshEntity(int me){meshEntity = me;}
  int getMeshEntity(){return meshEntity;}
};

class PostProGmshScalarTriangle : public PostProGmshShape
{
 public:
  PostProGmshScalarTriangle (int nbFields)
    {
      Values.reserve(nbFields*3);
	  Values.resize(nbFields*3);
    }

  void addValue (int i,double v1, double v2, double v3)
    {
	  i *= 3;
      Values[i]     += v1;
      Values[i + 1] += v2;
      Values[i + 2] += v3;
    }

  void setValue (int i,double v1, double v2, double v3)
    {
	  i *= 3;
      Values[i]     = v1;
      Values[i + 1] = v2;
      Values[i + 2] = v3;
    }

  void getValue (int i,double &v1, double &v2, double &v3)
    {
	  i *= 3;
      v1 = Values[i];
      v2 = Values[i + 1];
      v3 = Values[i + 2];
    }

  void print (ofstream & ofs, int q, int ITER) 
	  {
      ofs << "ST (" 
	  << X[0] << "," << Y[0] << "," << Z[0] << ","
	  << X[1]<< ","  << Y[1] << "," << Z[1] << "," 
	  << X[2] << "," << Y[2] << "," << Z[2] << ") {";
	  q *=3;
      for (int i=0; i <3; i++)
	{
	  if(i)ofs << "," << Values[q+i]/(double) ITER;
	  else ofs << " " << Values[q+i]/(double) ITER;
	} 
      ofs << "};" << endl;
    } 

  void print  (ofstream & ofs, int q)
    {
      ofs << "ST (" 
	  << X[0] << "," << Y[0] << "," << Z[0] << ","
	  << X[1]<< ","  << Y[1] << "," << Z[1] << "," 
	  << X[2] << "," << Y[2] << "," << Z[2] << ") {";
	  q *=3;
      for (int i=0; i <3; i++)
	{
	  if(i)ofs << "," << Values[q+i];
	  else ofs << " " << Values[q+i];
	} 
      ofs << "};" << endl;
    } 
};


class PostProGmshVectorTriangle : public PostProGmshShape
{
 public:
  PostProGmshVectorTriangle (int nbFields)
    {
      Values.reserve( nbFields*9);
	  Values.resize( nbFields*9);

    }
  void setValue (int i,double v1, double v2, double v3
		     ,double v4, double v5, double v6
		     ,double v7, double v8, double v9)
    {
      i *= 9;
      Values[i]     = v1;
      Values[i + 1] = v2;
      Values[i + 2] = v3;
      Values[i + 3] = v4;
      Values[i + 4] = v5;
      Values[i + 5] = v6;
      Values[i + 6] = v7;
      Values[i + 7] = v8;
      Values[i + 8] = v9;
    }
  void print (ofstream & ofs, int q)
    {
      ofs << "VT (" 
	  << X[0] << "," << Y[0] << "," << Z[0] << ","
	  << X[1]<< ","  << Y[1] << "," << Z[1] << "," 
	  << X[2] << "," << Y[2] << "," << Z[2] << ") {";
	  
	  for (int i=0; i <9; i++)
	{
	  q *= 9;
	  if(i)ofs << ","<< Values[q+i];
	  else ofs << " " << Values[q+i];
	  // if(i%10 == 9)ofs << endl;
	} 
      ofs << "};" << endl;
    } 
};

class PostProGmshScalarTetrahedron : public PostProGmshShape
{
 public:
  PostProGmshScalarTetrahedron (int nbFields = 1)
    {
      Values.reserve(nbFields*4);
	   Values.resize(nbFields*4);
    }
   void setValue (int i,double v1, double v2, double v3, double v4)
    {
	   i *=4;
      Values[i] = v1;
      Values[i + 1] = v2;
      Values[i + 2] = v3;
      Values[i + 3] = v4;
    }
  void print (ofstream & ofs, int q)
    {
      ofs << "SS ("
	  << X[0] << "," << Y[0] << "," << Z[0] << ","
	  << X[1]<< ","  << Y[1] << "," << Z[1] << "," 
	  << X[2]<< ","  << Y[2] << "," << Z[2] << "," 
	  << X[3] << "," << Y[3] << "," << Z[3] << ") {";
	  for (int i=0; i <1*4; i++)
	{
	  if(i)
	    {
	      ofs << ", " << Values[q*4+i];
	    }
	  else
	    { 
	      ofs << " " << Values[q*4+i];
	    }	  
	  if(i%7 == 6)ofs << endl;
	} 
      ofs << "};" << endl;
    } 
};

class PostProGmshVectorTetrahedron : public PostProGmshShape
{
 public:
  PostProGmshVectorTetrahedron (int nbFields)
    {
      Values.reserve( nbFields*12);
	  Values.resize( nbFields*12);

    }
  void setValue (int i,double v1, double v2, double v3
		     ,double v4, double v5, double v6
		     ,double v7, double v8, double v9,
			  double v10, double v11, double v12)
    {
      i *= 12;
      Values[i]     = v1;
      Values[i + 1] = v2;
      Values[i + 2] = v3;
      Values[i + 3] = v4;
      Values[i + 4] = v5;
      Values[i + 5] = v6;
      Values[i + 6] = v7;
      Values[i + 7] = v8;
      Values[i + 8] = v9;
	  Values[i + 9] = v10;
	  Values[i + 10] = v11;
	  Values[i + 11] = v12;

    }

  void print (ofstream & ofs, int q)
    {
      ofs << "VS (" 
	  << X[0] << "," << Y[0] << "," << Z[0] << ","
	  << X[1]<< ","  << Y[1] << "," << Z[1] << "," 
	  << X[2] << "," << Y[2] << "," << Z[2] << ","
	  << X[3] << "," << Y[3] << "," << Z[3] << ") {";
	  
	  for (int i=0; i <12; i++)
	{
	  q *= 12;
	  if(i)ofs << ","<< Values[q+i];
	  else ofs << " " << Values[q+i];
	} 
      ofs << "};" << endl;
    } 
};
class PostProGmshView
{
	// used to collect all postprocessing data in one array for printing
 private:
  int NbShapes;  //number of triangles or tets passed to the gmsh
  int CurrentShape; //index of the current triangle or tet
 public:
  PostProGmshShape **AllShapes; 
  int getCurrentShape(){return CurrentShape;}
  int getNbShapes(){return NbShapes;}

  PostProGmshView (int N):NbShapes(N)
    {
      AllShapes = new PostProGmshShape*[NbShapes];
      CurrentShape = 0;
    }
  ~PostProGmshView ()
    {
      for(int i=0;i<CurrentShape;i++)delete AllShapes[i];
      delete AllShapes;
    }
  PostProGmshShape * getShape(int i)
    {
      return AllShapes[i];
    }
  PostProGmshShape * NextInViewIsAScalarTriangle(int Nb = 1)
    {
      AllShapes[CurrentShape] = new PostProGmshScalarTriangle(Nb);
      return AllShapes[CurrentShape++];      
    }
  PostProGmshShape * NextInViewIsAVectorTriangle(int Nb = 1)
    {
      AllShapes[CurrentShape] = new PostProGmshVectorTriangle(Nb);
      return AllShapes[CurrentShape++];      
    }
  PostProGmshShape * NextInViewIsAScalarTetrahedron(int Nb = 1)
    {
      AllShapes[CurrentShape] = new PostProGmshScalarTetrahedron(Nb);
      return AllShapes[CurrentShape++];      
    }

  PostProGmshShape * NextInViewIsAVectorTetrahedron(int Nb = 1)
    {
      AllShapes[CurrentShape] = new PostProGmshVectorTetrahedron(Nb);
      return AllShapes[CurrentShape++];      
    }

  
  friend ofstream & operator << (ofstream & ofs, const PostProGmshView &p)
    {
      ofs << "View \"Exported field\" {" << endl; 
      for(int i=0;i<p.CurrentShape;i++)p.AllShapes[i]->print(ofs,i);
      ofs << "};" << endl;
      return ofs;
    }
};

#endif
