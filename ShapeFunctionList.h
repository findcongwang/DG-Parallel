#ifndef H_SHAPEFCTLIST
#define H_SHAPEFCTLIST

class Edge;
class Face;
class mRegion;

const int SHF_MAXORDER = 20;

class SHF
{
  double PascalScalarProduct (int i, int j);
  double orthogonalizationMatrix[SHF_MAXORDER][SHF_MAXORDER];  // max order 20
  void computeOrthogonalizationMatrix(int order);
 public:
  double getOrthElem(int,int);
  typedef enum _SHF_T {SHF_Pascal,SHF_Fourier,SHF_Wavelet,SHF_Legendre,SHF_PascalXY,SHF_Pascal_Orthogonal} SHF_T;
  static int getSize(int order, int level);
  static int getQSize(int order, int level);
  static void init(SHF_T);
// 1D
  static double getFct(int i, double u, Edge *e);
  static double getdFctdu(int i, double u, Edge *e);
  static double getdFctduu(int i, double u, Edge *e);
// Tet
  static double getTetFct(int i, double u, double v, double w );
  static double getTetdFctdu(int i, double u, double v, double w);
  static double getTetdFctdv(int i, double u, double v, double w);
  static double getTetdFctdw(int i, double u, double v, double w);
// Hex
  static double getHexFct(int i, double u, double v, double w );
  static double getHexdFctdu(int i, double u, double v, double w);
  static double getHexdFctdv(int i, double u, double v, double w);
  static double getHexdFctdw(int i, double u, double v, double w);
// Tet Orthogonal
  static double getTetFctOrthogonal(int i, double u, double v, double w );
  static double getTetdFctduOrthogonal(int i, double u, double v, double w);
  static double getTetdFctdvOrthogonal(int i, double u, double v, double w);
  static double getTetdFctdwOrthogonal(int i, double u, double v, double w);
// Triangle
  static double getFct(int i, double u, double v, Face * );
  static double getdFctdu(int i, double u, double v, Face *);
  static double getdFctdv(int i, double u, double v, Face *);
  static double getdFctduu(int i, double u, double v, Face *);
  static double getdFctdvv(int i, double u, double v, Face *);
  static double getdFctduv(int i, double u, double v, Face *);
  static double getdFctGeneral(int i, int derivOrder, int jth , double u, double v, Face *);
// QUADS...
  static double getQFct(int i, double u, double v, Face * );
  static double getQdFctdu(int i, double u, double v, Face *);
  static double getQdFctdv(int i, double u, double v, Face *);
  static double getQdFctduu(int i, double u, double v, Face *);
  static double getQdFctdvv(int i, double u, double v, Face *);
  static double getQdFctduv(int i, double u, double v, Face *);
  static double getQdFctGeneral(int i, int derivOrder, int jth , double u, double v, Face *);
 
 // 1D
  static double PascalgetFct(int i, double u);
  static double PascalgetdFctdu(int i, double u);
  static double PascalgetdFctduu(int i, double u);
// Triangle
  static double PascalgetFct(int i, double u, double v);
  static double PascalgetdFctdu(int i, double u, double v);
  static double PascalgetdFctdv(int i, double u, double v);
  static double PascalgetdFctduu(int i, double u, double v);
  static double PascalgetdFctdvv(int i, double u, double v);
  static double PascalgetdFctduv(int i, double u, double v);
  static double PascalgetdFctGeneral(int iFct, int ithDeriv, int jthTerm, double u, double v);

  static double PascalgetFctOrthogonal(int i, double u, double v);
  static double PascalgetdFctduOrthogonal(int i, double u, double v);
  static double PascalgetdFctdvOrthogonal(int i, double u, double v);
  static double PascalgetdFctGeneralOrthogonal(int iFct, int ithDeriv, int jthTerm, double u, double v);

// Tetraedron
  static double PascalgetFct(int i, double u, double v, double w);
  static double PascalgetdFctduu(int i, double u, double v, double w);
  static double PascalgetdFctdvv(int i, double u, double v, double w);
  static double PascalgetdFctdww(int i, double u, double v, double w);
  static double PascalgetdFctduv(int i, double u, double v, double w);
  static double PascalgetdFctduw(int i, double u, double v, double w);
  static double PascalgetdFctdvw(int i, double u, double v, double w);
protected:
  static SHF_T FunctionSpace;
	// 1D
  static double FouriergetFct(int i, double u);
  static double FouriergetdFctdu(int i, double u);
  static double FouriergetdFctduu(int i, double u);
};
#endif

