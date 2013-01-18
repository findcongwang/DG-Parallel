#ifndef _DGLIMITER_H_
#define _DGLIMITER_H_

class DGCell;

class DGLimiter 
{
  protected:
   void computeMinMaxEdge(DGCell*);
   void computeMinMaxVertex(DGCell*);
  public:
   virtual void limit (DGCell *,double) = 0;
};

class DGExpFilterLimiter : public DGLimiter
{
  public:
   virtual void limit (DGCell *,double);
};

class DGBarthLimiter : public DGLimiter
{
  public:
   virtual void limit (DGCell *,double);
};

class DGBarthLimiterEuler : public DGLimiter
{
  public:
   virtual void limit (DGCell *, double);
};

class DGVertexLimiter : public DGLimiter
{
  public:
   virtual void limit (DGCell *,double);
};

class VertLimiter : public DGLimiter
{
  public:
   virtual void limit (DGCell *,double);
};

class DGVertexLimiterEuler : public DGLimiter
{
  public:
   virtual void limit (DGCell *,double);
};

class DGPhysicalLimiter : public DGLimiter
{
  public:
   virtual void limit (DGCell *,double);
};

class DGSuperBee : public DGLimiter
{
  public:
   virtual void limit (DGCell *,double);
};

class DGMomentLimiter : public DGLimiter
{
  public:
   virtual void limit (DGCell *,double);
};

class DGMomentLimiterEuler : public DGLimiter
{
  public:
   virtual void limit (DGCell *,double);
};


#endif
