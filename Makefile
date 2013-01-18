#
# Makefile for "libmeshLib.a"
#

.IGNORE:
#ifdef $x86_64
CPP		= g++ -O3 
#endif
#CPP              =g++ -g -O3
#CPP             = CC
#CC_FLAGS        =   -p
#CC_FLAGS       = -library=stlport4 -dalign -xO3 -p
#CC_FLAGS        =   -dalign -xtarget=ultra2 -xO3 -p
#CC_FLAGS        =  -xO3 
#CFLAGS         = -Wno-unused  -Wno-deprecated \
#	-mcpu=supersparc -mtune=ultrasparc\
#	-fomit-frame-pointer  -Winline 
#VERSION_FLAGS = 

RM       = rm
RMFLAGS  = -f

LIB      = lib/libmeshLib.a

SRC = mEntity.cc \
      Edge.cc \
      Vertex.cc \
      mIdGenerator.cc \
      Face.cc \
      mTet.cc \
      mHex.cc \
      mImportExport.cc \
      mEntityContainer.cc \
      mAttachableDataContainer.cc \
      mException.cc \
      mDGMesh.cc \
      mMesh.cc \
      mPoint.cc \
      mTensor2.cc \
      mVector.cc \
      Integrator.cc \
      EulerLaw.cc \
      ScalarLaw.cc \
      Burgers.cc\
      Burgers2D.cc\
      Suli.cc\
      Examples.cc \
      Lilia.cc\
      LinearSystem.cc\
      DGCell.cc \
      DGError.cc \
      DGMyError.cc\
      DGPref.cc \
      DGAnalysis.cc \
      DGLimiter.cc \
      DGSensor.cc \
      DGOrthogonalize.cc \
      Mapping.cc \
      FieldEvaluator.cc \
      FunctionSpace.cc \
      matrixUtil.cc \
      ShapeFunctionList.cc \
      GaussLegendre1D.cc \
      GaussLegendreSimplex.cc \
      GaussQuadratureTet.cc \
      GaussQuadratureHex.cc \
      GaussQuadratureQuad.cc \
      GaussQuadratureTri.cc  \
      TimeIntegrator.cc	\
      BasisFunctions.cc \
      CutCell.cc\
      Geometry.cc
OBJ = $(SRC:.cc=.o)

.SUFFIXES: .o .cc

$(LIB): $(OBJ) 
	 ar r $(LIB) $(OBJ) 

.cc.o:
	$(CPP) $(CC_FLAGS) -c $< 

clean:
	$(RM) $(RMFLAGS) $(OBJ)

lint:
	$(LINT) $(CC_FLAGS) $(SRC)

exe:	
	$(CPP) $(CC_FLAGS)  main.cc lib/libmeshLib.a -o bin/dg -lc -lm

depend:
	(sed '/^# DO NOT DELETE THIS LINE/q' Makefile && \
	$(CPP) -MM  $(CFLAGS) ${SRC} \
	) >Makefile.new
	cp Makefile Makefile.bak
	cp Makefile.new Makefile
	$(RM) $(RMFLAGS) Makefile.new

# DO NOT DELETE THIS LINE
mEntity.o: mEntity.cc mEntity.h mAttachableDataContainer.h mException.h \
  gEntity.h Vertex.h mPoint.h Edge.h
Edge.o: Edge.cc Edge.h mEntity.h mAttachableDataContainer.h Vertex.h \
  mPoint.h
Vertex.o: Vertex.cc Vertex.h mPoint.h mEntity.h \
  mAttachableDataContainer.h mException.h mVector.h mIdGenerator.h
mIdGenerator.o: mIdGenerator.cc mIdGenerator.h
Face.o: Face.cc Face.h mEntity.h mAttachableDataContainer.h Vertex.h \
  mPoint.h Edge.h
mTet.o: mTet.cc mTet.h mRegion.h mEntity.h mAttachableDataContainer.h \
  Edge.h Face.h mException.h Vertex.h mPoint.h
mHex.o: mHex.cc mHex.h mRegion.h mEntity.h mAttachableDataContainer.h \
  Edge.h Face.h mException.h Vertex.h mPoint.h
mImportExport.o: mImportExport.cc mImportExport.h mMesh.h mEntity.h \
  mAttachableDataContainer.h mIdGenerator.h gEntity.h Vertex.h mPoint.h \
  Edge.h Face.h mTet.h mRegion.h mException.h
mEntityContainer.o: mEntityContainer.cc mEntity.h \
  mAttachableDataContainer.h mException.h
mAttachableDataContainer.o: mAttachableDataContainer.cc \
  mAttachableDataContainer.h
mException.o: mException.cc mException.h
mDGMesh.o: mDGMesh.cc mDGMesh.h mMesh.h mEntity.h \
  mAttachableDataContainer.h mIdGenerator.h gEntity.h Vertex.h mPoint.h \
  Edge.h Face.h
mMesh.o: mMesh.cc mMesh.h mEntity.h mAttachableDataContainer.h \
  mIdGenerator.h gEntity.h Vertex.h mPoint.h Edge.h Face.h mTet.h \
  mRegion.h mHex.h mException.h
mPoint.o: mPoint.cc mPoint.h
mTensor2.o: mTensor2.cc mTensor2.h mVector.h
mVector.o: mVector.cc mVector.h mTensor2.h mPoint.h mException.h \
  Constants.h
Integrator.o: Integrator.cc Integrator.h mEntity.h \
  mAttachableDataContainer.h IntPt.h BasisFunctions.h FunctionSpace.h \
  mVector.h ShapeFunctionList.h
EulerLaw.o: EulerLaw.cc EulerLaw.h ConservationLaw.h mPoint.h mVector.h \
  FieldEvaluator.h mCompiler.h
ScalarLaw.o: ScalarLaw.cc ScalarLaw.h ConservationLaw.h mVector.h \
  mPoint.h FieldEvaluator.h mCompiler.h
Burgers.o: Burgers.cc DGAnalysis.h Integrator.h FieldEvaluator.h \
  mCompiler.h mPoint.h Burgers.h ConservationLaw.h DGLimiter.h DGSensor.h \
  mDGMesh.h mMesh.h mEntity.h mAttachableDataContainer.h mIdGenerator.h \
  gEntity.h Constants.h mVector.h
Burgers2D.o: Burgers2D.cc DGAnalysis.h Integrator.h FieldEvaluator.h \
  mCompiler.h mPoint.h Burgers.h ConservationLaw.h DGLimiter.h DGSensor.h \
  mDGMesh.h mMesh.h mEntity.h mAttachableDataContainer.h mIdGenerator.h \
  gEntity.h Constants.h mVector.h
Suli.o: Suli.cc DGAnalysis.h Integrator.h FieldEvaluator.h mCompiler.h \
  mPoint.h ScalarLaw.h ConservationLaw.h DGLimiter.h DGSensor.h mDGMesh.h \
  mMesh.h mEntity.h mAttachableDataContainer.h mIdGenerator.h gEntity.h \
  mVector.h Constants.h
Examples.o: Examples.cc DGAnalysis.h Integrator.h FieldEvaluator.h \
  mCompiler.h mPoint.h EulerLaw.h ConservationLaw.h ScalarLaw.h \
  DGLimiter.h DGSensor.h mDGMesh.h mMesh.h mEntity.h \
  mAttachableDataContainer.h mIdGenerator.h gEntity.h Constants.h
Lilia.o: Lilia.cc DGAnalysis.h Integrator.h FieldEvaluator.h mCompiler.h \
  mPoint.h ScalarLaw.h ConservationLaw.h LinearSystem.h DGLimiter.h \
  DGSensor.h mDGMesh.h mMesh.h mEntity.h mAttachableDataContainer.h \
  mIdGenerator.h gEntity.h mVector.h Constants.h
LinearSystem.o: LinearSystem.cc LinearSystem.h ConservationLaw.h \
  mVector.h mPoint.h FieldEvaluator.h mCompiler.h
DGCell.o: DGCell.cc DGCell.h mPoint.h mVector.h \
  mAttachableDataContainer.h Constants.h FunctionSpace.h \
  ShapeFunctionList.h ConservationLaw.h Integrator.h Mapping.h mTensor2.h \
  mCompiler.h mEntity.h Face.h Vertex.h FieldEvaluator.h mMesh.h \
  mIdGenerator.h gEntity.h Edge.h
DGError.o: DGError.cc DGCell.h mPoint.h mVector.h \
  mAttachableDataContainer.h Constants.h FunctionSpace.h \
  ShapeFunctionList.h ConservationLaw.h Integrator.h Mapping.h mTensor2.h \
  mCompiler.h mEntity.h gEntity.h Edge.h Face.h Vertex.h FieldEvaluator.h \
  mMesh.h mIdGenerator.h
DGMyError.o: DGMyError.cc DGAnalysis.h Integrator.h DGMyError.h Mapping.h \
  mTensor2.h mCompiler.h FunctionSpace.h mVector.h ShapeFunctionList.h \
  ConservationLaw.h mDGMesh.h mMesh.h mEntity.h \
  mAttachableDataContainer.h mIdGenerator.h gEntity.h DGLimiter.h \
  DGCell.h mPoint.h Constants.h FieldEvaluator.h mImportExport.h
DGPref.o: DGPref.cc mEntity.h mAttachableDataContainer.h mDGMesh.h \
  mMesh.h mIdGenerator.h gEntity.h DGCell.h mPoint.h mVector.h \
  Constants.h FunctionSpace.h ShapeFunctionList.h ConservationLaw.h \
  DGAnalysis.h Integrator.h mImportExport.h FieldEvaluator.h mCompiler.h \
  Face.h Edge.h Vertex.h
DGAnalysis.o: DGAnalysis.cc DGAnalysis.h Integrator.h Mapping.h \
  mTensor2.h mCompiler.h FunctionSpace.h mVector.h ShapeFunctionList.h \
  ConservationLaw.h mDGMesh.h mMesh.h mEntity.h \
  mAttachableDataContainer.h mIdGenerator.h gEntity.h DGLimiter.h \
  DGCell.h mPoint.h Constants.h FieldEvaluator.h mImportExport.h \
  TimeIntegrator.h Vertex.h Edge.h
DGLimiter.o: DGLimiter.cc mEntity.h mAttachableDataContainer.h Mapping.h \
  mTensor2.h mCompiler.h Integrator.h ConservationLaw.h FunctionSpace.h \
  mVector.h ShapeFunctionList.h DGLimiter.h DGCell.h mPoint.h Constants.h \
  FieldEvaluator.h gEntity.h
DGSensor.o: DGSensor.cc DGSensor.h mPoint.h DGCell.h mVector.h \
  mAttachableDataContainer.h Constants.h FunctionSpace.h \
  ShapeFunctionList.h ConservationLaw.h mDGMesh.h mMesh.h mEntity.h \
  mIdGenerator.h gEntity.h DGAnalysis.h Integrator.h Mapping.h mTensor2.h \
  mCompiler.h postProGmsh.h
DGOrthogonalize.o: DGOrthogonalize.cc FunctionSpace.h mVector.h \
  ShapeFunctionList.h Integrator.h mEntity.h mAttachableDataContainer.h \
  Mapping.h mTensor2.h mCompiler.h
Mapping.o: Mapping.cc Mapping.h mTensor2.h mCompiler.h mEntity.h \
  mAttachableDataContainer.h mVector.h mPoint.h Vertex.h
FieldEvaluator.o: FieldEvaluator.cc FieldEvaluator.h mCompiler.h DGCell.h \
  mPoint.h mVector.h mAttachableDataContainer.h Constants.h \
  FunctionSpace.h ShapeFunctionList.h ConservationLaw.h
FunctionSpace.o: FunctionSpace.cc Mapping.h mTensor2.h mCompiler.h \
  FunctionSpace.h mVector.h ShapeFunctionList.h
matrixUtil.o: matrixUtil.cc
ShapeFunctionList.o: ShapeFunctionList.cc mEntity.h \
  mAttachableDataContainer.h Edge.h Face.h mRegion.h ShapeFunctionList.h \
  OrthoSFTri.h TetOrtho.h
GaussLegendre1D.o: GaussLegendre1D.cc
GaussLegendreSimplex.o: GaussLegendreSimplex.cc IntPt.h \
  GaussLegendreSimplex.h
GaussQuadratureTet.o: GaussQuadratureTet.cc IntPt.h \
  GaussLegendreSimplex.h
GaussQuadratureHex.o: GaussQuadratureHex.cc IntPt.h \
  GaussLegendreSimplex.h
GaussQuadratureQuad.o: GaussQuadratureQuad.cc IntPt.h \
  GaussLegendreSimplex.h
GaussQuadratureTri.o: GaussQuadratureTri.cc IntPt.h \
  GaussLegendreSimplex.h
TimeIntegrator.o: TimeIntegrator.cc TimeIntegrator.h DGLimiter.h \
  mDGMesh.h mMesh.h mEntity.h mAttachableDataContainer.h mIdGenerator.h \
  gEntity.h DGCell.h mPoint.h mVector.h Constants.h FunctionSpace.h \
  ShapeFunctionList.h ConservationLaw.h FieldEvaluator.h mCompiler.h \
  Mapping.h mTensor2.h
BasisFunctions.o: BasisFunctions.cc BasisFunctions.h mEntity.h \
  mAttachableDataContainer.h FunctionSpace.h mVector.h \
  ShapeFunctionList.h Integrator.h
CutCell.o: CutCell.cc CutCell.h DGCell.h mPoint.h mVector.h \
  mAttachableDataContainer.h Constants.h FunctionSpace.h \
  ShapeFunctionList.h ConservationLaw.h Integrator.h Mapping.h mTensor2.h \
  mCompiler.h mEntity.h Face.h Vertex.h FieldEvaluator.h mMesh.h \
  mIdGenerator.h gEntity.h Edge.h Geometry.h
Geometry.o: Geometry.cc Geometry.h mPoint.h DGCell.h Vertex.h mEntity.h
main.o: main.cc mImportExport.h mDGMesh.h mMesh.h mEntity.h \
  mAttachableDataContainer.h mIdGenerator.h gEntity.h
