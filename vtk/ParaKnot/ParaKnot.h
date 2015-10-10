// .NAME myParametricParaKnot - Generate a trefoil knot
// .SECTION Description
// myParametricParaKnot generates a knot.
//
#ifndef __myParametricParaKnot_h
#define __myParametricParaKnot_h

#include <vtkVersion.h>
#include "vtkParametricFunction.h"

class myParametricParaKnot : public vtkParametricFunction
{

public:
#if VTK_MAJOR_VERSION <= 5
  vtkTypeRevisionMacro(myParametricParaKnot,vtkParametricFunction);
#else
  vtkTypeMacro(myParametricParaKnot,vtkParametricFunction);
#endif

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct a trefoil knot with the following parameters:
  // MinimumU = -Pi, MaximumU = Pi,
  // MinimumV = =Pi, MaximumV = Pi,
  // JoinU = 1, JoinV = 1,
  // TwistU = 0, TwistV = 0,
  // ClockwiseOrdering = 1,
  // DerivativesAvailable = 1,
  static myParametricParaKnot *New();

  vtkSetMacro(SpeedP,double);
  vtkSetMacro(SpeedQ,double);
  vtkSetMacro(MultP,double);
  vtkSetMacro(MultQ,double);
  vtkSetMacro(WaveAmplitude,double);
  vtkSetMacro(TubeRadius,double);

  vtkGetMacro(SpeedP,double);
  vtkGetMacro(SpeedQ,double);
  vtkGetMacro(MultP,double);
  vtkGetMacro(MultQ,double);
  vtkGetMacro(WaveAmplitude,double);
  vtkGetMacro(TubeRadius,double);

  // Description
  // Return the parametric dimension of the class.
  virtual int GetDimension() {return 2;}

  // Description:
  // A Para Knot
  //
  // This function performs the mapping \f$ f(u,v) \rightarrow (x,y,x) \f$, returning it
  // as Pt. It also returns the partial derivatives Du and Dv.
  // \f$ Pt = (x, y, z), Du = (dx/du, dy/du, dz/du), Dv = (dx/dv, dy/dv, dz/dv) \f$.
  // Then the normal is \f$ N = Du X Dv \f$.
  virtual void Evaluate(double uvw[3], double Pt[3], double Duvw[9]);

  // Description:
  // Calculate a user defined scalar using one or all of uvw,Pt,Duvw.
  //
  // uvw are the parameters with Pt being the the Cartesian point,
  // Duvw are the derivatives of this point with respect to u, v and w.
  // Pt, Duvw are obtained from Evaluate().
  //
  // This function is only called if the ScalarMode has the value
  // vtkParametricFunctionSource::SCALAR_FUNCTION_DEFINED
  //
  // If the user does not need to calculate a scalar, then the
  // instantiated function should return zero.
  //
  virtual double EvaluateScalar(double uvw[3], double Pt[3], double Duvw[9]);

protected:
  myParametricParaKnot();
  ~myParametricParaKnot();

  // Variables
  double SpeedP;
  double SpeedQ;
  double MultP;
  double MultQ;
  double WaveAmplitude;
  double TubeRadius;

private:
  myParametricParaKnot(const myParametricParaKnot&);  // Not implemented.
  void operator=(const myParametricParaKnot&);  // Not implemented.
};

#endif
