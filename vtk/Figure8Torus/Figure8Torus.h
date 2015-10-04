// .NAME myParametricFigure8Torus - Generate a figure-8 torus.
// .SECTION Description
// myParametricFigure8Torus generates a torus.
//
#ifndef __myParametricFigure8Torus_h
#define __myParametricFigure8Torus_h

#include <vtkVersion.h>
#include "vtkParametricFunction.h"

class myParametricFigure8Torus : public vtkParametricFunction
{

public:
#if VTK_MAJOR_VERSION <= 5
  vtkTypeRevisionMacro(myParametricFigure8Torus,vtkParametricFunction);
#else
  vtkTypeMacro(myParametricFigure8Torus,vtkParametricFunction);
#endif

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct a figure-8 torus with the following parameters:
  // MinimumU = -Pi, MaximumU = Pi,
  // MinimumV = =Pi, MaximumV = Pi,
  // JoinU = 1, JoinV = 1,
  // TwistU = 0, TwistV = 0,
  // ClockwiseOrdering = 1,
  // DerivativesAvailable = 1,
  // Radius = 1
  static myParametricFigure8Torus *New();

  // Description:
  // Set/Get the radius from the center to the middle of the ring of the
  // figure-8 torus.  The default value is 1.0.
  vtkSetMacro(Radius,double);
  vtkGetMacro(Radius,double);

  // Description
  // Return the parametric dimension of the class.
  virtual int GetDimension() {return 2;}

  // Description:
  // A Figure-8 torus.
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
  myParametricFigure8Torus();
  ~myParametricFigure8Torus();

  // Variables
  double Radius;

private:
  myParametricFigure8Torus(const myParametricFigure8Torus&);  // Not implemented.
  void operator=(const myParametricFigure8Torus&);  // Not implemented.
};

#endif
