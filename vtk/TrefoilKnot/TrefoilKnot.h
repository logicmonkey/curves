// .NAME myParametricTrefoilKnot - Generate a trefoil knot
// .SECTION Description
// myParametricTrefoilKnot generates a knot.
//
#ifndef __myParametricTrefoilKnot_h
#define __myParametricTrefoilKnot_h

#include <vtkVersion.h>
#include "vtkParametricFunction.h"

class myParametricTrefoilKnot : public vtkParametricFunction
{

public:
#if VTK_MAJOR_VERSION <= 5
  vtkTypeRevisionMacro(myParametricTrefoilKnot,vtkParametricFunction);
#else
  vtkTypeMacro(myParametricTrefoilKnot,vtkParametricFunction);
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
  // Radius = 1
  static myParametricTrefoilKnot *New();

  // Description:
  // Set/Get the radius from the center to the middle of the knot.
  // The default value is 1.0.
  vtkSetMacro(Radius,double);
  vtkGetMacro(Radius,double);

  // Description
  // Return the parametric dimension of the class.
  virtual int GetDimension() {return 2;}

  // Description:
  // A Trefoil Knot
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
  myParametricTrefoilKnot();
  ~myParametricTrefoilKnot();

  // Variables
  double Radius;

private:
  myParametricTrefoilKnot(const myParametricTrefoilKnot&);  // Not implemented.
  void operator=(const myParametricTrefoilKnot&);  // Not implemented.
};

#endif
