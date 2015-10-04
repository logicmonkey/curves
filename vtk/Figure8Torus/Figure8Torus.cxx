#include "Figure8Torus.h"
#include "vtkObjectFactory.h"
#include "vtkMath.h"

vtkStandardNewMacro(myParametricFigure8Torus);

//----------------------------------------------------------------------------
myParametricFigure8Torus::myParametricFigure8Torus() :
  Radius(1.0)
{
  this->MinimumU = -vtkMath::Pi();
  this->MinimumV = -vtkMath::Pi();
  this->MaximumU = vtkMath::Pi();
  this->MaximumV = vtkMath::Pi();

  this->JoinU = 1;
  this->JoinV = 1;
  this->TwistU = 0;
  this->TwistV = 0;
  this->ClockwiseOrdering = 1;
  this->DerivativesAvailable = 1;
}

//----------------------------------------------------------------------------
myParametricFigure8Torus::~myParametricFigure8Torus()
{
}

//----------------------------------------------------------------------------
void myParametricFigure8Torus::Evaluate(double uvw[3], double Pt[3], double Duvw[9])
{
  double u = uvw[0];
  double v = uvw[1];
  double *Du = Duvw;
  double *Dv = Duvw + 3;

  // Precompute the terms.
  double cu = cos(u);
  double su = sin(u);
  double cv = cos(v);
  double c2v = cos(2*v);
  double sv = sin(v);
  double s2v = sin(2*v);
  double term2 = sv*cu - s2v*su/2;
  double term1 = this->Radius + term2;
  double term3 = cv*cu - c2v*su;
  double term4 = -sv*su-sin(2*v)*cu/2;

  // The point
  Pt[0] = cu*term1; //cos(u)*(this->RingRadius + sin(v)*cos(u) - sin(2*v)*sin(u)/2);
  Pt[1] = su*term1;//sin(u)*(this->RingRadius + sin(v)*cos(u) - sin(2*v)*sin(u)/2);
  Pt[2] = su*sv + cu*s2v/2; //sin(u)*sin(v) + cos(u)*sin(2*v)/2;

  //The derivatives are:
  Du[0] = -Pt[1] + cu*term4;
  Du[1] = Pt[0] + su*term4;
  Du[2] = term2;
  Dv[0] = cu*term3;
  Dv[1] = su*term3;
  Dv[2] = cv*su+c2v*cu;
}

//----------------------------------------------------------------------------
double myParametricFigure8Torus::EvaluateScalar(double* vtkNotUsed(uv[3]),
                                          double* vtkNotUsed(Pt[3]),
                                          double* vtkNotUsed(Duv[9]))
{
  return 0;
}

//----------------------------------------------------------------------------
void myParametricFigure8Torus::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Radius: " << this->Radius << "\n";
}
