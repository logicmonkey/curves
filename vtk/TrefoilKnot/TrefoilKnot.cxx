#include "TrefoilKnot.h"
#include "vtkObjectFactory.h"
#include "vtkMath.h"

vtkStandardNewMacro(myParametricTrefoilKnot);

//----------------------------------------------------------------------------
myParametricTrefoilKnot::myParametricTrefoilKnot() :
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
myParametricTrefoilKnot::~myParametricTrefoilKnot()
{
}

//----------------------------------------------------------------------------
void myParametricTrefoilKnot::Evaluate(double uvw[3], double Pt[3], double Duvw[9])
{
  double u = uvw[0];
  double v = uvw[1];
  double *Du = Duvw;
  double *Dv = Duvw + 3;

  // Precompute the terms.
  double cu  = cos(u);
  double c2u = cos(2*u);
  double c3u = cos(3*u);
  double su  = sin(u);
  double s2u = sin(2*u);
  double s3u = sin(3*u);

  double cv = cos(v);
  double sv = sin(v);

  // The point
  Pt[0] = (6*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*s3u*sv)/(sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*pow(s3u,2)+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*pow(s3u,2)+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2)))+((4*c2u*c3u-8*s2u*(s3u/3+1))*cv)/sqrt(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2))+4*s2u*(s3u/3+1);
  Pt[1] = -(6*(4*c2u*c3u-8*s2u*(s3u/3+1))*s3u*sv)/(sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*pow(s3u,2)+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*pow(s3u,2)+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2)))+((-8*c2u*(s3u/3+1)-4*s2u*c3u)*cv)/sqrt(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2))+4*c2u*(s3u/3+1);
  Pt[2] = (((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2))*sv)/(sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*pow(s3u,2)+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*pow(s3u,2)+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2)))+2*c3u;

  //The derivatives are:
  Du[0] = 6*s3u*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)*sv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2))+18*c3u*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*sv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2))-3*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*s3u*(2*((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2))*((8*c2u*(s3u/3+1)+4*s2u*c3u)*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)+(-8*c2u*(s3u/3+1)-4*s2u*c3u)*(-12*s2u*s3u-16*s2u*(s3u/3+1)+16*c2u*c3u)-2*(4*c2u*c3u-8*s2u*(s3u/3+1))*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u))+72*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*s3u*s3u*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)+72*(4*c2u*c3u-8*s2u*(s3u/3+1))*s3u*s3u*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u)+216*c3u*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u+216*c3u*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u)*sv/pow(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2),3/2)-(4*c2u*c3u-8*s2u*(s3u/3+1))*(2*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)+2*(4*c2u*c3u-8*s2u*(s3u/3+1))*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u))*cv/(2*pow(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2),3/2))+(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u)*cv/sqrt(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2))+8*c2u*(s3u/3+1)+4*s2u*c3u;
  Du[1] = -6*s3u*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u)*sv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2))-18*c3u*(4*c2u*c3u-8*s2u*(s3u/3+1))*sv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2))+3*(4*c2u*c3u-8*s2u*(s3u/3+1))*s3u*(2*((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2))*((8*c2u*(s3u/3+1)+4*s2u*c3u)*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)+(-8*c2u*(s3u/3+1)-4*s2u*c3u)*(-12*s2u*s3u-16*s2u*(s3u/3+1)+16*c2u*c3u)-2*(4*c2u*c3u-8*s2u*(s3u/3+1))*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u))+72*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*s3u*s3u*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)+72*(4*c2u*c3u-8*s2u*(s3u/3+1))*s3u*s3u*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u)+216*c3u*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u+216*c3u*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u)*sv/pow(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2),3/2)-(-8*c2u*(s3u/3+1)-4*s2u*c3u)*(2*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)+2*(4*c2u*c3u-8*s2u*(s3u/3+1))*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u))*cv/(2*pow(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2),3/2))+(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)*cv/sqrt(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2))-8*s2u*(s3u/3+1)+4*c2u*c3u;
  Du[2] = ((8*c2u*(s3u/3+1)+4*s2u*c3u)*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)+(-8*c2u*(s3u/3+1)-4*s2u*c3u)*(-12*s2u*s3u-16*s2u*(s3u/3+1)+16*c2u*c3u)-2*(4*c2u*c3u-8*s2u*(s3u/3+1))*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u))*sv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2))-((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2))*(2*((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2))*((8*c2u*(s3u/3+1)+4*s2u*c3u)*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)+(-8*c2u*(s3u/3+1)-4*s2u*c3u)*(-12*s2u*s3u-16*s2u*(s3u/3+1)+16*c2u*c3u)-2*(4*c2u*c3u-8*s2u*(s3u/3+1))*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u))+72*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*s3u*s3u*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)+72*(4*c2u*c3u-8*s2u*(s3u/3+1))*s3u*s3u*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u)+216*c3u*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u+216*c3u*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u)*sv/(2*pow(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2),3/2))-6*s3u;
  Dv[0] = 6*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*s3u*cv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2))-(4*c2u*c3u-8*s2u*(s3u/3+1))*sv/sqrt(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2));
  Dv[1] = -(-8*c2u*(s3u/3+1)-4*s2u*c3u)*sv/sqrt(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2))-6*(4*c2u*c3u-8*s2u*(s3u/3+1))*s3u*cv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2));
  Dv[2] = ((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2))*cv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2));
}

//----------------------------------------------------------------------------
double myParametricTrefoilKnot::EvaluateScalar(double* vtkNotUsed(uv[3]),
                                          double* vtkNotUsed(Pt[3]),
                                          double* vtkNotUsed(Duv[9]))
{
  return 0;
}

//----------------------------------------------------------------------------
void myParametricTrefoilKnot::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Radius: " << this->Radius << "\n";
}
