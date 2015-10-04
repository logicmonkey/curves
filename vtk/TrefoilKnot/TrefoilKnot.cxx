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

  // Precompute the common trenscendentals
  double cu  = cos(u);
  double c2u = cos(2*u);
  double c3u = cos(3*u);
  double su  = sin(u);
  double s2u = sin(2*u);
  double s3u = sin(3*u);

  double cv = cos(v);
  double sv = sin(v);

  // The Trefoil Knot is generated from a simple plane curve C2 with the
  // parametrization in u:
  //
  //     C2(u): [4*sin(2*u)*(sin(3*u)/3+1),
  //             4*cos(2*u)*(sin(3*u)/3+1)]
  //
  // Adding a sinusoidal third coordinate turns the plane curve into a
  // space curve C3(u):
  //
  //     C3(u): [4*sin(2*u)*(sin(3*u)/3+1),
  //          4*cos(2*u)*(sin(3*u)/3+1),
  //          2*cos(3*u)]
  //
  // Tracing the space curve with a circle centered on, and perpendicular
  // to its tangent introduces a second parameter v and sweeps out a
  // surface S. u takes us along the space curve, v takes us to a point
  // on the surface centred around C3. The calculation steps are as
  // follows.
  //
  // The tangent vector T(x,y,z) to C3 is its first derivative.
  // The normal vector to T in the xy plane is trivially calculated as:
  //
  //     N = (T.y, -T.x, 0)
  //
  // For a consistent tube radius along C3, the surface sweeping circle
  // needs to be in a plane through N and perpendicular to T. The binormal
  // vector B to T and N is given by the cross product:
  //
  //     B = TxN
  //
  // S is then the point on C3 plus an offset on the sweeping circle of
  // radius r in 3D space:
  //
  //     S(u,v) = C3(u) + r*cos(v)*NN + r*sin(v)*BB
  //
  //     where NN is normalised unit vector N
  //     and BB is normalised binormal B
  //
  // Pt is a point on S.

  Pt[0] = (6*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*s3u*sv)/(sqrt(36*pow(4*c2u*c3u\
    -8*s2u*(s3u/3+1),2)*pow(s3u,2)+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*\
    pow(s3u,2)+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*\
    c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2)))+((4*c2u*c3u-8*s2u*(s3u/3+\
    1))*cv)/sqrt(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+pow(-8*c2u*(s3u/3+1)-\
    4*s2u*c3u,2))+4*s2u*(s3u/3+1);

  Pt[1] = -(6*(4*c2u*c3u-8*s2u*(s3u/3+1))*s3u*sv)/(sqrt(36*pow(4*c2u*c3u\
    -8*s2u*(s3u/3+1),2)*pow(s3u,2)+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*\
    pow(s3u,2)+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*\
    c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2)))+((-8*c2u*(s3u/3+1)-4*s2u*\
    c3u)*cv)/sqrt(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+pow(-8*c2u*(s3u/3+1)-\
    4*s2u*c3u,2))+4*c2u*(s3u/3+1);

  Pt[2] = (((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow\
    (4*c2u*c3u-8*s2u*(s3u/3+1),2))*sv)/(sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u\
    /3+1),2)*pow(s3u,2)+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*pow(s3u,2)+\
    pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*\
    c2u*c3u-8*s2u*(s3u/3+1),2),2)))+2*c3u;

  // Now calculate the partial derivatives with respect to u and v for
  // each component of Pt. A computer algebra system is useful here.

  // Use of the partial derivatives is a compile option.

  Du[0] = 6*s3u*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)*sv/sqrt(36*pow(\
    4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*\
    c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*\
    s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2))+18*c3u*(-8*c2u*(s3u/3+\
    1)-4*s2u*c3u)*sv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36\
    *pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*\
    s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2\
    ),2))-3*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*s3u*(2*((-8*c2u*(s3u/3+1)-4*s2u\
    *c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2))*\
    ((8*c2u*(s3u/3+1)+4*s2u*c3u)*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u\
    )+(-8*c2u*(s3u/3+1)-4*s2u*c3u)*(-12*s2u*s3u-16*s2u*(s3u/3+1)+16*c2u*\
    c3u)-2*(4*c2u*c3u-8*s2u*(s3u/3+1))*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*\
    s2u*c3u))+72*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*s3u*s3u*(12*s2u*s3u+16*s2u\
    *(s3u/3+1)-16*c2u*c3u)+72*(4*c2u*c3u-8*s2u*(s3u/3+1))*s3u*s3u*(-12*\
    c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u)+216*c3u*pow(4*c2u*c3u-8*s2u*(\
    s3u/3+1),2)*s3u+216*c3u*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u)*sv/\
    pow(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3\
    +1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(\
    s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2),3/2)-(4*c2u*\
    c3u-8*s2u*(s3u/3+1))*(2*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*(12*s2u*s3u+16*\
    s2u*(s3u/3+1)-16*c2u*c3u)+2*(4*c2u*c3u-8*s2u*(s3u/3+1))*(-12*c2u*s3u\
    -16*c2u*(s3u/3+1)-16*s2u*c3u))*cv/(2*pow(pow(4*c2u*c3u-8*s2u*(s3u/3+\
    1),2)+pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2),3/2))+(-12*c2u*s3u-16*c2u*(\
    s3u/3+1)-16*s2u*c3u)*cv/sqrt(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+pow(-8\
    *c2u*(s3u/3+1)-4*s2u*c3u,2))+8*c2u*(s3u/3+1)+4*s2u*c3u;

  Du[1] = -6*s3u*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u)*sv/sqrt(36*\
    pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*\
    s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)\
    +4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2))-18*c3u*(4*c2u*c3u-8\
    *s2u*(s3u/3+1))*sv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+\
    36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4\
    *s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),\
    2),2))+3*(4*c2u*c3u-8*s2u*(s3u/3+1))*s3u*(2*((-8*c2u*(s3u/3+1)-4*s2u\
    *c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2))*\
    ((8*c2u*(s3u/3+1)+4*s2u*c3u)*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u\
    )+(-8*c2u*(s3u/3+1)-4*s2u*c3u)*(-12*s2u*s3u-16*s2u*(s3u/3+1)+16*c2u*\
    c3u)-2*(4*c2u*c3u-8*s2u*(s3u/3+1))*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*\
    s2u*c3u))+72*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*s3u*s3u*(12*s2u*s3u+16*s2u\
    *(s3u/3+1)-16*c2u*c3u)+72*(4*c2u*c3u-8*s2u*(s3u/3+1))*s3u*s3u*(-12*\
    c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*c3u)+216*c3u*pow(4*c2u*c3u-8*s2u*(\
    s3u/3+1),2)*s3u+216*c3u*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u)*sv/\
    pow(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3\
    +1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(\
    s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2),3/2)-(-8*\
    c2u*(s3u/3+1)-4*s2u*c3u)*(2*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*(12*s2u*s3u\
    +16*s2u*(s3u/3+1)-16*c2u*c3u)+2*(4*c2u*c3u-8*s2u*(s3u/3+1))*(-12*c2u\
    *s3u-16*c2u*(s3u/3+1)-16*s2u*c3u))*cv/(2*pow(pow(4*c2u*c3u-8*s2u*(\
    s3u/3+1),2)+pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2),3/2))+(12*s2u*s3u+16*\
    s2u*(s3u/3+1)-16*c2u*c3u)*cv/sqrt(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+\
    pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2))-8*s2u*(s3u/3+1)+4*c2u*c3u;

  Du[2] = ((8*c2u*(s3u/3+1)+4*s2u*c3u)*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*\
    c2u*c3u)+(-8*c2u*(s3u/3+1)-4*s2u*c3u)*(-12*s2u*s3u-16*s2u*(s3u/3+1)+\
    16*c2u*c3u)-2*(4*c2u*c3u-8*s2u*(s3u/3+1))*(-12*c2u*s3u-16*c2u*(s3u/3\
    +1)-16*s2u*c3u))*sv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u\
    +36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-\
    4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1)\
    ,2),2))-((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-\
    pow(4*c2u*c3u-8*s2u*(s3u/3+1),2))*(2*((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(\
    8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2))*((8*c2u\
    *(s3u/3+1)+4*s2u*c3u)*(12*s2u*s3u+16*s2u*(s3u/3+1)-16*c2u*c3u)+(-8*\
    c2u*(s3u/3+1)-4*s2u*c3u)*(-12*s2u*s3u-16*s2u*(s3u/3+1)+16*c2u*c3u)-2\
    *(4*c2u*c3u-8*s2u*(s3u/3+1))*(-12*c2u*s3u-16*c2u*(s3u/3+1)-16*s2u*\
    c3u))+72*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*s3u*s3u*(12*s2u*s3u+16*s2u*(\
    s3u/3+1)-16*c2u*c3u)+72*(4*c2u*c3u-8*s2u*(s3u/3+1))*s3u*s3u*(-12*c2u\
    *s3u-16*c2u*(s3u/3+1)-16*s2u*c3u)+216*c3u*pow(4*c2u*c3u-8*s2u*(s3u/3\
    +1),2)*s3u+216*c3u*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u)*sv/(2*pow(\
    36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-\
    4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+\
    1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2),3/2))-6*s3u;

  Dv[0] = 6*(-8*c2u*(s3u/3+1)-4*s2u*c3u)*s3u*cv/sqrt(36*pow(4*c2u*c3u-8*\
    s2u*(s3u/3+1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*\
    s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow\
    (4*c2u*c3u-8*s2u*(s3u/3+1),2),2))-(4*c2u*c3u-8*s2u*(s3u/3+1))*sv/\
    sqrt(pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)+pow(-8*c2u*(s3u/3+1)-4*s2u*c3u\
    ,2));

  Dv[1] = -(-8*c2u*(s3u/3+1)-4*s2u*c3u)*sv/sqrt(pow(4*c2u*c3u-8*s2u*(s3u\
    /3+1),2)+pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2))-6*(4*c2u*c3u-8*s2u*(s3u/\
    3+1))*s3u*cv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+1),2)*s3u*s3u+36*pow\
    (-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*c2u*(s3u/3+1)-4*s2u*\
    c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8*s2u*(s3u/3+1),2),2)\
    );

  Dv[2] = ((-8*c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(\
    4*c2u*c3u-8*s2u*(s3u/3+1),2))*cv/sqrt(36*pow(4*c2u*c3u-8*s2u*(s3u/3+\
    1),2)*s3u*s3u+36*pow(-8*c2u*(s3u/3+1)-4*s2u*c3u,2)*s3u*s3u+pow((-8*\
    c2u*(s3u/3+1)-4*s2u*c3u)*(8*c2u*(s3u/3+1)+4*s2u*c3u)-pow(4*c2u*c3u-8\
    *s2u*(s3u/3+1),2),2));
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
