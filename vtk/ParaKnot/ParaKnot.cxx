#include "ParaKnot.h"
#include "vtkObjectFactory.h"
#include "vtkMath.h"

vtkStandardNewMacro(myParametricParaKnot);

//----------------------------------------------------------------------------
myParametricParaKnot::myParametricParaKnot()
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
myParametricParaKnot::~myParametricParaKnot()
{
}

//----------------------------------------------------------------------------
void myParametricParaKnot::Evaluate(double uvw[3], double Pt[3], double Duvw[9])
{
  double u = uvw[0];
  double v = uvw[1];
  double *Du = Duvw;
  double *Dv = Duvw + 3;

  double p = this->SpeedP;
  double q = this->SpeedQ;
  double m = this->MultP;
  double n = this->MultQ;
  double r = this->TubeRadius;
  double w = this->WaveAmplitude;

Pt[0] = -(abs(q)+p)*r*(n*q*sin(q*u)+m*p*sin(p*u))*cos((abs(q)+p)*u)*sin(v)*w
                   *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                                       *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                         +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                                         *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                         +pow((-n*q*sin(q*u)-m*p*sin(p*u))
                               *(n*q*sin(q*u)+m*p*sin(p*u))
                               -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-1.0/2.0)
      +r*(n*q*cos(q*u)+m*p*cos(p*u))
        *pow(pow(n*q*sin(q*u)+m*p*sin(p*u),2)
              +pow(n*q*cos(q*u)+m*p*cos(p*u),2),-1.0/2.0)*cos(v)+n*cos(q*u)
      +m*cos(p*u);

Pt[1] = (abs(q)+p)*r*(n*q*cos(q*u)+m*p*cos(p*u))*cos((abs(q)+p)*u)*sin(v)*w
                  *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                                      *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                        +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                                        *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                        +pow((-n*q*sin(q*u)-m*p*sin(p*u))
                              *(n*q*sin(q*u)+m*p*sin(p*u))
                              -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-1.0/2.0)
      +r*(n*q*sin(q*u)+m*p*sin(p*u))
        *pow(pow(n*q*sin(q*u)+m*p*sin(p*u),2)
              +pow(n*q*cos(q*u)+m*p*cos(p*u),2),-1.0/2.0)*cos(v)+n*sin(q*u)
      +m*sin(p*u);

Pt[2] = r*((-n*q*sin(q*u)-m*p*sin(p*u))*(n*q*sin(q*u)+m*p*sin(p*u))
          -pow(n*q*cos(q*u)+m*p*cos(p*u),2))*sin(v)
         *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                             *pow(cos((abs(q)+p)*u),2)*pow(w,2)
               +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                               *pow(cos((abs(q)+p)*u),2)*pow(w,2)
               +pow((-n*q*sin(q*u)-m*p*sin(p*u))*(n*q*sin(q*u)+m*p*sin(p*u))
                     -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-1.0/2.0)
      +sin((abs(q)+p)*u)*w;

Du[0] = pow(abs(q)+p,2)*r*(n*q*sin(q*u)+m*p*sin(p*u))*sin((abs(q)+p)*u)*sin(v)
                       *w
                       *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                                           *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                             +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                                             *pow(cos((abs(q)+p)*u),2)
                                             *pow(w,2)
                             +pow((-n*q*sin(q*u)-m*p*sin(p*u))
                                   *(n*q*sin(q*u)+m*p*sin(p*u))
                                   -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-1.0/2.0)
      -(abs(q)+p)*r*(n*pow(q,2)*cos(q*u)+m*pow(p,2)*cos(p*u))
                 *cos((abs(q)+p)*u)*sin(v)*w
                 *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                                     *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                       +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                                       *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                       +pow((-n*q*sin(q*u)-m*p*sin(p*u))
                             *(n*q*sin(q*u)+m*p*sin(p*u))
                             -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-1.0/2.0)
      +(abs(q)+p)*r*(n*q*sin(q*u)+m*p*sin(p*u))*cos((abs(q)+p)*u)*sin(v)*w
                 *(2*pow(abs(q)+p,2)*(n*q*cos(q*u)+m*p*cos(p*u))
                    *(-n*pow(q,2)*sin(q*u)-m*pow(p,2)*sin(p*u))
                    *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                  +2*pow(abs(q)+p,2)*(n*pow(q,2)*cos(q*u)+m*pow(p,2)*cos(p*u))
                    *(n*q*sin(q*u)+m*p*sin(p*u))*pow(cos((abs(q)+p)*u),2)
                    *pow(w,2)
                  -2*pow(abs(q)+p,3)*cos((abs(q)+p)*u)*sin((abs(q)+p)*u)
                    *pow(n*q*sin(q*u)+m*p*sin(p*u),2)*pow(w,2)
                  -2*pow(abs(q)+p,3)*cos((abs(q)+p)*u)*sin((abs(q)+p)*u)
                    *pow(n*q*cos(q*u)+m*p*cos(p*u),2)*pow(w,2)
                  +2*(-2*(n*q*cos(q*u)+m*p*cos(p*u))
                        *(-n*pow(q,2)*sin(q*u)-m*pow(p,2)*sin(p*u))
                     +(-n*pow(q,2)*cos(q*u)-m*pow(p,2)*cos(p*u))
                      *(n*q*sin(q*u)+m*p*sin(p*u))
                     +(n*pow(q,2)*cos(q*u)+m*pow(p,2)*cos(p*u))
                      *(-n*q*sin(q*u)-m*p*sin(p*u)))
                    *((-n*q*sin(q*u)-m*p*sin(p*u))*(n*q*sin(q*u)+m*p*sin(p*u))
                     -pow(n*q*cos(q*u)+m*p*cos(p*u),2)))
                 *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                                     *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                       +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                                       *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                       +pow((-n*q*sin(q*u)-m*p*sin(p*u))
                             *(n*q*sin(q*u)+m*p*sin(p*u))
                             -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-3/2)
       /2
      +r*(-n*pow(q,2)*sin(q*u)-m*pow(p,2)*sin(p*u))
        *pow(pow(n*q*sin(q*u)+m*p*sin(p*u),2)
              +pow(n*q*cos(q*u)+m*p*cos(p*u),2),-1.0/2.0)*cos(v)
      -r*(n*q*cos(q*u)+m*p*cos(p*u))
        *(2*(n*q*cos(q*u)+m*p*cos(p*u))
           *(-n*pow(q,2)*sin(q*u)-m*pow(p,2)*sin(p*u))
         +2*(n*pow(q,2)*cos(q*u)+m*pow(p,2)*cos(p*u))
           *(n*q*sin(q*u)+m*p*sin(p*u)))
        *pow(pow(n*q*sin(q*u)+m*p*sin(p*u),2)
              +pow(n*q*cos(q*u)+m*p*cos(p*u),2),-3/2)*cos(v)
       /2-n*q*sin(q*u)-m*p*sin(p*u);

Du[1] = -pow(abs(q)+p,2)*r*(n*q*cos(q*u)+m*p*cos(p*u))*sin((abs(q)+p)*u)
                        *sin(v)*w
                        *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                                            *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                              +pow(abs(q)+p,2)
                               *pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                               *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                              +pow((-n*q*sin(q*u)-m*p*sin(p*u))
                                    *(n*q*sin(q*u)+m*p*sin(p*u))
                                    -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-1.0/2.0)
      +(abs(q)+p)*r*(-n*pow(q,2)*sin(q*u)-m*pow(p,2)*sin(p*u))
                 *cos((abs(q)+p)*u)*sin(v)*w
                 *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                                     *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                       +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                                       *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                       +pow((-n*q*sin(q*u)-m*p*sin(p*u))
                             *(n*q*sin(q*u)+m*p*sin(p*u))
                             -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-1.0/2.0)
      -(abs(q)+p)*r*(n*q*cos(q*u)+m*p*cos(p*u))*cos((abs(q)+p)*u)*sin(v)*w
                 *(2*pow(abs(q)+p,2)*(n*q*cos(q*u)+m*p*cos(p*u))
                    *(-n*pow(q,2)*sin(q*u)-m*pow(p,2)*sin(p*u))
                    *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                  +2*pow(abs(q)+p,2)*(n*pow(q,2)*cos(q*u)+m*pow(p,2)*cos(p*u))
                    *(n*q*sin(q*u)+m*p*sin(p*u))*pow(cos((abs(q)+p)*u),2)
                    *pow(w,2)
                  -2*pow(abs(q)+p,3)*cos((abs(q)+p)*u)*sin((abs(q)+p)*u)
                    *pow(n*q*sin(q*u)+m*p*sin(p*u),2)*pow(w,2)
                  -2*pow(abs(q)+p,3)*cos((abs(q)+p)*u)*sin((abs(q)+p)*u)
                    *pow(n*q*cos(q*u)+m*p*cos(p*u),2)*pow(w,2)
                  +2*(-2*(n*q*cos(q*u)+m*p*cos(p*u))
                        *(-n*pow(q,2)*sin(q*u)-m*pow(p,2)*sin(p*u))
                     +(-n*pow(q,2)*cos(q*u)-m*pow(p,2)*cos(p*u))
                      *(n*q*sin(q*u)+m*p*sin(p*u))
                     +(n*pow(q,2)*cos(q*u)+m*pow(p,2)*cos(p*u))
                      *(-n*q*sin(q*u)-m*p*sin(p*u)))
                    *((-n*q*sin(q*u)-m*p*sin(p*u))*(n*q*sin(q*u)+m*p*sin(p*u))
                     -pow(n*q*cos(q*u)+m*p*cos(p*u),2)))
                 *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                                     *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                       +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                                       *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                       +pow((-n*q*sin(q*u)-m*p*sin(p*u))
                             *(n*q*sin(q*u)+m*p*sin(p*u))
                             -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-3/2)
       /2
      +r*(n*pow(q,2)*cos(q*u)+m*pow(p,2)*cos(p*u))
        *pow(pow(n*q*sin(q*u)+m*p*sin(p*u),2)
              +pow(n*q*cos(q*u)+m*p*cos(p*u),2),-1.0/2.0)*cos(v)
      -r*(n*q*sin(q*u)+m*p*sin(p*u))
        *(2*(n*q*cos(q*u)+m*p*cos(p*u))
           *(-n*pow(q,2)*sin(q*u)-m*pow(p,2)*sin(p*u))
         +2*(n*pow(q,2)*cos(q*u)+m*pow(p,2)*cos(p*u))
           *(n*q*sin(q*u)+m*p*sin(p*u)))
        *pow(pow(n*q*sin(q*u)+m*p*sin(p*u),2)
              +pow(n*q*cos(q*u)+m*p*cos(p*u),2),-3/2)*cos(v)
       /2+n*q*cos(q*u)+m*p*cos(p*u);

Du[2] = r*(-2*(n*q*cos(q*u)+m*p*cos(p*u))
             *(-n*pow(q,2)*sin(q*u)-m*pow(p,2)*sin(p*u))
          +(-n*pow(q,2)*cos(q*u)-m*pow(p,2)*cos(p*u))
           *(n*q*sin(q*u)+m*p*sin(p*u))
          +(n*pow(q,2)*cos(q*u)+m*pow(p,2)*cos(p*u))
           *(-n*q*sin(q*u)-m*p*sin(p*u)))*sin(v)
         *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                             *pow(cos((abs(q)+p)*u),2)*pow(w,2)
               +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                               *pow(cos((abs(q)+p)*u),2)*pow(w,2)
               +pow((-n*q*sin(q*u)-m*p*sin(p*u))*(n*q*sin(q*u)+m*p*sin(p*u))
                     -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-1.0/2.0)
      -r*((-n*q*sin(q*u)-m*p*sin(p*u))*(n*q*sin(q*u)+m*p*sin(p*u))
         -pow(n*q*cos(q*u)+m*p*cos(p*u),2))*sin(v)
        *(2*pow(abs(q)+p,2)*(n*q*cos(q*u)+m*p*cos(p*u))
           *(-n*pow(q,2)*sin(q*u)-m*pow(p,2)*sin(p*u))
           *pow(cos((abs(q)+p)*u),2)*pow(w,2)
         +2*pow(abs(q)+p,2)*(n*pow(q,2)*cos(q*u)+m*pow(p,2)*cos(p*u))
           *(n*q*sin(q*u)+m*p*sin(p*u))*pow(cos((abs(q)+p)*u),2)*pow(w,2)
         -2*pow(abs(q)+p,3)*cos((abs(q)+p)*u)*sin((abs(q)+p)*u)
           *pow(n*q*sin(q*u)+m*p*sin(p*u),2)*pow(w,2)
         -2*pow(abs(q)+p,3)*cos((abs(q)+p)*u)*sin((abs(q)+p)*u)
           *pow(n*q*cos(q*u)+m*p*cos(p*u),2)*pow(w,2)
         +2*(-2*(n*q*cos(q*u)+m*p*cos(p*u))
               *(-n*pow(q,2)*sin(q*u)-m*pow(p,2)*sin(p*u))
            +(-n*pow(q,2)*cos(q*u)-m*pow(p,2)*cos(p*u))
             *(n*q*sin(q*u)+m*p*sin(p*u))
            +(n*pow(q,2)*cos(q*u)+m*pow(p,2)*cos(p*u))
             *(-n*q*sin(q*u)-m*p*sin(p*u)))
           *((-n*q*sin(q*u)-m*p*sin(p*u))*(n*q*sin(q*u)+m*p*sin(p*u))
            -pow(n*q*cos(q*u)+m*p*cos(p*u),2)))
        *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                            *pow(cos((abs(q)+p)*u),2)*pow(w,2)
              +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                              *pow(cos((abs(q)+p)*u),2)*pow(w,2)
              +pow((-n*q*sin(q*u)-m*p*sin(p*u))*(n*q*sin(q*u)+m*p*sin(p*u))
                    -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-3/2)
       /2+(abs(q)+p)*cos((abs(q)+p)*u)*w;

Dv[0] = -(abs(q)+p)*r*(n*q*sin(q*u)+m*p*sin(p*u))*cos((abs(q)+p)*u)*cos(v)*w
                   *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                                       *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                         +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                                         *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                         +pow((-n*q*sin(q*u)-m*p*sin(p*u))
                               *(n*q*sin(q*u)+m*p*sin(p*u))
                               -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-1.0/2.0)
      -r*(n*q*cos(q*u)+m*p*cos(p*u))
        *pow(pow(n*q*sin(q*u)+m*p*sin(p*u),2)
              +pow(n*q*cos(q*u)+m*p*cos(p*u),2),-1.0/2.0)*sin(v);

Dv[1] = (abs(q)+p)*r*(n*q*cos(q*u)+m*p*cos(p*u))*cos((abs(q)+p)*u)*cos(v)*w
                  *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                                      *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                        +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                                        *pow(cos((abs(q)+p)*u),2)*pow(w,2)
                        +pow((-n*q*sin(q*u)-m*p*sin(p*u))
                              *(n*q*sin(q*u)+m*p*sin(p*u))
                              -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-1.0/2.0)
      -r*(n*q*sin(q*u)+m*p*sin(p*u))
        *pow(pow(n*q*sin(q*u)+m*p*sin(p*u),2)
              +pow(n*q*cos(q*u)+m*p*cos(p*u),2),-1.0/2.0)*sin(v);

Dv[2] = r*((-n*q*sin(q*u)-m*p*sin(p*u))*(n*q*sin(q*u)+m*p*sin(p*u))
          -pow(n*q*cos(q*u)+m*p*cos(p*u),2))*cos(v)
         *pow(pow(abs(q)+p,2)*pow(n*q*sin(q*u)+m*p*sin(p*u),2)
                             *pow(cos((abs(q)+p)*u),2)*pow(w,2)
               +pow(abs(q)+p,2)*pow(n*q*cos(q*u)+m*p*cos(p*u),2)
                               *pow(cos((abs(q)+p)*u),2)*pow(w,2)
               +pow((-n*q*sin(q*u)-m*p*sin(p*u))*(n*q*sin(q*u)+m*p*sin(p*u))
                     -pow(n*q*cos(q*u)+m*p*cos(p*u),2),2),-1.0/2.0);
}

//----------------------------------------------------------------------------
double myParametricParaKnot::EvaluateScalar(double* vtkNotUsed(uv[3]),
                                          double* vtkNotUsed(Pt[3]),
                                          double* vtkNotUsed(Duv[9]))
{
  return 0;
}

//----------------------------------------------------------------------------
void myParametricParaKnot::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
