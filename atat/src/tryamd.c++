#include "version.h"
#include <fstream.h>
#include "linalg.h"

Real myfr(Real r) {
  return -r*exp(-r);
}

Real myf(Real x, Real y) {
  rVector2d z(x,y);
  rVector2d z0;
  Real acc=0.;
  for (z0(0)=-0.5; z0(0)<=2.6; z0(0)+=10.) {
    for (z0(1)=-0.5; z0(1)<=2.6; z0(1)+=10.) {
      acc+=myfr(norm(z-z0));
    }
  }
  return  acc;
}

int main(int argc, char *argv[]) {
  Real dx=0.01;
  ofstream file("ftmp.out");
  Array<Real> x(2);
  Array<Real> f(2);
  Array2d<Real> H(2,2);
  for (x(0)=-2.; x(0)<=2.; x(0)+=0.05) {
    for (x(1)=-2.; x(1)<=2.; x(1)+=0.05) {
      f(0)=(myf(x(0)+dx,x(1))-myf(x(0)-dx,x(1)))/(2.*dx);
      f(1)=(myf(x(0),x(1)+dx)-myf(x(0),x(1)-dx))/(2.*dx);
      H(0,0)=(myf(x(0)+dx,x(1))+myf(x(0)-dx,x(1))-2.*myf(x(0),x(1)))/sqr(dx);
      H(1,1)=(myf(x(0),x(1)+dx)+myf(x(0),x(1)-dx)-2.*myf(x(0),x(1)))/sqr(dx);
      H(1,0)=(myf(x(0)+dx,x(1)+dx)-myf(x(0)+dx,x(1)-dx) - ( myf(x(0)-dx,x(1)+dx)-myf(x(0)-dx,x(1)-dx) ) )/sqr(dx)/4.0;
      H(0,1)=H(1,0);
      Real t=(H(0,0)+H(1,1))/2.;
      Real d=H(0,0)*H(1,1)-H(1,0)*H(0,1); // 1*x^2-2*t*x+d;
      Real mincurv=t-sqrt(sqr(t)-d);
      Real maxcurv=t+sqrt(sqr(t)-d);
      file << x(0) << " " << x(1) << " " << myf(x(0),x(1)) << " " << mincurv << " " << maxcurv << endl;
    }
    file << endl;
  }
}
