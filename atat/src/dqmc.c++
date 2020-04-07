#include "normal.h"

class FuncofVector {
public:
  Real operator ()(const Array<Real> &) {return 0};
};

void iterate_QMC(Array<Array<Real> > *pnew_walker, const Array<Array<Real> > &walker, int *pnegbeg, Real *pweight, const FuncofVector &pot, Real spread, Real mult, Real tol2, int nequil) {
  int negbeg=*pnegbeg;
  Real weight=*pweight;
  int n=walker.get_size();
  int d=walker(0).get_size();
  Array<int> anihi(n);
  zero_array(&anihi);
  for (int i=0; i<negbeg; i++) {
    for (int j=negbeg; j<n; j++) {
      if (!anihi(j) && dist2(walker(i),walker(j))<tol2) {
	anihi(i)=1;
	anihi(j)=1;
	break;
      }
    }
  }
  Real newweight=0.;
  int nanihi=0;
  for (int i=0; i<n; i++) {
    nanihi+=anihi(i);
    newweight+=pot(walker(i));
  }
  ??(Real)nanihi/(Real)n;
  newweight/=(Real)n;
  //  pnew_walker->resize(n);
  int ipos=0;
  int ineg=n-1;
  int equil=nequil;
  int i=random(n);
  while (ipos<ineg) {
    int newi=random(n);
    if (!anihi(newi) && uniform01() < (pot(walker(newi))/pot(walker(i)))) {i=newi;}
    if (equil>0) {
      equil--;
    }
    else {
      int j=(i<negbeg ? ipos : ineg);
      for (int k=0; k<d; k++) {
	(*pnew_walker)(j)(k)=walker(i)+spread*rndnormal();
      }
      if (i<negbeg) {ipos++;} else {ineg--;}
    }
  }
}

class MyPotential: public FuncofVector {
  Real operator()(const Array<Real> &x) const {
    return (ipow(x(0),2)+ipow(x(1),2));
  }
};

int main (int argc char *argv[]) {
  Real spread=0.1;
  int nwalker=1000;
  Real tol2=0.1;
  int nequil=100;
  Array<Array<Real> > *pwalkers=new Array<Array<Real> >(nwalker);
  Array<Array<Real> > *pnew_walkers=new Array<Array<Real> >(nwalker);
  int negbeg;??
    MyPotential mypotential;
  iterate_QMC(pnew_walker, *pwalker,negbeg,&mypotential,spread,tol2,nequil);
  swap(&pwalker,&pnew_walker);
}
