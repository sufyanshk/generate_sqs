#include "calccorr.h"

// This is a skeleton to show how to include user-defined cluster functions;
// look for XXXX and replace it by the (internal) name of this cluster function type;
// look for xxxx abd replace it by the user-visible name of this cluster function type (it will be invoked by the -crf=xxxx option);
// look for myexpression and replace it by your own definition of the cluster functions;
// save this file under a new name, say corrXXXX.c++ ;
// add corrXXXX.o to the makefile on the line PLUGINCORR ;
// type "make";

class HALFCorrFuncTable: public CorrFuncTable {
 public:
  HALFCorrFuncTable(): CorrFuncTable() {}
  void init(int comp) {
    Array<Array<Array<Real> > > &table=*this;
    table.resize(comp-1);
    for (int m=2; m<=comp; m++) {
      table(m-2).resize(m-1);
      for (int t=0; t<m-1; t++) {
	table(m-2)(t).resize(m);
      }
      for (int s=0; s<m; s++) {
        for (int t=1; t<=(m/2); t++) {
	  table(m-2)(2*t-2)(s)=double(s)/double((m-1));//-cos(2*M_PI*s*t/m);
        }
        for (int t=1; t<=((m+1)/2-1); t++) {
	  table(m-2)(2*t-1)(s)=sin((s%2)*M_PI/2)+(ipow((Real)(s-1),s))*s/(2*(m-1));//-sin(2*M_PI*s*t/m);
        }
      }
    }
  }
};

SpecificPlugIn<CorrFuncTable,HALFCorrFuncTable> HALFCorrFuncPlugIn("half");
