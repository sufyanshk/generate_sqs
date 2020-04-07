#include <fstream.h>
#include <sys/stat.h>
#include "parse.h"
#include "getvalue.h"
#include "version.h"

char *helpstring="Insert more help here";

int find_inflect_pt(Real *pinflecty, Real *ppos, const Array<Real> &y) {
  Real tiny=1e-6;
  Real c1=y(0)+y(2)-2.*y(1);
  Real c2=y(1)+y(3)-2.*y(2);
  if (c1*c2>0) {return 0;}
  Real x=0.5;
  if ((c2-c1)>tiny*max(fabs(c1),fabs(c2))) {
    x=-c1/(c2-c1);
  }
  *pinflecty=y(1)*(1-x)+y(2)*x;
  if (ppos) {
    *ppos=x+1.;
  }
  return 1;
}

int find_local_min(Real *pminy, Real *ppos, const Array<Real> &y) {
  if (y(1)<y(0) && y(1)<y(2)) {
    Real b=(y(2)-y(0))/2.;
    Real a=(y(0)+y(2)-2.*y(1))/2.;
    Real c=y(1);
    Real x=-b/(2.*a);
    Real y=a*sqr(x)+b*x+c;
    (*pminy)=y;
    if (ppos) {
      (*ppos)=(x+1.);
    }
    return 1;
  }
  else {
    return 0;
  }
}

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  char *str1filename="str1.out";
  char *str2filename="str2.out";
  int nimage=0;
  Real midf=0.5;
  Real dimerf=0.0;
  int sigdig=5;
  int dohelp=0;
  int dummy=0;
  int doinflect=0;
  int doinflectpos=0;
  char *dirstem="00";
  AskStruct options[]={
    {"","STRucture PATH generator " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-s1","Input file defining structure 1 (default: str1.out)",STRINGVAL,&str1filename},
    {"-s2","Input file defining structure 2 (default: str2.out)",STRINGVAL,&str2filename},
    {"-ni","Number of images",INTVAL,&nimage},
    {"-f","Interpolation fraction in [0,1]",REALVAL,&midf},
    {"-df","Dimer distance [0,1]",REALVAL,&dimerf},
    {"-dp","Directory name pattern (default 00)",STRINGVAL,&dirstem},
    {"-ci","Calculate inflection point",BOOLVAL,&doinflect},
    {"-cil","Calculate inflection point location",BOOLVAL,&doinflectpos},
    {"-sig","Number of significant digits printed (default: 5)",INTVAL,&sigdig},
    {"-h","Display more help",BOOLVAL,&dohelp},
    {"-d","Use all default values",BOOLVAL,&dummy}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }

  cout.setf(ios::fixed);
  cout.precision(sigdig);

  if (doinflect || doinflectpos) {
    LinkedList<Real> energylist;
    while (1) {
      Real tmp=MAXFLOAT;
      cin >> tmp;
      if (tmp==MAXFLOAT) break;
      energylist << new Real(tmp);
    }
    Array<Real> energy;
    LinkedList_to_Array(&energy,energylist);
    int i;
    for (i=0; i<energy.get_size()-2; i++) {
      Real ei=0.;
      Real xi=0.;
      Array<Real> e3;
      extract_elements(&e3,energy,i,i+3);
      int found=0;
      found=find_local_min(&ei,&xi,e3);
      if (found) {cerr << "Local minimum found" << endl;}
      if (!found && i<energy.get_size()-3) {
	Array<Real> e4;
	extract_elements(&e4,energy,i,i+4);
	found=find_inflect_pt(&ei,&xi,e4);
      }
      if (found) {
	if (doinflectpos) {
	  cout << ((Real)i+xi)/(Real)(energy.get_size()-1) << endl;
	}
	else {
	  cout << ei << endl;
	}
	break;
      }
    }
    if (i==energy.get_size()-2) {
      if (doinflectpos) {
	cout << 1. << endl;
      }
      else {
	cout << energy(energy.get_size()-1) << endl;
      }
    }
    exit(0);
  }

  // parsing lattice and structure files. See parse.hh for detail;
  Structure str1;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    ifstream str1file(str1filename);
    if (!str1file) ERRORQUIT("Unable to open structure 1 file");
    parse_lattice_file(&str1.cell, &str1.atom_pos, &str1.atom_type, &labellookup, &label, str1file, &axes);
    wrap_inside_cell(&str1.atom_pos,str1.atom_pos,str1.cell);
    fix_atom_type(&str1, labellookup);
  }
  Structure str2;
  {
    ifstream str2file(str2filename);
    if (!str2file) ERRORQUIT("Unable to open structure 2 file");
    parse_structure_file(&str2.cell, &str2.atom_pos, &str2.atom_type, label, str2file, NULL);
    wrap_inside_cell(&str2.atom_pos,str2.atom_pos,str2.cell);
    reorder_atoms(&str2,str1);
  }

  Real f0,f1,df;
  if (nimage==0) {
    f0=midf;
    df=1;
    f1=midf+df/2;
  }
  else {
    f0=0;
    df=1./(Real)(nimage-1);
    f1=1+df/2;
  }
  int im=0;
  for (Real f=f0; f<f1; f+=df, im++) {
    Structure str;
    str.cell=(1-f)*str1.cell+f*str2.cell;
    str.atom_type=str1.atom_type;
    str.atom_pos.resize(str1.atom_pos.get_size());
    Array<rVector3d> dir;
    dir.resize(str1.atom_pos.get_size());
    Array<rVector3d> frac1;
    Array<rVector3d> frac2;
    rMatrix3d icell1=!str1.cell;
    rMatrix3d icell2=!str2.cell;
    for (int at=0; at<str1.atom_pos.get_size(); at++) {
      rVector3d frac1=icell1*str1.atom_pos(at);
      rVector3d frac2=icell2*str2.atom_pos(at);
      rVector3d cellsh=round(frac2-frac1);
      frac2=frac2-cellsh;
      str.atom_pos(at)=str.cell*((1.-f)*frac1+f*frac2);
      dir(at)=str.cell*(frac2-frac1)*dimerf;
    }
    
    if (nimage==0) {
      write_structure(str,label,axes,cout,0);
      if (dimerf>0) {
	for (int at=0; at<dir.get_size(); at++) {
	  cout << (!axes)*dir(at) << " ZZdimer" << endl;
	}
      }
    }
    else {
      AutoString dirname(dirstem);
      ostrstream num;
      num << im << '\0';
      const char *pnum=num.str();
      char *pdirname=(char *)dirname;
      strcpy(pdirname+strlen(pdirname)-strlen(pnum),pnum);
      mkdir(dirname,S_IRWXU | S_IRWXG | S_IRWXO);
      dirname+=AutoString("/str.out");
      ofstream file(dirname);
      file.setf(ios::fixed);
      file.precision(sigdig);
      write_structure(str,label,axes,file,0);
    }
  }
}
