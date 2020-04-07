#include <fstream.h>
#include "parse.h"
#include "getvalue.h"
#include "version.h"

char *helpstring="Insert more help here";

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  char *latfilename="lat.in";
  char *strfilename="str.out";
  Real r=0;
  int n=0;
  int b=0;
  int sigdig=5;
  int dohelp=0;
  int dummy=0;
  AskStruct options[]={
    {"","Skeleton for atat-like codes " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-l","Input file defining the lattice (defaults: lat.in)",STRINGVAL,&latfilename},
    {"-s","Input file defining the structure (default: str.out)",STRINGVAL,&strfilename},
    {"-r","A real parameter",REALVAL,&r},
    {"-n","An integer parameter",INTVAL,&n},
    {"-b","A boolean parameter",BOOLVAL,&b},
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

  // parsing lattice and structure files. See parse.hh for detail;
  Structure lat;
  Array<Arrayint> labellookup;
  Array<AutoString> label;
  rMatrix3d axes;
  {
    ifstream latfile(latfilename);
    if (!latfile) ERRORQUIT("Unable to open lattice file");
    parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &labellookup, &label, latfile, &axes);
    wrap_inside_cell(&lat.atom_pos,lat.atom_pos,lat.cell);
  }

  rMatrix3d supcel;
  find_smallest_supercell_enclosing_sphere(&supcel,lat.cell,r);
  rMatrix3d supcela=(!axes)*supcel;
  write_axes(axes,cout,0);
  for (int i=0; i<3; i++) {
    cout << supcela.get_column(i) << endl;
  }

}
