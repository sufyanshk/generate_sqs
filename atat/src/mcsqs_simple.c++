#include <fstream.h>
#include "clus_str.h"
#include "parse.h"
#include "getvalue.h"
#include "version.h"
#include "calccorr.h"

extern char *helpstring;


void calc_all_corr(Array<Real> *pcorr, const Structure &str, const Array<Array<MultiCluster> > &eqclus, const rMatrix3d unitcell, const Array<Array<Array<Real> > > &corrfunc) {
  pcorr->resize(eqclus.get_size());
  for (int t=0; t<eqclus.get_size(); t++) {
    (*pcorr)(t)=calc_correlation(str,eqclus(t),unitcell,corrfunc);
  }
}

Real calc_objective_func(const Array<Real> &corr, const Array<Real> &tcorr, Real mysqstol, Real weightdist, const Array<Real> &diam) {
  Array<Real> dcorr(tcorr.get_size());
  for (int t=0; t<tcorr.get_size(); t++) {
    dcorr(t)=fabs(corr(t)-tcorr(t));
  }
  Real objdist=MAXFLOAT;
  Real objdev=0.;
  for (int t=0; t<tcorr.get_size(); t++) {
    if (dcorr(t)>mysqstol) {
      objdist=min(diam(t),objdist);
      objdev+=dcorr(t);
    }
  }
  return -objdist*weightdist+objdev;
}

/*
int match_point(const Array<Real> &corr, const Array<Real> &tcorr, Real mysqstol, const Array<Real> &diam) {
  for (int t=0; t<tcorr.get_size(); t++) {
    if (diam(t)==0) {
      if (fabs(corr(t)-tcorr(t))>mysqstol) return 0;
    }
  }
  return 1
}
*/

void open_numbered_file(ofstream &file, char *prefix, int num, char *suffix) {
  ostrstream filename;
  if (num>-1) {
    filename << prefix << num << suffix << '\0';
  }
  else {
    filename << prefix << suffix << '\0';
  }
  file.open(filename.str());
  if (!file) {
    cerr << "Unable to open " << filename.str() << endl;
    ERRORQUIT("Aborting");
  }
}

class SupercellData {
public:
  SupercellData(void): str(), nbcomp(), corr(), obj() {};
  Structure str;
  Array<int> nbcomp;
  Array<Real> corr;
  Real obj;
public:
  void operator=(const SupercellData &mc) {
    str=mc.str;
    nbcomp=mc.nbcomp;
    corr=mc.corr;
    obj=mc.obj;
  }
};

int main(int argc, char *argv[]) {
  int maxvol=0; 
  char *latticefilename="lat.in";
  char *clusterfilename="clusters.out";
  char *corrfilename="tcorr.out";
  int readcell=0;
  int ip=-1;
  Real mysqstol=zero_tolerance;
  Real weightdist=1.;
  int seed=0;
  int maxtic=1000;
  int sigdig=6;
  char *corrfunc_label="trigo";
  int dohelp=0;
  AskStruct options[]={
    {"","GENerate Special Quasirandom Structures " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-n","nb of atom/unit cell",INTVAL,&maxvol},
    {"-sig","Number of significant digits to print in output files",INTVAL,&sigdig},
    {"-cf","Input file defining the clusters (default: clusters.out)",STRINGVAL,&clusterfilename},
    {"-tc","Input file defining the target correlations (default: tcorr.out)",STRINGVAL,&corrfilename},
    {"-tol","Tolerance for matching correlations (default: 1e-3)",REALVAL,&mysqstol},
    {"-wr","Weight assigned to range of perfect correlation match is objective function",REALVAL,&weightdist},
    {"-l","Input file defining the lattice (default: lat.in)",STRINGVAL,&latticefilename},
    {"-rc","Read unit cells from file",BOOLVAL,&readcell},
    {"-ip","Index of current process (for parallel operation)",INTVAL,&ip},
    {"-crf","Select correlation functions (default: trigo)",STRINGVAL,&corrfunc_label},
    {"-sd","Seed for random number generation (default: use clock)",INTVAL,&seed},
    {"-h","Display more help",BOOLVAL,&dohelp}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (dohelp) {
    cout << helpstring;
    return 1;
  }

  // read in lattice (see parse.h);
  Structure lat;
  Array<Arrayint> site_type_list;
  Array<AutoString> atom_label;
  rMatrix3d axes;
  ifstream file(latticefilename);
  if (!file) ERRORQUIT("Unable to open lattice file.");
  parse_lattice_file(&lat.cell, &lat.atom_pos, &lat.atom_type, &site_type_list, &atom_label, file, &axes);

  Array<int> comp(lat.atom_type.get_size());
  for (int at=0; at<comp.get_size(); at++) {
    comp(at)=site_type_list(lat.atom_type(at)).get_size();
  }

  // initialize a table of correlation functions;
  if (!check_plug_in(CorrFuncTable(),corrfunc_label)) {
    ERRORQUIT("Aborting");
  }
  CorrFuncTable *pcorrfunc=GenericPlugIn<CorrFuncTable>::create(corrfunc_label);
  pcorrfunc->init(max(comp));

  LinkedList<MultiCluster> clusterlist;
  LinkedList<Real> corrlist;
  {
    ifstream clusterfile(clusterfilename);
    if (!clusterfile) ERRORQUIT("Unable to open cluster file.");
    ifstream corrfile(corrfilename);
    if (!corrfile) ERRORQUIT("Unable to open target correlation file.");
    read_clusters_and_eci(&clusterlist, &corrlist, clusterfile, corrfile, axes);
  }
  Array<Real> tcorr;
  LinkedList_to_Array(&tcorr,corrlist);

  cout.setf(ios::fixed);
  cout.precision(sigdig);

  Array<rMatrix3d> pointgroup;
  find_pointgroup(&pointgroup,lat.cell);
  SpaceGroup spacegroup;
  find_spacegroup(&spacegroup.point_op,&spacegroup.trans,lat.cell,lat.atom_pos,lat.atom_type);
  spacegroup.cell=lat.cell;

  Array<Array<MultiCluster> > eqclus(clusterlist.get_size());
  Array<Real> diam(clusterlist.get_size());
  LinkedListIterator<MultiCluster> ic(clusterlist);
  for (int t=0; ic; t++, ic++) {
    find_equivalent_clusters(&eqclus(t),*ic,lat.cell,spacegroup.point_op,spacegroup.trans);
    diam(t)=get_length_quick(ic->clus);
  }

  ofstream logfile;
  open_numbered_file(logfile, "mcsqs",ip,".log");

  logfile << "Generating supercells..." << endl;
  Array<rMatrix3d> supercell;
  int v=(int)(maxvol/MAX(1,lat.atom_pos.get_size()));
  if (readcell) {
    ifstream cellfile("sqscell.out");
    int nc=0;
    cellfile >> nc;
    supercell.resize(nc);
    for (int i=0; i<supercell.get_size(); i++) {
      read_cell(&(supercell(i)),cellfile);
    }
  }
  else {
    find_supercells(&supercell, v, v, lat.cell, pointgroup);
    ofstream cellfile("sqscell.out");
    cellfile.setf(ios::fixed);
    cellfile.precision(sigdig);
    cellfile << supercell.get_size() << endl << endl;
    for (int i=0; i<supercell.get_size(); i++) {
      supercell(i)=find_symmetric_cell(supercell(i));
      write_axes(supercell(i),cellfile,0);
      cellfile << endl;
    }
  }

  rndseed(seed);

  Array<SupercellData> mc(supercell.get_size());
  SupercellData best;
  best.obj=MAXFLOAT;
  int cc=0;

  logfile << "Initializing random supercells..." << endl;
  for (int c=0; c<supercell.get_size(); c++) {
    mc(c).str.cell=supercell(c);
    find_all_atom_in_supercell(&(mc(c).str.atom_pos),
		 &(mc(c).str.atom_type),lat.atom_pos,
		 comp,
		 lat.cell, mc(c).str.cell);
    mc(c).nbcomp=mc(c).str.atom_type;
    for (int i=0; i<mc(c).str.atom_type.get_size(); i++) {
      if (mc(c).str.atom_type(i)==1) {
	mc(c).str.atom_type(i)=0;
      }
      else {
	mc(c).str.atom_type(i)=random(mc(c).nbcomp(i));
      }
    }
    calc_all_corr(&(mc(c).corr),mc(c).str,eqclus,lat.cell,*pcorrfunc);
    mc(c).obj=calc_objective_func(mc(c).corr,tcorr,mysqstol,weightdist,diam);
    if (mc(c).obj<best.obj) {
      best=mc(c);
      cc=c;
    }
  }
  logfile << "Initialization done." << endl;

  int tic=0;
  Real T=1.;
  Real obj=best.obj;
  best.obj=obj+1.;
  while (1) {
    if (obj<best.obj) {
 cerr << "BBBBBBBBBBBBBBBBBBBBBBBBBBBBB" << endl;
      best=mc(cc);
      ofstream strfile;
      open_numbered_file(strfile, "bestsqs",ip,".out");
      strfile.setf(ios::fixed);
      strfile.precision(sigdig);
      write_structure(best.str,lat,site_type_list,atom_label,axes,strfile);
      ofstream corrfile;
      open_numbered_file(corrfile, "bestcorr",ip,".out");
      corrfile.setf(ios::fixed);
      corrfile.precision(sigdig);
      for (int t=0; t<best.corr.get_size(); t++) {
	corrfile << best.corr(t) << endl;
      }
      logfile << best.obj << endl;
    }

    if (tic==0) {
      ifstream tempfile("temp.in");
      if (!tempfile) {ERRORQUIT("Unable to open temp.in file containing the current temperature.");}
      tempfile >> T;
      if (file_exists("stopsqs")) {break;}
      tic=maxtic;
    }
    int newcc=random(mc.get_size());
    //    int doswap=0;
    //    if (match_point(mc(newcc).corr,tcorr,mysqstol,diam)) {doswap=1;}
    int at;
    do {
      at=random(mc(newcc).str.atom_type.get_size());
    } while (mc(newcc).nbcomp(at)<2);
    int saveat=mc(newcc).str.atom_type(at);
    mc(newcc).str.atom_type(at)=(saveat+random(mc(newcc).nbcomp(at)-1)) % mc(newcc).nbcomp(at);
    
    Array<Real> newcorr;
    calc_all_corr(&newcorr,mc(newcc).str,eqclus,lat.cell,*pcorrfunc);
    Real newobj=calc_objective_func(newcorr,tcorr,mysqstol,weightdist,diam);
  cerr << newcc << " " << best.obj << " " << newobj << " ";
  for (int t=0; t<newcorr.get_size(); t++) {cerr << newcorr(t) << " ";}
    if (uniform01() < exp((obj-newobj)/T) ) {
      mc(newcc).corr=newcorr;
      mc(newcc).obj=newobj;
      cc=newcc;
      obj=newobj;
  cerr << "A";
    }
    else {
      mc(newcc).str.atom_type(at)=saveat;
  cerr << "r"; 
    }
  cerr << endl << flush;
    tic--;
  }
}
