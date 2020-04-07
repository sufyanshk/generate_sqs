#include <fstream.h>
#include "getvalue.h"
#include "version.h"
#include "integer.h"
#include "linalg.h"
#include "meshutil.h"

template<class T>
ostream& operator << (ostream &s, const LinkedList<T> &a) {
  s << a.get_size() << endl;
  LinkedListIterator<T> i(a);
  for (; i; i++) {
    s << *i << endl;
  }
  return s;
}

void debug_mesh(const LinkedList<Array<int> > & polylist, const Array<Array<Real> > &pts) {
  ofstream file("mesh.out");
  LinkedListIterator<Array<int> > it(polylist);
  for (;it; it++) {
    ofstream last("lastpoly.out");
    for (int i=0; i<it->get_size(); i++) {
      for (int ii=0; ii<i; ii++) {
	for (int j=0; j<pts((*it)(0)).get_size(); j++) {
	  file << pts((*it)(i))(j) << " ";
	  last << pts((*it)(i))(j) << " ";
	}
	file << endl;
	last << endl;
	for (int j=0; j<pts((*it)(0)).get_size(); j++) {
	  file << pts((*it)(ii))(j) << " ";
	  last << pts((*it)(ii))(j) << " ";
	}
	file << endl;
	file << endl;
	file << endl;
	last << endl;
	last << endl;
	last << endl;
      }
    }
  }
}

void find_perp_to_vector_array(Array<Real> *pperp, const Array<Array<Real> > & v) {
  int dim=v(0).get_size();
  pperp->resize(dim);
  Array2d<Real> B(dim,dim),iB;
  for (int i=0; i<dim-1; i++) {
    for (int j=0; j<dim; j++) {
      B(i,j)=v(i)(j);
    }
  }
  for (int j=0; j<dim; j++) {
    B(dim-1,j)=uniform01()-0.5;
  }
  invert_matrix(&iB,B);
  for (int i=0; i<dim; i++) {
    (*pperp)(i)=iB(i,dim-1);
  }
  normalize(pperp);
}

int calc_circumsphere(Array<Real> *pcenter, Real *pradius2, const Array<Array<Real> > &pts) {
  int dim=pts(0).get_size();
  int nbpts=pts.get_size();
  Array2d<Real> A(dim,dim);
  Array<Real> v(dim);
  for (int i=0; i<nbpts-1; i++) {
    Array<Real> delta;
    diff(&delta,pts(i+1),pts(i));
    set_row(&A,delta,i);
    v(i)=(norm2(pts(i+1))-norm2(pts(i)))/2.;
  }
  if (nbpts<dim+1) {
    Array2d<Real> B,iB;
    Array<Real> u(dim);
    B=A;
    for (int i=0; i<dim; i++) {
      B(dim-1,i)=uniform01()-0.5;
    }
    invert_matrix(&iB,B);
    for (int i=0; i<dim; i++) {
      u(i)=iB(i,dim-1);
      A(dim-1,i)=u(i);
    }
    v(dim-1)=inner_product(u,pts(0));
  }
  else if (nbpts != dim+1) {
    ERRORQUIT("wrong dimensions in calc_circumsphere");
  }
  solve_linsys(&A,&v);
  *pcenter=v;
  Array<Real> r;
  diff(&r,pts(0),*pcenter);
  *pradius2=norm2(r);
  return 1;
}

class SideFlag {
public:
  int done;
  Array<int> ipts;
  Array<Real> forward;
  SideFlag(void): done(0),ipts(0),forward(0) {}
  SideFlag(int _done, const Array<int>& _ipts): done(_done),ipts(_ipts),forward(0) {}
};

void debug_side(const LinkedList<SideFlag > & edgelist, const Array<Array<Real> > &pts) {
  ofstream file("normal.out");
  ofstream filedone("done.out");
  LinkedListIterator<SideFlag> it(edgelist);
  for (; it; it++) {
    Array<Real> center;
    center=pts(it->ipts(0));
    for (int i=1; i<it->ipts.get_size(); i++) {
      sum(&center,center,pts(it->ipts(i)));
    }
    product(&center,center,1./(Real)(it->ipts.get_size()));
    for (int j=0; j<center.get_size(); j++) {
      filedone << center(j) << " ";
    }
    filedone << it->done << endl;
    if (it->forward.get_size()!=0) {
      Array<Real> endarrow;
      sum(&endarrow,center,it->forward);
      for (int j=0; j<center.get_size(); j++) {
	file << center(j) << " ";
      }
      file << endl;
      for (int j=0; j<center.get_size(); j++) {
	file << endarrow(j) << " ";
      }
      file << endl;
      file << endl;
      file << endl;
    }
  }
}

class RealIntPair {
public:
  Real r;
  int i;
  RealIntPair(void): r(0.), i(0) {}
  RealIntPair(Real _r, int _i): r(_r), i(_i) {}
};

template<class T>
class OrderBy_r {
 public:
  int operator () (const T& a, const T& b) const {
    return (a.r<b.r);
  }
};

class PointFromCenterIterator {
  const Array<Array<Real> > &pts;
  Array<Real> center;
  Array<RealIntPair> distindex;
  int curindex;
public:
  PointFromCenterIterator(const Array<Array<Real> > &_pts, const Array<Real> &_center): pts(_pts), center(_center), distindex(_pts.get_size()), curindex(0) {
    for (int i=0; i<pts.get_size(); i++) {
      Array<Real> d;
      diff(&d,pts(i),center);
      distindex(i).r=norm(d);
      distindex(i).i=i;
    }
    sort_array(&distindex,OrderBy_r<RealIntPair>());
  }
  void init(void) {curindex=0;}
  operator void * () {return (void *)(curindex<distindex.get_size() ? 1 : 0);}
  int operator () (void)  {
    return distindex(curindex).i;
  }
  void operator++(int) {
    curindex++;
  }
};

int find_point_to_add(const SideFlag &side, const Array<Array<Real> > &pts, Real maxr=MAXFLOAT) {
  Array<Real> midside(pts(side.ipts(0)));
  for (int i=1; i<side.ipts.get_size(); i++) {
    sum(&midside,midside,pts(side.ipts(i)));
  }
  product(&midside,midside,1./(Real)(side.ipts.get_size()));
  Array<Array<Real> > polytope(side.ipts.get_size()+1);
  for (int i=0; i<side.ipts.get_size(); i++) {
    polytope(i+1)=pts(side.ipts(i));
    // cerr << "side[" << i << "]=" << side.ipts(i) << endl;
  }
  PointFromCenterIterator itpts(pts,midside);
  for (; itpts; itpts++) {
    // cerr << "here\n";
    // cerr << "itpts=" << itpts() << endl;
    Array<Real> delta;
    diff(&delta,pts(itpts()),midside);
    if (norm(delta)>maxr) {return -1;}
    if ( inner_product(side.forward,delta)>0 && !is_in_array(side.ipts,itpts()) ) {
      polytope(0)=pts(itpts());
      // cerr << polytope << endl;
      Array<Real> center;
      Real radius2;
      calc_circumsphere(&center,&radius2, polytope);
      // cerr << "center=" << center << endl;
      PointFromCenterIterator jtpts(pts,center);
      int success=0;
      for (; jtpts; jtpts++) {
	if (is_in_array(side.ipts,jtpts()) || jtpts()==itpts()) {
	  success=1;
	  break;
	}
	Array<Real> delta;
	diff(&delta,pts(jtpts()),center);
	// cerr << norm(delta) << "<>" << sqrt(radius2) << endl;
	if (norm2(delta)<radius2) break;
      }
      if (success) {
	// cerr << "picked" << endl;
	break;
      }
    }
  }
  if (!itpts) {return -1;}
  return itpts();
}

void calc_forward_vector(Array<Real> *pforward, const Array<int> &sideipts, int ptipts, const Array<Array<Real> > &pts) {
  int dim=pts(0).get_size();
  Array<Array<Real> > edge(dim-1);
  for (int j=0; j<sideipts.get_size()-1; j++) {
    diff(&(edge(j)),pts(sideipts(j+1)),pts(sideipts(0)));
  }
  // cerr << "edge=" << edge << endl;
  if (sideipts.get_size()<dim) {
    Array<Array<Real> > planeedge(dim-1);
    for (int j=0; j<dim-1; j++) {
      // cerr << "j=" << j << endl;
      // cerr << sideipts(j) << endl;
      diff(&(planeedge(j)),pts(sideipts(j)),pts(ptipts));
    }
    find_perp_to_vector_array(&(edge(dim-2)),planeedge);
  }
  // cerr << "here2\n";
  find_perp_to_vector_array(pforward,edge);
  // cerr << "here3\n";
  Array<Real> testv(dim);
  diff(&testv,pts(ptipts),pts(sideipts(0)));
  if (inner_product(testv,*pforward)>0) {
    product(pforward,*pforward,-1.);
  }
}

int add_point_to_mesh(LinkedList<Array<int> > *ppolylist, LinkedList<SideFlag> *psidelist, int pttoadd, SideFlag *pside, const Array<Array<Real> > &pts) {
  int good=1;
  LinkedList<SideFlag> newsidelist;
  for (int i=0; i<pside->ipts.get_size(); i++) {
    SideFlag *pnewside=new SideFlag();
    pnewside->ipts=pside->ipts;
    pnewside->ipts(i)=pttoadd;
    sort_array(&(pnewside->ipts));
    pnewside->done=1;
    calc_forward_vector(&pnewside->forward,pnewside->ipts,pside->ipts(i),pts);
    LinkedListIterator<SideFlag> itexistside(*psidelist);
    for (; itexistside; itexistside++) {
      if ( itexistside->ipts == pnewside->ipts ) break;
    }
    newsidelist << pnewside;
    if (itexistside) {
      if (itexistside->done!=1) {
	good=0;
	// cerr << "Not adding side: done=" << itexistside->done << endl;
      }
    }
  }
  if (good) {
    LinkedListIterator<SideFlag> it(newsidelist);
    while (it) {
      LinkedListIterator<SideFlag> itexistside(*psidelist);
      for (; itexistside; itexistside++) {
	if ( itexistside->ipts == it->ipts ) break;
      }
      if (!itexistside) {
	// cerr << "side added\n";
	(*psidelist) << newsidelist.detach(it);
      }
      else {
	itexistside->done=2;
	itexistside->forward.resize(0);
	it++;
      }
    }
    pside->done=2;
    pside->forward.resize(0);
    Array<int> *ppoly=new Array<int>(pside->ipts.get_size()+1);
    (*ppoly)(0)=pttoadd;
    for (int i=0; i<pside->ipts.get_size(); i++) {
      (*ppoly)(i+1)=pside->ipts(i);
    }
    // cerr << "POLY" << (*ppoly) << endl;
    sort_array(ppoly);
    (*ppolylist) << ppoly;
    return 1;
  }
  else {
    pside->done=-1;
    pside->forward.resize(0);
    return 0;
  }
}

void create_mesh(LinkedList<Array<int> > *ppolylist, const Array<Array<Real> > &pts, int fulldim, Real maxr) {
  LinkedList<SideFlag> sidelist;
  int dim=pts(0).get_size();
  if (pts.get_size()<dim+fulldim) return;
  Array<int> *ppoly=new Array<int>(dim+fulldim);
  Array<Real> acc(pts(0));
  (*ppoly)(0)=0;
  for (int n=1; n<ppoly->get_size(); n++) {
    Array<Real> avg;
    product(&avg,acc,1./(Real)n);
    PointFromCenterIterator itpts(pts,avg);
    while (1) {
      itpts++;
      int i;
      for (i=0; i<n; i++) {
	if ((*ppoly)(i)==itpts()) break;
      }
      if (i==n) break;
    }
    (*ppoly)(n)=itpts();
    sum(&acc,acc,pts(itpts()));
  }
  {
    Array<Real> center;
    Real radius2;
    Array<Array<Real> > mypts(ppoly->get_size());
    for (int i=0; i<mypts.get_size(); i++) {
      mypts(i)=pts((*ppoly)(i));
    }
    calc_circumsphere(&center,&radius2,mypts);
    PointFromCenterIterator jtpts(pts,center);
    for (;jtpts; jtpts++) {
      Array<Real> delta;
      diff(&delta,pts(jtpts()),center);
      if (norm2(delta)>radius2) break;
      // cerr << "bad point: " << jtpts() << " " << norm(delta) << " " << sqrt(radius2) << endl;
    }
  }
  (*ppolylist) << ppoly;

  // cerr << "poly= ";
  // cerr << *ppoly << endl;

  for (int n=0; n<ppoly->get_size(); n++) {
    SideFlag *pside=new SideFlag();
    pside->ipts.resize(ppoly->get_size()-1);
    for (int i=0; i<ppoly->get_size()-1; i++) {
      pside->ipts(i)=(*ppoly)(i+(i>=n ? 1:0));
    }
    sort_array(&(pside->ipts));
    calc_forward_vector(&(pside->forward),pside->ipts,(*ppoly)(n),pts);
    pside->done=1;
    sidelist << pside;
    // cerr << "side= ";
    // cerr << pside->ipts << endl;
    // cerr << pside->forward << endl;
  }
  // debug_mesh(*ppolylist,pts);
  // debug_side(sidelist,pts);
  // cerr << "done 1st poly\n";

  int didsome;
  do {
    LinkedListIterator<SideFlag> itside(sidelist);
    didsome=0;
    for (; itside; itside++) {
      if (itside->done==1) {
	// cerr << itside->ipts << endl;
	int pttoadd=find_point_to_add(*itside,pts,maxr);
	// cerr << "add=" << pttoadd << endl;
	if (pttoadd==-1) {
	  itside->done=-1;
	}
	else {
	  add_point_to_mesh(ppolylist,&sidelist, pttoadd,itside,pts);
	  // debug_mesh(*ppolylist,pts);
	  // debug_side(sidelist,pts);
	  Real tmp;
	  didsome=1;
	}
      }
    }
  } while (didsome);
}

int remove_nearby(Array<Array<Real> > *ppts, const Array<Array<Real> > &orgpts, Real cutoff) {
  int ret=1;
  Real cutoff2=sqr(cutoff);
  LinkedList<Array<Real> > list;
  for (int i=0; i<orgpts.get_size(); i++) {
    int j;
    for (j=0; j<i; j++) {
      Array<Real> delta;
      diff(&delta,orgpts(i),orgpts(j));
      if (inner_product(delta,delta)<cutoff2) {
	ret=0;
	break;
      }
    }
    if (j==i) {
      list << new Array<Real>(orgpts(i));
    }
  }
  LinkedList_to_Array(ppts,list);
  return ret;
}

char *helpstring="Insert more help here";

int main(int argc, char *argv[]) {
  // parsing command line. See getvalue.hh for details;
  char *tabfilename="tab.in";
  //  Real maxd=0.;
  Real color=0;
  Real rmax=0.1;
  Real rmin=1e-3;
  int sigdig=5;
  char *exporter_label="vtk";
  int doline=0;
  int nofixnormal=0;
  int backface=0;
  int dohelp=0;
  int dummy=0;
  AskStruct options[]={
    {"","Simplexize " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-t","Input table file",STRINGVAL,&tabfilename},
    //{"-md","Max distance between points of the same surface",REALVAL,&maxd},
    {"-col","Color",REALVAL,&color},
    {"-rmax","Maximum radius of curvature",REALVAL,&rmax},
    {"-rmin","Minimum distance between points kept",REALVAL,&rmin},
    {"-l","Plot lines",BOOLVAL,&doline},
    //{"-nfn","do Not fix Normal vector",BOOLVAL,&nofixnormal},
    //{"-bf","generate back face",BOOLVAL,&backface},
    {"-sig","Number of significant digits printed (default: 5)",INTVAL,&sigdig},
    {"-ef","Export format (default: vtk)",STRINGVAL,&exporter_label},
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

  if (!check_plug_in(MeshExporter(),exporter_label)) {
    // cerr << exporter_label << endl;
    ERRORQUIT("Unknown export format");
  }
  MeshExporter *pexporter=GenericPlugIn<MeshExporter>::create(exporter_label);


  cout.setf(ios::fixed);
  cout.precision(sigdig);

  Array<Array<Real> > rawpts;
  if (strcmp(tabfilename,"-")==0) {
    read_table(&rawpts,cin,1);
  }
  else {
    ifstream tabfile(tabfilename);
    if (!tabfile) ERRORQUIT("Unable to open table file");
    read_table(&rawpts,tabfile,1);
  }

  // Array<Array<Real> > allpts;
  LinkedList<rVector3d> listpts3;
  LinkedList<Array<int> > allpoly;

  if (doline) {
    int idx=0;
    for (int begi=0; begi<rawpts.get_size(); begi++) {
      int endi=begi;
      int reali=idx;
      while (endi<rawpts.get_size() && rawpts(endi).get_size()>0) endi++;
      for (int i=begi; i<endi; i++) {
	rVector3d *pv=new rVector3d();
	Array_to_FixedVector(pv,rawpts(i));
	listpts3 << pv;
	idx++;
      }
      int numi=endi-begi;
      for (int i=0; i<numi; i++) {
	for (int j=i+1; j<numi; j++) {
	  Array<int> *pline=new  Array<int>(2);
	  (*pline)(0)=reali+i;
	  (*pline)(1)=reali+j;
	  allpoly << pline;
	}
      }
      begi=endi;
    }
  }
  else {
    int surfbeg=0;
    while (surfbeg<rawpts.get_size()) {
      int surfend=surfbeg;
      while (surfend<rawpts.get_size() && rawpts(surfend).get_size()>0) surfend++;
      Array<Array<Real> > pts;
      extract_elements(&pts, rawpts,surfbeg,surfend-1);
      
      remove_nearby(&pts, pts,rmin);
      LinkedList<Array<int> > poly;
      create_mesh(&poly, pts,0,rmax);
      Array<rVector3d> pts3;
      ArrayArray_to_ArrayFixedVector(&pts3,pts);
      
      combine_mesh(&listpts3,&allpoly, pts3,poly);
      // cerr << "surfend= " << surfend << endl;
      // cerr << "meshptssize= " << listpts3.get_size() << endl;
      // cerr << "meshpolysize= " << allpoly.get_size() << endl;
      surfbeg=surfend+1;
    }
  }
  Array<Real> colors(1);
  colors(0)=color;
  Array<rVector3d> allpts3;
  LinkedList_to_Array(&allpts3,listpts3);
  pexporter->write(cout,allpts3,allpoly,colors);
}
