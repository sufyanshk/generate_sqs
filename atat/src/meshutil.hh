#ifndef _MESHUTIL_H_
#define _MESHUTIL_H_

#include <fstream.h>
#include "arraylist.h"
#include "stringo.h"
#include "plugin.h"

class MeshExporter {
public:
  virtual void write(ostream &file, const Array<rVector3d> &pts, const LinkedList<Array<int> > &poly, const Array<Real> &ptsdata) {}
  void write(ostream &file, const Array<rVector3d> &pts, const LinkedList<iVector3d> &poly, const Array<Real> &ptsdata);
  void write(ostream &file, const Array<Array<Real> > &pts, const LinkedList<Array<int> > &poly, const Array<Real> &ptsdata);
};

class VTKMeshExporter: public MeshExporter {
  virtual void write(ostream &file, const Array<rVector3d> &pts, const LinkedList<Array<int> > &poly, const Array<Real> &ptsdata);
};

void convert_to_rVector3d(Array<rVector3d> *ppts, const Array<Array<Real> > pts);
void combine_mesh(LinkedList<rVector3d> *ppts, LinkedList<Array<int> > *ppoly, const Array<rVector3d> &newpts, const LinkedList<Array<int> > &newpoly);
void transform(Array<rVector3d> *ptpts, const Array<rVector3d> pts, const rMatrix3d &op, const rVector3d &trans);

template<class T>
void ArrayArray_to_ArrayFixedVector(Array<Vector3d<T> > *paf, const Array<Array<T> > &aa) {
  paf->resize(aa.get_size());
  for (int i=0; i<aa.get_size(); i++) {
    Array_to_FixedVector(&((*paf)(i)),aa(i));
  }
}

template<class T>
void LinkedListFixedVector_to_LinkedListArray(LinkedList<Array<T> > *pla, const LinkedList<Vector3d<T> > &lf) {
  LinkedListIterator<Vector3d<T> > it(lf);
  for (;it; it++) {
    Array<T> *pa=new Array<T>();
    FixedVector_to_Array(pa,*it);
    (*pla) << pa;
  }
}


class PolyFont3D {
  Array<Array<rVector3d> > charpts;
  Array<LinkedList<Array<int> > > charpoly;
  Array<Real> charwidth;
  Real spacing;
public:
  PolyFont3D(void): charpts(256),charpoly(256),charwidth(256),spacing(0.1) {}
  PolyFont3D(ifstream &file): charpts(256),charpoly(256),charwidth(256),spacing(0.1) {init(file);}
  int init(ifstream &file);
  void write(LinkedList<rVector3d> *plistpts, LinkedList<Array<int> > *plistpoly, const rVector3d &where, const rVector3d &right, const rVector3d &up, char *string);
  Real get_length(char *string, int includelast=0);
};

#endif
