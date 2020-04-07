#ifndef _CLUS_STR_H_
#define _CLUS_STR_H_

#include "calccorr.h"
#include "stringo.h"
#include "lattype.h"
//#include <fstream.h>
#include "mpiinterf.h"

extern Real complexity_exp;

typedef Array<rVector3d> Cluster;
typedef Array<Cluster> ArrayCluster;
typedef Array<int> Arrayint;

inline Real get_length_quick(const Array<rVector3d> &c) {
  if (c.get_size()==0) return 0;
  return norm(c(c.get_size()-1)-c(0));
}

class ClusterBank {
  LinkedList<Cluster> cluster_list;
  LinkedListIterator<Cluster> current_cluster;
  int current_index;
  Real previous_length;
  AtomMultipletIterator multiplet;
  const SpaceGroup &equivalent_by_symmetry;
 public:
  ClusterBank(const rMatrix3d &cell,
              const Array<rVector3d> &atom_pos,
              int ntuple,
              const SpaceGroup &_equivalent_by_symmetry);
  ClusterBank(const ClusterBank &clusterbank,int ntuple);
  void reset(void);
  void operator++(int);
  int get_current_index(void) {
    return current_index;
  }
  Real get_previous_length(void) {
    return previous_length;
  }
  Real get_current_length(void) {
    return get_length_quick(*current_cluster);
  }
  operator Cluster& (void) {
    return *current_cluster;
  }
  LinkedList<Cluster>& get_cluster_list(void) {
    return cluster_list;
  }
};

class MultiClusterBank {
  Structure lat;
  ClusterBank cluster_bank;
  LinkedList<MultiCluster> cluster_list;
  LinkedListIterator<MultiCluster> current_cluster;
  int current_index;
  Real previous_length;
  const SpaceGroup &equivalent_by_symmetry;

  void make_new_clusters(void);
 public:
  MultiClusterBank(const Structure &_lat,
              int ntuple,
              const SpaceGroup &_equivalent_by_symmetry);
  MultiClusterBank(const MultiClusterBank &clusterbank,int ntuple);
  void reset(void);
  void operator++(int);
  int get_current_index(void) {
    return current_index;
  }
  Real get_previous_length(void) {
    return previous_length;
  }
  Real get_current_length(void) {
    return get_length_quick(current_cluster->clus);
  }
  operator MultiCluster& (void) {
    return *current_cluster;
  }
  MultiCluster* operator -> () {
    return current_cluster;
  }
  LinkedList<MultiCluster>& get_cluster_list(void) {
    return cluster_list;
  }
};

template<class T>
class StructureBank {
  int current_index;
  const SpaceGroup &equivalent_by_symmetry;
  Structure basic_structure;
  int current_volume;
  LinkedList<T> structure_list;
  LinkedListIterator<T> current_structure;
  int do2D;

  MultiDimIterator<Arrayint> curnew_config;
  Array<rMatrix3d> supercells;
  int curnew_supercell;
  LinkedListIterator<T> beginning_of_cell;

  void init(void) {
    current_volume=0;
    current_index=0;
    curnew_supercell=0;

    // special case to put pure structures first;
    LinkedListIterator<T> insert_at(structure_list);
    Structure blank_superstructure;
    blank_superstructure.cell=find_symmetric_cell(basic_structure.cell);
    find_all_atom_in_supercell(&blank_superstructure.atom_pos,
                               &blank_superstructure.atom_type,basic_structure.atom_pos,
                               basic_structure.atom_type,
                               basic_structure.cell, blank_superstructure.cell);
    Array<int> savetype=blank_superstructure.atom_type;
    for (int t=0; t<max(savetype); t++) {
      for (int s=0; s<blank_superstructure.atom_type.get_size(); s++) {
	blank_superstructure.atom_type(s)=MIN(t,savetype(s)-1);
      }
      LinkedListIterator<T> i(structure_list);
      for ( ; i ; i++) {
	if (equivalent_by_symmetry(*i,blank_superstructure)) break;
      }
      if (!i) {
	structure_list.add(new T(blank_superstructure),insert_at);
	insert_at++;
      }
    }
    reset();
  }
  int find_new_structure(void);
 public:
  StructureBank(const Structure &_basic_structure, const SpaceGroup &_equivalent_by_symmetry, int _do2D=0):
    structure_list(), current_structure(),
    equivalent_by_symmetry(_equivalent_by_symmetry),
    basic_structure(_basic_structure), do2D(_do2D), supercells(), curnew_config(), beginning_of_cell() {init();}
  StructureBank(const Structure &_basic_structure,
		const Array<Arrayint> &site_type_list,
		const SpaceGroup &_equivalent_by_symmetry, int _do2D=0):
    structure_list(), current_structure(),
    equivalent_by_symmetry(_equivalent_by_symmetry),
    basic_structure(_basic_structure), do2D(_do2D), supercells(), curnew_config(), beginning_of_cell() {
      for (int i=0; i<basic_structure.atom_type.get_size(); i++) {
        basic_structure.atom_type(i)=site_type_list(basic_structure.atom_type(i)).get_size();
      }
      init();
  }
  void set_2D_mode(int _do2D) {
    do2D=_do2D;
  }
  void reset(void) {
    current_index=0;
    current_structure.init(structure_list);
  }
  T& get_current_structure(void) {
    find_new_structure();
    return *current_structure;
  }
  operator T& (void) {
    find_new_structure();
    return *current_structure;
  }
  void operator++(int) {
    if (find_new_structure()) curend_config++;
    current_structure++;
    current_index++;
  }
  LinkedList<T>& get_structure_list(void) {
    return structure_list;
  }
  int get_current_index(void) {
    return current_index;
  }
  int get_max_volume(void) {
    return current_volume;
  }
  int end_of_this_volume(void) { // caution: this routine only works if add_structure is not called.
    if (supercells.get_size()==0) {return 0;}
    if (curnew_config) {return 0;} else {return 1;}
  }
  int find_new_structures(int new_volume) {
    while (1) {
      if (end_of_this_volume()) {
	if (current_volume==new_volume) break;
      }
      (*this)++;
    }
  } 
  int add_structure(const Structure &str, T **pp_str=NULL) {
    LinkedListIterator<T> i(structure_list);
    for ( ; i ; i++) {
      if (equivalent_by_symmetry(*i,str)) break;
    }
    if (!i) {
      LinkedListIterator<T> insert_at(structure_list);
      while (insert_at && insert_at->atom_pos.get_size() < str.atom_pos.get_size()+1) insert_at++;
      T *p_str=new T(str);
      structure_list.add(p_str,insert_at);
      if (pp_str) *pp_str=p_str;
      return 1;
    }
    else {
      if (pp_str) *pp_str=i;
      return 0;
    }
  }
};

template<class T>
int StructureBank<T>::find_new_structure(void) {
  int need_to_gen=0;
  if (!current_structure) {
    need_to_gen=1;
  }
  else if (current_structure->atom_pos.get_size()/basic_structure.atom_pos.get_size() > MAX(1,current_volume)) {
    need_to_gen=1;
  }
  if (!need_to_gen) return 0;

  while (1) {
    if (curnew_config) {
      for (int s=0; s<((Arrayint)curnew_config).get_size(); s++) {
	blank_superstructure.atom_type(s)=((Arrayint)curnew_config)(s);
      }
      if (!contains_pure_translations_or_lexico_successor(blank_superstructure,basic_structure.cell)) {
	LinkedListIterator<T> i=beginning_of_cell;
	for ( ; i ; i++) {
	  if (equivalent_by_symmetry(*i,blank_superstructure)) break;
	}
	if (!i) {
	  structure_list.add(current_structure,new T(blank_superstructure));
	  break; // exit routine;
	}
      }
    }
    else {
      curnew_supercell++;
      if (curnew_supercell>=supercells.get_size()) {
	current_volume=current_volume+1;
	Array<rMatrix3d> pointgroup;
	pointgroup_from_spacegroup(&pointgroup,equivalent_by_symmetry.point_op);
	if (do2D) {
	  find_supercells_2D(&supercells, current_volume, current_volume, basic_structure.cell, pointgroup);
	}
	else {
	  find_supercells(&supercells, current_volume, current_volume, basic_structure.cell, pointgroup);
	}
	curnew_supercell=0;
      }
      beginning_of_cell=current_structure;
      supercells(curnew_supercell)=find_symmetric_cell(supercells(curnew_supercell));
      Structure blank_superstructure;
      blank_superstructure.cell=supercells(curnew_supercell);
      find_all_atom_in_supercell(&blank_superstructure.atom_pos,
				 &blank_superstructure.atom_type,basic_structure.atom_pos,
				 basic_structure.atom_type,
				 basic_structure.cell, blank_superstructure.cell);
      curnew_config.init(blank_superstructure.atom_type);
    }
    curnew_config++;
  }
  return 1;
}

#endif
