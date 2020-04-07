#ifndef _MCLIB_H_
#define _MCLIB_H_

#include <fstream.h>
#include <time.h>
#include "clus_str.h"
#include "linalg.h"
#include "equil.h"
#include "kmeci.h"

#define SPIN_TYPE signed char

class FlipInfo {
 public:
  iVector3d cell;
  int incell;
  int spin;
  FlipInfo(void): cell() {}
};

class MultiMonteCarloState {
 public:
  Real cur_energy;
  Real cur_disorder_param;
  Array<Real> cur_rho;
  LinkedList<FlipInfo> saved_spins;
 public:
  MultiMonteCarloState(void): cur_rho(), saved_spins() {cur_energy=0; cur_disorder_param=0;}
};

class MultiMonteCarlo {
 protected:
  Structure lattice;
  Array<Array<int> > site_type_list;
  iVector3d supercell;
  iVector3d margin;
  iVector3d total_box;
  SPIN_TYPE *spin;
  SPIN_TYPE *spin_orig;
  Real rspin_size;
  int site_in_cell;
  int active_site_in_cell;
  int *which_site;
  int *nb_spin_val;
  int total_clusters;
  int *nb_clusters;
  //  int *nb_x_clusters;
  int **cluster_size;
  Real *rcluster_mult_per_atom;
  int ***site_offset;
  Real ****spin_val_clus;
  Real **eci;
  int **which_cluster;
  Real E_ref;
  int which_is_empty;
  int nb_point;
  int *which_is_point;
  Real *point_mult;
  Real cur_energy;
  Array<Real> cur_rho;
  Real *pcur_rho;
  Real *new_rho;
  Real cur_disorder_param;
  Array<Real> cur_conc;

  Real T;
  Array<Real> mu;
  Real *pmu;

  Array<Array<Array<Real> > > corrfunc;

  Array<Array<int> > allowed_flip_site;
  Array<Array<int> > allowed_flip_before;
  Array<Array<Array<int> > > allowed_flip_after;
  iVector3d flip_span;

#ifdef ATAT_MPI
  iVector3d parbox_subdiv;
  int parbox_step;
  Array<Array<iVector3d> > parbox_org;
  Array<Array<iVector3d> > parbox_size;
#endif

 protected:
  void calc_delta_point_corr(Array<Real> *pcorr, int site, int type);
  void calc_delta_point_corr(Array<Real> *pcorr, const Array<int> &sites,const Array<int> &types);
 public:
  void find_all_allowed_flips(const Array2d<Real> &corr_constraints, const iVector3d &_flip_span);
  int force_spin_flip(int *cell, int incell, int newspin);
  void save_state(MultiMonteCarloState *pstate);
  void save_spin(MultiMonteCarloState *pstate, int *cell, int incell);
  void restore_state(const MultiMonteCarloState &state);
  int access(int *cell, int incell);
  void multi_spin_flip(void);
 public:
  void calc_from_scratch(void);
  void spin_flip(void);
  void run(int mc_passes, int mode);
  MultiMonteCarlo(const Structure &_lattice, const Array<Array<int> > &_site_type_list, const iVector3d &_supercell, 
		const SpaceGroup &space_group, const LinkedList<MultiCluster> &cluster_list, 
		const Array<Array<Array<Real> > > &_corrfunc);
  ~MultiMonteCarlo(void);
  void set_eci(const Array<Real> &eci);
  void init_random(const Array<Array<Real> > &conc);
  void init_structure(const Structure &str);
  //  void set_concentration(Real concentration);
  void set_T_mu(Real _T, const Array<Real> &_mu);
  void view(const Array<Arrayint> &labellookup, const Array<AutoString> &atom_label, ofstream &file, const rMatrix3d &axes);

  const Array<Real> & get_cur_corr(void) const {
    return cur_rho;
  }
  Real get_cur_energy(void) {
    return E_ref+cur_energy+extension_get_energy();
  }
  const Array<Real> & get_cur_concentration(void);
  Real get_cur_disorder_param(void) const {
    return cur_disorder_param;
  }
  int get_total_clusters(void) const {
    return total_clusters;
  }
  void get_cluster_mult(Array<Real> *mult) const {
    mult->resize(total_clusters);
    for (int i=0; i<total_clusters; i++) {(*mult)(i)=rcluster_mult_per_atom[i];}
  }
  Real get_T(void) const {return T;}
  const Array<Real> &get_mu(void) const {return mu;}
  const iVector3d & get_cell_size() const {
    return supercell;
  }
  void get_thermo_data(Array<Real> *pdata);
#ifdef ATAT_MPI
  void set_parallel(const iVector3d &_parbox_subdiv, int _parbox_step) {
    parbox_subdiv=_parbox_subdiv;
    parbox_step=_parbox_step;
    iVector3d region_size=supercell;
    for (int i=0; i<3; i++) {
      region_size(i)/=parbox_subdiv(i);
    }
    parbox_org.resize(16);
    parbox_size.resize(16);
    int pass=0;
    for (int halfshift=0; halfshift<2; halfshift++) {
      MultiDimIterator<iVector3d> subbox(iVector3d(2,2,2));
      for (; subbox; subbox++) {
	parbox_org(pass).resize(parbox_subdiv);
	iVector3d vsubbox=subbox;
	for (int i=0; i<3; i++) {
	  vsubbox(i);
	}
	pass++;
      }
    }
  }
#endif
 protected:
  virtual void extension_calc_from_scratch(void) {}
  virtual void extension_save_state(void) {}
  virtual void extension_forget_state(void) {}
  virtual void extension_update_spin_flip(int *cell, int incell, int oldspin, int newspin) {}
  virtual void extension_undo_spin_flip(void) {}
  virtual Real extension_get_energy(void) {return 0;}
};


void run_mc(GenericAccumulator *accum, MultiMonteCarlo *pmc, int mode);

class FlippedSpin {
 public:
  FlippedSpin(int *_cell, int _incell, const Array<Real> &_dspin) {
	for (int i=0; i<3; i++) {cell[i]=_cell[i];}
	incell=_incell; dspin=_dspin;
  }
  int cell[3];
  int incell;
  Array<Real> dspin;
};

class KSpaceMultiMonteCarlo: public MultiMonteCarlo {
 protected:
  int nsite;
  int size;
  Real rsize;
  Array<Array<Array<Complex> > > ft_spin;
  Array<Array<Array<Array<Array<Complex> > > > > ft_eci;
  Array<Array<Array<Array<Array<Complex> > > > > dir_eci;
  Array<Array<Array<Complex> > > convol;
  KSpaceECI *p_kspace_eci;
  Real cur_recip_E;
  clock_t fft_time;
  clock_t flip_time;
  Array<Real> ref_x;
  Real threshold_dx;

  LinkedList<FlippedSpin> flipped_spins;
  FlippedSpin *psave_spins;
  Real save_recip_E;

 public:  
  KSpaceMultiMonteCarlo(const Structure &_lattice, const Array<Array<int> > &_site_type_list, const iVector3d &_supercell,
		const SpaceGroup &space_group, const LinkedList<MultiCluster> &cluster_list,
		const Array<Array<Array<Real> > > &_corrfunc, KSpaceECI *_p_kspace_eci);
  ~KSpaceMultiMonteCarlo();
 protected:
  void set_k_space_eci(void);
  void extension_calc_from_scratch(void);
  void extension_save_state(void);
  void extension_forget_state(void);
  void extension_update_spin_flip(int *cell, int incell, int oldspin, int newspin);
  void extension_undo_spin_flip(void);
  Real extension_get_energy(void) {return cur_recip_E;}
};

#endif
