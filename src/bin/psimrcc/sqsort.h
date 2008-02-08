#ifndef _psi_src_bin_psimrcc_sqsort_h
#define _psi_src_bin_psimrcc_sqsort_h
/***************************************************************************
 *  SQ
 *  Copyright (C) 2007 by Francesco Evangelista (frank@ccc.uga.edu)
 *  A second quantization code
 ***************************************************************************/

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define four(i,j,k,l) INDEX(INDEX(i,j),INDEX(k,l))

#include <vector>
#include <utility>

namespace psi{ namespace psimrcc{

/**
	@author Francesco Evangelista <frank@ccc.uga.edu>
*/
class SQSort{
public:
  SQSort();
  ~SQSort();
  double get_frozen_core_energy() {return(frozen_core_energy);}
  double get_h(int p,int q) {return(h[p][q]);}
  double get_tei(int p,int q,int r,int s) {return(tei[four(p,q,r,s)]);}
private:
  void init();
  void cleanup();
  void read_one_electron_integrals();
  void read_two_electron_integrals();
  void frozen_core_energy_oei_contribution(double**& h_pitzer);
  void frozen_core_energy_tei_contribution(double*&  tei_pitzer);
  double frozen_core_energy;
  int*     ioff;
  double** h;
  double** f;
  double*  tei;
};

extern SQSort* integrals;

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_sqsort_h
