#ifndef _psi_src_bin_psimrcc_sqproduct_h
#define _psi_src_bin_psimrcc_sqproduct_h
/***************************************************************************
 *  SQ
 *  Copyright (C) 2007 by Francesco Evangelista (frank@ccc.uga.edu)
 *  A second quantization code
 ***************************************************************************/

#include <vector>
#include <utility>

namespace psi{ namespace psimrcc{

/**
	@author Francesco Evangelista <frank@ccc.uga.edu>
*/
class SQProduct{
public:
  SQProduct();
  ~SQProduct();
  void    print();
  void    add_alpha_creator(int mo);
  void    add_beta_creator(int mo);
  void    add_alpha_annihilator(int mo);
  void    add_beta_annihilator(int mo);
  void    zero() {operators.clear();}
private:
  std::vector<std::pair<int, bool> > operators;
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_sqproduct_h
