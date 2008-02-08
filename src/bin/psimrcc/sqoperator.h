#ifndef _psi_src_bin_psimrcc_sqoperator_h
#define _psi_src_bin_psimrcc_sqoperator_h
/***************************************************************************
 *  SQ
 *  Copyright (C) 2007 by Francesco Evangelista (frank@ccc.uga.edu)
 *  A second quantization code
 ***************************************************************************/

#include <vector>
#include "sqproduct.h"

namespace psi{ namespace psimrcc{

/**
	@author Francesco Evangelista <frank@ccc.uga.edu>
*/
class SQOperator{
public:
  SQOperator();
  ~SQOperator();
  void Hamiltonian();
  void print();
  void add_element(double value, SQProduct& sqp) {matrix_element.push_back(value);sq_product.push_back(sqp);}
  double get_matrix_element(int n) {return(matrix_element[n]);}
  SQProduct& get_sq_product(int n) {return(sq_product[n]);}
  int  size() {return(matrix_element.size());}
private:
  std::vector<double>    matrix_element;
  std::vector<SQProduct> sq_product;
  double treshold;
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_sqoperator_h
