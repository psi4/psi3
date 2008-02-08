#ifndef _psi_src_bin_psimrcc_squlli_h
#define _psi_src_bin_psimrcc_squlli_h
/***************************************************************************
 *  SQ
 *  Copyright (C) 2007 by Francesco Evangelista (frank@ccc.uga.edu)
 *  A second quantization code
 ***************************************************************************/

#include "sqproduct.h"

#include <vector>
#include <utility>

namespace psi{ namespace psimrcc{

/**
	@author Francesco Evangelista <frank@ccc.uga.edu>
*/
class SQULL{
  typedef unsigned long long ull;
public:
  SQULL();
  ~SQULL();
  void print();
  ull bits;
  void reset()      {bits=0;}
  void reset(int n) {bits = bits & (~flags[n]);}
  void set(int n)   {bits = bits | flags[n];}
  bool test(int n)  {return (bits & flags[n] != 0);}
  SQULL operate(SQProduct& sqp) {SQULL result;return(result); // TODO Fix this
                                }
private:
  void init();
  void binary(ull number);
  static int maxbits;
  static bool initialized;
  static ull flags[128];
  static ull flods[128];
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_squlli_h
