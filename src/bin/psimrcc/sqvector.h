#ifndef _psi_src_bin_psimrcc_sqvector_h
#define _psi_src_bin_psimrcc_sqvector_h
/***************************************************************************
 *  SQ
 *  Copyright (C) 2007 by Francesco Evangelista (frank@ccc.uga.edu)
 *  A second quantization code
 ***************************************************************************/

#include <vector>
#include <map>
#include "squlli.h"
#include "sqoperator.h"

namespace psi{ namespace psimrcc{

typedef std::vector<SQULL>  DetsVec;
typedef std::vector<double> DoubleVec;
typedef std::map<SQULL,int> DetsMap;

/**
	@author Francesco Evangelista <frank@ccc.uga.edu>
*/
class SQVector{
public:
                    SQVector();
                    SQVector(SQULL& reference);
                    ~SQVector();
          void      apply_linear(double factor,SQOperator& op);
          void      apply_linear(double factor,SQOperator& op,int minexc,int maxexc);
          void      apply_exp(double factor,SQOperator& op);
          void      apply_exp(double factor,SQOperator& op,int minexc,int maxexc);
private:
          void      init();
          DoubleVec coefficients;
          DetsVec   determinants;
//   static  DetsMap   address;
  static  bool      initialized;
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_sqvector_h
