/***************************************************************************
 *  SQ
 *  Copyright (C) 2007 by Francesco Evangelista (frank@ccc.uga.edu)
 *  A second quantization code
 ***************************************************************************/

#include "sq.h"
#include "sqoperator.h"
#include "sqsort.h"

namespace psi{ namespace psimrcc{

SQSort* integrals;

void run_sq()
{
  integrals = new SQSort();

  SQOperator H;
  H.Hamiltonian();
  H.print();
  
  delete integrals;
}

}} /* End Namespaces */