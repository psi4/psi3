/***************************************************************************
 *  SQ
 *  Copyright (C) 2007 by Francesco Evangelista (frank@ccc.uga.edu)
 *  A second quantization code
 ***************************************************************************/

#include "sqvector.h"
#include "moinfo.h"

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

// DetsMap SQVector::address;
bool    SQVector::initialized;

using namespace std;

SQVector::SQVector()
{
  if(!initialized)
    init();
}

SQVector::SQVector(SQULL& reference)
{
  if(!initialized)
    init();
  coefficients.clear();
  determinants.clear();
  coefficients.push_back(1.0);
  determinants.push_back(reference);
}

void SQVector::init()
{
  // TODO Initialized the addressing;
}

SQVector::~SQVector()
{
}

void SQVector::apply_linear(double factor,SQOperator& op)
{
  // Copy the coefficients to a temp vector
  DoubleVec original;
  original = coefficients;
  for(int i=0;i<coefficients.size();++i)
    coefficients[i]=0.0;

  for(int i=0;i<coefficients.size();++i){
    for(int j=0;j<op.size();++j){
//       address[determinants[i].operate(op.get_sq_product(j))];
      op.get_matrix_element(j);
      op.get_sq_product(j);
    }
  }
}

//   void apply_linear(double factor,SQOperator& op,int minexc,int maxexc);
//   void apply_exp(double factor,SQOperator& op);
//   void apply_exp(double factor,SQOperator& op,int minexc,int maxexc);
}} /* End Namespaces */