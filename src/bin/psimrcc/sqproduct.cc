/***************************************************************************
 *  SQ
 *  Copyright (C) 2007 by Francesco Evangelista (frank@ccc.uga.edu)
 *  A second quantization code
 ***************************************************************************/

#include "moinfo.h"
#include "sqproduct.h"

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

SQProduct::SQProduct()
{
}

SQProduct::~SQProduct()
{
}

void SQProduct::add_alpha_creator(int mo)
{
  operators.push_back(make_pair(mo,true));
}

void SQProduct::add_beta_creator(int mo)
{
  operators.push_back(make_pair(mo,true));
}

void SQProduct::add_alpha_annihilator(int mo)
{
  operators.push_back(make_pair(mo,false));
}

void SQProduct::add_beta_annihilator(int mo)
{
  operators.push_back(make_pair(mo,false));
}

void SQProduct::print()
{
  if(operators.size()>0){
    for(int i = operators.size()-1;i>=0;i--){
      fprintf(outfile,"a%s(%d)",(operators[i].second ? "+" : "-" ),operators[i].first);
    }
  }else{
    fprintf(outfile,"I");
  }
  

}
}} /* End Namespaces */