/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <cmath>
#include "moinfo.h"
#include "sqoperator.h"
#include "sqsort.h"

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

SQOperator::SQOperator()
{
  treshold = 1.0e-20;
}

SQOperator::~SQOperator()
{
}

void SQOperator::Hamiltonian()
{
  SQProduct sqp;

  // H0
  sqp.zero();
  sq_product.push_back(sqp);
  matrix_element.push_back(integrals->get_frozen_core_energy());

  // H1 alpha-alpha
  for(int p = 0;p<moinfo->get_nmo();p++){
      for(int q = 0;q<moinfo->get_nmo();q++){
          if(fabs(integrals->get_h(p,q))>treshold){
              sqp.zero();
              sqp.add_alpha_annihilator(q);
              sqp.add_alpha_creator(p);
              sq_product.push_back(sqp);
              matrix_element.push_back(integrals->get_h(p,q));
          }
      }
  }
  // H1 beta-beta
  for(int p = 0;p<moinfo->get_nmo();p++){
      for(int q = 0;q<moinfo->get_nmo();q++){
          if(fabs(integrals->get_h(p,q))>treshold){
              sqp.zero();
              sqp.add_beta_annihilator(q);
              sqp.add_beta_creator(p);
              sq_product.push_back(sqp);
              matrix_element.push_back(integrals->get_h(p,q));
          }
      }
  }

  // H2
  for(int p = 0;p<moinfo->get_nmo();p++){
    for(int q = 0;q<moinfo->get_nmo();q++){
      for(int r = 0;r<moinfo->get_nmo();r++){
        for(int s = 0;s<moinfo->get_nmo();s++){
          if(fabs(integrals->get_tei(p,q,r,s))>treshold){
            if((s!=q)&&(r!=p)){
              sqp.zero();
              sqp.add_alpha_annihilator(q);
              sqp.add_alpha_annihilator(s);
              sqp.add_alpha_creator(r);
              sqp.add_alpha_creator(p);

              sq_product.push_back(sqp);
              matrix_element.push_back(0.5*integrals->get_tei(p,q,r,s));
            }

            sqp.zero();
            sqp.add_alpha_annihilator(s);
            sqp.add_beta_annihilator(q);
            sqp.add_beta_creator(p);
            sqp.add_alpha_creator(r);

            sq_product.push_back(sqp);
            matrix_element.push_back(integrals->get_tei(p,q,r,s));

            if((s!=q)&&(r!=p)){
              sqp.zero();
              sqp.add_beta_annihilator(q);
              sqp.add_beta_annihilator(s);
              sqp.add_beta_creator(r);
              sqp.add_beta_creator(p);

              sq_product.push_back(sqp);
              matrix_element.push_back(0.5*integrals->get_tei(p,q,r,s));
            }
          }
        }
      }
    }
  }
};


void SQOperator::print()
{
  fprintf(outfile,"\n\n  Printing Operator");
  for(int element = 0; element < matrix_element.size() ; element++){
    fprintf(outfile,"\n    [%5d] %16.12f * ",element,matrix_element[element]);
    sq_product[element].print();
  }
  fflush(outfile);
}

// void SQOperator::H()
// {
//   
//   int p,q,r,s;
//   long int *ioff;
//   double **h,*tei;
//   Operator    op;
// 
//   ioff    =ints.getIoff();
//   tei     =ints.getTei();
//   h       =ints.getH();
//   // H0
//   op.zero();
//   set.push_back(op);
//   t.push_back(ints.getCoreEnergy());
//   nset++;
// 
//   // H1 alpha-alpha
//   for(p=0;p<moinfo.nmo;p++){
//       for(q=0;q<moinfo.nmo;q++){
//           if(fabs(h[p][q])>t_treshold){
//               op.zero();
//               op.addAlphaAnnihilator(q);
//               op.addAlphaCreator(p);
//               set.push_back(op);
//               t.push_back(h[p][q]);
//               nset++;
//           }
//       }
//   }
//   // H1 beta-beta
//   for(p=0;p<moinfo.nmo;p++){
//       for(q=0;q<moinfo.nmo;q++){
//           if(fabs(h[p][q])>t_treshold){
//               op.zero();
//               op.addBetaAnnihilator(q);
//               op.addBetaCreator(p);    
//               set.push_back(op);
//               t.push_back(h[p][q]);
//               nset++;
//           }
//       }
//   }
//   // H2 alpha-alpha
//   for(p=0;p<moinfo.nmo;p++){
//       for(q=0;q<moinfo.nmo;q++){
//           for(r=0;r<moinfo.nmo;r++){
//               for(s=0;s<moinfo.nmo;s++){
//                   if(fabs(tei[four(p,q,r,s)])>t_treshold){
//                       if((s!=q)&&(r!=p)){
//                       op.zero();
//                       op.addAlphaAnnihilator(q); 
//                       op.addAlphaAnnihilator(s); 
//                       op.addAlphaCreator(r);     
//                       op.addAlphaCreator(p);     
// 
//                       set.push_back(op);
//                       t.push_back(0.5*tei[four(p,q,r,s)]);
//                       nset++;}
// 
//                       op.zero();
//                       op.addBetaAnnihilator(q); 
//                       op.addAlphaAnnihilator(s); 
//                       op.addAlphaCreator(r);     
//                       op.addBetaCreator(p);     
// 
//                       set.push_back(op);
//                       t.push_back(1.0*tei[four(p,q,r,s)]);
//                       nset++;
// 
//                       // The ollowing code is shaded because I introduced a factor of two in the previous expression. This simpplification is allowed when using RHF ,ROHF and CAS orbitals.
//                       /*op.zero();
//                       op.addAlphaAnnihilator(q); 
//                       op.addBetaAnnihilator(s); 
//                       op.addBetaCreator(r);     
//                       op.addAlphaCreator(p);     
// 
//                       set.push_back(op);
//                       t.push_back(0.5*tei[four(p,q,r,s)]);
//                       nset++;*/
// 
//                       if((s!=q)&&(r!=p)){
//                       op.zero();
//                       op.addBetaAnnihilator(q); 
//                       op.addBetaAnnihilator(s); 
//                       op.addBetaCreator(r);     
//                       op.addBetaCreator(p);     
// 
//                       set.push_back(op);
//                       t.push_back(0.5*tei[four(p,q,r,s)]);
//                       nset++;}
// 
//                   }
//               }
//           }
//       }
//   }
// }
}} /* End Namespaces */