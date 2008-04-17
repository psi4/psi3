#include <iostream>
#include <cmath>

#include <libmoinfo/libmoinfo.h>
#include <libciomr/libciomr.h>

#include "scf.h"
#include "memory_manager.h"
#include "sblock_matrix.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

double SCF::energy(int cycle,double old_energy)
{
  double electronic_energy = 0.0;

  T  = H;
  T += Fc;
  electronic_energy += dot(Dc,T);

  if(reference == rohf){
    T  = H;
    T.scale(0.5);
    T += Fo;
    electronic_energy += dot(Do,T);
  }

  total_energy = electronic_energy + moinfo_scf->get_nuclear_energy();

  if(reference == tcscf){
//     SBlockMatrix Dtc_sum("Dtc sum",nirreps,sopi,sopi);
// 
//     // Compute diagonal elements of H
//     for(int I = 0 ; I < nci; ++I){
//       Dtc_sum  = Dc;
//       Dtc_sum += Dtc[I];
//       construct_G(Dtc_sum,G,"PK");
//       T  = H;
//       T.scale(2.0);
//       T += G;
//       H_tcscf[I][I] = dot(Dtc_sum,T) + moinfo_scf->get_nuclear_energy();
//     }
// 
//     // Compute off-diagonal elements of H
//     for(int I = 0 ; I < nci; ++I){
//       for(int J = I + 1; J < nci; ++J){
//         construct_G(Dtc[I],G,"K");
//         H_tcscf[I][J] = H_tcscf[J][I] = - dot(Dtc[J],G);
//       }
//     }

//     fprintf(outfile,"\n  Hamiltonian");
//     for(int I = 0 ; I < nci; ++I){
//       fprintf(outfile,"\n    ");
//       for(int J = 0 ; J < nci; ++J)
//         fprintf(outfile," %11.8f ",H_tcscf[I][J]);
//     }
    

    // Compute the CI gradient
    norm_ci_grad = 0.0;
    for(int I = 0 ; I < nci; ++I){
      ci_grad[I] = 0.0;
      for(int J = 0; J < nci; ++J){
        ci_grad[I] += H_tcscf[I][J] * ci[J];
      }
      ci_grad[I] -= old_energy * ci[I];
      norm_ci_grad += fabs(ci_grad[I]);
    }
    
    double*  eigenvalues;
    double** eigenvectors;
    allocate1(double,eigenvalues,nci);
    allocate2(double,eigenvectors,nci,nci);

    sq_rsp(nci,nci,H_tcscf,eigenvalues,1,eigenvectors,1.0e-14);
  
    total_energy = eigenvalues[root];

    for(int I = 0 ; I < nci; ++I)
      ci[I] = eigenvectors[I][root]; 

    release1(eigenvalues);
    release2(eigenvectors);
  }

  fprintf(outfile,"\n  @SCF %4d  %20.12f  %20.12f",cycle,total_energy,total_energy-old_energy);

  if(reference == tcscf){
    fprintf(outfile,"\n    ci      = [");
    for(int I = 0 ; I < nci; ++I)
      fprintf(outfile,"%11.8f%s",ci[I], I != nci -1 ? "," : "");
    fprintf(outfile,"]");

    fprintf(outfile,"\n    ci_grad = [");
    for(int I = 0 ; I < nci; ++I)
      fprintf(outfile,"%11.8f%s",ci_grad[I], I != nci -1 ? "," : "");
    fprintf(outfile,"]");
  }

  fflush(outfile);

  return(total_energy); 
}

// double SCF::compute_energy_rohf(int cycle,double old_energy)
// {
//   double electronic_energy = 0.0;
//   double total_energy = 0.0;
//   for(int h =0; h < nirreps; ++h){
//     double** d_closed = Dc->get_block(h);
//     double** f_closed = Fc->get_block(h);
//     double** d_open   = Do->get_block(h);
//     double** f_open   = Fo->get_block(h);
//     double** h0 = H->get_block(h);
//     for(int m = 0; m < sopi[h]; ++m){
//       for(int n = 0; n < sopi[h]; ++n){
//         electronic_energy +=   d_closed[m][n] * (       h0[m][n] + f_closed[m][n])
//                              + d_open[m][n] *   ( 0.5 * h0[m][n] + f_open[m][n]);
//       }
//     }
//   }
//   total_energy = electronic_energy + moinfo_scf->get_nuclear_energy();
//   fprintf(outfile,"\n  @SCF %4d  %20.12f  %20.12f",cycle,total_energy,total_energy-old_energy);
//   fflush(outfile);
//   return(total_energy); 
// }

// double SCF::compute_energy_tcscf(int cycle,double old_energy)
// {
//   double electronic_energy = 0.0;
//   double total_energy = 0.0;
//   
//   
// 
//   for(int k = 0; k < 2; ++k){
//     double matrix_element = 0.0;
//     *Dtc_sum  = *Dc;
//     *Dtc_sum += *Dtc[k];
//     construct_G(Dtc_sum,G,PK);
//     for(int h =0; h < nirreps; ++h){
//       double** d  = Dtc_sum->get_block(h);
//       double** g  = G->get_block(h);
//       double** h0 = H->get_block(h);
//       for(int m = 0; m < sopi[h]; ++m){
//         for(int n = 0; n < sopi[h]; ++n){
//           matrix_element +=   d[m][n] * ( 2.0 * h0[m][n] + g[m][n]);
//         }
//       }
//       H_tcscf[k][k] = matrix_element + moinfo_scf->get_nuclear_energy();
//     }
//   }
//   G->zero();
//   construct_G(Dtc[0],G,K);
//   double matrix_element = 0.0;
//   for(int h =0; h < nirreps; ++h){
//       double** d  = Dtc[1]->get_block(h);
//       double** g  = G->get_block(h);
//       for(int m = 0; m < sopi[h]; ++m){
//         for(int n = 0; n < sopi[h]; ++n){
//           matrix_element +=   d[m][n] * g[m][n];
//         }
//       }
//       H_tcscf[0][1] = H_tcscf[1][0] = - matrix_element;
//     }
// 
//    fprintf(outfile,"\n  Hamiltonian");
//    fprintf(outfile,"\n  %20.12f  %20.12f",H_tcscf[0][0],H_tcscf[0][1]);
//    fprintf(outfile,"\n  %20.12f  %20.12f",H_tcscf[1][0],H_tcscf[1][1]);
// 
// 
// //   total_energy = electronic_energy + moinfo_scf->get_nuclear_energy();
//   fprintf(outfile,"\n  @SCF %4d  %20.12f  %20.12f  %15.12f  %15.12f",cycle,total_energy,total_energy-old_energy,ci[0],ci[1]);
//   fflush(outfile);
// 
//   delete Dtc_sum;
//   return(total_energy); 
// }

}} /* End Namespaces */

