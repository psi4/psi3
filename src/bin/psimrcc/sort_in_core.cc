/**
 *  @file ccsort_in_core.cpp
 *  @ingroup (PSIMRCC)
*/

#include <libmoinfo/libmoinfo.h>

#include "blas.h"
#include "matrix.h"
#include "sort.h"
#include "transform.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCSort::build_integrals_in_core()
{
  trans->read_integrals_from_transqt();
  frozen_core_energy_in_core();
  sort_integrals_in_core();
  trans->free_memory();
}

void CCSort::frozen_core_energy_in_core()
{
  // One electron contribution to the frozen core energy from each irrep
  efzc=0.0;
  for(int i=0;i<nfzc;i++){
    int ii=frozen_core[i];
    efzc+=2.0*trans->oei(ii,ii);
  }
  // Two electron contribution to the frozen core energy
  for(int i=0;i<nfzc;i++){
    for(int j=0;j<nfzc;j++){
      int ii=frozen_core[i];
      int jj=frozen_core[j];
      efzc+=2.0*trans->tei(ii,ii,jj,jj);
      efzc-=trans->tei(ii,jj,ii,jj);
    }
  }
}

void CCSort::sort_integrals_in_core()
{
  // Sort the TEI for CC computations
  MatrixMap matrix_map = blas->get_MatrixMap();
  for(MatrixMap::iterator iter = matrix_map.begin(); iter!=matrix_map.end(); ++iter){
    form_fock_in_core(iter);
    form_two_electron_integrals_in_core(iter);
  }
}

void CCSort::form_fock_in_core(MatrixMap::iterator& iter)
{
  CCMatrix* Matrix = iter->second;
  if(Matrix->is_fock()){
    string label     = Matrix->get_label();
    double*** matrix = Matrix->get_matrix();
    short* pq = new short[2];
    const intvec& oa2p = moinfo->get_occ_to_mo();

    bool alpha = true;
    if((label.find("O")!=string::npos) || (label.find("V")!=string::npos) || (label.find("A")!=string::npos) || (label.find("F")!=string::npos))
      alpha = false;

    // N.B. Never introduce Matrices/Vectors with O or V in the name before you compute the Fock matrix elements
    vector<int> aocc = moinfo->get_aocc(Matrix->get_reference(),AllRefs);
    vector<int> bocc = moinfo->get_bocc(Matrix->get_reference(),AllRefs);

    for(int n=0;n<moinfo->get_nirreps();n++)
      for(int i = 0;i<Matrix->get_left_pairpi(n);i++)
        for(int j = 0;j<Matrix->get_right_pairpi(n);j++){
          // Find p and q from the pairs
          Matrix->get_two_indices_pitzer(pq,n,i,j);
          // Add the h(p,q) contribution
          matrix[n][i][j]=trans->oei(pq[0],pq[1]);
          // Add the core contribution//
          for(int k=0;k<nfzc;k++){
            int kk=frozen_core[k];
            matrix[n][i][j]+=add_fock_two_in_core(pq[0],pq[1],kk,true);
            matrix[n][i][j]+=add_fock_two_in_core(pq[0],pq[1],kk,false);
          }
          for(int k=0;k<aocc.size();k++){
            int kk=oa2p[aocc[k]];
            if(alpha)
              matrix[n][i][j]+=add_fock_two_in_core(pq[0],pq[1],kk,true);
            else
              matrix[n][i][j]+=add_fock_two_in_core(pq[0],pq[1],kk,false);
          }
          for(int k=0;k<bocc.size();k++){
            int kk=oa2p[bocc[k]];
            if(!alpha)
              matrix[n][i][j]+=add_fock_two_in_core(pq[0],pq[1],kk,true);
            else
              matrix[n][i][j]+=add_fock_two_in_core(pq[0],pq[1],kk,false);
          }
        }
    delete[] pq;
  }
}

double CCSort::add_fock_two_in_core(int p, int q, int k, bool exchange)
{
  // Add the (pq|kk) contribution
  double term = trans->tei(p,q,k,k);
  // Add the -(pk|qk) contribution
  if(exchange)
    term -= trans->tei(p,k,q,k);
  return(term);
}

void CCSort::form_two_electron_integrals_in_core(MatrixMap::iterator& iter)
{
  CCMatrix* Matrix = iter->second;
  if(Matrix->is_integral()){
    short*      pqrs = new short[4];
    double*** matrix = Matrix->get_matrix();
    bool antisymmetric = Matrix->is_antisymmetric();
    if(Matrix->is_chemist()){
      for(int n=0;n<moinfo->get_nirreps();n++)
        for(int i = 0;i<Matrix->get_left_pairpi(n);i++)
          for(int j = 0;j<Matrix->get_right_pairpi(n);j++){
            Matrix->get_four_indices_pitzer(pqrs,n,i,j);
            // From (pq|rs) = <pr|qs> we define
            // (pq:rs) = <pr:qs> = (pq|rs) - (ps|qr)

            // Add the +<pr|qs> = (pq|rs) contribution
            matrix[n][i][j] += trans->tei(pqrs[0],pqrs[1],pqrs[2],pqrs[3]);

            // Add the -<pq|sr> = -(ps|qr) contribution
            if(antisymmetric)
              matrix[n][i][j] -= trans->tei(pqrs[0],pqrs[3],pqrs[1],pqrs[2]);
          }
    }else{
      for(int n=0;n<moinfo->get_nirreps();n++)
        for(int i = 0;i<Matrix->get_left_pairpi(n);i++)
          for(int j = 0;j<Matrix->get_right_pairpi(n);j++){
            Matrix->get_four_indices_pitzer(pqrs,n,i,j);
            // Add the +<pq|rs> = (pr|qs) contribution
            matrix[n][i][j] += trans->tei(pqrs[0],pqrs[2],pqrs[1],pqrs[3]);

            // Add the -<pq|sr> = -(ps|qr) contribution
            if(antisymmetric)
              matrix[n][i][j] -= trans->tei(pqrs[0],pqrs[3],pqrs[1],pqrs[2]);
          }
    }
    delete[] pqrs;
  }
}

}} /* End Namespaces */
