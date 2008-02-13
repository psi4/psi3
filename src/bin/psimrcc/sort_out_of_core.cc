/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
/**
 *  @file ccsort_out_of_core.cpp
 *  @ingroup PSIMRCC
*/

#include "memory_manager.h"
#include "transform.h"
#include "sort.h"
#include "matrix.h"
#include "moinfo.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

/**
 * Builds the integral matrices on disk using an out-of-core algorithm
 */
void CCSort::build_integrals_out_of_core()
{
  trans->read_oei_from_transqt();
  MatrixMap matrix_map= blas->get_MatrixMap();
  MatMapIt mat_it     = matrix_map.begin();
  MatMapIt mat_end    = matrix_map.end();
  int mat_irrep = 0;
  int cycle     = 0;
  while(mat_it!=mat_end){
    fprintf(outfile,"\n\n  CCSort::build_integrals_out_of_core(): cycle = %d",cycle);
    // Find how many matrices blocks we can store in 95% of the free memory and allocate them
    MatrixBlks to_be_processed;
    setup_out_of_core_list(mat_it,mat_irrep,mat_end,to_be_processed);
//     fprintf(outfile,"\nCCSort::the following %d matrix blocks will be processed:",to_be_processed.size());
//     for(MatBlksIt block_it=to_be_processed.begin();block_it!=to_be_processed.end();++block_it)
//       fprintf(outfile,"\n%s[%s]",block_it->first->get_label().c_str(),moinfo->get_irr_labs(block_it->second));
    int first_irrep = 0;
    int last_irrep  = 0;
    while(last_irrep<moinfo->get_nirreps()){
      last_irrep = trans->read_tei_mo_integrals_block(first_irrep);
      fprintf(outfile,"\n    CCSort: sorting block(s) %d->%d",first_irrep,last_irrep);
      if(cycle==0) frozen_core_energy_out_of_core();
      sort_integrals_out_of_core(first_irrep,last_irrep,to_be_processed);
      trans->free_tei_mo_integrals_block(first_irrep,last_irrep);
    }
    // Write to disk and free memory
    dump_integrals_to_disk(to_be_processed);
    cycle++;
  }
}

void CCSort::setup_out_of_core_list(MatMapIt& mat_it,int& mat_irrep,MatMapIt& mat_end,MatrixBlks& to_be_processed)
{
  double ccintegrals_memory = mem->get_free_memory()*0.5;
  bool out_of_memory = false;
  while((mat_it != mat_end) && !out_of_memory){
    if(mat_it->second->is_integral() || mat_it->second->is_fock()){
      CCMatrix* Matrix = mat_it->second;
      while(mat_irrep < moinfo->get_nirreps() && !out_of_memory){
        size_t block_size = Matrix->get_block_sizepi(mat_irrep);
        if(to_MB(block_size) < ccintegrals_memory){
          to_be_processed.push_back(make_pair(Matrix,mat_irrep));
          // Allocate the matrix, this will also take care of MOInfo::allocated_memory
          Matrix->allocate_block(mat_irrep);
          ccintegrals_memory -= to_MB(block_size);
          mat_irrep++;
        }else{
            out_of_memory = true;
        }
      }
      if(!out_of_memory)
        mat_irrep=0;
    }
    if(!out_of_memory)
      ++mat_it;
  }
}

void CCSort::sort_integrals_out_of_core(int first_irrep, int last_irrep, MatrixBlks& to_be_processed)
{
  for(MatBlksIt block_it=to_be_processed.begin();block_it!=to_be_processed.end();++block_it){
    form_fock_out_of_core(block_it->first,block_it->second);
    form_two_electron_integrals_out_of_core(block_it->first,block_it->second);
  }
}


void CCSort::form_two_electron_integrals_out_of_core(CCMatrix* Matrix, int h)
{
  if(Matrix->is_integral()){
    short*      pqrs = new short[4];
    double*** matrix = Matrix->get_matrix();
    bool antisymmetric = Matrix->is_antisymmetric();
    if(Matrix->is_chemist()){
      for(int i = 0;i<Matrix->get_left_pairpi(h);i++)
        for(int j = 0;j<Matrix->get_right_pairpi(h);j++){
          Matrix->get_four_indices_pitzer(pqrs,h,i,j);
          // From (pq|rs) = <pr|qs> we define
          // (pq:rs) = <pr:qs> = (pq|rs) - (ps|qr)

          // Add the +<pr|qs> = (pq|rs) contribution
          matrix[h][i][j] += trans->tei_block(pqrs[0],pqrs[1],pqrs[2],pqrs[3]);

          // Add the -<pq|sr> = -(ps|qr) contribution
          if(antisymmetric)
            matrix[h][i][j] -= trans->tei_block(pqrs[0],pqrs[3],pqrs[1],pqrs[2]);
        }
    }else{
      for(int i = 0;i<Matrix->get_left_pairpi(h);i++)
        for(int j = 0;j<Matrix->get_right_pairpi(h);j++){
          Matrix->get_four_indices_pitzer(pqrs,h,i,j);
          // Add the +<pq|rs> = (pr|qs) contribution
          matrix[h][i][j] += trans->tei_block(pqrs[0],pqrs[2],pqrs[1],pqrs[3]);

          // Add the -<pq|sr> = -(ps|qr) contribution
          if(antisymmetric)
            matrix[h][i][j] -= trans->tei_block(pqrs[0],pqrs[3],pqrs[1],pqrs[2]);
        }
    }
    delete[] pqrs;
  }
}

void CCSort::form_fock_out_of_core(CCMatrix* Matrix, int h)
{
  if(Matrix->is_fock()){
    string label     = Matrix->get_label();
    double*** matrix = Matrix->get_matrix();
    short* pq = new short[2];
    int* oa2p = moinfo->get_occ_to_pitzer();

    bool alpha = true;
    if((label.find("O")!=string::npos) || (label.find("V")!=string::npos) || (label.find("A")!=string::npos)) // NB This was missing the last bit, this might be a problem
      alpha = false;

    // N.B. Never introduce Matrices/Vectors with O or V in the name before you compute the Fock matrix elements
    vector<int> aocc = moinfo->get_aocc("a",Matrix->get_reference());
    vector<int> bocc = moinfo->get_bocc("a",Matrix->get_reference());

    for(int i = 0;i<Matrix->get_left_pairpi(h);i++)
      for(int j = 0;j<Matrix->get_right_pairpi(h);j++){
        // Find p and q from the pairs
        Matrix->get_two_indices_pitzer(pq,h,i,j);
        // Add the h(p,q) contribution
        matrix[h][i][j]=trans->oei(pq[0],pq[1]);
        // Add the core contribution//
        for(int k=0;k<nfzc;k++){
          int kk=frozen_core[k];
          matrix[h][i][j]+=add_fock_two_in_core(pq[0],pq[1],kk,true);
          matrix[h][i][j]+=add_fock_two_in_core(pq[0],pq[1],kk,false);
        }
        for(int k=0;k<aocc.size();k++){
          int kk=oa2p[aocc[k]];
          if(alpha)
            matrix[h][i][j]+=add_fock_two_in_core(pq[0],pq[1],kk,true);
          else
            matrix[h][i][j]+=add_fock_two_in_core(pq[0],pq[1],kk,false);
        }
        for(int k=0;k<bocc.size();k++){
          int kk=oa2p[bocc[k]];
          if(!alpha)
            matrix[h][i][j]+=add_fock_two_in_core(pq[0],pq[1],kk,true);
          else
            matrix[h][i][j]+=add_fock_two_in_core(pq[0],pq[1],kk,false);
        }
      }
    delete[] pq;
  }
}

double CCSort::add_fock_two_out_of_core(int p, int q, int k, bool exchange)
{
  // Add the (pq|kk) contribution
  double term = trans->tei_block(p,q,k,k);
  // Add the -(pk|qk) contribution
  if(exchange)
    term -= trans->tei_block(p,k,q,k);
  return(term);
}

void CCSort::frozen_core_energy_out_of_core()
{
  // One electron contribution to the frozen core energy from each irrep
  efzc=0.0;
  for(int i=0;i<nfzc;i++){
    int ii=frozen_core[i];
    efzc+=2*trans->oei(ii,ii);
  }
  // Two electron contribution to the frozen core energy
  for(int i=0;i<nfzc;i++){
    for(int j=0;j<nfzc;j++){
      int ii=frozen_core[i];
      int jj=frozen_core[j];
      efzc+=2.0*trans->tei_block(ii,ii,jj,jj);
      efzc-=trans->tei_block(ii,jj,ii,jj);
    }
  }
}

/**
 * Dump (write + free memory) the two electron integral matrix blocks contained in the list to_be_processed
 * @param to_be_processed
 */
void CCSort::dump_integrals_to_disk(MatrixBlks& to_be_processed)
{
  for(MatBlksIt block_it=to_be_processed.begin();block_it!=to_be_processed.end();++block_it){
    CCMatrix* Matrix = block_it->first;
    Matrix->dump_block_to_disk(block_it->second);
  }
}

}} /* End Namespaces */
