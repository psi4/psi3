#include <cmath>
#include <algorithm>

#include "memory_manager.h"
#include <libmoinfo/libmoinfo.h>
#include "transform.h"
#include "matrix.h"
#include <libutil/libutil.h>
#include "algebra_interface.h"
#include "blas.h"

#define CCTRANSFORM_USE_BLAS

#define MAX(i,j) ((i>j) ? i : j)
#define MIN(i,j) ((i>j) ? j : i)
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define four(i,j,k,l) INDEX(INDEX(i,j),INDEX(k,l))

#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include "psifiles.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

/**
 * Read at least one block of the two electron MO integrals from an iwl buffer assuming Pitzer ordering and store them in the packed array tei_mo
 */
int CCTransform::read_tei_mo_integrals_block(int first_irrep)
{
  // Read all the (frozen + non-frozen) TEI in Pitzer order
  // and store them in a in-core block-matrix
  CCIndex* mo_indexing = blas->get_index("[n>=n]");

  int last_irrep = allocate_tei_mo_block(first_irrep);

  double value;
  size_t p,q,r,s,pq,rs,pqrs,irrep;
  size_t ilsti,nbuf,fi,index,elements;
  elements = 0;
  struct iwlbuf ERIIN;
  iwl_buf_init(&ERIIN,PSIF_MO_TEI,0.0,1,1);
    do{
      ilsti = ERIIN.lastbuf;
      nbuf  = ERIIN.inbuf;
      fi    = 0;
      for(index = 0; index < nbuf; index++){
        // Compute the [pq] index for this pqrs combination
        p = abs(ERIIN.labels[fi]);
        q = ERIIN.labels[fi+1];
        r = ERIIN.labels[fi+2];
        s = ERIIN.labels[fi+3];
        value = ERIIN.values[index];
        irrep = mo_indexing->get_tuple_irrep(p,q);
        // Fill in only the blocks that fit
        if((first_irrep<=irrep) && (irrep<last_irrep)){
          pq    = mo_indexing->get_tuple_index(p,q);
          rs    = mo_indexing->get_tuple_index(r,s);
          pqrs  = INDEX(pq,rs);
          tei_mo[irrep][pqrs]=value;
        }
        fi += 4;
        elements++;
      }
      if(!ilsti)
        iwl_buf_fetch(&ERIIN);
    } while(!ilsti);
  fprintf(outfile,"\n    CCTransform: read %d non-zero integrals", elements);
  fflush(outfile);
  iwl_buf_close(&ERIIN,1);
  return(last_irrep);
}

/**
 * Allocate as many blocks of the tei_mo array and exit(EXIT_FAILURE) if there is not enough space
 */
int CCTransform::allocate_tei_mo_block(int first_irrep)
{
  if(first_irrep>moinfo->get_nirreps()){
    fprintf(outfile,"\n    CCTransform: allocate_tei_mo_block() was called with first_irrep > nirreps !");
    fflush(outfile);
    exit(EXIT_FAILURE);
  }
  
  int last_irrep = first_irrep;

  if(tei_mo==NULL){
    // Allocate the tei_mo matrix blocks
    tei_mo = new double*[moinfo->get_nirreps()];
    for(int h=0;h<moinfo->get_nirreps();h++)
      tei_mo[h] = NULL;
  }

  // Find how many irreps we can store in 95% of the free memory
  double cctransform_memory = mem->get_free_memory()*0.95;
  size_t matrix_size = 0;
  for(int h=first_irrep;h<moinfo->get_nirreps();h++){
    if(tei_mo_indexing->get_pairpi(h)>0){
      size_t block_size = INDEX(tei_mo_indexing->get_pairpi(h)-1,tei_mo_indexing->get_pairpi(h)-1)+1;
      if(to_MB(block_size)<cctransform_memory){
        matrix_size  += block_size;
        tei_mo[h]     = new double[block_size];
        zero_arr(tei_mo[h],block_size);
        cctransform_memory-=to_MB(block_size);
        mem->add_allocated_memory(to_MB(block_size));
        last_irrep++;
      }
    }else{
      last_irrep++;
    }
  }
  fprintf(outfile,"\n    CCTransform: irrep %d->%d will be read in core",first_irrep,last_irrep);
  if(first_irrep==last_irrep){
    fprintf(outfile,"\n    CCTransform: allocate_tei_mo_block() has not enough memory!");
    fflush(outfile);
    exit(EXIT_FAILURE);
  }
  first_irrep_in_core=first_irrep;
  last_irrep_in_core=last_irrep;
  return(last_irrep);
}

/**
 * Free the blocks included in the first_irrep->last_irrep range
 */
void CCTransform::free_tei_mo_integrals_block(int first_irrep, int last_irrep)
{
  for(int h=first_irrep;h<last_irrep;h++){
    if(tei_mo[h] != NULL){
      size_t block_size = INDEX(tei_mo_indexing->get_pairpi(h)-1,tei_mo_indexing->get_pairpi(h)-1)+1;
      delete[] tei_mo[h];
      mem->add_allocated_memory(-to_MB(block_size));
    }
  }
  if(last_irrep>=moinfo->get_nirreps()){
    delete[] tei_mo;
    tei_mo = NULL;
  }
}

double CCTransform::tei_block(int p, int q, int r, int s)
{
  // Get the (pq|rs) integral
  int irrep(tei_mo_indexing->get_tuple_irrep(MAX(p,q),MIN(p,q)));
  if((first_irrep_in_core <= irrep) && (irrep < last_irrep_in_core))
    return(tei_mo[tei_mo_indexing->get_tuple_irrep(MAX(p,q),MIN(p,q))][INDEX(tei_mo_indexing->get_tuple_index(MAX(p,q),MIN(p,q)),tei_mo_indexing->get_tuple_index(MAX(r,s),MIN(r,s)))]);
  else
    return(0.0);
}

}} /* End Namespaces */
