/*
** INDPAIRS.C
** 
** Contains code pertaining to the "independent pairs" of orbitals for which
** the energy is not invariant.  Only these pairs need to be considered in
** the MCSCF orbital rotation procedure.
**
** C. David Sherrill
** University of California, Berkeley
** May 1998
*/

extern "C" {
   #include <stdio.h>
   #include <stdlib.h>
}

#include "indpairs.h"

IndepPairs::IndepPairs() // Default constructor
{
  npairs = 0;
  nirreps = 0;
  p = NULL;
  q = NULL;
  map_pair_ir = NULL;
  map_pair_rel = NULL;
  ir_npairs = NULL;
  ir_p = NULL;
  ir_q = NULL;
  ir_p_rel = NULL;
  ir_q_rel = NULL;
  ir_map_pair = NULL;
}


// Regular constructor
IndepPairs::IndepPairs(int nirr, int num_ras, int **ras_opi, int ***ras_orbs,
                       int **fzc_orbs, int **fzv_orbs, int *frozen_docc, 
                       int *frozen_uocc, int *ci2relpitz, int fciflag) 
{

  set(nirr, num_ras, ras_opi, ras_orbs, fzc_orbs, fzv_orbs, 
      frozen_docc, frozen_uocc, ci2relpitz, fciflag);

}

void IndepPairs::set(int nirr, int num_ras, int **ras_opi, int ***ras_orbs,
                    int **fzc_orbs, int **fzv_orbs, int *frozen_docc, 
                    int *frozen_uocc, int *ci2relpitz, int fciflag) 
{
  int count, ir_count, *ir_cnt;
  int rasi, rasj, irrep;
  int i, j;

  nirreps = nirr;
  ir_npairs = new int[nirreps];

  // count everything up!

  for (irrep=0; irrep<nirreps; irrep++) ir_npairs[irrep] = 0;

  for (rasj=0; rasj<num_ras; rasj++)
    for (irrep=0; irrep<nirreps; irrep++) 
      ir_npairs[irrep] += ras_opi[rasj][irrep] * frozen_docc[irrep];

  if (!fciflag) {
    for (rasj=0; rasj<num_ras; rasj++)
      for (rasi=rasj+1; rasi<num_ras; rasi++)
        for (irrep=0; irrep<nirreps; irrep++)
          ir_npairs[irrep] += ras_opi[rasi][irrep] * ras_opi[rasj][irrep];
  }

  for (rasj=0; rasj<num_ras; rasj++)
    for (irrep=0; irrep<nirreps; irrep++)
      ir_npairs[irrep] += ras_opi[rasj][irrep] * frozen_uocc[irrep];

  for (irrep=0; irrep<nirreps; irrep++)
    ir_npairs[irrep] += frozen_docc[irrep] * frozen_uocc[irrep];

  npairs = 0;
  for (irrep=0; irrep<nirreps; irrep++) 
    npairs += ir_npairs[irrep];

  if (npairs==0)  {
    printf("(IndepPairs): Constructor called but no pairs!!\n");
    fprintf(stderr, "(IndepPairs): Constructor called but no pairs!!\n");
  }

  else {
    p = new int[npairs];
    q = new int[npairs];
    map_pair_ir  = new int[npairs];
    map_pair_rel = new int[npairs]; 

    ir_p = new int*[nirreps];
    ir_q = new int*[nirreps];
    ir_p_rel = new int*[nirreps];
    ir_q_rel = new int*[nirreps];
    ir_map_pair = new int*[nirreps];

    for (irrep=0; irrep<nirreps; irrep++) {
      i = ir_npairs[irrep];
      if (!i) continue;
      ir_p[irrep] = new int[i];
      ir_q[irrep] = new int[i];
      ir_p_rel[irrep] = new int[i];
      ir_q_rel[irrep] = new int[i];
      ir_map_pair[irrep] = new int[i];
    }

    ir_cnt = new int[nirreps];
    for (irrep=0; irrep<nirreps; irrep++) ir_cnt[irrep] = 0; 
    count = 0;

    // Now put everything in the proper arrays

    // Independent pairs for fzc with rasj 
    for (rasj=0; rasj<num_ras; rasj++) {
      for (irrep=0; irrep<nirreps; irrep++) {
        ir_count = ir_cnt[irrep];
        for (i=0; i<ras_opi[rasj][irrep]; i++) {
          for (j=0; j<frozen_docc[irrep]; j++) {
            p[count] = ras_orbs[rasj][irrep][i];
            q[count] = fzc_orbs[irrep][j];
            map_pair_ir[count] = irrep;
            map_pair_rel[count] = ir_count;
            ir_p[irrep][ir_count] = p[count];
            ir_q[irrep][ir_count] = q[count];
            ir_p_rel[irrep][ir_count] = ci2relpitz[p[count]];
            ir_q_rel[irrep][ir_count] = ci2relpitz[q[count]];
            ir_map_pair[irrep][ir_count] = count;
            count++;
            ir_count++;
          }
        }
      ir_cnt[irrep] = ir_count;
      }
    }


    // Independent pairs for rasi with rasj 
    if (!fciflag) {
      for (rasj=0; rasj<num_ras; rasj++) {
        for (rasi=rasj+1; rasi<num_ras; rasi++) {
          for (irrep=0; irrep<nirreps; irrep++) {
            ir_count = ir_cnt[irrep];
            for (j=0; j<ras_opi[rasj][irrep]; j++) {
              for (i=0; i<ras_opi[rasi][irrep]; i++) {
                p[count] = ras_orbs[rasi][irrep][i];
                q[count] = ras_orbs[rasj][irrep][j];
                map_pair_ir[count] = irrep;
                map_pair_rel[count] = ir_count;
                ir_p[irrep][ir_count] = p[count];
                ir_q[irrep][ir_count] = q[count];
                ir_p_rel[irrep][ir_count] = ci2relpitz[p[count]];
                ir_q_rel[irrep][ir_count] = ci2relpitz[q[count]];
                ir_map_pair[irrep][ir_count] = count;
                count++;
                ir_count++;
              }
            }
          ir_cnt[irrep] = ir_count;
          }
        }
      }
    }

    // Independent pairs for fzv with rasj 
    for (rasj=0; rasj<num_ras; rasj++) {
      for (irrep=0; irrep<nirreps; irrep++) {
        ir_count = ir_cnt[irrep];
        for (j=0; j<ras_opi[rasj][irrep]; j++) {
          for (i=0; i<frozen_uocc[irrep]; i++) {
            p[count] = fzv_orbs[irrep][i];
            q[count] = ras_orbs[rasj][irrep][j];
            map_pair_ir[count] = irrep;
            map_pair_rel[count] = ir_count;
            ir_p[irrep][ir_count] = p[count];
            ir_q[irrep][ir_count] = q[count];
            ir_p_rel[irrep][ir_count] = ci2relpitz[p[count]];
            ir_q_rel[irrep][ir_count] = ci2relpitz[q[count]];
            ir_map_pair[irrep][ir_count] = count;
            count++;
            ir_count++;
          }
        }
        ir_cnt[irrep] = ir_count;
      }
    }

    // Independent pairs for fzc with fzv 
    for (irrep=0; irrep<nirreps; irrep++) {
      ir_count = ir_cnt[irrep];
      for (i=0; i<frozen_uocc[irrep]; i++) {
        for (j=0; j<frozen_docc[irrep]; j++) {
          p[count] = fzv_orbs[irrep][i];
          q[count] = fzc_orbs[irrep][j];
          map_pair_ir[count] = irrep;
          map_pair_rel[count] = ir_count;
          ir_p[irrep][ir_count] = p[count];
          ir_q[irrep][ir_count] = q[count];
          ir_p_rel[irrep][ir_count] = ci2relpitz[p[count]];
          ir_q_rel[irrep][ir_count] = ci2relpitz[q[count]];
          ir_map_pair[irrep][ir_count] = count;
          count++;
          ir_count++;
        }
      }
      ir_cnt[irrep] = ir_count;
    }

    // check things
    if (count != npairs) {
      printf("(IndepPairs::set): mismatch in counted pairs!\n");
      fprintf(stderr, "(IndepPairs::set): mismatch in counted pairs!\n");
    }

  } 

}



IndepPairs::~IndepPairs() // Destructor
{
  int h;

  if (npairs) {
    delete [] p;
    delete [] q;
    delete [] map_pair_ir;
    delete [] map_pair_rel;
  }

  if (ir_npairs != NULL) {
    for (h=0; h<nirreps; h++) {
      if (ir_npairs[h]) {
        delete [] ir_p[h];
        delete [] ir_q[h];
        delete [] ir_p_rel[h];
        delete [] ir_q_rel[h];
        delete [] ir_map_pair[h];
      }  
    }
  }

  if (npairs) {
    delete [] ir_npairs;
    delete [] ir_p;
    delete [] ir_q;
    delete [] ir_p_rel;
    delete [] ir_q_rel;
    delete [] ir_map_pair;
  }

}


void IndepPairs::print(FILE *outfile)
{
  int h;

  fprintf(outfile, "\nList of all independent pairs:\n");
  print_selected(npairs, p, q, outfile);
  fprintf(outfile, "\nLists of independent pairs by irrep:\n");
  
  for (h=0; h<nirreps; h++) {
    if (!ir_npairs[h]) continue;
    fprintf(outfile, "\n\t Irrep %d:\n", h);
    print_selected(ir_npairs[h], ir_p[h], ir_q[h], 
                   ir_p_rel[h], ir_q_rel[h], outfile);
  }

}


void IndepPairs::print_selected(int num, int *parr, int *qarr, FILE *outfile)
{
  int ii;

  fprintf(outfile, "\n  %4d Independent Pairs\n", num);
  fprintf(outfile, "\t p\t q\n");
  fprintf(outfile,   "    -------------------\n");

  for (ii=0; ii<num; ii++) {
    fprintf(outfile,"\t %2d\t",parr[ii]);
    fprintf(outfile,"%2d\n",qarr[ii]);
  }
  fflush(outfile);

}

void IndepPairs::print_selected(int num, int *parr, int *qarr, 
                                int *prel, int *qrel, FILE *outfile)
{
  int ii;

  fprintf(outfile, "\n\t %4d Independent Pairs\n", num);
  fprintf(outfile, "\t p\t q\t P\t Q\n");
  fprintf(outfile,   "    ----------------------------------\n");

  for (ii=0; ii<num; ii++) {
    fprintf(outfile,"\t %2d\t%2d\t%2d\t%2d\n",
            parr[ii],qarr[ii],prel[ii],qrel[ii]);
  }
  fflush(outfile);

}


void IndepPairs::print_vec(double *arr, char *label, FILE *outfile) 
{
  int pair;

  fprintf(outfile, "%s\n", label);
  for (pair=0; pair<npairs; pair++) {
    fprintf(outfile, "Pair (%2d,%2d) = %12.6lf\n",
            p[pair], q[pair], arr[pair]);
  }
  fprintf(outfile, "\n");
  fflush(outfile);

}

/*
** get_irrep_vec()
**
** This function gets the piece of a vector belonging to a given irrep.
**  Assumes that the ordering of the input array is consistent with the
**  ordering of independent pairs in the total list of pairs, and has
**  the same length as the total number of independent pairs.
*/
double * IndepPairs::get_irrep_vec(int irrep, double *arr)
{

  int pair, target_len;
  double *newarr;

  target_len = ir_npairs[irrep];
  if (target_len==0)  return (NULL);

  newarr = new double[target_len];
  for (pair=0; pair<target_len; pair++) {
    newarr[pair] = arr[ir_map_pair[irrep][pair]];
  }

  return(newarr);

}


/*
** put_irrep_vec()
**
** This function scatters a vector for an irrep to the overall vector
**  using the mapping arrays in the independent pairs class.  This is
**  basically the reverse of get_irrep_vec().
*/
void IndepPairs::put_irrep_vec(int irrep, double *ir_vec, double *tot_vec)
{
  int pair;

  for (pair=0; pair<ir_npairs[irrep]; pair++) {
    tot_vec[ir_map_pair[irrep][pair]] = ir_vec[pair];
  }

} 



