/******************************************************************************

	IRREP.CC

	finds reduction coefficients of a delocalized internal coordinate

        arguments:
          &simples     -- address of an object of class internals
          **di_coord   -- pointer to delocalized internal corrdinate matrix 
		          (eigenvectors of G)

        returns:
          -- pointer to array of pointers, each of which points to a 
             properly symmetrized delocalized internal coordinate (salc)
             with no linear dependence on any other salc

	note: also determines irrep of salc and places correponding number in
              'irr' array
	                                          significant modifications by
						  J. Kenny June '00
******************************************************************************/

extern "C" {
  #include <stdio.h>
  #include <file30.h>
  #include <stdlib.h>
  #include <string.h>
  #include <physconst.h>
  #include <math.h>
  #include <libciomr.h>
  #include <ip_libv1.h>
}

#define EXTERN
#include "opt.h"
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

double **irrep(internals &simples, double **di_coord) {
  int i, j, a, ops, coord_num,                       /* counter variables */
      **ct_tors,
      **ict,
      order, num_irreps = syminfo.num_irreps;        /* order of and number of irreps in molecular point group */
  double sum,
         *tmp_evect,                                 /* stores effect of operation on salc */
         *evect_character,
         char_irrep,
         **symm_coord,
         **irreps_spanned;                          /* each row corresponds to one coordinate, 1 in position corresponding
						     to irrep (first entry first irrep etc.) indicates that irrep is
						     spanned, zero otherwise*/
  char buffer[MAX_LINELENGTH], *err;
  int h,
      *spanned_arr,                                    /* holds numbers irreps spanned by all salcs */
      num_spanned;                                    /* number of irreps spanned by one salc */
     
  
  /*** allocate some memory ***/
  /* note: for 'symm_coord' since we can just swap pointers to coordinates, we only need enough memory to hold the pointers */
  symm_coord = (double **) malloc(num_nonzero*sizeof(double *));
  evect_character = (double *) malloc(num_irreps*sizeof(double));
  tmp_evect = (double *) malloc(simples.get_num()*sizeof(double));
  irreps_spanned = init_matrix(num_nonzero,num_irreps);
  spanned_arr = init_int_array(num_nonzero);
  
  /*determine order of point group*/
  order = 0;
  for(i=0;i<num_irreps;i++)
    order += ops_in_class[i]; 
    
  
  /*** loop over number of salcs ***/
  for(coord_num=0;coord_num<num_nonzero;coord_num++) { 
      
    num_spanned = 0;
    zero_arr(evect_character,num_irreps);

    /* Determine the characters of the salc under different operations */
    for (ops = 0;ops<num_irreps;++ops) {
      zero_arr(tmp_evect, simples.get_num());
      for (i=0;i<simples.get_num();++i) {
        a = simples.id_to_index(syminfo.ict_ops[i][ops]);
        tmp_evect[a] += syminfo.ict_ops_sign[i][ops] * di_coord[coord_num][i];
      }
      for (i=0;i<simples.get_num();++i) {
             evect_character[ops] += tmp_evect[i]*di_coord[coord_num][i];
       }
    }
  
    /* Determine the irrep of the salc */
    for (i=0;i<num_irreps;++i) {
      char_irrep=0.0;
      for (j=0;j<ops;++j) {
        char_irrep += evect_character[j] * syminfo.ct[i][j] * ops_in_class[j];
      }
      char_irrep = char_irrep/ (float)order;
      if( char_irrep > 0.1) {
        ++num_spanned;
	irreps_spanned[coord_num][i] = char_irrep;
      }
    }
    spanned_arr[coord_num] = num_spanned;
  }


  
  /*** now the fun stuff ***
    need to project coordinates that span more than one irrep onto those irreps and put everything together
    without creating linear dependencies in the coordinates,  this should be implemented soon.  in the
    meantime molecules with eigenvectors that span only one irrep (I use this loosely for now) should
    optimize fine*/

  /* first count number of irreps spanned in total (reuse 'num_spanned') */ 
  num_spanned = 0;
  for(coord_num=0;coord_num<num_nonzero;coord_num++) {
      num_spanned += spanned_arr[coord_num];
    }
  
  /* if each coordinate spans only one irrep we won't bother with projection operator,
     this saves a little effort with simple cases (can't right now anyways) */
  if(num_spanned == num_nonzero) {
    for(coord_num=0;coord_num<num_nonzero;coord_num++) { 
	for(i=0;i<num_irreps;++i) {
          if(irreps_spanned[coord_num][i] > 0.1) irr[coord_num] = i; }
      symm_coord[coord_num] =  di_coord[coord_num];
    }
  }

  /* if more than one irrep spanned per coordinate, big problems */
  else if(num_spanned > num_nonzero) {
      fprintf(outfile,"\n-- delocalized internal coordinates span more than one irrep ... ");
      fprintf(outfile,"\n      ... currently unable to handle this situation, sorry");
      fprintf(outfile,"\n-- stopping execution \n\n");
      fflush(outfile);
      exit(1);
     
    }

return symm_coord;
}






























