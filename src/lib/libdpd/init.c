#include <stdio.h>
#include <stdarg.h>
#include <libciomr.h>
#include "dpd.h"
#include "dpd.gbl"

struct dpdpair {
    int *left_orbspi;
    int *right_orbspi;
    int *left_orbsym;
    int *right_orbsym;
    int *left_orboff;
    int *right_orboff;
    int permlr;
    int ler;
};


int dpd_init(int nirreps, int memory, int num_subspaces, ...)
{
  int h,h0,h1,cnt,***dp,l_irrep,r_irrep,p,q;
  int i,j,k,l,*count,offset1,offset2;
  int num_pairs, **orbspi, **orbsym, *numorbs, **orboff;
  int ***pairidx, **pairtot, ****pairorb;
  int **oe_orbidx, ***oe_orbs;
  struct dpdpair *pairs;
  va_list ap;

  dpd_nirreps = nirreps;
  dpd_memory = memory;
  dpd_num_subspaces = num_subspaces;

  /* Construct binary direct product array */
  dp = (int ***) malloc(nirreps * sizeof(int **));
  for(h=0; h < nirreps; h++) {
      dp[h] = init_int_matrix(nirreps,2);
      cnt=0;
      for(h0=0; h0 < nirreps; h0++) {
          for(h1=0; h1 < nirreps; h1++) {
              if((h0^h1)==h) {
                  dp[h][cnt][0] = h0;
                  dp[h][cnt++][1] = h1;
                }
            }
        }
    }

  /* Grab the irrep population and orbital symmetry arrays from the arg list */
  va_start(ap, num_subspaces);
  orbspi = (int **) malloc(sizeof(int *) * num_subspaces);
  orbsym = (int **) malloc(sizeof(int *) * num_subspaces);
  for(i=0; i < num_subspaces; i++) {
      orbspi[i] = va_arg(ap, int *);
      orbsym[i] = va_arg(ap, int *);
    }
  va_end(ap);
  dpd_orbspi = orbspi;
  dpd_orbsym = orbsym;

  /* Compute the number of orbitals in each subspace */
  numorbs = (int *) malloc(num_subspaces * sizeof(int));
  for(i=0; i < num_subspaces; i++) {
      numorbs[i] = 0;
      for(h=0; h < nirreps; h++)
	  numorbs[i] += orbspi[i][h];
    }
  dpd_numorbs = numorbs;

  /* Compute the orbital offset arrays */
  orboff = (int **) malloc(num_subspaces * sizeof(int *));
  for(i=0; i < num_subspaces; i++) {
      orboff[i] = init_int_array(nirreps);
      for(j=1; j < nirreps; j++)
	  orboff[i][j] = orboff[i][j-1] + orbspi[i][j-1];
    }
  dpd_orboff = orboff;

  /* Compute the number of bra or ket index combinations */
  num_pairs = (num_subspaces * (num_subspaces - 1)) + (5 * num_subspaces);
  dpd_num_pairs = num_pairs;

  /* Set up the pair structs for later use */
  pairs = (struct dpdpair *) malloc(num_pairs * sizeof(struct dpdpair));

  /* Build the row/column dimension arrays */
  pairtot = init_int_matrix(num_pairs, nirreps);
  dpd_pairtot = pairtot;
  /* Loop over the groups of three "diagonal" pairs */
  for(i=0; i < num_subspaces; i++) {

      pairs[5*i].left_orbspi = orbspi[i];
      pairs[5*i].left_orbsym = orbsym[i];
      pairs[5*i].left_orboff = orboff[i];
      pairs[5*i].right_orbspi = orbspi[i];
      pairs[5*i].right_orbsym = orbsym[i];
      pairs[5*i].right_orboff = orboff[i];
      pairs[5*i].permlr = 0;
      pairs[5*i].ler = 0;

      pairs[5*i+1].left_orbspi = orbspi[i];
      pairs[5*i+1].left_orbsym = orbsym[i];
      pairs[5*i+1].left_orboff = orboff[i];
      pairs[5*i+1].right_orbspi = orbspi[i];
      pairs[5*i+1].right_orbsym = orbsym[i];
      pairs[5*i+1].right_orboff = orboff[i];
      pairs[5*i+1].permlr = 1;
      pairs[5*i+1].ler = 0;

      pairs[5*i+2].left_orbspi = orbspi[i];
      pairs[5*i+2].left_orbsym = orbsym[i];
      pairs[5*i+2].left_orboff = orboff[i];
      pairs[5*i+2].right_orbspi = orbspi[i];
      pairs[5*i+2].right_orbsym = orbsym[i];
      pairs[5*i+2].right_orboff = orboff[i];
      pairs[5*i+2].permlr = -1;
      pairs[5*i+2].ler = 0;

      pairs[5*i+3].left_orbspi = orbspi[i];
      pairs[5*i+3].left_orbsym = orbsym[i];
      pairs[5*i+3].left_orboff = orboff[i];
      pairs[5*i+3].right_orbspi = orbspi[i];
      pairs[5*i+3].right_orbsym = orbsym[i];
      pairs[5*i+3].right_orboff = orboff[i];
      pairs[5*i+3].permlr = 1;
      pairs[5*i+3].ler = 1;

      pairs[5*i+4].left_orbspi = orbspi[i];
      pairs[5*i+4].left_orbsym = orbsym[i];
      pairs[5*i+4].left_orboff = orboff[i];
      pairs[5*i+4].right_orbspi = orbspi[i];
      pairs[5*i+4].right_orbsym = orbsym[i];
      pairs[5*i+4].right_orboff = orboff[i];
      pairs[5*i+4].permlr = -1;
      pairs[5*i+4].ler = 1;
  		  
      for(j=0; j < nirreps; j++) 
	  for(k=0; k < nirreps; k++) {
	      l_irrep = dp[j][k][0]; r_irrep = dp[j][k][1];

	      /* orbspi,orbspi */
	      pairtot[5*i][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];

	      if(l_irrep > r_irrep) {
		  /* orbspi < orbspi, +1 */
		  pairtot[5*i+1][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];
		  /* orbspi < orbspi, -1 */
		  pairtot[5*i+2][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];
		  /* orbspi <= orbspi, +1 */
		  pairtot[5*i+3][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];
		  /* orbspi <= orbspi, -1 */
		  pairtot[5*i+4][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];
		}
	      else if(l_irrep == r_irrep) {
		  /* orbspi < orbspi, +1 */
		  pairtot[5*i+1][j] +=
		     (orbspi[i][l_irrep] * (orbspi[i][l_irrep]-1))/2;
		  /* orbspi < orbspi, -1 */
		  pairtot[5*i+2][j] +=
		     (orbspi[i][l_irrep] * (orbspi[i][l_irrep]-1))/2;
		  /* orbspi <= orbspi, +1 */
		  pairtot[5*i+3][j] +=
		     (orbspi[i][l_irrep] * (orbspi[i][l_irrep]+1))/2;
		  /* orbspi <= orbspi, -1 */
		  pairtot[5*i+4][j] +=
		     (orbspi[i][l_irrep] * (orbspi[i][l_irrep]+1))/2;
		}
	    }
    }
	  
   /* Loop over the remaining "off diagonal" pairs */
  for(i=0,cnt=5*num_subspaces; i < num_subspaces; i++) 
      for(j=i+1; j < num_subspaces; j++,cnt+=2) {

	  pairs[cnt].left_orbspi = orbspi[i];
	  pairs[cnt].left_orbsym = orbsym[i];
	  pairs[cnt].left_orboff = orboff[i];
	  pairs[cnt].right_orbspi = orbspi[j];
	  pairs[cnt].right_orbsym = orbsym[j];
	  pairs[cnt].right_orboff = orboff[j];
	  pairs[cnt].permlr = 0;
	  pairs[cnt].ler = 0;

	  pairs[cnt+1].left_orbspi = orbspi[j];
	  pairs[cnt+1].left_orbsym = orbsym[j];
	  pairs[cnt+1].left_orboff = orboff[j];
	  pairs[cnt+1].right_orbspi = orbspi[i];
	  pairs[cnt+1].right_orbsym = orbsym[i];
	  pairs[cnt+1].right_orboff = orboff[i];
	  pairs[cnt+1].permlr = 0;
	  pairs[cnt+1].ler = 0;
  	  
	  for(k=0; k < nirreps; k++)
	      for(l=0; l < nirreps; l++) {
	      
		  l_irrep = dp[k][l][0]; r_irrep = dp[k][l][1];

		  /* orbspi[i],orbspi[j] */
		  pairtot[cnt][k] += orbspi[i][l_irrep] * orbspi[j][r_irrep];
		  /* orbspi[j],orbspi[i] */
		  pairtot[cnt+1][k] += orbspi[j][l_irrep] * orbspi[i][r_irrep];

		}
	}

  /* Temporary check until I'm sure I'm doing this right */
  if(num_pairs != cnt) { printf("Error in dpd_init()!\n"); exit(2); }

  /* Build the row/column index lookup arrays */
  pairidx = (int ***) malloc(num_pairs * sizeof(int **));
  pairorb = (int ****) malloc(num_pairs * sizeof(int ***));
  dpd_pairidx = pairidx;
  dpd_pairorb = pairorb;
  count = init_int_array(nirreps);
  /* Loop over the groups of three "diagonal" pairs */
  for(i=0; i < num_subspaces; i++) {

      for(l=0; l < 5; l++) {
	  pairidx[5*i+l] = init_int_matrix(numorbs[i],numorbs[i]);
	  for(j=0; j < numorbs[i]; j++)
	      for(k=0; k < numorbs[i]; k++) 
		  pairidx[5*i+l][j][k] = -1;
      
	  pairorb[5*i+l] = (int ***) malloc(nirreps * sizeof(int **));
	  for(j=0; j < nirreps; j++) {
	      pairorb[5*i+l][j] =
	       pairtot[5*i+l][j] ? init_int_matrix(pairtot[5*i+l][j],2) : NULL;
	      for(k=0; k < pairtot[5*i+l][j]; k++) {
		  pairorb[5*i+l][j][k][0] = -1;
		  pairorb[5*i+l][j][k][1] = -1;
		}
	    }
	}

      zero_int_array(count,nirreps);

      /* orbspi[i],orbspi[i] */
      for(j=0; j < nirreps; j++)
	  for(k=0; k < nirreps; k++) {
	      h0 = dp[j][k][0]; h1 = dp[j][k][1];
	      offset1 = orboff[i][h0];  offset2 = orboff[i][h1];
	      for(p=0; p < orbspi[i][h0]; p++)
		  for(q=0; q < orbspi[i][h1]; q++) {
		      pairorb[5*i][j][count[j]][0] = p+offset1;
		      pairorb[5*i][j][count[j]][1] = q+offset2;
		      pairidx[5*i][p+offset1][q+offset2] = count[j]++;
		    }
	    }

      zero_int_array(count, nirreps);

      /* orbspi[i] < orbspi[i], +1 */
      for(j=0; j < nirreps; j++)
	  for(k=0; k < nirreps; k++) {
	      h0 = dp[j][k][0]; h1 = dp[j][k][1];
	      offset1 = orboff[i][h0]; offset2 = orboff[i][h1];
	      if(h0 == h1) {
		  for(p=0; p < orbspi[i][h0]; p++)
		      for(q=0; q < p; q++) {
			  pairorb[5*i+1][j][count[j]][0] = p+offset1;
			  pairorb[5*i+1][j][count[j]][1] = q+offset2;
			  pairidx[5*i+1][p+offset1][q+offset2] = count[j];
			  pairidx[5*i+1][q+offset2][p+offset1] = count[j]++;
			}
		}
	      else if(h0 > h1) {
		  for(p=0; p < orbspi[i][h0]; p++)
		      for(q=0; q < orbspi[i][h1]; q++) {
			  pairorb[5*i+1][j][count[j]][0] = p+offset1;
			  pairorb[5*i+1][j][count[j]][1] = q+offset2;
			  pairidx[5*i+1][p+offset1][q+offset2] = count[j];
			  pairidx[5*i+1][q+offset2][p+offset1] = count[j]++;
			}
		}
	    }

      zero_int_array(count, nirreps);

      /* orbspi[i] < orbspi[i], -1 */
      for(j=0; j < nirreps; j++)
	  for(k=0; k < nirreps; k++) {
	      h0 = dp[j][k][0]; h1 = dp[j][k][1];
	      offset1 = orboff[i][h0]; offset2 = orboff[i][h1];
	      if(h0 == h1) {
		  for(p=0; p < orbspi[i][h0]; p++)
		      for(q=0; q < p; q++) {
			  pairorb[5*i+2][j][count[j]][0] = p+offset1;
			  pairorb[5*i+2][j][count[j]][1] = q+offset2;
			  pairidx[5*i+2][p+offset1][q+offset2] = count[j];
			  pairidx[5*i+2][q+offset2][p+offset1] = count[j]++;
			}
		}
	      else if(h0 > h1) {
		  for(p=0; p < orbspi[i][h0]; p++)
		      for(q=0; q < orbspi[i][h1]; q++) {
			  pairorb[5*i+2][j][count[j]][0] = p+offset1;
			  pairorb[5*i+2][j][count[j]][1] = q+offset2;
			  pairidx[5*i+2][p+offset1][q+offset2] = count[j];
			  pairidx[5*i+2][q+offset2][p+offset1] = count[j]++;
			}
		}
	    }

      zero_int_array(count, nirreps);

      /* orbspi[i] <= orbspi[i], +1 */
      for(j=0; j < nirreps; j++)
	  for(k=0; k < nirreps; k++) {
	      h0 = dp[j][k][0]; h1 = dp[j][k][1];
	      offset1 = orboff[i][h0]; offset2 = orboff[i][h1];
	      if(h0 == h1) {
		  for(p=0; p < orbspi[i][h0]; p++)
		      for(q=0; q <= p; q++) {
			  pairorb[5*i+3][j][count[j]][0] = p+offset1;
			  pairorb[5*i+3][j][count[j]][1] = q+offset2;
			  pairidx[5*i+3][p+offset1][q+offset2] = count[j];
			  pairidx[5*i+3][q+offset2][p+offset1] = count[j]++;
			}
		}
	      else if(h0 > h1) {
		  for(p=0; p < orbspi[i][h0]; p++)
		      for(q=0; q < orbspi[i][h1]; q++) {
			  pairorb[5*i+3][j][count[j]][0] = p+offset1;
			  pairorb[5*i+3][j][count[j]][1] = q+offset2;
			  pairidx[5*i+3][p+offset1][q+offset2] = count[j];
			  pairidx[5*i+3][q+offset2][p+offset1] = count[j]++;
			}
		}
	    }

      zero_int_array(count, nirreps);

      /* orbspi[i] <= orbspi[i], -1 */
      for(j=0; j < nirreps; j++)
	  for(k=0; k < nirreps; k++) {
	      h0 = dp[j][k][0]; h1 = dp[j][k][1];
	      offset1 = orboff[i][h0]; offset2 = orboff[i][h1];
	      if(h0 == h1) {
		  for(p=0; p < orbspi[i][h0]; p++)
		      for(q=0; q <= p; q++) {
			  pairorb[5*i+4][j][count[j]][0] = p+offset1;
			  pairorb[5*i+4][j][count[j]][1] = q+offset2;
			  pairidx[5*i+4][p+offset1][q+offset2] = count[j];
			  pairidx[5*i+4][q+offset2][p+offset1] = count[j]++;
			}
		}
	      else if(h0 > h1) {
		  for(p=0; p < orbspi[i][h0]; p++)
		      for(q=0; q < orbspi[i][h1]; q++) {
			  pairorb[5*i+4][j][count[j]][0] = p+offset1;
			  pairorb[5*i+4][j][count[j]][1] = q+offset2;
			  pairidx[5*i+4][p+offset1][q+offset2] = count[j];
			  pairidx[5*i+4][q+offset2][p+offset1] = count[j]++;
			}
		}
	    }

    }

  /* Loop over the remaining "off diagonal" pairs */
  for(i=0,cnt=5*num_subspaces; i < num_subspaces; i++) {
      for(j=i+1; j < num_subspaces; j++,cnt+=2) {

	  pairidx[cnt] = init_int_matrix(numorbs[i],numorbs[j]);
	  for(k=0; k < numorbs[i]; k++)
	      for(l=0; l < numorbs[j]; l++)
		  pairidx[cnt][k][l] = -1;
	  pairidx[cnt+1] = init_int_matrix(numorbs[j],numorbs[i]);
	  for(k=0; k < numorbs[j]; k++)
	      for(l=0; l < numorbs[i]; l++)
		  pairidx[cnt+1][k][l] = -1;

	  pairorb[cnt] = (int ***) malloc(nirreps * sizeof(int **));
	  pairorb[cnt+1] = (int ***) malloc(nirreps * sizeof(int **));
	  for(k=0; k < nirreps; k++) {
	      pairorb[cnt][k] =
               pairtot[cnt][k] ? init_int_matrix(pairtot[cnt][k],2) : NULL;
	      pairorb[cnt+1][k] =
               pairtot[cnt+1][k] ? init_int_matrix(pairtot[cnt+1][k],2) : NULL;
	      for(l=0; l < pairtot[cnt][k]; l++) {
		  pairorb[cnt][k][l][0] = -1;
		  pairorb[cnt][k][l][1] = -1;
		}
	      for(l=0; l < pairtot[cnt+1][k]; l++) {
		  pairorb[cnt+1][k][l][0] = -1;
		  pairorb[cnt+1][k][l][1] = -1;
		}

	    }

	  zero_int_array(count, nirreps);
	  
	  for(k=0; k < nirreps; k++)
	      for(l=0; l < nirreps; l++) {
		  h0 = dp[k][l][0]; h1 = dp[k][l][1];
		  offset1 = orboff[i][h0];  offset2 = orboff[j][h1];
		  for(p=0; p < orbspi[i][h0]; p++)
		      for(q=0; q < orbspi[j][h1]; q++) {
			  pairorb[cnt][k][count[k]][0] = p+offset1;
			  pairorb[cnt][k][count[k]][1] = q+offset2;
			  pairidx[cnt][p+offset1][q+offset2] = count[k]++;
			}
		}

	  zero_int_array(count, nirreps);

	  for(k=0; k < nirreps; k++)
	      for(l=0; l < nirreps; l++) {
		  h0 = dp[k][l][0]; h1 = dp[k][l][1];
		  offset1 = orboff[j][h0];  offset2 = orboff[i][h1];
		  for(p=0; p < orbspi[j][h0]; p++)
		      for(q=0; q < orbspi[i][h1]; q++) {
			  pairorb[cnt+1][k][count[k]][0] = p+offset1;
			  pairorb[cnt+1][k][count[k]][1] = q+offset2;
			  pairidx[cnt+1][p+offset1][q+offset2] = count[k]++;
			}
		}
	}
    }

  /* Temporary check until I'm sure I'm doing this right */
  if(num_pairs != cnt) { printf("Error in dpd_init()!\n"); exit(3); }

  for(h=0; h < nirreps; h++) free_int_matrix(dp[h], nirreps);
  free(dp);

  /* Now generate the global list of DPD parameters */
  dpd_params =
	  (struct dpdparams **) malloc(num_pairs*sizeof(struct dpdparams *));
  for(i=0; i < num_pairs; i++)
      dpd_params[i] =
	  (struct dpdparams *) malloc(num_pairs*sizeof(struct dpdparams));

  for(i=0; i < num_pairs; i++) {
      for(j=0; j < num_pairs; j++) {
	  dpd_params[i][j].nirreps = nirreps;

	  dpd_params[i][j].pqnum = i;
	  dpd_params[i][j].rsnum = j;

	  dpd_params[i][j].rowtot = pairtot[i];
	  dpd_params[i][j].coltot = pairtot[j];

	  dpd_params[i][j].rowidx = pairidx[i];
	  dpd_params[i][j].colidx = pairidx[j];

	  dpd_params[i][j].roworb = pairorb[i];
	  dpd_params[i][j].colorb = pairorb[j];

	  dpd_params[i][j].ppi = pairs[i].left_orbspi;
	  dpd_params[i][j].qpi = pairs[i].right_orbspi;

	  dpd_params[i][j].rpi = pairs[j].left_orbspi;
	  dpd_params[i][j].spi = pairs[j].right_orbspi;

	  dpd_params[i][j].psym = pairs[i].left_orbsym;
	  dpd_params[i][j].qsym = pairs[i].right_orbsym;

	  dpd_params[i][j].rsym = pairs[j].left_orbsym;
	  dpd_params[i][j].ssym = pairs[j].right_orbsym;

	  dpd_params[i][j].poff = pairs[i].left_orboff;
	  dpd_params[i][j].qoff = pairs[i].right_orboff;

	  dpd_params[i][j].roff = pairs[j].left_orboff;
	  dpd_params[i][j].soff = pairs[j].right_orboff;

	  dpd_params[i][j].perm_pq = pairs[i].permlr;
	  dpd_params[i][j].perm_rs = pairs[j].permlr;
	  dpd_params[i][j].peq = pairs[i].ler;
	  dpd_params[i][j].res = pairs[j].ler;
	}
    }

  /* Now generate the global list of one-electron DPD parameters */
  oe_orbidx = (int **) malloc(num_subspaces*sizeof(int *));
  for(i=0; i < num_subspaces; i++) {
      oe_orbidx[i] = init_int_array(numorbs[i]);
      for(j=0; j < numorbs[i]; j++)
	  oe_orbidx[i][j] = -1;
    }
  dpd_oe_orbidx = oe_orbidx;

  oe_orbs = (int ***) malloc(num_subspaces*sizeof(int **));
  for(i=0;i < num_subspaces; i++) {
      oe_orbs[i] = (int **) malloc(nirreps*sizeof(int *));
      for(j=0; j < nirreps; j++) {
	  oe_orbs[i][j] = orbspi[i][j] ? init_int_array(orbspi[i][j]) : NULL;
	  for(k=0; k < orbspi[i][j]; k++)
	      oe_orbs[i][j][k] = -1;
	}
    }
  dpd_oe_orbs = oe_orbs;

  for(i=0; i < num_subspaces; i++) {
      zero_int_array(count, nirreps);
      for(j=0; j < nirreps; j++) {
	  offset1 = orboff[i][j];
	  for(p=0; p < orbspi[i][j]; p++) {
	      oe_orbs[i][j][count[j]] = p+offset1;
	      oe_orbidx[i][p+offset1] = count[j]++;
	    }
	}
    }

  dpd_oe_params =
   (struct oe_dpdparams **) malloc(num_subspaces*sizeof(struct oe_dpdparams *));
  for(i=0; i < num_subspaces; i++)
      dpd_oe_params[i] =
       (struct oe_dpdparams *) malloc(num_subspaces*sizeof(struct oe_dpdparams));

  for(i=0,cnt=0; i < num_subspaces; i++) {
      for(j=0; j < num_subspaces; j++,cnt++) {
	  dpd_oe_params[i][j].nirreps = nirreps;

	  dpd_oe_params[i][j].pnum = i;
	  dpd_oe_params[i][j].qnum = j;

	  dpd_oe_params[i][j].rowtot = orbspi[i];
	  dpd_oe_params[i][j].coltot = orbspi[j];

	  dpd_oe_params[i][j].rowidx = oe_orbidx[i];
	  dpd_oe_params[i][j].colidx = oe_orbidx[j];

	  dpd_oe_params[i][j].roworb = oe_orbs[i];
	  dpd_oe_params[i][j].colorb = oe_orbs[j];

	  dpd_oe_params[i][j].ppi = orbspi[i];
	  dpd_oe_params[i][j].qpi = orbspi[j];

	  dpd_oe_params[i][j].poff = orboff[i];
	  dpd_oe_params[i][j].qoff = orboff[j];

	  dpd_oe_params[i][j].psym = orbsym[i];
	  dpd_oe_params[i][j].qsym = orbsym[j];
	}
    }

  free(count);
  
  free(pairs);

  return 0;
}

