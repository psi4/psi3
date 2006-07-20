#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <ccfiles.h>
#include <psifiles.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void idx_permute_presort(dpdfile4 *File, int this_bucket,
			   int **bucket_map, int **bucket_offset,
			   int p, int q, int r, int s,
			   double value, FILE *outfile);

int file_build_presort(dpdfile4 *File, int inputfile, double tolerance, 
		       long int memoryb, int keep, int fzc, 
		       double **Da, double **Db, double **fock_a, double **fock_b, int ref)
{
  struct iwlbuf InBuf;
  int lastbuf;
  long int memoryd, core_left, row_length;
  int h, nirreps, n, row, nump, numq, nbuckets;
  int **bucket_map, **bucket_offset, **bucket_rowdim;
  long int **bucket_size;
  Value *valptr;
  Label *lblptr;
  int idx, p, q, r, s, pq, rs;
  double value;
  struct iwlbuf *SortBuf;
  psio_address next;
  double **D, **fock;

  nirreps = File->params->nirreps;

  memoryd = memoryb/sizeof(double);

  /* It's annoying that I have to compute this here */
  for(h=0,nump=0,numq=0; h < nirreps; h++) {
    nump += File->params->ppi[h];
    numq += File->params->qpi[h];
  }
  bucket_map = init_int_matrix(nump,numq);

  /* Room for one bucket to begin with */
  bucket_offset = (int **) malloc(sizeof(int *));
  bucket_offset[0] = init_int_array(nirreps);
  bucket_rowdim = (int **) malloc(sizeof(int *));
  bucket_rowdim[0] = init_int_array(nirreps);
  bucket_size = (long int **) malloc(sizeof(long int *));
  bucket_size[0] = init_long_int_array(nirreps);
    
  /* Figure out how many passes we need and where each p,q goes */
  for(h=0,core_left=memoryd,nbuckets=1; h < nirreps; h++) {

    row_length = (long int) File->params->coltot[h^(File->my_irrep)];
	       
    for(row=0; row < File->params->rowtot[h]; row++) {

      if((core_left - row_length) >= 0) {
	core_left -= row_length;
	bucket_rowdim[nbuckets-1][h]++;
	bucket_size[nbuckets-1][h] += row_length;
      }
      else {
	nbuckets++;
	core_left = memoryd - row_length;

	/* Make room for another bucket */
	bucket_offset = (int **) realloc((void *) bucket_offset,
					 nbuckets * sizeof(int *));
	bucket_offset[nbuckets-1] = init_int_array(nirreps);
	bucket_offset[nbuckets-1][h] = row;

	bucket_rowdim = (int **) realloc((void *) bucket_rowdim,
					 nbuckets * sizeof(int *));
	bucket_rowdim[nbuckets-1] = init_int_array(nirreps);
	bucket_rowdim[nbuckets-1][h] = 1;

	bucket_size = (long int **) realloc((void *) bucket_size,
					    nbuckets * sizeof(long int *));
	bucket_size[nbuckets-1] = init_long_int_array(nirreps);
	bucket_size[nbuckets-1][h] = row_length;
      }

      p = File->params->roworb[h][row][0];
      q = File->params->roworb[h][row][1];
      bucket_map[p][q] = nbuckets - 1;
    }
  }

  if(fzc && (ref==0 || ref==1)) { /* for convenience in RHF/ROHF cases */
    D = Da; 
    fock = fock_a;
  }

  fprintf(outfile, "\tSorting File: %s nbuckets = %d\n", File->label, nbuckets);
  fflush(outfile);

  next = PSIO_ZERO;
  for(n=0; n < nbuckets; n++) { /* nbuckets = number of passes */

    /* Prepare target matrix */
    for(h=0; h < nirreps; h++) {
      File->matrix[h] = block_matrix(bucket_rowdim[n][h], File->params->coltot[h]);
    }

    iwl_buf_init(&InBuf, inputfile, tolerance, 1, 1);

    lblptr = InBuf.labels;
    valptr = InBuf.values;
    lastbuf = InBuf.lastbuf;

    for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = abs((int) lblptr[idx++]);
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      value = (double) valptr[InBuf.idx];

      idx_permute_presort(File,n,bucket_map,bucket_offset,p,q,r,s,value,outfile);

      if(fzc && !n) { /* build frozen-core operator only on first pass*/
	if(ref==0 || ref==1) { /* RHF/ROHF */
	  pq = INDEX(p,q);
	  rs = INDEX(r,s);
	  /* (pq|rs) */
	  fock[p][q] += 2.0 * D[r][s] * value;
	  fock[p][r] -= D[q][s] * value;
	  if(pq!=rs && p!=q && r!=s) {
	    /* (pq|sr) */
	    fock[p][q] += 2.0 * D[s][r] * value;
	    fock[p][s] -= D[q][r] * value;
	    /* (qp|rs) */
	    fock[q][p] += 2.0 * D[r][s] * value;
	    fock[q][r] -= D[p][s] * value;
	    /* (qp|sr) */
	    fock[q][p] += 2.0 * D[s][r] * value;
	    fock[q][s] -= D[p][r] * value;
	    /* (rs|pq) */
	    fock[r][s] += 2.0 * D[p][q] * value;
	    fock[r][p] -= D[s][q] * value;
	    /* (rs|qp) */
	    fock[r][s] += 2.0 * D[q][p] * value;
	    fock[r][q] -= D[s][p] * value;
	    /* (sr|pq) */
	    fock[s][r] += 2.0 * D[p][q] * value;
	    fock[s][p] -= D[r][q] * value;
	    /* (sr|qp) */
	    fock[s][r] += 2.0 * D[q][p] * value;
	    fock[s][q] -= D[r][p] * value;
	  }
	  else if(p!=q && r!=s && pq==rs) {
	    /* (pq|sr) */
	    fock[p][q] += 2.0 * D[s][r] * value;
	    fock[p][s] -= D[q][r] * value;
	    /* (qp|rs) */
	    fock[q][p] += 2.0 * D[r][s] * value;
	    fock[q][r] -= D[p][s] * value;
	    /* (qp|sr) */
	    fock[q][p] += 2.0 * D[s][r] * value;
	    fock[q][s] -= D[p][r] * value;
	  }
	  else if(p!=q && r==s) {
	    /* (qp|rs) */
	    fock[q][p] += 2.0 * D[r][s] * value;
	    fock[q][r] -= D[p][s] * value;
	    /* (rs|pq) */
	    fock[r][s] += 2.0 * D[p][q] * value;
	    fock[r][p] -= D[s][q] * value;
	    /* (rs|qp) */
	    fock[r][s] += 2.0 * D[q][p] * value;
	    fock[r][q] -= D[s][p] * value;
	  }
	  else if(p==q && r!=s) {
	    /* (pq|sr) */
	    fock[p][q] += 2.0 * D[s][r] * value;
	    fock[p][s] -= D[q][r] * value;
	    /* (rs|pq) */
	    fock[r][s] += 2.0 * D[p][q] * value;
	    fock[r][p] -= D[s][q] * value;
	    /* (sr|pq) */
	    fock[s][r] += 2.0 * D[p][q] * value;
	    fock[s][p] -= D[r][q] * value;
	  }
	  else if(p==q && r==s && pq!=rs) {
	    /* (rs|pq) */
	    fock[r][s] += 2.0 * D[p][q] * value;
	    fock[r][p] -= D[s][q] * value;
	  }
	} /* if(ref==0||ref==1) */
      } /* if(fzc) */

    } /* end loop through current buffer */

    /* Now run through the rest of the buffers in the file */
    while (!lastbuf) {
      iwl_buf_fetch(&InBuf);
      lastbuf = InBuf.lastbuf;

      for (idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
	p = abs((int) lblptr[idx++]);
	q = (int) lblptr[idx++];
	r = (int) lblptr[idx++];
	s = (int) lblptr[idx++];

	value = (double) valptr[InBuf.idx];

	idx_permute_presort(File,n,bucket_map,bucket_offset,p,q,r,s,value,outfile);
      
	if(fzc && !n) { /* build frozen-core operator only on first pass */
	  if(ref==0 || ref==1) { /* RHF/ROHF */
	    pq = INDEX(p,q);
	    rs = INDEX(r,s);
	    /* (pq|rs) */
	    fock[p][q] += 2.0 * D[r][s] * value;
	    fock[p][r] -= D[q][s] * value;
	    if(pq!=rs && p!=q && r!=s) {
	      /* (pq|sr) */
	      fock[p][q] += 2.0 * D[s][r] * value;
	      fock[p][s] -= D[q][r] * value;
	      /* (qp|rs) */
	      fock[q][p] += 2.0 * D[r][s] * value;
	      fock[q][r] -= D[p][s] * value;
	      /* (qp|sr) */
	      fock[q][p] += 2.0 * D[s][r] * value;
	      fock[q][s] -= D[p][r] * value;
	      /* (rs|pq) */
	      fock[r][s] += 2.0 * D[p][q] * value;
	      fock[r][p] -= D[s][q] * value;
	      /* (rs|qp) */
	      fock[r][s] += 2.0 * D[q][p] * value;
	      fock[r][q] -= D[s][p] * value;
	      /* (sr|pq) */
	      fock[s][r] += 2.0 * D[p][q] * value;
	      fock[s][p] -= D[r][q] * value;
	      /* (sr|qp) */
	      fock[s][r] += 2.0 * D[q][p] * value;
	      fock[s][q] -= D[r][p] * value;
	    }
	    else if(p!=q && r!=s && pq==rs) {
	      /* (pq|sr) */
	      fock[p][q] += 2.0 * D[s][r] * value;
	      fock[p][s] -= D[q][r] * value;
	      /* (qp|rs) */
	      fock[q][p] += 2.0 * D[r][s] * value;
	      fock[q][r] -= D[p][s] * value;
	      /* (qp|sr) */
	      fock[q][p] += 2.0 * D[s][r] * value;
	      fock[q][s] -= D[p][r] * value;
	    }
	    else if(p!=q && r==s) {
	      /* (qp|rs) */
	      fock[q][p] += 2.0 * D[r][s] * value;
	      fock[q][r] -= D[p][s] * value;
	      /* (rs|pq) */
	      fock[r][s] += 2.0 * D[p][q] * value;
	      fock[r][p] -= D[s][q] * value;
	      /* (rs|qp) */
	      fock[r][s] += 2.0 * D[q][p] * value;
	      fock[r][q] -= D[s][p] * value;
	    }
	    else if(p==q && r!=s) {
	      /* (pq|sr) */
	      fock[p][q] += 2.0 * D[s][r] * value;
	      fock[p][s] -= D[q][r] * value;
	      /* (rs|pq) */
	      fock[r][s] += 2.0 * D[p][q] * value;
	      fock[r][p] -= D[s][q] * value;
	      /* (sr|pq) */
	      fock[s][r] += 2.0 * D[p][q] * value;
	      fock[s][p] -= D[r][q] * value;
	    }
	    else if(p==q && r==s && pq!=rs) {
	      /* (rs|pq) */
	      fock[r][s] += 2.0 * D[p][q] * value;
	      fock[r][p] -= D[s][q] * value;
	    }
	  } /* if(ref==0||ref==1) */
	} /* if(fzc) */

      } /* end loop through current buffer */
    } /* end loop over reading buffers */

    iwl_buf_close(&InBuf, 1); /* close buffer for next pass */

    for(h=0; h < nirreps;h++) {
      if(bucket_size[n][h])
	psio_write(File->filenum, File->label, (char *) File->matrix[h][0],
		   bucket_size[n][h]*((long int) sizeof(double)), next, &next);
      free_block(File->matrix[h]);
    }

  } /* end loop over buckets/passes */

  /* Get rid of the input integral file */
  psio_open(inputfile, PSIO_OPEN_OLD);
  psio_close(inputfile, keep);

  free_int_matrix(bucket_map,nump);

  for(n=0; n < nbuckets; n++) {
    free(bucket_offset[n]);
    free(bucket_rowdim[n]);
    free(bucket_size[n]);
  }
  free(bucket_offset);
  free(bucket_rowdim);
  free(bucket_size);

  return 0;
}
