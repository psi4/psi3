#include <stdio.h>
#include <stdlib.h>
#include <psifiles.h>
#include <libciomr.h>
#include <iwl.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))

/* Globals needed for the post-backtransform sort */
int nbuckets;        /* number of sorting buckets   */
int *shell;          /* AO -> shell                 */
int *shell_size;     /* AOs in shell                */
int *pair_size;      /* Length of each shell pair   */
int *bucket_map;     /* shell-pair -> sort bucket   */
int *bucket_offset;  /* bucket -> quartet offset    */
int *bucket_quarts;  /* no. of quartets in a bucket */
int *bucket_firstpq; /* First pq in bucket          */
int *bucket_lastpq;  /* Last pq in bucket           */

void backsort_prep(void)
{
  int i, j, k, l, lend, size_i, size_j, size_k, size_l, size_ij;
  int this_pair_size, this_pair_sizeb, this_pair_quarts,core_left;
  int nshell, *snuc, num_pairs;
  int ii, ij, ijkl, imin, imax;

  snuc = moinfo.snuc;
  nshell = moinfo.nshell;
  num_pairs = ioff[nshell];

  shell_size = init_int_array(nshell);
  for(i=0; i < nshell; i++) shell_size[i] = ioff[moinfo.stype[i]];

  pair_size = init_int_array(num_pairs);
  for(i=0,ij=0; i < nshell; i++) {
      size_i = shell_size[i];
      for(j=0; j <= i; j++,ij++) {
          size_j = shell_size[j];

          pair_size[ij] = size_i * size_j;
	}
    }

  bucket_map = init_int_array(nshell*(nshell+1)/2);
  shell = init_int_array(moinfo.nao);

  /* Make room for one bucket to begin with */
  bucket_offset = init_int_array(1);
  bucket_quarts = init_int_array(1);
  bucket_firstpq = init_int_array(1);
  bucket_lastpq = init_int_array(1);
  
  /* Figure out how many buckets we need and where each */
  /* shell-pair (in canonical shell-ordering) goes.     */
  core_left = params.maxcor;
  nbuckets = 1;
  for(i=0,ij=0,ijkl=0; i < nshell; i++) {
      size_i = shell_size[i];

      /* Generate the orbital -> shell lookup while we're here */
      imin = moinfo.sloc[i] - 1;
      imax = imin + shell_size[i];
      for(ii=imin; ii < imax; ii++) shell[ii] = i;
      
      for(j=0; j <= i; j++,ij++) {
	  size_j = shell_size[j];
	  size_ij = pair_size[ij];

	  /* Number of density elements in this pair */
	  this_pair_size = 0;

	  /* Number of quartets in this pair */
	  this_pair_quarts = 0;
		   
	  for(k=0; k <= i; k++) {
	      size_k = shell_size[k];
	      lend = (k==i) ? j : k;
	      for(l=0; l <= lend; l++,ijkl++) {
		  size_l = shell_size[l];

		  this_pair_size += size_ij * size_k * size_l;
		  this_pair_quarts++;

		} /* l loop */
	    } /* k loop */

	  /* Add 4 indices to size and convert to bytes */
	  this_pair_sizeb = this_pair_size*(4*sizeof(int) + sizeof(double));
	  
	  if((core_left - this_pair_sizeb) >= 0) {
	      core_left -= this_pair_sizeb;
	      bucket_quarts[nbuckets-1] += this_pair_quarts;
	      bucket_lastpq[nbuckets-1] = ij;
	    }
	  else {
	      nbuckets++;
	      core_left = params.maxcor - this_pair_sizeb;
	      bucket_offset = (int *) realloc((void *) bucket_offset,
					      nbuckets * sizeof(int));
	      bucket_offset[nbuckets-1] = bucket_offset[nbuckets-2] +
					  bucket_quarts[nbuckets-2];

	      bucket_firstpq = (int *) realloc((void *) bucket_firstpq,
			 		       nbuckets * sizeof(int));
	      bucket_firstpq[nbuckets-1] = ij;
	      
	      bucket_lastpq = (int *) realloc((void *) bucket_lastpq,
					       nbuckets * sizeof(int));
	      bucket_lastpq[nbuckets-1] = ij;
	      
	      bucket_quarts = (int *) realloc((void *) bucket_quarts,
					      nbuckets * sizeof(int));
	      bucket_quarts[nbuckets-1] = this_pair_quarts;
	    }

	  bucket_map[ij] = nbuckets - 1;

	} /* j loop */
    } /* i loop */

}

void backsort(int first_tmp_file, double tolerance)
{
  int i, n, lastbuf, idx, nquarts, this_quartet, this_counter, *counter;
  int *snuc, sp, sq, sr, ss, send, size_r, size_s, size_pq;
  int p, q, r, s, pshell, qshell, rshell, sshell, pqshell, rsshell, pqrs;
  int pq, rs;
  int **pidx, **qidx, **ridx, **sidx;
  double **gamma, value;
  struct iwlbuf InBuf, OutBuf;
  int num_tpdm, quartet_size, lastq;

  snuc = moinfo.snuc;

  iwl_buf_init(&OutBuf, params.mfile, tolerance, 0, 0);

  num_tpdm = 0;
  for(n=0; n < nbuckets; n++) {
      iwl_buf_init(&InBuf, first_tmp_file+n, tolerance, 1, 0);
      lastbuf = 0;

      nquarts = bucket_quarts[n];
      
      pidx = (int **) malloc(nquarts * sizeof(int *));
      qidx = (int **) malloc(nquarts * sizeof(int *));
      ridx = (int **) malloc(nquarts * sizeof(int *));
      sidx = (int **) malloc(nquarts * sizeof(int *));
      gamma = (double **) malloc(nquarts * sizeof(double *));
      counter = init_int_array(nquarts);

      /* Compute quartet sizes for this bucket and allocate space */
      for(pq=bucket_firstpq[n],nquarts=0; pq <= bucket_lastpq[n]; pq++) {

	  size_pq = pair_size[pq];

	  for(r=0,rs=0; r < moinfo.nshell; r++) {
	      size_r = shell_size[r];
	      for(s=0; s <= r; s++,rs++) {
		  size_s = shell_size[s];

		  if(rs > pq) break;

		  quartet_size = size_pq * size_r * size_s;

		  pidx[nquarts] = init_int_array(quartet_size);
		  qidx[nquarts] = init_int_array(quartet_size);
		  ridx[nquarts] = init_int_array(quartet_size);
		  sidx[nquarts] = init_int_array(quartet_size);
		  gamma[nquarts] = init_array(quartet_size);
		  nquarts++;
		}
	    }
	}

      if(nquarts != bucket_quarts[n]) {
	  printf("Quartet error: nquarts = %d; bucket_quarts[%d] = %d\n",
		 nquarts, n, bucket_quarts[n]);
	  exit(20);
	}

      while(!lastbuf) {
	  iwl_buf_fetch(&InBuf);
	  lastbuf = InBuf.lastbuf;

	  for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
	      p = (int) InBuf.labels[idx++];
	      q = (int) InBuf.labels[idx++];
	      r = (int) InBuf.labels[idx++];
	      s = (int) InBuf.labels[idx++];

	      value = (double) InBuf.values[InBuf.idx];

	      pshell = shell[p]; qshell = shell[q];
	      rshell = shell[r]; sshell = shell[s];

	      /* Skip this quartet if on one center */
	      if(snuc[pshell] == snuc[qshell] &&
		 snuc[pshell] == snuc[rshell] &&
		 snuc[pshell] == snuc[sshell])
		  continue;

	      pqshell = INDEX(pshell,qshell); rsshell = INDEX(rshell,sshell);
	      pqrs = INDEX(pqshell,rsshell);

	      pq = INDEX(p,q);  rs = INDEX(r,s);

	      this_quartet = pqrs - bucket_offset[n];
	      this_counter = counter[this_quartet];

	      pidx[this_quartet][this_counter] = p;
	      qidx[this_quartet][this_counter] = q;
	      ridx[this_quartet][this_counter] = r;
	      sidx[this_quartet][this_counter] = s;
	      gamma[this_quartet][this_counter] = value;

	      counter[this_quartet]++;

	      /* Now run through the appropriate permutations */
	      if(pqshell != rsshell) {
		  if(pshell != qshell) {
		      if(rshell == sshell) {
			  if(qshell != rshell) { /* (pq|rr) */
			      if(r!=s) {
				  this_counter = counter[this_quartet];
				  pidx[this_quartet][this_counter] = p;
				  qidx[this_quartet][this_counter] = q;
				  ridx[this_quartet][this_counter] = s;
				  sidx[this_quartet][this_counter] = r;
				  gamma[this_quartet][this_counter] = value;
				  counter[this_quartet]++;
				}
			    }
			  else {  /* (pq|qq) */
			      if(r!=s) {
				  this_counter = counter[this_quartet];
				  pidx[this_quartet][this_counter] = p;
				  qidx[this_quartet][this_counter] = q;
				  ridx[this_quartet][this_counter] = s;
				  sidx[this_quartet][this_counter] = r;
				  gamma[this_quartet][this_counter] = value;
				  counter[this_quartet]++;
				}
			    }
			}
		    }
		  else {
		      if(rshell != sshell) {
			  if(qshell != rshell) { /* (pp|rs) */
			      if(p!=q) {
				  this_counter = counter[this_quartet];
				  pidx[this_quartet][this_counter] = q;
				  qidx[this_quartet][this_counter] = p;
				  ridx[this_quartet][this_counter] = r;
				  sidx[this_quartet][this_counter] = s;
				  gamma[this_quartet][this_counter] = value;
				  counter[this_quartet]++;
				}
			    }
			  else { /* (pp|ps) */
			      if(p!=q) {
				  this_counter = counter[this_quartet];
				  pidx[this_quartet][this_counter] = q;
				  qidx[this_quartet][this_counter] = p;
				  ridx[this_quartet][this_counter] = r;
				  sidx[this_quartet][this_counter] = s;
				  gamma[this_quartet][this_counter] = value;
				  counter[this_quartet]++;
				}
			    }
			}
		      else {
			  if(qshell != rshell) { /* (pp|rr) */
			      if(p!=q && r!=s) {
				  this_counter = counter[this_quartet];
				  pidx[this_quartet][this_counter] = q;
				  qidx[this_quartet][this_counter] = p;
				  ridx[this_quartet][this_counter] = r;
				  sidx[this_quartet][this_counter] = s;
				  gamma[this_quartet][this_counter] = value;
				  counter[this_quartet]++;
				  
				  this_counter = counter[this_quartet];
				  pidx[this_quartet][this_counter] = p;
				  qidx[this_quartet][this_counter] = q;
				  ridx[this_quartet][this_counter] = s;
				  sidx[this_quartet][this_counter] = r;
				  gamma[this_quartet][this_counter] = value;
				  counter[this_quartet]++;
				  
				  this_counter = counter[this_quartet];
				  pidx[this_quartet][this_counter] = q;
				  qidx[this_quartet][this_counter] = p;
				  ridx[this_quartet][this_counter] = s;
				  sidx[this_quartet][this_counter] = r;
				  gamma[this_quartet][this_counter] = value;
				  counter[this_quartet]++;
				}
			      else if(p!=q && r==s) {
				  this_counter = counter[this_quartet];
				  pidx[this_quartet][this_counter] = q;
				  qidx[this_quartet][this_counter] = p;
				  ridx[this_quartet][this_counter] = r;
				  sidx[this_quartet][this_counter] = s;
				  gamma[this_quartet][this_counter] = value;
				  counter[this_quartet]++;
				}
			      else if(p==q && r!=s) {
				  this_counter = counter[this_quartet];
				  pidx[this_quartet][this_counter] = p;
				  qidx[this_quartet][this_counter] = q;
				  ridx[this_quartet][this_counter] = s;
				  sidx[this_quartet][this_counter] = r;
				  gamma[this_quartet][this_counter] = value;
				  counter[this_quartet]++;
				}
			    }
			  else {
			      exit(1);
			    }
			}
		    }
		}
	      else { /* pqshell == rsshell */
		  if(pshell != qshell) { /* (pq|pq) */
		      if(pq != rs) {
			  this_counter = counter[this_quartet];
			  pidx[this_quartet][this_counter] = r;
			  qidx[this_quartet][this_counter] = s;
			  ridx[this_quartet][this_counter] = p;
			  sidx[this_quartet][this_counter] = q;
			  gamma[this_quartet][this_counter] = value;
			  counter[this_quartet]++;
			}
		    }
		  else { /* (pp|pp) */
		      /* This shouldn't actually occur because of the
			 snuc[] filter above */
		      exit(1);
		    }
		}
	    }
	}

      iwl_buf_close(&InBuf, 0);

      /* Now flush each quartet in order */
      for(pq=bucket_firstpq[n],nquarts=0; pq <= bucket_lastpq[n]; pq++) {
	  size_pq = pair_size[pq];

	  for(r=0,rs=0; r < moinfo.nshell; r++) {
	      size_r = shell_size[r];
	      for(s=0; s <= r; s++,rs++) {
		  size_s = shell_size[s];

		  if(rs > pq) break;

		  iwl_buf_wrt_arr(&OutBuf, gamma[nquarts],
				  pidx[nquarts], qidx[nquarts],
				  ridx[nquarts], sidx[nquarts],
				  counter[nquarts]);

		  /*
		  for(i=0; i < counter[nquarts]; i++) {
		      fprintf(outfile, "%d %d %d %d gamma = %20.12f\n", 
			      pidx[nquarts][i], qidx[nquarts][i],
			      ridx[nquarts][i], sidx[nquarts][i],
			      gamma[nquarts][i]);
		    }
		    */

		  /* Mark the end of the shell quartet */
		  if(counter[nquarts]) {
		      iwl_buf_wrt_val(&OutBuf, -1, -1, -1, -1,
				      9.9999999999, 0, outfile, 0);
		      /*
		      fprintf(outfile, "-1 -1 -1 -1 gamma = 9.9999999999\n");
		      */
		    }

/*		  num_tpdm += counter[nquarts]; */

		  nquarts++;
		}
	    }
	}

      /* Free the sort arrays */
      for(i=0; i < bucket_quarts[n]; i++) {
	  free(pidx[i]);  free(qidx[i]); free(ridx[i]);  free(sidx[i]);
	  free(gamma[i]);
	}
      free(pidx); free(qidx); free(ridx); free(sidx);
      free(gamma); free(counter);
      
    } /* end of bucket loop */

  iwl_buf_flush(&OutBuf, 1);
  iwl_buf_close(&OutBuf, 1);

/*  fprintf(outfile, "num_tpdm = %d\n", num_tpdm); */

  /* Free up the global data */
  free(shell);
  free(shell_size);
  free(pair_size);
  free(bucket_map);
  free(bucket_offset);
  free(bucket_quarts);
  free(bucket_firstpq);
  free(bucket_lastpq);
}
