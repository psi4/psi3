#include <stdio.h>
#include "dpd.h"

/* dpd_buf_mat_irrep_rd(): Reads an entire irrep from disk into a dpd
** two-electron buffer using the "rules" specified when the buffer was
** initialized by dpd_buf_init().
**
** Arguments:
**   struct dpdbuf *Buf: A pointer to the dpdbuf where the data will
**                       be stored.
**   int irrep: The irrep number to be read.
**   int print_flag: Boolean which indicates whether the irrep will be
**                   printed.
**   FILE *outfile: The formatted output file stream.
**
** Tested methods: 11,12,22,31,44. (12/3/96)
**
** There are definitely some problems with the read routines here.
** Mainly in the case of unpacking bra indices, there is a danger that
** the input row buffer won't be filled with new values or zeros.  See,
** e.g., method 21.
**
** T. Daniel Crawford
** December, 1996
**
** Minor modifications for newest dpd version.
** TDC
** October, 1997
*/

int dpd_buf_mat_irrep_rd(struct dpdbuf *Buf, int irrep, int print_flag,
			 FILE *outfile)
{
  int method;
  int pq, rs;  /* dpdbuf row and column indices */
  int p, q, r, s;  /* orbital indices */
  int filepq, filers, filesr;  /* Input dpdfile row and column indices */
  int rowtot, coltot;  /* dpdbuf row and column dimensions */
  int b_perm_pq, b_perm_rs, b_peq, b_res;
  int f_perm_pq, f_perm_rs, f_peq, f_res;
  int pq_permute, permute;
  double value; 

  rowtot = Buf->params->rowtot[irrep];
  coltot = Buf->params->coltot[irrep];

  b_perm_pq = Buf->params->perm_pq; b_perm_rs = Buf->params->perm_rs;
  f_perm_pq = Buf->file.params->perm_pq; f_perm_rs = Buf->file.params->perm_rs;
  b_peq = Buf->params->peq; b_res = Buf->params->res;
  f_peq = Buf->file.params->peq; f_res = Buf->file.params->res;

  if((b_perm_pq == f_perm_pq) && (b_perm_rs == f_perm_rs) &&
     (b_peq == f_peq) && (b_res == f_res)) {
      if(Buf->anti) method = 11;
      else method = 12;
      }
  else if((b_perm_pq != f_perm_pq) && (b_perm_rs == f_perm_rs) &&
	  (b_res == f_res)) {
      if(f_perm_pq && !b_perm_pq) {
	  if(Buf->anti) {
	      fprintf(outfile, "\n\tUnpack pq and antisymmetrize?\n");
	      exit(2);
	    }
	  method = 21;
	}
      else if(!f_perm_pq && b_perm_pq) {
	  if(Buf->anti) method = 22;
	  else method = 23;
	}
      else {
	  fprintf(outfile, "\n\tInvalid second-level method!\n");
	  exit(2);
	}
    }
  else if((b_perm_pq == f_perm_pq) && (b_perm_rs != f_perm_rs) &&
	  (b_peq == f_peq)) {
      if(f_perm_rs && !b_perm_rs) {
	  if(Buf->anti) {
	      fprintf(outfile, "\n\tUnpack rs and antisymmetrize?\n");
	      exit(3);
	    }
	  method = 31;
	}
      else if(!f_perm_rs && b_perm_rs) {
	  if(Buf->anti) method = 32;
	  else method = 33;
	}
      else {
	  fprintf(outfile, "\n\tInvalid third-level method!\n");
	  exit(3);
	}
    }
  else if((b_perm_pq != f_perm_pq) && (b_perm_rs != f_perm_rs)) {
      if(f_perm_pq && !b_perm_pq) {
	  if(f_perm_rs && !b_perm_rs) {
	      if(Buf->anti) {
		  fprintf(outfile, "\n\tUnpack pq and rs and antisymmetrize?\n");
		  exit(4);
		}
	      else method = 41;
	    }
	  else if(!f_perm_rs && b_perm_rs) {
	      if(Buf->anti) {
		  fprintf(outfile, "\n\tUnpack pq and antisymmetrize?\n");
		  exit(4);
		}
	      else method = 42;
	    }
	}
      else if(!f_perm_pq && b_perm_pq) {
	  if(f_perm_rs && !b_perm_rs) {
	      if(Buf->anti) {
		  fprintf(outfile, "\n\tUnpack rs and antisymmetrize?\n");
		  exit(4);
		}
	      else method = 43;
	    }
	  else if(!f_perm_rs && b_perm_rs) {
	      if(Buf->anti) method = 44;
	      else method = 45;
	    }
	}
      else {
	  fprintf(outfile, "\n\tInvalid fourth-level method!\n");
	  exit(4);
	}
    }
  else {
      fprintf(outfile, "\n\tInvalid method in dpd_buf_mat_irrep_rd!\n");
      exit(5);
    }


  if(print_flag)
      fprintf(outfile, "\n\tMethod for irrep %d, dpdbuf \"%s\": %d\n",
	      irrep, Buf->file.label, method);

  switch(method) {
  case 11: /* No change in pq or rs; antisymmetrize */
      
      /* Prepare the input buffer from the input file */
      dpd_file_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf/dpdfile */
      for(pq=0; pq < rowtot; pq++) {

	  /* Fill the buffer */
	  dpd_file_mat_irrep_row_rd(&(Buf->file), irrep, pq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep][rs][0];
	      s = Buf->params->colorb[irrep][rs][1];

	      /* Column indices in the dpdfile */
	      filers = rs;
	      filesr = Buf->file.params->colidx[s][r];

	      value = Buf->file.matrix[irrep][0][filers];

	      value -= Buf->file.matrix[irrep][0][filesr];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      if(print_flag) dpd_buf_mat_irrep_print(Buf, irrep, outfile);
  
      /* Close the input buffer */
      dpd_file_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  
  case 12: /* No change in pq or rs */
      
      Buf->file.matrix[irrep] = Buf->matrix[irrep];
      dpd_file_mat_irrep_rd(&(Buf->file), irrep, print_flag, outfile);
      
      break;
  case 21: /* Unpack pq; no change in rs */
      /* Prepare the input buffer from the input file */
      dpd_file_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf */
      for(pq=0; pq < rowtot; pq++) {
	  p = Buf->params->roworb[irrep][pq][0];
	  q = Buf->params->roworb[irrep][pq][1];
	  filepq = Buf->file.params->rowidx[p][q];

	  /* Set the permutation operator's value */
	  permute = ((p < q) && (f_perm_pq < 0) ? -1 : 1);

	  /* Fill the buffer */
	  if(filepq >= 0)
	      dpd_file_mat_irrep_row_rd(&(Buf->file), irrep, filepq);
	  else
	      dpd_file_mat_irrep_row_zero(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      filers = rs;

	      value = Buf->file.matrix[irrep][0][filers];

	      /* Assign the value, keeping track of the sign */
	      Buf->matrix[irrep][pq][rs] = permute*value;
	    }
	}

      if(print_flag) dpd_buf_mat_irrep_print(Buf, irrep, outfile);
  
      /* Close the input buffer */
      dpd_file_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  case 22: /* Pack pq; no change in rs; antisymmetrize */
      /* Prepare the input buffer from the input file */
      dpd_file_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf */
      for(pq=0; pq < rowtot; pq++) {
	  p = Buf->params->roworb[irrep][pq][0];
	  q = Buf->params->roworb[irrep][pq][1];
	  filepq = Buf->file.params->rowidx[p][q];

	  /* Fill the buffer */
	  dpd_file_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep][rs][0];
	      s = Buf->params->colorb[irrep][rs][1];

	      /* Column indices in the dpdfile */
	      filers = rs;
	      filesr = Buf->file.params->colidx[s][r];

	      value = Buf->file.matrix[irrep][0][filers];

	      value -= Buf->file.matrix[irrep][0][filesr];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      if(print_flag) dpd_buf_mat_irrep_print(Buf, irrep, outfile);
  
      /* Close the input buffer */
      dpd_file_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  case 23: /* Pack pq; no change in rs */
      /* Prepare the input buffer from the input file */
      dpd_file_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf */
      for(pq=0; pq < rowtot; pq++) {
	  p = Buf->params->roworb[irrep][pq][0];
	  q = Buf->params->roworb[irrep][pq][1];
	  filepq = Buf->file.params->rowidx[p][q];

	  dpd_file_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      filers = rs;

	      value = Buf->file.matrix[irrep][0][filers];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      if(print_flag) dpd_buf_mat_irrep_print(Buf, irrep, outfile);
  
      /* Close the input buffer */
      dpd_file_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  case 31: /* No change in pq; unpack rs */
      /* Prepare the input buffer from the input file */
      dpd_file_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf/dpdfile */
      for(pq=0; pq < rowtot; pq++) {
	  filepq = pq;

	  /* Fill the buffer */
	  dpd_file_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep][rs][0];
	      s = Buf->params->colorb[irrep][rs][1];
	      filers = Buf->file.params->colidx[r][s];

	      /* rs permutation operator */
	      permute = ((r < s) && (f_perm_rs < 0) ? -1 : 1);

	      /* Is this fast enough? */
	      value = ((filers < 0) ? 0 :
		       Buf->file.matrix[irrep][0][filers]);

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = permute*value;
	    }
	}

      if(print_flag) dpd_buf_mat_irrep_print(Buf, irrep, outfile);
  
      /* Close the input buffer */
      dpd_file_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  case 32: /* No change in pq; pack rs; antisymmetrize */
      /* Prepare the input buffer from the input file */
      dpd_file_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf/dpdfile */
      for(pq=0; pq < rowtot; pq++) {
	  filepq = pq;

	  /* Fill the buffer */
	  dpd_file_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep][rs][0];
	      s = Buf->params->colorb[irrep][rs][1];

	      /* Column indices in the dpdfile */
	      filers = Buf->file.params->colidx[r][s];
	      filesr = Buf->file.params->colidx[s][r];

	      value = Buf->file.matrix[irrep][0][filers];
	      value -= Buf->file.matrix[irrep][0][filesr];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      if(print_flag) dpd_buf_mat_irrep_print(Buf, irrep, outfile);
  
      /* Close the input buffer */
      dpd_file_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  case 33: /* No change in pq; pack rs */
      /* Prepare the input buffer from the input file */
      dpd_file_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf/dpdfile */
      for(pq=0; pq < rowtot; pq++) {
	  filepq = pq;

	  /* Fill the buffer */
	  dpd_file_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep][rs][0];
	      s = Buf->params->colorb[irrep][rs][1];
	      filers = Buf->file.params->colidx[r][s];

	      value = Buf->file.matrix[irrep][0][filers];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      if(print_flag) dpd_buf_mat_irrep_print(Buf, irrep, outfile);
  
      /* Close the input buffer */
      dpd_file_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  case 41: /* Unpack pq and rs */
      /* Prepare the input buffer from the input file */
      dpd_file_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf */
      for(pq=0; pq < rowtot; pq++) {
	  p = Buf->params->roworb[irrep][pq][0];
	  q = Buf->params->roworb[irrep][pq][1];
	  filepq = Buf->file.params->rowidx[p][q];

	  /* Set the value of the pq permutation operator */
	  pq_permute = ((p < q) && (f_perm_pq) ? -1 : 1);

	  /* Fill the buffer */
	  if(filepq >= 0)
	      dpd_file_mat_irrep_row_rd(&(Buf->file), irrep, filepq);
	  else
	      dpd_file_mat_irrep_row_zero(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep][rs][0];
	      s = Buf->params->colorb[irrep][rs][1];
	      filers = Buf->file.params->colidx[r][s];

	      /* Set the value of the pqrs permutation operator */
	      permute = ((r < s) && (f_perm_rs) ? -1 : 1)*pq_permute;

              value = 0;

	      if(filers >= 0)
		  value = Buf->file.matrix[irrep][0][filers];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = permute*value;
	    }
	}

      if(print_flag) dpd_buf_mat_irrep_print(Buf, irrep, outfile);
  
      /* Close the input buffer */
      dpd_file_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  case 42: /* Pack pq; unpack rs */
      fprintf(outfile, "\n\tHaven't programmed method 42 yet!\n");
      exit(42);

      break;
  case 43: /* Unpack pq; pack rs */
      fprintf(outfile, "\n\tHaven't programmed method 43 yet!\n");
      exit(43);

      break;
  case 44: /* Pack pq; pack rs; antisymmetrize */
      /* Prepare the input buffer from the input file */
      dpd_file_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf */
      for(pq=0; pq < rowtot; pq++) {
	  p = Buf->params->roworb[irrep][pq][0];
	  q = Buf->params->roworb[irrep][pq][1];
	  filepq = Buf->file.params->rowidx[p][q];

	  /* Fill the buffer */
	  dpd_file_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep][rs][0];
	      s = Buf->params->colorb[irrep][rs][1];

	      /* Column indices in the dpdfile */
	      filers = Buf->file.params->colidx[r][s];
	      filesr = Buf->file.params->colidx[s][r];

	      value = Buf->file.matrix[irrep][0][filers];
	      value -= Buf->file.matrix[irrep][0][filesr];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      if(print_flag) dpd_buf_mat_irrep_print(Buf, irrep, outfile);
  
      /* Close the input buffer */
      dpd_file_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  case 45: /* Pack pq and rs */
      /* Prepare the input buffer from the input file */
      dpd_file_mat_irrep_row_init(&(Buf->file), irrep);

      /* Loop over rows in the dpdbuf */
      for(pq=0; pq < rowtot; pq++) {
	  p = Buf->params->roworb[irrep][pq][0];
	  q = Buf->params->roworb[irrep][pq][1];
	  filepq = Buf->file.params->rowidx[p][q];

	  dpd_file_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

	  /* Loop over the columns in the dpdbuf */
	  for(rs=0; rs < coltot; rs++) {
	      r = Buf->params->colorb[irrep][rs][0];
	      s = Buf->params->colorb[irrep][rs][1];
	      filers = Buf->file.params->colidx[r][s];

	      if(filers < 0) {
		  fprintf(outfile, "\n\tNegative colidx in method 44?\n");
		  exit(44);
		}

	      value = Buf->file.matrix[irrep][0][filers];

	      /* Assign the value */
	      Buf->matrix[irrep][pq][rs] = value;
	    }
	}

      if(print_flag) dpd_buf_mat_irrep_print(Buf, irrep, outfile);
  
      /* Close the input buffer */
      dpd_file_mat_irrep_row_close(&(Buf->file), irrep);

      break;
  default:  /* Error trapping */
      fprintf(outfile, "\n\tInvalid switch case in dpd_buf_mat_irrep_rd!\n");
      exit(6);
      break;
    }
  
  return 0;

}
