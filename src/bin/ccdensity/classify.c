#include <stdio.h>
#include <iwl.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

/*
** This is the ugliest piece of code I have written in my life.
** T. Daniel Crawford, September, 1996
*/

void classify(int p, int q, int r, int s, double value,
	      struct iwlbuf *ABuf, struct iwlbuf *BBuf,
	      struct iwlbuf *CBuf, struct iwlbuf *DBuf,
	      struct iwlbuf *EBuf, struct iwlbuf *FBuf)
{
  int *occ, *vir, *socc;
  int *occ_rel, *vir_rel;
  int dirac=1;
  int soccs;

  occ = moinfo.occ; vir = moinfo.vir; socc = moinfo.socc;
  occ_rel = moinfo.occ_rel; vir_rel = moinfo.vir_rel;

  soccs = socc[p] + socc[q] + socc[r] + socc[s];
 
  /* A (oo|oo) integrals */
  if((occ[p] && occ[q] && occ[r] && occ[s]))
      iwl_buf_wrt_val(ABuf, occ_rel[p], occ_rel[q], occ_rel[r], occ_rel[s],
		      value, 0, outfile, dirac);

  /* B (vv|vv) integrals */
  if((vir[p] && vir[q] && vir[r] && vir[s]))
      iwl_buf_wrt_val(BBuf, vir_rel[p], vir_rel[q], vir_rel[r], vir_rel[s],
		      value, 0, outfile, dirac);

  /* C (oo|vv) integrals */
  if(soccs > 1) {
      if((occ[p] && occ[q] && vir[r] && vir[s]))
	  iwl_buf_wrt_val(CBuf, occ_rel[p], occ_rel[q], vir_rel[r], vir_rel[s],
			  value, 0, outfile, dirac);
      if((occ[r] && occ[s] && vir[p] && vir[q]))
	  iwl_buf_wrt_val(CBuf, occ_rel[r], occ_rel[s], vir_rel[p], vir_rel[q],
			  value, 0, outfile, dirac);
    }
  else if((occ[p] && occ[q] && vir[r] && vir[s]))
      iwl_buf_wrt_val(CBuf, occ_rel[p], occ_rel[q], vir_rel[r], vir_rel[s],
		      value, 0, outfile, dirac);
  else if((occ[r] && occ[s] && vir[p] && vir[q]))
      iwl_buf_wrt_val(CBuf, occ_rel[r], occ_rel[s], vir_rel[p], vir_rel[q],
		      value, 0, outfile, dirac);

  /* D (ov|ov) integrals */
  if(soccs > 1) {
      if((occ[p] && vir[q] && occ[r] && vir[s]))
	  iwl_buf_wrt_val(DBuf, occ_rel[p], vir_rel[q], occ_rel[r], vir_rel[s],
			  value, 0, outfile, dirac);
      if((occ[q] && vir[p] && occ[r] && vir[s]))
	  iwl_buf_wrt_val(DBuf, occ_rel[q], vir_rel[p], occ_rel[r], vir_rel[s],
			  value, 0, outfile, dirac);
      if((occ[p] && vir[q] && occ[s] && vir[r]))
	  iwl_buf_wrt_val(DBuf, occ_rel[p], vir_rel[q], occ_rel[s], vir_rel[r],
			  value, 0, outfile, dirac);
      if((occ[q] && vir[p] && occ[s] && vir[r]))
	  iwl_buf_wrt_val(DBuf, occ_rel[q], vir_rel[p], occ_rel[s], vir_rel[r],
			  value, 0, outfile, dirac);
    }
  else if((occ[p] && vir[q] && occ[r] && vir[s]))
      iwl_buf_wrt_val(DBuf, occ_rel[p], vir_rel[q], occ_rel[r], vir_rel[s],
		      value, 0, outfile, dirac);
  else if((occ[q] && vir[p] && occ[r] && vir[s]))
      iwl_buf_wrt_val(DBuf, occ_rel[q], vir_rel[p], occ_rel[r], vir_rel[s],
		      value, 0, outfile, dirac);
  else if((occ[p] && vir[q] && occ[s] && vir[r]))
      iwl_buf_wrt_val(DBuf, occ_rel[p], vir_rel[q], occ_rel[s], vir_rel[r],
		      value, 0, outfile, dirac);
  else if((occ[q] && vir[p] && occ[s] && vir[r]))
      iwl_buf_wrt_val(DBuf, occ_rel[q], vir_rel[p], occ_rel[s], vir_rel[r],
		      value, 0, outfile, dirac);

  /* E (vo|oo) integrals */
  if(soccs > 1) {
      if((vir[p] && occ[q] && occ[r] && occ[s]))
	  iwl_buf_wrt_val(EBuf, vir_rel[p], occ_rel[q], occ_rel[r], occ_rel[s],
			  value, 0, outfile, dirac);
      if((vir[q] && occ[p] && occ[r] && occ[s]))
	  iwl_buf_wrt_val(EBuf, vir_rel[q], occ_rel[p], occ_rel[r], occ_rel[s],
			  value, 0, outfile, dirac);
      if((vir[r] && occ[s] && occ[p] && occ[q]))
	  iwl_buf_wrt_val(EBuf, vir_rel[r], occ_rel[s], occ_rel[p], occ_rel[q],
			  value, 0, outfile, dirac);
      if((vir[s] && occ[r] && occ[p] && occ[q]))
	  iwl_buf_wrt_val(EBuf, vir_rel[s], occ_rel[r], occ_rel[p], occ_rel[q],
			  value, 0, outfile, dirac);
    } 
  else if((vir[p] && occ[q] && occ[r] && occ[s]))
      iwl_buf_wrt_val(EBuf, vir_rel[p], occ_rel[q], occ_rel[r], occ_rel[s],
		      value, 0, outfile, dirac);
  else if((vir[p] && occ[q] && occ[s] && occ[r]))
      iwl_buf_wrt_val(EBuf, vir_rel[q], occ_rel[p], occ_rel[r], occ_rel[s],
		      value, 0, outfile, dirac);
  else if((vir[r] && occ[s] && occ[p] && occ[q]))
      iwl_buf_wrt_val(EBuf, vir_rel[r], occ_rel[s], occ_rel[p], occ_rel[q],
		      value, 0, outfile, dirac);
  else if((vir[r] && occ[s] && occ[q] && occ[p]))
      iwl_buf_wrt_val(EBuf, vir_rel[s], occ_rel[r], occ_rel[p], occ_rel[q], 
		      value, 0, outfile, dirac);

  /* F (ov|vv) integrals */
  if(soccs > 1) {
      if((occ[p] && vir[q] && vir[r] && vir[s]))
	  iwl_buf_wrt_val(FBuf, occ_rel[p], vir_rel[q], vir_rel[r], vir_rel[s],
			  value, 0, outfile, dirac);
      if((occ[q] && vir[p] && vir[r] && vir[s]))
	  iwl_buf_wrt_val(FBuf, occ_rel[q], vir_rel[p], vir_rel[r], vir_rel[s],
			  value, 0, outfile, dirac);
      if((occ[r] && vir[s] && vir[p] && vir[q]))
	  iwl_buf_wrt_val(FBuf, occ_rel[r], vir_rel[s], vir_rel[p], vir_rel[q],
			  value, 0, outfile, dirac);
      if((occ[s] && vir[r] && vir[p] && vir[q]))
	  iwl_buf_wrt_val(FBuf, occ_rel[s], vir_rel[r], vir_rel[p], vir_rel[q],
			  value, 0, outfile, dirac);
    }
  else if((occ[p] && vir[q] && vir[r] && vir[s]))
      iwl_buf_wrt_val(FBuf, occ_rel[p], vir_rel[q], vir_rel[r], vir_rel[s],
		      value, 0, outfile, dirac);
  else if((occ[q] && vir[p] && vir[r] && vir[s]))
      iwl_buf_wrt_val(FBuf, occ_rel[q], vir_rel[p], vir_rel[r], vir_rel[s],
		      value, 0, outfile, dirac);
  else if((occ[r] && vir[s] && vir[p] && vir[q]))
      iwl_buf_wrt_val(FBuf, occ_rel[r], vir_rel[s], vir_rel[p], vir_rel[q],
		      value, 0, outfile, dirac);
  else if((occ[s] && vir[r] && vir[p] && vir[q]))
      iwl_buf_wrt_val(FBuf, occ_rel[s], vir_rel[r], vir_rel[p], vir_rel[q],
		      value, 0, outfile, dirac);
}

