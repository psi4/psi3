#include <stdio.h>
#include <iwl.h>
#include <ip_libv1.h>
#include <libciomr.h>
#define EXTERN
#include "globals.h"

void classify(int p, int q, int r, int s, double value,
	      struct iwlbuf *ABuf, struct iwlbuf *CBuf,
	      struct iwlbuf *DBuf, struct iwlbuf *EBuf)
{
  int *occ, *vir, *socc;
  int *cc_occ, *cc_vir;
  int dirac=1;
  int soccs;

  occ = moinfo.occ; vir = moinfo.vir; socc = moinfo.socc;
  cc_occ = moinfo.cc_occ; cc_vir = moinfo.cc_vir;

  soccs = socc[p] + socc[q] + socc[r] + socc[s];
 
  /* A (oo|oo) integrals */
  if((occ[p] && occ[q] && occ[r] && occ[s]))
      iwl_buf_wrt_val(ABuf, cc_occ[p], cc_occ[q], cc_occ[r], cc_occ[s],
		      value, 0, outfile, dirac);

  /* C (oo|vv) integrals */
  if(soccs > 1) {
      if((occ[p] && occ[q] && vir[r] && vir[s]))
	  iwl_buf_wrt_val(CBuf, cc_occ[p], cc_occ[q], cc_vir[r], cc_vir[s],
			  value, 0, outfile, dirac);
      if((occ[r] && occ[s] && vir[p] && vir[q]))
	  iwl_buf_wrt_val(CBuf, cc_occ[r], cc_occ[s], cc_vir[p], cc_vir[q],
			  value, 0, outfile, dirac);
    }
  else if((occ[p] && occ[q] && vir[r] && vir[s]))
      iwl_buf_wrt_val(CBuf, cc_occ[p], cc_occ[q], cc_vir[r], cc_vir[s],
		      value, 0, outfile, dirac);
  else if((occ[r] && occ[s] && vir[p] && vir[q]))
      iwl_buf_wrt_val(CBuf, cc_occ[r], cc_occ[s], cc_vir[p], cc_vir[q],
		      value, 0, outfile, dirac);

  /* D (ov|ov) integrals */
  if(soccs > 1) {
      if((occ[p] && vir[q] && occ[r] && vir[s]))
	  iwl_buf_wrt_val(DBuf, cc_occ[p], cc_vir[q], cc_occ[r], cc_vir[s],
			  value, 0, outfile, dirac);
      if((occ[q] && vir[p] && occ[r] && vir[s]))
	  iwl_buf_wrt_val(DBuf, cc_occ[q], cc_vir[p], cc_occ[r], cc_vir[s],
			  value, 0, outfile, dirac);
      if((occ[p] && vir[q] && occ[s] && vir[r]))
	  iwl_buf_wrt_val(DBuf, cc_occ[p], cc_vir[q], cc_occ[s], cc_vir[r],
			  value, 0, outfile, dirac);
      if((occ[q] && vir[p] && occ[s] && vir[r]))
	  iwl_buf_wrt_val(DBuf, cc_occ[q], cc_vir[p], cc_occ[s], cc_vir[r],
			  value, 0, outfile, dirac);
    }
  else if((occ[p] && vir[q] && occ[r] && vir[s]))
      iwl_buf_wrt_val(DBuf, cc_occ[p], cc_vir[q], cc_occ[r], cc_vir[s],
		      value, 0, outfile, dirac);
  else if((occ[q] && vir[p] && occ[r] && vir[s]))
      iwl_buf_wrt_val(DBuf, cc_occ[q], cc_vir[p], cc_occ[r], cc_vir[s],
		      value, 0, outfile, dirac);
  else if((occ[p] && vir[q] && occ[s] && vir[r]))
      iwl_buf_wrt_val(DBuf, cc_occ[p], cc_vir[q], cc_occ[s], cc_vir[r],
		      value, 0, outfile, dirac);
  else if((occ[q] && vir[p] && occ[s] && vir[r]))
      iwl_buf_wrt_val(DBuf, cc_occ[q], cc_vir[p], cc_occ[s], cc_vir[r],
		      value, 0, outfile, dirac);

  /* E (vo|oo) integrals */
  if(soccs > 1) {
      if((vir[p] && occ[q] && occ[r] && occ[s]))
	  iwl_buf_wrt_val(EBuf, cc_vir[p], cc_occ[q], cc_occ[r], cc_occ[s],
			  value, 0, outfile, dirac);
      if((vir[q] && occ[p] && occ[r] && occ[s]))
	  iwl_buf_wrt_val(EBuf, cc_vir[q], cc_occ[p], cc_occ[r], cc_occ[s],
			  value, 0, outfile, dirac);
      if((vir[r] && occ[s] && occ[p] && occ[q]))
	  iwl_buf_wrt_val(EBuf, cc_vir[r], cc_occ[s], cc_occ[p], cc_occ[q],
			  value, 0, outfile, dirac);
      if((vir[s] && occ[r] && occ[p] && occ[q]))
	  iwl_buf_wrt_val(EBuf, cc_vir[s], cc_occ[r], cc_occ[p], cc_occ[q],
			  value, 0, outfile, dirac);
    } 
  else if((vir[p] && occ[q] && occ[r] && occ[s]))
      iwl_buf_wrt_val(EBuf, cc_vir[p], cc_occ[q], cc_occ[r], cc_occ[s],
		      value, 0, outfile, dirac);
  else if((vir[p] && occ[q] && occ[s] && occ[r]))
      iwl_buf_wrt_val(EBuf, cc_vir[q], cc_occ[p], cc_occ[r], cc_occ[s],
		      value, 0, outfile, dirac);
  else if((vir[r] && occ[s] && occ[p] && occ[q]))
      iwl_buf_wrt_val(EBuf, cc_vir[r], cc_occ[s], cc_occ[p], cc_occ[q],
		      value, 0, outfile, dirac);
  else if((vir[r] && occ[s] && occ[q] && occ[p]))
      iwl_buf_wrt_val(EBuf, cc_vir[s], cc_occ[r], cc_occ[p], cc_occ[q], 
		      value, 0, outfile, dirac);

}

