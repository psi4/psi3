#include <stdio.h>
#include <psio.h>
#include <libciomr.h>
#include "iwl.h"
extern FILE *outfile;

int iwl_rdone_all(int itap, int ntri, double *ints, double *e_fzc, int erase);
int iwl_rdone_all_act(int itap, double *ints, double *e_fzc, int *ioff,
                      int norbs, int nfzc, int nfzv, int erase);


/*
** IWL_RDONE()
**
** This function reads the one-electron integrals in the MO basis
**  from disk and stores them in core.  Substantially revised on
**  29 April 1998 to filter out frozen orbitals if requested.
**  This change requires a very different argument list from the 
**  previous version of this code.
**
** David Sherrill, January 1994
** Revised by David Sherrill, April 1998
**
** Arguments:
**   itap       = tape to read ints from
**   ints       = buffer (already allocated) to store the integrals
**   e_fzc      = pointer to hold the frozen core energy
**   ioff       = standard offset array
**   norbs      = number of molecular orbitals
**   nfzc       = number of frozen core orbs...set to 0 if all ints desired
**   nfzv       = number of frozen virt orbs...set to 0 if all ints desired
**   erase      = erase itap (1=yes, 0=no)
**   printflg   = printing flag.  Set to 1 to print ints; otherwise, set to 0
**   outfile    = file pointer for output of ints or error messages
**
** Returns:
**   1 if reading was successful, 0 otherwise
*/
int iwl_rdone(int itap, double *ints, double *e_fzc, int *ioff,
              int norbs, int nfzc, int nfzv, int erase, 
              int printflg, FILE *outfile)
{

  int nact, ntri, stat;

  nact = norbs - nfzc - nfzv;
  ntri = (nact * (nact + 1)) / 2;

  if ((nfzc == 0) && (nfzv == 0))
    stat = iwl_rdone_all(itap,ntri,ints,e_fzc,erase);

  else 
    stat = iwl_rdone_all_act(itap,ints,e_fzc,ioff,norbs,nfzc,nfzv,erase);
 
  if (printflg) print_array(ints, nact, outfile);

  return(stat);
}



int iwl_rdone_all(int itap, int ntri, double *ints, double *e_fzc, int erase)
{
  psio_open(itap,PSIO_OPEN_OLD);
  if (psio_toclen(itap)==0) {
    fprintf(outfile, "iwl_rdone_all: Can't open one-elec integral file %d\n",
      itap);
    psio_close(itap,0);
    return(0);
  }
  
  psio_read_entry(itap, IWL_KEY_EFZC, (char *) e_fzc, sizeof(double));
  psio_read_entry(itap, IWL_KEY_ONEL, (char *) ints, ntri*sizeof(double));
  psio_close(itap, !erase);
  return(1);
}


int iwl_rdone_all_act(int itap, double *ints, double *e_fzc, int *ioff,
                      int norbs, int nfzc, int nfzv, int erase)
{
  int nact, ntri_full, ntri_act;
  int i, j, ij, ij2;
  double *tmpbuf;

  psio_open(itap, PSIO_OPEN_OLD);
  if (psio_toclen(itap)==0) {
    fprintf(outfile, "iwl_rdone_all: Can't open one-elec integral file %d\n",
      itap);
    psio_close(itap,0);
    return(0);
  }

  nact = norbs - nfzc - nfzv;
  ntri_full = (norbs * (norbs + 1)) / 2;
  ntri_act =  (nact  * (nact  + 1)) / 2;

  /* read the full-size array */
  tmpbuf = init_array(ntri_full);
  psio_read_entry(itap, IWL_KEY_EFZC, (char *) e_fzc, sizeof(double));
  psio_read_entry(itap, IWL_KEY_ONEL, (char *) tmpbuf,ntri_full*sizeof(double));
  psio_close(itap, !erase);

  /* filter out the frozen orbitals, put result in target array "ints" */
  for (i=0,ij=0; i<nact; i++) {
    for (j=0; j<=i; j++,ij++) { 
      ij2 = ioff[i+nfzc] + (j + nfzc);
      ints[ij] = tmpbuf[ij2];
    }
  }

  free(tmpbuf);

  return(1);
}

