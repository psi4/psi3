#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include <qt.h>
#include <iwl.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"

#define MAXIOFF3 255
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void make_arrays(double *evals, double **evals_docc, double **evals_virt,
                  int ndocc, int nvirt, int **ioff3, int **focact,
                  int **locact, int **fvract, int **lvract);


void energy(void)
{
  double energy=0.0;
  double escf;
  int nfile;
  double tolerance;
  struct iwlbuf NBuff;
  double *te;
  int ndocc, nvirt, ndv, nte;
  int *ioff3;
  int *clsdpi, *orbspi;
  int nirreps;
  int i,h;
  double *evals_docc, *evals_virt;
  int a,ia,j,ja,b,jb,ib,iajb,ibja;
  int asym,isym,bsym,jsym,iasym;
  int ifirst,jfirst,afirst,bfirst;
  int ilast,jlast,alast,blast;
  int print_lvl,keep_integrals;
  int *focact,*locact,*fvract,*lvract;

  escf=moinfo.escf;
  tolerance = params.tolerance;
  nfile = params.nfile;
  nirreps = moinfo.nirreps;
  clsdpi = moinfo.clsdpi;
  orbspi = moinfo.orbspi;
  print_lvl = params.print_lvl;
  keep_integrals = params.keep_integrals;

  ndocc = 0;
  nvirt = 0;
  for(h=0; h < nirreps; h++) {
      ndocc += clsdpi[h];
      nvirt += orbspi[h] - clsdpi[h];
    }

  ndv = ndocc*nvirt;
  nte = ndv*(ndv+1)/2;
  te = init_array(nte);

  fprintf(outfile, "\n\tStarting MP2 Energy Calculation...\n");
  fflush(outfile);

  make_arrays(moinfo.evals,&evals_docc,&evals_virt,ndocc,nvirt,&ioff3,
              &focact, &locact, &fvract, &lvract);
  
  iwl_buf_init(&NBuff, nfile, tolerance, 1, 1);

  iwl_buf_rd_all(&NBuff,te, ioff3, ioff3, 1, ioff, (print_lvl > 3), outfile);

  for(isym=0; isym < nirreps; isym++) {
      ifirst = focact[isym];
      ilast = locact[isym];
      for(i=ifirst; i <= ilast; i++) {
          for(asym=0; asym < nirreps; asym++) {
              afirst = fvract[asym];
              alast = lvract[asym];
              iasym = isym^asym;
              for(a=afirst; a <= alast; a++) {
                  ia = ioff3[i] + a;
                  for(jsym=0; jsym < nirreps; jsym++) {
                      jfirst = focact[jsym];
                      jlast = locact[jsym];
                      bsym = iasym^jsym;
                      bfirst = fvract[bsym];
                      blast = lvract[bsym];
                      for(j=jfirst; j <= jlast; j++) {
                          ja = ioff3[j] + a;
                          for(b=bfirst; b <= blast; b++) {
                              jb = ioff3[j] + b;
                              ib = ioff3[i] + b;
                              iajb = INDEX(ia,jb);
                              ibja = INDEX(ib,ja);
                              energy += te[iajb]*(2.0*te[iajb]-te[ibja])/
                                        (evals_docc[i] + evals_docc[j] -
                                         evals_virt[a] - evals_virt[b]);
                            }
                        }
                    }
                }
            }
        }
    }

  fprintf(outfile, "\n\tMBPT(2) Energy = %20.10lf\n", energy);
  fprintf(outfile,   "\tTotal Energy   = %20.10lf\n", escf+energy);
  fflush(outfile);

  iwl_buf_close(&NBuff, keep_integrals);
  free(te);

}

void make_arrays(double *evals, double **evals_docc, double **evals_virt,
                 int ndocc, int nvirt, int **ioff3, int **focact,
                 int **locact, int **fvract, int **lvract)
{

  int h,i,j;
  int num_docc;
  int count, offset;
  int *clsdpi;
  int *orbspi;
  int *virtpi;
  int *frdocc, *fruocc;
  int nirreps;
  int first_offset, last_offset;

  clsdpi = moinfo.clsdpi;
  orbspi = moinfo.orbspi;
  nirreps = moinfo.nirreps;
  virtpi = moinfo.virtpi;
  frdocc = moinfo.frdocc;
  fruocc = moinfo.fruocc;

  /* Re-order eigenvalues from Pitzer to QTP Relative ordering */
  (*evals_docc) = init_array(ndocc);
  count=0;
  offset=0;
  for(h=0; h < nirreps; h++) {
      if(h) offset += orbspi[h-1];
      num_docc = clsdpi[h];
      for(j=offset; j < (offset + num_docc); j++) {
          (*evals_docc)[count] = evals[j];
          count++;
        }
    }

  (*evals_virt) = init_array(nvirt);
  count=0;
  offset=0;
  for(h=0; h < nirreps; h++) {
      if(h) offset += orbspi[h-1];
      for(j=offset+clsdpi[h]; j < (offset+orbspi[h]); j++) {
          (*evals_virt)[count] = evals[j];
          count++;
        }
    }

  /* Generate ioff3 array.  This array gives the row offset for an
     ndocc x nvirt matrix */
  (*ioff3) = init_int_array(MAXIOFF3);
  for(i=0; i < MAXIOFF3; i++) {
      (*ioff3)[i] = i*nvirt;
    }


  /* Construct first and last indexing arrays for occupied and virtual
     orbitals in QTP ordering */

  *focact = init_int_array(nirreps);
  *locact = init_int_array(nirreps);
  for(h=0; h < nirreps; h++) {
      (*focact)[h] = -1;
      (*locact)[h] = -2;
    }
  first_offset=frdocc[0];
  last_offset=clsdpi[0]-1;
  (*focact)[0] = first_offset;
  (*locact)[0] = last_offset;
  for(h=1; h < nirreps; h++) {
      first_offset += frdocc[h] + clsdpi[h-1] - frdocc[h-1];
      last_offset += clsdpi[h];
      if(clsdpi[h]-frdocc[h]) {
          (*focact)[h] = first_offset;
          (*locact)[h] = last_offset;
        }
    }

  *fvract = init_int_array(nirreps);
  *lvract = init_int_array(nirreps);
  for(h=0; h < nirreps; h++) {
      (*fvract)[h] = -1;
      (*lvract)[h] = -2;
    }
  first_offset=0;
  last_offset=virtpi[0]-fruocc[0]-1;
  (*fvract)[0] = first_offset;
  (*lvract)[0] = last_offset;
  for(h=1; h < nirreps; h++) {
      first_offset += virtpi[h-1];
      last_offset += fruocc[h-1] + virtpi[h] - fruocc[h];
      if(virtpi[h]-fruocc[h]) {
          (*fvract)[h] = first_offset;
          (*lvract)[h] = last_offset;
        }
    }
}


