#include <stdio.h>
#include <stdlib.h>
extern "C" {
#include <libciomr/libciomr.h>
#include <libfile30/file30.h>
#include <libpsio/psio.h>
#include <psifiles.h>
}
#include "molecule.h"
#include "params.h"

extern Molecule_t Molecule;
extern Params_t Params;
extern FILE *outfile;

extern void done(const char *);

void eval_dboc_rhf()
{
  int disp = 1;  // Displacement counter
  const int ndisp = Params.disp_per_coord * Molecule.natom * 3;
  char *inputcmd = new char[40];

  for(disp=1; disp<=ndisp; disp+=2) {
    sprintf(inputcmd,"input --getgeom %d",disp);
    int errcod = system(inputcmd);
    if (errcod) {
      done("input failed");
    }
    errcod = system("psi");
    if (errcod) {
      done("psi failed");
    }
    disp++;
    sprintf(inputcmd,"input --getgeom %d --savemos",disp);
    errcod = system(inputcmd);
    if (errcod) {
      done("input failed");
    }
    errcod = system("psi");
    if (errcod) {
      done("psi failed");
    }

    file30_init();
    int num_mo = file30_rd_nmo();
    int num_so = file30_rd_nso();
    double **rhf_evec_p = file30_rd_scf();
    file30_close();
    
    psio_open(PSIF_OLD_CHKPT, PSIO_OPEN_OLD);
    double **rhf_evec_m = block_matrix(num_so,num_mo);
    psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:MOs alpha", (char *)&(rhf_evec_m[0][0]),
		    num_so*num_mo*sizeof(double));
    double **Spm = block_matrix(num_so,num_so);
    psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:New-Old basis overlap", (char *)&(Spm[0][0]),
		    num_so*num_so*sizeof(double));
    psio_close(PSIF_OLD_CHKPT, 1);
    
    double **tmpmat1 = block_matrix(num_mo,num_so);
    mmult(rhf_evec_p,1,Spm,0,tmpmat1,0,
	  num_mo,num_so,num_so,0);
    double **CSC = block_matrix(num_mo,num_mo);
    mmult(tmpmat1,0,rhf_evec_m,0,CSC,0,
	  num_mo,num_so,num_mo,0);
    free_block(tmpmat1);
    
    fprintf(outfile,"  -Cp*Spm*Cm :\n");
    print_mat(CSC,num_mo,num_mo,outfile);
    free_block(CSC);
  }

}
