#include <stdio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#define EXTERN
#include "globals.h"

/* get_moinfo(): Routine to obtain basic orbital information from
** file30 and compute the associated lookup arrays.
**
** T. Daniel Crawford, October 1996
** Modified for STABLE by TDC March, 2000 */

void get_moinfo(void)
{
  int i, j, h, count, ocount, vcount, errcod, col;
  int cl_offset, op_offset, vr_offset;
  int nopen, nclsd;
  double enuc, escf;
  
  chkpt_init();
  moinfo.nirreps = chkpt_rd_nirreps();
  moinfo.nmo = chkpt_rd_nmo();
  moinfo.iopen = chkpt_rd_iopen();
  moinfo.labels = chkpt_rd_irr_labs();
  enuc = chkpt_rd_enuc();
  escf = chkpt_rd_escf();
  moinfo.orbspi = chkpt_rd_orbspi();
  moinfo.clsdpi = chkpt_rd_clsdpi();
  moinfo.openpi = chkpt_rd_openpi();
  chkpt_close();

  /* Number of occupied and virtual orbirals per irrep including
     open-shells */
  moinfo.occpi = init_int_array(moinfo.nirreps);
  moinfo.virtpi = init_int_array(moinfo.nirreps);
  moinfo.uoccpi = init_int_array(moinfo.nirreps);
  nclsd = 0; nopen = 0;
  for(h=0; h < moinfo.nirreps; h++) {
      moinfo.occpi[h] = moinfo.clsdpi[h] + moinfo.openpi[h];
      moinfo.virtpi[h] = moinfo.orbspi[h] - moinfo.clsdpi[h];
      moinfo.uoccpi[h] = moinfo.orbspi[h] - moinfo.clsdpi[h] - moinfo.openpi[h];
      nclsd += moinfo.clsdpi[h];
      nopen += moinfo.openpi[h];
    }
 
  moinfo.frdocc = init_int_array(moinfo.nirreps);
  moinfo.fruocc = init_int_array(moinfo.nirreps);
  errcod = ip_int_array("FROZEN_DOCC", moinfo.frdocc, moinfo.nirreps);
  errcod = ip_int_array("FROZEN_UOCC", moinfo.fruocc, moinfo.nirreps);

  moinfo.nfzc = moinfo.nfzv = 0;
  for(i=0; i < moinfo.nirreps; i++) {
      moinfo.nfzc += moinfo.frdocc[i];
      moinfo.nfzv += moinfo.fruocc[i];
    }

  /* Build the boolean arrays for the orbital classification routines.
     The occ and vir arrays identify active occupied or virtual
     orbitals (including open shells).  The socc array identifies only
     the singly occupied orbitals.  The argument for all arrays much
     be a QT-ordered index. Similar arrays are used in the integral
     sorting routines found in CCSORT and CCDENSITY. */

  moinfo.occ = init_int_array(moinfo.nmo);
  moinfo.vir = init_int_array(moinfo.nmo);
  moinfo.socc = init_int_array(moinfo.nmo);

  count = 0;
  for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.clsdpi[i]; j++, count++) {
	  moinfo.occ[count] = 1;
	}
    }
  for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.openpi[i]; j++, count++) {
	  moinfo.occ[count] = 1;
	  moinfo.vir[count] = 1;
	  moinfo.socc[count] = 1;
	}
    }
  for(i=0; i < moinfo.nirreps; i++) {
      for(j=0; j < moinfo.uoccpi[i]; j++, count++) {
	  moinfo.vir[count] = 1;
	}
    }

  /* Build relative index arrays which map a QT Standard orbital index
     into the CC occupied or virtual relative index.  In the CC
     occupied ordering (given by "cc_occ"), the orbitals are organized
     with all frozen core orbitals first, followed by doubly occupied
     orbitals, followed by singly occupied orbitals.  In the CC
     virtual ordering (given by "cc_vir"), the orbitals are organized
     with all frozen virtual orbitals first, followed by active
     unoccupied orbitals, followed by the singly occupied orbitals.
     These arrays are needed by classification routines used in
     sorting two-electron integrals where integrals are divided up by
     their occupied and virtual index patterns (see classify.c and
     distribute.c).
    
     For example, for 3B1 CH2 in a DZ basis with one frozen core
     orbital and three frozen virtual orbitals:
    
     QT:      0   1   2   3   4   5   6   7   8   9  10  11  12  13
     Space:  fc   d   d   s   s   u   u   u   u   u   u  fv  fv  fv
     Symm:   a1  a1  b2  a1  b1  a1  a1  a1  b1  b2  b2  a1  a1  b2
     --------------------------------------------------------------
     cc_occ:  0   1   4   2   3  -1  -1  -1  -1  -1  -1  -1  -1  -1
     cc_vir: -1  -1  -1   5   7   2   3   4   6   9  10   0   1   8

     A "-1" is used to mark QT indices irrelevant to the given CC
     ordering array.

     Note that I also compute the CC->QT versions of the cc_occ and
     cc_vir arrays (qt_occ and qt_vir, respectively) simultaneously.
     
     */

  /* CC ordering and symmetry arrays */
  moinfo.cc_occ = init_int_array(moinfo.nmo);
  moinfo.cc_vir = init_int_array(moinfo.nmo);
  moinfo.qt_occ = init_int_array(moinfo.nmo);
  moinfo.qt_vir = init_int_array(moinfo.nmo);
  for(i=0; i < moinfo.nmo; i++) {
      moinfo.cc_occ[i] = moinfo.cc_vir[i] = -1;
      moinfo.qt_occ[i] = moinfo.qt_vir[i] = -1;
    }

  moinfo.occ_sym = init_int_array(moinfo.nmo);
  moinfo.vir_sym = init_int_array(moinfo.nmo);

  count = 0;
  cl_offset = 0;
  op_offset = nclsd;
  for(h=0; h < moinfo.nirreps; h++) {
      if(h) cl_offset += moinfo.clsdpi[h-1];
      for(i=0; i < moinfo.clsdpi[h]; i++,count++) {
	  moinfo.cc_occ[cl_offset+i] = count;
	  moinfo.qt_occ[count] = cl_offset+i;
	  moinfo.occ_sym[count] = h;
	}
      if(h) op_offset += moinfo.openpi[h-1];
      for(i=0; i < moinfo.openpi[h]; i++,count++) {
	  moinfo.cc_occ[op_offset+i] = count;
	  moinfo.qt_occ[count] = op_offset+i;
	  moinfo.occ_sym[count] = h;
	}
    }

  count = 0;
  vr_offset = nclsd + nopen;
  op_offset = nclsd;
  for(h=0; h < moinfo.nirreps; h++) {
      if(h) vr_offset += moinfo.uoccpi[h-1];
      for(i=0; i < moinfo.uoccpi[h]; i++,count++) {
	  moinfo.cc_vir[vr_offset+i] = count;
	  moinfo.qt_vir[count] = vr_offset+i;
	  moinfo.vir_sym[count] = h;
	}
      if(h) op_offset += moinfo.openpi[h-1];
      for(i=0; i < moinfo.openpi[h]; i++,count++) {
	  moinfo.cc_vir[op_offset+i] = count;
	  moinfo.qt_vir[count] = op_offset+i;
	  moinfo.vir_sym[count] = h;
	}
    }

  /* Calculate occupied and virtual orbital offsets within each irrep */
  moinfo.occ_off = init_int_array(moinfo.nirreps);
  moinfo.vir_off = init_int_array(moinfo.nirreps);
  ocount = moinfo.occpi[0]; vcount = moinfo.virtpi[0];
  for(h=1; h < moinfo.nirreps; h++) {
      moinfo.occ_off[h] = ocount;
      ocount += moinfo.occpi[h];
      moinfo.vir_off[h] = vcount;
      vcount += moinfo.virtpi[h];
    }

  /*
  fprintf(outfile,"\n\tChkpt Parameters:\n");
  fprintf(outfile,"\t------------------\n");
  fprintf(outfile,"\tNumber of irreps     = %d\n",moinfo.nirreps);
  fprintf(outfile,"\tNumber of MOs        = %d\n",moinfo.nmo);
  fprintf(outfile,"\tNumber of active MOs = %d\n\n",moinfo.nactive);
  fprintf(outfile,
      "\tLabel\t# MOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
  fprintf(outfile,
      "\t-----\t-----\t------\t------\t------\t------\t------\n");
  for(i=0; i < moinfo.nirreps; i++) {
      fprintf(outfile,
              "\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n",
              moinfo.labels[i],moinfo.orbspi[i],moinfo.frdocc[i],
              moinfo.clsdpi[i],moinfo.openpi[i],moinfo.uoccpi[i],
              moinfo.fruocc[i]);
    }
  fprintf(outfile,"\n\tNuclear Rep. energy (chkpt) = %20.15f\n", moinfo.enuc);
  fprintf(outfile,  "\tSCF energy          (chkpt) = %20.15f\n", escf);
  */
}
