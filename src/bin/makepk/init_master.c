
/* $Log$
 * Revision 1.1  2000/02/04 22:51:33  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:52:27  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:06:29  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

/* allocate space on disk for master file  */

/* the following is from fortran MASTER *********************************
C
C
C  NO.  :  CONTENTS
C
C   TRIA. = LOWER TRIANGULAR MATRIX
C   SQUA. = SQUARE MATRIX
C   RECT. = RECTANGULAR MATRIX
C   1   : LOCATION
C   2   : PARAMETERS
C   3   : ENERGIES AND REAL CONSTANTS
C   4   : MO ORDERING IN SCF
C   5   : MO ORDERING IN DRT
C   6   : MO CODE FOR CI
C   7   : MO CODE FOR MCSCF
C   8   : NUCLEAR CHARGE AND GEOMETRY
C   9   : AO-SO TRANSFORMATION MATRIX (RECT.)
C  10   : EIGENVALUES IN PITZER'S SCF (SQUA.)
C  11   : OCCUPATION IN PITZER'S SCF
C  12   : SO-MO EIGENVECTORS IN PITZER'S SCF (SQUA.)
C  13   : AO-MO EIGENVECTORS IN PITZER'S SCF (RECT.)
C  14   : SO OVERLAP INTEGRALS IN PITZER'S SCF (TRIA.)
C  15   : ONE ELECTRON SO INTEGRALS IN PITZER'S SCF (TRIA.)
C  16   : ONE ELECTRON MO INTEGRALS IN PITZER'S SCF (TRIA.)
C  17   : EIGENVALUES IN SORTED SCF
C  18   : OCCUPATION IN SORTED SCF
C  19   : SO-MO EIGENVECTORS IN SORTED SCF (SQUA.)
C  20   : AO-MO EIGENVECTORS IN SORTED SCF (RECT.)
C  21   : SO OVERLAP INTEGRALS IN SORTED SCF (TRIA.)
C  22   : ONE ELECTRON SO INTEGRALS IN SORTED SCF (TRIA.)
C  23   : ONE ELECTRON MO INTEGRLAS IN SORTED SCF (TRIA.)
C  24   : LAGRANGIAN MATRIX FOR SCF IN AO BASIS (TRIA.)
C  25   : LAGRANGIAN MATRIX FOR SCF IN MO BASIS (TRIA.)
C  26   : K MATRIX FOR HIGH SPIN OPEN-SHELL SCF (TRIA.)
C  27   : FIRST ZETA MATRIX  FOR GENERALIZED OPEN-SHELL SCF (TRIA.)
C  28   : SECOND ZETA MATRIX FOR GENERALIZED OPEN-SHELL SCF (TRIA.)
C  29   : THIRD ZETA MATRIX  FOR GENERALIZED OPEN-SHELL SCF (TRIA.)
C  30   : FOURTH ZETA MATRIX FOR GENERALIZED OPEN-SHELL SCF (TRIA.)
C  31   : FIFTH ZETA MATRIX  FOR GENERALIZED OPEN-SHELL SCF (TRIA.)
C  32   : LAGRANGIAN MATRIX FOR CI IN AO BASIS (TRIA.)
C  33   : LAGRANGIAN MATRIX FOR CI IN MO BASIS (SQUA.)
C  34   : ONE PDM IN AO BASIS FOR CI (TRIA.)
C  35   : ONE PDM IN MO BASIS FOR CI (TRIA.)
C  36   : LAGRANGIAN MATRIX FOR MCSCF IN AO BASIS (TRIA.)
C  37   : LAGRANGIAN MATRIX FOR MCSCF IN MO BASIS (SQUA.)
C  38   : ONE PDM IN AO BASIS FOR MCSCF (TRIA.)
C  39   : ONE PDM IN MO BASIS FOR MCSCF (TRIA.)
c  40   : alp matrix in scf order
c  41   : bet matrix in scf order
c  42   : ao-mo eigenvectors in drt order

*********************************************************************/

init_master()

   {
      int i,j,nsecth,size;
      int nsect2,nbasl,ngeom,ngeoml,nbas2,nbas2l;
      int nbasq,nbasql,nbasr,nbasrl,ntri2,ntril;
      int nbatr2,nbatrl,nbasa,nbasal;
      int *dum;

      rfile(itap40);

      nsect=1024;

      size = MAX0(nsect,nbfao*nbfao*2);

      block_locs = (int *) init_array(nsect);
      dum = (int *) init_array(size);

      nsecth=nsect/2;
      nsect2=nsect*2;
      nbasl = (nbfso-1)/nsect + 1;
      ngeom = natom*8;
      ngeoml = (ngeom-1)/nsect + 1;
      nbas2 = nbfso*2;
      nbas2l = (nbas2-1)/nsect + 1;
      nbasq = nbfso*nbfso*2;
      nbasql = (nbasq-1)/nsect + 1;
      nbasr = nbfao*nbfso*2;
      nbasrl = (nbasr-1)/nsect + 1;
      ntri2 = nbstri*2;
      ntril = (ntri2-1)/nsect + 1;
      nbatr2 = nbatri*2;
      nbatrl = (nbatr2-1)/nsect + 1;
      nbasa = nbfao*nbfao*2;
      nbasal = (nbasa-1)/nsect + 1;

      dum[0] = 1;
      dum[nsecth] = nsect;
      dum[1] = dum[0]+1;
      dum[1+nsecth] = nsect;
      dum[2] = dum[1]+1;
      dum[2+nsecth] = nsect2;
      dum[3] = dum[2]+2;
      dum[3+nsecth] = nbfso;
      dum[4] = dum[3]+nbasl;
      dum[4+nsecth] = nbfso;
      dum[5] = dum[4]+nbasl;
      dum[5+nsecth] = nbfso;
      dum[6] = dum[5]+nbasl;
      dum[6+nsecth] = nbfso;
      dum[7] = dum[6]+nbasl;
      dum[7+nsecth] = ngeom;
      dum[8] = dum[7]+ngeoml;
      dum[8+nsecth] = nbasr;
      dum[9] = dum[8]+nbasrl;
      dum[9+nsecth] = nbas2;
      dum[10] = dum[9]+nbas2l;
      dum[10+nsecth] = nbas2;
      dum[11] = dum[10]+nbas2l;
      dum[11+nsecth]= nbasq;
      dum[12] = dum[11]+nbasql;
      dum[12+nsecth] = nbasr;
      dum[13] = dum[12]+nbasrl;
      dum[13+nsecth]= ntri2;
      dum[14] = dum[13]+ntril;
      dum[14+nsecth] = ntri2;
      dum[15] = dum[14]+ntril;
      dum[15+nsecth] = ntri2;
      dum[16] = dum[15]+ntril;
      dum[16+nsecth] = nbas2;
      dum[17] = dum[16]+nbas2l;
      dum[17+nsecth] = nbas2;
      dum[18] = dum[17]+nbas2l;
      dum[18+nsecth] = nbasq;
      dum[19] = dum[18]+nbasql;
      dum[19+nsecth] = nbasr;
      dum[20] = dum[19]+nbasrl;
      dum[20+nsecth] = ntri2;
      dum[21] = dum[20]+ntril;
      dum[21+nsecth] = ntri2;
      dum[22] = dum[21]+ntril;
      dum[22+nsecth] = ntri2;
      dum[23] = dum[22]+ntril;
      dum[23+nsecth] = nbatr2;
      dum[24] = dum[23]+nbatrl;
      dum[24+nsecth] = ntri2;
      dum[25] = dum[24]+ntril;
      dum[25+nsecth] = ntri2;
      dum[26] = dum[25]+ntril;
      dum[26+nsecth] = ntri2;
      dum[27] = dum[26]+ntril;
      dum[27+nsecth] = ntri2;
      dum[28] = dum[27]+ntril;
      dum[28+nsecth] = ntri2;
      dum[29] = dum[28]+ntril;
      dum[29+nsecth] = ntri2;
      dum[30] = dum[29]+ntril;
      dum[30+nsecth] = ntri2;
      dum[31] = dum[30]+ntril;
      dum[31+nsecth] = nbatr2;
      dum[32] = dum[31]+nbatrl;
      dum[32+nsecth] = nbasq;
      dum[33] = dum[32]+nbasql;
      dum[33+nsecth] = nbatr2;
      dum[34] = dum[33]+nbatrl;
      dum[34+nsecth] = ntri2;
      dum[35] = dum[34]+ntril;
      dum[35+nsecth] = nbatr2;
      dum[36] = dum[35]+nbatrl;
      dum[36+nsecth] = nbasq;
      dum[37] = dum[36]+nbasql;
      dum[37+nsecth] = nbatr2;
      dum[38] = dum[37]+nbatrl;
      dum[38+nsecth] = ntri2;
      dum[39] = dum[38]+nbatrl;
      dum[39+nsecth] = ntri2;
      dum[40] = dum[39]+nbatrl;
      dum[40+nsecth] = ntri2;
      dum[41] = dum[40]+nbasrl;
      dum[41+nsecth] = nbasr;
 
      for(i=0; i < nsect ; i++) block_locs[i]=dum[i];

      rwrit(itap40,(char *) dum,sizeof(int)*nsect,block_locs[0]);

      dum[0] = nsect;
      dum[1] = 0;
      dum[2] = nbfso;
      dum[3] = nbstri;
      dum[5] = (twocon) ? -iopen : iopen;
      dum[9] = natom;
      dum[10] = nbfao;
      dum[11] = nbfso;
      dum[12] = n_so_typs;
      dum[16] = natom*3;

      rwrit(itap40,(char *) dum,sizeof(int)*nsect,block_locs[1]);
      rwrit(itap40,(char *) dum,sizeof(int)*nsect2,block_locs[2]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbfso,block_locs[3]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbfso,block_locs[4]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbfso,block_locs[5]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbfso,block_locs[6]);
      rwrit(itap40,(char *) dum,sizeof(int)*ngeom,block_locs[7]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbasr,block_locs[8]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbas2,block_locs[9]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbas2,block_locs[10]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbasq,block_locs[11]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbasr,block_locs[12]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[13]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[14]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[15]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbas2,block_locs[16]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbas2,block_locs[17]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbasq,block_locs[18]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbasr,block_locs[19]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[20]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[21]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[22]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbatr2,block_locs[23]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[24]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[25]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[26]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[27]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[28]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[29]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[30]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbatr2,block_locs[31]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbasq,block_locs[32]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbatr2,block_locs[33]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[34]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbatr2,block_locs[35]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbasq,block_locs[36]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbatr2,block_locs[37]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[38]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[39]);
      rwrit(itap40,(char *) dum,sizeof(int)*ntri2,block_locs[40]);
      rwrit(itap40,(char *) dum,sizeof(int)*nbasr,block_locs[41]);

      }
