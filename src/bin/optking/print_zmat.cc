/*** PRINT_ZMAT() printout z-matrix ***/ 

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
}

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"

void compute_zmat(cartesians &carts, int *unique_zvars);
void print_zmat(FILE *outfile, int *unique_zvars);

/* recompute values for z-matrix and write back to chkpt */
void compute_zmat(cartesians &carts, int *unique_zvars) {
  int i, a, cnt;
  int nallatom, natom, *to_nodummy, *atom_dummy;
  struct internals *simples;
  struct z_entry *zmat;
  char sym[30];
  double *opt_fgeom;

  nallatom = optinfo.nallatom;
  chkpt_init(PSIO_OPEN_OLD);
  zmat = chkpt_rd_zmat();
  chkpt_close();

  /* determine and save the unique variables */
  cnt = -1;
  for (i=1; i<nallatom; ++i) {
    unique_zvars[++cnt] = 1;
    strcpy(sym, zmat[i].bond_label);
    for (a=0;a<i;++a) {
      if (strcmp(sym, zmat[a].bond_label) == 0) {
        unique_zvars[cnt] = 0;
        break;
      }
    }
    if (i>1) {
      unique_zvars[++cnt] = 1;
      strcpy(sym, zmat[i].angle_label);
      for (a=0;a<i;++a) {
        if (strcmp(sym, zmat[a].angle_label) == 0) {
          unique_zvars[cnt] = 0;
          break;
        }
      }
    }
    if (i>2) {
      unique_zvars[++cnt] = 1;
      strcpy(sym, zmat[i].tors_label);
      for (a=0;a<i;++a) {
        if (strcmp(sym, zmat[a].tors_label) == 0) {
          unique_zvars[cnt] = 0;
          break;
        }
      }
    }
  }

  int *nints;
  nints = (int *) malloc(5*sizeof(int));
  nints[0] = nallatom-1;
  nints[1] = nallatom-2;
  nints[2] = nallatom-3;
  nints[3] = 0;
  internals zints(nints);
  /* compute the value of the unique variables */
  zints.stre.set_num(nallatom-1);
  zints.bend.set_num(nallatom-2);
  zints.tors.set_num(nallatom-3);
  zints.out.set_num(0);

  cnt = 0;
  for (i=0;i<nallatom;++i) {
    if (i>0) {
      zints.stre.set_id(i-1,++cnt);
      zints.stre.set_A(i-1,i);
      zints.stre.set_B(i-1,zmat[i].bond_atom-1);
    }
    if (i>1) {
      zints.bend.set_id(i-2,++cnt);
      zints.bend.set_A(i-2,i);
      zints.bend.set_B(i-2,zmat[i].bond_atom-1);
      zints.bend.set_C(i-2,zmat[i].angle_atom-1);
    }
    if (i>2) {
      zints.tors.set_id(i-3,++cnt);
      zints.tors.set_A(i-3,i);
      zints.tors.set_B(i-3,zmat[i].bond_atom-1);
      zints.tors.set_C(i-3,zmat[i].angle_atom-1);
      zints.tors.set_D(i-3,zmat[i].tors_atom-1);
    }
  }
  // compute value of the zmatrix coordinates
  opt_fgeom = carts.get_fcoord();
  zints.compute_internals(nallatom, opt_fgeom);
  // zints.print(outfile, 1);

  // insert computed values into zmatrix object
  for (i=0;i<nallatom;++i) {
    if (i>0) {
      zmat[i].bond_val = zints.stre.get_val(i-1);
    }
    if (i>1) {
      zmat[i].angle_val = zints.bend.get_val(i-2);
    }
    if (i>2) {
      zmat[i].tors_val = zints.tors.get_val(i-3);
    }
  }

  // write recomputed z-matrix to file30
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_zmat(zmat);
  chkpt_close();
}

/*** ZMAT_TO_INTCO if simples are not already there and zmat_simples
* is not turned off, then generate simple internals from zmatrix ***/
void print_zmat(FILE *outfile, int *unique_zvars) {
  int i, a, b, c, d, cnt = 0;
  int nallatom, natom, *to_nodummy, *atom_dummy;
  char **felement;
  char buf[2], sym[3];
  double *zvals;
  struct z_entry *zmat;

  nallatom = optinfo.nallatom;
  natom = optinfo.natom;
  atom_dummy = optinfo.atom_dummy;
  to_nodummy = optinfo.to_nodummy;

  chkpt_init(PSIO_OPEN_OLD);
  zmat = chkpt_rd_zmat();
  felement = chkpt_rd_felement();
  zvals = chkpt_rd_zvals();
  chkpt_close();

  fprintf(outfile,"  zmat = ( \n");
  for (i=0; i<nallatom; ++i) {
    if (atom_dummy[i]) strcpy(sym,"X");
    else zval_to_symbol(zvals[to_nodummy[i]], sym);
    fprintf(outfile,"    ( %s ", sym);
    if (i > 0) {
      fprintf(outfile," %d", zmat[i].bond_atom);
      fprintf(outfile," %s", zmat[i].bond_label);
      if (zmat[i].bond_opt) fprintf(outfile,"$");
    }
    if (i > 1) {
      fprintf(outfile," %d", zmat[i].angle_atom);
      fprintf(outfile," %s", zmat[i].angle_label);
      if (zmat[i].angle_opt) fprintf(outfile,"$");
    }
    if (i > 2) {
      fprintf(outfile," %d", zmat[i].tors_atom);
      fprintf(outfile," %s$", zmat[i].tors_label);
    }
    fprintf(outfile,")\n");
  }
  fprintf(outfile,"  )\n");
  cnt = -1;
  fprintf(outfile,"  zvars = ( \n");
  for (i=0; i<nallatom; ++i) {
    if (i > 0) {
      if (unique_zvars[++cnt])
        fprintf(outfile,"    ( %s %10.5lf )\n",
            zmat[i].bond_label, zmat[i].bond_val);
    }
    if (i > 1) {
      if (unique_zvars[++cnt])
        fprintf(outfile,"    ( %s %10.5lf )\n",
            zmat[i].angle_label, zmat[i].angle_val);
    }
    if (i > 2) {
      if (unique_zvars[++cnt])
        fprintf(outfile,"    ( %s %10.5lf )\n",
            zmat[i].tors_label, zmat[i].tors_val);
    }
  }
  fprintf(outfile,"  )\n");
  return;
}

/*
void zval_to_symbol(double zval, char *sym) {
  int z;
  z = (int) zval;

  if (z==0) strcpy(sym,"G");
  else if (z==1) strcpy(sym,"H"); 
  else if (z==2) strcpy(sym,"HE"); 
  else if (z==3) strcpy(sym,"LI"); 
  else if (z==4) strcpy(sym,"BE"); 
  else if (z==5) strcpy(sym,"B"); 
  else if (z==6) strcpy(sym,"C"); 
  else if (z==7) strcpy(sym,"N"); 
  else if (z==8) strcpy(sym,"O"); 
  else if (z==9) strcpy(sym,"F"); 
  else if (z==10) strcpy(sym,"NE"); 
  else if (z==11) strcpy(sym,"NA"); 
  else if (z==12) strcpy(sym,"MG"); 
  else if (z==13) strcpy(sym,"AL"); 
  else if (z==14) strcpy(sym,"SI"); 
  else if (z==15) strcpy(sym,"P"); 
  else if (z==16) strcpy(sym,"S"); 
  else if (z==17) strcpy(sym,"CL"); 
  else if (z==18) strcpy(sym,"AR"); 
  else if (z==19) strcpy(sym,"K"); 
  else if (z==20) strcpy(sym,"CA"); 
  else if (z==21) strcpy(sym,"SC"); 
  else if (z==22) strcpy(sym,"TI"); 
  else if (z==23) strcpy(sym,"V"); 
  else if (z==24) strcpy(sym,"CR"); 
  else if (z==25) strcpy(sym,"MN"); 
  else if (z==26) strcpy(sym,"FE"); 
  else if (z==27) strcpy(sym,"CO"); 
  else if (z==28) strcpy(sym,"NI"); 
  else if (z==29) strcpy(sym,"CU"); 
  else if (z==30) strcpy(sym,"ZN"); 
  else if (z==31) strcpy(sym,"GA"); 
  else if (z==32) strcpy(sym,"GE"); 
  else if (z==33) strcpy(sym,"AS"); 
  else if (z==34) strcpy(sym,"SE"); 
  else if (z==35) strcpy(sym,"BR"); 
  else if (z==36) strcpy(sym,"KR"); 
  else if (z==37) strcpy(sym,"RB"); 
  else if (z==38) strcpy(sym,"SR"); 
  else if (z==39) strcpy(sym,"Y"); 
  else if (z==40) strcpy(sym,"ZR"); 
  else if (z==41) strcpy(sym,"NB");
  else if (z==42) strcpy(sym,"MO");
  else if (z==43) strcpy(sym,"TC");
  else if (z==44) strcpy(sym,"RU");
  else if (z==45) strcpy(sym,"RH");
  else if (z==46) strcpy(sym,"PD");
  else if (z==47) strcpy(sym,"AG");
  else if (z==48) strcpy(sym,"CD");
  else if (z==49) strcpy(sym,"IN"); 
  else if (z==50) strcpy(sym,"SN");
  else if (z==51) strcpy(sym,"SB");
  else if (z==52) strcpy(sym,"TE");
  else if (z==53) strcpy(sym,"I");
  else if (z==54) strcpy(sym,"XE");
  else if (z==55) strcpy(sym,"CS");
  else if (z==56) strcpy(sym,"BA");
  else if (z==57) strcpy(sym,"LA");
  else if (z==58) strcpy(sym,"CE");
  else if (z==59) strcpy(sym,"PR");
  else if (z==60) strcpy(sym,"ND");
  else if (z==61) strcpy(sym,"PM");
  else if (z==62) strcpy(sym,"SM");
  else if (z==63) strcpy(sym,"EU");
  else if (z==64) strcpy(sym,"GD");
  else if (z==65) strcpy(sym,"TB");
  else if (z==66) strcpy(sym,"DY");
  else if (z==67) strcpy(sym,"HO");
  else if (z==68) strcpy(sym,"ER");
  else if (z==69) strcpy(sym,"TM");
  else if (z==70) strcpy(sym,"TY");
  else if (z==71) strcpy(sym,"LU");
  else if (z==72) strcpy(sym,"HF");
  else if (z==73) strcpy(sym,"TA");
  else if (z==74) strcpy(sym,"W");
  else if (z==75) strcpy(sym,"RE");
  else if (z==76) strcpy(sym,"OS");
  else if (z==77) strcpy(sym,"IR");
  else if (z==78) strcpy(sym,"PT");
  else if (z==79) strcpy(sym,"AU");
  else if (z==80) strcpy(sym,"HG");
  else if (z==81) strcpy(sym,"TL");
  else if (z==82) strcpy(sym,"PB");
  else if (z==83) strcpy(sym,"BI");
  else if (z==84) strcpy(sym,"PO");
  else if (z==85) strcpy(sym,"AT");
  else if (z==86) strcpy(sym,"RN");
  else if (z==87) strcpy(sym,"FR");
  else if (z==88) strcpy(sym,"RA");
  else if (z==89) strcpy(sym,"AC");
  else if (z==90) strcpy(sym,"TH");
  else if (z==91) strcpy(sym,"PA");
  else if (z==92) strcpy(sym,"U");
  else if (z==93) strcpy(sym,"NP");
  else if (z==94) strcpy(sym,"PU");
  else if (z==95) strcpy(sym,"AM");
  else if (z==96) strcpy(sym,"CM");
  else if (z==97) strcpy(sym,"BK");
  else if (z==98) strcpy(sym,"CF");
  else if (z==99) strcpy(sym,"ES");
  else if (z==100) strcpy(sym,"FM");
  else if (z==101) strcpy(sym,"MD");
  else if (z==102) strcpy(sym,"NO");
  else if (z==103) strcpy(sym,"UNQ");
  else if (z==104) strcpy(sym,"UNP");
  else if (z==105) strcpy(sym,"UNH");
  else if (z==106) strcpy(sym,"UNS");
  return ;
}
*/
