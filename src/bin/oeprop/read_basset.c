#define EXTERN
#include "includes.h"
#include "oeprop.gbl"
#include "oeprop.h"

void read_basset_info()
{
  int i;

  exps  = file30_rd_exps();
  contr = file30_rd_contr();
  sprim = file30_rd_sprim();
  snuc  = file30_rd_snuc();
  stype = file30_rd_stype();
  snumg = file30_rd_snumg();
  sloc  = file30_rd_sloc();
  
  lmax = 0;
  for(i=0;i<nshell;i++)
    lmax = (lmax >= stype[i]) ? lmax : stype[i];
  lmax--;

  if (print_lvl >= PRINTBASSETLEVEL) {
    fprintf(outfile,"    Prim#  \t   Exponent \tContr. coeff.\n");
    for(i=0;i<nprim;i++)
      fprintf(outfile,"  %5d  \t%12.5lf\t%13.11lf\n",i+1,exps[i],contr[i]);
    fprintf(outfile,"\n\n  Shell# Nuc#\t  L\tSPRIM\t SLOC\tSNUMG\n");
    for(i=0;i<nshell;i++)
      fprintf(outfile,"  %3d  \t%3d\t%3d\t%3d\t%3d\t%3d\t%3d\t%3d\n",i+1,
              snuc[i],stype[i]-1,sprim[i],sloc[i],snumg[i]);
    fprintf(outfile,"\n\n");
  }

}
