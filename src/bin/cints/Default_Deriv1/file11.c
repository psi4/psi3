#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libint.h>

#include "defines.h"
#define EXTERN
#include "global.h"

void file11()
{
  int i;
  char wfnstring[80];
  FILE *fp11;

  sprintf(wfnstring,"%s forces (a.u.)",UserOptions.wfn);
  print_atomvec(wfnstring,Grad);
  ffile(&fp11, "file11.dat", 1);
  fprintf(fp11,"%-59.59s %-10.10s%-8.8s\n",Molecule.label,UserOptions.wfn,UserOptions.dertype);
  fprintf(fp11,"%5d%20.10lf\n",Molecule.num_atoms,MOInfo.Escf);
  for(i=0;i<Molecule.num_atoms;i++)
    fprintf(fp11,"%20.10lf%20.10lf%20.10lf%20.10lf\n",
	    Molecule.centers[i].Z_nuc,Molecule.centers[i].x,Molecule.centers[i].y,Molecule.centers[i].z);
  for(i=0;i<Molecule.num_atoms;i++)
    fprintf(fp11,"                    %20.10lf%20.10lf%20.10lf\n",
	    Grad[i][0],Grad[i][1],Grad[i][2]);
  fclose(fp11);

  return;
}
  
