#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libciomr/libciomr.h>
#include <libfile30/file30.h>
#include <libint/libint.h>

#include "defines.h"
#define EXTERN
#include "global.h"

void file11()
{
  int i;
  char wfnstring[80];
  double **Geom, **GradRef, **GeomRef;
  FILE *fp11;

  /*--- Ok, have to be sure we get these energies, reload if necessary ---*/
  /*--- should have Escf around at least ---*/
  if (strcmp(UserOptions.wfn,"SCF")==0) {
    MOInfo.Eref = MOInfo.Escf;
  }
  else {
    MOInfo.Eref = file30_rd_eref();
    MOInfo.Ecorr = file30_rd_ecorr();
  }

  /*--- Geometry in the canonical frame ---*/
  Geom = block_matrix(Molecule.num_atoms,3);
  /*--- Geometry and gradient in the reference frame ---*/
  GeomRef = block_matrix(Molecule.num_atoms,3);
  GradRef = block_matrix(Molecule.num_atoms,3);
  
  /*--- Rotate back to the reference frame ---*/
  for(i=0;i<Molecule.num_atoms;i++) {
      Geom[i][0] = Molecule.centers[i].x;
      Geom[i][1] = Molecule.centers[i].y;
      Geom[i][2] = Molecule.centers[i].z;
  }
  mmult(Geom,0,Molecule.Rref,0,GeomRef,0,Molecule.num_atoms,3,3,0);
  mmult(Grad,0,Molecule.Rref,0,GradRef,0,Molecule.num_atoms,3,3,0);

  
  sprintf(wfnstring,"%s forces in the reference frame (a.u.)",UserOptions.wfn);
  print_atomvec(wfnstring,GradRef);
  ffile(&fp11, "file11.dat", 1);

  fprintf(fp11,"%-59.59s %-10.10s%-8.8s\n",Molecule.label,UserOptions.wfn,UserOptions.dertype);

  fprintf(fp11,"%5d",Molecule.num_atoms);
  if (strcmp(UserOptions.wfn,"SCF") == 0)
    fprintf(fp11,"%20.10lf\n",MOInfo.Escf);
  else
    fprintf(fp11,"%20.10lf\n",MOInfo.Eref+MOInfo.Ecorr);
  
  for(i=0;i<Molecule.num_atoms;i++)
    fprintf(fp11,"%20.10lf%20.10lf%20.10lf%20.10lf\n",
	    Molecule.centers[i].Z_nuc,GeomRef[i][0],GeomRef[i][1],GeomRef[i][2]);
  for(i=0;i<Molecule.num_atoms;i++)
    fprintf(fp11,"                    %20.10lf%20.10lf%20.10lf\n",
	    GradRef[i][0],GradRef[i][1],GradRef[i][2]);
  fclose(fp11);

  free_block(Geom);
  free_block(GeomRef);
  free_block(GradRef);
  
  return;
}
  
