// This function forms a set of delocalized, symmetry-adapted
// internal coordinates, given a set of simples and their
// s-vectors.  The new coordinates will not mix coordinates of
// different types if mix_types == 1

#if HAVE_CMATH
# include <cmath>
#else
# include <math.h>
#endif

extern "C" {
#include <stdio.h>
#include <libchkpt/chkpt.h>
#include <stdlib.h>
#include <string.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
}

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

extern double **irrep(internals &simples, double **evectst);

void rm_rotations(internals &simples, cartesians &carts, int &num_nonzero, double **evects);

void delocalize(internals &simples, cartesians &carts) {
  int error,i,j,k,a,b,c,d,id,count,intco_type,sub_index,row[4],dim[4];
  int rotor_type, degrees_of_freedom, col, natom;
  double **stre_mat, **bend_mat, **tors_mat, **out_mat;
  double **evectst, **evectst_symm, **coord_symm, *fmass;
  double *evals, **evects, **B, **uBt, **BBt, **u, **temp_mat;

  natom = optinfo.natom;

  // Build B matrix for simples
  B = block_matrix(simples.get_num(),natom*3);
  count = -1;
  for (i=0;i<simples.stre.get_num();++i) {
    a = simples.stre.get_A(i);
    b = simples.stre.get_B(i);
    ++count;
    for (k=0;k<3;++k) {
      B[count][3*a+k] += simples.stre.get_s_A(i,k);
      B[count][3*b+k] += simples.stre.get_s_B(i,k);
    }
  }
  for (i=0;i<simples.bend.get_num();++i) {
    a = simples.bend.get_A(i);
    b = simples.bend.get_B(i);
    c = simples.bend.get_C(i);
    ++count;
    for (k=0;k<3;++k) {
      B[count][3*a+k] += simples.bend.get_s_A(i,k);
      B[count][3*b+k] += simples.bend.get_s_B(i,k);
      B[count][3*c+k] += simples.bend.get_s_C(i,k);
    }
  }
  for (i=0;i<simples.tors.get_num();++i) {
    a = simples.tors.get_A(i);
    b = simples.tors.get_B(i);
    c = simples.tors.get_C(i);
    d = simples.tors.get_D(i);
    ++count;
    for (k=0;k<3;++k) {
      B[count][3*a+k] += simples.tors.get_s_A(i,k);
      B[count][3*b+k] += simples.tors.get_s_B(i,k);
      B[count][3*c+k] += simples.tors.get_s_C(i,k);
      B[count][3*d+k] += simples.tors.get_s_D(i,k);
    }
  }
  for (i=0;i<simples.out.get_num();++i) {
    a = simples.out.get_A(i);
    b = simples.out.get_B(i);
    c = simples.out.get_C(i);
    d = simples.out.get_D(i);
    ++count;
    for (k=0;k<3;++k) {
      B[count][3*a+k] += simples.out.get_s_A(i,k);
      B[count][3*b+k] += simples.out.get_s_B(i,k);
      B[count][3*c+k] += simples.out.get_s_C(i,k);
      B[count][3*d+k] += simples.out.get_s_D(i,k);
    }
  }
  //  print_mat2(B,simples.get_num(),natom*3,outfile);

  uBt = block_matrix(3*natom, simples.get_num());
  fmass = carts.get_fmass();
  u = mass_mat(fmass);
  free(fmass);

  // Form BBt matrix
  BBt = block_matrix(simples.get_num(),simples.get_num());
  if (optinfo.mix_types) {
    mmult(u,0,B,1,uBt,0,3*natom,3*natom,simples.get_num(),0);

    mmult(B,0,uBt,0,BBt,0,simples.get_num(),3*natom,simples.get_num(),0);

    // mmult(B,0,B,1,BBt,0,simples.get_num(),natom*3,simples.get_num(),0);
  }
  else {
    // Make BBt block diagonal by multiplying only stre*stre, etc.
    dim[0] = simples.stre.get_num();
    dim[1] = simples.bend.get_num();
    dim[2] = simples.tors.get_num();
    dim[3] = simples.out.get_num();
    row[0] = 0;
    row[1] = dim[0];
    row[2] = row[1]+dim[1];
    row[3] = row[2]+dim[2];
    double **ptr;
    ptr = (double **) malloc(simples.get_num()*sizeof(double *));
    for (i=0;i<4;++i) {
      for (j=0;j<dim[i];++j) {
        ptr[j] = BBt[row[i]+j] + row[i];
      }
      if (dim[i] != 0)

        mmult(u,0,&(B[row[i]]),1,uBt,0,3*natom,3*natom,dim[i],0);
      mmult(&(B[row[i]]),0,uBt,0,ptr,0,dim[i],3*natom,dim[i],0);

      //  mmult(&(B[row[i]]),0,&(B[row[i]]),1,ptr,0,dim[i],natom*3,dim[i],0);
    }
    free(ptr);
  }
  free_block(u);
  free_block(B);
  free_block(uBt);

  //fprintf(outfile,"The BB^t Matrix:\n");
  //print_mat2(BBt,simples.get_num(),simples.get_num(),outfile);

  // Diagonalize BBt
  evals = init_array(simples.get_num());
  evects = block_matrix(simples.get_num(),simples.get_num());
  sq_rsp(simples.get_num(),simples.get_num(), BBt, evals, 1, evects, 1.0E-14);

  /* check eigenvectors of BBt */
  if (optinfo.print_delocalize) {
    fprintf(outfile,"\n\n+++Delocalized Coordinate Formation+++\n");
    fprintf(outfile,"\n\nChecking eigenvectors of BBt...\n");
    eivout(evects, evals, simples.get_num(), simples.get_num(), outfile);
  }
  temp_mat = block_matrix(simples.get_num(),simples.get_num());
  mmult(BBt,0,evects,0,temp_mat,0,simples.get_num(),simples.get_num(),simples.get_num(),0);
  for (j=0;j<simples.get_num();++j) {
    error = 0;
    for (i=0;i<simples.get_num();++i) {
      if ( fabs(temp_mat[i][j] - evals[j]*evects[i][j]) > 1.0E-13)  error = 1;
      if (error == 1) { fprintf(outfile,"Error in evect %d\n",j); error = 0;}
    }
  }
  free_block(temp_mat);
  free_block(BBt);

  /* check for proper number of non-redundant coordinates (3n-6) */ 
  num_nonzero = 0;
  for (i=0;i<simples.get_num();++i)
    if( evals[i] > optinfo.ev_tol ) ++num_nonzero;

  chkpt_init(PSIO_OPEN_OLD);
  rotor_type = chkpt_rd_rottype();
  chkpt_close();

  switch (rotor_type) {
    case 3:
      degrees_of_freedom = 3 * natom - 5;
      break;
    case 6:
      degrees_of_freedom = 0;
      break;
    default:
      degrees_of_freedom = 3 * natom - 6;
      break;
  }
  if (num_nonzero == degrees_of_freedom) {
    fprintf(outfile,"\nGood: # of delocalized coordinates = # of degrees");
    fprintf(outfile," of freedom.\n");
  }
  else if (num_nonzero < degrees_of_freedom) {
    fprintf(outfile,"# of delocalized coordinates < # of degrees of ");
    fprintf(outfile,"freedom!\n You may need more simple internals.\n");
  }
  else if (num_nonzero > degrees_of_freedom) {
    fprintf(outfile,"# of delocalized coordinates > # of degrees of ");
    fprintf(outfile,"freedom!\n You might try increasing EV_TOL.\n");
    rm_rotations(simples, carts, num_nonzero, evects);
    fprintf(outfile,"Eigenvectors of BBt after removal of rotations\n");
    eivout(evects, evals, simples.get_num(), simples.get_num(), outfile);
  }


  /* transpose evects matrix, also throw out redundant coordinates
     (eigenvectors corresponding to zero eigenvalues*/

  char aline[MAX_LINELENGTH], **buffer, *err;
  int h;
  irr = init_int_array(simples.get_num());

  evectst = block_matrix(num_nonzero,simples.get_num());

  for (i=0;i<simples.get_num();++i) {
    for (j=simples.get_num()-num_nonzero;j<simples.get_num();++j) {
      evectst[j-simples.get_num()+num_nonzero][i] = evects[i][j];
    }
  }
  if (optinfo.print_delocalize == 1) {
    fprintf(outfile,"\nNon-redundant delocalized internal coordinates");
    fprintf(outfile,"(each row is a coordinate):\n");
    print_mat(evectst,num_nonzero,simples.get_num(),outfile);
  }
  free(evals);
  free_block(evects);

  /* send evectst to irrep.cc to be properly symmetrized */
  evectst_symm = irrep(simples, evectst);
  free_block(evectst);

  if (optinfo.print_delocalize == 1) {
    fprintf(outfile,"\nSymmetrized evects\n");
    print_mat(evectst_symm,num_nonzero,simples.get_num(),outfile);
  }

  // print out coordinates to intco.dat
  ffile(&fp_intco, "intco.dat", 2);
  count = 0;
  for( ; ; ) {
    err = fgets(aline, MAX_LINELENGTH, fp_intco);
    if (err == NULL) break;
    ++count;
  }
  rewind(fp_intco);

  /* read all but the last line of the file into memory... */
  buffer = (char **) malloc((count-1) * sizeof(char *));
  for(i=0; i < count-1; i++) {
    buffer[i] = (char *) malloc(MAX_LINELENGTH * sizeof(char));
    err = fgets(buffer[i], MAX_LINELENGTH, fp_intco);
  }
  rewind(fp_intco);

  /* ...and overwite */
  for(i=0; i < count-1; i++) {
    fprintf(fp_intco, "%s", buffer[i]);
    free(buffer[i]);
  }
  free(buffer);

  // Print out coordinates to intco.dat
  fprintf(fp_intco,"  symm = ( \n");
  for (h=0; h<syminfo.nirreps; ++h) {
    if (h==1) fprintf(fp_intco, "  )\n  asymm = (\n");
    for (i=0;i<num_nonzero;++i) {
      if (irr[i] == h) {

        // print irrep label
        fprintf(fp_intco,"    (");
        fprintf(fp_intco,"\"%s\"",syminfo.clean_irrep_lbls[irr[i]]);

        // print simples ids
        fprintf(fp_intco," (");
        for (col=0, j=0;j<simples.get_num();++j) {
          if ( fabs(evectst_symm[i][j]) > 1.0E-10 ) {
            fprintf(fp_intco," %d",simples.index_to_id(j));
            if (col == 20) { fprintf(fp_intco,"\n    "); col = -1; }
            ++col;
          }
        }
        fprintf(fp_intco,")\n");

        // print coefficients of salc
        fprintf(fp_intco,"   (");
        for (col=0, j=0;j<simples.get_num();++j) {
          if ( fabs(evectst_symm[i][j]) > 1.0E-10 ) {
           fprintf(fp_intco,"%15.7lf",evectst_symm[i][j]);
          //fprintf(fp_intco,"%15.10lf",evectst_symm[i][j]);
          if (col == 7) { fprintf(fp_intco,"\n    "); col = -1; }
          ++col;
         }
        }
        fprintf(fp_intco,"))\n");
      }
    }
  }

  fprintf(fp_intco,"  )\n)\n"); 
  fflush(fp_intco);
  fclose(fp_intco);

  free_block(evectst_symm);
  return;
}

double repulsion(double *fatomic_num, double *geom);

// removes rotations - returns number of rotations removed
void rm_rotations(internals &simples, cartesians &carts, int &num_nonzero,
    double **evects) {

  int i, j, k, a, b, c, d, ivect, cnt;
  double *disp_coord, scale=0.001, *fatomic_num, energy, disp_energy;
  double *coord, rot_tol = 1.0E-10;

  fatomic_num = new double[optinfo.natom];
  for (i=0;i<optinfo.natom;++i)
    fatomic_num[i] = carts.get_fatomic_num(i);

  coord = carts.get_coord();

  energy = repulsion(fatomic_num,coord);

  disp_coord = new double [3*optinfo.natom];
  for (ivect=0; ivect<num_nonzero; ++ivect) {

    for (i=0;i<(3*optinfo.natom);++i)
      disp_coord[i] = coord[i];

    cnt = -1;
    for (i=0;i<simples.stre.get_num();++i) {
      a = simples.stre.get_A(i);
      b = simples.stre.get_B(i);
      ++cnt;
      for (k=0;k<3;++k) {
        disp_coord[3*a+k] += scale * evects[cnt][ivect] * simples.stre.get_s_A(i,k);
        disp_coord[3*b+k] += scale * evects[cnt][ivect] * simples.stre.get_s_B(i,k);
      }
    }
    for (i=0;i<simples.bend.get_num();++i) {
      a = simples.bend.get_A(i);
      b = simples.bend.get_B(i);
      c = simples.bend.get_C(i);
      ++cnt;
      for (k=0;k<3;++k) {
        disp_coord[3*a+k] += scale * evects[cnt][ivect] * simples.bend.get_s_A(i,k);
        disp_coord[3*b+k] += scale * evects[cnt][ivect] * simples.bend.get_s_B(i,k);
        disp_coord[3*c+k] += scale * evects[cnt][ivect] * simples.bend.get_s_C(i,k);
      }
    }
    for (i=0;i<simples.tors.get_num();++i) {
      a = simples.tors.get_A(i);
      b = simples.tors.get_B(i);
      c = simples.tors.get_C(i);
      d = simples.tors.get_D(i);
      ++cnt;
      for (k=0;k<3;++k) {
        disp_coord[3*a+k] += scale * evects[cnt][ivect] * simples.tors.get_s_A(i,k);
        disp_coord[3*b+k] += scale * evects[cnt][ivect] * simples.tors.get_s_B(i,k);
        disp_coord[3*c+k] += scale * evects[cnt][ivect] * simples.tors.get_s_C(i,k);
        disp_coord[3*d+k] += scale * evects[cnt][ivect] * simples.tors.get_s_D(i,k);
      }
    }
    for (i=0;i<simples.out.get_num();++i) {
      a = simples.out.get_A(i);
      b = simples.out.get_B(i);
      c = simples.out.get_C(i);
      d = simples.out.get_D(i);
      ++cnt;
      for (k=0;k<3;++k) {
        disp_coord[3*a+k] += scale * evects[cnt][ivect] * simples.out.get_A(i);
        disp_coord[3*b+k] += scale * evects[cnt][ivect] * simples.out.get_B(i);
        disp_coord[3*c+k] += scale * evects[cnt][ivect] * simples.out.get_C(i);
        disp_coord[3*d+k] += scale * evects[cnt][ivect] * simples.out.get_D(i);
      }
    }
    disp_energy = repulsion(fatomic_num, disp_coord);
    // fprintf(outfile,"dispenergy - energy: %16.12lf\n",disp_energy - energy);
    if (fabs(disp_energy - energy) < rot_tol) {
      fprintf(outfile,"rotational coordinate eliminated");
      for (i=ivect+1; i<num_nonzero; ++i) 
        for (j=0; j<simples.get_num(); ++j) 
          evects[i-1][j] = evects[i][j];

      --ivect;
      --num_nonzero;
    }
  }
  free(coord);
  delete [] disp_coord;
  return;
}

double repulsion(double *fatomic_num, double *coord) {
  int i, j, dim;
  double dist, tval = 0.0;

  dim = optinfo.natom;
  for (i=0; i<dim; ++i)
    for (j=0; j<i; ++j) {
      dist = sqrt(
          SQR(coord[3*i+0]-coord[3*j+0])
          + SQR(coord[3*i+1]-coord[3*j+1])
          + SQR(coord[3*i+2]-coord[3*j+2]) );

      tval += fatomic_num[i]*fatomic_num[j] / dist;
    }
  fprintf(outfile,"returning repulsion: %15.10lf \n", tval);
  return tval;
}

