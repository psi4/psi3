
// This function forms a set of delocalized, symmetry-adapted
// internal coordinates, given a set of simples and their
// s-vectors.  The new coordinates will not mix stretches and
// bends unless MIX_TYPES is set to 1.

extern "C" {
  #include <stdio.h>
  #include <file30.h>
  #include <stdlib.h>
  #include <string.h>
  #include <physconst.h>
  #include <math.h>
  #include <libciomr.h>
  #include <ip_libv1.h>
}

#define EXTERN
#include "opt.h"
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

extern double **irrep(internals &simples, double **salc);

void delocalize(int num_atoms,internals &simples) {
  int error,i,j,k,a,b,c,d,id,count,intco_type,sub_index,row[4],dim[4];
  int rotor_type, degrees_of_freedom;
  double **stre_mat, **bend_mat, **tors_mat, **out_mat;
  double *evals, **evects, **B, **BBt, **u, **temp_mat;

 // Build B matrix for simples
  B = init_matrix(simples.get_num(),num_atoms*3);
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
  //  print_mat2(B,simples.get_num(),num_atoms*3,outfile);

 // Form BBt matrix
  BBt = init_matrix(simples.get_num(),simples.get_num());
  if (optinfo.mix_types) {
    mmult(B,0,B,1,BBt,0,simples.get_num(),num_atoms*3,simples.get_num(),0);
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
        mmult(&(B[row[i]]),0,&(B[row[i]]),1,ptr,0,dim[i],num_atoms*3,dim[i],0);
    }
    free(ptr);
  }
  //  fprintf(outfile,"The BB^t Matrix:\n");
  //  print_mat2(BBt,simples.get_num(),simples.get_num(),outfile);

 // Diagonalize BBt
  evals = init_array(simples.get_num());
  evects = init_matrix(simples.get_num(),simples.get_num());
  sq_rsp(simples.get_num(),simples.get_num(), BBt, evals, 1, evects, 1.0E-14);

  /* check eigenvectors of BBt */
  if (optinfo.print_delocalize) {
     fprintf(outfile,"\n\n+++Delocalized Coordinate Formation+++\n");
     fprintf(outfile,"\n\nChecking eigenvectors of BBt...\n");
     eivout(evects, evals, simples.get_num(), simples.get_num(), outfile);
   }
  temp_mat = init_matrix(simples.get_num(),simples.get_num());
  mmult(BBt,0,evects,0,temp_mat,0,simples.get_num(),simples.get_num(),simples.get_num(),0);
  for (j=0;j<simples.get_num();++j) {
     error = 0;
     for (i=0;i<simples.get_num();++i) {
        if ( fabs(temp_mat[i][j] - evals[j]*evects[i][j]) > 1.0E-13)  error = 1;
        if (error == 1) { fprintf(outfile,"Error in evect %d\n",j); error = 0;}
       }
    }

    fflush(outfile);
    free_matrix(temp_mat,simples.get_num());
  

  
  /* check for proper number of non-redundant coordinates (3n-6) */ 
  num_nonzero = 0;
  for (i=0;i<simples.get_num();++i)
    if( evals[i] > optinfo.ev_tol ) ++num_nonzero;

  rewind(fp_input);
  ip_initialize(fp_input,outfile);
  file30_init();
  rotor_type = file30_rd_rottype();
  file30_close();
  ip_done();
  switch (rotor_type) {
  case 3:
      degrees_of_freedom = 3 * num_atoms - 5;
      break;
  case 6:
      degrees_of_freedom = 0;
      break;
  default:
      degrees_of_freedom = 3 * num_atoms - 6;
      break;
  }
  if (num_nonzero < degrees_of_freedom) {
    if (num_nonzero < degrees_of_freedom) {
       punt("too few non-redundant coordinates, try reducing eigenvalue tolerance\n\tand check atom connectivity");
      }
    if (num_nonzero > degrees_of_freedom) {
       punt("too many non-redundant coordinates, try increasing eigenvalue tolerance");
      }
  }


  /* transpose evects matrix, also throw out redundant coordinates
     (eigenvectors corresponding to zero eigenvalues*/

  double **evectst;
  char buffer[MAX_LINELENGTH], *err;
  int h;
  irr = init_int_array(simples.get_num());

  evectst = init_matrix(num_nonzero,simples.get_num());
  
  for (i=0;i<simples.get_num();++i) {
      for (j=simples.get_num()-num_nonzero;j<simples.get_num();++j) {
	  evectst[j-simples.get_num()+num_nonzero][i] = evects[i][j];
	}
    }
  if (optinfo.print_delocalize == 1) {
    fprintf(outfile,"\n\nTransposed non-redundant delocalized internal coordinate matrix (each row is a coordinate):\n");
    print_mat(evectst,num_nonzero,simples.get_num(),outfile);
  }
  free_matrix(evects,simples.get_num());



  fp_intco = fopen("intco.dat", "r+");
  count = 0;
  for( ; ; ) {
    err = fgets(buffer, MAX_LINELENGTH, fp_intco);
    if (err == NULL) break;
    ++count;
  }
  rewind(fp_intco);
  for(i=0;i<(count-1);++i)
    err = fgets(buffer, MAX_LINELENGTH, fp_intco);
  fflush(fp_intco);



  /* send evectst to irrep.cc to be properly symmetrized */
    evectst = irrep(simples, evectst);

double** pointy = evectst;


 // Print out coordinates to intco.dat
  int col;
  for (h=0;h<num_nonzero;++h) {
    if (h==0) fprintf(fp_intco,"  symm = ( \n");
    for (i=0;i<num_nonzero;++i) {
      if (h == irr[i]) {
        fprintf(fp_intco,"    (");
        fprintf(fp_intco,"\"%s\"",syminfo.clean_irrep_lbls[irr[i]]);
        fprintf(fp_intco," (");
        for (col=0, j=0;j<simples.get_num();++j) {
          if ( fabs(evectst[i][j]) > 1.0E-10 ) {
            fprintf(fp_intco," %d",simples.index_to_id(j));
            if (col == 20) { fprintf(fp_intco,"\n    "); col = -1; }
            ++col;
          }
        }
        fprintf(fp_intco,")\n");
        fprintf(fp_intco,"   (");
        for (col=0, j=0;j<simples.get_num();++j) {
          if ( fabs(evectst[i][j]) > 1.0E-10 ) {
            fprintf(fp_intco,"%10.6lf",evectst[i][j]);
            if (col == 6) { fprintf(fp_intco,"\n    "); col = -1; }
            ++col;
          }
        }
        fprintf(fp_intco,"))\n");
      }
    }
    if (h==0) fprintf(fp_intco, "  )\n  asymm = (\n");
  }

  fflush(fp_intco);

  free_matrix(B,simples.get_num());
  free_matrix(BBt,simples.get_num());
  free(evals);
  free_matrix(pointy,num_nonzero); 

  fprintf(fp_intco,"  )\n)\n"); 
  fflush(fp_intco);
  fclose(fp_intco);
                          
  return;
}

