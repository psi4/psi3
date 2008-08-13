/***************************************************************************
 *   Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *   frank@ccc.uga.edu
 *   SR/MRCC Code
 ***************************************************************************/
#if 0
#include <libmoinfo/libmoinfo.h>
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include <libutil/libutil.h>

namespace psi{ namespace psimrcc{

extern FILE* outfile;

using namespace std;

void CCMRCC::canonicalize()
{
  for(int n=0;n<moinfo->get_nrefs();n++){
    string n_str  = to_string(moinfo->get_ref_number("a",n));
    string factor = to_string(pow(eigenvector[moinfo->get_ref_number("a",n)],2.0));
    solve("avg_fock[o][o] += " + factor + " fock[o][o]{" + n_str + "}");
    fprintf(outfile,"\n\nFactor = %s ",factor.c_str());
    blas->print("fock[o][o]{" + n_str + "}");
  }
  blas->print("avg_fock[o][o]");

  double** scf = moinfo->get_scf_mos();

  for(int i=0;i<moinfo->get_nirreps();i++){
    int matrix_size = moinfo->get_actv()[i];
    if(matrix_size>1){
      int matrix_offset = moinfo->get_docc()[i];
      double*  evalues = new double[matrix_size];
      double** avg_fock = CCMatrix::get_matrix("avg_fock[o][o]")[i];
      double** matrix;
      double** evectors;
      double** newMOs;
      init_dmatrix(matrix,matrix_size,matrix_size);
      init_dmatrix(evectors,matrix_size,matrix_size);
      init_dmatrix(newMOs,moinfo->get_nso(),matrix_size);

      for(int j=0;j<matrix_size;j++)
        for(int k=0;k<matrix_size;k++){
          matrix[j][k]=avg_fock[j+matrix_offset][k+matrix_offset];
        }
      fprintf(outfile,"\n\nIrrep = %s",moinfo->get_irr_labs(i));
      print_dmatrix(matrix,matrix_size,matrix_size,outfile,"Averaged Fock Matrix");

      // TBM use blas!!
      sq_rsp(matrix_size,matrix_size,matrix,evalues,1,evectors,1.0E-14);
      print_dmatrix(evectors,matrix_size,matrix_size,outfile,"Eigenvectors");

      int index=0;
      for(int j=moinfo->get_first_occupied_pitzer(alpha)[i]+matrix_offset;j<moinfo->get_last_occupied_pitzer(alpha)[i];j++){
        fprintf(outfile,"\nYou need to mix orbital %d",j);
        for(int k=0;k<moinfo->get_nso();k++)
          newMOs[k][index]=scf[k][j];
        index++;
      }
      index=0;
      for(int j=moinfo->get_first_occupied_pitzer(alpha)[i]+matrix_offset;j<moinfo->get_last_occupied_pitzer(alpha)[i];j++){
        for(int k=0;k<moinfo->get_nso();k++){
          scf[k][j]=0.0;
          for(int l=0;l<matrix_size;l++)
            scf[k][j]+=newMOs[k][l]*evectors[l][index];
        }
        index++;
      }

      print_dmatrix(evectors,matrix_size,matrix_size,outfile,"Eigenvectors");

      free_dmatrix(matrix,matrix_size,matrix_size);
      free_dmatrix(evectors,matrix_size,matrix_size);
      free_dmatrix(newMOs,moinfo->get_nso(),matrix_size);
      delete[] evalues;
    }
  }
  moinfo->write_scf_mos();
}

}} /* End Namespaces */

#endif
