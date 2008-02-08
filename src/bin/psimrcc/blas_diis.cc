/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include "algebra_interface.h"
#include "calculation_options.h"
#include "blas.h"
#include "moinfo.h"
#include "memory_manager.h"
#include "utilities.h"
#include "error.h"
#include <limits>
#include <cmath>

#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

vector<pair<string,string> > diis_matrices;
const double diis_singular_tollerance = 1.0e-12;

void CCBLAS::diis_add(string amps, string delta_amps)
{
  vector<string> amps_names = moinfo->get_matrix_names(amps);
  vector<string> delta_amps_names = moinfo->get_matrix_names(delta_amps);
  for(int n=0;n<amps_names.size();n++){
    diis_matrices.push_back(make_pair(amps_names[n],delta_amps_names[n]));
  }
}

void CCBLAS::diis_save_t_amps(int cycle)
{
  int diis_step = cycle % options->get_int_option("MAXDIIS");
  for(vector<pair<string,string> >::iterator it=diis_matrices.begin();it!=diis_matrices.end();++it){
    for(int h=0;h<moinfo->get_nirreps();h++){
      CCMatIrTmp Amps = get_MatIrTmp(it->first,h,none);
      double** matrix = Amps->get_matrix()[h];
      size_t   block_sizepi = Amps->get_block_sizepi(h);
      if(block_sizepi>0){
        char data_label[80];
        sprintf(data_label,"%s_%s_%d_%d",(it->first).c_str(),"DIIS",h,diis_step);
        psio_write_entry(MRCC_ON_DISK,data_label,(char*)&(matrix[0][0]),block_sizepi*sizeof(double));
      }
    }
  }
}

void CCBLAS::diis(int cycle)
{
  int diis_step = cycle % options->get_int_option("MAXDIIS");

  for(vector<pair<string,string> >::iterator it=diis_matrices.begin();it!=diis_matrices.end();++it){
    if(it->second.find("t3_delta")==string::npos){
      for(int h=0;h<moinfo->get_nirreps();h++){
        CCMatIrTmp DeltaAmps = get_MatIrTmp(it->second,h,none);
        double** matrix = DeltaAmps->get_matrix()[h];
        size_t   block_sizepi = DeltaAmps->get_block_sizepi(h);
        if(block_sizepi>0){
          char data_label[80];
          sprintf(data_label,"%s_%s_%d_%d",(it->second).c_str(),"DIIS",h,diis_step);
          psio_write_entry(MRCC_ON_DISK,data_label,(char*)&(matrix[0][0]),block_sizepi*sizeof(double));
        }
      }
    }
  }

  fprintf(outfile,"   S");

  // Do a DIIS step 
  if(diis_step==options->get_int_option("MAXDIIS")-1){
    double** diis_B;
    double*  diis_A;
    allocate1(double,diis_A,options->get_int_option("MAXDIIS")+1);
    allocate2(double,diis_B,options->get_int_option("MAXDIIS")+1,options->get_int_option("MAXDIIS")+1);
    bool singularities_found = false;
    for(vector<pair<string,string> >::iterator it=diis_matrices.begin();it!=diis_matrices.end();++it){
      // Zero A and B
      for(int i=0;i<options->get_int_option("MAXDIIS");i++){
        diis_A[i]=0.0;
        diis_B[i][options->get_int_option("MAXDIIS")]=diis_B[options->get_int_option("MAXDIIS")][i]=-1.0;
        for(int j=0;j<options->get_int_option("MAXDIIS");j++)
          diis_B[i][j]=0.0;
      }
      diis_B[options->get_int_option("MAXDIIS")][options->get_int_option("MAXDIIS")]=0.0;
      diis_A[options->get_int_option("MAXDIIS")]=-1.0;

      // Build B
      for(int h=0;h<moinfo->get_nirreps();h++){
        CCMatIrTmp Amps       = get_MatIrTmp(it->first,h,none);
        size_t   block_sizepi = Amps->get_block_sizepi(h);
        if(block_sizepi>0){
          double*  i_matrix;
          double*  j_matrix;
          allocate1(double,i_matrix,block_sizepi);
          allocate1(double,j_matrix,block_sizepi);

          // Build the diis_B matrix
          for(int i=0;i<options->get_int_option("MAXDIIS");i++){
            // Load vector i irrep h
            char i_data_label[80];
            sprintf(i_data_label,"%s_%s_%d_%d",(it->second).c_str(),"DIIS",h,i);
            psio_read_entry(MRCC_ON_DISK,i_data_label,(char*)&(i_matrix[0]),block_sizepi*sizeof(double));

            for(int j=i;j<options->get_int_option("MAXDIIS");j++){
              // Load vector j irrep h
              char j_data_label[80];
              sprintf(j_data_label,"%s_%s_%d_%d",(it->second).c_str(),"DIIS",h,j);
              psio_read_entry(MRCC_ON_DISK,j_data_label,(char*)&(j_matrix[0]),block_sizepi*sizeof(double));

              int dx = 1;
              int lenght = block_sizepi;
              if( block_sizepi < static_cast<size_t>(numeric_limits<int>::max()) ){
                diis_B[i][j] += F_DDOT(&lenght,i_matrix,&dx,j_matrix,&dx);
                diis_B[j][i] = diis_B[i][j];
              }else{
                print_error("The numeric limits for int was reached for F_DDOT",__FILE__,__LINE__);
              }
            }
          }
          release1(i_matrix);
          release1(j_matrix);
        }
      }

      // Solve B x = A
      int  matrix_size = options->get_int_option("MAXDIIS") + 1;          
      int* IPIV = new int[matrix_size];
      int nrhs = 1;
      int info = 0;
      F_DGESV(&matrix_size, &nrhs, &(diis_B[0][0]),&matrix_size, &(IPIV[0]), &(diis_A[0]),&matrix_size, &info);         
      delete[] IPIV;

      // Update T = sum t(i) * A(i);
      if(!info){
        for(int h=0;h<moinfo->get_nirreps();h++){
          CCMatIrTmp Amps       = get_MatIrTmp(it->first,h,none);
          size_t   block_sizepi = Amps->get_block_sizepi(h);
          if(block_sizepi>0){
            // Update the amplitudes
            double*  i_matrix;
            double*  j_matrix;
            allocate1(double,i_matrix,block_sizepi);
            allocate1(double,j_matrix,block_sizepi);
            double* t_matrix = &(Amps->get_matrix()[h][0][0]);
            Amps->zero_matrix_block(h);
            for(int i=0;i<options->get_int_option("MAXDIIS");i++){
              char i_data_label[80];
              sprintf(i_data_label,"%s_%s_%d_%d",(it->first).c_str(),"DIIS",h,i);
              psio_read_entry(MRCC_ON_DISK,i_data_label,(char*)&(i_matrix[0]),block_sizepi*sizeof(double));
              for(size_t n=0;n<block_sizepi;n++){
                t_matrix[n] += diis_A[i]*i_matrix[n];
              }
            }
            release1(i_matrix);
            release1(j_matrix);
          }
        }
      }else{
        singularities_found = true;
      }

    }
    fprintf(outfile,"/E");
    if(singularities_found)
      fprintf(outfile," (singularities found)");
    release1(diis_A);
    release2(diis_B);
  }
}


/*
SOME OLD CODE THAT WAS SUPPOSED TO BE FANCIER, IT DOESN'T REALLY IMPROVE

//           // This is a stable matrix inversion routine
// 

//           double** eigenvectors;
//           double** original_matrix;
//           double*  eigenvalues;
// 
//           int  matrix_size = options->get_int_option("MAXDIIS") + 1;          
//           init_matrix<double>(eigenvalues,matrix_size);
//           init_matrix<double>(eigenvectors,matrix_size,matrix_size);
//           init_matrix<double>(original_matrix,matrix_size,matrix_size);
// 
//           for(int i = 0; i < matrix_size;i++)
//             for(int j = 0; j < matrix_size;j++)
//               original_matrix[i][j] = diis_B[i][j];
// 
//           sq_rsp(matrix_size,matrix_size,diis_B,eigenvalues,1,eigenvectors,1e-14);
// 
//           fprintf(outfile,"\n DIIS %s[%d]",it->first.c_str(),h);
//           for(int i=0;i<matrix_size;i++)
//             fprintf(outfile,"\n eigenvalues[%d] = %20.12f",i,eigenvalues[i]);
// 
//           for(int k = 0; k < matrix_size;k++){
//             if(fabs(eigenvalues[k]) < diis_singular_tollerance){
//               eigenvalues[k]=0.0;
//             }else{
//               eigenvalues[k]=1/eigenvalues[k];
//             }
//           }
//           for(int i=0;i<matrix_size;i++)
//             fprintf(outfile,"\n eigenvalues[%d] = %20.12f",i,eigenvalues[i]);
//           
//           for(int i = 0; i < matrix_size;i++){
//             for(int j = 0; j < matrix_size;j++){
//               diis_B[i][j] = 0.0;
//               for(int k = 0; k < matrix_size;k++){
//                 diis_B[i][j] += eigenvectors[i][k]*eigenvectors[j][k]*eigenvalues[k];
//               }
//             }
//           }
//           print_dmatrix(diis_B,matrix_size,matrix_size,outfile,"diis_B");
// 
//           for(int i = 0; i < matrix_size;i++){
//             for(int j = 0; j < matrix_size;j++){
//               eigenvectors[i][j] = 0.0;
//               for(int k = 0; k < matrix_size;k++){
//                 eigenvectors[i][j] += original_matrix[i][k]*diis_B[k][j];
//               }
//             }
//           }
//           fprintf(outfile,"\n DIIS %s[%d] Testing inverse",it->first.c_str(),h);
//           print_dmatrix(eigenvectors,matrix_size,matrix_size,outfile);
// 
//           for(int i = 0; i < matrix_size;i++){
//             eigenvalues[i] = 0.0;
//             for(int j = 0; j < matrix_size;j++){
//               eigenvalues[i] += diis_B[i][j] * diis_A[j];
//             }
//           }
//           for(int i = 0; i < matrix_size;i++)
//             diis_A[i] = eigenvalues[i];
//           for(int i=0;i<matrix_size;i++)
//             fprintf(outfile,"\n A[%d] = %20.12f",i,diis_A[i]);
// 
//           free_matrix<double>(eigenvalues,matrix_size);
//           free_matrix<double>(eigenvectors,matrix_size,matrix_size);
//           free_matrix<double>(original_matrix,matrix_size,matrix_size);
// 
// 
//           fprintf(outfile,"\n DIIS %s[%d]",it->first.c_str(),h);
//           print_dmatrix(diis_B,matrix_size,matrix_size,outfile);
//           for(int i=0;i<matrix_size;i++){
//             fprintf(outfile,"\n A[%d] = %20.12f",i,diis_A[i]);
//           }

*/

}} /* End Namespaces */