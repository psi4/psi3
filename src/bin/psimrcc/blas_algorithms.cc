/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include "blas.h"
#include "debugging.h"
#include "error.h"
#include "moinfo.h"
#include "memory_manager.h"
#include "utilities.h"

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCBLAS::zero_right_four_diagonal(char* cstr)
{
  string str(cstr);
  // To zero diagonals of things like "Fae[v][v]{u}"
  vector<string> names = moinfo->get_matrix_names(str);
  for(int n=0;n<names.size();n++){
    CCMatrix* Matrix = get_Matrix(names[n]);
    Matrix->zero_right_four_diagonal();
    DEBUGGING(5,
      fprintf(outfile,"\n...setting the right diagonal terms of %s to zero",names[n].c_str());
    );
  }
}

void CCBLAS::zero_non_doubly_occupied(char* cstr)
{
  string str(cstr);
  // To zero non-doubly occupied MOs of things like "Fae[v][v]{u}"
  vector<string> names = moinfo->get_matrix_names(str);
  for(int n=0;n<names.size();n++){
    CCMatrix* Matrix = get_Matrix(names[n]);
    Matrix->zero_non_doubly_occupied();
    DEBUGGING(5,
      fprintf(outfile,"\n...setting the right diagonal terms of %s to zero",names[n].c_str());
    );
  }
}

void CCBLAS::zero_non_external(char* cstr)
{
  string str(cstr);
  // To zero non-external MOs of things like "Fae[v][v]{u}"
  vector<string> names = moinfo->get_matrix_names(str);
  for(int n=0;n<names.size();n++){
    CCMatrix* Matrix = get_Matrix(names[n]);
    Matrix->zero_non_external();
    DEBUGGING(5,
      fprintf(outfile,"\n...setting the right diagonal terms of %s to zero",names[n].c_str());
    );
  }
}

void CCBLAS::reduce_spaces(char* out,char* in)
{
  string  in_str(in);
  string out_str(out);
  // To zero diagonals of things like "Fae[v][v]{u}"
  vector<string>  in_names = moinfo->get_matrix_names(in_str);
  vector<string> out_names = moinfo->get_matrix_names(out_str);
  if(in_names.size()!=out_names.size())
    print_error("CCBLAS::map_spaces, number of references mismatch",__FILE__,__LINE__);
  for(int n=0;n<in_names.size();n++){
    CCMatrix*  in_Matrix = get_Matrix(in_names[n]);
    CCMatrix* out_Matrix = get_Matrix(out_names[n]);
    process_reduce_spaces(out_Matrix,in_Matrix);
  }
}

void CCBLAS::process_reduce_spaces(CCMatrix* out_Matrix,CCMatrix* in_Matrix)
{
  double*** out_matrix = out_Matrix->get_matrix();
  int*      act_to_occ = moinfo->get_actv_to_occ();
  int*      act_to_vir = moinfo->get_actv_to_vir();

  string& out_index_label = out_Matrix->get_index_label();
  string&  in_index_label =  in_Matrix->get_index_label();

  int index_label_size = out_index_label.size();

  int** map;
  allocate2(int,map,index_label_size,moinfo->get_nmo());

  for(int k=0;k<index_label_size;k++){
    if(out_index_label[k]=='a' && in_index_label[k]=='o'){
      for(int l=0;l<moinfo->get_nactv();l++){
        map[k][l] = act_to_occ[l];
      }
    }else if(out_index_label[k]=='a' && in_index_label[k]=='v'){
      for(int l=0;l<moinfo->get_nactv();l++){
        map[k][l] = act_to_vir[l];
      }
    } else {
      for(int l=0;l<moinfo->get_nmo();l++){
        map[k][l] = l;
      }
    }
  }

  if(index_label_size==2){
    short* pq = new short[2];

    for(int h=0;h<moinfo->get_nirreps();h++){
      for(int i=0;i<out_Matrix->get_left_pairpi(h);i++){
        for(int j=0;j<out_Matrix->get_right_pairpi(h);j++){
          out_Matrix->get_two_indices(pq,h,i,j);
          out_matrix[h][i][j] = in_Matrix->get_two_address_element(map[0][pq[0]],map[1][pq[1]]);
        }
      }
    }

    delete[] pq;
  }else if(index_label_size==4){
    short* pqrs = new short[4];
    for(int h=0;h<moinfo->get_nirreps();h++){
      for(int i=0;i<out_Matrix->get_left_pairpi(h);i++){
        for(int j=0;j<out_Matrix->get_right_pairpi(h);j++){
          out_Matrix->get_four_indices(pqrs,h,i,j);
          out_matrix[h][i][j] = in_Matrix->get_four_address_element(map[0][pqrs[0]],map[1][pqrs[1]],map[2][pqrs[2]],map[3][pqrs[3]]);
        }
      }
    }
    delete[] pqrs;
  }
  release2(map);
}

}} /* End Namespaces */