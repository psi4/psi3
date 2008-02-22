/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include "blas.h"
#include "debugging.h"
#include "memory_manager.h"
#include "moinfo.h"
#include "utilities.h"
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <cstdlib>

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCBLAS::add_index(char* cstr)
{
  // Make sure that the element that we are adding is not present
  string str(cstr);
  to_lower(str);
  if(indices.find(str)==indices.end()){
    indices.insert(make_pair(str,new CCIndex(str)));
  }
}

void CCBLAS::add_Matrix(char* cstr)
{
  string str(cstr);
  vector<string> names = moinfo->get_matrix_names(str);
  for(int n=0;n<names.size();n++)
    add_Matrix_ref(names[n]);
}

void CCBLAS::add_Matrix(string str)
{
  vector<string> names = moinfo->get_matrix_names(str);
  for(int n=0;n<names.size();n++)
    add_Matrix_ref(names[n]);
}

void CCBLAS::add_Matrix_ref(std::string& str)
{
  // Make sure that the element that we are adding is not present
  MatrixMap::iterator iter = matrices.find(str);
  if(iter==matrices.end()){
    CCIndex* index_pointer[2];
    // Default: assume the [] indexing
    index_pointer[0]=get_index("[]");
    index_pointer[1]=get_index("[]");
    vector<string> index_string_vec = split_indices(str);
    for(int i=0;i<index_string_vec.size();++i)
      index_pointer[i]=get_index(index_string_vec[i]);
    matrices.insert(make_pair(str,new CCMatrix(str,index_pointer[0],index_pointer[1])));
  }
}

CCIndex* CCBLAS::get_index(char* cstr)
{
  string str(cstr);
  to_lower(str);
  // Make sure that the element that we are retrieving is present
  IndexMap::iterator iter = indices.find(str);
  if(iter!=indices.end()){
    return(indices[str]);
  }
  string err("\nCCBLAS::get_index() couldn't find index " + str);
  print_error(err,__FILE__,__LINE__);
  return(NULL);
}

CCIndex* CCBLAS::get_index(string& str)
{
  to_lower(str);
  // Make sure that the element that we are retrieving is present
  IndexMap::iterator iter = indices.find(str);
  if(iter!=indices.end()){
    return(indices[str]);
  }
  string err("\nCCBLAS::get_index() couldn't find index " + str);
  print_error(err,__FILE__,__LINE__);
  return(NULL);
}

CCMatTmp CCBLAS::get_MatTmp(std::string str, int reference, DiskOpt disk_option)
{
  append_reference(str,reference);
  load(get_Matrix(str));
  return(CCMatTmp(get_Matrix(str),disk_option)); 
}

CCMatTmp CCBLAS::get_MatTmp(std::string str, DiskOpt disk_option)
{
  load(get_Matrix(str));
  return(CCMatTmp(get_Matrix(str),disk_option)); 
}

CCMatTmp CCBLAS::get_MatTmp(CCMatrix* Matrix, DiskOpt disk_option)
{
  load(Matrix);
  return(CCMatTmp(Matrix,disk_option)); 
}

CCMatIrTmp CCBLAS::get_MatIrTmp(std::string str, int reference, int irrep,  DiskOpt disk_option)
{
  append_reference(str,reference);
  load_irrep(get_Matrix(str),irrep);
  return(CCMatIrTmp(get_Matrix(str),irrep,disk_option)); 
}

CCMatIrTmp CCBLAS::get_MatIrTmp(std::string str, int irrep, DiskOpt disk_option)
{
  load_irrep(get_Matrix(str),irrep);
  return(CCMatIrTmp(get_Matrix(str),irrep,disk_option));
}

CCMatIrTmp CCBLAS::get_MatIrTmp(CCMatrix* Matrix, int irrep, DiskOpt disk_option)
{
  load_irrep(Matrix,irrep);
  return(CCMatIrTmp(Matrix,irrep,disk_option));
}

CCMatrix* CCBLAS::get_Matrix(char* cstr, int reference)
{
  string str(cstr);
  append_reference(str,reference);
  return(get_Matrix(str));
}

CCMatrix* CCBLAS::get_Matrix(char* cstr)
{
  string str(cstr);
  return(get_Matrix(str));
}

CCMatrix* CCBLAS::get_Matrix(string& str)
{
  // Make sure that the element that we are retrieving is present
  MatrixMap::iterator iter = matrices.find(str);
  if(iter!=matrices.end())
    return(matrices[str]);
  string err("\nCCBLAS::get_matrix() couldn't find matrix " + str);
  print_error(err,__FILE__,__LINE__);
  return(NULL);
}

CCMatrix* CCBLAS::get_Matrix(string& str, string& expression)
{
  // Make sure that the element that we are retrieving is present
  MatrixMap::iterator iter = matrices.find(str);
  if(iter!=matrices.end()){
    return(matrices[str]);
  }
  fprintf(outfile,"\n\nCCBLAS::parse() couldn't find matrix %s in the CCMatrix list\n\nwhile parsing the string:\n\t%s\n\n",str.c_str(),expression.c_str());
  fflush(outfile);
  exit(1);
}

void CCBLAS::set_scalar(char* cstr,int reference,double value)
{
  string str(cstr);
  set_scalar(str,reference,value);
}

void CCBLAS::set_scalar(string& str,int reference,double value)
{
  string matrix_str = add_reference(str,reference);
  // Make sure that the element that we are retrieving is present
  MatrixMap::iterator iter = matrices.find(matrix_str);
  if(iter!=matrices.end()){
    load(iter->second);
    iter->second->set_scalar(value);
    return;
  }
  string err("\nCCBLAS::set_scalar() couldn't find matrix " + matrix_str);
  print_error(err.c_str(),__FILE__,__LINE__);
}

double CCBLAS::get_scalar(char* cstr,int reference)
{
  string str(cstr);
  return(get_scalar(str,reference));
}

double CCBLAS::get_scalar(string& str,int reference)
{
  string matrix_str(str);
  append_reference(matrix_str,reference);
  // Make sure that the element that we are retrieving is present
  MatrixMap::iterator iter = matrices.find(matrix_str);
  if(iter!=matrices.end()){
    load(iter->second);
    return(iter->second->get_scalar());
  }
  string err("\nCCBLAS::get_scalar() couldn't find matrix " + matrix_str);
  print_error(err.c_str(),__FILE__,__LINE__);
  return (0.0);
}

double CCBLAS::get_scalar(string str)
{
  // Make sure that the element that we are retrieving is present
  MatrixMap::iterator iter = matrices.find(str);
  if(iter!=matrices.end()){
    load(iter->second);
    return(iter->second->get_scalar());
  }
  string err("\nCCBLAS::get_scalar() couldn't find matrix " + str);
  print_error(err.c_str(),__FILE__,__LINE__);
  return (0.0);
}

void CCBLAS::load(CCMatrix* Matrix)
{
  if(Matrix->is_allocated()){
    DEBUGGING(2,
      fprintf(outfile,"\nCCBLAS::load(%s): matrix is in core.",Matrix->get_label().c_str());
    );
  }else{
    DEBUGGING(2,
      fprintf(outfile,"\nCCBLAS::load(%s): matrix is not in core. Loading it :[",Matrix->get_label().c_str());
    );
    // Do we have enough memory to fit the entire matrix in core?
    double memory_required = Matrix->get_memory();
    make_space(memory_required);
    Matrix->load();
    DEBUGGING(2,
      fprintf(outfile,"\n] <- done.");
    );
  }
}

void CCBLAS::load_irrep(CCMatrix* Matrix,int h)
{
  if(Matrix->is_block_allocated(h)){
    DEBUGGING(2,
      fprintf(outfile,"\nCCBLAS::load_irrep(%s,%d): matrix block is in core.",Matrix->get_label().c_str(),h);
    )
  }else{
    DEBUGGING(2,
      fprintf(outfile,"\nCCBLAS::load_irrep(%s,%d): matrix block is not in core. Loading it : [",Matrix->get_label().c_str(),h);
    )
    // Do we have enough memory to fit the entire matrix in core?
    double memory_required = Matrix->get_memorypi(h);
    make_space(memory_required);
    Matrix->load_irrep(h);
    DEBUGGING(2,
      fprintf(outfile,"\n] <- done.");
    )
  }
}

void CCBLAS::make_space(double memory_required)
{
  if(memory_required < mem->get_free_memory())
    return;
  else{
    fprintf(outfile,"\nCCBLAS::make_space() not implemented yet!!!");
    // Attempt #1
    

  }
}

// double*** CCBLAS::get_sortmap(CCIndex* T_left,CCIndex* T_right,int thread)
// {
//   string sortstr = T_left->get_label() + T_right->get_label() + to_string(thread);
//   double*** T_matrix=NULL;
//   SortMap::iterator iter = sortmap.find(sortstr);
//   if(iter==sortmap.end()){
//     int T_matrix_offset = 0;
//     T_matrix = new double**[moinfo->get_nirreps()];
//     for(int irrep=0;irrep<moinfo->get_nirreps();irrep++){
//       T_matrix[irrep] = new double*[T_left->get_pairpi(irrep)];
//       for(int i=0;i<T_left->get_pairpi(irrep);i++){
//         T_matrix[irrep][i]=&(work[thread][T_matrix_offset+i*T_right->get_pairpi(irrep)]);
//       }
//       zero_arr(&(T_matrix[irrep][0][0]),T_left->get_pairpi(irrep)*T_right->get_pairpi(irrep));
//       T_matrix_offset+=T_left->get_pairpi(irrep)*T_right->get_pairpi(irrep);
//     }
//     sortmap.insert(make_pair(sortstr,T_matrix));
//   }else{
//     T_matrix = sortmap[sortstr];
//   }
//   return(T_matrix);
// }

}} /* End Namespaces */
