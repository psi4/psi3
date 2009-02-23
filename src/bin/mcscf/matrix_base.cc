#include <cstring>
#include <iostream>

#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libutil/libutil.h>

#include "matrix_base.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

MatrixBase::MatrixBase(int rows, int cols) : rows_(rows),cols_(cols),elements_(rows*cols),matrix_(NULL)
{
  allocate2(double,matrix_,rows_,cols_);
}

MatrixBase::~MatrixBase()
{
  release2(matrix_);
}

void MatrixBase::print()
{
  for(int i=0 ; i < rows_; i++){
    fprintf(outfile,"\n  ");
    for(int j=0 ; j < cols_; j++)
      fprintf(outfile,"%10.6f",matrix_[i][j]);
  }
  fprintf(outfile,"\n");
}

void MatrixBase::scale(double factor)
{
  if(elements_>0)
    C_DSCAL(elements_,
            factor,
            &(matrix_[0][0]),
            1);
}

void MatrixBase::transpose()
{
  if(elements_>0){
    double temp;
    for(size_t i = 0; i < rows_; ++i){
      for(size_t j = i + 1; j < cols_; ++j){
        temp          = matrix_[i][j];
        matrix_[i][j] = matrix_[j][i];
        matrix_[j][i] = temp;
      }
    }
  }
}

void MatrixBase::zero()
{
  if(elements_>0)
    memset(&(matrix_[0][0]),'\0', sizeof(double) * elements_);
}

void MatrixBase::zero_diagonal()
{
  if(elements_>0 && (rows_ == cols_))
    for(int i=0 ; i < rows_; i++)
      matrix_[i][i] = 0.0;
}


void MatrixBase::multiply(bool transpose_A, bool transpose_B, MatrixBase* A, MatrixBase* B)
{
  char transa = transpose_A ? 't' : 'n';
  char transb = transpose_B ? 't' : 'n';
  if(elements_>0){
    // Multiply A and B
    int m = rows_;       // TODO This only works for square matrices!
    int n = rows_;
    int k = rows_;
    int nca = rows_;
    int ncb = rows_;
    int ncc = rows_;
    C_DGEMM(transa,
            transb,
            m,
            n,
            k,
            1.0,
            A->get_matrix()[0],nca,
            B->get_matrix()[0],ncb,
            0.0,
            get_matrix()[0],ncb);
  }
}

void MatrixBase::diagonalize(MatrixBase* eigenmatrix, VectorBase* eigenvalues)
{
  // Diagonalize the block
  if(elements_>0 && (rows_ == cols_)){
    sq_rsp(rows_,
           cols_,
           matrix_,
           eigenvalues->get_vector(),
           1,
           eigenmatrix->get_matrix(),
           1.0e-14);
  }
}

double dot(MatrixBase* A, MatrixBase* B)
{
  double value = 0.0;
  if(A->rows_ * A->cols_>0){
    for(int i = 0; i < A->rows_; ++i)
      for(int j = 0; j < A->cols_; ++j)
        value += A->matrix_[i][j] * B->matrix_[i][j];
  }
  return(value);
}

MatrixBase& MatrixBase::operator+=(const MatrixBase& rhs)
{
  if(elements_>0){
    for(size_t i = 0; i < rows_; ++i)
      for(size_t j = 0; j < cols_; ++j)
        matrix_[i][j] += rhs.matrix_[i][j];
  }
  return(*this);
}

}}

// void BlockMatrix::diagonalize(BlockMatrix* eigenvectors, double* eigenvalues)
// {
//   for(int h=0; h < nirreps; ++h){
//     if(block_size[h]>0){

//   }
// }
/*
void BlockMatrix::minus(BlockMatrix* B)
{
  for(int h=0; h < nirreps; ++h){
    double** A_matrix_block = matrix[h];
    double** B_matrix_block = B->get_block(h);
    if(block_size[h]>0){
      for(int i = 0; i < block_size[h]; ++i)
        for(int j = 0; j < block_size[h]; ++j)
          A_matrix_block[i][j] -= B_matrix_block[i][j];
    }
  }
}

*/

/*
BlockMatrix& BlockMatrix::operator=(const BlockMatrix& rhs)
{
  if(this == &rhs){
    return(*this);
  }
  for(int h=0; h < nirreps; ++h){
    double** lhs_matrix_block = matrix[h];
    const double** rhs_matrix_block = rhs.get_block(h);
    if(block_size[h]>0){
      for(int i = 0; i < block_size[h]; ++i)
        for(int j = 0; j < block_size[h]; ++j)
          lhs_matrix_block[i][j] = rhs_matrix_block[i][j];
    }
  }
  return(*this);
}

double operator^(const BlockMatrix& rhs,const BlockMatrix& lhs)
{
  double value = 0.0;
  int nirreps = rhs.get_nirreps();
  for(int h=0; h < nirreps; ++h){
    const double** rhs_matrix_block = rhs.get_block(h);
    const double** lhs_matrix_block = lhs.get_block(h);
    int block_size = rhs.get_block_size(h);
    if(block_size>0){
      for(int i = 0; i < block_size; ++i)
        for(int j = 0; j < block_size; ++j)
          value += lhs_matrix_block[i][j] * rhs_matrix_block[i][j];
    }
  }
  return(value);
}

*/
