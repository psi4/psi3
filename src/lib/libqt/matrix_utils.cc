/*! \file
    \ingroup QT
    \brief misc matrix operations
*/
#include <stdio.h>
#include <stdlib.h>
#include <libqt/qt.h>
#include <math.h>

extern "C" {

  /// Canonicalize phases of matrix A
  /// phases are canonical when the largest-magnitude coefficient in each column is positive
  void canonicalize_column_phases(double** A, int nrow, int ncol) {

    const double soft_zero = 1e-12;

    // in each column
    for(int c=0; c<ncol; ++c) {
      // find the largest-magnitude element
      int rmax = -1;
      int valmax = 0.0;
      for(int r=0; r<nrow; ++r) {
        const double value = fabs(A[r][c]);
        if (value - valmax > soft_zero) {
          valmax = value;
          rmax = r;
        }
      }
      // if the largest-magnitude element is negative, scale the column by -1
      if (A[rmax][c] < -soft_zero) {
        for(int r=0; r<nrow; ++r) {
          A[r][c] *= -1.0;
        }
      }
    }

  }

  /// Canonicalize phases of matrix A
  /// phases are canonical when the largest-magnitude coefficient in each row is positive
  void canonicalize_row_phases(double** A, int nrow, int ncol) {

    const double soft_zero = 1e-12;

    // in each row
    for(int r=0; r<nrow; ++r) {
      // find the largest-magnitude element
      int cmax = -1;
      int valmax = 0.0;
      for(int c=0; c<ncol; ++c) {
        const double value = fabs(A[r][c]);
        if (value - valmax > soft_zero) {
          valmax = value;
          cmax = c;
        }
      }
      // if the largest-magnitude element is negative, scale the row by -1
      if (A[r][cmax] < -soft_zero) {
        for(int c=0; c<ncol; ++c) {
          A[r][c] *= -1.0;
        }
      }
    }

  }

} /* extern "C" */

