
#ifndef _psi3_dboc_linalg_h_
#define _psi3_dboc_linalg_h_

#include "float.h"
#include <stdio.h>

FLOAT** create_matrix(int a, int b);
void delete_matrix(FLOAT** M);
FLOAT** convert_matrix(double **M, int a, int b, int transpose);
void print_mat(FLOAT** a, int m, int n, FILE* out);
int matrix_mult(FLOAT** A, int arow, int acol, FLOAT** B, int brow, int bcol, FLOAT** C);
void lu_decom(FLOAT** a, int n, int* indx, FLOAT* d);

#endif
