/*! \file shell_block_matrix.h
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/

double ****init_shell_block_matrix(void);
void free_shell_block_matrix(double****);
void shell_block_to_block(double****, double**);
void GplusGt(double****, double****);
