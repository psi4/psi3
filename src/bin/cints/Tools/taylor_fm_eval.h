/*! \file taylor_fm_eval.h
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/

void init_Taylor_Fm_Eval(unsigned int mmax, double epsilon);
void taylor_compute_fm(double *F, double T, unsigned int l);
void free_Taylor_Fm_Eval();
