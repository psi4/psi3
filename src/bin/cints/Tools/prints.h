#ifndef _psi_src_bin_cints_Tools_prints_h
#define _psi_src_bin_cints_Tools_prints_h

/*! \file prints.h
    \ingroup CINTS
*/
namespace psi { namespace CINTS {
void print_intro();
void print_scalars();
void print_basisset();
void print_basisset_aobasis();
void print_quote();
void print_opdm();
void print_atomvec(const char *quantity, double **vecs);
void print_atommat(const char *quantity, double **mat);
 void print_moinfo_corr();
}}
#endif
