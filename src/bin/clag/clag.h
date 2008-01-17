/*! \file 
    \ingroup (CLAG)
    \brief Enter brief description of file here 
*/

#ifndef _psi3_bin_clag_clag_h_
#define _psi3_bin_clag_clag_h_

namespace psi { namespace clag {

double **rdopdm(int nbf, int print_lvl, int opdm_file) ;  
double *rdtpdm(int nbf, int print_lvl, int tpdm_file) ; 
void init_io(int argc, char **argv) ; 
void close_io(void) ;
void trace_opdm(double **opdm, int nbf);
void trace_tpdm(double *tpdm, int nbf);
double lagcalc(double **OPDM, double *TPDM, double *h, double *TwoElec, 
             double **lag, int nmo, int npop, int print_lvl, int lag_file); 
void ci_energy(double **OPDM, double *TPDM, double *h, double *TwoElec,
               int nbf, double enuc, double eci_30, double lagtr); 

}} // end namespace psi::clag

#endif // header guard
