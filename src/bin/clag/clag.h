/****************************************************************************/
/* lag.c header file including standard definitions and prototypes          */
/****************************************************************************/
 
double **rdopdm(int nbf, int print_lvl, int opdm_file) ;  
double *rdtpdm(int nbf, int print_lvl, int tpdm_file) ; 
void init_io(void) ; 
void close_io(void) ;
void trace_opdm(double **opdm, int nbf);
void trace_tpdm(double *tpdm, int nbf);
double lagcalc(double **OPDM, double *TPDM, double *h, double *TwoElec, 
             double **lag, int nmo, int npop, int print_lvl, int lag_file); 
void ci_energy(double **OPDM, double *TPDM, double *h, double *TwoElec,
               int nbf, double enuc, double eci_30, double lagtr); 
void onel_to_cas(double *onel_ints, int *corr_to_pitz, int nso,
                  int print_lvl, int itap);
void twoel_to_cas(double *twoel_ints, int *corr_to_pitz, int nso,
                  int print_lvl, int itap);
void onepdm_to_cas(double **onepdm, int *corr_to_pitz, int nso, int npop,
                   int print_lvl, int itap);
void twopdm_to_cas(double *tpdm, int *corr_to_pitz, int nso,
                   int npop, int print_lvl, int itap);
void lag_to_cas(double **lag, int *corr_to_pitz, int nso,
                int print_lvl, int itap);

