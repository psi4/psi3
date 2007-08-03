/*! \file 
    \ingroup (MP2R12)
    \brief Enter brief description of file here 
*/

/*! \defgroup MP2R12 Add a description of the group MP2R12 */

/* Struct for PSIF_CHKPT molecular orbital information */
struct MOInfo {
    int nmo;
    int nirreps;
    int iopen;
    int *orbspi;
    int *clsdpi;
    int *openpi;
    int *virtpi;
    int *frdocc;
    int *fruocc;
    int *first;
    int *last;
    int noeints;
    int nteints;
    char **labels;
    double enuc;
    double escf;
    double *evals;
    double *te_ints;
};
    
