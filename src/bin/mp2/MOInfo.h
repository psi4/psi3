/* Struct for file30 molecular orbital information */
struct MOInfo {
    int nmo;
    int nirreps;
    int iopen;
    int *orbspi;
    int *clsdpi;
    int *openpi;
    int *virtpi;
    int *fruocc;
    int *frdocc;
    int noeints;
    int nteints;
    char **labels;
    double enuc;
    double escf;
    double *evals;
    double *te_ints;
};
    
