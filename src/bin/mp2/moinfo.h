struct moinfo {
  int nmo;               /* no. of MOs */
  int nso;               /* no. of symmetry adapted MOs */
  int nirreps;           /* no. of irreducible representations */
  int nfzdocc;           /* no. of frozen occupied MOs */
  int nfzvirt;           /* no. of frozen virtual MOs */
  int nactmo;
  int nactdocc;          /* no. of active occupied MOs */
  int nactvirt;          /* no. of active virtual MOs */
  int ndocc;             /* no. of all occupied MOs */  
  int nvirt;             /* no. of all virtual MOs */   
  int *fzdoccpi;         /* frozen occupied MOs per irrep */
  int *fzvirtpi;         /* frozen virtual MOs per irrep */
  int *actdoccpi;        /* active occupied MOs per irrep */
  int *actdoccsym;
  int *actvirtpi;        /* active virtual MOs per irrep */
  int *actvirtsym;
  int *doccpi;           /* all occupied MOs per irrep */
  int *virtpi;           /* all virtual MOs per irrep */
  int *mopi;             /* all MOs per irrep */
  int *ioff;             /* ioff array */
  double Enuc;           /* Nuclear repulsion energy */
  double Escf;           /* SCF energy */
  double *scfevals;      /* SCF eigenvalues */
  char **irreplabels;
  int *docc_off;
  int *virt_off;
  int *qt_docc;
  int *qt_virt;
  double **opdm;
};
