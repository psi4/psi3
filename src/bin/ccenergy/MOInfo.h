struct MOInfo {
  int nirreps;           /* no. of irreducible representations */
  int nmo;               /* no. of molecular orbitals */
  int nso;               /* no. of symmetry orbitals */
  int iopen;             /* 0=closed shell; >0=open shell */
  int phase;             /* Boolean for consistency of orbital phases */
  int aobasis;
  int *orbspi;           /* no. of MOs per irrep */
  int *orbsym;           /* orbital symmetry (Pitzer/SO) */
  int *clsdpi;           /* no. of closed-shells per irrep excl. frdocc */
  int *openpi;           /* no. of open-shells per irrep */
  int *uoccpi;           /* no. of unoccupied orbitals per irr. ex. fruocc */
  int *frdocc;           /* no. of frozen core orbitals per irrep */
  int *fruocc;           /* no. of frozen unoccupied orbitals per irrep */
  char **labels;         /* irrep labels */
  int *occpi;            /* no. of occupied orbs. (incl. open) per irrep */
  int *aoccpi;           /* no. of alpha occupied orbs. (incl. open) per irrep */
  int *boccpi;           /* no. of beta occupied orbs. (incl. open) per irrep */
  int *virtpi;           /* no. of virtual orbs. (incl. open) per irrep */
  int *avirtpi;          /* no. of alpha virtual orbs. (incl. open) per irrep */
  int *bvirtpi;          /* no. of beta virtual orbs. (incl. open) per irrep */
  int *occ_sym;          /* relative occupied index symmetry */
  int *aocc_sym;         /* relative alpha occupied index symmetry */
  int *bocc_sym;         /* relative beta occupied index symmetry */
  int *vir_sym;          /* relative virtual index symmetry */
  int *avir_sym;         /* relative alpha virtual index symmetry */
  int *bvir_sym;         /* relative beta virtual index symmetry */
  int iter;              /* Current CCSD iteration */
  double conv;           /* Current convergence level */
  double enuc;           /* Nuclear repulsion energy */
  double emp2;           /* MP2 energy */
  double escf;           /* SCF energy (from file30) */
  double eref;           /* Reference energy (file100) */
  double ecc;            /* Current coupled cluster energy */
  double t1diag;         /* Standard open- or closed-shell T1 diagnostic */
  double d1diag;         /* Janssen and Nielsen's D1 Diagnostic */
  double new_d1diag;     /* Lee's modified D1 Diagnostic */
  double ***C;           /* Virtual orbital transformation matrix (for AO-basis B terms) */
  double ***Ca;          /* UHF alpha virtual orbital transformation matrix (for AO-basis B terms) */
  double ***Cb;          /* UHF beta virtual orbital transformation matrix (for AO-basis B terms) */
};
