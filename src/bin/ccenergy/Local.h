struct Local {
  int natom;
  int nao;
  int nocc;
  int nvir;
  int *aostart;
  int *aostop;
  int **pairdomain;
  int *pairdom_len;
  int *pairdom_nrlen;
  double ***V;
  double ***W;
  double *eps_occ;
  double **eps_vir;
  double cutoff;
};
