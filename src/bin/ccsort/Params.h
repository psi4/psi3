struct Params {
    int print_lvl;         /* Output level control */
    int TEIFile;           /* Unit no. for input two-electron integrals */
    int keep_TEIFile;      /* Should we keep the input two-elec. integrals? */
    int OEIFile;           /* Unit no. for input one-electron integrals */
    int keep_OEIFile;      /* Should we keep the input one-elec. integrals? */
    int FZCFile;           /* Unit no. for input frozen core operator */
    double tolerance;      /* Cutoff value for integrals in IWL Buffers */
    int memory;            /* Memory available (in bytes) */
};
