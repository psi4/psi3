/* $Id$ */
/* $Log$
 * Revision 1.4  2000/06/22 22:14:58  evaleev
 * Modifications for KS DFT. Reading in XC Fock matrices and XC energy in formg_direct need to be uncommented (at present those are not produced by CINTS yet).
 *
/* Revision 1.3  2000/06/02 13:32:14  kenny
/*
/*
/* Added dynamic integral accuracy cutoffs for direct scf.  Added a few global
/* variables.  Added keyword 'dyn_acc'; true--use dynamic cutoffs.  Use of
/* 'dconv' and 'delta' to keep track of density convergence somewhat awkward,
/* but avoids problems when accuracy is switched and we have to wipe out density
/* matrices.  Also added error message and exit if direct rohf singlet is
/* attempted since it doesn't work.
/* --Joe Kenny
/*
/* Revision 1.2  2000/03/28 15:45:31  evaleev
/* Increased the MAX_BASIS and MAXIOFF to 4096
/*
 * Revision 1.1.1.1  2000/02/04  22:52:29  evaleev
 * Started PSI 3 repository
 *
/* Revision 1.10  1999/11/11 21:15:13  localpsi
/* Altered cscf to do some guess at the multiplicity from SOCC. -STB (11/11/99)
/*
/* OH and in case your wondering who localpsi is, it is the superuser on my pc
/* that contains my psi files.
/*
/* Revision 1.9  1999/11/04 19:24:28  localpsi
/* STB (11/4/99) - Added the orb_mix feature which is equivalent to guess = mix
/* in G94 and also fixed restarting so that if you have different wavefuntions,
/* everything works.  Also if you specify no DOCC and SOCC and restart, if the
/* wavefunctions are different, it will guess again.
/*
/* Revision 1.8  1999/11/02 23:55:55  localpsi
/* Shawn Brown - (11/2/99) Modified to the code in a few major ways.
/*
/* 1.  Added the capability to do UHF.  All of the features available with the
/* other refrences have been added for UHF.
/*
/* 2.  For UHF, I had to alter the structure of file30. (See cleanup.c for a
/* map)  This entailed adding a pointer array right after the header in the SCF
/* section of file30 that pointed to all of the data for the SCF caclulation.
/* Functions were added to libfile30 to account for this and they are
/* incorporated in this code.
/*
/* 3.  Updated and fixed all of the problems associated with my previous
/* guessing code.  The code no longer uses OPENTYPE to specify the type of
/* occupation.  The keword REFERENCE and MULTP can now be used to indicate any
/* type of calculation.  (e.g. ROHF with MULTP of 1 is an open shell singlet
/* ROHF calculation)  This code was moved to occ_fun.c.  The code can also
/* guess at any multplicity in a highspin case, provided enough electrons.
/*
/* Revision 1.7  1999/11/02 18:10:12  evaleev
/* Direct SCF improved
/*
/* Revision 1.6  1999/10/22 19:47:17  evaleev
/* A direct SCF-enabled version (set DIRECT_SCF=TRUE in input.dat).
/*
/* Revision 1.5  1999/08/17 19:04:13  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.4  1999/08/11 19:24:53  evaleev
/* Unhardwired the size of the ioff array (set it to 1024 for now) and increased MAX_BASIS to 1024.
/*
/* Revision 1.3  1999/08/11 18:39:03  evaleev
/* Added some checks on the lowest eigenvalue of the overlap matrix.
/*
/* Revision 1.2  1999/07/24 18:13:49  crawdad
/* Renamed variable "nint" to "cscf_nint" to avoid DEC compiler type conflict.
/* -Daniel
/*
 * Revision 1.1.1.1  1999/04/12  16:59:25  evaleev
 * Added a version of CSCF that can work with CINTS.
 * -Ed
 * */

#include <libciomr.h>

#define MAX_BASIS 4096
#define MAXIOFF 4096
#define SEVAL_CUTOFF 1.0E-6

#define SMAT 0
#define TMAT 1
#define VMAT 2

#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

EXTERN FILE *infile, *outfile, *JK,*gmat,*diis_out;

EXTERN double dampsv;           /* scale factor in diis */
EXTERN double repnuc;           /* nuclear repulsion */
EXTERN double etot;             /* electronic and total energies */
EXTERN double exc;              /* KS DFT exchange-correlation energy */
EXTERN double lshift;           /* levelshift */
EXTERN double diiser;           /* max off-diag. element in MO fock mat. */
EXTERN double save_ci1,save_ci2; /* ci coefficients for tcscf */
EXTERN double dampd;
EXTERN double dampo;
EXTERN double eri_cutoff;       /* accuracy of integrals to request from cints if doing direct */

EXTERN int direct_scf;          /* 1 to request direct formation of the Fock matrices */
EXTERN int diisflg;             /* 0 for diis, 1 disables diis */
EXTERN int iconv;               /* dmat convg. criterion */
EXTERN int iopen;               /* 0 for closed, 1 for open, 2 for twocon */
EXTERN int inflg;               /* 0 default, 1 use old guess, 2 use core H */
EXTERN int print;               /* print flag */
EXTERN int fock_typ;            /* 0 for default, 1 for simpler op sh fock m */
EXTERN int ndiis;               /* # of error matrices to keep in diis */
EXTERN int it_diis;             /* iteration to begin diis extrapolation */
EXTERN int itmax;               /* max iterations */
EXTERN int use_iwl;             /* use IWL format */
EXTERN int delete_ints;         /* delete ints? */
EXTERN int delete_1e;           /* delete one-electron ints? */
EXTERN int delete_2e;           /* delete two-electron ints? */

EXTERN int multp;		/* multiplicity of the molecule */
EXTERN int mflag;               /* 1 if multp specified */
EXTERN int charge;		/* charge of the molecule */
EXTERN int natom;		/* number of atoms in the molecule */
EXTERN int nelec;		/* number of electrons in the molecule */
EXTERN int nbfso;		/* total number of symmetry-adapted basis functions */
EXTERN int nmo;                 /* total number of molecular orbitals */
EXTERN char *reference;         /* RHF,UHF,ROHF,TCSCF,RKS,UKS */
EXTERN char *functional;        /* KS DFT functional name, just to print out */

typedef enum {ref_rhf = 0, ref_uhf = 1, ref_rohf = 2, ref_tcscf = 3,
      ref_rks = 4, ref_uks = 5} reftype;
EXTERN reftype refnum;

EXTERN int mo_out;              /* 1 if display orbitals in new format at end*/
EXTERN int n_so_typs;           /* number of irreps w/ non-zero num of so's */
EXTERN int nbasis;              /* # basis functions */
EXTERN int nsfmax;              /* max # of so's per irrep */
EXTERN int n_open;              /* # open shells */
EXTERN int num_ir;              /* # of symmetry types */
EXTERN int mxcoef;              /* sum of (# so's per irrep)**2 */
EXTERN int mxcoef2;             /* sum of ioff[# so's per irrep] */
EXTERN int readflg;             /* 1 if using buffered io */
EXTERN int maxbuf;              /* number of integrals per buffer */
EXTERN int num_bufs;            /* number of buffers used */
EXTERN int num_ints;            /* total integrals written to supermatrix */
EXTERN int iter;                /* iteration */
EXTERN int converged;           /* 1 if converged */
EXTERN int hsos;                /* 1 if high spin open shell */
EXTERN int singlet;             /* 1 if open shell singlet */
EXTERN int uhf;                 /* 1 if uhf 0 if RHF or ROHF */
EXTERN int special;             /* 1 if OPENTYPE=special */
EXTERN int twocon;              /* 1 if tcscf */
EXTERN int ksdft;               /* 1 if Kohn-Sham DFT */
EXTERN int mixing;              /* 1 if mixing for UHF, default is 0 */
EXTERN PSI_FPTR pos34;               /* pointer to location in file34 */
EXTERN int cscf_nint;                /* number of pki ints in present batch */
EXTERN int opshl1,opshl2;
EXTERN int opblk1,opblk2;
EXTERN int second_root;         /* get the second root of the MCSCF */
EXTERN int icheck_rot;          /* check orbital rotations? */
EXTERN int ediff;

EXTERN int itap30,itap34,itapS,itapT,itapV,itap33,itap92,itap93,itapDSCF;
EXTERN double alpha1,alpha2,alpha3;  /* two configuration things */

EXTERN int ioff[MAXIOFF];       /* matrix offsets */

EXTERN double *alpha, *beta;    /* arrays with energy coupling coeffs */
EXTERN double *zvals;		/* array for nuclear charges */
EXTERN int *symm_tot;           /* array containing the orbital symmetries in order of energy */
EXTERN double *ener_tot;        /* array containing the orbital energies in order */

EXTERN int *i10;  

EXTERN union psi_buffer {
          int *lbli;
          double *pki;
          double **pki_p;
          } oubuf;

EXTERN struct symm {
    double *smat;
    double *tmat;
    double *hmat;
    double *fock_pac;
    double *fock_open;
    double *fock_eff;
    double *fock_evals;
    double *gmat;
    double *gmato;
    double *xcmat;  /* Exchange-correlation Fock matrix for KS DFT */
    double *pmat;   /* Closed-shell density matrixin RHF, alpha or beta in UHF) */
    double *pmato;
    double *pmat2;
    double *pmato2;
    double *dpmat;
    double *dpmato;
    double **cmat;
    /* STB(4/1/98) - Added array for saving evalues of core H */
    double *hevals;
    /* TDC(6/19/96) - Added array for saving original MO vector */
    double **cmat_orig;
    double **sahalf;
    double *occ_num;
    int nclosed;
    int nopen;
    int nhalf;
    /* who in the hell needs to know the degeneracy of irreps in Abelian subgroups??? */
/*  int degeneracy;   */
    int num_so;          /* Number of SOs in this symmetry block */
    int num_mo;          /* Number of MOs in this symmetry block,
			    may be different from num_so */
    int os_num;
    int ideg;
    char *irrep_label;
    /* STB -7/2/99 I know this is a little redundant but it is for 
       UHF */
    int noccup;
} *scf_info;

/* STB - 10/11/99 - structure added to handle spin */
EXTERN struct spin {
    struct symm *scf_spin;
    char *spinlabel;
} *spin_info;

/* TDC(6/19/96) - Added flag for success or failure of phase checking
   routine */
EXTERN int phase_check;

/* EFV(10/24/98) - Added an array that maps an SO number to the symmetry block number */
EXTERN int *so2symblk;

/* JPK(6/1/00) added variables for dynamic integral accuracy in direct scf*/
EXTERN int tight_ints, ok_ints,    /*keeps track of acccuracy being used*/
       dyn_acc,                    /*1 for dynamic integral accuracy, else 0*/    
       acc_switch;                 /*accuracy switch:  1 -> accuracy has been 
                                     switched*/
EXTERN double delta;               /*just another density convergence
                                     variable*/
