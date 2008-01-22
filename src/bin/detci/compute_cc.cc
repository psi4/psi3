/*
** COMPUTE_CC.CC
**
** This will be a prototype implementation of arbitrary-order 
** coupled-cluster code.
** 
** C. David Sherrill
** Center for Computational Molecular Science and Technology
** Georgia Institute of Technology
** March 2005
**
** Note: I think I need onel ints as g for formation of sigma
** in non-FCI cases, but make sure any CC parts don't try to get
** h and actually get g instead...
*/

#define EXTERN

#include <math.h>

extern "C" {
   #include <stdio.h>
   #include <stdlib.h>
   #include <psifiles.h>
   #include <libipv1/ip_lib.h>
   #include <libciomr/libciomr.h>
   #include <libqt/qt.h>
   #include <libqt/slaterdset.h>
   #include <libpsio/psio.h>
   #include "structs.h"
   #include "globals.h"
   #include "cc.h"

   extern unsigned char ***Occs;
   extern struct olsen_graph *AlphaG;
   extern struct olsen_graph *BetaG;
   extern struct stringwr **alplist;
   extern struct stringwr **betlist;

   extern BIGINT strings2det(int alp_code, int alp_idx,
     int bet_code, int bet_idx);
                                                         
   extern int og_lex_addr(struct olsen_graph *Graph, int *occs, int nel,
     int *listnum);

   extern void print_vec(unsigned int nprint, int *Iacode, int *Ibcode,
     int *Iaidx, int *Ibidx, double *coeff,
     struct olsen_graph *AlphaG, struct olsen_graph *BetaG,
     struct stringwr **alplist, struct stringwr **betlist,
     FILE *outfile);

   extern void xpeay(double *x, double a, double *y, int size);

   extern void orb2lbl(int orbnum, char *label);
   extern int  lbl2orb(char *label);
   
   struct T_blks Tblks; 
   struct diis_obj diis;

   extern void parse_import_vector(SlaterDetSet *sdset, int *i_alplist,
      int *i_alpidx, int *i_betlist, int *i_betidx, int *i_blknums);

}

#include <iostream>
#define ODOMETER_DECL_ONLY
#include "odometer.h"
                                                                                
#include "civect.h"

int cc_reqd_sblocks[CI_BLK_MAX];

/*
** compute_cc()
**
** This is the top-level function that controls the coupled-cluster
** computation
**
*/
void compute_cc(void)
{
  printf("compute_cc: Not yet available\n");
}
