/* 
** OLSENGRAPH.CC
**
** Routines needed to maintain the Olsen Graph Object.  The Olsen graph
** for alpha/beta strings has a different subgraph for each value of
** irrep and RAS I hole/RAS III particle combination.  
**
** C. David Sherrill, May 1995
** Based on previous code by David Sherrill, 1994
**
*/

#define EXTERN

extern "C" {
   #include <stdlib.h>
   #include <stdio.h>
   #include <libciomr/libciomr.h>
   #include <libqt/qt.h>
   #include "structs.h"
   #include "globals.h"
   extern struct stringwr **alplist;  
   extern struct stringwr **betlist; 

   extern void stringlist(struct olsen_graph *Graph, struct stringwr **slist);
   extern void set_ciblks(struct olsen_graph *AG, struct olsen_graph *BG);
   extern void print_ci_space(struct stringwr *strlist, int num_strings,
      int nirreps, int strtypes, int nel, FILE *outfile);
   extern void str_abs2rel(int absidx, int *relidx, int *listnum,
      struct olsen_graph *Graph);
} ;

#include <iostream>
#define ODOMETER_DECL_ONLY
#include "odometer.h"


//#define DEBUG  

/* FUNCTION PROTOTYPES for this module */
void olsengraph(struct olsen_graph *Graph, int ci_orbs, int num_el, 
      int nirreps, int *orbsym, int ras1_lvl, int ras1_min, int ras1_max, 
      int ras3_lvl, int ras3_max, int num_fzc_orbs, int num_cor_orbs,
      int ras4_lvl, int ras4_max, int ras34_max);
void og_add_walk(int ras1_idx, int ras3_num, int ras4_num, int *occs, 
      int nel_expl, int norb, int nirreps, int num_fzc_orbs,
      struct olsen_graph *Graph);
int og_calc_y(struct level *lvl, int ci_orbs);
void og_fill(int num_el, int norb, int nirreps, int num_fzc_orbs,
      struct olsen_graph *Graph);
int subgr_lex_addr(struct level *head, int *occs, int nel, int norb);
int og_lex_addr(struct olsen_graph *Graph, int *occs, int nel,
      int *listnum);
void og_print(struct olsen_graph *Graph, FILE *outfile);



/*
** form_strings(): Form the alpha and beta string graphs.  Use a different
**    graph for every irrep * (num el in RAS I) * (num el in RAS III)
**
*/
void form_strings(void)
{
   int i, nlists, nirreps, ncodes;
   int irrep, code, listnum;
   int *occs;

   AlphaG = new olsen_graph[1];

   // Make the graph
   olsengraph(AlphaG, CalcInfo.num_ci_orbs, CalcInfo.num_alp,
      CalcInfo.nirreps, CalcInfo.orbsym,
      Parameters.a_ras1_lvl, Parameters.a_ras1_min, Parameters.a_ras1_max,
      Parameters.ras3_lvl, Parameters.a_ras3_max,
      CalcInfo.num_fzc_orbs, CalcInfo.num_cor_orbs,
      Parameters.ras4_lvl, Parameters.a_ras4_max, Parameters.a_ras34_max);

   if (Parameters.print_lvl > 3)
      og_print(AlphaG, outfile) ;

   ncodes = AlphaG->subgr_per_irrep;
   nirreps = AlphaG->nirreps;
   nlists = nirreps * ncodes;
   /* alplist = new stringwr*[nlists]; */
   alplist = (struct stringwr **) malloc(nlists * sizeof(struct stringwr *));
   for (i=0; i<nlists; i++) alplist[i] = NULL;

   stringlist(AlphaG, alplist);

   if (Parameters.print_lvl>=4) {
      for (irrep=0,listnum=0; irrep < nirreps; irrep++) {
         for (code=0; code < ncodes; code++, listnum++) {
            fprintf(outfile, "Strings for irrep %d code %2d (list %2d)\n", 
               irrep, code, listnum);
            print_ci_space(alplist[irrep * ncodes + code],
               AlphaG->sg[irrep][code].num_strings,
               nirreps, nlists, AlphaG->num_el_expl, outfile) ;
            }
         }
      }

   /* for beta string graph if necessary */
   if (CalcInfo.iopen && !(Parameters.Ms0)) {

      BetaG = new olsen_graph[1];

      olsengraph(BetaG, CalcInfo.num_ci_orbs, CalcInfo.num_bet,
          CalcInfo.nirreps, CalcInfo.orbsym,
          Parameters.b_ras1_lvl,  Parameters.b_ras1_min, Parameters.b_ras1_max,
          Parameters.ras3_lvl, Parameters.b_ras3_max,
          CalcInfo.num_fzc_orbs, CalcInfo.num_cor_orbs,
          Parameters.ras4_lvl, Parameters.b_ras4_max, Parameters.b_ras3_max);

      if (Parameters.print_lvl > 3)
         og_print(BetaG, outfile) ;

      ncodes = BetaG->subgr_per_irrep;
      nirreps = BetaG->nirreps;
      nlists = nirreps * ncodes;
      betlist = (struct stringwr **) malloc(nlists * sizeof(struct stringwr *));
      for (i=0; i<nlists; i++) betlist[i] = NULL;


      stringlist(BetaG, betlist);

      if (Parameters.print_lvl>=4) {
         for (irrep=0; irrep < nirreps; irrep++) {
            for (code=0; code < ncodes; code++) {
               fprintf(outfile, "Strings for irrep %d code %2d\n", irrep,
                  code);
               print_ci_space(betlist[irrep * ncodes + code],
                  BetaG->sg[irrep][code].num_strings,
                  nirreps, nlists, BetaG->num_el_expl, outfile) ;
               }
            }
         }
      }

   else {
      betlist = alplist;
      BetaG = AlphaG;
      }

   /* get number of alpha/beta strings, ref symmetry, etc */
   set_ciblks(AlphaG, BetaG) ;

   /* if the user wants to filter out some initial guesses based on
      phases of two determinants, we need to convert their absolute
      string numbers into string lists (subgraphs) and relative 
      addresses 
    */
   if (Parameters.filter_guess) {
     str_abs2rel(Parameters.filter_guess_Ia, &Parameters.filter_guess_Iaridx, 
                 &Parameters.filter_guess_Iac, AlphaG);
     str_abs2rel(Parameters.filter_guess_Ib, &Parameters.filter_guess_Ibridx, 
                 &Parameters.filter_guess_Ibc, BetaG);
     str_abs2rel(Parameters.filter_guess_Ja, &Parameters.filter_guess_Jaridx, 
                 &Parameters.filter_guess_Jac, AlphaG);
     str_abs2rel(Parameters.filter_guess_Jb, &Parameters.filter_guess_Jbridx, 
                 &Parameters.filter_guess_Jbc, BetaG);
   }
}



/*
** olsengraph(): Form an Olsen/Roos graph for the alpha/beta strings.
** 
** Unfortunately, I found it necessary to split the DRT graph into many
**    subgraphs, as described in Olsen & Roos 1988.  Each graph must 
**    contain walks of only one irrep, and with only one value for the 
**    number of electrons at the RAS I level.  All these graphs are 
**    stored in a collective container, struct olsen_graph.  For open
**    shell systems, there will be two of these, one for alpha and one
**    for beta.  For closed shell systems only the alpha graph is stored.
**
** Even more unfortunately, each subgraph must also be labeled according
**    to the number of electrons in RAS III, or so I believe I figured
**    out at one time.  That's how I'll do it now.  [CDS 5/95]
**
** Even more more unfortunately, I'm going to split it according to 
**    number of electrons in RAS IV!!!  That's just crazy.  [CDS 8/95]
**
** For a full CI, make all strings with a given irrep belong to the same
**    graph.  (6/19/95)
**
** Parameters:
**    Graph        =  struct olsen_graph to hold all subgraphs
**    ci_orbs      =  number of orbitals explicitly in CI 
**    num_el       =  number of electrons for the relevant string
**    nirreps      =  number of irreducible representations
**    orbsym       =  orbital symmetry array
**    ras1_lvl     =  last level in RAS I
**    ras1_min     =  min number of electrons at RAS I level for the string 
**                    (nb this definition includes core electrons!)
**    ras1_max     =  maximum number of electrons at RAS I level 
**    ras3_lvl     =  first level in RAS III
**    ras3_max     =  max number of electrons in RAS III _for the string_
**    num_fzc_orbs = number of frozen core orbitals
**    num_cor_orbs = number of restricted core orbitals
**    ras4_lvl     =  first level of the new RAS IV 
**    ras4_max     =  max number of electrons in RAS IV for the string
**    ras34_max    =  max number of electrons in RAS III and IV
**
** Returns: none
*/
void olsengraph(struct olsen_graph *Graph, int ci_orbs, int num_el, 
      int nirreps, int *orbsym, int ras1_lvl, int ras1_min, int ras1_max, 
      int ras3_lvl, int ras3_max, int num_fzc_orbs, int num_cor_orbs,
      int ras4_lvl, int ras4_max, int ras34_max)
{
   Odometer Ras1, Ras2, Ras3, Ras4 ;
   int n1, n2, n3, n4 ;
   int n1max, n1min ;
   int max_el_ras1 ;
   int *occs, *array1, *array2, *array3, *array4, **encode_tmp;
   int i, j, k ;
   int orbs_frozen, maxj, fzc_sym=0, code=0;
   struct stringgraph *sgptr;


   #ifdef DEBUG
   fprintf(outfile, "ras1_lvl = %d   ras1_min = %d  ras1_max = %d\n",
      ras1_lvl, ras1_min, ras1_max) ;
   fprintf(outfile, "ras3_lvl = %d   ras3_max = %d\n", ras3_lvl, ras3_max) ;
   #endif

   // Go ahead and set the occupations of the frozen orbs 
   occs = init_int_array(num_el) ;
   for (i=0; i<num_cor_orbs; i++) occs[i] = i ;
 
   orbs_frozen = num_fzc_orbs + num_cor_orbs;

   for (i=0; i<num_fzc_orbs; i++) fzc_sym ^= orbsym[i];

   maxj = nirreps * (num_el - num_fzc_orbs + 1);

   // go ahead and make room for the occupations of RAS I, II, and III
   array1 = init_int_array(num_el) ;
   array2 = init_int_array(num_el) ;
   array3 = init_int_array(num_el) ;
   array4 = init_int_array(num_el) ;

   // Initialize the Graph data structure
   Graph->num_el = num_el;
   Graph->num_el_expl = num_el - num_fzc_orbs;
   Graph->num_orb = ci_orbs ;
   Graph->num_fzc_orbs = num_fzc_orbs;
   Graph->num_cor_orbs = num_cor_orbs;
   Graph->fzc_sym = fzc_sym;
   Graph->orbsym = orbsym;
   Graph->ras1_lvl = ras1_lvl ;
   Graph->ras1_min = ras1_min ;
   Graph->ras1_max = ras1_max ;
   Graph->ras3_lvl = ras3_lvl ; 
   Graph->ras3_max = ras3_max ;
   Graph->ras4_lvl = ras4_lvl ;
   Graph->ras4_max = ras4_max ;
   Graph->nirreps = nirreps;
   Graph->str_per_irrep = init_int_array(nirreps);
   n1max = ras1_max - orbs_frozen; 
   n1min = ras1_min - orbs_frozen;

   Graph->decode = (int ***) malloc ((ras1_max - ras1_min + 1) *
      sizeof(int **));
   for (i=0; i<(ras1_max - ras1_min + 1); i++) {
      Graph->decode[i] = init_int_matrix(ras3_max + 1, ras4_max + 1);
      }

   encode_tmp = init_int_matrix(3, (ras1_max - ras1_min + 1) * (ras3_max + 1)
                   * (ras4_max + 1));

   /* need to know how many possible RAS I/RAS III/RAS IV combinations */
   if (!Parameters.fci_strings) {
      for (i=ras1_max; i>=ras1_min; i--) {
         for (j=0; j<=ras3_max; j++) {
            for (k=0; k<=ras4_max; k++) {
               if ((i+j+k<=num_el) && (num_el-i-j-k<=ras3_lvl-ras1_lvl-1)
                   && (j+k <= ras34_max) && (!(Parameters.r4s && k>=2 &&
                   ras1_max - i > Parameters.ex_lvl))) {
                  Graph->decode[i-ras1_min][j][k] = code;
                  encode_tmp[0][code] = i - num_fzc_orbs; 
                  encode_tmp[1][code] = j;
                  encode_tmp[2][code] = k;
                  code++;
                  }
               else Graph->decode[i-ras1_min][j][k] = -1;
               }
            }
         }
      }
   else { /* all strings w/ given irrep belong to same graph for FCI */
      for (i=ras1_max; i>=ras1_min; i--) {
         for (j=0; j<=ras3_max; j++) {
            for (k=0; k<=ras4_max; k++) {
               if ((i+j+k<=num_el) && (num_el-i-j-k<=ras3_lvl-ras1_lvl-1)
                   && (j+k <= ras34_max)) {
                  Graph->decode[i-ras1_min][j][k] = 0;
                  }
               else Graph->decode[i-ras1_min][j][k] = -1;
               }
            }
         }
      code = 1;
      }

   Graph->encode = init_int_matrix(3,code);
   for (i=0; i<code; i++) {
      Graph->encode[0][i] = encode_tmp[0][i];
      Graph->encode[1][i] = encode_tmp[1][i];
      Graph->encode[2][i] = encode_tmp[2][i];
      } 
   free_int_matrix(encode_tmp, 3);
 
 
   Graph->subgr_per_irrep = code;
   Graph->sg = (struct stringgraph **) malloc (nirreps * 
      sizeof(struct stringgraph *));

   Graph->list_offset = init_int_array(nirreps * Graph->subgr_per_irrep);

   for (i=0; i<nirreps; i++) {

      (Graph->sg)[i] = (struct stringgraph *) malloc (code *
         sizeof(struct stringgraph));

      for (j=0; j<code; j++) {

         sgptr = (Graph->sg)[i] + j;
         sgptr->lvl = (struct level *) malloc ((ci_orbs+1) *
            sizeof(struct level));
         sgptr->ktmp = (int ***) malloc (2 * sizeof (int **));
         sgptr->ktmp[0] = init_int_matrix(maxj, ci_orbs);
         sgptr->ktmp[1] = init_int_matrix(maxj, ci_orbs);

         for (k=0; k < ci_orbs+1; k++) {
            (sgptr->lvl)[k].num_j = 0;
            }
         }
      }

   // loop over the possible number of e- in RAS I (n1) and III (n3).
   // and now IV (n4)
   // the number of electrons in RAS II is defined via 
   //    n2 = num_el_this_spin - n1 - n3 - n4
   //
   // Employ the very useful Generalized Odometer
   //

   for (n1 = n1max; n1 >= n1min; n1--) {
      Ras1.resize(n1) ;
      Ras1.set_min_lex(num_cor_orbs) ;
      Ras1.set_max_lex(ras1_lvl) ;

      for (n3 = 0; n3 <= ras3_max; n3++) {

         Ras3.resize(n3) ;
         Ras3.set_min_lex(ras3_lvl) ;
         /* Ras3.set_max_lex(ci_orbs-1) ; */
         Ras3.set_max_lex(ras4_lvl-1);

         for (n4 = 0; n4 <= ras4_max && n4 <= ras34_max - n3; n4++) {

            n2 = num_el - orbs_frozen - n1 - n3 - n4;
            if (n2 < 0 || n2 > ras3_lvl - ras1_lvl - 1) continue ; 

            /* CDS 8/24/95 */
            if (Parameters.r4s && n4 >= 2 && n1max - n1 > Parameters.ex_lvl) 
               continue;

            #ifdef DEBUG
            fprintf(outfile, "n1 = %d, n2 = %d, n3 = %d, n4 = %d\n", 
               n1, n2, n3, n4) ;
            if (n2 < 0) printf("Error: n2 < 0 in form_strings()\n") ;
            #endif

            Ras2.resize(n2) ;
            Ras4.resize(n4) ; 
            Ras2.set_min_lex(ras1_lvl+1) ;
            Ras2.set_max_lex(ras3_lvl-1) ;
            Ras4.set_min_lex(ras4_lvl);
            Ras4.set_max_lex(ci_orbs-1);

            Ras1.reset() ; Ras2.reset() ; Ras3.reset() ;  Ras4.reset() ;

            do {
               Ras1.get_value(array1) ;
               do {
                  Ras2.get_value(array2) ;
                  do {
                     Ras3.get_value(array3) ;
                     do {
                        Ras4.get_value(array4) ;
                        for (i=n1-1, j = num_cor_orbs; i>=0; i--)
                           occs[j++] = array1[i] ;
                        for (i=n2-1; i>=0; i--)
                           occs[j++] = array2[i] ;
                        for (i=n3-1; i>=0; i--)
                           occs[j++] = array3[i] ;
                        for (i=n4-1; i>=0; i--)
                           occs[j++] = array4[i] ;
                        
                        // print out occupations for debugging
                        #ifdef DEBUG
                        for (i=0; i<num_el - num_fzc_orbs; i++) 
                           fprintf(outfile, "%2d ", occs[i]) ;
                        fprintf(outfile, "\n") ;
                        #endif
                  
                        // add this walk to the graph
                        og_add_walk(n1-n1min, n3, n4, occs, 
                           num_el-num_fzc_orbs, 
                           ci_orbs, nirreps, num_fzc_orbs, Graph) ;

                        Ras4.increment_lex() ;
                     } while (!Ras4.at_min()) ;
                     Ras3.increment_lex() ;
                     } while (!Ras3.at_min()) ;
                  Ras2.increment_lex() ;
                  } while (!Ras2.at_min()) ;
               Ras1.increment_lex() ;
               } while (!Ras1.at_min()) ;

            } /* end loop over n4 */

         }  /* end loop over n3 */

      } /* end loop over n1 */

   free(array1) ;
   free(array2) ;
   free(array3) ;
   free(array4) ;
   free(occs) ;

   /* fill up the olsen graph from the ki's */
   og_fill(num_el - num_fzc_orbs, ci_orbs, nirreps, num_fzc_orbs, Graph);

   /* free the temporary storage for the links */
   for (i=0; i<nirreps; i++) {
      for (j=0; j<code; j++) {
         sgptr = Graph->sg[i] + j;
         free_int_matrix(sgptr->ktmp[0], maxj);
         free_int_matrix(sgptr->ktmp[1], maxj);
         }
      }

   return; 
}



/*
** og_add_walk(): Add a walk to a subgraph within the Olsen/Roos scheme
**    of subgraphs.  Uses struct olsen_graph.
**
** Parameters:
**    ras1_idx     = Number of e- in RAS I - minimum # of e- in RAS I
**                   for the given _string_
**    ras3_num     = Number of electrons in RAS III
**    ras4_num     = Number of electrons in RAS IV
**    occs         = array listing orbital each electron occupies 
**    nel_expl     = number of electrons explicitly treated (i.e. minus
**                   all implicitly treated frozen core electrons)
**    norb         = number of orbitals _explicitly_ included 
**    nirreps      = number of irreps
**    num_fzc_orbs = number of frozen core orbitals
**    Graph        = Olsen Graph structure containing all subgraphs for a 
**                   given electron spin (alpha or beta)
**
** Returns: none
**
** Note: The newidx code excludes frozen core electrons from the num_el
**       factor
*/
void og_add_walk(int ras1_idx, int ras3_num, int ras4_num, int *occs, 
      int nel_expl, int norb, int nirreps, int num_fzc_orbs,
      struct olsen_graph *Graph)
{
   int i;
   int irrep;
   int orb = 0;
   int cur_el = 0, idx, newidx, cur_b = 0;
   int ***ktp;
   struct stringgraph *subgraph;
   int code;
   int *orbsym;

   orbsym = Graph->orbsym;
   orbsym += num_fzc_orbs;
   irrep = Graph->fzc_sym;
   idx = irrep + 1;

   /* figure out the irrep for this walk */
   for (i=0; i<nel_expl; i++) {
      orb = occs[i] ;
      irrep ^= orbsym[orb] ;
      }

   /* get a pointer to the right subgraph */
   code = Graph->decode[ras1_idx][ras3_num][ras4_num];
   subgraph = Graph->sg[irrep] + code;

   if (subgraph == NULL) {
      printf("Error (og_add_walk): NULL subgraph pointer\n");
      return;
      }
   if (code < 0) {
      printf("Error (og_add_walk): negative RAS code\n");
      return;
      }

   /* loop over all (explicitly included) orbitals */
   cur_b = Graph->fzc_sym;
   ktp = subgraph->ktmp;

   for (i=0; i<norb; i++) {

      /* if the current orbital is occupied */
      if (cur_el < nel_expl && occs[cur_el] == i) { 
         cur_el++ ;
         cur_b ^= orbsym[i];
         newidx = cur_el * nirreps + cur_b + 1;
         ktp[1][idx-1][i] = newidx;
         idx = newidx;
         }
      
      else {
         ktp[0][idx-1][i] = idx;
         }

      } /* end loop over i */

}



/*
** og_fill(): Fill out the Olsen-Roos subgraph DRT's.  So far, all we have
**    is a list of preliminary links (level.ki's).  The values of the
**    preliminary links are given according to the formula
**
**          node_index = a * irreps + b + 1 
**
**   (recall that 1 is added so that a 0 can mean no link).  Once again,
**   my 'b' is the symmetry species (0 to nirrep-1) and 'a' is the number
**   of electrons accounted for so far in the partial walk.
**
**   Uses struct olsen_graph.
**
** Parameters:
**    num_el   = total of explicit electrons for string
**    norb     = number of CI orbitals
**    nirreps  = number of irreps
**    Graph    = Olsen Graph structure containing all subgraphs for a given
**               electron spin (alpha or beta)
**
** Returns: none
**
** David Sherrill, July 1994
*/
void og_fill(int num_el, int norb, int nirreps, int num_fzc_orbs,
      struct olsen_graph *Graph) 
{
   int maxj, i, j, s, m, a, b, code ;
   int newa, newb, newj ;
   int ras1_idx, irrep ;
   int idx, newidx ;
   int *temp_a, *temp_b, **temp_kbar;
   struct level *curr, *next ;
   struct stringgraph *subgraph;
   int offset, absoffset, tot_num_str = 0;
   int max_str=0;
   int ncodes;

   // maxj is the max number of distinct rows per level
   maxj = nirreps * (num_el + 1);

   temp_a = init_int_array(maxj) ;
   temp_b = init_int_array(maxj) ;
   temp_kbar = init_int_matrix(2,maxj);

   ncodes = Graph->subgr_per_irrep;

   /* loop over all the subgraphs */
   for (irrep=0,absoffset=0; irrep<nirreps; irrep++) {
      for (code=0,offset=0; code < ncodes; code++) {
         Graph->list_offset[irrep*ncodes+code] = absoffset;

         subgraph = (Graph->sg)[irrep] + code;
         if (subgraph == NULL) {
            printf("Error (fill_og): NULL subgraph pointer\n") ;
            return;
            }

         /* set up the first row, make sure this subgraph is used */
         a = 0;
         b = Graph->fzc_sym;
         if (subgraph->ktmp[0][b][0] == 0 && subgraph->ktmp[1][b][0] == 0) {
            subgraph->num_strings = 0;
            free(subgraph->lvl);
            continue;
            }
         else {
           curr = subgraph->lvl;
           curr->num_j = 1;
           curr->a = init_int_array(1);
           curr->b = init_int_array(1);
           curr->k = init_int_matrix(2,1);
           curr->kbar = init_int_matrix(2,1);
           curr->x = init_int_array(1);
           curr->y = init_int_array(1);
           (curr->a)[0] = num_fzc_orbs;
           (curr->b)[0] = Graph->fzc_sym;
           (curr->x)[0] = 1;
           }  

         for (i=0; i<norb; i++) {
            curr = subgraph->lvl + i ; /* pointer to current row */
            next = subgraph->lvl + i + 1 ;
            zero_int_matrix(temp_kbar, 2, maxj);
 
            for (j=0; j < curr->num_j; j++) {
               a = curr->a[j];
               b = curr->b[j];
               idx = (a - num_fzc_orbs) * nirreps + b + 1;

               for (s=0; s<2; s++) {
                  newidx = subgraph->ktmp[s][idx-1][i];
                  if (newidx) {
                     newa = (newidx - 1) / nirreps + num_fzc_orbs ; 
                     newb = (newidx - 1) % nirreps ;

                     // now get a j for this new row if it's already
                     // been seen, else assign it the next ava j value
                     for (m=0; m<(next->num_j); m++) {
                        if (temp_a[m] == newa && temp_b[m] == newb) {
                           newj = m+1;
                           break;
                           }
                        }

                     if (m == next->num_j) {
                        temp_a[m] = newa ;
                        temp_b[m] = newb ;
                        next->num_j++ ;
                        newj = m+1 ;
                        }

                     (curr->k)[s][j] = newj ;   /* j may be wrong var here */
                     temp_kbar[s][newj-1] = j+1 ; /* "" */
                     }                     
                  }
               }

            next->a = init_int_array(next->num_j);
            next->b = init_int_array(next->num_j);
            next->k = init_int_matrix(2, next->num_j);
            next->y = init_int_array(next->num_j);
            next->x = init_int_array(next->num_j);
            next->kbar = init_int_matrix(2, next->num_j);

            for (j=0; j<(next->num_j); j++) {
               (next->a)[j] = temp_a[j];
               (next->b)[j] = temp_b[j];
               for (s=0; s<2; s++) 
                  (next->kbar)[s][j] = temp_kbar[s][j];
               }

            } /* end loop over i (orbitals) */

         /* probably best to do x,y's of subgraph here because we know
            the current one has >= 1 string */
         subgraph->num_strings = og_calc_y(subgraph->lvl, norb);
         subgraph->offset = offset;
         offset += subgraph->num_strings;
         absoffset += subgraph->num_strings;
         } /* end loop over code */

      Graph->str_per_irrep[irrep] = offset;
      if (offset > max_str) max_str = offset;
      tot_num_str += offset;
      } /* end loop over irreps */ 

   Graph->num_str = tot_num_str;
   Graph->max_str_per_irrep = max_str;

   free(temp_a);
   free(temp_b);
   free(temp_kbar);
}


/*
** og_calc_y() : Function calculates the lexical indices y
**
**
** Returns: The number of strings in the subgraph
*/
int og_calc_y(struct level *lvl, int ci_orbs)
{
   struct level *curr, *next ;
   int i, j, s ;
   int xcur ;
   int ksij, ksijbar ;
   int num_strings = 0;


   /* first get the x's (vertex weights) */
   for (i=0; i<ci_orbs; i++) {

      curr = lvl + i;
      next = lvl + i + 1;

      for (j=0; j<(curr->num_j); j++) {

         xcur = (curr->x)[j] ;

         for (s=0; s<2; s++) {
            ksij = (curr->k)[s][j] - 1 ;
            if (ksij >= 0) (next->x)[ksij] += xcur ;
            }
         }

      }

   /* count up how many strings */
   for (i=0; i<next->num_j; i++) 
      num_strings += next->x[i];


   /* now get the y's (arc weights) */
   for (i=0; i<ci_orbs; i++) {

      curr = lvl + i;
      next = lvl + i + 1;

      for (j=0; j<(curr->num_j); j++) {
         ksij = (curr->k)[1][j] - 1 ;
         if (ksij < 0) (curr->y)[j] = 0 ;
         else {
            ksijbar = (next->kbar)[0][ksij] - 1  ;
            if (ksijbar >= 0)
               (curr->y)[j] = (curr->x)[ksijbar] ;
            }
         }
      }

   return(num_strings);
}



void og_print(struct olsen_graph *Graph, FILE *outfile)
{

   struct level *curr;
   int ras1_min, ras1_max, ras3_max, ras4_max, code;
   int i,j,k,l,a,b;
   struct stringgraph *subgraph;

   ras1_min = Graph->ras1_min;
   ras1_max = Graph->ras1_max;
   ras3_max = Graph->ras3_max;
   ras4_max = Graph->ras4_max;

   fprintf(outfile,"\nOlsen Graph:\n");
   fprintf(outfile,"%3c%2d Electrons\n",' ',Graph->num_el);
   fprintf(outfile,"%3c%2d Frozen core orbitals\n",' ',Graph->num_fzc_orbs);
   fprintf(outfile,"%3c%2d Restricted core orbs\n",' ',Graph->num_cor_orbs);
   fprintf(outfile,"%3c%2d Explicit electrons\n",' ',Graph->num_el_expl);
   fprintf(outfile,"%3c%2d Explicit Orbitals\n",' ',Graph->num_orb);
   fprintf(outfile,"%3c%2d RAS I level\n",' ',Graph->ras1_lvl);
   fprintf(outfile,"%3c%2d RAS I minimum\n",' ',ras1_min);
   fprintf(outfile,"%3c%2d RAS I maximum\n",' ',ras1_max);
   fprintf(outfile,"%3c%2d RAS III level\n",' ',Graph->ras3_lvl);
   fprintf(outfile,"%3c%2d RAS III maximum\n",' ',ras3_max);
   fprintf(outfile,"%3c%2d RAS IV maximum\n",' ',ras4_max);
   fprintf(outfile,"%3c%2d Number of irreps\n",' ',Graph->nirreps);
   fprintf(outfile,"%3c%2d Subgraphs per irrep\n",' ',
      Graph->subgr_per_irrep);
   fprintf(outfile,"%3c%2d Max strings in irrep\n", ' ', 
      Graph->max_str_per_irrep);
   fprintf(outfile,"%3c%2d Strings in total\n\n", ' ', Graph->num_str);

   fprintf(outfile, "\n");
   for (i=ras1_min; i<=ras1_max; i++) {
      for (j=0; j<=ras3_max; j++) {
         for (k=0; k<=ras4_max; k++) {
            if ((code = Graph->decode[i-ras1_min][j][k]) >= 0) {
               fprintf(outfile, "%5cDecode (%2d,%2d,%2d) = %3d\n",' ',
                  i,j,k,code);
               } 
            }
         }
      }

   fprintf(outfile, "\n%4cString Distinct Row Tables\n", ' ');
   fprintf(outfile, "%7c%3s %3s %3s %3s %3s %3s %3s %3s %3s %3s\n", ' ',
      "i", "j", "a", "b", "k0", "k1", "k0b", "k1b", "x", "y");
   for (i=0; i<Graph->nirreps; i++) {
      fprintf(outfile, "\n%4cIrrep %2d has %d strings\n", ' ', i, 
         Graph->str_per_irrep[i]);
      for (j=0; j<Graph->subgr_per_irrep; j++) {
         subgraph = Graph->sg[i] + j;
         if (subgraph->num_strings) {
            fprintf(outfile, "%6cCode(%3d) : %4d strings, offset = %4d\n",
               ' ', j, subgraph->num_strings, subgraph->offset);
            curr = subgraph->lvl;
            for (k=0; k<Graph->num_orb+1; k++,curr++) {
               for (l=0; l<curr->num_j; l++) {
                  a = (curr->a)[l];
                  b = (curr->b)[l];
                  fprintf(outfile, 
                     "%7c%3d %3d %3d %3d %3d %3d %3d %3d %3d %3d\n", ' ',
                      k, l+1, a, b, (curr->k)[0][l], (curr->k)[1][l],
                      (curr->kbar)[0][l], (curr->kbar)[1][l],
                      (curr->x)[l], (curr->y)[l]);
                  }
               }
            }
         }
      }

   fprintf(outfile, "\n");
   fflush(outfile);
}


