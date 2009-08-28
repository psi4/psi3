/*! \file
    \ingroup OPTKING
    \brief Class declaration for simple internal coordinate set
*/

#ifndef _psi3_bin_optking_simples_h_
#define _psi3_bin_optking_simples_h_

#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <cov_radii.h>
#include <exception>

#include "stre.h"
#include "bend.h"
#include "tors.h"
#include "out.h"
#include "linb.h"
#include "frag.h"

#include <cctype>
#include <vector>

namespace psi { namespace optking {

class salc_set; // forward declaration for friend functions

using std::vector;

class simples_class {

  vector<stre_class> stre;
  vector<bend_class> bend;
  vector<tors_class> tors;
  vector<out_class> out;
  vector<linb_class> linb;
  vector<frag_class> frag;

  public:

   friend double **compute_B(const simples_class &, const salc_set &);
   friend double *compute_q(const simples_class &, const salc_set &);
   friend void empirical_H(const simples_class &, const salc_set &, const cartesians &);
   friend void delocalize(const simples_class &, const cartesians &);
   friend void rm_rotations(const simples_class & simples, const cartesians & carts,
       int &num_nonzero, double **evects);
   friend void get_syminfo(const simples_class & );
   friend int opt_step(cartesians &carts, simples_class &simples, const salc_set &symm);
   friend int *read_constraints(const simples_class & simples);
   friend void step_limit(const simples_class & simples, const salc_set &symm, double *dq);
   friend void check_zero_angles(const simples_class & simples, const salc_set & symm, double *dq);

   // constructor in frag.cc
   // user_intcos = 1; read in simple coordinates from intco.dat
   //               0; generate simples from geometry
   simples_class(cartesians& carts, int user_intcos);

   simples_class() : stre(), bend(), tors(), out(), linb(), frag() { };

   ~simples_class() {
     stre.clear();
     bend.clear();
     tors.clear();
     out.clear();
     linb.clear();
     frag.clear();
   }

   // print_vals: false => print definitions to a file in intco.dat format
   //             true  => print intcos and their values
   // print_frag_weights: whether to print weights for fragments
   void print(FILE *fp_out, bool print_vals, bool print_frag_weights = false) const;

   // print s vectors
   void print_s(FILE *fp_out) const {
     int i;
     fprintf(fp_out,"S vectors for simple internal coordinates\n");
     for (i=0; i<stre.size(); ++i) stre.at(i).print_s(fp_out);
     for (i=0; i<bend.size(); ++i) bend.at(i).print_s(fp_out);
     for (i=0; i<tors.size(); ++i) tors.at(i).print_s(fp_out);
     for (i=0; i<out.size(); ++i)  out.at(i).print_s(fp_out);
     for (i=0; i<linb.size(); ++i) linb.at(i).print_s(fp_out);
     for (i=0; i<frag.size(); ++i) frag.at(i).print_s(fp_out);
     fprintf(fp_out,"\n");
   }

   // compute values of simple internal coordinates
   void compute(double *geom) {
     int i;
     for (i=0; i<stre.size(); ++i) stre.at(i).compute(geom);
     for (i=0; i<bend.size(); ++i) bend.at(i).compute(geom);
     for (i=0; i<tors.size(); ++i) tors.at(i).compute(geom);
     for (i=0; i<out.size(); ++i)  out.at(i).compute(geom);
     for (i=0; i<linb.size(); ++i) linb.at(i).compute(geom);
     for (i=0; i<frag.size(); ++i) frag.at(i).compute(geom);
     return;
   }

   // compute and store s sectors 
   void compute_s(double *geom) {
     int i;
     for (i=0; i<stre.size(); ++i) stre.at(i).compute_s(geom);
     for (i=0; i<bend.size(); ++i) bend.at(i).compute_s(geom);
     for (i=0; i<tors.size(); ++i) tors.at(i).compute_s(geom);
     for (i=0; i<out.size(); ++i)  out.at(i).compute_s(geom);
     for (i=0; i<linb.size(); ++i) linb.at(i).compute_s(geom);
     for (i=0; i<frag.size(); ++i) frag.at(i).compute_s(geom);
     return;
   }

   // get number of simple internal coordinates
   int get_num(void) const {
     int i, n;
     n = stre.size() + bend.size() + tors.size() + out.size() + linb.size();
     for (i=0; i<frag.size(); ++i)
       n += frag[i].get_dim();
     return n;
   }

   void fix_near_180(void) {
     int i;
     for (i=0; i<tors.size(); ++i)
       tors[i].fix_near_180();
     for (i=0; i<frag.size(); ++i)
       frag[i].fix_near_180();
   }

   // ** following functions in frag.cc **

   // given id number, return type of coordinate, index within the type, and the
   // sub_index2 (I=1-6) for interfragment coordinates
   void locate_id(int id, Intco_type *itype, int *sub_index, int *sub_index2) const;

   // returns double ** bond connectivity matrix
   double **bond_connectivity_matrix(int natoms) const;

   // given the overall optking index (0-N), returns the user-assigned id number
   // optking index skips over interfragment coordinates turned "off"
   int index_to_id(int index) const;

   // given the id number, returns the overall optking index (0-N)
   int id_to_index(int id) const;

   // return id number that cooresponds to given simple type and atom numbers
   int get_id_from_atoms_stre(int a, int b) const;
   int get_id_from_atoms_bend(int a, int b, int c) const;
   int get_id_from_atoms_tors(int a, int b, int c, int d) const;
   int get_id_from_atoms_out(int a, int b, int c, int d, int *sign) const;
   int get_id_from_atoms_linb(int a, int b, int c, int linval) const;
   int get_id_from_atoms_frag(int a_natom, int b_natom, int *a_atom, int *b_atom) const;

   //int *atom2fragment(int natom) { return frag.atom2fragment(natom); } 

   bool is_unique (tors_class & t1) const {
     int i;
     bool unique = true;

     for (i=0; i<tors.size(); ++i)
       if (t1 == tors[i])
         unique = false;
     return unique;
   }


};

}} /* namespace psi::optking */

#endif
