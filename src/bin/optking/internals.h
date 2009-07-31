/*! \file
    \ingroup OPTKING
    \brief Class declaration for internals
*/

#ifndef _psi3_bin_optking_internals_h_
#define _psi3_bin_optking_internals_h_

#include "stretch.h"
#include "bend.h"
#include "torsion.h"
#include "out_of_plane.h"
#include "lin_bend.h"
#include "fragment.h"

namespace psi { namespace optking {

class internals {

  public:
   
   stretch_set stre;
   bend_set bend;
   torsion_set tors;
   out_set out;
   lin_bend_set lin_bend;
   fragment_set frag;
  
    void compute_internals(int natom, double *geom);
    void compute_s(int natom, double *geom);
    void print_s();
    int get_num() {
      int num = 0;
      num += stre.get_num();
      num += bend.get_num();
      num += tors.get_num();
      num += out.get_num();
      num += lin_bend.get_num();
      num += frag.get_dim();
      return num;
    }

    ~internals();

   // user_intcos = 1; read in internals
   //               0; generate internals
    internals(cartesians& carts, int user_intcos);
    internals(int *num_arr);

   // print_flag: 0   print to a file in intco.dat format
   //             1   print intcos and their values
    void print(FILE *outfile, int print_flag);

    void locate_id(int id, int *intco_type, int *sub_index, int *sub_index2);
    int index_to_id(int index);
    int id_to_index(int id);
    void fix_near_180(void) {
      tors.fix_near_180();
      frag.fix_near_180();
    }
    double **bond_connectivity_matrix(int natom) const {
      return stre.bond_connectivity_matrix(natom);
    }
    int *atom2fragment(int natom) {
      return frag.atom2fragment(natom);
    } 

};

}} /* namespace psi::optking */

#endif
