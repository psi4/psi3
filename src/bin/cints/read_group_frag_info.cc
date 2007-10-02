/*! \file read_group_frag_info.cc
    \ingroup (CINTS)
    \brief Read the matrix of group and fragment assignments.
*/
/*! \ingroup (OPTIONS)@{
  GROUPS: Matrix of group and fragment memberships for probabilistic combination of atomic potentials (formerly LCAP).
  The first integer denotes the group
  and the second denotes the fragment within the group.@}
*/
#define EXTERN

#include <cstdio>
#include <cstdlib>
#include "defines.h"
#include "global.h"
#include "data_structs.h"

#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <vector>
#include <stdexcept>

using namespace std;

namespace psi {
  namespace CINTS {
    /*! Read in the group and fragment memberships of each atom.*/
    void read_group_frag_info()
    {
      int i, j, errcod;
      int atomcount;
      double Z = 0.0;
      int tmp = 0;
      int num_elem;
      vector<vector<int> > groups_and_fragments(Molecule.num_atoms);
      
      ip_count("GROUPS", &num_elem, 0);
    
      if(num_elem%2 != Molecule.num_atoms) 
	throw domain_error("Problem with number of elements in GROUPS. num_elem%2 != Molecule.num_atoms");
      
      if (num_elem == 0)
	throw domain_error("GROUPS is empty!");
      
      atomcount = 0;
      for(i=0;i<Molecule.num_atoms;i++) {
	for(j=0; j<2;j++){
	  errcod = ip_data("GROUPS","%ld", &tmp,1,2*i+j);
	  if (errcod != IPE_OK)
	    throw domain_error("Problem with the GROUPS array.");
	  else
	    groups_and_fragments[i].push_back(tmp);
	}
      }
      
      return;
    }
  }
}
