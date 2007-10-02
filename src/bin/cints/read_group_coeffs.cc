/*! \file read_group_coeffs.cc
    \ingroup (CINTS)
    \brief This function reads in a matrix of coefficients 
    for the fragments in a set of groups as defined by the GROUPS keyword.
*/
/*! \external \ingroup (OPTIONS)
  \param GROUP_COEFFS A matrix of non-negative doubles.
  Each row has to add up to one.
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include "global.h"
#include "defines.h"

#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>

#include "data_structs.h"

using namespace std;

namespace psi { namespace CINTS {
  
  //! Read the fragment coefficients in each group.
  void read_groups_coeffs()
  {
    int i, j, errcod, ngroups,nfrags;
    vector<vector<double> > group_coeffs;
    
    try {
      if( ip_exist("GROUP_COEFFS",0) ) {
	ip_count("GROUP_COEFFS", &ngroups, 0) ;
	
	group_coeffs.resize(ngroups);
	
	for(i=0;i<ngroups;i++) {
	  ip_count("GROUP_COEFFS",&nfrags,1,i);
	  if (errcod != IPE_OK) {
	    throw domain_error("Problem reading the GROUP_COEFFS matrix.");
	  }
	  group_coeffs[i].resize(nfrags);
	  for(j=0;j<nfrags;j++)
	    ip_data("GROUP_COEFFS","%lf",&(group_coeffs[i][j]),2,i,j);
	  double sum=0;
	  for(j=0;j<nfrags;j++) 
	    {
	      if(group_coeffs[i][j]<0.0)
		throw domain_error(" negative fragment coefficient.");
	      sum+=group_coeffs[i][j];
	    }
	  if(fabs(1.0-sum)>1e-16)
	    throw domain_error("fragment does not add to one.");
	}
      }
      /* IF USER DOES NOT SPECIFY GROUPS, POINT TO DEFAULT GROUPS */
      else {
	group_coeffs.resize(Molecule.num_atoms);
	for(i=0;i<Molecule.num_atoms;i++)
	  group_coeffs[i].push_back(1.0);
      }
      
      //      Molecule_cpp.group_coeffs=group_coeffs;
      
    } catch (domain_error e) {
      cerr << e.what() << endl;
      throw domain_error("error in read_group_coeffs");
    }
    return;
  }
};

};
