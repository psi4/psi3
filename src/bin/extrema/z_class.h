/*#############################################################################

  z_class.h

  declaration of the z-matrix derived class

  ###########################################################################*/

#include <file30.h>

class z_class : public coord_base<simple_internal> {
  double **B_mat;
  public:
    struct z_entry *z_geom;
    z_class(int num_coord);
    ~z_class();
};



