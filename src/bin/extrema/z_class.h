/*#############################################################################

  z_class.h

  declaration of the z-matrix derived class

  ###########################################################################*/

class z_class : public coord_base<simple_internal> {
  double **B_mat;
  public:
    z_class(int num_coord);
};



