/*#############################################################################

  z_class.h

  declaration of the z-matrix derived class

  ###########################################################################*/

#include <file30.h>

class z_class : public coord_base<simple_internal> {
  public:
    struct z_entry *z_geom;
    z_class(int num_coord);
    ~z_class();
    int get_bond_atom(int i){ return z_geom[i].bond_atom;};
    int get_angle_atom(int i){ return z_geom[i].angle_atom;};
    int get_tors_atom(int i){ return z_geom[i].tors_atom;};
    void set_bond(int i, double val){ 
	z_geom[i].bond_val = val;
	return;}
    void set_angle(int i, double val){
	z_geom[i].angle_val = val;
	return;}
    void set_tors(int i, double val){
	z_geom[i].tors_val = val;
	return;}
};



