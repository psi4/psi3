/*##################################################################################
  simple.h

  simple coordinate concrete class declaration and definition
  a type which coord_base can hold an array of

  we need a consistent definition for the internal type:
     bond length = 0
     valence angle = 1;
     torsional angle = 2;
  ################################################################################*/

class simple {
    int type;
    double val;
    int atom, bond, angle, tors;
  public:
     void set_simple(int ty, double value, int at, int bd, int an, int tr) {
	 type=ty; val=value; atom=at; bond=bd; angle=an; tors=tr;
	return;
    }    
    void set_val(double value) { val = value; return; }
    double get_val() { return val; }
    int get_type() { return type; }
    int get_atom() { return atom; }
    int get_bond() { return bond; }
    int get_angle() { 
	if(type>0)
	    return angle;
	else
	    punt("class simple error: non-existent angle reference atom asked for");
    }
    int get_tors() {
	if(type>1)
	    return tors;
	else
	    punt("class simple error: non-existent torsion reference atom asked for");
	    }  
};

 

























































































