/*##################################################################################
  simple_internal.h

  simple coordinate concrete class declaration and definition
  a type which coord_base can hold an array of

  we need a consistent definition for the internal type:
     bond length = 0
     valence angle = 1;
     torsional angle = 2;
  ################################################################################*/

class simple_internal {
    int type;
    double val;
  public:
    void simple_internal::simple_internal(double value, int A, int B) {
    void set(int i) { type = i; return;}
    void set(double a) { val = a; return;}
    void set(int i, double a) { type = i; val = a; return;}
    void set(double a, int i) { type = i; val = a; return;}
    int get_type() { return type;}
    double get_val() { return val;}
};

    const class simple {
	int type;
	mutable double val;
	int *ref;
      public:
	void simple::simple(double value, int atom, int bond); //bond constructor
        void simple::simple(double value, int atom, int bond, int angle); //valence angle contructor
        void set_val(double value) { val = value; return; }
        double get_val() { return val; }
        int get_atom() { return ref[0]; }
        int get_bond() { return ref[1]; }
        int get_angle() { 
	    if(type>0)
		return ref[2];
            else
                punt("class simple error: non-existent angle reference atom asked for");
	}
	int get_tors() {
	    if(type>1)
		return ref[3];
            else
                punt("class simple error: non-existent torsion reference atom asked for");
	}     
    };

    //bond constructor
    void simple::simple(double value, int atom, int bond) {  
	try {
	    ref = new int[2];
	}
	catch(bad_alloc) {
	    punt("malloc error: no memory left");
	}
	val = value;
	ref[0] = atom;
	ref[1] = bond;
	return;
    }

    //angle constructor
    void simple::simple(double value, int atom, int bond, int angle) {
	try {
	    ref = new int[3];
	}
	catch(bad_alloc) {
	    punt("malloc error: no memory left");
	}
	val = value;
	ref[0] = atom;
	ref[1] = bond;
        ref[2] = angle;
	return;
    }    
    
