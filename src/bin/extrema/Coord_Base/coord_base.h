/*###########################################################################*/
/*! \file coord_base.h
  \brief Coordinate base class declaration. */

/*! \class coord_base
  \brief First level of abstract coordinate classes.

  The <b>coord_base</b> class contains data and functions common to 
  all coordinate types.  All coordinate types derive from this class.
  Member data includes generic coordinates and gradients, cartesians and 
  cartesian gradients, a generic Hessian, information form previous iterations
  ("old" variables), and basic user suppied parameters.  Generic functions
  for coordinate data manipulations are members of this class.

  "generic" variables hold values and no information regarding
  coordinates to which they correspond.  Classes deriving from this 
  class determine the actual coordinate type and are responsible for 
  proper handling of these variables.  */
/*  						Joseph P. Kenny 11/29/01
  ###########################################################################*/

class coord_base : public math_tools {
  
  protected:

    int iteration, /*!< current iteration */ 
	num_atoms, /*!< number of atoms */
	num_entries, /*!< number of entries (including dummy atoms) */
	num_coords, /*!< number of coordinates which are actually optimized */
	grad_max, /*!< max allowable gradient is 10^-(grad_max) */
	print_lvl; /*!< print level */

    double *coords, /*!< generic coordinate array */ 
	*carts, /*!< cartesian coordinate array */
	*grads, /*!< generic gradient array */
	*c_grads, /*!< cartesian gradient array */
	**Hi, /*!< generic inverse hessian matrix */
	**u, /*!< the u matrix */
	*masses, /*!< atomic masses array */
	*coords_old, /*!< generic coordinates from previous iteration */
	*grads_old, /*!< generic gradients from previous iteration */
	**Hi_old, /*!< generic hessian inverse */
	*coord_write; /*!< holds coordinate values prior to optimization step
			until opt.dat is written */

    /*! \note "generic" variables hold values and no information regarding
       coordinates to which they correspond.  Classes deriving from this 
       class determine the actual coordinate type and are responsible for 
       proper handling of these variables */ 

    char *update, /*!< the hessian inverse update method */
	**felement; /*!< the full list of element names 
		      (including dummy atoms) */

  public:

    coord_base() : math_tools() {return;}
    void construct_coord_base(int natm, int nent, int nopt);
    ~coord_base(){

	free(carts);
	free(c_grads);
	free(coords);
	free(grads);
	free_matrix(Hi, num_coords);
	free_matrix(u,3*num_entries);
	free(masses);
	free(coords_old);
	free(grads_old);
	free_matrix(Hi_old,num_coords);
	free(coord_write);
	int i;
        for(i=0;i<num_entries;++i)
	    free(felement[i]);

	ip_done();
	tstop(outfile);
	fclose(infile);
	fclose(outfile);

	return;
    }
    void parse_input();
    void print_carts(double conv);
    void print_c_grads();
    void print_Hi();
    void read_opt();
    void write_opt();
    void update_Hi();
    void read_file11();
    virtual void write_file30();
    void grad_test();
    virtual void initial_Hi()=0;
    void H_test();
};

/*! \fn coord_base::coord_base()
  \brief Dummy constructor.  
  
  This is a dummy constructor which exists to make compilers happy.  It does
  nothing.  The actual constructor may depend on data which is not know
  when a derived class is initialized and may be called later in that
  class's constructor. */

/*! \fn coord_base::~coord_base()
  \brief Destructor.  Memory is freed and io is stopped */

/*! \fn coord_base::initial_Hi()
  \brief Derived classes must provide a inverse Hessian guess function. */

/*! \fn coord_base::write_file30()
  \brief Derived classes may override this function and write necessary 
  information to file30 */
