/*###########################################################################*/
/*! \file math_tools.h
  \brief Math tools concrete class header file. */

/*! \class math_tools
   \brief A general mathematics routine class.

   General implementations of basic mathematical algorithms.  The only 
   dependency other than standard libraries is the <b>PSI 3.0</b>
   library <b>libciomr</b>. */
/*                                                  Joseph P. Kenny 11/29/01
  ##########################################################################*/

class math_tools{

  protected:
    math_tools() { return; }
    ~math_tools() { return; }

    double* newton_step(int dim, double **_Hi, double *g);
    double** update_bfgs(int dim, double *_var_dif, 
			 double *grad_dif, double **_Hi_old);
    double** update_ms(int dim, double *_var_dif, 
		       double *_grad_dif, double **_Hi_old);
};

/*! \fn math_tools::math_tools()
  \brief Default constructor, does nothing. */
/*! \fn math_tools::~math_tools()
  \brief Default destructor, does nothing. */
