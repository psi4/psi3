/*###########################################################################*/
/*! \file defines.h
  \brief #define parameters */
/*###########################################################################*/

#define NORMAL_PRINT 1 /*!< normal print level */
#define RIDICULOUS_PRINT 3 /*!< level at which printing becomes ridiculous */
#define BT_CONV 1.0e-10 /*!< convergence level for internals->cartesians
			  iterative back transformation */
#define POS_NEG_TORS 1.0e-6 /*!< tolerance for pos/neg torsion pairs */
#define BT_LOOP 10 /*!< max iterations for internals->cartesians 
		     iterative back transformation */
#define MAX_LINELENGTH 133 /*!< max length of lines in file11 */
#define EQUIV_GRAD 1.0e-6 /*!< tolerance for equivalent gradients */
#define ALMOST_ONE 1.0e-8 /*!< tolerance for numbers near 1.0 in 
			    cart_to_internals */
#define BOND_LIM 0.1 /*!< bond limit in angstroms */
#define ANGLE_LIM 5.0 /*!< angle limit in degrees */
