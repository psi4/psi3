/*
** GLOBALDEFS.H
**
** Global defines for DETCAS
**
** C. David Sherrill
** University of California, Berkeley
** April 1998
**
*/

#define MAX_RAS_SPACES 4
#define IOFF_MAX       10302
#define INDEX(i,j) ( (i>j) ? (ioff[(i)] + (j)): (ioff[(j)] + (i)) )
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define MAX_COMMENT 10
