/*!
** \file get_file_ptr.c
*/

#include "includes.h"
#include "pointers.h"
 
/*!
** GET_FILE_PTR():  Simply return the bytewise global file pointer for the
** specified unit.  Allows one to get this pointer without having to go to
** the trouble of declaring things like ptr, etc.
**
**   \param unit = file unit number
**
** Returns:
**   bytewise file pointer
**
** Daniel, January 1996
*/
PSI_FPTR get_file_ptr(int unit)
{
  return(ptr.wptr[unit]);
}
