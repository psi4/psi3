/*!
** \file rclose.c
** \ingroup (CIOMR)
*/

extern void ioclos_(int *, int *);

/*!
** rclose: close a binary file
**
** \param unit = file number
** \param status = 3 to keep file, 4 to erase it
**
** \ingroup (CIOMR)
*/
void rclose(int unit, int status)
{
      ioclos_(&unit,&status);
   }
