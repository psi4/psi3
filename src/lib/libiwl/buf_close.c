/*!
  \file buf_close.c
*/
#include <stdio.h>
#include <psio.h>
#include "iwl.h"

/*!
** IWL_BUF_CLOSE()
** 
**	\param Buf      Buffer to be closed
**	\param keep    Do not delete if keep==1
**
** Close a Integrals With Labels Buffer
*/
void iwl_buf_close(struct iwlbuf *Buf, int keep)
{

   psio_close(Buf->itap, keep ? 1 : 0);
   free(Buf->labels);
   free(Buf->values);
}

