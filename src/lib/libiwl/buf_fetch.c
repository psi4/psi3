#include <stdio.h>
#include <psio.h>
#include "iwl.h"

/*
** iwl_buf_fetch()
**
** Fetch an IWL buffer from disk
** David Sherrill, 26 June 1996
**
*/
void iwl_buf_fetch(struct iwlbuf *Buf)
{
  psio_read(Buf->itap, IWL_KEY_BUF, (char *) &(Buf->lastbuf), sizeof(int),
	    Buf->bufpos, &Buf->bufpos);
  psio_read(Buf->itap, IWL_KEY_BUF, (char *) &(Buf->inbuf), sizeof(int),
	    Buf->bufpos, &Buf->bufpos);
  psio_read(Buf->itap, IWL_KEY_BUF, (char *) Buf->labels, Buf->ints_per_buf * 
	    4 * sizeof(Label), Buf->bufpos, &Buf->bufpos);
  psio_read(Buf->itap, IWL_KEY_BUF, (char *) Buf->values, Buf->ints_per_buf *
	    sizeof(Value), Buf->bufpos, &Buf->bufpos);
  Buf->idx = 0;
}
