/*! \file iwl_tebuf.h
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/

void iwl_buf_wrt_struct_nocut(struct iwlbuf *Buf, struct tebuf *Tebuf, int size);
void iwl_buf_wrt_struct(struct iwlbuf *Buf, struct tebuf *Tebuf, int size, double cutoff);

