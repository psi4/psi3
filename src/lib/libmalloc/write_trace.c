#ifdef	MALLOCTRACE


/* NOT copyright by SoftQuad. - msb, 1988 */
#ifndef lint
static char *SQ_SccsId = "@(#)write_trace.c	1.4 88/08/24";
#endif
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  _mal_write_trace				See ./ORIGINS for credits
 *
 * 	_mal_write_trace writes a line to the malloc trace file.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include	<stdio.h>
#include	"malloc.h"

/* Machine-dependent stuff:  the Sun's stack format. */

struct	stackframe {
	struct stackframe	*link;
	char			*ret_addr;
};

typedef struct stackframe	*FRAMEPTR;

#define	MYFRAME(myfirstarg)	 ((FRAMEPTR) (((Size *) &(myfirstarg)) - 2))

#define	NEXTFRAME(myframe)	((myframe)->link)
#define	RET_ADDR(myframe)	((myframe)->ret_addr)

/* Test whether a stdio file is still writable */

#define	STDIO_WRITABLE(fp)	(((fp)->_flag & _IOWRT) != NULL)


void
_mal_write_trace (event_label, req_nbytes, act_nbytes, mem)
char *event_label;
Size req_nbytes, act_nbytes;
char *mem;
{
	extern enum _malstate _malstate;
	extern FILE *_malfp;
	extern int _malbuff;
	FRAMEPTR frame;

	if (_malstate != S_TRACING || !STDIO_WRITABLE(_malfp))
		return;

	/* Reset the state in case stdio does any mallocing */
	_malstate = S_IN_STDIO;

	if (req_nbytes == act_nbytes)
		fprintf (_malfp, "%s of %ld at %ld\n",
			event_label, (long) req_nbytes, (long) mem);
	else
		fprintf (_malfp, "%s of %ld gets %ld at %ld\n",
			event_label, (long) req_nbytes, (long) act_nbytes,
								(long) mem);

	for (frame = NEXTFRAME (MYFRAME (event_label));
			frame;
			frame = NEXTFRAME (frame))
			
		fprintf (_malfp, "caller %8.8lx\n", (long) RET_ADDR (frame));

	fprintf (_malfp, "\n");

	/* Flush if desired */
	if (!_malbuff) fflush (_malfp);

	/* And back to normal */
	_malstate = S_TRACING;
}
#endif
