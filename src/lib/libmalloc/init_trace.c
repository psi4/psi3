#ifdef	MALLOCTRACE

#include <stdio.h>
#include "malloc.h"


/* NOT copyright by SoftQuad. - msb, 1988 */
#ifndef lint
static char *SQ_SccsId = "@(#)init_trace.c	1.3 88/08/24";
#endif
/*
 * Open the trace output file.  The _malstate serves to prevent
 * any infinite recursion if stdio itself calls malloc.
 */

void
_mal_init_trace()
{
	char *filename, *bufflag;
	char *getenv();
	extern enum _malstate _malstate;
	extern int _malbuff;
	extern FILE *_malfp;

	if (_malstate != S_INITIAL)	/* Can't happen, but what the heck */
		return;

	_malstate = S_IN_STDIO;

	filename = getenv (TRACEENVVAR);
	if (filename == NULL)
		filename = TRACEFILE;

	_malfp = fopen (filename, "w");
	if (_malfp == NULL) {

		perror (filename);
		fprintf(stderr, "(will run without malloc tracing)\n");

	} else {
		/*
		 * Nonportable kludge that should trigger a malloc
		 * in stdio while doing no harm, and also reduce the
		 * likelihood of a future malloc in stdio.
		 */

		putc (0, _malfp);
		fflush (_malfp);
		fseek (_malfp, 0L, 0);

		/* Were we requested not to flush after each entry? */

		bufflag = getenv (TRACEBUFVAR);
		_malbuff = (bufflag != NULL && *bufflag != '\0');

		/* Initialization successful */

		_malstate = S_TRACING;
	}
}
#endif
