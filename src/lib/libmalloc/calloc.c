#include "malloc.h"
#define NULL 0


/* NOT copyright by SoftQuad. - msb, 1988 */
#ifndef lint
static char *SQ_SccsId = "@(#)calloc.c	1.2 88/08/24";
#endif
/*
** And added by ado. . .
*/

char *
calloc(n, s)
unsigned	n;
unsigned	s;
{
	unsigned	cnt;
	char *		cp;

	cnt = n * s;
	cp = malloc(cnt);
	if (cp != NULL)
		bzero(cp, (int) cnt);
	return cp;
}

void
cfree(mem)
char *	mem;
{
	free(mem);
}

