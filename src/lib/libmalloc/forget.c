#include "malloc.h"
#define NULL 0

/* NOT copyright by SoftQuad. - msb, 1988 */
#ifndef lint
static char *SQ_SccsId = "@(#)forget.c	1.4 88/08/24";
#endif
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *  forget				William L. Sebok
 * A "smarter" malloc v1.0		Sept. 24, 1984 rev. June 30,1986
 *			Then modified by Arthur David Olsen
 *
 *	forget returns to the malloc arena all memory allocated by sbrk()
 *	 above "bpnt".
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
forget(bpnt)
char *bpnt;
{
	register struct overhead *p, *q, *r, *b;
	register Size l;
	struct overhead *crbrk;
	char pinvalid, oendfree;

	/*
	 * b = forget point
	 * p = beginning of entry
	 * q = end of entry, beginning of gap
	 * r = end of gap, beginning of next entry (or the break)
	 * pinvalid = true when groveling at forget point
	 */

	pinvalid = 0;
	oendfree = endfree;	endfree = 0;
	b = (struct overhead *)bpnt;
	p = FROMADJ(adjhead.q_back);
	r = crbrk = (struct overhead *)CURBRK;

	for (;pinvalid == 0 && b < r; p = FROMADJ(TOADJ(r = p)->q_back)) {
		if ( p == FROMADJ(&adjhead)
		 || (q = (struct overhead *)((char *)p + p->ov_length)) < b
		) {
			pinvalid = 1;
			q = b;
		}

		if (q == r)
			continue;

		ASSERT(q < r,
"\nforget: addresses in adjacency chain are out of order!\n");

		/* end of gap is at break */
		if (oendfree && r == crbrk) {
			(void)BRK((char *)q);	/* free it yourself */
			crbrk = r = q;
			continue;
		}

		if (pinvalid)
			q = (struct overhead *) /* align q pointer */
				(((long)q + (NALIGN-1)) & (~(NALIGN-1)));

		l = (char *)r - (char *)q;
		/*
		 * note: unless something is screwy: (l%NALIGN) == 0
		 * as r should be aligned by this point
		 */

		if (l >= (int) sizeof(struct overhead)) {
			/* construct busy entry and free it */
			q->ov_magic = MAGIC_BUSY;
			q->ov_length = l;
			insque(TOADJ(q),TOADJ(p));
			free((char *)q + sizeof(struct overhead));
		} else if (pinvalid == 0) {
			/* append it to previous entry */
			p->ov_length += l;
			if (p->ov_magic == MAGIC_FREE) {
				remque(TOBUK(p));
				{
					register struct qelem *	bp;

					bp = &buckets[mlindx(p->ov_length)];
					if (bp > hifreebp)
						hifreebp = bp;
					insque(TOBUK(p),bp);
				}
			}
		}
	}
	endfree = oendfree;
	if (endfree)
		mlfree_end();
	return;
}
