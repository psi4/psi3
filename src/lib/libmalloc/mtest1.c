
/* NOT copyright by SoftQuad. - msb, 1988 */
#ifndef lint
static char *SQ_SccsId = "@(#)mtest1.c	1.2 88/08/24";
#endif
/*
 * tstmalloc	this routine tests and exercizes the malloc/free package
 */
#include <stdio.h>
#include <setjmp.h>
#include "malloc.h"

jmp_buf env;

main(argc,argv)
int argc; char *argv[];
{
	char lin[100], arg1[20], arg2[20], arg3[20];
	char *res, *malloc(), *realloc();
	register struct overhead *p, *q;
	register struct qelem *qp;
	int arg, nargs, argn;
	int l;
	void malerror();

	mlabort = &malerror;
	setjmp(env);

	for (;;) {
		printf("*  ");
		if (fgets(lin,sizeof lin, stdin)== NULL)
			exit(0);
		nargs = sscanf(lin,"%s%s%s",arg1,arg2,arg3);
		switch (arg1[0]) {

		case 'b':
			if (nargs == 2) {
				arg = atoi(arg2);
				if (arg<0 ||  arg>=NBUCKETS)
					goto bad;

				qp = &buckets[arg];
				printf("Bucket %2d\t\t\t  buk=%08lx %08lx\n",
					arg, qp->q_forw,qp->q_back);
				qp = qp->q_forw;
				for (; qp != &buckets[arg]; qp = qp->q_forw) {
					p = FROMBUK(qp);
					if (dump(p))
						break;
				}
			} else {
				printf("Buckets:");
				for (qp=buckets; qp<&buckets[NBUCKETS];qp++){
					if (qp->q_forw != qp)
						printf(" %d", qp-buckets);
				}
				printf("\n");
			}
			break;
		case 'e':
			endfree = 1;
			break;
		case 'E':
			endfree = 0;
			break;
		case 'f':
			if (nargs != 2) goto bad;
			sscanf(arg2,"%lx",&arg);
			printf("free(%x)\n",arg);
			free(arg);
			break;
		case 'F':
			if (nargs != 2) goto bad;
			sscanf(arg2,"%lx",&arg);
			printf("forget(%x)\n",arg);
			forget(arg);
			break;
		case 'h':
			printf("\
b	print bucket chains that aren't empty\n\
b [n]	trace through bucket chains\n\
e	turn on freeing of end of memory\n\
E	turn off freeing of end of memory\n\
f addr		free addr\n\
F addr		forget below addr\n\
h	print this help file\n\
m bytes		malloc bytes\n\
q	quit\n\
r addr bytes	realloc\n\
s	print break addr\n\
S	sbrk count\n\
t	trace through adjacency chain\n\
");
			break;
		case 'm':
			if (nargs != 2) goto bad;
			arg = atoi(arg2);
			res = malloc(arg);
			printf("malloc(%d) = %lx\n",arg, res);
			break;
		case 'r':
			if (nargs != 3) goto bad;
			sscanf(arg2,"%lx",&arg);
			argn = atoi(arg3);
			res = realloc(arg,argn);
			printf("realloc(%lx,%d) = %lx\n",arg,argn,res);
			break;
		case 'q':
			exit(0);
			break;
		case 's':
			printf("brk = %08x\n",sbrk(0));
			break;
		case 'S':
			if (nargs != 2) goto bad;
			sscanf(arg2,"%ld",&arg);
			printf("sbrk(%d)\n",arg);
			sbrk(arg);
			break;
		case 't':
			printf("\t\t\t\t\t\thead    adj=%08lx %08lx\n",
				adjhead.q_forw,adjhead.q_back);
			for (qp = adjhead.q_forw; qp!=&adjhead; qp=qp->q_forw) { 
				p = FROMADJ(qp);
				if (dump(p))
					break;
				q = FROMADJ(qp->q_forw);
				if (q==FROMADJ(&adjhead))
					q = (struct overhead *)sbrk(0);
				l = (char *)q - (char *)p - p->ov_length;
				if (l>0)
					printf("%08x free space  len=%8d\n",
						(char *)p + p->ov_length, l);
			}
			break;
		default:
	bad:		printf("Bad command\n");
		}
	}
}

dump(p)
register struct overhead *p;
{
	register char *s;
	int stat = 0;

	if (p->ov_magic == MAGIC_FREE)
		s = "MAGIC_FREE ";
	else if (p->ov_magic == MAGIC_BUSY) 
		s = "MAGIC_BUSY ";
	else {
		s = "BAD MAGIC  ";
		stat = 1;
	}
		
	printf( "%08x %s len=%8d buk=%08x %08x adj=%08x %08x\n",
		(&p[1]),s,p->ov_length,p->ov_buk.q_forw,p->ov_buk.q_back,
		p->ov_adj.q_forw,p->ov_adj.q_back
	);
	return(stat);
}

void
malerror()
{
	write(2,"malloc error\n",13);
	longjmp(env,1);
}
