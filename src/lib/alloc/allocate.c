
#include <stdio.h>

#if defined(APOLLO)
main()
#elif defined(SUN3)
MAIN_()
#elif defined(DEC)
MAIN__()
#elif defined(MIPS)
MAIN__()
#elif defined(AIX)
main()
#else
MAIN__()
#endif
{
	int maxcor,maxcri,sizeofreal=8;
	char *core,*malloc();
#if defined(APOLLO) || defined(AIX)
	fndcro(&maxcor); 
#else
	fndcro_(&maxcor); 
#endif

        if (!maxcor) {
#if defined(APOLLO) || defined(AIX)
	  fndcor(&maxcor,&maxcri,&sizeofreal);
#else
	  fndcor_(&maxcor,&maxcri,&sizeofreal);
#endif
          }
        else {
          fprintf(stdout,"NOTE: using old style memory request because\n");
          fprintf(stdout,"      # MAXCOR # or # CORE ### was found.\n");
          }

	core = malloc((unsigned)((maxcor)*sizeofreal));
	if (core == NULL) {
		perror("Fortran memory allocation");
		exit(1);
		}

#if defined(__i386__) && defined(__GNUC__)
        /* make floating point errors cause an exception (except for
         * denormalized operands, since small numbers are denormalized) */
        asm("fldcw %0" : : "o" (0x372));
#endif

#if defined(APOLLO) || defined(AIX)
	fentry(core, core, &maxcor);
#else
	fentry_(core, core, &maxcor);
#endif
	}
