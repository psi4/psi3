#include "iomrparam.h"

#ifdef ALLOC_GLOBALS
#define EXTERN
#else
#define EXTERN extern
#endif

struct pointer {
     PSI_FPTR *wptr;
     };

EXTERN struct pointer ptr;
EXTERN int sector;
EXTERN time_t time_start, time_end;
#if !defined(SGI)
EXTERN struct tms total_tmstime;
#endif

#undef EXTERN
