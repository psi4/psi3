
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <sys/param.h>
#include <sys/mount.h>

#if (defined(DEC)||defined(SUN)||defined(MIPS)||defined(SGI))
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

