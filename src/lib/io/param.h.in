/* Change this macro to unsigned long long int on systems that allow it */
#define PSI_FPTR unsigned long int

#define MAX_UNIT 300
#define MAX_VOLUME 8
#define MAX_STRING 512

#define PSIIO_UNOPENED 0
#define PSIIO_SEQUENTIAL 1
#define PSIIO_S_ASYNC 2
#define PSIIO_R_ASYNC 3
#define PSIIO_RAM 4
#define PSIIO_FORTRAN 5

#define IOOP_READ 1
#define IOOP_WRITE 2

#if defined(AIX)
# define C_FFRDBO ffrdbo
# define C_FFRDI ffrdi
# define C_FFRDC ffrdc
# define C_CPRGID cprgid
# define C_IOFO iofo
# define C_IOFC iofc
# define C_IOFW iofw
# define C_IOFR iofr
# define C_IOOPEN ioopen
# define C_IOCLOS ioclos
# define C_IOWRR  iowrr
# define C_IORDR  iordr
# define C_LOCATE locate
# define C_IOGTLN iogtln
# define C_MABORT mabort
#else
# define C_FFRDBO ffrdbo_
# define C_FFRDI ffrdi_
# define C_FFRDC ffrdc_
# define C_CPRGID cprgid_
# define C_IOFO iofo_
# define C_IOFC iofc_
# define C_IOFW iofw_
# define C_IOFR iofr_
# define C_IOOPEN ioopen_
# define C_IOCLOS ioclos_
# define C_IOWRR  iowrr_
# define C_IORDR  iordr_
# define C_LOCATE locate_
# define C_IOGTLN iogtln_
# define C_MABORT mabort_
#endif
