#ifdef APOLLO
      subroutine ioerr(ierr)
C   BASE.INS.FTN, /sys/ins, jrw, 07/08/81
C
C   CHANGES:
C     87/02/10  mishkin Added ios_$stdin, etc.
C     86/07/24  geiser IADDR is now an intrinsic.
C     85/03/26  RpS    Added lib_$member_of_set declaration
C     81/09/17  CAS    added declaration of IADDR
C
C   base include file for fortran customers
                                               
C---System defined stream id's:

       integer*2
     +     stream_$stdin, stream_$stdout, stream_$errin, stream_$errout

       parameter (
     +     stream_$stdin   = 0,
     +     stream_$stdout  = 1,
     +     stream_$errin   = 2,
     +     stream_$errout  = 3)

       integer*2
     +     ios_$stdin, ios_$stdout, ios_$errin, ios_$errout

       parameter (
     +     ios_$stdin   = 0,
     +     ios_$stdout  = 1,
     +     ios_$errin   = 2,
     +     ios_$errout  = 3)


C---Status "ok" definition:

       integer*4
     +     status_$ok

       parameter (
     +     status_$ok = 0)


C---Functions:

       logical
     +     lib_$member_of_set
C   ERROR.INS.FTN, us>ins, jrw, 06/15/81
C   error function declarations for fortran

C   Changes:
C      06/15/81 jrw  original coding.


       integer*2
     +     error_$subsys, error_$module, error_$code
       logical
     +     error_$fail

      call error_$print(ierr)
      return
      end
#else
      subroutine ioerr(ierr)
      write(6,*) 'ioerr: error ', ierr, ' encountered.'
      return
      end
#endif

