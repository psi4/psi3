      SUBROUTINE MABORT                                                 
#ifdef APOLLO
      integer fault
#endif
C                                                                       
clj   call lnkerr('MAbort called')
C                                                                       
      write(6,*) 'mabort: ERROR, aborting'
      write(*,*) 'mabort: ERROR, aborting'
#ifdef AIX
      write(0,*) 'mabort: ERROR, aborting'
#endif
clj      call exit()
#ifdef APOLLO
      fault = 16#00120010
      call pfm_$error_trap(fault)
#endif
#ifndef NO_EXIT_ROUTINE
      call exit(3)
#endif
      stop 876
C                                                                       
      END                                                               
