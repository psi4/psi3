c The main abort sequence
      subroutine mabort
#if defined(SUN)
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif

      write(6,*) 'mabort: main abort sequence begins'
      write(6,*) '       This routine will try to print out all'
      write(6,*) '       available information at abort time.  Not'
      write(6,*) '       all this will be relevant.  Output from C'
      write(6,*) '       library abort sequence routines might not'
      write(6,*) '       appear in this file.  Check other files'
      write(6,*) '       produced by your calculation for other error'
      write(6,*) '       messages.'

      call jabort
      call pabort
      call ioabort

      write(6,*) 'mabort: main abort sequence done'

      call fioflu(6)

clj Try to clean up
      call psdone()

#if !defined(NO_EXIT_ROUTINE)
      call exit(1)
#else
      stop 1
#endif
      end

c----------------------------------------------------------------
c The quiet abort sequence
      subroutine qabort
#ifdef sun
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif

      call fioflu(6)

#if !defined(NO_EXIT_ROUTINE)
      call exit(1)
#else
      stop 1
#endif
      end
c
