c If an application does not have a ioabort then this will be linked in as
c the default
c This only be used if the io library is being used.
      subroutine ioabort
      write(6,*) 'no io abort sequence, skipping'
      return
      end
