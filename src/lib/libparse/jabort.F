c If an application does not have a jabort then this will be linked in as
c the default
      subroutine jabort
      write(6,*) 'no job abort sequence, skipping'
      return
      end
