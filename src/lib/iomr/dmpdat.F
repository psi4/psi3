C   *** (MEMBER:DMPDAT  ) CONVERTED IN POTOPS-00-00. CONV.DATE:89-03-01
      subroutine dmpdat(dat,ndat)
      implicit integer (a-z)
      integer dat(ndat)
      write(6,*) 'dump of data:'
      do 101 i = 1, ndat
        write(6,1) i, dat(i)
    1   format(1X,I7,1X,I15)
  101   continue
      return
      end
