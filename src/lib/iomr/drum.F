*Deck Drum
      subroutine drum
#if defined(SUN4)
      call abrupt_underflow()
#elif defined(SUN3)
      call abrupt_underflow()
#endif
c
c     close (unit=5,err=1)
c   1 continue
c     close (unit=6,err=2)
c   2 continue
c     open (unit=5,file='input.dat',status='old')
c     rewind 5
c     open (unit=6,file='output.dat',status='unknown')
c     call wind(6)
c     open (unit=3,file='tape3.dat',status='unknown')
c     call wind(3)
      call ffile(5,'input',0)
      call ffile(6,'output',1)
      call ffile(3,'tape3',0)
c
      return
      end
