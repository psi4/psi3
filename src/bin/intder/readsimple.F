      subroutine rdsimp(type,ia,s,ns,nst,nbnd,ntor,nlin,nout,nspf)
      implicit integer(a-z)
#include <error.h>
      character*5 type(ns)
      character*80 keywrd
      real*8 s(ns)
      dimension ia(ns,6)

      if(nst.ne.0) then
        keywrd = ':INTCO:STRE'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 10 i=1,nst
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)=' STRE'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdint(ia(ij,6))
            if(errcod.ne.Eevect) errcod = frdrpn()
  10      continue
        endif
      endif

      if(nbnd.ne.0) then
        keywrd = ':INTCO:BEND'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 11 i=1,nbnd
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)=' BEND'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdint(ia(ij,3))
            errcod=frdint(ia(ij,6))
            if(errcod.ne.Eevect) errcod = frdrpn()
  11      continue
        endif
      endif

      if(ntor.ne.0) then
        keywrd = ':INTCO:TORS'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 12 i=1,ntor
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)=' TORS'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdint(ia(ij,3))
            errcod=frdint(ia(ij,4))
            errcod=frdint(ia(ij,6))
            if(errcod.ne.Eevect) errcod = frdrpn()
  12      continue
        endif
      endif

      if(nlin.ne.0) then
        keywrd = ':INTCO:LIN2'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 13 i=1,nlin
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)=' LIN1'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdint(ia(ij,3))
            errcod=frdint(ia(ij,4))
            errcod=frdint(ia(ij,6))
            if(errcod.ne.Eevect) errcod = frdrpn()
  13      continue
        endif
      endif

      keywrd = ':INTCO:OUT'
      errcod = frdpos(keywrd)
      if(errcod.eq.EOK) then
        errcod = frdlpn()
        if(errcod.ne.EOK) goto 901
        do 14 i=1,nout
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          errcod = frdint(ij)
          type(ij)='  OUT'
          errcod=frdint(ia(ij,1))
          errcod=frdint(ia(ij,2))
          errcod=frdint(ia(ij,3))
          errcod=frdint(ia(ij,4))
          errcod=frdint(ia(ij,6))
          if(errcod.ne.Eevect) errcod = frdrpn()
  14    continue
      endif

      if(nspf.ne.0) then
        keywrd = ':INTCO:SPF'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 15 i=1,nspf
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)='  SPF'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdflt(s(ij))
            if(errcod.eq.Eevect) goto 15
            errcod=frdint(ia(ij,6))
            if(errcod.ne.Eevect) errcod = frdrpn()
  15      continue
        endif
      endif

      return

  901 continue
      write(6,*) 'trouble parsing simple internal coordinates'
      call pabort
      call qabort
      end
