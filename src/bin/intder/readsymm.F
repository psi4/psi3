      subroutine rdsymm(lbls,iu,u,nsym,nasym,ns)
      implicit integer(a-z)
#include <error.h>
      character*80 lbls(nsym)
      real*8 u(ns,100),cff
      dimension iu(ns,0:100)

      errcod = frdpos(':INTCO:SYMM')
      if(errcod.ne.EOK) goto 901

      do 115  k=1,nsym
 115  iu(k,0)=0

      errcod = frdlpn()
      do 10 symco=1,nsym-nasym
        errcod = frdlpn()
        errcod = frdchr(lbls(symco))
        if(errcod.ne.EOK) goto 901
        errcod = frdlpn()
        if(errcod.ne.EOK) goto 901
        nsim=0
 201    errcod = frdint(simno)
        if(errcod.ne.Eevect) then
          nsim=nsim+1
          if(simno.lt.0) then
            iu(symco,nsim) = -simno
            u(symco,nsim) = -1.0d0
          else
            iu(symco,nsim) = simno
            u(symco,nsim) = 1.0d0
          endif
          goto 201
        endif
        errcod = frdlpn()
        if(errcod.eq.EOK) then
          do 100 sim=1,nsim
            errcod = frdflt(cff)
            u(symco,sim)=u(symco,sim)*cff
 100      continue
          errcod=frdrpn()
          errcod=frdrpn()
        endif
        iu(symco,0)=nsim
  10  continue

      errcod = frdpos(':INTCO:ASYMM')

      errcod = frdlpn()
      do 20 symco=nsym-nasym+1,nsym
        errcod = frdlpn()
        errcod = frdchr(lbls(symco))
        if(errcod.ne.EOK) goto 902
        errcod = frdlpn()
        if(errcod.ne.EOK) goto 902
        nsim=0
 221    errcod = frdint(simno)
        if(errcod.ne.Eevect) then
          nsim=nsim+1
          if(simno.lt.0) then
            iu(symco,nsim) = -simno
            u(symco,nsim) = -1.0d0
          else
            iu(symco,nsim) = simno
            u(symco,nsim) = 1.0d0
          endif
          goto 221
        endif
        errcod = frdlpn()
        if(errcod.eq.EOK) then
          do 200 sim=1,nsim
            errcod = frdflt(cff)
            u(symco,sim)=u(symco,sim)*cff
 200      continue
          errcod=frdrpn()
          errcod=frdrpn()
        endif
        iu(symco,0)=nsim
  20  continue
      return 

 901  write(6,*) 'error parsing SYMM vector'
      call pabort
      call qabort
 902  write(6,*) 'error parsing SYMM vector'
      call pabort
      call qabort
      end
