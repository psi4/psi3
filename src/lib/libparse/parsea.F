c These routines provide a more efficient method to parse arrays.
c First the user calls frdpos to position the input.
c If need '(' call frdlpn
c If need ')' call frdrpn
c If need integer call frdint
c If need real call frdrl
c If need character call frdchr
c If the array dimensions are unknown then watch the return codes.
c This allows an array to be read in without parsing the entire input
c (frdpos - the time consumer - is called only once)
c---------------------------------------------------------------------------
c This expects the next token to be a character datum, if it is
c a ')' then Eevect is returned if is is a '(' then Ebvect is return.
      integer function frdchr(result)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) result

      character*(lentok) token
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer gtoken

      frdchr = gtoken(token)
      if (frdchr .ne. EOK) return

      if (token .eq. ')') then
        errcod = Eevect
        frdchr = errcod
        emess = 'Reached the end of an array (frdchr)'
        return
      else if (token .eq. '(') then
        errcod = Ebvect
        frdchr = errcod
        emess = 'Found the beginning of a vector (frdchr)'
        return
        endif

      result = token
      if (result .ne. token) then
        errcod = Etrunc
        frdchr = errcod
        emess = 'Character data was truncated (frdchr)'
        endif

      return
      end
c---------------------------------------------------------------------------
c This expects the next token to be a integer datum, if it is
c a ')' then Eevect is returned if is is a '(' then Ebvect is return.
      integer function frdint(result)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer result

      character*(lentok) token
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer sign
      integer gtoken, chrint

      frdint = gtoken(token)
      if (frdint .ne. EOK) return

      if (token .eq. ')') then
        errcod = Eevect
        frdint = errcod
        emess = 'Reached the end of an array (frdint)'
        return
      else if (token .eq. '(') then
        errcod = Ebvect
        frdint = errcod
        emess = 'Found the beginning of a vector (frdint)'
        return
        endif

      sign = 0
      if (token .eq. '-') then
        sign = -1
      else if (token .eq. '+') then
        sign = 1
        endif

      if (sign .ne. 0) then
        frdint = gtoken(token)
        if (frdint .ne. EOK) return

        if (token .eq. ')') then
          errcod = Eevect
          frdint = errcod
          emess = 'Reached the end of an array (frdint)'
          return
        else if (token .eq. '(') then
          errcod = Ebvect
          frdint = errcod
          emess = 'Found the beginning of a vector (frdint)'
          return
          endif
       endif

      frdint = chrint(token,result)
      if (frdint .ne. EOK) return
      if (sign.ne.0) result = result*sign

      return
      end
c---------------------------------------------------------------------------
c This expects the next token to be a floating point datum, if it is
c a ')' then Eevect is returned (end of vector) if is is a '(' then
c Ebvect is return.
      integer function frdflt(result)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      real*8 result

      character*(lentok) token
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer sign
      integer gtoken, chrflt

      frdflt = gtoken(token)
      if (frdflt .ne. EOK) return

      if (token .eq. ')') then
        errcod = Eevect
        frdflt = errcod
        emess = 'Reached the end of an array (frdflt)'
        return
      else if (token .eq. '(') then
        errcod = Ebvect
        frdflt = errcod
        emess = 'Found the beginning of a vector (frdflt)'
        return
        endif

      sign = 0
      if (token .eq. '-') then
        sign = -1
      else if (token .eq. '+') then
        sign = 1
        endif

      if (sign .ne. 0) then
        frdflt = gtoken(token)
        if (frdflt .ne. EOK) return

        if (token .eq. ')') then
          errcod = Eevect
          frdflt = errcod
          emess = 'Reached the end of an array (frdflt)'
          return
        else if (token .eq. '(') then
          errcod = Ebvect
          frdflt = errcod
          emess = 'Found the beginning of a vector (frdflt)'
          return
          endif
        endif

      frdflt = chrflt(token,result)
      if (frdflt .ne. EOK) return
      if (sign .ne. 0) result = result*sign

      return
      end
c---------------------------------------------------------------------------
c This expects the next token to be char, if it is not, Esyntx is returned.
      integer function fexpct(char)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) char

      character*(lentok) token
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer gtoken

      fexpct = gtoken(token)
      if (fexpct .ne. EOK) return

      if (token .ne. char) then
        errcod = Esyntx
        fexpct = errcod
        emess = 'The expected token, '''
        call scat(emess,char)
        call scat(emess,''', was not found (fexpct)')
        return
        endif

      fexpct = errcod
      return
      end
c---------------------------------------------------------------------------
c This expects the next token to be a '('
      integer function frdlpn()
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"

      integer fexpct

      frdlpn = fexpct('(')
      return
      end
c---------------------------------------------------------------------------
c This expects the next token to be a ')'
      integer function frdrpn()
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"

      integer fexpct

      frdrpn = fexpct(')')
      return
      end
