c Here are some character routines.
c---------------------------------------------------------------------------
c concat second string into first string
      subroutine scat(s1,s2)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer (a-z)
#endif
      character*(*) s1,s2

      integer i,next

      next = 1
      do 101 i = len(s1), 1, -1
        if (s1(i:i) .ne. ' ') then
          next = i + 1
          goto 901
          endif
 101    continue

 901  continue

      do 102 i = 1, len(s2)
        if (next .gt. len(s1)) goto 902
        s1(next:next) = s2(i:i)
        next = next + 1
 102    continue

 902  continue

      return
      end
c---------------------------------------------------------------------------
c concat second string into first string adding a blank between the two
      subroutine scatb(s1,s2)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer (a-z)
#endif
      character*(*) s1,s2

      integer i,next

c     write(6,*) ' s1,s2 = ', s1,s2
c     write(6,*) ' lengths = ', len(s1),len(s2)

      next = 1
      do 101 i = len(s1), 1, -1
        if (s1(i:i) .ne. ' ') then
          next = i + 2
          goto 901
          endif
 101    continue

 901  continue

      do 102 i = 1, len(s2)
c       write(6,*) 'next,i,len(s1)=',next,i,len(s1)
        if (next .gt. len(s1)) goto 902
        s1(next:next) = s2(i:i)
        next = next + 1
 102    continue

 902  continue

c     write(6,*) ' s1 = ',s1

      return
      end
c----------------------------------------------------------------------------
c Convert an integer into a character string
      subroutine intchr(int,chr)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer (a-z)
#endif
      integer int
      character*(*) chr

      integer pos,curint,rem,nxtint

      character*1 digchr

      curint = int
      pos = len(chr)

      chr = ' '

 801  continue
      if (pos .eq. 0) then
        write(6,*) 'ERROR: intchr: chr too small, aborting'
        call mabort
        endif
      nxtint = curint/10
      rem = curint - nxtint*10
      chr(pos:pos) = digchr(rem)
      curint = nxtint
      if (curint .eq. 0) goto 902
      pos = pos - 1
      goto 801

 902  continue
      call leftju(chr)
c     write(6,*) 'intchr: int,chr = ',int,chr
      return
      end
c------------------------------------------------------------------------------
c Left justify a string
      subroutine leftju(string)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer (a-z)
#endif
      character*(*) string

      integer i,pos

c Find out where the nonblank characters start
      do 101 i = 1, len(string)
        if (string(i:i) .ne. ' ') then
          pos = i
          goto 902
          endif
 101    continue
      return

 902  continue
c     write(6,*) 'leftju: pos = ',pos
c Copy the nonblank characters to the beginning of the string
      do 102 i = 1, len(string)
        if (pos .le. len(string)) then
          string(i:i) = string(pos:pos)
        else
          string(i:i) = ' '
          endif
        pos = pos + 1
 102    continue

      return
      end
c---------------------------------------------------------------------------
c Convert a digit to a character
      character*1 function digchr(dig)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer (a-z)
#endif
      integer dig

      character*10 digits
      data digits/'0123456789'/

      if ((dig .lt. 0).or.(dig .gt. 9)) then
        write(6,*) 'digchr: argument is out of range'
        call mabort
        endif

      digchr = digits(dig+1:dig+1)
      return
      end
c---------------------------------------------------------------------------
c If the regular expression re matches string then return true.
c NOTE: The only re accepted is '.'
      logical function remtch(re,string)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer (a-z)
#endif
      character*(*) re,string

      integer i, lrex, lstr, length

      lrex = len(re)
      lstr = len(string)

c If one of the strings is longer than the other then make sure the
c longer one has nothing past the end of the shorter one.
      if (lrex .gt. lstr) then
        if (re(lstr+1:lrex).ne.' ') then
          remtch = .false.
          return
          endif
        length = lstr
      else if (lstr .gt. lrex) then
        if (string(lrex+1:lstr).ne.' ') then
          remtch = .false.
          return
          endif
        length = lrex
      else
        length = lrex
        endif
c Now only need to work within the first length characters

      do 101 i = 1, length
        if (re(i:i) .eq. '.') then
          if (string(i:i) .eq. ' ') then
            remtch = .false.
            return
            endif
        else if (re(i:i) .ne. string(i:i)) then
          remtch = .false.
          return
          endif
 101    continue

      remtch = .true.
      return
      end
c---------------------------------------------------------------------------
c  Takes a character datum and converts it into a floating point number.
      integer function chrflt(char,result)
#ifdef UNDEFINED
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) char
      real*8 result

      logical num
      character*(lentok) cnum, cden
      integer i, start, sign, lden, lnum, ipos, digit
      real*8 fnum, fden
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer c1int, chrlen

      if (char(1:1).eq.'+') then
        sign = +1
        start = 2
      else if (char(1:1).eq.'-') then
        sign = -1
        start = 2
      else
        sign = +1
        start = 1
        endif

      ipos = 1
      num = .true.
      cnum = ' '
      cden = ' '
      do 101 i = start, len(char)
        if (char(i:i) .eq. '.') then
          num = .false.
          ipos = 1
        else if ((char(i:i).eq.'+').or.(char(i:i).eq.'-')) then
c         Signs have to be check for here because chrint will allow
c         signs in certain places.
          errcod = Etype
          chrflt = errcod
          emess = 'Bad floating point number (chrflt)'
          return
        else if (char(i:i) .eq. ' ') then
          if (char(i:len(char)) .ne. ' ') then
            errcod = Esyntx
            chrflt = errcod
            emess = 'Blank in the middle of the data (chrflt)'
            return
          else
            goto 901
            endif
        else if (num) then
          cnum(ipos:ipos) = char(i:i)
          ipos = ipos + 1
        else
          cden(ipos:ipos) = char(i:i)
          ipos = ipos + 1
          endif
  101   continue

  901 continue

c The first character of the denonimator cannot be a sign
      if ((cden(1:1).eq.'+').or.(cden(1:1).eq.'-')) then
        errcod = Etype
        chrflt = errcod
        emess = 'Bad floating point number (chrflt)'
        return
        endif

c Compute the numerator
c The length of the numerator to the first blank:
      lnum = chrlen(cnum)
      fnum = 0.D0
      do 104 i = 1, lnum
        chrflt = c1int(cnum(i:i),digit)
        if (chrflt.ne.EOK) return
        fnum = fnum * 10.D0 + digit
 104    continue

c Compute the denominator
c The length of the denonimator to the first blank:
      lden = chrlen(cden)
      fden = 0.D0
      do 105 i = lden, 1, -1
        chrflt = c1int(cden(i:i),digit)
        if (chrflt.ne.EOK) return
        fden = (fden + digit)*0.1D0
 105    continue

c If we have neither a num or den then it is an error
      if ((lden.eq.0).and.(lnum.eq.0)) then
        errcod = Etype
        chrflt = errcod
        emess = 'Bad floating point number (chrflt)'
        return
        endif

c Get the length of the denominator

  902 continue

      result = sign * (fnum + fden)

c     write(6,*) 'chrflt: char =',char
c     write(6,*) ' cnum,cden = ',cnum,cden
c     write(6,*) ' lnum,lden = ', lnum, lden
c     write(6,*) ' sign, fnum, fden =',sign,fnum,fden
c     write(6,*) ' result = ',result

      errcod = EOK
      chrflt = errcod
      emess = ' '
      return
      end
c---------------------------------------------------------------------------
c  Take a character and returns the integer it represents in result.
c The return value is an error code.
      integer function chrint(char,result)
#ifdef UNDEFINED
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) char
      integer result

      integer lchar, sign, starti, nchar, ipower, itmp, i, itmp2
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer c1int

c Compute the length of char
      nchar = 0
      lchar = len(char)
  801 continue
          nchar = nchar + 1
          if (nchar.gt.lchar) goto 901
          if (char(nchar:nchar).ne.' ') goto 801
c Once we get to the first blank, everything else must be blank
          if (char(nchar:lchar).ne.' ') then
            errcod = Esyntx
            chrint = errcod
            emess = 'Blank in the middle of the data (chrint)'
            return
            endif
  901 continue
      nchar = nchar - 1

      ipower = nchar
      itmp = 0
      sign = +1
      starti = 1
      if (char(1:1) .eq. '-') then
         sign = -1
         ipower = ipower - 1
         starti = 2
      else if (char(1:1) .eq. '+') then
         ipower = ipower - 1
         starti = 2
         endif
      if (lchar - starti + 1 .le. 0) then
        errcod = Esyntx
        chrint = errcod
        emess = 'Invalid integer (chrint)'
        return
        endif
      do 101 i = starti, nchar
          ipower = ipower - 1
          chrint = c1int(char(i:i),itmp2)
          if (chrint .ne. EOK) return
          itmp = itmp2*10**ipower + itmp
  101     continue
      result = itmp * sign

      errcod = EOK
      chrint = errcod
      emess = ' '

      return
      end
c---------------------------------------------------------------------------
c Given a character*1 return the integer value of the digit in result.
c The return value is an error code.
      integer function c1int(char,result)
#ifdef UNDEFINED
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*1 char
      integer result

      integer i
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      character*1 digit(0:9)
      data digit/'0','1','2','3','4','5','6','7','8','9'/

      do 101 i = 0, 9
          if (char .eq. digit(i)) then
              errcod = EOK
              c1int = errcod
              emess = ' '
              result=i
              return
              endif
  101     continue

      errcod = Etype
      c1int = Etype
      emess = 'A character other than ''0''-''9'' was found (c1int)'

      end
c---------------------------------------------------------------------------
c  This returns the number of characters preceeding and including the last
c  nonblank character in char.
      integer function chrend(char)
#ifdef UNDEFINED
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
      character*(*) char

      integer i

      do 104 i = len(char), 1, -1
        if (char(i:i) .ne. ' ') then
          chrend = i
          return
          endif
  104   continue

      chrend = 0
      return
      end
c---------------------------------------------------------------------------
c  This returns the number of characters preceeding the first blank in char.
      integer function chrlen(char)
#ifdef UNDEFINED
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
      character*(*) char

      integer i

      do 104 i = 1, len(char)
        if (char(i:i) .eq. ' ') then
          chrlen = i-1
          return
          endif
  104   continue

      chrlen = len(char)
      return
      end
c-----------------------------------------------------------------------------
c Converts a string to lower case
      subroutine tolowr(string)
#ifdef UNDEFINED
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
      character*(*) string

      integer i
      character*(1) uctolc

      do 101 i = 1, len(string)
        string(i:i) = uctolc(string(i:i))
  101   continue

      return
      end
c-----------------------------------------------------------------------------
c Converts a string to upper case
      subroutine touppr(string)
#ifdef UNDEFINED
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
      character*(*) string

      character*(1) lctouc
      integer i

      do 101 i = 1, len(string)
        string(i:i) = lctouc(string(i:i))
  101   continue

      return
      end
c--------------------------------------------------------------
c If given a lower case letter this will return the uppercase equivalent
c Otherwise the given character is returned.
      character*1 function lctouc(char)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      character*1 char
      integer i
      character*26 uplet, lolet
      data lolet/'abcdefghijklmnopqrstuvwxyz'/
      data uplet/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      do 101 i = 1, 26
        if (char .eq. lolet(i:i)) then
          lctouc = uplet(i:i)
          return
          endif
 101    continue
      lctouc = char
      return
      end
c--------------------------------------------------------------
c If given a upper case letter this will return the lowercase equivalent
c Otherwise the given character is returned.
      character*1 function uctolc(char)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      character*1 char
      integer i
      character*26 uplet, lolet
      data lolet/'abcdefghijklmnopqrstuvwxyz'/
      data uplet/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      do 101 i = 1, 26
        if (char .eq. uplet(i:i)) then
          uctolc = lolet(i:i)
          return
          endif
 101    continue
      uctolc = char
      return
      end
