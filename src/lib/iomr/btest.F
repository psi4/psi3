      logical function btest(flag,bit)
      integer flag
      integer bit
      integer flag2
      flag2 = flag
      do 101 i = 1, bit
        flag2 = flag2/2
 101  continue
      btest = .not. (((flag2/2)*flag2) .eq. flag2)
      return
      end
