      integer EOK, Etype, Ekey, Elngth, Esyntx, Einflp, Efind, EEOF
      integer Etrunc, Eseg, Eparen, Eevect, Erdrct, EEOD, Earg, Ebvect
      integer Eopen, Eclose, Elimit, Emrpar
      parameter (
c                    OK, no error.
     & EOK = 0,
c                    Bad type in a formatted read.
     & Etype = 1,
c                    Could not find the keyword.
     & Ekey  = 2,
c                    Bad array length.
     & Elngth  = 3,
c                    Syntax error.
     & Esyntx  = 4,
c                    Infinite loop
     & Einflp  = 5,
c                    Efind (not used, 6 is available for replacement)
     & Efind = 6,
c                    End of File
     & EEOF = 7,
c                    Character value has been truncated
     & Etrunc = 8,
c                    Too many segments in a keyword
     & Eseg = 9,
c                    Too many levels of parenthesis
     & Eparen = 10,
c                    End of vector (not really an error)
     & Eevect = 11,
c                    Redirection depth is too deep.
     & Erdrct = 12,
c                    End of Data is like EOF but expected.
     & EEOD = 13,
c                    Bad value in argument list
     & Earg = 14,
c                    Could not open a file.
     & Eopen = 15,
c                    Beginning of vector (not really an error)
     & Ebvect = 16,
c                    Some error occurred while closing a file.
     & Eclose = 17,
c                    An internal limit has been exceeded.
     & Elimit = 18,
c                    A right parenthesis is missing
     & Emrpar = 19 )
