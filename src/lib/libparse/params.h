c Length of the error message buffer
c (this is currently hardwired to 80 in some places)
      integer lemess
      parameter (lemess = 80)
c Maximum length of a token
      integer lentok
      parameter (lentok = 80)
c Maximum number of segments in a keyword
      integer maxseg
      parameter (maxseg = 20)
c Maximum number of current working keywords
      integer maxcwk
      parameter (maxcwk = 10)
c Maximum length of a filename
      integer maxnam
      parameter (maxnam = 100)
c Maximum number of units (if this is less than 100 then problems may occur,
c if this is more than one hundred then some portability problems may
c arise)
      integer maxunt
      parameter (maxunt = 100)
