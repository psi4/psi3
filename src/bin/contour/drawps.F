c The values in the drawps common block
c real*8  conval = contour value
c integer npoint = the number of points drawn so far 
c integer copen = 1 if a contour is started
c integer usecol = 1 if color is to be used

c-------------------------------------------------------------------------
c Initialize the postscript file.
      subroutine initps(x0,y0,x1,y1,label,color,cmin,cmax,cstep,
     &            orbtyp,orbnum,rotate)
      real*8 x0,y0,x1,y1,cmin,cmax,cstep
      character*(*) label,orbtyp
      integer color,orbnum,rotate
c The width and height of the drawing
      real*8 xlen,ylen
c The center of the drawing
      real*8 xcen,ycen
c Scale factors
      real*8 xscale, yscale, scale
c The page length, page height, text area size, hash mark size, and page border
      real*8 px,py,pt,pm,pb
      real*8 tmp
      integer start, end, mark
      integer chrend
      character*100 msg,chrnum

      real*8 conval
      integer npoint,copen,usecol
      common/drawps/ conval,npoint,copen,usecol

      copen = 0
      usecol = color

      px = 612
      py = 792
      pt = 30
      pm = 26
      pb = 72

      xlen = x1-x0
      ylen = y1-y0

      if (xlen.lt.0.) then
        write(6,*) 'ERROR: xlen < 0'
        call qabort
        endif

      if (ylen.lt.0.) then
        write(6,*) 'ERROR: ylen < 0'
        call qabort
        endif

      xcen = x0 + xlen/2
      ycen = y0 + ylen/2

c Open the postscript output file.
      call fioopn(7,'contour.ps','USER','WRITE')
      write(7,1) '%! Contour plotting program for molecular orbitals'
      write(7,1) '%  Program by Yaoming Xie'
      write(7,1) '%  Postscript output by Curtis Janssen'
    1 format(10a)

c Define some functions to reduce the length of the output
      write(7,1) '/L {lineto} def'
      write(7,1) '/M {moveto} def'
      write(7,1) '/R {rlineto} def'

      write(7,2) '% Box ranges:',x0,x1,y0,y1
    2 format(a,2(' (',f10.4,', ',f10.4,')'))

c Move the origin to the center of the page
      write(7,1) '% Move the origin to the center of the page'
      write(7,3) px/2, py/2, 'translate'
    3 format(f10.4,1x,f10.4,1x,a)

c See if we should rotate the page
      if (    (rotate.eq.1)
     &    .or.((rotate.eq.-1).and.(xlen.gt.ylen))) then
        write(7,1) '% Rotate the page by 90 degrees'
        write(7,1) '90 rotate'
        tmp = px
        px = py
        py = tmp
        endif

c Adjust the page to make room for the text
      write(7,1) '% Adjust the origin to make room for text'
      write(7,3) 0.,pt,'translate'
      py = py - pt

c Adjust the page to make room for the hash marks
      write(7,1) '% Adjust the page to make room for the hash marks'
      write(7,3) pm,pm,'translate'
      py = py - pm
      px = px - pm

c Get a font
      write(7,1) '% Make a label'
      write(7,4) '/Helvetica findfont',12.0,'scalefont setfont'
    4 format(a,1x,f10.4,1x,a)

c Display the label
      write(7,1) '(',label(1:chrend(label)),')'
      write(7,1) 'dup stringwidth pop 2 div neg'
      write(7,5) - (py-2*pb)/2 - pm - 14, 'M show'
    5 format(f10.4,1x,a)

c Display other information about the plot.
      msg = ' '
      call scatb(msg,orbtyp)
      call scatb(msg,'number')
      call intchr(orbnum,chrnum)
      call scatb(msg,chrnum)
      write(7,1) '(',msg(1:chrend(msg)),')'
      write(7,1) 'dup stringwidth pop 2 div neg'
      write(7,5) - (py-2*pb)/2 - pm - 30, 'M show'

c Determine the scale factor
      write(7,1) '% Scale the page to the units are now Bohr'
      xscale = (px - 2*pb)/xlen
      yscale = (py - 2*pb)/ylen
      if (xscale.gt.yscale) then
        scale = yscale
      else
        scale = xscale
        endif
      write(7,3) scale,scale,'scale'
c Converts points into bohr
      write(7,4) '/point {',scale,'div} def'

c Translate to make the picture fit within the page
      write(7,1) '% Translate to the plot will fit in the drawing area'
      write(7,3) -xcen,-ycen,'translate'

c Draw a border
      write(7,1) '% Draw a border'
      write(7,1) '0.4 point setlinewidth'
      write(7,3) x0,y0,'M'
      write(7,3) x1,y0,'L'
      write(7,3) x1,y1,'L'
      write(7,3) x0,y1,'L'
      write(7,3) x0,y0,'L'
      write(7,1) 'stroke'

c Get a font
      write(7,1) '% Get a font for the hash marks labels'
      write(7,4) '/Helvetica findfont',12.0,'point scalefont setfont'

c Draw some marks around the border

      write(7,1) '% Draw and label the hash marks'
c ........ X marks
      start = dint(x0*10 + 1)
      end = dint(x1*10 - 1)
      do 101 mark = start, end

c Make the mark
        write(7,3) mark/10.0d0,y0,'M'
        if (10*(mark/10).eq.mark) then
          tmp = -6.0
        else if (5*(mark/5).eq.mark) then
          tmp = -4.0
        else
          tmp = -2.0
          endif
        write(7,3) 0.,tmp,'point R'
        write(7,1) 'stroke'

c Show some text
        if (10*(mark/10).eq.mark) then
          call genstr(mark/10)
          write(7,1) 'dup stringwidth pop 2 div neg'
          write(7,5) mark/10.0d0, 'add'
          write(7,5) y0, '-22 point add M show'
          endif
  101   continue

c ........ Y marks
      start = dint(y0*10 + 1)
      end = dint(y1*10 - 1)
      do 102 mark = start, end

c Make the mark
        write(7,3) x0, mark/10.0d0,'M'
        if (10*(mark/10).eq.mark) then
          tmp = -6.0
        else if (5*(mark/5).eq.mark) then
          tmp = -4.0
        else
          tmp = -2.0
          endif
        write(7,3) tmp,0.,'exch point exch R'
        write(7,1) 'stroke'

c Show some text
        if (10*(mark/10).eq.mark) then
          call genstr(mark/10)
          write(7,1) 'dup stringwidth pop -4 point'
          write(7,5) mark/10.0d0, 'add exch neg'
          write(7,5) x0, '-10 point add add exch M show'
          endif
  102   continue

c Set up the clipping path
clj This is acting wierd and I don't have time to fix it:
clj   write(7,1) '% Clip off anything drawn outside the picture'
clj   write(7,1) 'newpath'
clj   write(7,3) x0,y0,'M'
clj   write(7,3) x1,y0,'L'
clj   write(7,3) x1,y1,'L'
clj   write(7,3) x0,y1,'L'
clj   write(7,3) x0,y0,'L'
clj   write(7,1) 'closepath clippath'

c Set the line width for the plot
      write(7,1) '% Set the line width for the plot'
      write(7,1) '0.7 point setlinewidth'

      return
      end

c-------------------------------------------------------------------------
c This draws atoms
      subroutine drawat(geom,chgs,natom,x0,y0,x1,y1,z,plane)
      integer natom
      real*8 geom(3,natom),chgs(natom),x0,y0,x1,y1,z
      character*(*) plane

      real*8 conval
      integer npoint,copen,usecol
      common/drawps/ conval,npoint,copen,usecol

      integer i
      character*(1) cx,cy
      real*8 ax,ay,az

      cx = plane(1:1)
      cy = plane(2:2)

      write(7,1) '% Mark where the atoms appear'
   1  format(a)
      do 101 i = 1, natom
        if (cx.eq.'X') then
          ax = geom(1,i)
          if (cy.eq.'Y') then
            ay = geom(2,i)
            az = geom(3,i)
          else
            ay = geom(3,i)
            az = geom(2,i)
            endif
        else if (cx.eq.'Y') then
          ax = geom(2,i)
          if (cy.eq.'X') then
            ay = geom(1,i)
            az = geom(3,i)
          else
            ay = geom(3,i)
            az = geom(1,i)
            endif
        else
          ax = geom(3,i)
          if (cy.eq.'X') then
            ay = geom(1,i)
            az = geom(2,i)
          else
            ay = geom(2,i)
            az = geom(1,i)
            endif
          endif
        if (     (ax.ge.x0).and.(ax.le.x1)
     &      .and.(ay.ge.y0).and.(ay.le.y1))
     &    call atomou(chgs(i),ax,ay,az,z)
  101   continue

      return
      end

c-------------------------------------------------------------------------
c This will output a single atom, given info about it.
      subroutine atomou(charge,x,y,z,c0)
      real*8 charge,x,y,z,c0

      real*8 conval
      integer npoint,copen,usecol
      common/drawps/ conval,npoint,copen,usecol

c Save the current graphics state
      write(7,1) 'gsave'

      write(7,1) '0 setgray'
      write(7,1) '0.4 point setlinewidth'
    1 format(10a)

c Make atoms appearing in the plane solid and those out of the plane
c dashed.
C  YX  12/19/91.  use big dot instead of cross.
c     if ( dabs(c0-z).gt.1.d-6 ) then
c       write(7,1) '[ 2 point 2 point ] 0 setdash'
c       endif

c     write(7,2) x,y,'M  6 point  6 point rlineto stroke'
c     write(7,2) x,y,'M -6 point  6 point rlineto stroke'
c     write(7,2) x,y,'M -6 point -6 point rlineto stroke'
c     write(7,2) x,y,'M  6 point -6 point rlineto stroke'
      write(7,2) x,y,' 2 point 0 360 arc '
      write(7,1) '  gsave'
      write(7,1) '  fill'
      write(7,1) '  grestore'
      write(7,1) 'stroke'
C  YX  12/19/91 end of change
    2 format(f10.4,1x,f10.4,1x,a)

      write(7,1) 'grestore'
      return
      end

c-------------------------------------------------------------------------
c This converts the mark to a string--not a pleasant task in Fortran 77
      subroutine genstr(mark)
      integer mark

   1  format('(',I1,'.0)')
   2  format('(',I2,'.0)')
   3  format('(',I3,'.0)')
   4  format('(',I4,'.0)')
   5  format('(',I5,'.0)')

      if (mark.lt.-1000) then
        write(7,5) mark
      else if (mark.lt.-100) then
        write(7,4) mark
      else if (mark.lt.-10) then
        write(7,3) mark
      else if (mark.lt.0) then
        write(7,2) mark
      else if (mark.lt.10) then
        write(7,1) mark
      else if (mark.lt.100) then
        write(7,2) mark
      else if (mark.lt.1000) then
        write(7,3) mark
      else
        write(7,4) mark
        endif

      return
      end

c-------------------------------------------------------------------------
c Finish up with the postscript output
      subroutine doneps()
      write(7,1) 'showpage'
    1 format(a)
      return
      end

c-------------------------------------------------------------------------
c Start a contour (val = the contour value)
      subroutine startc(val)
      real*8 val
      real*8 conval
      integer npoint,copen,usecol
      common/drawps/ conval,npoint,copen,usecol

      if (copen.eq.1) then
        write(6,*) 'ERROR: copen is 1 in startc'
        call qabort
        endif
      copen = 1
      conval = val
      npoint = 0
      write(7,*)
      write(7,1) val
   1  format('% Starting contour with value ',f10.4)

c Draw negative parts:
      if (val.lt.-1.d-6) then
        if (usecol.ne.0) then
          write(7,2) '1 0 0 setrgbcolor'
   2      format(a)
        else
          write(7,2) '[ .03 .03 ] 0 setdash'
          endif
c Draw nodes:
      else if (val.lt.1.d-6) then
        if (usecol.ne.0) then
          write(7,2) '0 setgray'
        else
CYX  12/19/91    change the length of dashes.
C         write(7,2) '[ .10 .03 ] 0 setdash'
          write(7,2) '[ .30 .06 ] 0 setdash'
          endif
c Draw positive parts:
      else
        if (usecol.ne.0) then
          write(7,2) '0 0 1 setrgbcolor'
        else
          write(7,2) '[ ] 0 setdash'
          endif
        endif
      return
      end

c-------------------------------------------------------------------------
c The next point for the contour
      subroutine point(x,y)
      real*8 x,y
      real*8 conval
      integer npoint,copen,usecol
      common/drawps/ conval,npoint,copen,usecol

      if (copen.ne.1) then
        write(6,*) 'ERROR: point was called but no contour is open'
        call qabort
        endif

      if (npoint.eq.0) then
        firstl = 0
        write(7,1) x,y,'M'
   1    format(f10.4,1x,f10.4,1x,a)
      else
        firstl = 0
        write(7,1) x,y,'L'
        endif

      if ((npoint.gt.0).and.(600*(npoint/600).eq.npoint)) then
        write(7,2) 'currentpoint'
        write(7,2) 'stroke'
        write(7,2) 'M'
   2    format(a)
        endif

      npoint = npoint + 1

      return
      end

c-------------------------------------------------------------------------
c For the whole set of the points of the contour line (either closed or open).
cyx 08/26/92
      subroutine points(xay,ipnt)
      real*8 xay(2,ipnt)
      real*8 conval
      integer npoint,copen,usecol
      common/drawps/ conval,npoint,copen,usecol

      if (copen.ne.1) then
        write(6,*) 'ERROR: point was called but no contour is open'
        call qabort
        endif

        firstl = 0
        write(7,1) xay(1,1),xay(2,1),'M'
   1    format(f10.4,1x,f10.4,1x,a)
        do 10 i=2,ipnt
        write(7,1) xay(1,i),xay(2,i),'L'
 10     continue

        return
        end
c-------------------------------------------------------------------------
c We reached a boundary
      subroutine bound
      call closec()
      return
      end

c-------------------------------------------------------------------------
c Close a contour
      subroutine closec
      real*8 conval
      integer npoint,copen,usecol
      common/drawps/ conval,npoint,copen,usecol

      if (copen.eq.0) then
        write(6,*) 'ERROR: closec was called but there is no contour'
        call qabort
        endif
      copen = 0
      write(7,1) 'stroke'
   1  format(a)
      return
      end
