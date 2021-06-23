C======================================================================|
      subroutine NEOCL2
C----------------------------------------------------------------------|
      integer	j
      save	j
      data	j/0/
C----------------------------------------------------------------------|
      if (j .eq. 0)	then
         write(*,*)" >>> Warning >>> the subroutine NEOCL2 is obsolete."
         write(*,*)
     >   "     Call will be ignored. Please use NEOCL(100) instead"
         j = 1
         pause
      endif
      end
C======================================================================|

