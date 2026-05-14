! Module with user functions (Parametric Slip Inversion)
!
! Author: Miroslav Hallo (10/2017)
! Charles University, Prague, Faculty of Mathematics and Physics
!
! This code is published under the GNU General Public License (GPL)
! for non-commercial purposes. To any licensee is given permission to 
! modify the work, as well as to copy and redistribute. Still
! we would like to kindly ask you to acknowledge the authors and don't
! remove their names from the code. This code is distributed in the
! hope that it will be useful, but WITHOUT ANY WARRANTY.
!---------------------------------------------------------------------

MODULE LibBox
!---------------------------------------------------------------------
!  Module containing global functions
!---------------------------------------------------------------------
CONTAINS
    
    
    
      FUNCTION GASDEV(idum)
!-------------------------------------------------------------------
! Numerical Recipes random number generator for 
! a Gaussian distribution (FORTRAN 77)
!-------------------------------------------------------------------
      integer          idum
      real GASDEV
      real v1,v2,r,fac
      real ran3
      if (idum.lt.0) iset=0
 10   v1=2*ran3(idum)-1
      v2=2*ran3(idum)-1
      r=v1**2+v2**2
      if(r.ge.1.or.r.eq.0) GOTO 10
      fac=sqrt(-2*log(r)/r)
      GASDEV=v2*fac
      RETURN
!-------------------------------------------------------------------
      END



FUNCTION binarysearch(length,array,value)
!---------------------------------------------------------------------
! Given an array and a value, returns the index of the element that
! is closest to, but less than, the given value.
! Uses a binary search algorithm.
! Used by interp2D function
!---------------------------------------------------------------------
    implicit none
    real(8),parameter:: dpr=1.d-9
    integer,intent(in):: length
    real(8),dimension(length),intent(in):: array
    real(8),intent(in):: value
    integer:: binarysearch
    integer:: left, middle, right
    !----------------------------------
    !  Search for index
    left = 1
    right = length
    do
      if (left > right) then
        exit
      endif
      middle = nint((left+right) / 2.0)
      if ( abs(array(middle) - value) <= dpr) then
        binarySearch = middle
        RETURN
      elseif (array(middle) > value) then
        right = middle - 1
      else
        left = middle + 1
      endif
    enddo
    binarysearch = right
    RETURN
!---------------------------------------------------------------------
END FUNCTION



!---------------------------------------------------------------------
END MODULE


