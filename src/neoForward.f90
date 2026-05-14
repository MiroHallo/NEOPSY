! Forward problem (Geopsy 3.X)
!
! Author: Miroslav Hallo, ETH Zurich, Swiss Seismological Service
! Version Fjord (3/2026)
! Hallo et al. (2021, https://doi.org/10.1093/gji/ggab116)
! 
! Copyright (C) 2019-2021 ETH Zurich
! This program is published under the GNU General Public License (GPL).
! This program is free software: you can modify it and/or redistribute it
! or any derivative version under the terms of the GNU General Public
! License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
! This code is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
! and don't remove their names from the code.
! You should have received copy of the GNU General Public License along
! with this program. If not, see <http://www.gnu.org/licenses/>.
!---------------------------------------------------------------------

MODULE ForwardBox
!---------------------------------------------------------------------
!  Module containing variables of the forward problem
!---------------------------------------------------------------------
    implicit none
    real(8),allocatable,dimension(:):: im_thick,im_vs,im_vp,im_rho
!---------------------------------------------------------------------
END MODULE



SUBROUTINE Forward(im,misfit,logPDF,error)
!---------------------------------------------------------------------
!  Calculate forward problem, evaluate misfit and -log(PDF)
!  In:  im - model identifer (integer)
!  Out: misfit - Misfit of data and synthetics
!       logPDF - Negative logrithm of the posterior PDF
!       error - error iostat (=0 if all OK)
!---------------------------------------------------------------------
    use InputBox
    use DataBox
    use GlobalBox
    use LibBox
    use ForwardBox
    use, NON_INTRINSIC :: IO
    implicit none
    integer,intent(in):: im
    real(8),intent(out):: logPDF
    real(8),intent(out):: misfit
    integer,intent(out):: error
    integer(4):: RorL,mode,group,nL,dL
    integer:: m,ci,i,Ind,ios
    real(8):: step,fmin,fmax
    character(len=80):: dir
    character(len=4):: cproc
    character(len=255):: path,file_NP,file_GP,file_DC,file_ELL
    character(len=10):: fminTXT,fmaxTXT,fnTXT
    character(len=9):: toTXT,toTXT1,toTXT2
    character(len=512):: comm
    
    !----------------------------------
    ! Interpret model parameters
    call InterpMod(im,error)
    if(error.ne.0)then
      return
    endif
    
    !----------------------------------
    ! Constants
    nL = m_Nlay(im)
    if(ExpertSwitch(2).eq.1)then
      toTXT='-5'
      toTXT1='-1'
      toTXT2='-2'
    else
      toTXT='-s 9 5'
      toTXT1='-s 9 0.1'
      toTXT2='-s 9 0.5'
    endif
    if(ExpertSwitch(3).eq.1)then
      group = 1
    else
      group = 0
    endif
    
    !----------------------------------
    ! Filenames
    dir = dirG    ! work directory from InputBox
    write(cproc,'(I0)') rank
    path=trim(dir)//'node'//trim(cproc)
    file_NP = trim(path)//'.npin'
    file_DC = trim(path)//'.npdc'
    file_GP = trim(path)//'.model'
    file_ELL = trim(path)//'.ell'
    
    !----------------------------------
    ! Compute Rayleigh and Love wave dispersion curves
    ci = 1
    do m=1,8
      if(f_n(m).gt.0)then
        if(m.le.4)then
          RorL = 0
          mode = m
        else
          RorL = 1
          mode = m-4
        endif
        dL = f_n(m)
        
        ! Write input binary file
        open(unit=31,form='unformatted',access='stream',file=file_NP,action='write',status='replace',iostat=error)
        if(error.ne.0)then
          comm = 'rm '//trim(file_NP)//' > /dev/null 2>&1'
          call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
          return
        endif
        write(31) RorL,mode,group ! 3X integer(4)
        write(31) nL,dL ! 2X integer(4)
        write(31) im_thick(1:nL) ! real(8)
        write(31) im_vs(1:nL) ! real(8)
        write(31) im_vp(1:nL) ! real(8)
        write(31) im_rho(1:nL) ! real(8)
        write(31) f_vec(ci:ci+f_n(m)-1) ! real(8)
        flush(31)
        close(31,iostat=error)
        if(error.ne.0)then
          comm = 'rm '//trim(file_NP)//' > /dev/null 2>&1'
          call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
          return
        endif
        
        ! Delete content of the DC file
        open(unit=32,form='unformatted',access='stream',file=file_DC,action='write',status='replace',iostat=error)
        if(error.ne.0)then
          comm = 'rm '//trim(file_DC)//' > /dev/null 2>&1'
          call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
          return
        endif
        write(32) error
        flush(32)
        close(32,iostat=error)
        if(error.ne.0)then
          comm = 'rm '//trim(file_DC)//' > /dev/null 2>&1'
          call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
          return
        endif
        
        ! Run npdc
        comm = 'timeout '//trim(toTXT)//' ./npdc '//trim(path)//' 2> /dev/null'
        call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
        
        ! Load results
        open(unit=32,form='unformatted',access='stream',file=file_DC,action='read',iostat=error)
        if(error.ne.0)then
          comm = 'rm '//trim(file_DC)//' > /dev/null 2>&1'
          call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
          return
        endif
        read(32,iostat=error) D_tmp(1:f_n(m))
        if(error.ne.0) return
        close(32,iostat=error)
        if(error.ne.0)then
          comm = 'rm '//trim(file_DC)//' > /dev/null 2>&1'
          call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
          return
        endif
        if(group.eq.0)then
          do i=1,f_n(m)
            if(D_tmp(i).le.0.d0)then
              error = 2
              return
            endif
          enddo
        else
          do i=2,f_n(m)-1
            if(D_tmp(i).le.0.d0)then
              error = 2
              return
            endif
          enddo
        endif
        ! Vector of synthetic data
        D_vec(ci:ci+f_n(m)-1) = D_tmp(1:f_n(m))
        ci = ci+f_n(m)
      endif
    enddo
    
    !----------------------------------
    ! Compute Rayleigh wave ellipticity curve (abs or angle)
    if(f_n(9).gt.0 .or. f_n(10).gt.0)then
      if(f_n(9).gt.0)then
        m = 9
      else
        m = 10
      endif
      ! Write the velocity model to ASCII file
      call INV2GEOPSY(file_GP,im_thick(1:nL),im_vp(1:nL),im_vs(1:nL),im_rho(1:nL),error)
      if(error.ne.0)then
        comm = 'rm '//trim(file_GP)//' > /dev/null 2>&1'
        call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
        return
      endif
      ! Data oversampling
      if((log10(f_lim(2,m))-log10(f_lim(1,m))).lt.0.2d0)then
        dL = max(50,f_n(m))
      elseif((log10(f_lim(2,m))-log10(f_lim(1,m))).lt.2.d0)then
        dL = max(100,f_n(m))
      else
        dL = max(200,f_n(m))
      endif
      step = (f_lim(2,m)/f_lim(1,m))**(1.d0/(dble(dL-1)))
      fmin = f_lim(1,m)*(step**(-1))
      fmax = f_lim(1,m)*(step**(dL))
      
      ! Delete content of the ELL file
      open(unit=33,form='formatted',file=file_ELL,status='replace',iostat=error)
      if(error.ne.0)then
        comm = 'rm '//trim(file_ELL)//' > /dev/null 2>&1'
        call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
        return
      endif
      write(33,'(I3)') error
      flush(33)
      close(33,iostat=error)
      if(error.ne.0)then
        comm = 'rm '//trim(file_ELL)//' > /dev/null 2>&1'
        call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
        return
      endif
      
      ! Run gpell
      write(fminTXT,'(F7.3)') fmin
      write(fmaxTXT,'(F7.3)') fmax
      write(fnTXT,'(I0)') dL
      ! First try
      comm = 'timeout '//trim(toTXT1)//' gpell -min '//trim(fminTXT)//' -max '//trim(fmaxTXT)//' -n '//&
             trim(fnTXT)//' < '//trim(file_GP)//' > '//trim(file_ELL)//' 2> /dev/null'
      call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
      call GEOPSY2INV(file_ELL, f_tmp(1:dL), D_tmp(1:dL), ios)
      
      if(ios.ne.0)then
        ! Second try
        comm = 'timeout '//trim(toTXT2)//' gpell -min '//trim(fminTXT)//' -max '//trim(fmaxTXT)//' -n '//&
             trim(fnTXT)//' < '//trim(file_GP)//' > '//trim(file_ELL)//' 2> /dev/null'
        call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
        call GEOPSY2INV(file_ELL, f_tmp(1:dL), D_tmp(1:dL), ios)
        
        if(ios.ne.0)then
          ! The last try
          comm = 'timeout '//trim(toTXT)//' gpell -min '//trim(fminTXT)//' -max '//trim(fmaxTXT)//' -n '//&
             trim(fnTXT)//' < '//trim(file_GP)//' > '//trim(file_ELL)//' 2> /dev/null'
          call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
          ! Load results
          call GEOPSY2INV(file_ELL, f_tmp(1:dL), D_tmp(1:dL), error)
          
          if(error.ne.0) return
          
        endif
      endif
      
      ! Vector of synthetic data
      if(m.eq.9)then
        D_tmp(1:dL) = dlog10(abs(D_tmp(1:dL)))
      else
        D_tmp(1:dL) = atan(D_tmp(1:dL))*57.29577951308d0
      endif
      
      ! Data interpolate
      do i=1,f_n(m)
        Ind = binarysearch(dL,f_tmp(1:dL),f_vec(ci+i-1))
        if(Ind.lt.1)then
          D_vec(ci+i-1) = D_tmp(1)
        elseif(Ind.ge.dL)then  
          D_vec(ci+i-1) = D_tmp(dL)
        elseif(m.eq.10 .and. abs(D_tmp(Ind+1)-D_tmp(Ind)).gt.45.d0)then
          D_vec(ci+i-1) = D_tmp(Ind)
        else
          D_vec(ci+i-1) = D_tmp(Ind) + (log10(f_vec(ci+i-1))-log10(f_tmp(Ind))) *&
                       (D_tmp(Ind+1)-D_tmp(Ind)) / (log10(f_tmp(Ind+1))-log10(f_tmp(Ind)))
        endif
      enddo
      ci = ci+f_n(m)
    endif
    
    !----------------------------------
    ! Prepare standardized synthetic data
    D_vec = D_vec*W_vec
    
    !----------------------------------
    ! Evaluate misfit function
    misfit=dot_product(D_vec-UW_vec,D_vec-UW_vec)
    
    !----------------------------------
    ! NaN test
    if(isnan(misfit))then
      error = 3
      return
    endif
    if(misfit.lt.0.d0)then
      error = 3
      return
    endif
    if(misfit.gt.1.d300)then
      error = 3
      return
    endif
    
    !----------------------------------
    ! Calculate minus logPDF
    ! Note: logPDF = -log( exp(-0.5d0*misfit) )
    logPDF = (0.5d0*misfit)
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE InterpMod(im,error)
!---------------------------------------------------------------------
!  Interprets the model with identifer im into a velocity profile
!  It interprets into im_thick, im_vs, im_vp, im_rho arrays from ForwardBox
!  error 0 - NO_ERROR
!---------------------------------------------------------------------
    use GlobalBox
    use ForwardBox
    implicit none
    integer,intent(in):: im
    integer,intent(out):: error
    integer:: i
    real(8):: boundary0,boundary
    
    error=0
    !----------------------------------
    ! Make it zero, for sure
    im_thick=0.d0
    im_vs=0.d0
    im_vp=0.d0
    im_rho=0.d0
    
    !----------------------------------
    ! Interprets all layeres
    boundary0 = 0.d0
    do i=1,m_Nlay(im)-1
      ! Find thickness
      boundary = m_dspr(1,i,im) + (m_dspr(1,i+1,im)-m_dspr(1,i,im))/2.d0
      im_thick(i) = boundary-boundary0
      boundary0 = boundary
      ! Assign the other parameters
      im_vs(i) = m_dspr(2,i,im)
      im_vp(i) = m_dspr(3,i,im)
      im_rho(i) = abs(m_dspr(4,i,im))
    enddo
    
    !----------------------------------
    ! Interprets bottom halfspace
    im_thick(m_Nlay(im)) = 0.d0
    im_vs(m_Nlay(im)) = m_dspr(2,m_Nlay(im),im)
    im_vp(m_Nlay(im)) = m_dspr(3,m_Nlay(im),im)
    im_rho(m_Nlay(im)) = abs(m_dspr(4,m_Nlay(im),im))
    
!---------------------------------------------------------------------
END SUBROUTINE


