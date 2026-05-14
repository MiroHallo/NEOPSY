! Initiation from DATA files - MTI
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

MODULE DataBox
!---------------------------------------------------------------------
!  Module containing data and their weights
!---------------------------------------------------------------------
    implicit none
    real(8),allocatable,dimension(:):: f_vec,U_vec,W_vec,UW_vec
    integer,allocatable,dimension(:):: f_n
    real(8),allocatable,dimension(:,:):: f_lim
    real(8),allocatable,dimension(:):: f_tmp,D_tmp,D_vec !for synthetics
    real(8),allocatable,dimension(:,:):: init_model
    integer:: maxMode,N_used,init_nL
!---------------------------------------------------------------------
END MODULE



SUBROUTINE InitData()
!---------------------------------------------------------------------
!  Read data files and create data and weight vectors
!---------------------------------------------------------------------
    use InputBox
    use DataBox
    implicit none
    integer:: maxN,ios,fcount,m,ci,i
    real(8),allocatable,dimension(:,:):: f_mat,U_mat,W_mat
    real(8),allocatable,dimension(:):: U_max
    character(len=255):: file_tmp
    character(len=80):: dir
    character(len=255):: filename
    real(8):: WlevelS,WlevelE,WlevelA
    
    !----------------------------------
    !  Prepare variables for data
    maxMode = 10  ! The last is for ellipticity
    maxN = 1000   ! Maximal number of discrete frequencies (of data samples)
    dir = dirG           ! work directory
    WlevelS = 1.d-4      ! water level of stddev of slowness [s/m]
    WlevelE = 0.14d0     ! water level of stddev of log10(ellipticity)
    WlevelA = 1.d0       ! water level of stddev of ellipticity angle
    
    !----------------------------------
    !  Allocate global variables
    allocate(f_n(maxMode))
    allocate(f_lim(2,maxMode))
    allocate(f_tmp(maxN))
    allocate(D_tmp(maxN))
    f_n = 0
    f_lim = 0.d0
    f_tmp = 0.d0
    D_tmp = 0.d0
    
    !----------------------------------
    !  Allocate local variables
    allocate(f_mat(maxN,maxMode))
    allocate(U_mat(maxN,maxMode))
    allocate(W_mat(maxN,maxMode))
    allocate(U_max(maxMode))
    f_mat = 0.d0
    U_mat = 0.d0
    W_mat = 0.d0
    U_max = 0.d0
    
    !----------------------------------
    !  Read initial model (if required)
    if(ExpertSwitch(4).eq.1)then
      file_tmp = dirINI
      ! Read the data
      open(169,form='formatted',file=file_tmp,action='read',status='old') ! Open the file
      read(169,*)
      read(169,*) init_nL
      read(169,*)
      allocate(init_model(4,init_nL))
      init_model = 0.d0
      do i=1,init_nL-1
        read(169,*,iostat=ios) init_model(1:4,i)
        if(ios.ne.0)then
          write(*,*) 'ERROR 214: File with the initial model is corrupted!'
          stop
        endif
      enddo
      read(169,*)
      read(169,*,iostat=ios) init_model(1:4,init_nL)
      if(ios.ne.0)then
        write(*,*) 'ERROR 214: File with the initial model is corrupted!'
        stop
      endif
      close(169)
    else
      allocate(init_model(4,1))
      init_model = 0.d0
      init_nL = 0
    endif
    
    !----------------------------------
    !  Read all data into data matrices
    do m=1,maxMode
    
      ! Select data type
      if(m.eq.1)then
        if(DCRm(1).le.0) cycle
        file_tmp = dirDCR0
      elseif(m.eq.2)then
        if(DCRm(2).le.0) cycle
        file_tmp = dirDCR1
      elseif(m.eq.3)then
        if(DCRm(3).le.0) cycle
        file_tmp = dirDCR2
      elseif(m.eq.4)then
        if(DCRm(4).le.0) cycle
        file_tmp = dirDCR3
      elseif(m.eq.5)then
        if(DCLm(1).le.0) cycle
        file_tmp = dirDCL0
      elseif(m.eq.6)then
        if(DCLm(2).le.0) cycle
        file_tmp = dirDCL1
      elseif(m.eq.7)then
        if(DCLm(3).le.0) cycle
        file_tmp = dirDCL2
      elseif(m.eq.8)then
        if(DCLm(4).le.0) cycle
        file_tmp = dirDCL3
      elseif(m.eq.9)then
        if(ELLm(1).le.0) cycle
        file_tmp = dirELL
      elseif(m.eq.10)then
        if(ELLm(2).le.0) cycle
        file_tmp = dirELA
      else
        write(*,*) 'ERROR: Unknown observed data!'
        stop
      endif
      
      ! Read the data
      open(170,form='formatted',file=file_tmp,action='read',status='old') ! Open the file
      read(170,*)
      fcount = 1
      do while(.true.) ! Loop over lines within the file
        read(170,*,iostat=ios) f_mat(fcount,m),U_mat(fcount,m),W_mat(fcount,m)
        if(ios.ne.0) exit
        fcount = fcount + 1
      enddo
      close(170)
      f_n(m) = fcount - 1
      if(fcount.gt.1)then
        f_lim(1,m) = f_mat(1,m)
        f_lim(2,m) = f_mat(f_n(m),m)
        U_max(m) = maxval(abs(U_mat(1:f_n(m),m)))
      endif
      
    enddo
    
    !----------------------------------
    !  Allocate global vectors
    allocate(f_vec(sum(f_n)))
    allocate(U_vec(sum(f_n)))
    allocate(W_vec(sum(f_n)))
    allocate(UW_vec(sum(f_n)))
    allocate(D_vec(sum(f_n)))
    f_vec = 0.d0
    U_vec = 0.d0
    W_vec = 1.d0
    UW_vec = 0.d0
    D_vec = 0.d0
    
    !----------------------------------
    !  Re-structulize into data vectors
    ci = 1
    N_used = 0
    ! Dispersion curves (R0,R1,R2,R3,L0,L1,L2,L3)
    do m=1,8
      if(f_n(m).eq.0) cycle
      f_vec(ci:ci+f_n(m)-1) = f_mat(1:f_n(m),m)
      U_vec(ci:ci+f_n(m)-1) = U_mat(1:f_n(m),m)
      !  std to weight
      do i=1,f_n(m)
        if(W_mat(i,m).gt.0.d0)then
          W_vec(ci+i-1) = 1.d0 / W_mat(i,m)
          N_used = N_used+1
        elseif(W_mat(i,m).eq.0.d0)then
          W_vec(ci+i-1) = 1.d0 / WlevelS
          N_used = N_used+1
        else
          W_vec(ci+i-1) = 0.d0
        endif
      enddo
      !  group velocity
      if(ExpertSwitch(3).eq.1)then
          W_vec(ci) = 0.d0
          W_vec(ci+f_n(m)-1) = 0.d0
      endif
      ci = ci+f_n(m)
    enddo
    ! Ellipticity (ELL)
    m = 9
    if(f_n(m).gt.0)then
      f_vec(ci:ci+f_n(m)-1) = f_mat(1:f_n(m),m)
      U_vec(ci:ci+f_n(m)-1) = dlog10(U_mat(1:f_n(m),m))
      do i=1,f_n(m)
        if(W_mat(i,m).gt.1.d0)then
          W_vec(ci+i-1) = 1.d0 / dlog10(W_mat(i,m))
          N_used = N_used+1
        elseif(W_mat(i,m).eq.1.d0)then
          W_vec(ci+i-1) = 1.d0 / WlevelE
          N_used = N_used+1
        else
          W_vec(ci+i-1) = 0.d0
        endif
      enddo
      ci = ci+f_n(m)
    endif
    ! Ellipticity angle (ELA)
    m = 10
    if(f_n(m).gt.0)then
      f_vec(ci:ci+f_n(m)-1) = f_mat(1:f_n(m),m)
      U_vec(ci:ci+f_n(m)-1) = U_mat(1:f_n(m),m)
      !  std to weight
      do i=1,f_n(m)
        if(W_mat(i,m).gt.0.d0)then
          W_vec(ci+i-1) = 1.d0 / W_mat(i,m)
          N_used = N_used+1
        elseif(W_mat(i,m).eq.0.d0)then
          W_vec(ci+i-1) = 1.d0 / WlevelA
          N_used = N_used+1
        else
          W_vec(ci+i-1) = 0.d0
        endif
      enddo
      ci = ci+f_n(m)
    endif
    
    !----------------------------------
    ! Prepare standardized data UW
    UW_vec = U_vec*W_vec
    
    !----------------------------------
    ! Write init data
    filename = trim(dir)//'in_data.txt'
    open(190,file=filename)
    write(190,'(A)') '# Number of frequency samples for each mode (R0,R1,R2,R3,L0,L1,L2,L3,ELL,ELA)'
    write(190,'(1000I5)') f_n
    do m=1,maxMode
      if(f_n(m).eq.0) cycle
      ci = sum(f_n(1:m-1))
      if(m.le.4)then
        write(190,'(A,I3)') '# f[Hz] U[s/m] W[m/s] UW[standardized] - Dispersion R',m-1
      elseif(m.le.8)then
        write(190,'(A,I3)') '# f[Hz] U[s/m] W[m/s] UW[standardized] - Dispersion L',m-5
      elseif(m.eq.9)then
        write(190,'(A)') '# f[Hz] log10(U(H/V)) W[unitless] UW[standardized] - Ellipticity'
      elseif(m.eq.10)then
        write(190,'(A)') '# f[Hz] alfa[deg] W[1/deg] UW[standardized] - Ellipticity angle'
      else
        write(190,'(A)') '# Unknown data'
      endif
      do i=1,f_n(m)
        write(190,'(4E17.7)') f_vec(ci+i),U_vec(ci+i),W_vec(ci+i),UW_vec(ci+i)
      enddo
    enddo
    close(190)
    
    !----------------------------------
    !write(*,*) '-> Data OK'
    
    deallocate(f_mat,U_mat,W_mat,U_max)
    
!---------------------------------------------------------------------
END SUBROUTINE


