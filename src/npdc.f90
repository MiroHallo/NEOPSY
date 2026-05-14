! Computes phase or group velocity dispersion curve for either 
! Rayleigh or Love waves given a layered 1D earth model.
!
! Author: Miroslav Hallo (4/2021)
! ETH Zurich, Swiss Seismological Service
!
! This code is published under the GNU General Public License (GPL)
! for non-commercial purposes. To any licensee is given permission to 
! modify the work, as well as to copy and redistribute. Still
! we would like to kindly ask you to acknowledge the authors and don't
! remove their names from the code. This code is distributed in the
! hope that it will be useful, but WITHOUT ANY WARRANTY.
!---------------------------------------------------------------------

PROGRAM  npdc
!---------------------------------------------------------------------
!  Main program
!---------------------------------------------------------------------
    implicit none
    integer(4):: RorL,mode,group,nL,dL,verbose
    real(8),dimension(100):: im_thick,im_vs,im_vp,im_rho
    real(8),dimension(1000):: freq,omega
    real(8),dimension(4000):: dc
    integer:: i,ifile,ios
    character(len=255):: path,filename
    
    !----------------------------------
    ! Read command-line arguments
    if(COMMAND_ARGUMENT_COUNT().NE.1)then
      write(*,*) 'ERROR npdc: Invalid command-line argument'
      stop
    endif
    call GET_COMMAND_ARGUMENT(1,path)
    
    !----------------------------------
    ! Write .npin binary file (for debugging)
    !if(.false.)then
    !! Settings
    !RorL = 1
    !mode = 1
    !group = 0
    !! Velocity model
    !nL = 4
    !im_thick(1:nL) = (/ 20.d0, 50.d0, 90.d0, 1000.d0 /)
    !im_vs(1:nL) = (/ 200.d0, 450.d0, 1000.d0, 2000.d0 /)
    !im_vp(1:nL) = (/ 360.d0, 810.d0, 1800.d0, 3600.d0 /)
    !im_rho(1:nL) = (/ 1800.d0, 1950.d0, 2000.d0, 2700.d0 /)
    !! Frequency settings
    !dL = 20
    !do i=1,dL
    !  freq(i)=0.5+((i-1)*0.5)
    !enddo
    !! Write file
    !ifile = 31
    !filename=trim(path)//'.npin'
    !open(unit=ifile,form='unformatted',access='stream',file=filename,action='write',status='replace')
    !write(ifile) RorL,mode,group ! 3X integer
    !write(ifile) nL,dL ! 2X integer
    !write(ifile) im_thick(1:nL) ! nL real(8)
    !write(ifile) im_vs(1:nL) ! nL real(8)
    !write(ifile) im_vp(1:nL) ! nL real(8)
    !write(ifile) im_rho(1:nL) ! nL real(8)
    !write(ifile) freq(1:dL) ! dL real(8)
    !flush(ifile)
    !close(ifile)
    !endif
    
    !----------------------------------
    ! Read .npin binary file
    ifile = 31
    filename=trim(path)//'.npin'
    open(unit=ifile,form='unformatted',access='stream',file=filename,action='read',iostat=ios)
    if(ios.ne.0)then
      write(*,*) 'ERROR npdc: Cannot open .npin file'
      stop
    endif
    read(ifile,iostat=ios) RorL,mode,group
    if(ios.ne.0)then
      write(*,*) 'ERROR npdc: Cannot read line1 of .npin file'
      stop
    endif
    read(ifile,iostat=ios) nL,dL
    if(ios.ne.0)then
      write(*,*) 'ERROR npdc: Cannot read line2 of .npin file'
      stop
    endif
    read(ifile,iostat=ios) im_thick(1:nL)
    if(ios.ne.0)then
      write(*,*) 'ERROR npdc: Cannot read line3 of .npin file'
      stop
    endif
    read(ifile,iostat=ios) im_vs(1:nL)
    if(ios.ne.0)then
      write(*,*) 'ERROR npdc: Cannot read line4 of .npin file'
      stop
    endif
    read(ifile,iostat=ios) im_vp(1:nL)
    if(ios.ne.0)then
      write(*,*) 'ERROR npdc: Cannot read line5 of .npin file'
      stop
    endif
    read(ifile,iostat=ios) im_rho(1:nL)
    if(ios.ne.0)then
      write(*,*) 'ERROR npdc: Cannot read line6 of .npin file'
      stop
    endif
    read(ifile,iostat=ios) freq(1:dL)
    if(ios.ne.0)then
      write(*,*) 'ERROR npdc: Cannot read line7 of .npin file'
      stop
    endif
    close(ifile)
    
    !----------------------------------
    ! Check values
    if(RorL.ne.0 .and. RorL.ne.1)then
      write(*,*) 'ERROR npdc: Select either Rayleigh (RorL=0) or Love (RorL=1) waves'
      stop
    endif
    if(mode.le.0 .or. mode.gt.4)then
      write(*,*) 'ERROR npdc: Mode number has to be 0<mode<=4'
      stop
    endif
    if(group.ne.0 .and. group.ne.1)then
      write(*,*) 'ERROR npdc: Select either phase (group=0) or group (group=1) velocities'
      stop
    endif
    if(nL.le.0 .or. nL.gt.100)then
      write(*,*) 'ERROR npdc: Number of layers has to be 0<nL<=100'
      stop
    endif
    if(dL.le.0 .or. dL.gt.1000)then
      write(*,*) 'ERROR npdc: Number of frequencies has to be 0<dL<=1000'
      stop
    endif
    
    !----------------------------------
    ! Dispersion curves allocation
    dc = 0.d0
    verbose = 0
    call dispersion_curve_init(verbose)
    
    !----------------------------------
    ! Dispersion curves (Geopsy wrapper)
    do i=1,dL
      omega(i)=6.28318530718d0*freq(i)
    enddo
    
    if(RorL.eq.0)then
      call dispersion_curve_rayleigh(nL,im_thick(1:nL),im_vp(1:nL),im_vs(1:nL),im_rho(1:nL),dL,omega(1:dL),mode,dc,group)
    else
      call dispersion_curve_love(nL,im_thick(1:nL),im_vs(1:nL),im_rho(1:nL),dL,omega(1:dL),mode,dc,group)
    endif
    
    !----------------------------------
    ! Write .npdc output binary file
    ifile = 32
    filename=trim(path)//'.npdc'
    open(unit=ifile,form='unformatted',access='stream',file=filename,action='write',status='replace',iostat=ios)
    if(ios.ne.0)then
      write(*,*) 'ERROR npdc: Cannot write .npdc file'
      stop
    endif
    write(ifile,iostat=ios) dc(((mode-1)*dL)+1:(mode*dL)) ! dL real(8)
    if(ios.ne.0)then
      write(*,*) 'ERROR npdc: Cannot write dc into .npdc file'
      stop
    endif
    flush(ifile)
    close(ifile)
    
!---------------------------------------------------------------------
END PROGRAM


