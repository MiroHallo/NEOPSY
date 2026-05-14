! Initiation from INPUT files - MTI
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

MODULE InputBox
!---------------------------------------------------------------------
!  Module containing parameters from input files
!---------------------------------------------------------------------
    implicit none
    integer:: nfix,nmax,nzone
    real(8):: ptrans,difd,thr_d(2),stp_d,lvz
    real(8),allocatable,dimension(:):: zone_d
    integer,allocatable,dimension(:):: zone_Ntmp
    real(8),allocatable,dimension(:,:):: thr_vs,thr_vp,thr_rho,thr_nu,stp_spr
    real(8),allocatable,dimension(:,:):: fix_dspr
    integer:: nchainsG,iburnG,nstepsG,fix_N
    character(len=80):: dirG
    integer:: DCRm(4),DCLm(4),ELLm(2),ExpertSwitch(4)
    character(len=255):: dirDCR0,dirDCR1,dirDCR2,dirDCR3,dirDCL0,dirDCL1,dirDCL2,dirDCL3
    character(len=255):: dirELL,dirELA,dirINI,dirREF
!---------------------------------------------------------------------
END MODULE



SUBROUTINE InitInput()
!---------------------------------------------------------------------
!  Read parameters from input.dat and and prepare parameters
!---------------------------------------------------------------------
    use InputBox
    use LibBox
    implicit none
    integer:: i,j,k
    real(8):: nu,tmpA,tmpB
    character(len=255):: filename,para1char,para2char
    logical:: file_exists
    character(len=1):: dir_tmp
    integer:: i_tmp
    
    !----------------------------------
    !  Input command-line arguments
    if(COMMAND_ARGUMENT_COUNT().NE.2)then
      write(*,*)'ERROR 101: Two command-line arguments required!'
      stop
    endif
       !read command-line arguments
    CALL GET_COMMAND_ARGUMENT(1,para1char)
    CALL GET_COMMAND_ARGUMENT(2,para2char)
    para1char=trim(para1char)
    para2char=trim(para2char)
    !  Check if input files OK
    inquire(file=para1char,exist=file_exists)
    if(.not.file_exists)then
      write(*,*) 'ERROR 102: The input file 1 does not exist!'
      stop
    endif
    inquire(file=para2char,exist=file_exists)
    if(.not.file_exists)then
      write(*,*) 'ERROR 103: The input file 2 does not exist!'
      stop
    endif
    
    !----------------------------------
    !  Read parameters from input file
    open(10,file=para1char,action='read',status='old')
    read(10,*); read(10,*); read(10,*); read(10,*)
    read(10,*) thr_d(1:2)
    read(10,*)
    read(10,*)
    read(10,*) nzone
    if(nzone.le.0)then
      write(*,*) 'ERROR 104: The number of zones must be a positive integer!'
      stop
    else
      allocate(zone_d(nzone))
      allocate(thr_vs(2,nzone))
      allocate(thr_vp(2,nzone))
      allocate(thr_rho(2,nzone))
      allocate(thr_nu(2,nzone))
      allocate(stp_spr(3,nzone))
    endif
    read(10,*)
    read(10,*)
    read(10,*) zone_d(1:nzone)
    read(10,*)
    read(10,*)
    do i=1,nzone
      read(10,*) thr_vs(1:2,i)
    enddo
    read(10,*)
    read(10,*)
    do i=1,nzone
      read(10,*) thr_vp(1:2,i)
    enddo
    read(10,*)
    read(10,*)
    do i=1,nzone
      read(10,*) thr_rho(1:2,i)
    enddo
    read(10,*)
    read(10,*)
    do i=1,nzone
      read(10,*) thr_nu(1:2,i)
    enddo
    read(10,*)
    read(10,*)
    do i=1,nzone
      read(10,*) stp_spr(1:3,i)
    enddo
    read(10,*)
    read(10,*)
    read(10,*) stp_d
    read(10,*)
    read(10,*)
    read(10,*) lvz
    read(10,*)
    read(10,*)
    read(10,*) nchainsG,iburnG,nstepsG
    read(10,*)
    read(10,*)
    read(10,'(A)') dirG
    read(10,*)
    read(10,*)
    read(10,*) nmax
    read(10,*)
    read(10,*)
    read(10,*) ptrans
    read(10,*)
    read(10,*)
    read(10,*) ExpertSwitch(1:4)
    read(10,*)
    read(10,*)
    read(10,*) fix_N
    read(10,*)
    read(10,*)
    read(10,*) nfix
    read(10,*)
    allocate(fix_dspr(4,nfix))
    do i=1,nfix
      read(10,*) fix_dspr(1:4,i)
    enddo
    close(10)
    
    !----------------------------------
    ! Read parameters from data input file
    open(20,file=para2char,action='read',status='old')
    read(20,*); read(20,*); read(20,*); read(20,*)
    read(20,*) DCRm(1:4)
    read(20,*)
    read(20,*)
    read(20,*) DCLm(1:4)
    read(20,*)
    read(20,*)
    read(20,*) ELLm(1:2)
    read(20,*)
    read(20,*)
    read(20,*)
    read(20,'(A)') dirDCR0
    read(20,'(A)') dirDCR1
    read(20,'(A)') dirDCR2
    read(20,'(A)') dirDCR3
    read(20,*)
    read(20,'(A)') dirDCL0
    read(20,'(A)') dirDCL1
    read(20,'(A)') dirDCL2
    read(20,'(A)') dirDCL3
    read(20,*)
    read(20,'(A)') dirELL
    read(20,'(A)') dirELA
    read(20,*)
    read(20,*)
    read(20,'(A)') dirINI
    read(20,*)
    read(20,*)
    read(20,'(A)') dirREF
    close(20)
    
    !----------------------------------
    ! Zones thickness to depth
    tmpA = zone_d(1)
    zone_d(1) = 0.d0
    do i=2,nzone
      tmpB = zone_d(i)
      zone_d(i) = zone_d(i-1) + tmpA
      tmpA = tmpB
    enddo
    allocate(zone_Ntmp(nzone))
    zone_Ntmp = 0
    
    !----------------------------------
    ! Minimum depth difference [m] between two initial nuclei (it is free during inverson)
    difd = 0.1d0
    
    !----------------------------------
    !  Input corrections/errors
    if(nzone.gt.nmax)then
      nmax = nzone+1
      write(*,'(A,I6)') 'AUTOC 105: The number of mandatory Voronoi nuclei exceeds nmax-1! NEW: ',nmax
    endif
    if(thr_d(1).le.0.d0)then
      thr_d(1) = difd
      write(*,'(A,F6.3)') 'AUTOC 106: The minimum layer-interface depth has to be positive! NEW: ',thr_d(1)
    endif
    if(thr_d(1).gt.thr_d(2))then
      write(*,*) 'ERROR 107: The minimum layer-interface depth exceeds max depth!'
      stop
    endif
    if(ptrans.gt.0.45d0)then
      ptrans = 0.45d0
      write(*,'(A,F6.3)') 'AUTOC 108: The trans-D proposal probability exceeds 0.45! NEW: ',ptrans
    endif
    if(nfix.gt.0 .and. ptrans.gt.0.d0)then
      ptrans = 0.d0
      write(*,*) 'AUTOC 125: Inversion with fixed Voronoi nuclei requires trans-D proposal probability = 0.0!'
    elseif(fix_N.gt.0 .and. ptrans.gt.0.d0)then
      ptrans = 0.d0
      write(*,*) 'AUTOC 126: Inversion with fixed number of layers requires trans-D proposal probability = 0.0!'
    endif
    if(ptrans.le.0.d0 .and. (nfix.le.0 .and. fix_N.le.0))then
      write(*,*) 'ERROR 127: Select inversion type (set trans-D proposal probability or fixed number of layers)!'
      stop
    endif
    if(nzone.ne.1 .and. (nfix.gt.0 .or. fix_N.gt.0))then
      write(*,*) 'ERROR 128: Inversion with a fixed number of layers can be performed in one zone only!'
      stop
    endif
    if(ptrans.le.0.d0 .and. (ExpertSwitch(4).eq.1))then
      write(*,*) 'ERROR 129: Input initial model can be set only for trans-D inversion!'
      stop
    endif
    if(nzone.gt.1 .and. (ExpertSwitch(4).eq.1))then
      write(*,*) 'ERROR 130: Input initial model can be set only for single-zonal inversion!'
      stop
    endif
    if(ptrans.le.0.d0 .and. (nfix.gt.nmax .or. fix_N.gt.nmax)  )then
      nmax = max(nfix,fix_N)+1
      write(*,'(A,I6)') 'AUTOC 109: The number of mandatory Voronoi nuclei exceeds nmax! NEW: ',nmax
    endif
    if(stp_d.gt.(log(thr_d(2))-log(thr_d(1)))/4.d0)then
      write(*,*) 'ERROR 110: Log-Depth step size too large for this profile!'
      stop
    endif
    if(thr_d(2).le.zone_d(nzone))then
      write(*,*) 'ERROR 111: There is defined a zone below the max depth of the profile!'
      stop
    endif
    
    do i=1,nzone
      if(thr_vs(1,i).gt.thr_vs(2,i))then
        tmpA = thr_vs(1,i)
        thr_vs(1,i) = thr_vs(2,i)
        thr_vs(2,i) = tmpA
        write(*,'(A,I6)') 'SWAP 112: The minimum S-wave velocity exceeds its maximum!',i
      endif
      if(thr_vp(1,i).gt.thr_vp(2,i))then
        tmpA = thr_vp(1,i)
        thr_vp(1,i) = thr_vp(2,i)
        thr_vp(2,i) = tmpA
        write(*,'(A,I6)') 'SWAP 113: The minimum P-wave velocity exceeds its maximum!',i
      endif
      if(thr_rho(1,i).gt.thr_rho(2,i))then
        tmpA = thr_rho(1,i)
        thr_rho(1,i) = thr_rho(2,i)
        thr_rho(2,i) = tmpA
        write(*,'(A,I6)') 'SWAP 114: The minimum density exceeds its maximum!',i
      endif
      if(thr_nu(1,i).gt.thr_nu(2,i))then
        tmpA = thr_nu(1,i)
        thr_nu(1,i) = thr_nu(2,i)
        thr_nu(2,i) = tmpA
        write(*,'(A,I6)') 'SWAP 115: The minimum Poisson ratio exceeds its maximum!',i
      endif
      if((thr_vs(2,i)-thr_vs(1,i)).gt.0.d0 .and. stp_spr(1,i).le.0.d0)then
        write(*,'(A,I6)') 'ERROR 116: Zero size of the S-wave velocity step!',i
        stop
      endif
      if((thr_vp(2,i)-thr_vp(1,i)).gt.0.d0 .and. stp_spr(2,i).le.0.d0)then
        write(*,'(A,I6)') 'ERROR 117: Zero size of the P-wave velocity step!',i
        stop
      endif
      if((thr_rho(2,i)-thr_rho(1,i)).gt.0.d0 .and. stp_spr(3,i).le.0.d0)then
        write(*,'(A,I6)') 'ERROR 118: Zero size of the Rho step!',i
        stop
      endif
      if(stp_spr(1,i).gt.(thr_vs(2,i)-thr_vs(1,i))/4.d0)then
        stp_spr(1,i) = (thr_vs(2,i)-thr_vs(1,i))/4.d0
        write(*,'(A,F7.3)') 'AUTOC 119: S-wave velocity step size too large! NEW: ',stp_spr(1,i)
      endif
      if(stp_spr(2,i).gt.(thr_vp(2,i)-thr_vp(1,i))/4.d0)then
        stp_spr(2,i) = (thr_vp(2,i)-thr_vp(1,i))/4.d0
        write(*,'(A,F7.3)') 'AUTOC 120: P-wave velocity step size too large! NEW: ',stp_spr(2,i)
      endif
      if(stp_spr(3,i).gt.(thr_rho(2,i)-thr_rho(1,i))/4.d0)then
        stp_spr(3,i) = (thr_rho(2,i)-thr_rho(1,i))/4.d0
        write(*,'(A,F7.3)') 'AUTOC 121: Density step size too large! NEW: ',stp_spr(3,i)
      endif
    enddo
    
    !----------------------------------
    !  Input warnings
    if(nchainsG.lt.5)then
      write(*,*) 'WARNING: Too few Markov chains! It may result into dependence of sampling models!'
    endif
    if(nchainsG.gt.100)then
      write(*,*) 'WARNING: Too many Markov chains! It may result into slow production of sampling models!'
    endif
    if(iburnG.lt.3000)then
      write(*,*) 'WARNING: Too few burn-in steps! It may result into biased posterior PDF!'
    endif
    if(nstepsG.lt.5000)then
      write(*,*) 'WARNING: Too few chain steps! It may result into biased posterior PDF!'
    endif
    if(stp_d.lt.(log(thr_d(2))-log(thr_d(1)))/1000.d0)then
      write(*,*) 'WARNING: Log-Depth step size is probably too small!'
    endif
    
    do i=1,nzone
      if(stp_spr(1,i).lt.(thr_vs(2,i)-thr_vs(1,i))/10000.d0)then
        write(*,'(A,I6)') 'WARNING: S-wave velocity step size is probably too small!',i
      endif
      if(stp_spr(2,i).lt.(thr_vp(2,i)-thr_vp(1,i))/10000.d0)then
        write(*,'(A,I6)') 'WARNING: P-wave velocity step size is probably too small!',i
      endif
      if(stp_spr(3,i).lt.(thr_rho(2,i)-thr_rho(1,i))/10000.d0)then
        write(*,'(A,I6)') 'WARNING: Rho step size is probably too small!',i
      endif
      if(thr_nu(1,i).lt.0.d0 .or. thr_nu(2,i).gt.0.5d0)then
        write(*,'(A,I6)') 'WARNING: Nu has some strange range. But I can deal with it.',i
      endif
    enddo
    
    do i=1,nfix
      k = binarysearch(nzone,zone_d(1:nzone),fix_dspr(1,i))
      if(k.le.0)then
        k = 1
      endif
      
      if(fix_dspr(4,i).gt.0.d0)then
        if((fix_dspr(2,i).lt.thr_vs(1,k)).or.(fix_dspr(2,i).gt.thr_vs(2,k)))then
          write(*,'(A,I6)') 'WARNING: S-wave velocity of a fixed nuclei exceeds limits!',i
        endif
        if((fix_dspr(3,i).lt.thr_vp(1,k)).or.(fix_dspr(3,i).gt.thr_vp(2,k)))then
          write(*,'(A,I6)') 'WARNING: P-wave velocity of a fixed nuclei exceeds limits!',i
        endif
        if((fix_dspr(4,i).lt.thr_rho(1,k)).or.(fix_dspr(4,i).gt.thr_rho(2,k)))then
          write(*,'(A,I6)') 'WARNING: Density of a fixed nuclei exceeds limits!',i
        endif
        nu = ((fix_dspr(3,i)**2) - 2.d0*(fix_dspr(2,i)**2)) /&
           (2.d0*((fix_dspr(3,i)**2) - (fix_dspr(2,i)**2)))
        if((nu.lt.thr_nu(1,k)).or.(nu.gt.thr_nu(2,k)))then
          write(*,'(A,I6,A1,F6.3)') 'WARNING: Poisson ratio of a fixed nuclei exceeds limits!',i,' /',nu
        endif
      endif
      
      do j=1,nfix
        if(i.eq.j)then
          cycle
        endif
        if(abs(fix_dspr(1,i)-fix_dspr(1,j)).lt.difd)then
          write(*,*) 'ERROR 122: Two of the fixed nuclei have depth difference less than the minimum allowed limit!'
          stop
        endif
      enddo
    enddo
    
    !----------------------------------
    !  Additional info for user
    !if(nfix.gt.0)then
    !  write(*,'(A,I6,A)') 'INFO: Inversion with',nfix,' fixed Voronoi nuclei'
    !endif
    !
    !do i=1,nzone
    !  if((thr_vs(2,i)-thr_vs(1,i)).eq.0.d0)then
    !    write(*,'(A,I6,A)') 'INFO: Fixed S-wave velocity in zone',i,'. I can deal with it, but are you sure?'
    !  endif
    !  if((thr_vp(2,i)-thr_vp(1,i)).eq.0.d0)then
    !    write(*,'(A,I6,A)') 'INFO: Fixed P-wave velocity in zone',i,'. I can deal with it, but are you sure?'
    !  endif
    !  if((thr_rho(2,i)-thr_rho(1,i)).eq.0.d0)then
    !    write(*,'(A,I6,A)') 'INFO: Fixed Rho in zone',i,'. OK then'
    !  endif
    !  if((thr_nu(2,i)-thr_nu(1,i)).eq.0.d0)then
    !    write(*,'(A,I6,A)') 'INFO: Fixed nu in zone',i,'. OK then'
    !  endif
    !enddo
    
    !----------------------------------
    !  Check working directory
    
    dirG = trim(dirG)
    i_tmp = len(trim(dirG(:)))
    dir_tmp = dirG(i_tmp:i_tmp)
    if(dir_tmp.ne.'/')then
      dirG = trim(dirG)//'/'
    endif
    !write(*,*) '-> Working directory: ',trim(dirG)
    inquire(file=dirG,exist=file_exists)
    if(.not.file_exists)then
      write(*,*) 'ERROR 123: The Working Directory does not exist, please create it!'
      stop
    endif
    filename = trim(dirG) // 'log/'
    !write(*,*) '-> Log directory: ',trim(filename)
    inquire(file=filename,exist=file_exists)
    if(.not.file_exists)then
      write(*,*) 'ERROR 124: The Log Directory does not exist, please create it!'
      stop
    endif
    
    !----------------------------------
    !  Prepare data directory
    dirDCR0 = trim(dirDCR0)
    dirDCR1 = trim(dirDCR1)
    dirDCR2 = trim(dirDCR2)
    dirDCR3 = trim(dirDCR3)
    dirDCL0 = trim(dirDCL0)
    dirDCL1 = trim(dirDCL1)
    dirDCL2 = trim(dirDCL2)
    dirDCL3 = trim(dirDCL3)
    dirELL = trim(dirELL)
    dirELA = trim(dirELA)
    dirINI = trim(dirINI)
    dirREF = trim(dirREF)
    
    !----------------------------------
    !  Check if input data are OK
    if((sum(DCRm(1:4))+sum(DCLm(1:4))).lt.1)then
      write(*,*) 'WARNING 201: Include at least one mode of a dispersion curve!'
      !stop
    endif
    if(ELLm(1).gt.0 .and. ELLm(2).gt.0)then
      write(*,*) 'ERROR 202: Include ellipticity OR ellipticity angle!'
      stop
    endif
    if(DCRm(1).gt.0)then
      inquire(file=dirDCR0,exist=file_exists)
      if(.not.file_exists)then
        write(*,*) 'ERROR 203: Datafile for the R0 mode does not exist!'
        stop
      endif
    endif
    if(DCRm(2).gt.0)then
      inquire(file=dirDCR1,exist=file_exists)
      if(.not.file_exists)then
        write(*,*) 'ERROR 204: Datafile for the R1 mode does not exist!'
        stop
      endif
    endif
    if(DCRm(3).gt.0)then
      inquire(file=dirDCR2,exist=file_exists)
      if(.not.file_exists)then
        write(*,*) 'ERROR 205: Datafile for the R2 mode does not exist!'
        stop
      endif
    endif
    if(DCRm(4).gt.0)then
      inquire(file=dirDCR3,exist=file_exists)
      if(.not.file_exists)then
        write(*,*) 'ERROR 206: Datafile for the R3 mode does not exist!'
        stop
      endif
    endif
    
    if(DCLm(1).gt.0)then
      inquire(file=dirDCL0,exist=file_exists)
      if(.not.file_exists)then
        write(*,*) 'ERROR 207: Datafile for the L0 mode does not exist!'
        stop
      endif
    endif
    if(DCLm(2).gt.0)then
      inquire(file=dirDCL1,exist=file_exists)
      if(.not.file_exists)then
        write(*,*) 'ERROR 208: Datafile for the L1 mode does not exist!'
        stop
      endif
    endif
    if(DCLm(3).gt.0)then
      inquire(file=dirDCL2,exist=file_exists)
      if(.not.file_exists)then
        write(*,*) 'ERROR 209: Datafile for the L2 mode does not exist!'
        stop
      endif
    endif
    if(DCLm(4).gt.0)then
      inquire(file=dirDCL3,exist=file_exists)
      if(.not.file_exists)then
        write(*,*) 'ERROR 210: Datafile for the L3 mode does not exist!'
        stop
      endif
    endif
    
    if(ELLm(1).gt.0)then
      inquire(file=dirELL,exist=file_exists)
      if(.not.file_exists)then
        write(*,*) 'ERROR 211: Datafile for the ellipticity does not exist!'
        stop
      endif
    endif
    if(ELLm(2).gt.0)then
      inquire(file=dirELA,exist=file_exists)
      if(.not.file_exists)then
        write(*,*) 'ERROR 212: Datafile for the ellipticity angle does not exist!'
        stop
      endif
    endif
    
    if(ExpertSwitch(4).eq.1)then
      inquire(file=dirINI,exist=file_exists)
      if(.not.file_exists)then
        write(*,*) 'ERROR 213: ASCII file with the initial model does not exist!'
        stop
      endif
    endif
    
    inquire(file=dirREF,exist=file_exists)
    if(.not.file_exists)then
      write(*,*) 'ERROR 215: File with the reference velocity model does not exist!'
      stop
    endif
	
    !----------------------------------
    !write(*,*) '-> Input OK'
    
!---------------------------------------------------------------------
END SUBROUTINE


