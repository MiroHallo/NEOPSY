! Multizonal Transdimensional Inversion - MTI
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

MODULE GlobalBox
!---------------------------------------------------------------------
!  Module containing global variables
!---------------------------------------------------------------------
    implicit none
    ! Global model parameters (last dim for chains)
    real(8),allocatable,dimension(:,:,:):: m_dspr
    integer,allocatable,dimension(:):: m_Nlay
    
    ! PT global control
    integer:: nproc,rank,iseed,pert,ifile,afile
    real(8),allocatable:: logPPDstore(:),misfitstore(:)
    integer:: abins,aaccept(2),aerrors(3),athr(7),aNnuclei(3)
    real(8):: aVR,sVR,sqUvec
    real,external:: ran3
!---------------------------------------------------------------------
END MODULE



PROGRAM  neoinv
!---------------------------------------------------------------------
!  Main program of the Transdimensional Inversion
!---------------------------------------------------------------------
    use GlobalBox
    implicit none
    
    !----------------------------------
    ! Read input files
    call InitInput()
    
    !----------------------------------
    ! Read data
    call InitData()
    
    !----------------------------------
    ! Execute of the inversion
    call InvExe()
    
!---------------------------------------------------------------------
END PROGRAM



SUBROUTINE InvExe()
!---------------------------------------------------------------------
!  Execute of the inversion
!---------------------------------------------------------------------
    use InputBox
    use GlobalBox
    implicit none
    integer:: ialg,nsteps,nchains,iburn,iseed0,nbins
    real(8):: swaprate,tlow,thigh
    character(len=80):: dir
    real(8),allocatable:: Chaintemp(:) ! Temperatures of each chain
    real(8):: Chaintemp_dummy(1)
    character(len=4):: cproc
    
    !----------------------------------
    ! Set up PT control parameters
    ialg = 0            ! Algorithm type:
                        ! ialg=0, PT with T swap between all levels
                        ! ialg=2, PT with T swap only between neighbouring
    nchains = nchainsG  ! Define number of chains
    swaprate = 1.0d0    ! Rate at which exchange swaps are proposed relative to within
                        ! chain steps. 1.0 = one exchage swap proposed for every McMC step.
                        ! Set this value to zero to turn off Parallel Tempering altogether.
    nsteps   = nstepsG  ! Number of chain steps per temperature
    iburn    = iburnG   ! Number of burn in samples
    iseed0 = 61527237   ! Random number seed
    tlow  = 1.d0        ! Lowest temperature of chains
    thigh = 50.d0       ! Highest temperature of chains
    nbins = 100         ! Number of temperature bins for diagnostics
    dir = dirG          ! Home directory for I/O files (some systems want full path names under MPI)
    
    ! Set global variables
    abins = 10         ! Number of steps (as one bin) for accept rate diagnostics
    
    !----------------------------------
    ! Initialize PT library and determine process ID (if parallel)
    ! If we are in serial mode, then rank=0,nproc=1
    call PT(0,0,0,0,0,Chaintemp_dummy,0.d0,0.d0,0,0.d0,0,dir,nproc,rank)
    write(cproc,'(I0)') rank
    
    if(rank.eq.0)then
      write(*,*) "------------------------------------------------------"
      write(*,*) "Multizonal Transdimensional Inversion (Fjord 3/2026)"
      write(*,*) "Miroslav Hallo, SED, ETH Zurich"
      write(*,*) "------------------------------------------------------"
    endif
    
    !----------------------------------
    ! Initial model parameters
    write(*,*)'-> Initial models (node'//trim(cproc)//')'
    call InitMod(nchains,dir)
    
    !----------------------------------
    ! Set temperatures for each chain (array Chaintemp)
    ! Random T-ladder with a log-uniform distribution between tlow and thigh
    write(*,*)'-> Set up temperatures (node'//trim(cproc)//')'
    allocate(Chaintemp(nchains))
    call Setuptempladder(nchains,tlow,thigh,dir,Chaintemp)
    
    !----------------------------------
    ! Calculate transition acceptance rates between chains
    write(*,*)'-> Parallel Tempering (node'//trim(cproc)//')'
    call PT(1,ialg,nchains,nsteps,iburn,Chaintemp,thigh,tlow,&
              nbins,swaprate,iseed0,dir,nproc,rank)
    
    !----------------------------------
    ! Finish MPI and clean up
    write(*,*)'-> Finish and clean up (node'//trim(cproc)//')'
    call PT(99,0,0,0,0,Chaintemp_dummy,0.d0,0.d0,0,0.d0,0,dir,nproc,rank)
    
    close(ifile)
    close(afile)
    deallocate(Chaintemp)
    
    !----------------------------------
    ! Done
    if(rank.eq.0)then
      write(*,*) "Done."
    endif
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE InitMod(nchains,dir)
!---------------------------------------------------------------------
!  Initial model parameters
!---------------------------------------------------------------------
    use InputBox
    use DataBox
    use GlobalBox
    use LibBox
    use ForwardBox
    implicit none
    integer,intent(in):: nchains
    character(len=80),intent(in):: dir
    character(len=4):: cproc
    character(len=255):: filename
    integer:: i,j,im,error,outThr,coun,ai,k
    real(8):: a,sig,temp(5),thr_tmp(2)
    real(8),allocatable,dimension(:):: nu_vec
    logical:: blby
    
    !----------------------------------
    ! Initialize random number generator
    iseed=-(181615141+1000*rank*rank)
    
    !----------------------------------
    ! Allocate model arrays
    ! The last chain (nchains+1) is for the perturbed model
    pert = nchains+1
    allocate(m_dspr(5,nmax,pert))
    allocate(m_Nlay(pert))
    allocate(logPPDstore(nchains))
    allocate(misfitstore(nchains))
    m_dspr=0.d0
    m_Nlay=0
    
    !----------------------------------
    ! Allocate forward im arrays
    allocate(im_thick(nmax))
    allocate(im_vs(nmax))
    allocate(im_vp(nmax))
    allocate(im_rho(nmax))
    im_thick=0.d0
    im_vs=0.d0
    im_vp=0.d0
    im_rho=0.d0
    
    !----------------------------------
    ! Allocate local arrays
    allocate(nu_vec(nmax))
    
    !----------------------------------
    ! Open output file for diagnostics
    write(cproc,'(I0)') rank
    filename=trim(dir)//'log/init_'//trim(cproc)//'.log'
    open(unit=48,file=filename,action='write',status='replace')
    write(48,'(A,I8)') 'Initial models (one node), number of Markov chains: ',nchains
    write(48,'(A)') 'Legend:'
    write(48,'(A)') '  Thr[1] - Initial number of layers reached a threshold'
    write(48,'(A)') '  Thr[2] - Threshold reached in depth of Voronoi nuclei'
    write(48,'(A)') '  Thr[3] - Threshold reached in Vs'
    write(48,'(A)') '  Thr[4] - Threshold reached in combination of [Vs,Vp,nu]'
    write(48,'(A)') '  Thr[5] - Threshold reached in Vp'
    write(48,'(A)') '  Thr[6] - Threshold reached in Rho'
    write(48,'(A)') '  Thr[7] - Low-velocity zone present where prohibited'
    write(48,'(A)') '  Err[1] - Geopsy does not produce results'
    write(48,'(A)') '  Err[2] - Geopsy reduced the number of samples'
    write(48,'(A)') '  Err[*] - Generic (undefined) forward problem error'
    write(48,'(A)') '------------------------------------------------------------'
    
    if(ExpertSwitch(4).eq.1)then
      !----------------------------------
      ! Initial model parameters from file
      write(48,'(A)') 'Initial models from file:'
      write(48,'(A)') '  '//trim(dirINI)
      if(ptrans.le.0.d0)then
        write(*,'(A)') 'ERROR 304: Input initial model can be set only for trans-D inversion!'
        stop
      endif
      if(nzone.gt.1)then
        write(*,'(A)') 'ERROR 305: Input initial model can be set only for single-zonal inversion!'
        stop
      endif
      if(init_nL.le.0)then
        write(*,'(A)') 'ERROR 306: Corrupted input initial model (<1 layers)!'
        stop
      endif
      if(init_nL.gt.1 .and. thr_d(1).gt.init_model(1,1))then
        write(*,'(A)') 'ERROR 307: Input initial model has too thin uppermost layer for the set depth thresholds!'
        stop
      endif
      
      im = 1
      outThr = 0
      m_dspr(:,:,im) = 0.d0
      
      ! Assign nuclei depth
      m_Nlay(im) = init_nL
      if(init_nL.eq.1)then
        sig = (log(thr_d(2))-log(thr_d(1)))*0.5d0
        a = log(thr_d(1)) + sig
        a = exp(a)
        if(a<(thr_d(1)))then
          a=thr_d(1)
          outThr = 2
        elseif(a>thr_d(2))then
          a=thr_d(2)
          outThr = 2
        endif
        m_dspr(1,1,im) = a
      else
        a = max((init_model(1,1)-(init_model(1,1)/2.d0)),thr_d(1))
        m_dspr(1,1,im) = a
        sig = 0.d0
        do i=2,m_Nlay(im)
          sig = sig + init_model(1,i-1)
          m_dspr(1,i,im) = sig + ( sig - m_dspr(1,i-1,im) )
          if(m_dspr(1,i-1,im)>sig)then
            outThr = 2
          endif
        enddo
      endif
      
      k = 1
      do i=1,m_Nlay(im) 
        ! Assign nuclei Vs
        a = init_model(3,i)
        if(a<thr_vs(1,k))then
          a=thr_vs(1,k)
          outThr = 3
        elseif(a>thr_vs(2,k))then
          a=thr_vs(2,k)
          outThr = 3
        endif
        m_dspr(2,i,im) = a
        ! Assign nuclei Vp
        a = init_model(2,i)
        if(a<thr_vp(1,k))then
          a=thr_vp(1,k)
          outThr = 5
        elseif(a>thr_vp(2,k))then
          a=thr_vp(2,k)
          outThr = 5
        endif
        m_dspr(3,i,im) = a
        ! Check nu
        a = ((m_dspr(3,i,im)**2) - 2.d0*(m_dspr(2,i,im)**2)) /&
                (2.d0*((m_dspr(3,i,im)**2) - (m_dspr(2,i,im)**2)))
        if(a<thr_nu(1,k))then
          a=thr_nu(1,k)
          outThr = 4
        elseif(a>thr_nu(2,k))then
          a=thr_nu(2,k)
          outThr = 4
        endif  
        ! Assign nuclei Rho
        a = init_model(4,i)
        if(a<thr_rho(1,k))then
          a=thr_rho(1,k)
          outThr = 6
        elseif(a>thr_rho(2,k))then
          a=thr_rho(2,k)
          outThr = 6
        endif
        m_dspr(4,i,im) = a
      enddo
        
      ! Check order
      if(m_Nlay(im).gt.1)then
        do i=2,m_Nlay(im)
          if((m_dspr(1,i,im)-m_dspr(1,i-1,im)).le.0.d0)then
            outThr = 2
          endif
        enddo
      endif
      
      ! Check low-velocity zones
      if(m_Nlay(im).gt.1)then
        do i=2,m_Nlay(im)
          if(m_dspr(1,i,im).gt.lvz)then
            if((m_dspr(2,i,im)-m_dspr(2,i-1,im)).lt.0.d0 .or. (m_dspr(3,i,im)-m_dspr(3,i-1,im)).lt.0.d0)then
              outThr = 7
            endif
          endif
        enddo
      endif
      
      ! Check if any model parameter out of range
      if(outThr.ne.0)then
        write(48,'(A6,I1,A)') 'Thr[',outThr,'] in the input initial model'
        flush(48)
        write(*,'(A)') 'ERROR 308: Input initial model is not compatible with parameterization!'
        stop
      endif
      
      ! Calculate -logPPD of initial model
      if(ExpertSwitch(1).eq.1)then
        call InterpMod(im,error)
        misfitstore(im) = 0.d0
        logPPDstore(im) = 0.d0
      else
        call Forward(im,misfitstore(im),logPPDstore(im),error)
      endif
      
      ! Check if any error
      if(error.ne.0)then
        write(48,'(A6,I1,A)') 'Err[',error,'] in the input initial model'
        flush(48)
        write(*,'(A)') 'ERROR 309: Input initial model is not compatible with parameterization or observed data!'
        stop
      endif
      
      ! Copy input initial model to all chains
      do im=2,nchains
        m_Nlay(im) = m_Nlay(1)
        m_dspr(:,:,im) = m_dspr(:,:,1)
        misfitstore(im) = misfitstore(1)
        logPPDstore(im) = logPPDstore(1)
      enddo
      
    else
      write(48,'(A)') 'Generation of random initial models:'
      !----------------------------------
      ! Initial random model parameters
      do im=1,nchains
        coun = 0
        do while(.true.) ! until no error
          
          outThr = 0
          blby = .FALSE.
        
          ! Avoid infinite cycle
          coun = coun+1
          if(coun.gt.9999)then
            write(*,*) 'ERROR 301: Intervals [min,max] are too tight or a problem in Geopsy forward problem!'
            stop
          endif
          
          ! Init random number of Voronoi nuclei (log-uniform)
          if(ptrans.gt.0.d0)then
            sig = log(dble(nmax-nzone))
            ai = int(exp(log(dble(nzone))+ran3(iseed)*sig) + 0.5d0)
            ! Cycle if out of range
            if((ai.lt.nzone).or.(ai.gt.nmax))then
              outThr = 1
              write(48,'(A6,I1,A11,I4)') 'Thr[',outThr,'] in chain:',im
              flush(48)
              cycle
            endif
          else
            if(nfix.gt.0)then
              ai = nfix
            elseif(fix_N.gt.0)then
              ai = fix_N
            else
              write(*,*) 'ERROR 302: The fixed number of layers is not set!'
              stop
            endif
            if(ai.gt.nmax)then
              write(*,*) 'ERROR 303: The fixed number of layers is greather than the nmax!'
              stop
            endif
          endif
          
          ! Init all random Voronoi nucleus
          m_Nlay(im) = ai
          
          do i=1,ai
          
            if(ptrans.gt.0.d0)then
              ! Assign zone Voronoi nuclei depth (uniform)
              if(i.le.nzone)then
                if(i.lt.nzone)then
                  sig = zone_d(i+1)-zone_d(i)
                  a = zone_d(i) + ran3(iseed)*sig
                  if(a.le.zone_d(i))then
                    a=zone_d(i)
                    outThr = 2
                  elseif(a.ge.zone_d(i+1))then
                    a=zone_d(i+1)
                    outThr = 2
                  endif
                else
                  sig = thr_d(2)-zone_d(i)
                  a = zone_d(i) + ran3(iseed)*sig
                  if(a.le.zone_d(i))then
                    a=zone_d(i)
                    outThr = 2
                  elseif(a.ge.thr_d(2))then
                    a=thr_d(2)
                    outThr = 2
                  endif
                endif
              ! Assign free Voronoi nuclei depth (log-uniform)
              else
                sig = (log(thr_d(2))-log(thr_d(1)))
                a = log(thr_d(1)) + ran3(iseed)*sig
                a = exp(a)
              endif
              if(a<(thr_d(1)))then
                a=thr_d(1)
                outThr = 2
              elseif(a>thr_d(2))then
                a=thr_d(2)
                outThr = 2
              endif
              m_dspr(1,i,im) = a
            
            else
              ! Assign fixed Voronoi nuclei
              if(nfix.gt.0)then
                if(fix_dspr(4,i).gt.0.d0)then
                  ! A-type
                  m_dspr(1:4,i,im) = fix_dspr(1:4,i)
                  m_dspr(5,i,im) = 123.d0 ! Flag
                  cycle
                else
                  ! D-type
                  m_dspr(1,i,im) = fix_dspr(1,i)
                  m_dspr(5,i,im) = -123.d0 ! Flag
                endif
              ! Assign depths of the fixed number of layers (log-uniform)
              elseif(fix_N.gt.0)then
                sig = (log(thr_d(2))-log(thr_d(1)))
                a = log(thr_d(1)) + ran3(iseed)*sig
                a = exp(a)
                if(a<(thr_d(1)))then
                  a=thr_d(1)
                  outThr = 2
                elseif(a>thr_d(2))then
                  a=thr_d(2)
                  outThr = 2
                endif
                m_dspr(1,i,im) = a
              else
                write(*,*) 'ERROR 302: The fixed number of layers is not set!'
                stop  
              endif
              
            endif
          
            ! Find the zone index of the nucleus
            k = binarysearch(nzone,zone_d(1:nzone),m_dspr(1,i,im))
            if(k.le.0)then
             k = 1
            endif
            
            ! Init random Vs (uniform)
            sig = (thr_vs(2,k)-thr_vs(1,k))
            a = thr_vs(1,k) + ran3(iseed)*sig
            if(a<thr_vs(1,k))then
              a=thr_vs(1,k)
              outThr = 3
            elseif(a>thr_vs(2,k))then
              a=thr_vs(2,k)
              outThr = 3
            endif
            m_dspr(2,i,im) = a
            
            ! Init random Vp (uniform but with fitting nu)
            thr_tmp(1) = (2.d0*(m_dspr(2,i,im)**2)*(1.d0-thr_nu(1,k))) / (1.d0-(2.d0*thr_nu(1,k)))
            thr_tmp(2) = (2.d0*(m_dspr(2,i,im)**2)*(1.d0-thr_nu(2,k))) / (1.d0-(2.d0*thr_nu(2,k)))
            thr_tmp(1) = sqrt( max(0.d0, thr_tmp(1)) )
            thr_tmp(2) = sqrt( max(0.d0, thr_tmp(2)) )
            if(thr_tmp(1).lt.thr_vp(1,k))then
              thr_tmp(1) = thr_vp(1,k)
            endif
            if(thr_tmp(2).gt.thr_vp(2,k))then
              thr_tmp(2) = thr_vp(2,k)
            endif
            if(thr_tmp(1).gt.thr_tmp(2))then
              outThr = 5
            endif
            a = thr_tmp(1) + ran3(iseed)*(thr_tmp(2)-thr_tmp(1))
            if(a<thr_vp(1,k))then
              a=thr_vp(1,k)
              outThr = 4
            elseif(a>thr_vp(2,k))then
              a=thr_vp(2,k)
              outThr = 4
            endif
            m_dspr(3,i,im) = a
            
            ! Init random Rho (uniform)
            sig = (thr_rho(2,k)-thr_rho(1,k))
            a = thr_rho(1,k) + ran3(iseed)*sig
            if(a<thr_rho(1,k))then
              a=thr_rho(1,k)
              outThr = 6
            elseif(a>thr_rho(2,k))then
              a=thr_rho(2,k)
              outThr = 6
            endif
            m_dspr(4,i,im) = a
          enddo
        
          ! Cycle if any model parameter out of range
          if(outThr.ne.0)then
            write(48,'(A6,I1,A11,I4)') 'Thr[',outThr,'] in chain:',im
            flush(48)
            m_dspr(:,:,im) = 0.d0
            cycle
          endif
          
          ! Bubble sorting by depth
          do i=1,m_Nlay(im)
            do j=m_Nlay(im),i+1,-1
              if(m_dspr(1,j-1,im)>m_dspr(1,j,im))then
                temp(1:5) = m_dspr(1:5,j-1,im)
                m_dspr(1:5,j-1,im) = m_dspr(1:5,j,im)
                m_dspr(1:5,j,im) = temp(1:5)
              endif
            enddo
          enddo
        
          ! Check distances
          if(m_Nlay(im).gt.1)then
            do i=2,m_Nlay(im)
              if(abs(m_dspr(1,i,im)-m_dspr(1,i-1,im)).lt.difd)then
                outThr = 2
                write(48,'(A6,I1,A11,I4)') 'Thr[',outThr,'] in chain:',im
                flush(48)
                m_dspr(:,:,im) = 0.d0
                blby = .TRUE.
                exit
              endif
            enddo
          endif
          ! It is bad model, cycle to the next
          if(blby)then
            cycle
          endif
        
          ! Check low-velocity zones (Desert)
          if(m_Nlay(im).gt.1)then
            do i=2,m_Nlay(im)
              if(m_dspr(1,i,im).gt.lvz)then
                if((m_dspr(2,i,im)-m_dspr(2,i-1,im)).lt.0.d0 .or. (m_dspr(3,i,im)-m_dspr(3,i-1,im)).lt.0.d0)then
                  outThr = 7
                  write(48,'(A6,I1,A11,I4)') 'Thr[',outThr,'] in chain:',im
                  flush(48)
                  m_dspr(:,:,im) = 0.d0
                  blby = .TRUE.
                  exit
                endif
              endif
            enddo
          endif
          ! It is bad model, cycle to the next
          if(blby)then
            cycle
          endif
          
          ! Calculate -logPPD of initial model
          if(ExpertSwitch(1).eq.1)then
            call InterpMod(im,error)
            misfitstore(im) = 0.d0
            logPPDstore(im) = 0.d0
          else
            call Forward(im,misfitstore(im),logPPDstore(im),error)
          endif
          
          ! Cycle if any error
          if(error.ne.0)then
            write(48,'(A6,I1,A11,I4)') 'Err[',error,'] in chain:',im
            flush(48)
            m_dspr(:,:,im) = 0.d0
            cycle
          endif
          
          ! If no error
          if(error.eq.0.and.outThr.eq.0)then
            exit
          endif
          
        enddo ! until no error
      enddo
    endif
    
    !----------------------------------
    ! Prepare vars for diagnostics
    aaccept = 0
    aerrors = 0
    athr = 0
    aNnuclei = 0
    aVR = -99.99d0
    sVR = -99.99d0
    sqUvec = N_used
    
    ! Write resultant initial models
    write(48,'(A)') '------------------------------------------------------------'
    write(48,'(A)') 'Final initial models:'
    do im=1,nchains
      call InterpMod(im,error)
      if(error.ne.0)then
        write(*,'(A,I4)') 'ERROR: Unable to interpret a final initial model in chain:',im
        stop
      endif
      write(48,'(A18,100000F10.1)') 'Nucleus Depth [m]:',m_dspr(1,1:m_Nlay(im),im)
      write(48,'(A18,100000F10.1)') 'Layer thick. [m]:',im_thick(1:m_Nlay(im))
      write(48,'(A18,100000F10.1)') 'Vs [m/s]:',m_dspr(2,1:m_Nlay(im),im)
      write(48,'(A18,100000F10.1)') 'Vp [m/s]:',m_dspr(3,1:m_Nlay(im),im)
      write(48,'(A18,100000F10.1)') 'Rho [kg/m3]:',abs(m_dspr(4,1:m_Nlay(im),im))
      do i=1,m_Nlay(im)
        nu_vec(i)=((m_dspr(3,i,im)**2) - 2.d0*(m_dspr(2,i,im)**2)) /&
                  (2.d0*((m_dspr(3,i,im)**2) - (m_dspr(2,i,im)**2)))
      enddo
      write(48,'(A18,100000F10.4)') 'Poisson ratio:',nu_vec(1:m_Nlay(im))
      write(48,'(A18,F10.1)') 'VR [%]:',(1-(misfitstore(im)/sqUvec))*100.d0
      write(48,'(A)') '----------------------------------------'
    enddo
    close(48)
    
    !----------------------------------
    ! Open output file for models
    ifile = 55
    filename=trim(dir)//'xmodels_'//trim(cproc)//'.dat'
    open(unit=ifile,form='unformatted',access='stream',file=filename,action='write',status='replace')
    flush(ifile)
    
    !----------------------------------
    ! Open output file for diagnostics
    afile = 69
    filename=trim(dir)//'log/stats_'//trim(cproc)//'.log'
    open(unit=afile,file=filename,action='write',status='replace')
    write(afile,'(A)') 'Runtime log from MCMC sampling (one node)'
    write(afile,'(A)') 'Legend:'
    write(afile,'(A)') '  MCstep - Number of already performed steps of Markov chains'
    write(afile,'(A)') '  AccT1 - Number of new accepted models on the sampling temperature'
    write(afile,'(A)') '  AccAll - Number of new accepted models on all temperatures'
    write(afile,'(A)') '  Err[1] - Number of cases when Geopsy did not produced results'
    write(afile,'(A)') '  Err[2] - Number of cases when Geopsy reduced the number of samples'
    write(afile,'(A)') '  Err[*] - Number of occurrences of an undefined forward problem error'
    write(afile,'(A)') '  aVR[%] - The highest data Variance Reduction of new accepted models from all chains'
    write(afile,'(A)') '  sVR[%] - The highest data Variance Reduction of models assigned to the sampling chain'
    write(afile,'(A)') '          (100% .. perfect fit, 0% .. 1sigma fit, -300% .. 2sigma fit, -9999 .. NaN)'
    write(afile,'(A)') '  DimAll - Average number of layers in proposed models on all temperatures'
    write(afile,'(A)') '  Dim[-] - Fraction of death move proposals (proportionate part of all moves)'
    write(afile,'(A)') '  Dim[+] - Fraction of birth move proposals (proportionate part of all moves)'
    write(afile,'(A)') '  Thr[1] - Fraction of cases when a threshold reached in proposed depth'
    write(afile,'(A)') '  Thr[2] - Fraction of cases when a threshold reached in proposed order of layers'
    write(afile,'(A)') '  Thr[3] - Fraction of cases when a threshold reached in proposed Vs'
    write(afile,'(A)') '  Thr[4] - Fraction of cases when a threshold reached in proposed Vp'
    write(afile,'(A)') '  Thr[5] - Fraction of cases when a threshold reached in proposed [Vs,Vp,nu] trinity!!!'
    write(afile,'(A)') '  Thr[6] - Fraction of cases when a threshold reached in proposed Rho'
    write(afile,'(A)') '  Thr[7] - Fraction of cases when a low-velocity zone occurred where prohibited'
    write(afile,'(A)') '----------------------------------------------------------------------------'
    write(afile,'(6A8,2A10,10A8)') 'MCstep','AccT1','AccAll','Err[1]','Err[2]','Err[*]','aVR[%]','sVR[%]', &
         'DimAll','Dim[-]','Dim[+]','Thr[1]','Thr[2]','Thr[3]','Thr[4]','Thr[5]','Thr[6]','Thr[7]'
    write(afile,'(A8,5I8,2A10,I8,9A8)')'max(val)',abins,abins*nchains,abins*nchains,abins*nchains,&
                           abins*nchains,'100%','100%',nmax,'100%','100%','100%','100%','100%','100%','100%','100%','100%'
    write(afile,'(A)') '----------------------------------------------------------------------------'
    flush(afile)
    
    !----------------------------------
    ! Deallocate local arrays
    deallocate(nu_vec)
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE AdvanceChain(ichain,is,Temper,logPPD)
!---------------------------------------------------------------------
!  User routine to advance a Metropolis-Hastings random 
!  walker using Markov chain Monte Carlo.
!  In:   ichain     Integer   : Random walk (chain) identifer
!        is         Integer   : Main loop ID over chain steps (iburn+nsteps)
!        Temper     Doubler   : Temperature of chain supplied by the calling routine
!  Out:  logPPD     Real*8    : Negative logrithm of the posterior
!                               probability density function of updated chain.
!---------------------------------------------------------------------
    use InputBox
    use GlobalBox
    use LibBox
    use ForwardBox
    implicit none
    integer,intent(in):: ichain,is
    real(8),intent(in):: Temper
    real(8),intent(out):: logPPD
    
    integer:: i,j,k,kz,kz2,error,errInt,freeSum
    real(8):: misfit,a,sig,nu,logPPD1,logPPD2,mai
    real(8):: logQ12,logQ21,q12,q21,tmp,temp(5),thr_tmp(2),thr_vptmp(2)
    logical:: yn,blby
    character(len=4):: aproc
    !character(len=7):: cproc
    !character(len=255):: filename
    
    !----------------------------------
    !  Get negative log PPD of current model
    logPPD1 = logPPDstore(ichain)
    logPPD = logPPD1
    
    m_dspr(:,:,pert) = 0.d0
    
    !----------------------------------
    ! Select move type (death/birth/perturb)
    ! Generic proposal distribution
    q12 = 1.d0  ! q(m1|m2)
    q21 = 1.d0  ! q(m2|m1)
    k = m_Nlay(ichain)
    
    !----------------------------------
    ! Count number of nuclei in zones
    zone_Ntmp = 0
    do i=1,m_Nlay(ichain)
      if(m_dspr(5,i,ichain).eq.0.d0)then ! Free nucleus
        kz = binarysearch(nzone,zone_d(1:nzone),m_dspr(1,i,ichain))
        if(kz.le.0)then
          kz = 1
        endif
        zone_Ntmp(kz) = zone_Ntmp(kz)+1
      endif
    enddo
    
    mai = ran3(iseed)
    if(mai.lt.ptrans)then
      !----------------------------------
      ! Dead move
      !----------------------------------
      m_Nlay(pert) = m_Nlay(ichain)
      m_dspr(:,:,pert) = m_dspr(:,:,ichain)
      
      if(k.gt.nzone .and. m_Nlay(pert).gt.1)then
        ! Proposal distribution
        q12 = dble(k)            ! q(m1|m2)
        q21 = dble(k-1)          ! q(m2|m1)
        ! Reciprocal distribution as an uninformed prior distribution p(k)
        ! The reciprocal distribution has a density function of the form
        ! f(x)~1/x for 0<a<x<b (~ means "is proportional to" -> norm. constant)
        ! The inverse distribution in this case is of the form
        ! p(k)~1/k for 0<=(1/b)<k<=(1/a)
        
        ! Count zones with multiple Voronoi nuclei
        freeSum = 0
        do kz=1,nzone
          if(zone_Ntmp(kz).gt.1)then
            freeSum = freeSum + zone_Ntmp(kz)
          endif
        enddo
        ! Select one free Voronoi nucleus
        j = 1 + int(ran3(iseed)*freeSum) ! j-th free nucleus
        j = min(j,freeSum) ! In very rare case
        do i=1,m_Nlay(pert)
          if(m_dspr(5,i,pert).eq.0.d0)then ! Free nucleus
            kz = binarysearch(nzone,zone_d(1:nzone),m_dspr(1,i,pert))
            if(kz.le.0)then
              kz = 1
            endif
            if(zone_Ntmp(kz).gt.1)then
              j = j - 1
              if(j.lt.1)then
                exit ! i-th nucleus is about to die
              endif
            endif
          endif
        enddo
        ! Dead of the i-th nucleus
        do j=i,m_Nlay(pert)-1
          m_dspr(1:5,j,pert) = m_dspr(1:5,j+1,pert)
        enddo
        m_dspr(1:5,m_Nlay(pert),pert) = 0.d0
        m_Nlay(pert) = m_Nlay(pert) - 1
        aNnuclei(2) = aNnuclei(2) + 1
      endif
      aNnuclei(1) = aNnuclei(1) + m_Nlay(pert)
      
    elseif(mai.lt.(2.d0*ptrans))then
      !----------------------------------
      ! Birth move
      !----------------------------------
      m_Nlay(pert) = m_Nlay(ichain)
      m_dspr(:,:,pert) = m_dspr(:,:,ichain)
      blby = .FALSE.
      
      if(m_Nlay(pert).lt.nmax)then
        ! Proposal distribution
        q12 = dble(k)            ! q(m1|m2)
        q21 = dble(k+1)          ! q(m2|m1)
        
        ! Birth of the i-th nucleus
        i = m_Nlay(pert)+1
        m_Nlay(pert) = i
          
        ! Init random depth (log-uniform)
        sig = (log(thr_d(2))-log(thr_d(1)))
        a = log(thr_d(1)) + ran3(iseed)*sig
        m_dspr(1,i,pert) = exp(a)
          
        ! Find the zone index of the nucleus
        kz = binarysearch(nzone,zone_d(1:nzone),m_dspr(1,i,pert))
        if(kz.le.0)then
          kz = 1
        endif
        
        ! Init random Vs (Desert)
        do j=1,i
          if(m_dspr(1,j,pert).ge.m_dspr(1,i,pert))then
            exit
          endif
        enddo
        if(j.eq.i)then
          if(m_dspr(1,i,pert).le.lvz)then ! Free velocity
            m_dspr(2,i,pert) = thr_vs(1,kz) + ran3(iseed)*(thr_vs(2,kz)-thr_vs(1,kz))
          else ! Velocity higher than in j-1
            m_dspr(2,i,pert) = m_dspr(2,j-1,pert) + ran3(iseed)*(thr_vs(2,kz)-m_dspr(2,j-1,pert))
          endif
        else
          if(m_dspr(1,j,pert).le.lvz)then ! Free velocity
            m_dspr(2,i,pert) = thr_vs(1,kz) + ran3(iseed)*(thr_vs(2,kz)-thr_vs(1,kz))
          elseif(m_dspr(1,i,pert).le.lvz .or. j.eq.1)then ! Velocity lower than in j
            sig = min(thr_vs(2,kz),m_dspr(2,j,pert)) - thr_vs(1,kz)
            if(sig.le.0.d0)then
              blby = .TRUE.
              athr(7)=athr(7)+1
            else
              m_dspr(2,i,pert) = thr_vs(1,kz) + ran3(iseed)*sig
            endif
          else
            sig = min(thr_vs(2,kz),m_dspr(2,j,pert)) - max(thr_vs(1,kz),m_dspr(2,j-1,pert))
            if(sig.le.0.d0)then
              blby = .TRUE.
              athr(7)=athr(7)+1
            else
              m_dspr(2,i,pert) = max(thr_vs(1,kz),m_dspr(2,j-1,pert)) + ran3(iseed)*sig
            endif
          endif
        endif
        
        ! Init random Vp (Desert)
        if(j.eq.i)then
          if(m_dspr(1,i,pert).le.lvz)then ! Free velocity
            thr_vptmp(1) = thr_vp(1,kz)
            thr_vptmp(2) = thr_vp(2,kz)
          else ! Velocity higher than in j-1
            thr_vptmp(1) = m_dspr(3,j-1,pert)
            thr_vptmp(2) = thr_vp(2,kz)
          endif
        else
          if(m_dspr(1,j,pert).le.lvz)then ! Free velocity
            thr_vptmp(1) = thr_vp(1,kz)
            thr_vptmp(2) = thr_vp(2,kz)
          elseif(m_dspr(1,i,pert).le.lvz .or. j.eq.1)then ! Velocity lower than in j
            thr_vptmp(1) = thr_vp(1,kz)
            thr_vptmp(2) = min(thr_vp(2,kz),m_dspr(3,j,pert))
          else
            thr_vptmp(1) = max(thr_vp(1,kz),m_dspr(3,j-1,pert))
            thr_vptmp(2) = min(thr_vp(2,kz),m_dspr(3,j,pert))
          endif
        endif
        ! Init random Vp (thr_nu and finish)
        thr_tmp(1) = (2.d0*(m_dspr(2,i,pert)**2)*(1.d0-thr_nu(1,kz))) / (1.d0-(2.d0*thr_nu(1,kz)))
        thr_tmp(2) = (2.d0*(m_dspr(2,i,pert)**2)*(1.d0-thr_nu(2,kz))) / (1.d0-(2.d0*thr_nu(2,kz)))
        thr_tmp(1) = sqrt( max(0.d0, thr_tmp(1)) )
        thr_tmp(2) = sqrt( max(0.d0, thr_tmp(2)) )
        if(thr_tmp(1).lt.thr_vptmp(1))then
          thr_tmp(1) = thr_vptmp(1)
        endif
        if(thr_tmp(2).gt.thr_vptmp(2))then
          thr_tmp(2) = thr_vptmp(2)
        endif
        if(thr_tmp(1).gt.thr_tmp(2))then
          blby = .TRUE.
          athr(5)=athr(5)+1
        endif
        a = thr_tmp(1) + ran3(iseed)*(thr_tmp(2)-thr_tmp(1))
        m_dspr(3,i,pert) = a
          
        ! Init random Rho (uniform)
        m_dspr(4,i,pert) = thr_rho(1,kz) + ran3(iseed)*(thr_rho(2,kz)-thr_rho(1,kz))
        
        ! Flag
        m_dspr(5,i,pert) = 0.d0
        
        ! Instant death
        if(blby)then
          m_dspr(:,:,pert) = m_dspr(:,:,ichain)
          m_Nlay(pert) = m_Nlay(ichain)
        else
          aNnuclei(3) = aNnuclei(3) + 1
        endif
        
      endif
      aNnuclei(1) = aNnuclei(1) + m_Nlay(pert)
      
    else
      !----------------------------------
      ! Perturb move
      !----------------------------------
      m_Nlay(pert) = m_Nlay(ichain)
      do i=1,m_Nlay(ichain)
        aNnuclei(1) = aNnuclei(1) + 1
      
        ! Find the zone index of the nucleus
        kz = binarysearch(nzone,zone_d(1:nzone),m_dspr(1,i,ichain))
        if(kz.le.0)then
          kz = 1
        endif
        
        if(m_dspr(5,i,ichain).gt.0.d0)then     ! Fixed nucleus of A-type
          m_dspr(1:5,i,pert) = m_dspr(1:5,i,ichain)
        else
          if(m_dspr(5,i,ichain).lt.0.d0)then   ! Fixed nucleus of D-type
            m_dspr(1,i,pert) = m_dspr(1,i,ichain)
            m_dspr(5,i,pert) = m_dspr(5,i,ichain)
          else                                 ! Free nucleus
            m_dspr(5,i,pert) = m_dspr(5,i,ichain)
            
            ! Random depth step (gauss in log-depth)
            a = log(m_dspr(1,i,ichain)) + dble(gasdev(iseed)) * stp_d
            if(a.lt.log(thr_d(1)))then
              a=(2.d0*log(thr_d(1)))-a
              if(a.gt.log(thr_d(2)))then ! For sure
                a=log(m_dspr(1,i,ichain))
              endif
              athr(1)=athr(1)+1
            elseif(a.gt.log(thr_d(2)))then
              a=(2.d0*log(thr_d(2)))-a
              if(a.lt.log(thr_d(1)))then ! For sure
                a=log(m_dspr(1,i,ichain))
              endif
              athr(1)=athr(1)+1
            endif
            
            ! Trans-zone move
            blby = .FALSE.
            kz2 = binarysearch(nzone,zone_d(1:nzone),exp(a))
            if(kz2.le.0)then
              kz2 = 1
            endif
            if(kz2.ne.kz)then
              ! Check parameters in the new zone
              sig = m_dspr(2,i,ichain)
              if(sig<thr_vs(1,kz2) .or. sig>thr_vs(2,kz2))then
                blby = .TRUE.
              endif
              sig = m_dspr(3,i,ichain)
              if(sig<thr_vp(1,kz2) .or. sig>thr_vp(2,kz2))then
                blby = .TRUE.
              endif
              sig = m_dspr(4,i,ichain)
              if(sig<thr_rho(1,kz2) .or. sig>thr_rho(2,kz2))then
                blby = .TRUE.
              endif
              nu = ((m_dspr(3,i,ichain)**2) - 2.d0*(m_dspr(2,i,ichain)**2)) /&
                   (2.d0*((m_dspr(3,i,ichain)**2) - (m_dspr(2,i,ichain)**2)))
              if(nu<thr_nu(1,kz2) .or. nu>thr_nu(2,kz2))then
                blby = .TRUE.
              endif
              ! Check if you can leave the old zone
              if(zone_Ntmp(kz).le.1)then
                blby = .TRUE.
              endif
            endif
            if(blby)then
              m_dspr(1,i,pert) = m_dspr(1,i,ichain)
            else
              if(kz2.ne.kz)then
                ! Multizonal prior PDF ratio (Desert+)
                sig = 1.d0/max(1.d0,(thr_vs(2,kz2)-thr_vs(1,kz2)))
                sig = sig * (1.d0/max(1.d0,(thr_vp(2,kz2)-thr_vp(1,kz2))))
                sig = sig * (1.d0/max(1.d0,(thr_rho(2,kz2)-thr_rho(1,kz2))))
                q12 = q12 * sig          ! p(m2)
                sig = 1.d0/max(1.d0,(thr_vs(2,kz)-thr_vs(1,kz)))
                sig = sig * (1.d0/max(1.d0,(thr_vp(2,kz)-thr_vp(1,kz))))
                sig = sig * (1.d0/max(1.d0,(thr_rho(2,kz)-thr_rho(1,kz))))
                q21 = q21 * sig          ! p(m1)
                zone_Ntmp(kz) = zone_Ntmp(kz)-1
                zone_Ntmp(kz2) = zone_Ntmp(kz2)+1
                kz = kz2
              endif
              m_dspr(1,i,pert) = exp(a)
            endif
          endif
          
          ! Random Vs step (gauss)
          a = m_dspr(2,i,ichain) + dble(gasdev(iseed)) * stp_spr(1,kz)
          if(a<thr_vs(1,kz))then
            a=(2.d0*thr_vs(1,kz))-a
            if(a>thr_vs(2,kz))then ! For sure
              a=m_dspr(2,i,ichain)
            endif
            athr(3)=athr(3)+1
          elseif(a>thr_vs(2,kz))then
            a=(2.d0*thr_vs(2,kz))-a
            if(a<thr_vs(1,kz))then ! For sure
              a=m_dspr(2,i,ichain)
            endif
            athr(3)=athr(3)+1
          endif
          m_dspr(2,i,pert) = a
        
          ! Random Vp step (gauss)
          thr_tmp(1) = (2.d0*(m_dspr(2,i,pert)**2)*(1.d0-thr_nu(1,kz))) / (1.d0-(2.d0*thr_nu(1,kz)))
          thr_tmp(2) = (2.d0*(m_dspr(2,i,pert)**2)*(1.d0-thr_nu(2,kz))) / (1.d0-(2.d0*thr_nu(2,kz)))
          thr_tmp(1) = sqrt( max(0.d0, thr_tmp(1)) )
          thr_tmp(2) = sqrt( max(0.d0, thr_tmp(2)) )
          if(thr_tmp(1).lt.thr_vp(1,kz))then
            thr_tmp(1) = thr_vp(1,kz)
          endif
          if(thr_tmp(2).gt.thr_vp(2,kz))then
            thr_tmp(2) = thr_vp(2,kz)
          endif
          if(thr_tmp(1).gt.thr_tmp(2))then
            m_dspr(2:3,i,pert) = m_dspr(2:3,i,ichain)
            athr(5)=athr(5)+1
          endif
          if( abs(thr_tmp(2)-thr_tmp(1)).lt.1.d0)then
            a = thr_tmp(1) + (thr_tmp(2)-thr_tmp(1))/2.d0
          else
            a = m_dspr(3,i,ichain) + dble(gasdev(iseed)) * stp_spr(2,kz)
          endif
          if(a<thr_tmp(1))then
            a=(2.d0*thr_tmp(1))-a
            if(a>thr_tmp(2))then ! For sure
              a=m_dspr(3,i,ichain)
              m_dspr(2:3,i,pert) = m_dspr(2:3,i,ichain)
            endif
            athr(4)=athr(4)+1
          elseif(a>thr_tmp(2))then
            a=(2.d0*thr_tmp(2))-a
            if(a<thr_tmp(1))then ! For sure
              a=m_dspr(3,i,ichain)
              m_dspr(2:3,i,pert) = m_dspr(2:3,i,ichain)
            endif
            athr(4)=athr(4)+1
          endif
          m_dspr(3,i,pert) = a
        
          ! Random Rho step (gauss)
          a = m_dspr(4,i,ichain) + dble(gasdev(iseed)) * stp_spr(3,kz)
          if(a<thr_rho(1,kz))then
            a=(2.d0*thr_rho(1,kz))-a
            if(a>thr_rho(2,kz))then ! For sure
              a=m_dspr(4,i,ichain)
            endif
            athr(6)=athr(6)+1
          elseif(a>thr_rho(2,kz))then
            a=(2.d0*thr_rho(2,kz))-a
            if(a<thr_rho(1,kz))then ! For sure
              a=m_dspr(4,i,ichain)
            endif
            athr(6)=athr(6)+1
          endif
          m_dspr(4,i,pert) = a
        
        endif
      enddo
      
    endif ! End of moves
    
    !----------------------------------
    ! Bubble sorting by depth
    do i=1,m_Nlay(pert)
      do j=m_Nlay(pert),i+1,-1
        if(m_dspr(1,j-1,pert)>m_dspr(1,j,pert))then
          temp(1:5) = m_dspr(1:5,j-1,pert)
          m_dspr(1:5,j-1,pert) = m_dspr(1:5,j,pert)
          m_dspr(1:5,j,pert) = temp(1:5)
        endif
      enddo
    enddo
    
    !----------------------------------
    ! Check distances between layers (over-parameterization)
    if(m_Nlay(pert).gt.1)then
      do i=2,m_Nlay(pert)
        if(abs(m_dspr(1,i,pert)-m_dspr(1,i-1,pert)).lt.difd)then
          athr(2)=athr(2)+1
          ! But it is still OK, do nothing to preserve prior
          !m_dspr(:,:,pert) = m_dspr(:,:,ichain)
          !m_Nlay(pert) = m_Nlay(ichain)
        endif
      enddo
    endif
    
    !----------------------------------
    ! Check low-velocity zones (Desert)
    if(m_Nlay(pert).gt.1)then
      do i=2,m_Nlay(pert)
        if(m_dspr(1,i,pert).gt.lvz)then
          if((m_dspr(2,i,pert)-m_dspr(2,i-1,pert)).lt.0.d0 .or. (m_dspr(3,i,pert)-m_dspr(3,i-1,pert)).lt.0.d0)then
            athr(7)=athr(7)+1
            m_dspr(:,:,pert) = m_dspr(:,:,ichain)
            m_Nlay(pert) = m_Nlay(ichain)
            exit
          endif
        endif    
      enddo
    endif
    
    !----------------------------------
    ! Calculate -logPPD of proposed model in ichain
    if(ExpertSwitch(1).eq.1)then
      call InterpMod(pert,error)
      misfit = 0.d0
      logPPD2 = 0.d0
    else
      call Forward(pert,misfit,logPPD2,error)
    endif
    
    !----------------------------------
    ! log of probability of moving state
    logQ21 = log(q21)
    logQ12 = log(q12)
    
    !----------------------------------
    ! Decide whether to accept perturbed model using M-H criterion.
    ! raised to the power of inverse temperature.
    if(error.eq.0)then
      call PT_McMC_accept(Temper,logPPD1,logQ12,logPPD2,logQ21,yn)
    else
      yn = .FALSE.
    endif
    
    !----------------------------------
    ! We accept the new model
    if(yn)then
      m_dspr(:,:,ichain) =  m_dspr(:,:,pert)
      m_Nlay(ichain) = m_Nlay(pert)
      logPPDstore(ichain) = logPPD2
      misfitstore(ichain) = misfit
      logPPD = logPPD2
    endif
    
    !----------------------------------
    ! Write out target model
    if(Temper.eq.1.d0 .and. is.gt.iburnG)then
      misfit = misfitstore(ichain)
      call InterpMod(ichain,errInt)
      write(ifile) is ! 1X integer (Eagle)
      write(ifile) logPPD,misfit  ! 2X real(8)
      write(ifile) m_Nlay(ichain) ! 1X integer
      write(ifile) im_thick(1:m_Nlay(ichain))
      write(ifile) im_vs(1:m_Nlay(ichain))
      write(ifile) im_vp(1:m_Nlay(ichain))
      write(ifile) im_rho(1:m_Nlay(ichain))
      flush(ifile)
    endif
    
    !----------------------------------
    ! Screen if start of production
    if(ichain.eq.1 .and. is.eq.iburnG)then
      write(aproc,'(I0)') rank
      write(*,*)'-> Production phase begins (node'//trim(aproc)//')'
    endif
    
    !----------------------------------
    ! Write out diagnostics
    if(yn)then
      aaccept(2)=aaccept(2)+1
      if(Temper.eq.1)then
        aaccept(1)=aaccept(1)+1
      endif
      if(aVR.lt.(1-(misfit/sqUvec)))then
        aVR = 1-(misfit/sqUvec)
      endif
    elseif(error.ne.0)then
       if(error.gt.0 .and. error.lt.3)then
         aerrors(error)=aerrors(error)+1
       else
         aerrors(3)=aerrors(3)+1
       endif
    endif
    if(Temper.eq.1.d0)then
      misfit = misfitstore(ichain)
      if(sVR.lt.(1-(misfit/sqUvec)))then
        sVR = 1-(misfit/sqUvec)
      endif
    endif
    if(mod(is,abins).eq.0 .and. ichain.eq.1)then
      tmp = dble(abins*nchainsG)
      write(afile,'(6I8,2F10.2,3F8.1,7F8.2)') is,aaccept,aerrors,aVR*100.d0,sVR*100.d0,&
         aNnuclei(1)/tmp,100.*aNnuclei(2:3)/tmp,100.*athr(1:7)/aNnuclei(1)
      flush(afile)
      aaccept=0
      aerrors=0
      athr=0
      aNnuclei=0
      aVR=-99.99d0
      sVR=-99.99d0
    endif
    
    !----------------------------------
    ! Write a time-slice from sampling chain
    !if(ichain.eq.1 .and. (rank.eq.1 .or. nproc.eq.1) )then
    !  if(is<10)then
    !    write(cproc,'(I1)') is
    !    filename=trim(dirG)//'mov/000000'//trim(cproc)//'.txt'
    !  elseif(is<100)then
    !    write(cproc,'(I2)') is
    !    filename=trim(dirG)//'mov/00000'//trim(cproc)//'.txt'
    !  elseif(is<1000)then
    !    write(cproc,'(I3)') is
    !    filename=trim(dirG)//'mov/0000'//trim(cproc)//'.txt'
    !  elseif(is<10000)then
    !    write(cproc,'(I4)') is
    !    filename=trim(dirG)//'mov/000'//trim(cproc)//'.txt'
    !  elseif(is<100000)then
    !    write(cproc,'(I5)') is
    !    filename=trim(dirG)//'mov/00'//trim(cproc)//'.txt'
    !  elseif(is<1000000)then
    !    write(cproc,'(I6)') is
    !    filename=trim(dirG)//'mov/0'//trim(cproc)//'.txt'
    !  else
    !    write(cproc,'(I7)') is
    !    filename=trim(dirG)//'mov/'//trim(cproc)//'.txt'
    !  endif
    !  ! Interpret model
    !  call InterpMod(ichain,error)
    !  ! Save model
    !  open(235,file=filename)
    !  write(235,'(I5)') m_Nlay(ichain)
    !  write(235,'(F10.3)') (1-(misfitstore(ichain)/sqUvec))*100.d0
    !  write(235,'(100000F8.1)') im_thick(1:m_Nlay(ichain))
    !  write(235,'(100000F8.1)') im_vs(1:m_Nlay(ichain))
    !  write(235,'(100000F8.1)') im_vp(1:m_Nlay(ichain))
    !  write(235,'(100000F8.1)') im_rho(1:m_Nlay(ichain))
    !  write(235,'(100000F8.4)') (((im_vp(i)**2) - 2.d0*(im_vs(i)**2))/&
    !         (2.d0*((im_vp(i)**2) - (im_vs(i)**2))),i=1,m_Nlay(ichain))
    !  close(235)
    !endif
    
    return
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE Setuptempladder(nchains,tlow,thigh,dir,Chaintemp)
!---------------------------------------------------------------------
!     Setuptempladder - User routine to set up temperatures for each chain
!     Input:
!        nchains    Integer   : Number of chains
!        Tlow       Double    : Lowest temperature
!        Thigh      Double    : Highest temperature
!        dir        Character(len=80)
!     Output:
!        Chaintemp  Double array (nchains)    : Temperatures of each chain
!                                               (on this processor)
!     Notes:
!           This utility routine calculates and returns the array Chaintemps.
!           It calculates a random log-uniform set of temperature values
!           between input bounds and put the results in Chaintemps.
!           Chaintemps is written out to file `tlevels'.
!---------------------------------------------------------------------
    use GlobalBox
    implicit none
#if defined MPI
    include "mpif.h"
    integer,dimension(MPI_STATUS_SIZE):: status
#endif
    integer,intent(in):: nchains
    real(8),intent(in):: tlow,thigh
    character(len=80),intent(in):: dir
    real(8):: Chaintemp(nchains)
#if defined MPI
    real(8),allocatable:: AllTemps(:)
    integer:: k,ierror
#endif
    real(8):: aval,bval
    integer:: it
    Character(len=255) filename
    
    !----------------------------------
    ! Selected temperatures randomly using log-uniform distribution
    aval = log(tlow)
    bval = log(thigh)
    do it=1,nchains
      Chaintemp(it) = exp(aval + ran3(iseed)*(bval-aval))
    enddo
    
    !----------------------------------
    ! Force first chain to be at tlow
    !if(rank==1)Chaintemp(1) = tlow       ! Force some chains to be at tlow
    Chaintemp(1) = tlow
    
    !----------------------------------
    ! Write to file
#if defined MPI
    allocate(AllTemps(nchains*nproc))
    AllTemps = 0.d0
    ! Send all Temperatures to master for output to file (for diagnostics)
    call MPI_GATHER(Chaintemp,nchains,MPI_DOUBLE_PRECISION,&
        AllTemps,nchains,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
    if(rank == 0)then
      filename = trim(dir)//'tlevels'
      open(15,file=filename,status='unknown')
      k = nchains*(nproc-1)
      if(k.gt.1)call HPSORT(k,AllTemps(nchains+1)) ! Order the temperatures for neat output
      do it = nchains+1,nproc*nchains
        write(15,*)it-nchains,AllTemps(it)
      enddo
      close(15)
    endif
    deallocate(AllTemps)
#else
    filename = trim(dir)//'log/tlevels.log'
    open(15,file=filename,status='unknown')
    do it=1,nchains
      write(15,*)it,Chaintemp(it)
    enddo
    close(15)
#endif
    RETURN
!---------------------------------------------------------------------
END SUBROUTINE



