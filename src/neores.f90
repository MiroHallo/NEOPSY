! Statistics from the ensemble of solutions - MTI
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

MODULE ResBox
!---------------------------------------------------------------------
!  Module containing global variables
!---------------------------------------------------------------------
    implicit none
    ! Global model parameters
    integer,allocatable,dimension(:):: m_Nlay
    real(8),allocatable,dimension(:,:):: im_thick,im_vs,im_vp,im_rho
    
    ! Global pop bins
    integer:: nBins,nDepth,nXBins
    real(8),allocatable,dimension(:,:):: allBins
    real(8),allocatable,dimension(:):: depBins,logBins,vs30Bins
    real(8),allocatable,dimension(:):: layBins
    integer,allocatable,dimension(:,:):: pop_vs,pop_vp,pop_nu,pop_rho
    integer,allocatable,dimension(:):: pop_dep,pop_lay,pop_vr
    
    ! Global MAP search
    integer,allocatable,dimension(:,:):: tmpop_vs,tmpop_vp
    real(8),allocatable,dimension(:):: max_vs,max_vp,tm_vs,tm_vp
    
    ! Model count
    integer:: mcount,fcount
    real(8):: VRlim(2)
    
    ! Global QWL
    integer:: nFreq
    real(8):: FreqLim(2),DepLim(2),ImpLim(2)
    real(8),allocatable,dimension(:):: Freq,DepVec,ImpVec
    real(8),allocatable,dimension(:,:,:):: resQWL
    real(8),allocatable,dimension(:,:):: Vs30
    integer,allocatable,dimension(:,:):: pop_qwl_vs,pop_qwl_dep,pop_qwl_imp
    integer,allocatable,dimension(:):: pop_vs30
    
    ! Global datafit
    real(8),allocatable,dimension(:):: DW_vec
	
    ! Global SH amplification
    complex(8),allocatable,dimension(:,:,:):: resSH
    real(8),allocatable,dimension(:):: AmpVec,fox
    real(8):: th,AmpLim(2)
    integer,allocatable,dimension(:,:):: pop_SH_unr,pop_SH_ref
    
    ! Reference model for amplification
    character(len=255):: ref_file
    real(8),allocatable,dimension(:,:):: refmodel,refSH
    real(8),allocatable,dimension(:):: ampcor
    
!---------------------------------------------------------------------
END MODULE



PROGRAM  neores
!---------------------------------------------------------------------
!  Afterprocess of the Multizonal Transdimensional Inversion
!---------------------------------------------------------------------
    use SHlib
    use QWLib
    use LibBox
    use InputBox
    use DataBox
    use ResBox
    implicit none
    integer:: nproc,p,run,b,ifile,ofile,ios,im,ib1,ib2,ib3,i,popInd,popD1,popD2,popD1l,popD2l,nL
    character(len=80):: dir
    character(len=255):: filename
    character(len=128):: fileprefix
    character(len=4):: cproc,mfnam
    real(8):: logPPD,misfit,depTmp,nu
    real(8):: best_misb1,best_misb2,best_misb3,sqUvec,best_norm2,best_norm3
    logical:: file_exists
    integer:: is,error
    
    !----------------------------------
    ! Settings
    nBins = 100 ! number of bins for all statistics
    nXBins  = 1000 ! number of bins for vs30
    nFreq = 300 ! number of frequencies for QWL search
    
    VRlim(1:2) = (/ 0.d0, 100.d0 /) ! MIN and MAX of variance reduction bins
    FreqLim(1:2) = (/ 0.5d0, 30.d0 /) ! MIN and MAX freq. liniths for QWL bins (log sampling)
    DepLim(1:2) = (/ 1.0d0, 10000.d0 /) ! MIN and MAX depth liniths for QWL bins (log sampling)
    ImpLim(1:2) = (/ 0.0d0, 8.d0 /) ! MIN and MAX impedance for QWL bins (linear sampling)
    AmpLim(1:2) = (/ 0.5d0, 10.d0 /) ! MIN and MAX amplification liniths for SH bins (log sampling)
	
    mfnam = 'out' ! text identificator of the output files
    nproc = 1000 ! maximal number of CPU procesors (xmodels files)
	
    write(*,*) "------------------------------------------------------"
    write(*,*) "Multizonal Transdimensional Inversion (Fjord 3/2026)"
    write(*,*) "Miroslav Hallo, SED, ETH Zurich"
    write(*,*) "------------------------------------------------------"
    
    !----------------------------------
    ! Read input files
    write(*,*) '-> Read inputs'
    call InitInput()
    dir = dirG        ! work directory with models and for outputs
    ref_file = dirREF ! path to ascii file with the reference model
    
    !----------------------------------
    ! Read data
    call InitData()
    
    !----------------------------------
    ! Prepare data arrays
    sqUvec = N_used
    allocate(DW_vec(sum(f_n)))
    DW_vec = 0.d0
    
    !----------------------------------
    ! Initial arrays for after processing
    write(*,*) '-> Initial arrays'
    call InitFin()
    
    ! Init param
    im = 1 ! The test model index
    ib1 = 2  ! The maximum likelyhood model index
    ib2 = 3  ! The MAP S+P model index
    ib3 = 4  ! The MAP S model index
    best_misb1 = 1d10 ! initial (very high) misfit
    best_misb2 = 1d10 ! initial (very high) misfit
    best_misb3 = 1d10 ! initial (very high) misfit
    best_norm2 = 1d10 ! initial (very high) L1-norm misfit
    best_norm3 = 1d10 ! initial (very high) L1-norm misfit
    
    !----------------------------------
    ! Open file for synthetic data ensemble
    ofile = 93
    filename=trim(dir)//trim(mfnam)//'_pop_data.bin'
    open(unit=ofile,form='unformatted',access='stream',file=filename,action='write',status='replace')
    flush(ofile)
    
    ! Loop over two runs (POP and MAP)
    do run=1,2
    
      mcount = 0
      fcount = 0
      
      if(run.eq.1)then
        write(*,*) '-> POP loop'
      else
        write(*,*) '-> MAP loop'
      endif
      
      !----------------------------------
      ! Loop over all files with results
      do p=0,nproc-1
        ! File
        write(cproc,'(I0)') p
        filename=trim(dir)//'xmodels_'//trim(cproc)//'.dat'
        ! Check if the file exists
        inquire(file=filename, exist=file_exists)
        if(.NOT.file_exists)then
          cycle
        endif
        fcount = fcount+1
        
        ! Open the file
        ifile = 73
        open(unit=ifile,form='unformatted',access='stream',file=filename,action='read')
        write(*,'(A3,A)') '   ',trim(filename)
        
        !----------------------------------
        ! Loop over models within the file
        do while(.true.)
          ! Clean previous model
          im_thick(:,im)=0.d0
          im_vs(:,im)=0.d0
          im_vp(:,im)=0.d0
          im_rho(:,im)=0.d0
          
          ! Read new model (Eagle)
          read(ifile,iostat=ios) is
          if(ios.ne.0) exit
          read(ifile,iostat=ios) logPPD,misfit
          if(ios.ne.0) exit
          read(ifile,iostat=ios) m_Nlay(im)
          if(ios.ne.0) exit
          read(ifile,iostat=ios) im_thick(1:m_Nlay(im),im)
          if(ios.ne.0) exit
          read(ifile,iostat=ios) im_vs(1:m_Nlay(im),im)
          if(ios.ne.0) exit
          read(ifile,iostat=ios) im_vp(1:m_Nlay(im),im)
          if(ios.ne.0) exit
          read(ifile,iostat=ios) im_rho(1:m_Nlay(im),im)
          if(ios.ne.0) exit
          if(m_Nlay(im).eq.0) exit
          
          ! Skip model if
          if(isnan((1.d0-(misfit/sqUvec))*100.d0)) cycle
          if(is.le.iburnG) cycle
          
          mcount = mcount+1
        
          ! Screen output
          if(mod(mcount,1000).eq.0)then
            write(*,'(I11,A17)') mcount,' models processed'
          endif
          
          if(run.eq.1)then
            ! Count layers population
            popInd = binarysearch(nmax,layBins(:),dble(m_Nlay(im)))
            pop_lay(popInd) = pop_lay(popInd)+1
            
            ! Count VR
            popInd = binarysearch(nBins,allBins(:,5),(1.d0-(misfit/sqUvec))*100.d0)
            if(popInd.gt.0)then
              pop_vr(popInd) = pop_vr(popInd)+1
            endif
            
            ! Count population
            depTmp = 0.d0
            do i=1,m_Nlay(im)
              ! Depth index
              if(i.eq.1)then
                popD1 = 1
                popD1l = 1
              else
                popD1 = popD2 + 1
                popD1l = popD2l + 1
              endif
              if(i.eq.m_Nlay(im))then
                popD2 = nDepth
                popD2l = nDepth
              else
                depTmp = depTmp + im_thick(i,im)
                ! Depth of interface
                popD2 = binarysearch(nDepth,depBins(:),depTmp)
                popD2 = popD2 - 1
                ! Pop of the Log-Depth of interfaces
                popD2l = binarysearch(nDepth,logBins(:),depTmp)
                pop_dep(popD2l) = pop_dep(popD2l) + 1
                popD2l = popD2l - 1
              endif
              ! Vs
              popInd = binarysearch(nBins,allBins(:,1),im_vs(i,im))
              pop_vs(popInd,popD1:popD2) = pop_vs(popInd,popD1:popD2) + 1
              ! Vp
              popInd = binarysearch(nBins,allBins(:,2),im_vp(i,im))
              pop_vp(popInd,popD1:popD2) = pop_vp(popInd,popD1:popD2) + 1
              ! nu
              nu = ((im_vp(i,im)**2) - 2.d0*(im_vs(i,im)**2)) /&
                   (2.d0*((im_vp(i,im)**2) - (im_vs(i,im)**2)))
              if(nu.lt.0)then
                nu=0.d0
              elseif(nu.gt.0.5d0)then
                nu=0.5d0
              endif
              popInd = binarysearch(nBins,allBins(:,3),nu)
              pop_nu(popInd,popD1:popD2) = pop_nu(popInd,popD1:popD2) + 1
              ! Rho
              popInd = binarysearch(nBins,allBins(:,4),im_rho(i,im))
              pop_rho(popInd,popD1:popD2) = pop_rho(popInd,popD1:popD2) + 1
              ! Vs log for MAP
              popInd = binarysearch(nBins,allBins(:,1),im_vs(i,im))
              tmpop_vs(popInd,popD1l:popD2l) = tmpop_vs(popInd,popD1l:popD2l) + 1
              ! Vp log for MAP
              popInd = binarysearch(nBins,allBins(:,2),im_vp(i,im))
              tmpop_vp(popInd,popD1l:popD2l) = tmpop_vp(popInd,popD1l:popD2l) + 1
            enddo
        
            ! Quarter-wavelength representation
            nL = m_Nlay(im)
            resQWL(1:nFreq,1:4,im) = qwlf(nL,nFreq,im_thick(1:nL,im),im_vs(1:nL,im),im_rho(1:nL,im),Freq(1:nFreq))
            
            ! Vs30 and respective quarter-wavelength frequency
            Vs30(1:2,im) = qwl30(nL,im_thick(1:nL,im),im_vs(1:nL,im),30.d0)
            
            ! Count QWL population
            do i=1,nFreq
              ! QWL-Vs
              popInd = binarysearch(nBins,allBins(:,1),resQWL(i,2,im))
              if(popInd.gt.0)then
                pop_qwl_vs(popInd,i) = pop_qwl_vs(popInd,i) + 1
              endif
              ! QWL-Depth
              popInd = binarysearch(nBins,DepVec,resQWL(i,1,im))
              if(popInd.gt.0)then
                pop_qwl_dep(popInd,i) = pop_qwl_dep(popInd,i) + 1
              endif
              ! QWL-Impedance
              popInd = binarysearch(nBins,ImpVec,1.d0/resQWL(i,4,im))
              if(popInd.gt.0)then
                pop_qwl_imp(popInd,i) = pop_qwl_imp(popInd,i) + 1
              endif
            enddo
            
            ! Count Vs30 population
            popInd = binarysearch(nXBins,vs30Bins,Vs30(2,im))
            pop_vs30(popInd) = pop_vs30(popInd) + 1
            
            ! SH-transfer
        resSH(1:nFreq,1:3,im)=respSH(nL,nFreq,im_thick(1:nL,im),im_vs(1:nL,im),im_rho(1:nL,im),fox(1:nL),Freq(1:nFreq),th)
            ! correction to the reference model
            ampcor(1:nFreq) = qwla(nFreq,im_vs(m_Nlay(im),im),2500.d0,refmodel(1:nFreq,2),refmodel(1:nFreq,3))
            do i=1,nFreq
              refSH(i,im) = abs(resSH(i,3,im))/ampcor(i)
            enddo
            
            ! Count SH amplification
            do i=1,nFreq
              ! SH, elastic, outcrop, unreferenced
              popInd = binarysearch(nBins,AmpVec,abs(resSH(i,3,im)))
              if(popInd.gt.0)then
                pop_SH_unr(popInd,i) = pop_SH_unr(popInd,i) + 1
              endif
              ! SH, elastic, outcrop, referenced
              popInd = binarysearch(nBins,AmpVec,refSH(i,im))
              if(popInd.gt.0)then
                pop_SH_ref(popInd,i) = pop_SH_ref(popInd,i) + 1
              endif
            enddo
        
            ! Compute and save synthetic data ensemble
            if(ExpertSwitch(1).eq.1)then
              if(mod(mcount,3000).eq.0)then
                ! Forward problem
                call ForwardRes(im,error)
                if(error.eq.0)then
                  ! Save
                  write(ofile) D_vec(1:sum(f_n))
                  write(ofile) DW_vec(1:sum(f_n))
                  flush(ofile)
                endif
              endif
            else
              if(mod(mcount,300).eq.0)then
                ! Forward problem
                call ForwardRes(im,error)
                if(error.eq.0)then
                  ! Save
                  write(ofile) D_vec(1:sum(f_n))
                  write(ofile) DW_vec(1:sum(f_n))
                  flush(ofile)
                endif
              endif
            endif
        
            ! The best model
            if(misfit.lt.best_misb1)then
              ! Rewrite 1st
              m_Nlay(ib1) = m_Nlay(im)
              im_thick(:,ib1) = im_thick(:,im)
              im_vs(:,ib1) = im_vs(:,im)
              im_vp(:,ib1) = im_vp(:,im)
              im_rho(:,ib1) = im_rho(:,im)
              resQWL(:,:,ib1) = resQWL(:,:,im)
              Vs30(:,ib1) = Vs30(:,im)
              resSH(:,:,ib1) = resSH(:,:,im)
              refSH(:,ib1) = refSH(:,im)
              best_misb1 = misfit
            endif
            
          else
            ! Prepare model in log-depths
            depTmp = 0.d0
            tm_vs = 0.d0
            tm_vp = 0.d0
            do i=1,m_Nlay(im)
              ! Depth index
              if(i.eq.1)then
                popD1l = 1
              else
                popD1l = popD2l + 1
              endif
              if(i.eq.m_Nlay(im))then
                popD2l = nDepth
              else
                depTmp = depTmp + im_thick(i,im)
                ! Pop of the Log-Depth of interfaces
                popD2l = binarysearch(nDepth,logBins(:),depTmp)
                popD2l = popD2l - 1
              endif
              ! the model in log-depths
              tm_vs(popD1l:popD2l) = im_vs(i,im)
              tm_vp(popD1l:popD2l) = im_vp(i,im)
            enddo
            
            ! Compute L1-norm to the max of PDF
            do b=1,nDepth
              tm_vs(b) = abs(max_vs(b)-tm_vs(b))
              tm_vp(b) = abs(max_vp(b)-tm_vp(b))
            enddo
            
            ! Compare with the saved MAP-SP model
            if((sum(tm_vs(1:nDepth))+(sum(tm_vp(1:nDepth))/2.d0)).lt.best_norm2)then
              ! Compute Quarter-wavelength representation
              nL = m_Nlay(im)
              resQWL(1:nFreq,1:4,im) = qwlf(nL,nFreq,im_thick(1:nL,im),im_vs(1:nL,im),im_rho(1:nL,im),Freq(1:nFreq))
              ! Vs30 and respective quarter-wavelength frequency
              Vs30(1:2,im) = qwl30(nL,im_thick(1:nL,im),im_vs(1:nL,im),30.d0)
              ! SH-transfer
        resSH(1:nFreq,1:3,im)=respSH(nL,nFreq,im_thick(1:nL,im),im_vs(1:nL,im),im_rho(1:nL,im),fox(1:nL),Freq(1:nFreq),th)
              ! correction to the reference model
              ampcor(1:nFreq) = qwla(nFreq,im_vs(m_Nlay(im),im),2500.d0,refmodel(1:nFreq,2),refmodel(1:nFreq,3))
              do i=1,nFreq
                refSH(i,im) = abs(resSH(i,3,im))/ampcor(i)
              enddo
              ! Rewrite
              m_Nlay(ib2) = m_Nlay(im)
              im_thick(:,ib2) = im_thick(:,im)
              im_vs(:,ib2) = im_vs(:,im)
              im_vp(:,ib2) = im_vp(:,im)
              im_rho(:,ib2) = im_rho(:,im)
              resQWL(:,:,ib2) = resQWL(:,:,im)
              Vs30(:,ib2) = Vs30(:,im)
              resSH(:,:,ib2) = resSH(:,:,im)
              refSH(:,ib2) = refSH(:,im)
              best_misb2 = misfit
              best_norm2 = sum(tm_vs(1:nDepth)) + (sum(tm_vp(1:nDepth))/2.d0)
            endif
            
            ! Compare with the saved MAP-S model
            if(sum(tm_vs(1:nDepth)).lt.best_norm3)then
              ! Compute Quarter-wavelength representation
              nL = m_Nlay(im)
              resQWL(1:nFreq,1:4,im) = qwlf(nL,nFreq,im_thick(1:nL,im),im_vs(1:nL,im),im_rho(1:nL,im),Freq(1:nFreq))
              ! Vs30 and respective quarter-wavelength frequency
              Vs30(1:2,im) = qwl30(nL,im_thick(1:nL,im),im_vs(1:nL,im),30.d0)
              ! SH-transfer
        resSH(1:nFreq,1:3,im)=respSH(nL,nFreq,im_thick(1:nL,im),im_vs(1:nL,im),im_rho(1:nL,im),fox(1:nL),Freq(1:nFreq),th)
              ! correction to the reference model
              ampcor(1:nFreq) = qwla(nFreq,im_vs(m_Nlay(im),im),2500.d0,refmodel(1:nFreq,2),refmodel(1:nFreq,3))
              do i=1,nFreq
                refSH(i,im) = abs(resSH(i,3,im))/ampcor(i)
              enddo
              ! Rewrite
              m_Nlay(ib3) = m_Nlay(im)
              im_thick(:,ib3) = im_thick(:,im)
              im_vs(:,ib3) = im_vs(:,im)
              im_vp(:,ib3) = im_vp(:,im)
              im_rho(:,ib3) = im_rho(:,im)
              resQWL(:,:,ib3) = resQWL(:,:,im)
              Vs30(:,ib3) = Vs30(:,im)
              resSH(:,:,ib3) = resSH(:,:,im)
              refSH(:,ib3) = refSH(:,im)
              best_misb3 = misfit
              best_norm3 = sum(tm_vs(1:nDepth))
            endif
            
          endif
          
        enddo ! End loop over models within the file
      
        !----------------------------------
        ! Close the file
        close(ifile)
      
      enddo  ! End loop over all files
      
      if(run.eq.1)then
        !----------------------------------
        ! Close the pop_data file
        close(ofile)
        
        !----------------------------------
        ! Write stats of after processing
        fileprefix = trim(dir)//trim(mfnam)//'_pop'
        call SaveStats(fileprefix)
        
        !----------------------------------
        ! Write ML model
        fileprefix = trim(dir)//trim(mfnam)//'_modelML'
        call SaveMod(ib1,fileprefix,1.d0-(best_misb1/sqUvec))
      
        !----------------------------------
        ! Prepare max for the MAP run
        do b=1,nDepth
          max_vs(b) = allBins(maxloc(tmpop_vs(1:nBins,b),1),1) + ((allBins(2,1)-allBins(1,1))/2.d0)
          max_vp(b) = allBins(maxloc(tmpop_vp(1:nBins,b),1),2) + ((allBins(2,2)-allBins(1,2))/2.d0)
        enddo
        
      else
        !----------------------------------
        ! Write MAP-SP model
        fileprefix = trim(dir)//trim(mfnam)//'_modelMAP-SP'
        call SaveMod(ib2,fileprefix,1.d0-(best_misb2/sqUvec))
        
        !----------------------------------
        ! Write MAP-S model
        fileprefix = trim(dir)//trim(mfnam)//'_modelMAP-S'
        call SaveMod(ib3,fileprefix,1.d0-(best_misb3/sqUvec))
      endif
      
    enddo  ! End loop over two runs
    
    write(*,*) "------------------------------------------------------"
    write(*,*) "Inversion Summary"
    write(*,'(A,I14)') ' Number of sampling models in total:',mcount
    write(*,'(A,I14)') ' Number of visited models in total: ',nchainsG*mcount
    write(*,'(A,I14)') ' Number of burn-in steps per node:  ',iburnG
    write(*,'(A,I14)') ' Average number of production steps per node: ',int(mcount/max(fcount-1,1))
    write(*,'(A,3I5)') ' Number of nodes (total/master/slave):  ',fcount,min(fcount-1,1),max(fcount-1,1)
    write(*,'(A,F6.1,A1)') ' Data VR of ML model:      ',(1.d0-(best_misb1/sqUvec))*100.d0,'%'
    write(*,'(A,F6.1,A1)') ' Data VR of MAP-SP model:  ',(1.d0-(best_misb2/sqUvec))*100.d0,'%'
    write(*,'(A,F6.1,A1)') ' Data VR of MAP-S model:   ',(1.d0-(best_misb3/sqUvec))*100.d0,'%'
    write(*,*) "Done."
    
!---------------------------------------------------------------------
END PROGRAM



SUBROUTINE InitFin()
!---------------------------------------------------------------------
!  Initial after processing
!---------------------------------------------------------------------
    use QWLib
    use InputBox
    use ResBox
    implicit none
    integer:: nm,b
    real(8):: minD
    
    nm = 4    ! 1=test, 2=best1, 3=best2, 4=best3
    
    !----------------------------------
    ! Discretization of pop in depth (min step)
    if(thr_d(2).le.100.d0)then
      minD = 0.1d0
    elseif(thr_d(2).le.1000.d0)then
      minD = 1.d0
    elseif(thr_d(2).le.10000.d0)then
      minD = 10.d0
    else
      minD = 100.d0
    endif
    
    !----------------------------------
    ! Allocate model arrays
    allocate(m_Nlay(nm))
    m_Nlay=0
    
    !----------------------------------
    ! Allocate forward im arrays
    allocate(im_thick(nmax,nm))
    allocate(im_vs(nmax,nm))
    allocate(im_vp(nmax,nm))
    allocate(im_rho(nmax,nm))
    im_thick=0.d0
    im_vs=0.d0
    im_vp=0.d0
    im_rho=0.d0
    
    !----------------------------------
    ! Alocate arrays for pop
    nDepth = int(thr_d(2)/minD + 0.5d0)
    allocate(allBins(nBins,5))
    allocate(depBins(nDepth))
    allocate(logBins(nDepth))
    allocate(vs30Bins(nXBins))
    
    allocate(layBins(nmax))
    allocate(pop_vs(nBins,nDepth))
    allocate(pop_vp(nBins,nDepth))
    allocate(pop_nu(nBins,nDepth))
    allocate(pop_rho(nBins,nDepth))
    allocate(pop_dep(nDepth))
    allocate(pop_lay(nmax))
    allocate(pop_vr(nBins))
    pop_vs=0
    pop_vp=0
    pop_nu=0
    pop_rho=0
    pop_dep=0
    pop_lay=0
    pop_vr=0
    
    !----------------------------------
    ! Alocate arrays for MAP search
    allocate(tmpop_vs(nBins,nDepth))
    allocate(tmpop_vp(nBins,nDepth))
    allocate(max_vs(nDepth))
    allocate(max_vp(nDepth))
    allocate(tm_vs(nDepth))
    allocate(tm_vp(nDepth))
    tmpop_vs=0
    tmpop_vp=0
    max_vs=0.d0
    max_vp=0.d0
    tm_vs=0.d0
    tm_vp=0.d0
    
    !----------------------------------
    ! Alocate arrays for QWL
    allocate(Freq(nFreq))
    allocate(DepVec(nBins))
    allocate(ImpVec(nBins))
    allocate(resQWL(nFreq,4,nm))
    allocate(Vs30(2,nm))
    Freq(1:nFreq) = loguni(nFreq,FreqLim(1),FreqLim(2))
    DepVec(1:nBins) = loguni(nBins,DepLim(1),DepLim(2))
    do b=1,nBins
      ImpVec(b) = ImpLim(1) + dble(b-1)*(ImpLim(2)-ImpLim(1))/dble(nBins)
    enddo
    resQWL=0.d0
    Vs30=0.d0
    
    !----------------------------------
    ! Alocate arrays for QWL pop
    allocate(pop_qwl_vs(nBins,nFreq))
    allocate(pop_qwl_dep(nBins,nFreq))
    allocate(pop_qwl_imp(nBins,nFreq))
    allocate(pop_vs30(nXBins))
    pop_qwl_vs=0
    pop_qwl_dep=0
    pop_qwl_imp=0
    pop_vs30=0
    
    !----------------------------------
    ! Alocate arrays for SH amplification
    allocate(resSH(nFreq,3,nm))
    allocate(fox(nmax))
    allocate(AmpVec(nBins))
    allocate(pop_SH_unr(nBins,nFreq))
    allocate(pop_SH_ref(nBins,nFreq))
    AmpVec(1:nBins) = loguni(nBins,AmpLim(1),AmpLim(2))
    resSH=0.d0
    fox=0.d0
    th=0.d0
    pop_SH_unr=0
    pop_SH_ref=0
    
	!----------------------------------
    ! Swiss reference amplification
    allocate(refmodel(nFreq,4))
    allocate(ampcor(nFreq))
    allocate(refSH(nFreq,nm))
    refmodel(1:nFreq,1:4) = qwlref(nFreq,Freq(1:nFreq),ref_file)
    ampcor=1.d0
    refSH=0.d0
    
    !----------------------------------
	! Prepare bins for parameters
    allBins=0.d0
    do b=1,nBins
      allBins(b,1) = minval(thr_vs(1,:)) + dble(b-1)*(maxval(thr_vs(2,:))-minval(thr_vs(1,:)))/dble(nBins) ! Vs
      allBins(b,2) = minval(thr_vp(1,:)) + dble(b-1)*(maxval(thr_vp(2,:))-minval(thr_vp(1,:)))/dble(nBins) ! Vp
      allBins(b,3) = minval(thr_nu(1,:)) + dble(b-1)*(maxval(thr_nu(2,:))-minval(thr_nu(1,:)))/dble(nBins) ! nu
      allBins(b,4) = minval(thr_rho(1,:)) + dble(b-1)*(maxval(thr_rho(2,:))-minval(thr_rho(1,:)))/dble(nBins) ! Rho
      allBins(b,5) = VRlim(1) + dble(b-1)*(VRlim(2)-VRlim(1))/dble(nBins) ! VR
    enddo
    depBins=0.d0
    do b=1,nDepth
      depBins(b) = dble(b-1)*minD
    enddo
    logBins=0.d0
    do b=2,nDepth
      logBins(b) = exp(log(thr_d(1)) + dble(b-1)*(log(thr_d(2))-log(thr_d(1)))/dble(nDepth-1))
    enddo
    layBins=0.d0
    do b=1,nmax
      layBins(b) = dble(b)
    enddo
    vs30Bins=0.d0
    do b=1,nXBins
      vs30Bins(b) = minval(thr_vs(1,:)) + dble(b-1)*(maxval(thr_vs(2,:))-minval(thr_vs(1,:)))/dble(nXBins)
    enddo
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE SaveStats(fileprefix)
!---------------------------------------------------------------------
!  Save pop statistics into .txt files
!---------------------------------------------------------------------
    use ResBox
    implicit none
    character(len=128):: fileprefix
    character(len=255):: filename
    integer:: b
    
    !  Headers of files with results
    filename=trim(fileprefix)//'_headers.txt'
    open(230,file=filename)
    write(230,'(A)') '#nBins / nDepth / models / dummy'
    write(230,'(I5)') nBins
    write(230,'(I5)') nDepth
    write(230,'(I11)') mcount
    write(230,'(F6.3)') 0.666d0
    write(230,'(1000000F8.1)') allBins(1:nBins,1) ! Vs
    write(230,'(1000000F8.1)') allBins(1:nBins,2) ! Vp
    write(230,'(1000000F8.4)') allBins(1:nBins,3) ! nu
    write(230,'(1000000F8.1)') allBins(1:nBins,4) ! Rho
    write(230,'(1000000F8.4)') allBins(1:nBins,5) ! VR
    write(230,'(1000000F9.1)') depBins(1:nDepth) ! Depth
    write(230,'(1000000E15.6)') logBins(1:nDepth) ! Log-Depth
    close(230)
    
    ! Vs
    filename=trim(fileprefix)//'_vs2D.txt'
    open(231,file=filename)
    do b=1,nDepth
      write(231,'(1000000I11)') pop_vs(1:nBins,b)
    enddo
    close(231)
    
    ! Vs
    filename=trim(fileprefix)//'_vp2D.txt'
    open(232,file=filename)
    do b=1,nDepth
      write(232,'(1000000I11)') pop_vp(1:nBins,b)
    enddo
    close(232)
    
    ! nu
    filename=trim(fileprefix)//'_nu2D.txt'
    open(233,file=filename)
    do b=1,nDepth
      write(233,'(1000000I11)') pop_nu(1:nBins,b)
    enddo
    close(233)
    
    ! Rho
    filename=trim(fileprefix)//'_rho2D.txt'
    open(234,file=filename)
    do b=1,nDepth
      write(234,'(1000000I11)') pop_rho(1:nBins,b)
    enddo
    close(234)
    
    ! Depth
    filename=trim(fileprefix)//'_dep1D.txt'
    open(235,file=filename)
    do b=1,nDepth
      write(235,'(1000000I11)') pop_dep(b)
    enddo
    close(235)
    
    ! Number of layeres
    filename=trim(fileprefix)//'_lay1D.txt'
    open(236,file=filename)
    write(236,'(1000000I11)') pop_lay(:)
    close(236)
    
    ! Number of layeres
    filename=trim(fileprefix)//'_vr1D.txt'
    open(237,file=filename)
    write(237,'(1000000I11)') pop_vr(:)
    close(237)
    
    !  QWL headers
    filename=trim(fileprefix)//'_QWL_head.txt'
    open(330,file=filename)
    write(330,'(A)') '#nBins / nFreq'
    write(330,'(I5)') nBins
    write(330,'(I5)') nFreq
    write(330,'(I11)') mcount
    write(330,'(1000000E15.6)') Freq(1:nFreq) ! Frequencies
    write(330,'(1000000F8.1)') allBins(1:nBins,1) ! QWL Vs
    write(330,'(1000000E15.6)') DepVec(1:nBins) ! QWL depths
    write(330,'(1000000E15.6)') ImpVec(1:nBins) ! QWL impedance
    write(330,'(1000000E15.6)') AmpVec(1:nBins) ! SH amplification
    write(330,'(I5)') nXBins
    write(330,'(1000000E15.6)') vs30Bins(1:nXBins) ! Vs30 bins
    close(330)
    
    ! QWL - Vs
    filename=trim(fileprefix)//'_QWL_vs2D.txt'
    open(331,file=filename)
    do b=1,nFreq
      write(331,'(1000000I11)') pop_qwl_vs(1:nBins,b)
    enddo
    close(331)
    
    ! QWL - Depth
    filename=trim(fileprefix)//'_QWL_dep2D.txt'
    open(332,file=filename)
    do b=1,nFreq
      write(332,'(1000000I11)') pop_qwl_dep(1:nBins,b)
    enddo
    close(332)
    
    ! QWL - Impedance
    filename=trim(fileprefix)//'_QWL_imp2D.txt'
    open(333,file=filename)
    do b=1,nFreq
      write(333,'(1000000I11)') pop_qwl_imp(1:nBins,b)
    enddo
    close(333)
    
    ! Vs30
    filename=trim(fileprefix)//'_vs30.txt'
    open(334,file=filename)
    write(334,'(1000000I11)') pop_vs30(:)
    close(334)
    
	! SH outcrop amplification unreferenced
    filename=trim(fileprefix)//'_SH_unr2D.txt'
    open(335,file=filename)
    do b=1,nFreq
      write(335,'(1000000I11)') pop_SH_unr(1:nBins,b)
    enddo
    close(335)
    
    ! SH outcrop amplification referenced
    filename=trim(fileprefix)//'_SH_ref2D.txt'
    open(336,file=filename)
    do b=1,nFreq
      write(336,'(1000000I11)') pop_SH_ref(1:nBins,b)
    enddo
    close(336)
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE SaveMod(im,fileprefix,VR)
!---------------------------------------------------------------------
!  Save a model into .txt file
!---------------------------------------------------------------------
    use InputBox
    use DataBox
    use ResBox
    implicit none
    integer:: im
    real(8):: VR
    character(len=128):: fileprefix
    character(len=255):: filename
    integer:: i,m,ci
    real(8),allocatable,dimension(:):: nu
    integer:: error
    
    ! Allocate and compute nu
    allocate(nu(m_Nlay(im)))
    do i=1,m_Nlay(im)
      nu(i)=((im_vp(i,im)**2) - 2.d0*(im_vs(i,im)**2)) /&
            (2.d0*((im_vp(i,im)**2) - (im_vs(i,im)**2)))
    enddo
    
    !----------------------------------
    ! Write model
    filename=trim(fileprefix)//'.txt'
    open(235,file=filename)
    write(235,'(I5)') m_Nlay(im)
    write(235,'(F10.3)') VR*100.d0
    write(235,'(100000F8.1)') im_thick(1:m_Nlay(im),im)
    write(235,'(100000F8.1)') im_vs(1:m_Nlay(im),im)
    write(235,'(100000F8.1)') im_vp(1:m_Nlay(im),im)
    write(235,'(100000F8.1)') im_rho(1:m_Nlay(im),im)
    write(235,'(100000F8.4)') nu(1:m_Nlay(im))
    close(235)
    
    deallocate(nu)
    
    !----------------------------------
    ! Write QWL of the model
    filename=trim(fileprefix)//'_QWL.txt'
    open(390,file=filename)
    write(390,'(I5)') nFreq
    write(390,'(2E15.6)') Vs30(1:2,im) ! Vs30
    write(390,'(1000000E15.6)') Freq(1:nFreq) ! Frequencies
    write(390,'(1000000E15.6)') resQWL(1:nFreq,1,im)  ! QWL - Depths
    write(390,'(1000000E15.6)') resQWL(1:nFreq,2,im)  ! QWL - Vs
    write(390,'(1000000E15.6)') resQWL(1:nFreq,3,im)  ! QWL - Densities
    write(390,'(1000000E15.6)') 1.d0/resQWL(1:nFreq,4,im)  ! QWL - Impedance
    write(390,'(1000000E15.6)') abs(resSH(1:nFreq,3,im))  ! SH amplification unreferenced
    write(390,'(1000000E15.6)') refSH(1:nFreq,im)  ! SH amplification referenced
    close(390)
    
    !----------------------------------
    ! Data fit (forward problem)
    call ForwardRes(im,error)
    if(error.ne.0)then
      if(ExpertSwitch(1).eq.1)then
        write(*,*) 'WARNING 401: A problem in GEOPSY!'
      else
        write(*,*) 'ERROR 401: A problem in GEOPSY!'
        stop
      endif
    else
      !----------------------------------
      ! Write model data
      filename=trim(fileprefix)//'_data.txt'
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
          write(190,'(4E17.7)') f_vec(ci+i),D_vec(ci+i),W_vec(ci+i),DW_vec(ci+i)
        enddo
      enddo
      close(190)
    endif
    
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE ForwardRes(im,error)
!---------------------------------------------------------------------
!  Calculate forward problem (Geopsy 3.X)
!  In:  im - model identifer (integer)
!  Out: error iostat (=0 if all OK)
!---------------------------------------------------------------------
    use InputBox
    use DataBox
    use LibBox
    use ResBox
    use, NON_INTRINSIC :: IO
    implicit none
    integer,intent(in):: im
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
    write(cproc,'(I0)') im-1
    path=trim(dir)//'best'//trim(cproc)
    file_NP = trim(path)//'.npin'
    file_DC = trim(path)//'.npdc'
    file_GP = trim(path)//'.model'
    file_ELL = trim(path)//'.ell'
    
    !----------------------------------
    ! Compute Rayleigh and Love wave dispersion curves (Eagle)
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
        write(31) im_thick(1:nL,im) ! real(8)
        write(31) im_vs(1:nL,im) ! real(8)
        write(31) im_vp(1:nL,im) ! real(8)
        write(31) im_rho(1:nL,im) ! real(8)
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
      call INV2GEOPSY(file_GP,im_thick(1:nL,im),im_vp(1:nL,im),im_vs(1:nL,im),im_rho(1:nL,im),error)
      if(error.ne.0)then
        comm = 'rm '//trim(file_GP)//' > /dev/null 2>&1'
        call EXECUTE_COMMAND_LINE(comm, wait=.TRUE.)
        return
      endif
      ! Data oversampling
      if((log10(f_lim(2,m))-log10(f_lim(1,m))).lt.0.2d0)then
        dL = max(100,f_n(m))
      elseif((log10(f_lim(2,m))-log10(f_lim(1,m))).lt.2.d0)then
        dL = max(200,f_n(m))
      else
        dL = max(300,f_n(m))
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
    ! Prepara standardized synthetic data
    DW_vec = D_vec*W_vec
    
!---------------------------------------------------------------------
END SUBROUTINE




