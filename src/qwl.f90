! Library (module) of quarter-wavelength representation
!
! Author: Miroslav Hallo (11/2019)
! ETH Zurich, Swiss Seismological Service
!
! This code is published under the GNU General Public License (GPL)
! for non-commercial purposes. To any licensee is given permission to 
! modify the work, as well as to copy and redistribute. Still
! we would like to kindly ask you to acknowledge the authors and don't
! remove their names from the code. This code is distributed in the
! hope that it will be useful, but WITHOUT ANY WARRANTY.
!---------------------------------------------------------------------

MODULE QWLib
!---------------------------------------------------------------------
!  QWL library of functions of the quarter-wavelength representation
!  loguni(), qwlf(), qwl30(), qwla() qwlref()
!---------------------------------------------------------------------
CONTAINS



FUNCTION qwlf(Nd,Nf,Thick,Vs,Rho,Freq)
!---------------------------------------------------------------------
!  Quarter-wavelength depth, velocity, and density from velocity profile
!  INPUT: Nd .. number of layers
!         Nf .. number of frequencies
!         Thick .. array of layer thicknesses [m]
!         Vs .. array of layer velocities [m/s]
!         Rho .. array of layer densities [kg/m^3]
!         Freq .. array of frequencies to compute [Hz]
!  OUTPUT:  qwlf(1:Nf,1) .. quarter-wavelength depths [m]
!           qwlf(1:Nf,2) .. quarter-wavelength velocities [m/s]
!           qwlf(1:Nf,3) .. quarter-wavelength densities [kg/m^3]
!           qwlf(1:Nf,4) .. quarter-wavelength Impedance contrast
!---------------------------------------------------------------------
    implicit none
    integer,intent(in):: Nd,Nf
    real(8),intent(in):: Thick(Nd), Vs(Nd), Rho(Nd), Freq(Nf)
    real(8):: qwlf(Nf,4)
    real(8):: Depth(Nd),T0(Nd),RhoZ0(Nd)
    integer,parameter:: Nz = 2048 ! Number of discrete depths
    !integer,parameter:: Nz = 512 ! Number of discrete depths
    real(8):: z(Nz),T(Nz),VsQWL(Nz),VsDiff(Nz)
    integer:: IZ(Nz)
    integer:: i,ni,f,IZ_tmp,minI,IZ_f,minC
    real(8):: T_tmp,z_tmp,VsQWL_tmp
    
    !----------------------------------
    ! Prepare the vector of depths z (grid search)
    z(1) = 0.d0
    z(2:Nz) = loguni(Nz-1,0.1d0,10000.d0)
    
    !----------------------------------
    ! Prepare vector of layer depths
    qwlf = 0.d0
    Depth(1) = 0.d0
    do i=2,Nd
      Depth(i) = Depth(i-1) + Thick(i-1)
    enddo
    
    !----------------------------------
    ! Find layer indexes
    IZ(1:Nz) = Nd
    IZ_tmp = 1
    ni = 0
    do while(.true.)
      ni = ni+1
      if(IZ_tmp.eq.Nd .or. ni.gt.Nz)then
        exit
      elseif(Depth(IZ_tmp+1).ge.z(ni))then
        IZ(ni) = IZ_tmp
      else
        IZ_tmp = IZ_tmp + 1
        ni = ni-1
      endif
    enddo
    
    !----------------------------------
    ! Compute vertical travel-times
    T0(1) = 0.d0
    do i=2,Nd
      T0(i) = T0(i-1) + (Thick(i-1) / Vs(i-1))
    enddo
    do i=1,Nz
      T(i) = T0(IZ(i)) + (abs(z(i) - Depth(IZ(i))) / Vs(IZ(i)))
    enddo
    
    !----------------------------------
    ! Compute reference densities
    RhoZ0(1) = 0.d0
    do i=2,Nd
      RhoZ0(i) = RhoZ0(i-1) + (Rho(i-1) * Thick(i-1))
    enddo
    
    !----------------------------------
    ! Find the quarter-wavelength depth, velocity and density
    VsQWL(1) = 999999.d0
    do i=2,Nz
      VsQWL(i) = z(i) / T(i)
    enddo
    
    do f=1,Nf
      ! Find depth with minimal difference
      do i=1,Nz
        VsDiff(i) = abs( z(i) - (VsQWL(i)/(4.d0*Freq(f))) )
      enddo
      if((1.d0/(4.d0*Freq(f))).ge.T(Nz))then
        minI = Nz
      else
        minI = minloc(VsDiff,1)
      endif
      ! Assign quarter-wavelength depth and velocity
      qwlf(f,1) = z(minI)
      qwlf(f,2) = VsQWL(minI)
      ! Find quarter-wavelength density
      IZ_f = IZ(minI)
      qwlf(f,3) = (RhoZ0(IZ_f) + abs(z(minI)-Depth(IZ_f)) * Rho(IZ_f)) / z(minI)
      
      !----------------------------------
      ! Find quarter-wavelength Impedance contrast
      do i=minI,Nz
        T_tmp = T(i) - T(minI)
        z_tmp = z(i) - z(minI)
        VsQWL_tmp = z_tmp/T_tmp
        VsDiff(i) = abs( z_tmp - (VsQWL_tmp/(4.d0*Freq(f))) )
      enddo
      if((1.d0/(4.d0*Freq(f))).ge.T_tmp)then
        minC = Nz
      else
        minC = (minI-1) + minloc(VsDiff(minI:Nz),1)
      endif
      ! Assign quarter-wavelength Impedance values
      VsQWL_tmp = (z(minC) - z(minI)) / (T(minC) - T(minI))
      qwlf(f,4) = VsQWL(minI) / VsQWL_tmp
    enddo
    
    RETURN
!---------------------------------------------------------------------
END FUNCTION



FUNCTION qwl30(Nd,Thick,Vs,d30)
!---------------------------------------------------------------------
!  Reference (d30) quarter-wavelength freqency and velocity
!  INPUT: Nd .. number of layers
!         Thick .. array of layer thicknesses [m]
!         Vs .. array of layer velocities [m/s]
!         d30 .. reference depth (standard is 30) [m]
!  OUTPUT:  qwl30(1) .. reference (d30) quarter-wavelength freqency [Hz]
!           qwl30(2) .. reference (d30) velocity [m/s]
!---------------------------------------------------------------------
    implicit none
    integer,intent(in):: Nd
    real(8),intent(in):: Thick(Nd), Vs(Nd), d30
    real(8):: qwl30(2)
    real(8):: Depth(Nd),T0(Nd),T
    integer:: IZ,IZ_tmp,i,ni
    
    !----------------------------------
    ! Prepare vector of layer depths
    qwl30 = 0.d0
    Depth(1) = 0.d0
    do i=2,Nd
      Depth(i) = Depth(i-1) + Thick(i-1)
    enddo
    
    !----------------------------------
    ! Find the layer index
    IZ = Nd
    IZ_tmp = 1
    ni = 0
    do while(.true.)
      ni = ni+1
      if(IZ_tmp.eq.Nd .or. ni.gt.1)then
        exit
      elseif(Depth(IZ_tmp+1).ge.d30)then
        IZ = IZ_tmp
      else
        IZ_tmp = IZ_tmp + 1
        ni = ni-1
      endif
    enddo
    
    !----------------------------------
    ! Compute vertical travel-times
    T0(1) = 0.d0
    do i=2,IZ
      T0(i) = T0(i-1) + (Thick(i-1) / Vs(i-1))
    enddo
    T = T0(IZ) + (abs(d30 - Depth(IZ)) / Vs(IZ))
    
    !----------------------------------
    ! Find the vs30 frequency and velocity
    qwl30(1) = 1.d0/(4.d0*T)
    qwl30(2) = d30 / T
    
    RETURN
!---------------------------------------------------------------------
END FUNCTION



FUNCTION qwla(Nf,VsC,RhoC,Vs_f,Rho_f)
!---------------------------------------------------------------------
!  Quarter-wavelength amplification factor
!  INPUT: Nf .. number of frequencies
!         VsC .. Reference wave velocity [m/s]
!         RhoC .. Reference density [kg/m^3]
!         Vs_f .. array of QWL velocities [m/s]
!         Rho_f .. array of QWL densities [kg/m^3]
!  OUTPUT:  qwla(1:Nf) .. The amplification factor
!---------------------------------------------------------------------
    implicit none
    integer,intent(in):: Nf
    real(8),intent(in):: VsC, RhoC, Vs_f(Nf), Rho_f(Nf)
    real(8):: qwla(Nf)
    integer:: f
    
    !----------------------------------
    ! Find the amplification factor
    do f=1,Nf
      qwla(f) = sqrt( (VsC*RhoC)/(Vs_f(f)*Rho_f(f))  )
    enddo
    
    RETURN
!---------------------------------------------------------------------
END FUNCTION



FUNCTION qwlref(Nf,Freq,ref_file)
!---------------------------------------------------------------------
!  Quarter-wavelength amplification factor in the reference ascii model
!  INPUT: Nf .. number of frequencies
!         Freq .. array of frequencies to compute [Hz]
!         ref_file .. filepath of an ascii file with the reference model
!  OUTPUT:  qwlref(1:Nf,4) .. The referenced d_f, Vs_f, Rho_f, IC_f
!---------------------------------------------------------------------
    implicit none
    integer,intent(in):: Nf
    real(8),intent(in):: Freq(Nf)
    character(len=255),intent(in):: ref_file
    real(8):: qwlref(Nf,4)
    real(8),allocatable,dimension(:,:):: ref_model
    integer:: maxN,fcount,ios,nL
	
	!----------------------------------
    !  Allocate local variables
    maxN = 10000   ! Maximal number of layers
    allocate(ref_model(maxN,3))
    ref_model = 2500.d0
	
	!----------------------------------
    ! Read from ascii input file
    open(170,form='formatted',file=ref_file,action='read',status='old')
    fcount = 1
    do while(.true.) ! Loop over lines within the file
      read(170,*,iostat=ios) ref_model(fcount,1),ref_model(fcount,2)
      if(ios.ne.0) exit
      fcount = fcount + 1
    enddo
    close(170)
    nL = fcount - 1
	
    !----------------------------------
    ! Compute QWL of the reference model
    qwlref(1:Nf,1:4) = qwlf(nL,Nf,ref_model(1:nL,1),ref_model(1:nL,2),ref_model(1:nL,3),Freq(1:Nf))
	
	!----------------------------------
    ! Clear
    deallocate(ref_model)
	
    RETURN
!---------------------------------------------------------------------
END FUNCTION



FUNCTION loguni(Nf,fmin,fmax)
!---------------------------------------------------------------------
!  Log-uniform sampling (Nf samples between [fmin,fmax] )
!---------------------------------------------------------------------
    implicit none
    integer,intent(in):: Nf
    real(8),intent(in):: fmin,fmax
    real(8):: loguni(Nf)
    real(8):: lfmin,lfmax,lfdif
    integer:: i
    
    !----------------------------------
    ! just to be sure
    loguni = 0.d0
    if(Nf.lt.2)then
      write(*,'(A)') 'ERROR1 in loguni'
      RETURN
    endif
    if(fmin.ge.fmax)then
      write(*,'(A)') 'ERROR2 in loguni'
      RETURN
    endif
    
    !----------------------------------
    ! Log-uniform sampling
    lfmin = dlog10(fmin)
    lfmax = dlog10(fmax)
    lfdif = (lfmax - lfmin)/dble(Nf-1)
    do i=1,Nf
      loguni(i) = 10.d0 ** (lfmin + dble(i-1)*lfdif)
    enddo
    
    RETURN
!---------------------------------------------------------------------
END FUNCTION



!---------------------------------------------------------------------
END MODULE




