! Library (module) of SH-transfer function
!
! Author: Miroslav Hallo (04/2020)
! ETH Zurich, Swiss Seismological Service
!
! This code is published under the GNU General Public License (GPL)
! for non-commercial purposes. To any licensee is given permission to 
! modify the work, as well as to copy and redistribute. Still
! we would like to kindly ask you to acknowledge the authors and don't
! remove their names from the code. This code is distributed in the
! hope that it will be useful, but WITHOUT ANY WARRANTY.
!---------------------------------------------------------------------

MODULE SHlib
!---------------------------------------------------------------------
!  Library with the seismic response of damped soil layers on rock.
!  Incident SH-wave, following pp. 268-270 in Kramer (1996)
!---------------------------------------------------------------------
CONTAINS



FUNCTION respSH(Nd,Nf,Thick,Vs,Rho,Fox,Freq,th)
!---------------------------------------------------------------------
!  Computing the seismic response of damped soil layers on rock
!  INPUT: Nd .. number of layers
!         Nf .. number of frequencies
!         Thick .. array of layer thicknesses [m]
!         Vs .. array of layer velocities [m/s]
!         Rho .. array of layer densities [kg/m^3]
!         Fox .. array of layer damping ratios
!         Freq .. array of frequencies to compute [Hz]
!         th .. wave incidence angle [deg] (0-60)
!  OUTPUT:  respSH(1:Nf,1) .. "surface to incident SH" transfer function
!           respSH(1:Nf,2) .. "surface to borehole" transfer function
!           respSH(1:Nf,3) .. "surface to rock outcrop" transfer function
!---------------------------------------------------------------------
    implicit none
    integer,intent(in):: Nd,Nf
    real(8),intent(in):: Thick(Nd), Vs(Nd), Rho(Nd), Fox(Nd), Freq(Nf), th
    complex(8):: respSH(Nf,3)
    complex(8):: VsS(Nd),Vslow(Nd),VsSimp(Nd),alphaS(Nd),A(Nd,Nf),B(Nd,Nf),ksH(Nf)
    real(8):: Hslow
    integer:: i,m
	
    !----------------------------------
    ! Prepare variables
    respSH=0.d0
    if(Nd.eq.1)then
      respSH(1:Nf,1) = 2.d0
      respSH(1:Nf,2) = 1.d0
      respSH(1:Nf,3) = 1.d0
      RETURN
    endif
    A=1.d0
    B=1.d0
	
    !----------------------------------
    ! Prepare the Complex shear wave velocity
    do i=1,Nd
      VsS(i) = complex(Vs(i),Vs(i)*Fox(i))
    enddo
    
    !----------------------------------
    ! Vertical inc. shear wave slowness
    Hslow = sin(th*0.01745329252d0)/Vs(Nd)
    do i=1,Nd
      Vslow(i) = sqrt((1.d0/(VsS(i)**2)) - (Hslow**2))
    enddo
	
	!----------------------------------
    ! Complex impedance ratio
    do i=1,Nd
      VsSimp(i) = (VsS(i)**2) * Vslow(i)
    enddo
    do i=1,Nd-1
      alphaS(i) = (Rho(i)*VsSimp(i)) / (Rho(i+1)*VsSimp(i+1))
    enddo
	
	!----------------------------------
    ! Recursive loop over layers
    do m=1,Nd-1
      ! Complex wave number x thickness
      do i=1,Nf
        ksH(i) = 6.2831853071796d0*Freq(i)*Thick(m)*Vslow(m)
      enddo
      ! Amplitudes of up-going and down-going waves
      do i=1,Nf
        A(m+1,i)=0.5d0*A(m,i)*(1.d0+alphaS(m))*exp(complex(0.d0,1.d0)*(ksH(i))) + &
                 0.5d0*B(m,i)*(1.d0-alphaS(m))*exp(complex(0.d0,-1.d0)*(ksH(i)))
        B(m+1,i)=0.5d0*A(m,i)*(1.d0-alphaS(m))*exp(complex(0.d0,1.d0)*(ksH(i))) + &
                 0.5d0*B(m,i)*(1.d0+alphaS(m))*exp(complex(0.d0,-1.d0)*(ksH(i)))
      enddo
    enddo
    
    !----------------------------------
    ! Transfer function surface to incident
    respSH(1:Nf,1) = 2.d0 / A(Nd,1:Nf)
	
    ! Transfer function surface to borehole
    respSH(1:Nf,2) = ( A(1,1:Nf) + B(1,1:Nf) ) / (  A(Nd,1:Nf) + B(Nd,1:Nf) )

    ! Transfer function surface to rock outcrop
    respSH(1:Nf,3) = 1.d0 / A(Nd,1:Nf)
    
    RETURN
!---------------------------------------------------------------------
END FUNCTION



!---------------------------------------------------------------------
END MODULE




