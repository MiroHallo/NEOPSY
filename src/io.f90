MODULE IO
  ! Architect: Walter Imperatori, 2019
  ! SET OF SUBROUTINES TO READ/WRITE OUTPUT/INPUT FILES FOR GEOPSY UTILITIES "GPDC", "GPELL"
  ! IMPORTANT NOTE: "GEOPSY2INV" MUST BE PROBABLY MODIFIED IF MORE THAN ONE MODE IS SOUGHT.

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL32, REAL64

  IMPLICIT NONE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! ALL SUBROUTINES ARE PUBLIC
  PUBLIC

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTEGER, PARAMETER, PRIVATE :: RSP = REAL32              !< REAL SINGLE PRECISION
  INTEGER, PARAMETER, PRIVATE :: RDP = REAL64              !< REAL DOUBLE PRECISION

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! OVERLOADING FOR SINGLE/DOUBLE PRECISION
  INTERFACE INV2GEOPSY
    MODULE PROCEDURE INV2GEOPSY_SP, INV2GEOPSY_DP
  END INTERFACE INV2GEOPSY

  INTERFACE GEOPSY2INV
    MODULE PROCEDURE GEOPSY2INV_SP, GEOPSY2INV_DP
  END INTERFACE GEOPSY2INV

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE INV2GEOPSY_SP(FNAME, THICKNESS, VP, VS, RHO, IERR)

      ! SUBPROGRAM TO PREPARE INPUT FILE FOR GEOPSY UTILITIES GPDC AND GPELL. I/O OPERATION CHECK ONLY WHEN OPENING/CLOSING FILE.

      CHARACTER(LEN=*),               INTENT(IN)  :: FNAME                               !< FILE NAME
      REAL(RSP),        DIMENSION(:), INTENT(IN)  :: THICKNESS                           !< THICKNESS (M)
      REAL(RSP),        DIMENSION(:), INTENT(IN)  :: VP                                  !< VP (M/S)
      REAL(RSP),        DIMENSION(:), INTENT(IN)  :: VS                                  !< VS (M/S)
      REAL(RSP),        DIMENSION(:), INTENT(IN)  :: RHO                                 !< DENSITY (KG/M3)
      INTEGER,                        INTENT(OUT) :: IERR                                !< ERROR FLAG
      INTEGER                                     :: I, N                                !< COUNTERS

      !-----------------------------------------------------------------------------------------------------------------------------

      ! GET NUMBER OF LAYERS
      N = SIZE(THICKNESS)

      OPEN(22, FILE = FNAME, STATUS = 'REPLACE', FORM = 'FORMATTED', IOSTAT = IERR)

      IF (IERR .NE. 0) RETURN

      WRITE(22, '(A)') '# First line: number of layers'
      WRITE(22, '(I0)') N

      WRITE(22, '(A)') '# Thickness (m), Vp (m/s), Vs (m/s) and density(kg/m3)'

      DO I = 1, N - 1
        WRITE(22, '(4F11.4)') THICKNESS(I), VP(I), VS(I), RHO(I)
      ENDDO

      WRITE(22, '(A)') '# Last line is the half-space'
      WRITE(22, '(4F11.4)') 0, VP(N), VS(N), RHO(N)

      CLOSE(22, IOSTAT = IERR)

    END SUBROUTINE INV2GEOPSY_SP

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE INV2GEOPSY_DP(FNAME, THICKNESS, VP, VS, RHO, IERR)

      ! SUBPROGRAM TO PREPARE INPUT FILE FOR GEOPSY UTILITIES GPDC AND GPELL. I/O OPERATION CHECK ONLY WHEN OPENING/CLOSING FILE.

      CHARACTER(LEN=*),               INTENT(IN)  :: FNAME                               !< FILE NAME
      REAL(RDP),        DIMENSION(:), INTENT(IN)  :: THICKNESS                           !< THICKNESS (M)
      REAL(RDP),        DIMENSION(:), INTENT(IN)  :: VP                                  !< VP (M/S)
      REAL(RDP),        DIMENSION(:), INTENT(IN)  :: VS                                  !< VS (M/S)
      REAL(RDP),        DIMENSION(:), INTENT(IN)  :: RHO                                 !< DENSITY (KG/M3)
      INTEGER,                        INTENT(OUT) :: IERR                                !< ERROR FLAG
      INTEGER                                     :: I, N                                !< COUNTERS

      !-----------------------------------------------------------------------------------------------------------------------------

      ! GET NUMBER OF LAYERS
      N = SIZE(THICKNESS)

      OPEN(22, FILE = FNAME, STATUS = 'REPLACE', FORM = 'FORMATTED', IOSTAT = IERR)

      IF (IERR .NE. 0) RETURN

      WRITE(22, '(A)') '# First line: number of layers'
      WRITE(22, '(I0)') N

      WRITE(22, '(A)') '# Thickness (m), Vp (m/s), Vs (m/s) and density(kg/m3)'

      DO I = 1, N - 1
        WRITE(22, '(4F11.4)') THICKNESS(I), VP(I), VS(I), RHO(I)
      ENDDO

      WRITE(22, '(A)') '# Last line is the half-space'
      WRITE(22, '(4F11.4)') 0., VP(N), VS(N), RHO(N)

      CLOSE(22, IOSTAT = IERR)

    END SUBROUTINE INV2GEOPSY_DP

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE GEOPSY2INV_SP(FNAME, FREQ, X, IERR)

      ! SUBPROGRAM TO READ OUTPUT FILE FROM GEOPSY UTILITIES GPDC AND GPELL. I/O OPERATION CHECK ONLY WHEN OPENING/CLOSING FILE.

      CHARACTER(LEN=*),               INTENT(IN)    :: FNAME                               !< FILE NAME
      REAL(RSP),        DIMENSION(:), INTENT(INOUT) :: FREQ                                !< FREQUENCY
      REAL(RSP),        DIMENSION(:), INTENT(INOUT) :: X                                   !< EITHER DISPERSION CURVE OR ELLIPTICITY
      INTEGER,                        INTENT(OUT)   :: IERR                                !< ERROR FLAG
      INTEGER                                       :: OK, I, L, N                         !< COUNTERS

      !-----------------------------------------------------------------------------------------------------------------------------

      ! NO ERROR BY DEFAUL
      IERR = 0

      ! GET NUMBER OF EXPECTED FREQUENCY POINTS
      L = SIZE(FREQ)

      OPEN(23, FILE = FNAME, STATUS = 'OLD', FORM = 'FORMATTED', IOSTAT = IERR)

      IF (IERR .NE. 0) RETURN

      OK = 0
      N  = 0

      ! CRUDE CHECK TO SEE IF GEOPSY CALCULATED SOMETHING...
      DO WHILE (OK .EQ. 0)
        READ(23, * , IOSTAT = OK)
        N = N + 1
      ENDDO

      ! EXIT WITH ERROR IF FILE HAS ONLY HEADER (OR EMPTY)
      IF (N .LE. 6) THEN
        IERR = 1
        CLOSE(23)
        RETURN
      ENDIF
      
      ! EXIT WITH ERROR IF FILE HAS NOT THE REQUIRED LENGTH
      IF (N .LT. L+5) THEN
        IERR = 2
        CLOSE(23)
        RETURN
      ENDIF

      ! ... OTHERWISE CONTINUE AND READ RESULTS
      REWIND(23)

      ! ALWAYS SKIP FIRST 5 LINES
      DO I = 1, 5
        READ(23, *)
      ENDDO

      ! FOR "*.dc" FILES WE NEED TO SKIP ONE MORE LINE
      IF (INDEX(FNAME, '.dc') .NE. 0) READ(23, *)

      ! START LOADING RESULTS (DISPERSION OR ELLIPTICITY CURVE)
      DO I = 1, L
        READ(23, * , IOSTAT = OK) FREQ(I), X(I)
      ENDDO
      IF (OK .NE. 0) THEN
        IERR = 2
        CLOSE(23)
        RETURN
      ENDIF

      CLOSE(23, IOSTAT = IERR)

    END SUBROUTINE GEOPSY2INV_SP

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE GEOPSY2INV_DP(FNAME, FREQ, X, IERR)

      ! SUBPROGRAM TO READ OUTPUT FILE FROM GEOPSY UTILITIES GPDC AND GPELL. I/O OPERATION CHECK ONLY WHEN OPENING/CLOSING FILE.

      CHARACTER(LEN=*),               INTENT(IN)    :: FNAME                               !< FILE NAME
      REAL(RDP),        DIMENSION(:), INTENT(INOUT) :: FREQ                                !< FREQUENCY
      REAL(RDP),        DIMENSION(:), INTENT(INOUT) :: X                                   !< EITHER DISPERSION CURVE OR ELLIPTICITY
      INTEGER,                        INTENT(OUT)   :: IERR                                !< ERROR FLAG
      INTEGER                                       :: OK, I, L, N                         !< COUNTERS

      !-----------------------------------------------------------------------------------------------------------------------------

      ! NO ERROR BY DEFAUL
      IERR = 0

      ! GET NUMBER OF EXPECTED FREQUENCY POINTS
      L = SIZE(FREQ)

      OPEN(23, FILE = FNAME, STATUS = 'OLD', FORM = 'FORMATTED', IOSTAT = IERR)

      IF (IERR .NE. 0) RETURN

      OK = 0
      N  = 0

      ! CRUDE CHECK TO SEE IF GEOPSY CALCULATED SOMETHING...
      DO WHILE (OK .EQ. 0)
        READ(23, * , IOSTAT = OK)
        N = N + 1
      ENDDO

      ! EXIT WITH ERROR IF FILE HAS ONLY HEADER (OR EMPTY)
      IF (N .LE. 6) THEN
        IERR = 1
        CLOSE(23)
        RETURN
      ENDIF
      
      ! EXIT WITH ERROR IF FILE HAS NOT THE REQUIRED LENGTH
      IF (N .LT. L+5) THEN
        IERR = 2
        CLOSE(23)
        RETURN
      ENDIF

      ! ... OTHERWISE CONTINUE AND READ RESULTS
      REWIND(23)

      ! ALWAYS SKIP FIRST 5 LINES
      DO I = 1, 5
        READ(23, *)
      ENDDO

      ! FOR "*.dc" FILES WE NEED TO SKIP ONE MORE LINE
      IF (INDEX(FNAME, '.dc') .NE. 0) READ(23, *)

      ! START LOADING RESULTS (DISPERSION OR ELLIPTICITY CURVE)
      DO I = 1, L
        READ(23, * , IOSTAT = OK) FREQ(I), X(I)
      ENDDO
      IF (OK .NE. 0) THEN
        IERR = 2
        CLOSE(23)
        RETURN
      ENDIF
      
      CLOSE(23, IOSTAT = IERR)

    END SUBROUTINE GEOPSY2INV_DP

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE IO
