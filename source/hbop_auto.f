      SUBROUTINE HBOP_AUTO(WAVE, NMAX, NH, NHE, NE, T, DOP, 
     *                     TOTAL, CONTIN,contonly,lineonly)             
C
C  wave   = interrogation wavelength in A               real*8
C  nmax   = highest principal quantum number to include integer
C  t      = kinetic temperature in K
C  ne     = electron number density in cm-3
C  nh     = number density of H I in cm-3
C  nhe    = number density of He I in cm-3
C  dop    = reduced Doppler width delta_lambda / lambda_0
C         = reduced Doppler width delta_nu / nu_0
C  total  = returns the total (line + continuous) absorption coefficient
C  contin = returns the continuous absorption coefficient
C
C  This code is a wrapper to HBOP, which will give the total bound-x
C  (bound-bound and bound-free) neutral H opacity at a given wavelength,
C  from the occupation probability formalism, WITHOUT the need for a 
C  line list.  This wrapper will simply produce a line list including 
C  all lines up to a highest principal quantum number NMAX, in the 
C  relevant spectral region (i.e. at the interrogation wavelength).
C
C  Note for efficiency one maybe wants to call HBOP directly if CONTIN 
C  is calculated separately.
C
C  At present this is limited to LTE; i.e. no possibility to pass non-
C  LTE populations.
C
C  All wavelengths will be calculated in air
C
C  **Important**: 
C     - Make sure NLEVELS in HBOP is at least as large as NMAX
C     - NMAX must be less than or equal to NMAXMAX, which is the maximum
C       n that can be handled, and must be set separately from 
C       NMAX to allow SAVE to work and avoid calculating lines
C       at every call
C     - NBF is the number of bound-free continua included and should be 
C       consistent with HBOP.  It is only here to catch if WAVE
C       is inconsistent with this.   Should be less or equal to 15 (see 
C       fucntion HBF
C
C  Paul Barklem, December 2024
C
      IMPLICIT NONE

      integer, parameter :: dp = selected_real_kind(15, 307)

      INTEGER I, J, NMAX, NBF, NMAXMAX, NLOW
      PARAMETER (NBF = 15, NMAXMAX = 100)
      integer  NLO(NMAXMAX), NUP(NMAXMAX), N
      real(dp) ::  NH, NHE, NE, T, CONTIN, TOTAL, DOP
      real(dp) ::  WAVE, WAVEH(NMAX), VACAIR
      real(dp) ::  EHYD(NMAXMAX), WCONTH(NMAXMAX)
      real(dp) ::  WLINEH(NMAXMAX,NMAXMAX)
      LOGICAL FIRST, contonly, lineonly
      SAVE FIRST, EHYD, WCONTH, WLINEH
      DATA FIRST/.TRUE./
C
C  On first call, compute and store:
C
C     energy levels in the H atom - EHYD(n)
C     bound-free continua break wavelengths in \AA\ in air - WCONTH(n)
C     line wavelengths in \AA in air - WLINEH(n_low,n_upper)
C
C        where n are relevant principal quantum numbers
C
      IF (NMAX.GT.NMAXMAX) THEN
        STOP 'HBOP_AUTO: NMAX exceeds NMAXMAX'
      ENDIF  

      IF (FIRST) THEN
        EHYD(1) =      0.000D0
        EHYD(2) =  82259.105D0
        EHYD(3) =  97492.302D0
        EHYD(4) = 102823.893D0
        EHYD(5) = 105291.651D0
        EHYD(6) = 106632.160D0
        EHYD(7) = 107440.444D0
        EHYD(8) = 107965.051D0
        DO 1 I = 9, NMAXMAX
        EHYD(I) = 109678.764D0 - 109677.576D0/(I*I)  
 1      CONTINUE 
        DO 2 I = 1, NMAXMAX
        WCONTH(I) = VACAIR(1.D8/(109678.764D0 - EHYD(I)))
 2      CONTINUE
        DO 4 I = 1, NMAXMAX-1
        DO 3 J = I+1, NMAXMAX
        WLINEH(I,J) = VACAIR(1.D8/(EHYD(J)-EHYD(I)))
 3      CONTINUE
 4      CONTINUE
        FIRST = .FALSE.
      ENDIF 
C
C  Find the relevant line series for this WAVE  
C
      NLOW = 1
!      DO 5 I = 2, NBF
      do 5 i=2,nmax
        IF (WAVE.GT.WCONTH(I)) NLOW = I
 5    CONTINUE

! BPz 27/02-2025 : do not stop. Instead hbop.f is modified to skip continuum calculation
!                  above NBF. However line opacity is computed for higher series.
!      IF (NLOW.EQ.NBF) THEN
!        STOP 'HBOP_AUTO: NLOW exceeds NBF' 
!      ENDIF
       if (nlow == nmax) then
         stop 'HBOP_AUTO: nlow exceeds nmax'
       endif
C
C  Extract the lines
C        
      N = NMAX - NLOW
      WAVEH(1:N) = WLINEH(NLOW,NLOW+1:NMAX)
      DO 7 I = 1, NMAX
        NLO(I) = NLOW
        NUP(I) = NLOW + I
 7    CONTINUE
C
C  Call HBOP (in LTE)
C
      CALL HBOP(WAVE, N, NLO(1:N), NUP(1:N), WAVEH(1:N), 
     &          NH, NHE, NE, T, DOP, 0., 0, NBF, TOTAL, CONTIN,
     &          contonly,lineonly)
C

      END      


