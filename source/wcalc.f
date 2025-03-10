*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
      REAL*8 FUNCTION WCALC(NH, NE, NHE, NS, TEMP)	
C  	
C  Computes the occupation probability for a H atom in 
C  state with effective principle quantum number NS in a plasma
C  enviroment with NH, NE, NHE number densities of H, ions,
C  and He atoms respectively.  This code assumes the perturbing
C  neutral H and He are in the ground state, (noting the hard
C  sphere approximation used is quite crude anyway) and that ions 
C  are predominantly singly ionized (ie. N(Z=1 ions) = Ne).
C
C  See eqn 4.71 Hummer & Milhalas (1988) 
C  Ions are now treated via Appendix A Hubeny, Hummer & Lanz (1994)
C  which is a fit to detailed calculation including correlations,
C  thus the temperature dependence
C
C  Sizes of neutral atoms adopted are sqrt(<r^2>)
C  
C  Coded by Paul Barklem and Kjell Eriksson, Aug 2003
C 	
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: sp = selected_real_kind(6, 37)
      real(dp) :: A, A0, BETAC, CHI, E, F, ION, IONH, K, NE, NEUTR
      real(dp) :: NH, NHE, NS, PI, RIH, TEMP, WION, WNEUTR, X
      DATA A0   / 5.29177d-9 /
      DATA E    / 4.803207d-10 /
      DATA IONH / 2.17991d-11 /
      DATA PI   / 3.1415927_dp /
C  
C  Neutral perturbers
C
      CHI=IONH/(NS*NS)
      RIH=SQRT(2.5*NS**4 + 0.5*NS*NS)*A0
      NEUTR=NH*(RIH + 1.73*A0)**3 + NHE*(RIH + 1.02*A0)**3
      WNEUTR = EXP(-4.*PI/3. * (NEUTR))
C
C  Charged perturbers
C      
      K=1.0_dp
      IF(NS .gt. 3.0_dp) THEN
        K = 16.0_dp/3.0_dp*(NS/(NS+1.0_dp))**2 *
     &      (NS + 7.0_dp/6.0_dp)/(NS*NS+NS+0.5_dp)
      ENDIF
c      ION  = NE*16.0_dp*(E**2/CHI/DSQRT(K))**3
c      WION = DEXP(-4.0_dp*PI/3.0_dp * (ION))
      IF ((NE.GT.10.0_dp).AND.(TEMP.GT.10.0_dp)) THEN  ! just in case!!! 
        A = 0.09_dp * NE**0.16667_dp / SQRT(TEMP)
        X = (1.0_dp + A)**3.15_dp
        BETAC = 8.3d14 * NE**(-0.667_dp) * K / (NS**4.0_dp)
        F = 0.14020_dp*X*BETAC*BETAC*BETAC /
     &      (1.0_dp+0.1285_dp*X*BETAC*DSQRT(BETAC))
        WION = F/(1.0_dp+F)
      ELSE 
        WION = 1.0_dp
      ENDIF 
C      
      WCALC = WION * WNEUTR 
      RETURN
      END
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
