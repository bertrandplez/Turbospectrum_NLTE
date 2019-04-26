C  INCLUDE FILE FOR SCMARCS   version 26-jan-89   BPz modifiee 200495 
C                                                     mkomp etc
C
C Maximal number of depth points in the model.
C
      INTEGER   NDP
      PARAMETER ( NDP=150 )
C
C Maximal number of rays for the solution of 
C the transfer equation.
C
      INTEGER   NRAYS
      PARAMETER ( NRAYS=150 ) 
*
*  equation of state parameters   (for tabs and absko)
*
*  mkomp  = max # of absorption coefficient contributors
*  mkompr = max # of temp. indep. absorption coefficient contributors
*
* mkomp=46, mkompr=17 before 22/09-1998. Changed to 47, 16 for the new FeI 
* and II opacities.. Not compatible with the old ones anymore.
*
      integer ifadim,kfadim,mkomp,mkompr,nabdim
      parameter (mkomp=52,mkompr=16,nabdim=15000,
     &  ifadim=(mkomp-mkompr)*ndp,kfadim=ifadim*3+mkompr*ndp)

