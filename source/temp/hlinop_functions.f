C***********************************************************************
      REAL*8 FUNCTION HPROFL(N,M,WAVE,WAVEH,T,XNE,H1FRC,HE1FRC,DOPPH)
C  
C hprofl = normalised line profile 
C          normalised to unity wrt integration over frequency
C n      = lower level principal quantum number
C m      = upper level principal quantum number
C wave   = interrogation wavelength in A           real*8
C waveh  = hydrogen line center wavelength in A    real*8
C t      = kinetic temperature in K
C xne    = electron number density in cm-3
C h1frc  = number density of H I in cm-3
C he1frc = number density of He I in cm-3
C dopph  = reduced Doppler width delta_lambda / lambda_0
C        = reduced Doppler width delta_nu / nu_0
C
C  Based on code by Deane Peterson and Bob Kurucz
C
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: NHaH2 = 49, NHaH2p = 74
      integer, parameter :: NLyaH2 = 95, NLyaH2p = 45, NLybH2 = 45
      integer, parameter :: NLybH2p = 26, NLygH2p = 27
      integer :: i, icut, ifins, ipos, k, m, m1, mmn, n, n1, nwid
      real(dp), PARAMETER :: CLIGHT = 2.99792458E18
      real(dp), PARAMETER :: CLIGHTCM = 2.99792458E10
      real(dp), PARAMETER :: PI = 3.14159265359, SQRTPI = 1.77245385
      real(dp) :: A, ALF(3), ASUM(100), ASUMLYMAN(100), BALFAC(3)
      real(dp) :: D, del, DELW, dlHaH2, dlHaH2p, dlLyaH2, dlLyaH2p
      real(dp) :: dlLybH2, dlLybH2p, dlLygH2p, DOP, dopph, EHYD(100)
      real(dp) :: F, FINEST(14), FINSWT(14), FO, FREQ, FREQNM, FREQSQ
      real(dp) :: gnm, grdwave, h1frc, HaH2(NHaH2), HaH2p(NHaH2p)
      real(dp) :: he1frc, hfnm, hfwid, hhw, HTOTAL, hwlor, hwstk
      real(dp) :: ISTAL(4), LNCOMP(4), LNGHAL(4), LyaH2(NLyaH2)
      real(dp) :: LyaH2p(NLyaH2p), LybH2(NLybH2)
      real(dp) :: LybH2p(NLybH2p), LygH2p(NLygH2p), radamp, RAYLCUT
      real(dp) :: resont, satb, SIGMA(3), STALPH(34), stark, stark1
      real(dp) :: STCOMP(5,4), STCPWT(5,4), STWTAL(34), t,t3nhe, t43
      real(dp) :: vbarh, vdw, WAVE, WAVEH, wavermp, wavermp0, wbrgd
      real(dp) :: xfac, xfacmax, xHaH2, xHaH2p, xknm, XKNMTB(4,3)
      real(dp) :: xlh1frc, xlHaH2, xlHaH2p, xlne, xlLyaH2
      real(dp) :: xlLyaH2p, xlLybH2, xlLybH2p, xlLygH2p, xlnHaH2
      real(dp) :: xlnHaH2p, xlnlyah2, xlnlyah2p, xlnlybh2, xlnlybh2p
      real(dp) :: xlnlygh2p, xLyaH2, xLyaH2p, xLybH2, xLybh2p, xLygH2p
      real(dp) :: xm, xm2, xm2mn2, xmn2, xn, XN2, xne, xself, xself0
      real(dp) :: xstark, xstark0, xstarkb
      DATA   xlLyaH2,  xlLyaH2p,   xlLybH2,  xlLybH2p,  xlLygH2p /
     &     1278.7_dp, 1290.3_dp, 1034.6_dp, 1032.4_dp, 975.72_dp /
      DATA    xlHaH2,   xlHaH2p /
     &     6649.3_dp, 6819.8_dp /
      DATA   dlLyaH2, dlLyaH2p, dlLybH2, dlLybH2p, dlLygH2p /
     &        4.0_dp,   3.0_dp,  3.0_dp,   2.0_dp,  1.00_dp /
      DATA    dlHaH2, dlHaH2p /
     &       38.0_dp, 42.4_dp /
      DATA   xLNLyaH2, xLNLyaH2p, xLNLybH2, xLNLybH2p, xLNLygH2p /
     &       -18.0_dp,  -17.0_dp, -18.0_dp,  -17.0_dp,  -17.0_dp /
      DATA    xLNHaH2,  xLNHaH2p /
     &       -17.0_dp,  -17.0_dp /
      LOGICAL LYMANALF
      SAVE
      COMMON / CLEVELS / EHYD
C
C  Einstein A-value sums for H lines
C
      DATA ASUM/
     & 0.000d+00, 4.696d+08, 9.980d+07, 3.017d+07, 1.155d+07, 5.189d+06,
     & 2.616d+06, 1.437d+06, 8.444d+05, 5.234d+05, 3.389d+05, 2.275d+05,
     & 1.575d+05, 1.120d+05, 8.142d+04, 6.040d+04, 4.560d+04, 3.496d+04,
     & 2.719d+04, 2.141d+04, 1.711d+04, 1.377d+04, 1.119d+04, 9.166d+03,
     & 7.572d+03, 6.341d+03, 5.338d+03, 4.523d+03, 3.854d+03, 3.302d+03,
     & 2.844d+03, 2.460d+03, 2.138d+03, 1.866d+03, 1.635d+03, 1.438d+03,
     & 1.269d+03, 1.124d+03, 9.983d+02, 8.894d+02, 7.947d+02, 7.120d+02,
     & 6.396d+02, 5.759d+02, 5.198d+02, 4.703d+02, 4.263d+02, 3.873d+02,
     & 3.526d+02, 3.215d+02, 2.938d+02, 2.689d+02, 2.465d+02, 2.264d+02,
     & 2.082d+02, 1.918d+02, 1.769d+02, 1.634d+02, 1.512d+02, 1.400d+02,
     & 1.298d+02, 1.206d+02, 1.121d+02, 1.043d+02, 9.720d+01, 9.066d+01,
     & 8.465d+01, 7.912d+01, 7.403d+01, 6.933d+01, 6.498d+01, 6.097d+01,
     & 5.725d+01, 5.381d+01, 5.061d+01, 4.765d+01, 4.489d+01, 4.232d+01,
     & 3.994d+01, 3.771d+01, 3.563d+01, 3.369d+01, 3.188d+01, 3.019d+01,
     & 2.860d+01, 2.712d+01, 2.572d+01, 2.442d+01, 2.319d+01, 2.204d+01,
     & 2.096d+01, 1.994d+01, 1.898d+01, 1.808d+01, 1.722d+01, 1.642d+01,
     & 1.566d+01, 1.495d+01, 1.427d+01, 1.363d+01/
C
C  For Lyman lines only the s-p transition is allowed.
C
      DATA ASUMLYMAN/
     & 0.000d+00, 6.265d+08, 1.897d+08, 8.126d+07, 4.203d+07, 2.450d+07,
     & 1.236d+07, 8.249d+06, 5.782d+06, 4.208d+06, 3.158d+06, 2.430d+06,
     & 1.910d+06, 1.567d+06, 1.274d+06, 1.050d+06, 8.752d+05, 7.373d+05,
     & 6.269d+05, 5.375d+05, 4.643d+05, 4.038d+05, 3.534d+05, 3.111d+05,
     & 2.752d+05, 2.447d+05, 2.185d+05, 1.959d+05, 1.763d+05, 1.593d+05,
     & 1.443d+05, 1.312d+05, 1.197d+05, 1.094d+05, 1.003d+05, 9.216d+04,
     & 8.489d+04, 7.836d+04, 7.249d+04, 6.719d+04, 6.239d+04, 5.804d+04,
     & 5.408d+04, 5.048d+04, 4.719d+04, 4.418d+04, 4.142d+04, 3.888d+04,
     & 3.655d+04, 3.440d+04, 3.242d+04, 3.058d+04, 2.888d+04, 2.731d+04,
     & 2.585d+04, 2.449d+04, 2.322d+04, 2.204d+04, 2.094d+04, 1.991d+04,
     & 1.894d+04, 1.804d+04, 1.720d+04, 1.640d+04, 1.566d+04, 1.496d+04,
     & 1.430d+04, 1.368d+04, 1.309d+04, 1.254d+04, 1.201d+04, 1.152d+04,
     & 1.105d+04, 1.061d+04, 1.019d+04, 9.796d+03, 9.419d+03, 9.061d+03,
     & 8.721d+03, 8.398d+03, 8.091d+03, 7.799d+03, 7.520d+03, 7.255d+03,
     & 7.002d+03, 6.760d+03, 6.530d+03, 6.310d+03, 6.100d+03, 5.898d+03,
     & 5.706d+03, 5.522d+03, 5.346d+03, 5.177d+03, 5.015d+03, 4.860d+03,
     & 4.711d+03, 4.569d+03, 4.432d+03, 4.300d+03 /
C
      DATA N1/0/, M1/0/
C
C  Fine structure components for alpha lines in FREQ*10**-7
C
      DATA STALPH / -730.0_dp,  370.0_dp,  188.0_dp,  515.0_dp,
     &    327.0_dp,  619.0_dp, -772.0_dp, -473.0_dp, -369.0_dp,
     &    120.0_dp,  256.0_dp,  162.0_dp,  285.0_dp, -161.0_dp,
     &    -38.3_dp,   6.82_dp, -174.0_dp, -147.0_dp, -101.0_dp,
     &    -77.5_dp,   55.0_dp,  126.0_dp,   75.0_dp,  139.0_dp,
     &    -60.0_dp,    3.7_dp,   27.0_dp,  -69.0_dp,  -42.0_dp,
     &    -18.0_dp,   -5.5_dp,   -9.1_dp,  -33.0_dp,   -24._dp /
C
C  Alpha component weights
C
      DATA STWTAL / 1.0_dp,2.0_dp,1.0_dp,2.0_dp,1.0_dp,2.0_dp,1.0_dp,
     &              2.0_dp,3.0_dp,1.0_dp,2.0_dp,1.0_dp,2.0_dp,1.0_dp,
     &              4.0_dp,6.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp,1.0_dp,
     &              2.0_dp,1.0_dp,2.0_dp,1.0_dp,4.0_dp,6.0_dp,1.0_dp,
     &              7.0_dp,6.0_dp,4.0_dp,4.0_dp,4.0_dp,5.0_dp /
      DATA ISTAL /1, 3, 10, 21/
      DATA LNGHAL/2, 7, 11, 14/
C
C  Fine structure for M.EQ.INFINITY IN FREQ*10**-7
C
      DATA STCOMP /  0.0_dp,   0.0_dp,    0.0_dp,    0.0_dp,   0.0_dp,
     &             468.0_dp, 576.0_dp, -522.0_dp,    0.0_dp,   0.0_dp,
     &             260.0_dp, 290.0_dp,  -33.0_dp, -140.0_dp,   0.0_dp,
     &             140.0_dp, 150.0_dp,   18.0_dp,  -27.0_dp,  -51.0_dp /
C
C  Weights for fine structure components
C
      DATA STCPWT/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,
     &            1.0_dp, 1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp,
     &            1.0_dp, 1.0_dp, 4.0_dp, 3.0_dp, 0.0_dp,
     &            1.0_dp, 1.0_dp, 4.0_dp, 6.0_dp, 4.0_dp /
C
      DATA LNCOMP/1, 3, 4, 5/
C
      DATA XKNMTB / 0.0001716_dp, 0.0090190_dp, 0.1001000_dp,
     &              0.5820000_dp, 0.0005235_dp, 0.0177200_dp,
     &              0.1710000_dp, 0.8660000_dp, 0.0008912_dp,
     &              0.0250700_dp, 0.2230000_dp, 1.0200000_dp /
C
C  Data for satellites in the wings of hydrogen lines arising from collisions
C  with H and H+ are taken from N. Allard's web-site:
C     http://mygepi.obspm.fr/~allard/lymantables.html
C  I have only including satellites from binary collisions <=> linear in N.
C
C  For, e.g., neutral perturbers on Lyman\alpha, the data are stored in LyaH2
C  with NLyaH2 elements, and cover the wavelength region
C           \lambda/[\AA] = [xlLyaH2; xlLyaH2+NLyaH2*dlLyaH2]
C  The data are given for a perturber density of 10^xLNLyaH2. The satellites
C  are assumed to be independent of temperature.
C
C  The transition to the self-broadened profile (for neutral perturbers)
C  and to the Stark profile (for charged ones) is accomplished with a
C  linear or sinusoidal bridging function between (the locally defined)
C  WAVERMP and the blue edge of the data. The self-broadened and Stark
C  profiles are allowed to increase (linearly from line-center to blue edge
C  of the satellite data) by a factor up to (the locally defined) XFACMAX
C  in order to join with the satellite profile - otherwise the satellite is
C  decreased to join smoothly with the former. 
C
C  The Ly\alpha and \beta lines are covered in
C  N.F Allard, A. Royer, J.F. Kielkopf and N. Feautrier,
C  1999, Phys. Rev. A, 60, 1021, 1033 
C  http://adsabs.harvard.edu/abs/1999PhRvA..60.1021A
C
C  Lyman \alpha satellites from H collisions  for N(H) = 1e18 cm^-3
C  wavelength  = 1278.7\AA  + (N-1)*4.    N=1,95   1654.7\AA
C  Major satellite at 1580\AA and minor ones at 1430 and 1345\AA.
C
      DATA LyaH2/
     & -19.079_dp, -19.127_dp, -19.172_dp, -19.213_dp, -19.247_dp,
     & -19.271_dp, -19.294_dp, -19.320_dp, -19.345_dp, -19.361_dp,
     & -19.367_dp, -19.367_dp, -19.369_dp, -19.371_dp, -19.370_dp,
     & -19.365_dp, -19.359_dp, -19.357_dp, -19.360_dp, -19.366_dp,
     & -19.372_dp, -19.378_dp, -19.383_dp, -19.388_dp, -19.392_dp,
     & -19.394_dp, -19.392_dp, -19.386_dp, -19.375_dp, -19.361_dp,
     & -19.347_dp, -19.332_dp, -19.317_dp, -19.303_dp, -19.289_dp,
     & -19.275_dp, -19.258_dp, -19.240_dp, -19.232_dp, -19.237_dp,
     & -19.248_dp, -19.264_dp, -19.281_dp, -19.294_dp, -19.301_dp,
     & -19.307_dp, -19.312_dp, -19.320_dp, -19.330_dp, -19.340_dp,
     & -19.352_dp, -19.361_dp, -19.366_dp, -19.362_dp, -19.347_dp,
     & -19.324_dp, -19.294_dp, -19.261_dp, -19.229_dp, -19.200_dp,
     & -19.175_dp, -19.154_dp, -19.135_dp, -19.118_dp, -19.101_dp,
     & -19.084_dp, -19.064_dp, -19.043_dp, -19.019_dp, -18.994_dp,
     & -18.969_dp, -18.946_dp, -18.927_dp, -18.911_dp, -18.900_dp,
     & -18.892_dp, -18.890_dp, -18.893_dp, -18.901_dp, -18.915_dp,
     & -18.933_dp, -18.957_dp, -18.985_dp, -19.018_dp, -19.056_dp,
     & -19.098_dp, -19.144_dp, -19.193_dp, -19.245_dp, -19.300_dp,
     & -19.358_dp, -19.420_dp, -19.483_dp, -19.547_dp, -19.609_dp /
C
C  Lyman \alpha satellites from H+ collisions  for N(H+) = 1e17 cm^-3
C  wavelength  = 1290.3\AA  + (N-1)*3.    N=1,45   1654.7\AA
C  Major satellite at 1390\AA and minor ones at 1295 and 1330\AA.
C
      DATA LyaH2p/
     & -18.553_dp, -18.583_dp, -18.613_dp, -18.641_dp, -18.672_dp,
     & -18.703_dp, -18.725_dp, -18.736_dp, -18.743_dp, -18.752_dp,
     & -18.760_dp, -18.764_dp, -18.769_dp, -18.785_dp, -18.814_dp,
     & -18.845_dp, -18.868_dp, -18.881_dp, -18.895_dp, -18.915_dp,
     & -18.935_dp, -18.934_dp, -18.906_dp, -18.860_dp, -18.816_dp,
     & -18.784_dp, -18.761_dp, -18.739_dp, -18.711_dp, -18.672_dp,
     & -18.633_dp, -18.595_dp, -18.567_dp, -18.564_dp, -18.582_dp,
     & -18.622_dp, -18.682_dp, -18.759_dp, -18.859_dp, -18.974_dp,
     &  -19.105_dp, -19.249_dp, -19.400_dp, -19.548_dp, -19.689_dp /

C  Lyman \beta satellites from H collisions  for N(H) = 1e18 cm^-3
C  wavelength  = 1034.6\AA  + (N-1)*3.    N=1,45   1166.6\AA
C  Major satellite at 1150 and minor ones at 1087 and 1110\AA.
C
      DATA LybH2/
     & -18.464_dp, -18.632_dp, -18.732_dp, -18.810_dp, -18.925_dp,
     & -19.006_dp, -19.055_dp, -19.074_dp, -19.094_dp, -19.099_dp,
     & -19.087_dp, -19.079_dp, -19.059_dp, -19.040_dp, -19.023_dp,
     & -19.004_dp, -18.988_dp, -19.002_dp, -19.063_dp, -19.116_dp,
     & -19.106_dp, -19.055_dp, -19.039_dp, -19.052_dp, -19.034_dp,
     & -18.969_dp, -18.887_dp, -18.839_dp, -18.834_dp, -18.846_dp,
     & -18.856_dp, -18.843_dp, -18.796_dp, -18.703_dp, -18.572_dp,
     & -18.449_dp, -18.360_dp, -18.255_dp, -18.180_dp, -18.199_dp,
     & -18.317_dp, -18.514_dp, -18.778_dp, -19.067_dp, -19.344_dp /
C
C  Lyman \beta satellites from H+ collisions  for N(H+) = 1e17 cm^-3
C  wavelength  = 1032.4\AA  + (N-1)*2.    N=1,26   1082.4\AA
C  Major satellites at 1060 and 1077\AA.
C
      DATA LybH2p/
     & -15.321_dp, -15.499_dp, -15.648_dp, -16.080_dp, -16.365_dp,
     & -16.470_dp, -16.540_dp, -16.599_dp, -16.648_dp, -16.702_dp,
     & -16.741_dp, -16.677_dp, -16.614_dp, -16.541_dp, -16.513_dp,
     & -16.636_dp, -16.862_dp, -17.092_dp, -17.170_dp, -17.123_dp,
     & -16.937_dp, -16.827_dp, -16.722_dp, -16.740_dp, -16.978_dp,
     & -17.356_dp /
C
C  Lyman \gamma satellites from H+ collisions  for N(H+) = 1e17 cm^-3
C  wavelength  = 975.72\AA  + (N-1)*1.    N=1,27   1001.72\AA
C  Major satellite at 992\AA and a minor one at 976\AA.
C  N.F. Allard_dp, I. Noselidze and J.W. Kruk_dp, 2009_dp, A&A_dp, 506_dp, 993_dp,
C  http://adsabs.harvard.edu/abs/2009A%26A...506..993A
C
      DATA LygH2p/
     & -14.173_dp, -14.380_dp, -14.553_dp, -14.751_dp, -14.944_dp,
     & -15.184_dp, -15.394_dp, -15.496_dp, -15.547_dp, -15.601_dp,
     & -15.676_dp, -15.748_dp, -15.807_dp, -15.810_dp, -15.812_dp,
     & -15.805_dp, -15.773_dp, -15.789_dp, -15.863_dp, -15.968_dp,
     & -16.020_dp, -16.039_dp, -16.152_dp, -16.395_dp, -16.742_dp,
     & -17.108_dp, -17.444_dp /
C
C  Balmer \alpha satellites from H collisions  for N(H) = 1e17 cm^-3
C  wavelength  = 6649.3\AA  + (N-1)*38.0    N=1,49    8473.3\AA
C  Broad satellite at 7500\AA.
C  From J.F. Kielkopf, N.F. Allard and A. Decrette, 2002,
C  Europ. Phys. J. D, 18(1),51
C  http://adsabs.harvard.edu/abs/2002EPJD...18...51K
C
      DATA HaH2/
     & -17.092_dp, -17.345_dp, -17.519_dp, -17.682_dp, -17.856_dp,
     & -17.961_dp, -18.086_dp, -18.186_dp, -18.291_dp, -18.370_dp,
     & -18.422_dp, -18.501_dp, -18.559_dp, -18.600_dp, -18.632_dp,
     & -18.638_dp, -18.630_dp, -18.606_dp, -18.585_dp, -18.562_dp,
     & -18.540_dp, -18.521_dp, -18.509_dp, -18.500_dp, -18.497_dp,
     & -18.503_dp, -18.515_dp, -18.529_dp, -18.548_dp, -18.573_dp,
     & -18.607_dp, -18.645_dp, -18.692_dp, -18.736_dp, -18.787_dp,
     & -18.840_dp, -18.897_dp, -18.961_dp, -19.025_dp, -19.099_dp,
     & -19.177_dp, -19.262_dp, -19.365_dp, -19.466_dp, -19.569_dp,
     & -19.685_dp, -19.784_dp, -19.931_dp, -20.145_dp /
C
C     Balmer \alpha satellites from H+ collisions  for N(H+) = 1e17 cm^-3
C     wavelength  = 6819.8\AA  + (N-1)*42.4    N=1,76    9999.8\AA
C     Major satellite at 8600\AA and smaller ones at 7050 and 9600\AA.
C
      DATA HaH2p/
     & -15.770_dp, -16.020_dp, -16.286_dp, -16.448_dp, -16.516_dp,
     & -16.551_dp, -16.622_dp, -16.725_dp, -16.841_dp, -16.963_dp,
     & -17.055_dp, -17.132_dp, -17.180_dp, -17.235_dp, -17.255_dp,
     & -17.278_dp, -17.305_dp, -17.294_dp, -17.300_dp, -17.301_dp,
     & -17.278_dp, -17.266_dp, -17.266_dp, -17.254_dp, -17.263_dp,
     & -17.286_dp, -17.290_dp, -17.301_dp, -17.332_dp, -17.346_dp,
     & -17.321_dp, -17.290_dp, -17.277_dp, -17.267_dp, -17.228_dp,
     & -17.164_dp, -17.109_dp, -17.076_dp, -17.055_dp, -17.030_dp,
     & -17.003_dp, -16.989_dp, -17.002_dp, -17.044_dp, -17.111_dp,
     & -17.204_dp, -17.321_dp, -17.456_dp, -17.597_dp, -17.730_dp,
     & -17.849_dp, -17.956_dp, -18.048_dp, -18.119_dp, -18.162_dp,
     & -18.182_dp, -18.194_dp, -18.206_dp, -18.216_dp, -18.222_dp,
     & -18.221_dp, -18.212_dp, -18.198_dp, -18.186_dp, -18.179_dp,
     & -18.178_dp, -18.184_dp, -18.197_dp, -18.217_dp, -18.242_dp,
     & -18.267_dp, -18.295_dp, -18.321_dp, -18.345_dp /
C
C  Most model atmosphere codes include Rayleigh scattering by H atoms 
C  elsewhere, eg. quantum mechanical calculations. This parameter cuts
C  the Lyman alpha natural absorption at this chosen point.  
C
C Changed from 1240 to 1400 by BPz 04/10-2007
      PARAMETER (RAYLCUT = 1400.0D0) ! in Angstroms
C
C  Data for self-broadening from calculations of Barklem, Piskunov and 
C  O'Mara (2000, A&A 363, 1091).
C  SIGMA is in m^2.
C  BALFAC = (4/Pi)**(alf/2) Gamma(2-2/alf) the leading factor for
C  computing the width.
C
      DATA ALF    / 0.677_dp, 0.455_dp, 0.380_dp /
      DATA BALFAC / 0.97875_dp, 0.97659_dp, 0.97794_dp /
      DATA SIGMA  / 3.304d-18, 6.497d-18, 1.178d-17 /
C
C  Precompute variables depending on current physical conditions
C
      T43  = (T/10000.0_dp)**0.3_dp
      T3NHE= T43*HE1FRC
      FO   = 1.25d-9*XNE**0.66667_dp ! Holtsmark normal field strength
      xLH1FRC=dlog10(H1FRC)
      XLNE =  dlog10(XNE)
C
C  Convert wavelengths to frequencies and compute detunings
C
      DELW = WAVE-WAVEH
      FREQ = CLIGHT/WAVE
      FREQNM = CLIGHT/WAVEH
      DEL = FREQ-FREQNM
C
C  Variables dependent on line - compute first time only
C
      IF((N.NE.N1).OR.(M.NE.M1)) THEN
c        print '("Stat.Weight:",2I3,4f10.6)',
c    &      HFNM(1,N),HFNM(1,M),HFNM(N,M),HFNM(M,N)
         N1 = N
         M1 = M
         MMN = M-N
         XN = N
         XN2 = XN*XN
         XM = M
         XM2 = XM*XM
         XMN2 = XM2*XN2
         XM2MN2 = XM2-XN2
         GNM = XM2MN2/XMN2
         IF ((MMN.LE.3).AND.(N.LE.4)) THEN
            XKNM = XKNMTB(N,MMN)
         ELSE
            XKNM = 5.5d-5/GNM*XMN2/(1.0_dp+0.13_dp/FLOAT(MMN))
         END IF
C
C  Lyman alpha wings require special treatment
C
         LYMANALF = .FALSE.
         IF ((N.EQ.1).AND.(M.EQ.2)) LYMANALF = .TRUE.
C
C  Natural damping taken from tables 
C  dnu/ nu_0 = 1/2 * 1/nu_0 * 1/2 pi * Gamma
C
         IF (N.EQ.1) THEN
            RADAMP = ASUMLYMAN(M)
         ELSE
            RADAMP = ASUM(N)+ASUM(M)
         ENDIF
         RADAMP = RADAMP/FREQNM/(4.0_dp*PI)
C
C  Resonance broadening following Ali & Griem (1966, Phys Rev 144, 366).
C  RESONT is dnu/nu per unit H density. Only the lower state is included 
C  (p-d approx for Balmer lines).  For N > 3 p-states play a very small 
C  role, and van der Waals might even dominate, there is however no 
C  available theory for this.  The lower state resonance broadening is
C  used as a guess.
C
         IF (N.NE.1) then
           RESONT = HFNM(1,N)/(1.0_dp-1.0_dp/ XN2)
         else
           resont = hfnm(1,m)/(1.0_dp-1.0_dp/ xm2)
         endif
         RESONT = RESONT * 2.07d-24/GNM
         VDW = 4.45d-26/GNM*(XM2*(7.0_dp*XM2+5.0_dp))**0.4_dp
         STARK = 1.6678d-18*FREQNM*XKNM
C
C  Approximately include fine structure.  Exact pattern for alpha lines, 
C  M infinite pattern used for others.
C
         IF (N.GT.4 .OR. M.GT.10) THEN
            IFINS = 1
            FINEST(1) = 0.0_dp
            FINSWT(1) = 1.0_dp
         ELSE IF (MMN.GT.1) THEN
            IFINS = LNCOMP(N)
            DO 1 I = 1,IFINS
            FINEST(I) = STCOMP(I, N)*1.0D7
            FINSWT(I) = STCPWT(I, N)/XN2
   1        CONTINUE
         ELSE
C
C  eg: Ly alpha IFINS=2, IPOS=1, FINEST=-7.3E9,3.7E9, FINSWT=1/3, 2/3
C
            IFINS = LNGHAL(N)
            IPOS = ISTAL(N)
            DO 2 I = 1,IFINS
            K = IPOS-1+I
            FINEST(I) = STALPH(K)*1.0D7
            FINSWT(I) = STWTAL(K)/XN2/3.0D0
   2        CONTINUE
         END IF
      END IF
C
C  Now compute the profile for the given physical conditions.
C
C  Firstly half-widths: DOPPH, HWSTK, HWLOR, are dnu/nu. DOP is the 
C  doppler width dnu.
C
      HWSTK = STARK*FO
      DOP = DOPPH*FREQNM
C      
C  For low Balmer lines use Barklem, Piskunov & O'Mara (2000, A&A 363, 
C  1091) instead of Ali & Griem for self-broadening. 1E6 factor is due 
C  to SI -> CGS.
C
      IF ((N.EQ.2) .AND. (M.LE.5)) THEN
         VBARH = SQRT(42008.0_dp*T)
         RESONT = BALFAC(M-N) * 1.0d4*SIGMA(M-N) 
     &                        * (VBARH/1.0d4)**(1.0_dp-ALF(M-N))
         RESONT = RESONT/FREQNM*1.0d6/2.0_dp/PI
      ENDIF
C
      HWLOR = RESONT*H1FRC + VDW*T3NHE + RADAMP
C
C  *Approximate* convolution is achieved by two cases:
C   (a) close to the core -> take the dominant mechanism.
C   (b) outside the core -> add all mechanisms.
C
C  Mechanisms are divided into three:
C   (1) Doppler core with fine structure.
C   (2) (Lorentz) damping parts - self, natural and helium.
C   (3) Stark broadening.
C

C
C  Determine the largest half-width 
C    NWID=1, Doppler
C        =2, Lorentz
C        =3, Stark
C
      IF (DOPPH.GE.HWLOR.AND.DOPPH.GE.HWSTK) THEN
        NWID = 1
        HFWID = DOPPH*FREQNM
      ELSE IF(HWLOR.GE.DOPPH.AND.HWLOR.GE.HWSTK) THEN
        NWID = 2
        HFWID = HWLOR*FREQNM
      ELSE
        NWID = 3
        HFWID = HWSTK*FREQNM
      END IF
C
C  HTOTAL is the profile.
C
      HTOTAL=0.0D0
C
C  Case 1: Doppler core including fine structure 
C
      IF ((ABS(DEL).GT.HFWID).OR.(NWID.EQ.1)) THEN
          DO 3 I = 1,IFINS
          D = DABS(FREQ-FREQNM-FINEST(I))/DOP
          IF(D.LE.7.) HTOTAL = HTOTAL + EXP(-D*D)*FINSWT(I)/(SQRTPI*DOP)
   3      CONTINUE
      ENDIF
C  
C  Case 2: Self-broadening section - includes natural and vdW due to 
C  helium for convenience.
C
      IF ((ABS(DEL).GT.HFWID).OR.(NWID.EQ.2)) THEN
C         
C  All lines except Lyman alpha:
C
C  The full classical expression for the line shape is used (eqn 4-73 in
C  Aller). This extends the profile to large detuning, but note the 
C  impact approximation is not valid.  
C
C  The classical damping constant is frequency dependent. Assume similar 
C  dependence of quantum mechanical damping (cf. eqn 4-116 Aller).
C
         IF (LYMANALF) THEN
C
C  Lyman alpha:
C
C  Natural broadening of Lyman alpha
C  This is cut at some specified detuning in the red.
C
            IF (WAVE.LT.RAYLCUT) THEN 
               HHW = RADAMP*FREQNM   
               HHW = HHW * FREQ*FREQ/(FREQNM*FREQNM)    
               FREQSQ = FREQ*FREQ
               HTOTAL = HTOTAL + 4.0_dp*FREQSQ*HHW/PI/
     &           ((DEL*(FREQNM+FREQ))**2.0_dp + FREQSQ*4.*HHW**2.0_dp)
            ENDIF
C
C  Self broadening of Lyman alpha:
C 
C  Red wing of Lyman Alpha before the satellite region and blue wing 
C  (ie. redward of ~1190 A).
C
C  Far red wing of Lyman Alpha in satellite region:
C
c  Changed by BPz on 25/09-2007 to stop extrapolation where Doyle's H+H 
c  CIA takes over in jonabs_vac.dat (detabs.f) : 1750A = 1d8/57142.8571
c
c  An analytical fit to Doyle (1968)'s results for 1750\AA:
c  log10[(H+H)(1750\AA)] =
c   -256.85471+xLT*(233.23529+xLT*(-76.960342+xLT*(11.304402-xLT*0.62315132)))
c  for logT > 4.7,  log10[(H+H)(1750\AA)] = 8.8761660 + (xLT-4.7d0)*0.16103775
c  with xLT = dlog10(T). From ~/convdat/{tables/subs.17elm.dat,tabcode/detabs.f}
c
cccc        IF      (FREQ.LT.20000.*CLIGHTCM) THEN 
            IF      (WAVE.GT.1750.0_dp) THEN                                 IGNORE
               XSELF = 0.0_dp
            ELSE IF (WAVE.GT.xlLyaH2) THEN                                INTERP
c                    wave > 1278.70\AA
               ICUT   = min0(NLyaH2-2,max0(0,int(
     &                      (WAVE-xlLyaH2)/dlLyaH2))) + 1
               GRDWAVE= xlLyaH2 + (ICUT-1)*dlLyaH2
               XLyaH2 = (LyaH2(ICUT+1)-LyaH2(ICUT))
     &               /dlLyaH2*(WAVE-GRDWAVE) + LyaH2(ICUT)
               XLyaH2 = 10.0_dp**(XLyaH2+xLH1FRC+xLNLyaH2)
               XSELF  = XLyaH2
            ELSE                                                          SCALE
C                    wave < 1278.70\AA
C  
C  The impact approximation breaks down quickly here validity is ~ 1 A 
C  at T~5000 K (Lortet & Roueff 1969, A&A 3, 462) but we use Ali & Griem 
C  anyway (~ Barklem et al for 2p).
C
C  A factor 1.17 is needed to match this region to the far red wing 
C  (below) at 1278.70\AA. We use a linear transition
C  from the line-center to the bluest point of the wing-table (LyaH2).
C  01.02.2012 Changed from 1.17 -> 1.0, instead adjusted the (H+H) data/RT
c              XFAC   = 1.+0.13*(DABS(DELW)/
c    &                     (1d8/(EHYD(M)-xlLyaH2) - 1d8/EHYD(M)))
               XFAC   = 1.0
c              IF (FREQ.GT.FREQNM) XFAC = 1.0
               HHW    = (XFAC*RESONT*H1FRC + VDW*T3NHE) * FREQNM
               XSELF  = HHW/PI/(DEL*DEL+HHW*HHW)
            ENDIF
C
C  Lyman beta satellites
C
         ELSE IF ((N.EQ.1).AND.(M.EQ.3)) THEN
            IF (WAVE.LT.RAYLCUT) THEN 
               HHW = RADAMP*FREQNM   
               HHW = HHW * FREQ*FREQ/(FREQNM*FREQNM)    
               FREQSQ = FREQ*FREQ
               HTOTAL = HTOTAL + 4.0_dp*FREQSQ*HHW/PI/
     &           ((DEL*(FREQNM+FREQ))**2.0_dp + FREQSQ*4.*HHW**2.0_dp)
            ENDIF
            IF      (WAVE.GT.1200.0_dp) THEN
               XSELF  = 0.
            ELSE IF (WAVE.GT.xlLybH2) THEN
               ICUT   = min0(NLybH2-2,max0(0,int(
     &                      (WAVE-xlLybH2)/dlLybH2))) + 1
               GRDWAVE= xlLybH2 + (ICUT-1)*dlLybH2
               XLybH2 = (LybH2(ICUT+1)-LybH2(ICUT))
     &               /dlLybH2*(WAVE-GRDWAVE) + LybH2(ICUT)
               XLybH2 = 10.0_dp**(XLybH2+xLH1FRC+xLNLybH2)
               XSELF  = XLybH2
            ELSE
               HHW    = (RESONT*H1FRC + VDW*T3NHE) * FREQNM
               XSELF  = HHW/PI/(DEL*DEL+HHW*HHW)
            ENDIF
C
C  Balmer alpha satellites
C
         ELSE IF ((N.EQ.2).AND.(M.EQ.3)) THEN
            WAVERMP = 7200.0_dp
c                  \lambda > 9999.8\AA
            IF      (WAVE.GT.(xlHaH2+dlHaH2*(NHaH2-1.0_dp))) THEN
c  Linear extrapolation in logarithmic absorption
               XHaH2  = HaH2(NHaH2) + (HaH2(NHaH2)-HaH2(NHaH2-1))/dlHaH2
     &                   *(WAVE - (xlHaH2+dlHaH2*(NHaH2-1.0_dp)))
               XHaH2  = 10.0_dp**(XHaH2+xLH1FRC+xLNHaH2)
               XSELF  = XHaH2
c                  \lambda > 6819.8\AA
            ELSE IF (WAVE.GT.xlHaH2) THEN
c  Linear interpolation in logarithmic absorption
               ICUT   = min0(NHaH2-2,max0(0,int(
     &                      (WAVE-xlHaH2)/dlHaH2))) + 1
               GRDWAVE= xlHaH2 + (ICUT-1)*dlHaH2
               XHaH2  = (HaH2(ICUT+1)-HaH2(ICUT))
     &                /dlHaH2*(WAVE-GRDWAVE) + HaH2(ICUT)
               XHaH2  = 10.0_dp**(XHaH2+xLH1FRC+xLNHaH2)
               XSELF  = XHaH2
c                  \lambda < 7200.0\AA
             IF     (WAVE.LT.WAVERMP) THEN
c
c  Ramp down that bluest part of the satellite
               HHW = HWLOR * FREQ*FREQ/FREQNM
               FREQSQ = FREQ*FREQ
               XSELF0 = 4.0_dp*FREQSQ*HHW/PI/
     &         ((DEL*(FREQNM+FREQ))**2.0_dp + FREQSQ*4.0_dp*HHW**2.0_dp)
               WBRGD  = (WAVERMP - WAVE)/(WAVERMP-xlHaH2)
               XSELF  = XSELF0*WBRGD + XSELF*(1.0_dp-WBRGD)
             ENDIF
c                  \lambda <= 6819.8\AA  - blue side of H+ satellites
            ELSE
c              HHW    = (RESONT*H1FRC + VDW*T3NHE) * FREQNM
c              XSELF  = HHW/PI/(DEL*DEL+HHW*HHW)
c  This is just copied from the default case, below.
               HHW    = HWLOR*FREQNM 
               HHW    = HHW * FREQ*FREQ/(FREQNM*FREQNM)    
               FREQSQ = FREQ*FREQ
               XSELF  = 4.0_dp*FREQSQ*HHW/PI/
     &           ((DEL*(FREQNM+FREQ))**2.0_dp + FREQSQ*4.*HHW**2.0_dp)
            ENDIF
         ELSE
C
c           HHW = HWLOR*FREQNM 
c           HHW = HHW * FREQ*FREQ/(FREQNM*FREQNM)    
            HHW = HWLOR * FREQ*FREQ/FREQNM
            FREQSQ = FREQ*FREQ
            XSELF = 4.0_dp*FREQSQ*HHW/PI/
     &        ((DEL*(FREQNM+FREQ))**2.0_dp + FREQSQ*4.*HHW**2.0_dp)
         ENDIF
C
         HTOTAL = HTOTAL + XSELF
C
C  End self-broadening section
C
      ENDIF
C
C  Case 3: Stark section
* Ignore if T<=3000K 2014-08-27
C
      IF (T.GE.3000.0_dp.AND.((ABS(DEL).GT.HFWID).OR.(NWID.EQ.3))) THEN
C
C  Stark wings due to ions and electrons.
C
C  Stark wings due to protons in the satellite regions of Lyman alpha,
C  beta and gamma, as well as Balmer alpha:
C
C  Assumed insensitive to T, and to scale linearly with N(H+).
C  We further assume N(H+) = N(e-).
C
C  There may be a contribution due to electrons also at long range. 
C  The static Holtsmark profile (1/2 the total STARK1) is used noting 
C  that the correction for quantum effects suggested by Stehle (1994, 
C  A&AS 104, 509 eqn 7) has been applied already.
C
         XSTARK0= STARK1(N,M,WAVE,WAVEH,T,XNE)
         IF (LYMANALF.AND.(WAVE.GT.(1.0d8/EHYD(M)))) THEN
C                    wave >  1215.6701\AA
            WAVERMP = 1325.0_dp
            XFACMAX = 5.00_dp
            IF      (WAVE.GT.5000.0_dp) THEN                              IGNORE
c                    wave > 1800.00\AA
               XSTARK = 0.0_dp
            ELSE IF (WAVE.GT.xlLyaH2p) THEN                               INTERP
c                    wave > 1290.30\AA
               ICUT = min0(NLyaH2p-2,max0(0,int(
     &                      (WAVE-xlLyaH2p)/dlLyaH2p))) + 1
               GRDWAVE      = xlLyaH2p + (ICUT-1)*dlLyaH2p
               XLyaH2p = (LyaH2p(ICUT+1)-LyaH2p(ICUT))
     &               /dlLyaH2p*(WAVE-GRDWAVE) + LyaH2p(ICUT)
               XLyaH2p = 10.0_dp**(XLyaH2p+XLNE+xLNLyaH2p)
             IF     (WAVE.LT.WAVERMP) THEN
               XSTARKb= 0.5_dp * STARK1(N,M,xlLyaH2p,WAVEH,T,XNE)
               SATb   = 10.0_dp**(LyaH2p(1)+XLNE+xLNLyaH2p)
c              if (WAVE.LT.1290.4_dp) print *,'XFAC: ',SATb/XSTARKb
               XFAC = 1.0_dp+(dmin1(SATb/XSTARKb,XFACMAX)*
     &                XSTARKb/SATb-1.0_dp)*(1.0_dp-cos(PI*
     &                (WAVE-WAVERMP)/(xlLyaH2p-WAVERMP)))*0.5_dp
               XLyaH2p= XFAC * XLyaH2p
             ENDIF
               XSTARK = XLyaH2p + 0.5_dp*XSTARK0
            ELSE                                                          SCALE
C                    wave < 1277.80\AA      > 78259/cm
C  If Lyman alpha we match the static ion part to the Allard et al 
C  value at the blue end of their data at 1290.30\AA in the red wing.
               XSTARKb = 0.5_dp * STARK1(N,M,xlLyaH2p,WAVEH,T,XNE)
               SATb = 10.0_dp**(LyaH2p(1)+XLNE+xLNLyaH2p)
               XFAC = 1.0_dp+(dmin1(SATb/XSTARKb,XFACMAX)-1.)*
     &                DABS(DELW)/(xlLyaH2p-1.0D8/EHYD(M))
               XSTARK = XSTARK0 * 0.5 * (1.0_dp + XFAC)
            ENDIF
C
C  Lyman beta satellites
         ELSE IF((N.EQ.1).AND.(M.EQ.3).AND.WAVE.GT.(1.0D8/EHYD(M)))THEN
            WAVERMP = 1052.5_dp
            XFACMAX = 2.39_dp
            IF      (WAVE.GT.1216.0_dp) THEN
               XSTARK = 0.0_dp
            ELSE IF (WAVE.GT.xlLybH2p) THEN
               ICUT = min0(NLybH2p-2,max0(0,int(
     &                      (WAVE-xlLybH2p)/dlLybH2p))) + 1
               GRDWAVE      = xlLybH2p + (ICUT-1)*dlLybH2p
               XLybH2p = (LybH2p(ICUT+1)-LybH2p(ICUT))
     &               /dlLybH2p*(WAVE-GRDWAVE) + LybH2p(ICUT)
               XLybH2p = 10.0_dp**(XLybH2p+XLNE+xLNLybH2p)
             IF     (WAVE.LT.WAVERMP) THEN
               XSTARKb= 0.5_dp * STARK1(N,M,xlLybH2p,WAVEH,T,XNE)
               SATb   = 10.0_dp**(LybH2p(1)+XLNE+xLNLybH2p)
c              if (WAVE.LT.1032.5_dp) print *,'XFAC: ',SATb/XSTARKb
               XFAC = 1.0_dp+(dmin1(SATb/XSTARKb,XFACMAX)*
     &                XSTARKb/SATb-1.0_dp)*
     &                (1.0_dp-cos(PI*(WAVE-WAVERMP)/
     &                (xlLybH2p-WAVERMP)))*0.5_dp
               XLybH2p= XFAC * XLybH2p
             ENDIF
               XSTARK = XLybH2p + 0.5_dp*XSTARK0
            ELSE
               XSTARKb = 0.5_dp * STARK1(N,M,xlLybH2p,WAVEH,T,XNE)
               SATb = 10.0_dp**(LybH2p(1)+XLNE+xLNLybH2p)
               XFAC = 1.0_dp+(dmin1(SATb/XSTARKb,XFACMAX)-1.0_dp)*
     &                DABS(DELW)/(xlLybH2p-1D8/EHYD(M))
               XSTARK = XSTARK0 * 0.5_dp * (1. + XFAC)
            ENDIF
C
C  Lyman gamma satellites
         ELSE IF ((N.EQ.1).AND.(M.EQ.4).AND.WAVE.GT.(1D8/EHYD(M))) THEN
            WAVERMP = 1002.0_dp
C  Limit to how much larger than XSTARK0 at xlLygH2p the satellites can be
            XFACMAX = 2.1_dp
            IF      (WAVE.GT.1026.0_dp) THEN
               XSTARK = 0.0_dp
            ELSE IF (WAVE.GT.xlLygH2p) THEN
               ICUT = min0(NLygH2p-2,max0(0,int(
     &                      (WAVE-xlLygH2p)/dlLygH2p))) + 1
               GRDWAVE      = xlLygH2p + (ICUT-1)*dlLygH2p
               XLygH2p = (LygH2p(ICUT+1)-LygH2p(ICUT))
     &               /dlLygH2p*(WAVE-GRDWAVE) + LygH2p(ICUT)
               XLygH2p = 10.0_dp**(XLygH2p+XLNE+xLNLygH2p)
             IF     (WAVE.LT.WAVERMP) THEN
               XSTARKb= 0.5_dp * STARK1(N,M,xlLygH2p,WAVEH,T,XNE)
               SATb   = 10.0_dp**(LygH2p(1)+XLNE+xLNLygH2p)
               XFAC   = 1.0_dp+(dmin1(SATb/XSTARKb,XFACMAX)*
     &                  XSTARKb/SATb-1.0_dp)*(1.0_dp-cos(PI*
     &                  (WAVE-WAVERMP)/(xlLygH2p-WAVERMP0)))*0.5_dp
c    &                            (WAVE-WAVERMP)/(xlLygH2p-WAVERMP)
               XLygH2p= XFAC * XLygH2p
             ENDIF
               XSTARK = XLygH2p + 0.5*XSTARK0
            ELSE
               XSTARKb= 0.5_dp * STARK1(N,M,xlLygH2p,WAVEH,T,XNE)
               SATb = 10.0_dp**(LygH2p(1)+XLNE+xLNLygH2p)
               XFAC = 1.0_dp+(dmin1(SATb/XSTARKb,XFACMAX)-1.)*
     &                DABS(DELW)/(xlLygH2p-1D8/EHYD(M))
               XSTARK = XSTARK0 * 0.5_dp * (1.0_dp + XFAC)
            ENDIF
C
C  Balmer alpha satellites
         ELSE IF ((N.EQ.2).AND.(M.EQ.3).AND.
     &            WAVE.GT.(1.0D8/(EHYD(M)-EHYD(N)))) THEN
            WAVERMP = 7275.0_dp
            XFACMAX = 2.34_dp
C  Beyond the satellites: lambda > 9915.\AA
            IF      (WAVE.GT.(xlHaH2p+dlHaH2p*(NHaH2p-1.))) THEN
               XSTARK = 0.0_dp
c                  \lambda > 6819.8\AA
            ELSE IF (WAVE.GT.xlHaH2p) THEN
               ICUT = min0(NHaH2p-2,max0(0,int(
     &                      (WAVE-xlHaH2p)/dlHaH2p))) + 1
               GRDWAVE      = xlHaH2p + (ICUT-1)*dlHaH2p
               XHaH2p = (HaH2p(ICUT+1)-HaH2p(ICUT))
     &               /dlHaH2p*(WAVE-GRDWAVE) + HaH2p(ICUT)
               XHaH2p = 10.0_dp**(XHaH2p+XLNE+xLNHaH2p)
c                  \lambda < 7275.\AA
             IF     (WAVE.LT.WAVERMP) THEN
               XSTARKb= 0.5_dp * STARK1(N,M,xlHaH2p,WAVEH,T,XNE)
               SATb   = 10.0_dp**(HaH2p(1)+XLNE+xLNHaH2p)
c              if (WAVE.LT.6820.5) print *,'XFAC: ',SATb/XSTARKb
               XFAC = 1.0_dp+(dmin1(SATb/XSTARKb,XFACMAX)*
     &                XSTARKb/SATb-1.0_dp)*(1.0_dp-cos(PI*
     &                (WAVE-WAVERMP)/(xlHaH2p-WAVERMP)))*0.5_dp
c    &                            (WAVE-WAVERMP)/(xlHaH2p-WAVERMP)
               XHaH2p = XFAC * XHaH2p
             ENDIF
               XSTARK = XHaH2p + 0.5_dp*XSTARK0
            ELSE
C
C  Ramp down from satellite to Stark wing over 80\AA from blue edge of
C  satellite data, xlHaH2p, using a differentiable cosine-bridge.
C  WARNING: This gives a rather prominent peak at log(Ne) >~ 19.
               XSTARKb= 0.5_dp * STARK1(N,M,xlHaH2p,WAVEH,T,XNE)
               SATb   = 10.0_dp**(HaH2p(1)+XLNE+xLNHaH2p)
               XFAC = 1.0_dp+(dmin1(SATb/XSTARKb,XFACMAX)-1.0_dp)*
     &                DABS(DELW)/(xlHaH2p-1D8/(EHYD(M)-EHYD(N)))
               XSTARK = XSTARK0 * 0.5_dp * (1.0_dp + XFAC)
            ENDIF
         ELSE
            XSTARK = XSTARK0
         ENDIF

         HTOTAL = HTOTAL + XSTARK
C
C  End Stark section
C
      END IF
C
      HPROFL = HTOTAL
C
      RETURN
      END

C***********************************************************************
      real*8 FUNCTION STARK1(N,M,WAVE,WAVEH,T,XNE)
C
C  Returns the Stark broadened line profile.  The code follows Griem's 
C  theories (mostly 1960 ApJ, 132, 883) with corrections to approximate 
C  the Vidal, Cooper & Smith (1973, ApJS 25, 37) profiles.
C
C  Area normalised to unity with frequency.
C
C  by Deane Peterson & Bob Kurucz.
C  (adapted, corrected and comments added by PB)
C
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer :: m, m1, mmn, mp, n, n1, np
      real(dp), PARAMETER :: K = 1.380649E-16  !Boltzmann in cgs
      real(dp), PARAMETER :: CLIGHT = 2.99792458E18
      real(dp), PARAMETER :: H = 6.6260715E-27  !Planck in cgs
      real(dp), PARAMETER :: PI = 3.14159265359, SQRTPI = 1.77245385
      real(dp) :: beta, c1, c1con, c1d, c2, c2con, c2d, dbeta, DEL
      real(dp) :: DELW, F, FO, fns, FREQ, FREQNM, g, g1, gam, gcon1
      real(dp) :: gcon2, gnm, he1frc, p1, pp, prqs, sofbet, stark1p
      real(dp) :: t, t3nhe, t4, t43, tp, vcse1f
      real(dp) :: WAVE, WAVEH, wavehp, wavep, wty1, xknm, xknmtb, xm
      real(dp) :: xm2, xm2mn2, xmn2, xn, xn2, xne, xne16,xnep, y1, y1b
      real(dp) :: y1num, y1s, y1scal, y1wht, y1wtm, y2
      DIMENSION Y1WTM(2,2),XKNMTB(4,3)
      LOGICAL LYMANALF
      SAVE
C
C  Knm constants as defined by Griem (1960, ApJ 132, 883) for the long 
C  range Holtsmark profile (due to ions only). Lyman and Balmer series 
C  are from VCS, higher series from elsewhere.
C
      DATA XKNMTB /
     &     0.0001716_dp, 0.0090190_dp, 0.1001000_dp, 0.5820000_dp,
     &     0.0005235_dp, 0.0177200_dp, 0.1710000_dp, 0.8660000_dp,
     &     0.0008912_dp, 0.0250700_dp, 0.2230000_dp, 1.0200000_dp /
C
      DATA Y1WTM / 1.0d18, 1.0d17, 1.0d16, 1.0d14 /
      DATA N1 / 0 /, M1 / 0 /
C
C
C  Variables depending on conditions
C

c we try to save time  BPz 12/10-2007
      if (n.eq.np .and. m.eq.mp .and. wave.eq.wavep .and. 
     &     waveh.eq.wavehp .and. t.eq.tp .and. xne.eq.xnep) then
         stark1=stark1p
         return
      endif

      T4 = T/10000.0_dp
      T43 = T4**0.3_dp
      T3NHE = T43*HE1FRC
      XNE16 = XNE**0.1666667_dp
      PP = XNE16*0.08989_dp/SQRT(T) ! the shielding parameter 
      FO = XNE16**4*1.25d-9      ! Holtsmark normal field strength
      Y1B = 2.0_dp/(1.0_dp+0.012_dp/T*SQRT(XNE/T))
      Y1S = T43/XNE16
      C1D = FO*78940.0_dp/ T
      C2D = FO**2/5.96d-23/XNE
      GCON1 = 0.2_dp+0.09_dp*SQRT(T4)/(1.0_dp+XNE/1.0d13)
      GCON2 = 0.2_dp/(1.0_dp+XNE/1.0d15)
C
      DELW = WAVE-WAVEH
      FREQNM = CLIGHT/WAVEH
      FREQ = CLIGHT/WAVE
      DEL = FREQ-FREQNM
C
C  Variables dependent on line - compute first time only
C
      IF((N.NE.N1).OR.(M.NE.M1)) THEN  
         N1 = N
         M1 = M
         MMN = M-N
         XN = N
         XN2 = XN*XN
         XM = M
         XM2 = XM*XM
         XMN2 = XM2*XN2
         XM2MN2 = XM2-XN2
         GNM = XM2MN2/XMN2
C
C  Knm constants not tabulated from approximate asymptotic expression 
C  (Griem 1960 eqn 33) where 1./(1.+.13/FLOAT(MMN)) appears to be a 
C  correction factor to match to the tables.
C
         IF ((MMN.LE.3).AND.(N.LE.4)) THEN
            XKNM = XKNMTB(N,MMN)
         ELSE
            XKNM = 5.5d-5/GNM*XMN2/(1.0_dp+0.13_dp/FLOAT(MMN))
         END IF
C
C  Some weighting factors which relate to y1, which is the velocity at 
C  which the minimum impact parameter (where second order perturbation 
C  theory breaks down) and the Lewis cutoff (the limit of validity of 
C  the impact approximation) are the same.
C
         IF(M.EQ.2) THEN
            Y1NUM = 550.0_dp
         ELSE IF (M.EQ.3) THEN
            Y1NUM = 380.0_dp
         ELSE
            Y1NUM = 320.0_dp
         END IF
         IF (MMN.LE.2 .AND. N.LE.2) THEN
            Y1WHT = Y1WTM(N,MMN)
         ELSE IF (MMN.LE.3) THEN
            Y1WHT = 1.0d14
         ELSE
            Y1WHT = 1.0d13
         END IF
C
         C1CON = XKNM/WAVEH*GNM*XM2MN2
         C2CON = (XKNM/WAVEH)**2
C
      ENDIF
C
C  Compute line profile
C
C  PRQS is the quasistatic ion contribution
C  FNS  is the quasistatic electron contribution rel to PRQS
C  F    is the impact electron contribution
C
C  First compute the width of the impact electron profile roughly Griem
C  (1967, ApJ 147, 1092) eqn for w.
C
      WTY1 = 1.0_dp/(1.0_dp+XNE/Y1WHT)
      Y1SCAL = Y1NUM*Y1S*WTY1+Y1B*(1.0_dp-WTY1)
      C1 = C1D*C1CON*Y1SCAL
      C2 = C2D*C2CON
      G1 = 6.77_dp*SQRT(C1)
      BETA = DABS(DELW)/FO/XKNM
      Y1 = C1*BETA
      Y2 = C2*BETA**2
      IF ((Y2.LE.1.0d-4).AND.(Y1.LE.1.0d-5)) THEN
         GAM = G1*dMAX1(0.0_dp,0.2114_dp+LOG(SQRT(C2)/C1))*
     &         (1.0_dp-GCON1-GCON2)
      ELSE
         GAM = G1*(0.5_dp*EXP(-MIN(80.0_dp,Y1))+VCSE1F(Y1)-0.5_dp*
     &         VCSE1F(Y2))*(1.0_dp-GCON1/(1.0_dp+(90.0_dp*Y1)**3)-
     &         GCON2/(1.0_dp+2000.0_dp*Y1))
         IF (GAM.LE.1.0d-20) GAM = 0.0_dp
      END IF
C
C  Compute individual quasistatic and impact profiles.
C
      PRQS = SOFBET(BETA,PP,N,M)
      IF (GAM.GT.0.) THEN
         F = GAM/PI/(GAM*GAM+BETA*BETA)
      ELSE
         F = 0.0D0
      ENDIF
C
C  Fraction of electrons which count as quasistatic. A fit to eqn 8 
C  (2nd term) of Griem (1967, ApJ 147, 1092).
C
      P1 = (0.9_dp*Y1)**2
      FNS = (P1+0.03_dp*SQRT(Y1))/(P1+1.0_dp)
C
C  DBETA (=dBeta/dfreq) changes the area normalisation. 
C  DSQRT(WAVE/WAVEH) corrects the long range part to dfreq**-5/2
C  asymptote, (see Stehle and Hutcheon 1999, A&AS 140, 93).
C
      DBETA = CLIGHT/FREQ/FREQ/XKNM/FO
      STARK1 = (PRQS*(1.0_dp+FNS)+F)*DBETA * DSQRT(WAVE/WAVEH)
C
C  The red wing is multiplied by the Boltzmann factor to roughly account
C  for quantum effects (Stehle 1994, A&AS 104, 509 eqn 7). Assume 
C  absorption case.  If emission do for DEL.GT.0.
C
      IF (DEL.LT.0.0d0) STARK1 = STARK1 * DEXP(-DABS(H*DEL)/K/T)
      
      np=n
      mp=m
      wavep=wave
      wavehp=waveh
      tp=t
      xnep=xne
      stark1p=stark1

      return
C
      END


C***********************************************************************
      real*8 FUNCTION HFNM(N,M)
C
C  HFNM calculates hydrogen oscillator strengths
C
C  From Kurucz codes.
C
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer :: m, mstr, n, nstr
      real(dp) :: fk, fkn, fnm, gca, ginf, wt, wtc
      real(dp) :: xm, xmn, xmn12, xn
      SAVE
      DATA NSTR / 0 /, MSTR / 0 /
C
      HFNM=0.0_dp
      IF(M.LE.N) RETURN
      IF(N.NE.NSTR) THEN
        XN=N
        GINF=0.2027_dp/XN**0.71_dp
        GCA=0.124_dp/XN
        FKN=XN*1.9603_dp
        WTC=0.45_dp-2.4_dp/XN**3*(XN-1.0_dp)
        NSTR=N
        XM=M
        XMN=M-N
        FK=FKN*(XM/(XMN*(XM+XN)))**3
        XMN12=XMN**1.2_dp
        WT=(XMN12-1.0_dp)/(XMN12+WTC)
        FNM=FK*(1.0_dp-WT*GINF-(0.222_dp+GCA/XM)*(1.0_dp-WT))
        MSTR=M
      ELSE IF(M.NE.MSTR) THEN
        XM=M
        XMN=M-N
        FK=FKN*(XM/(XMN*(XM+XN)))**3
        XMN12=XMN**1.2_dp
        WT=(XMN12-1.0_dp)/(XMN12+WTC)
        FNM=FK*(1.0_dp-WT*GINF-(0.222_dp+GCA/XM)*(1.0_dp-WT))
        MSTR=M
      END IF
      HFNM=FNM
C
      RETURN
      END



C***********************************************************************
      real*8 FUNCTION VCSE1F(X)
C
C  E1 function calculator for VCS approximation. It's rough, but 
C  arranged to be fast. X must be >=0.
C
C  From Kurucz codes.
C
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: sp = selected_real_kind(6, 37)
      real(dp) :: x
      VCSE1F = 0.0_dp
      IF(X.LE.0.0_dp) RETURN
      IF(X.LE.0.01_dp) THEN
        VCSE1F = -LOG(X)-0.577215_dp+X
      ELSE IF(X.LE.1.0_dp) THEN
        VCSE1F = -LOG(X)-0.57721566_dp+
     &           X*(0.99999193_dp+X*(-0.24991055_dp+
     &           X*(0.05519968_dp+X*(-0.00976004_dp+
     &           X*0.00107857_dp))))
      ELSE IF(X.LE.300.0_dp) THEN
        VCSE1F = (X*(X+2.334733_dp)+0.25062_dp)/(X*(X+3.330657_dp)+
     &           1.681534_dp)/X*dEXP(-X)
      END IF
C
      RETURN
      END


C***********************************************************************
      real*8 FUNCTION SOFBET(B,P,N,M)
C
C  Calculates S(BETA,P) for Hydrogen lines ie. the Holtsmark profile for
C  quasistatic charged particles.  The alpha and beta lines of the first
C  three series are explicitly included. All other cases use the H18 
C  profile. Profiles are normalised to full oscillator strength. Method 
C  is based on Griem (1960, ApJ 132, 883).
C
C  By Deane Peterson and Bob Kurucz.
C
C  STORAGE FOR CORRECTIONS (P,BETA,IND),(P,IND),(P,IND)
C
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer :: im, indx, ip, jm, jp, m, mmn, n
      real(dp) :: b, b2, BETA(15), C(5,7), C1(5), C2(5), C3(5), C4(5)
      real(dp) :: C5(5), C6(5), C7(5), cbm, cbp, cc, corr, D(5,7), D1(5)
      real(dp) :: D2(5), D3(5), D4(5), D5(5), D6(5), D7(5), dd, p
      real(dp) :: PP(5), pr1, pr2, PROB1(75), PROB2(75), PROB3(75)
      real(dp) :: PROB4(75), PROB5(75), PROB6(75), PROB7(75)
      real(dp) :: PROPBM(5,15,7), sb, wt, wtbm, wtbp, wtpm, wtpp
      EQUIVALENCE (PROPBM(1,1,1),PROB1(1)),(PROPBM(1,1,2),PROB2(1))
      EQUIVALENCE (PROPBM(1,1,3),PROB3(1)),(PROPBM(1,1,4),PROB4(1))
      EQUIVALENCE (PROPBM(1,1,5),PROB5(1)),(PROPBM(1,1,6),PROB6(1))
      EQUIVALENCE (PROPBM(1,1,7),PROB7(1))
      EQUIVALENCE (C(1,1),C1(1)),(C(1,2),C2(1)),(C(1,3),C3(1))
      EQUIVALENCE (C(1,4),C4(1)),(C(1,5),C5(1)),(C(1,6),C6(1))
      EQUIVALENCE (C(1,7),C7(1))
      EQUIVALENCE (D(1,1),D1(1)),(D(1,2),D2(1)),(D(1,3),D3(1))
      EQUIVALENCE (D(1,4),D4(1)),(D(1,5),D5(1)),(D(1,6),D6(1))
      EQUIVALENCE (D(1,7),D7(1))
      SAVE PROPBM,C,D,PP,BETA
C
C  Lyman alpha
C
      DATA PROB1/
     & -0.980_dp, -0.967_dp, -0.948_dp, -0.918_dp, -0.873_dp,
     & -0.968_dp, -0.949_dp, -0.921_dp, -0.879_dp, -0.821_dp,
     & -0.950_dp, -0.922_dp, -0.883_dp, -0.830_dp, -0.764_dp,
     & -0.922_dp, -0.881_dp, -0.830_dp, -0.770_dp, -0.706_dp,
     & -0.877_dp, -0.823_dp, -0.763_dp, -0.706_dp, -0.660_dp,
     & -0.806_dp, -0.741_dp, -0.682_dp, -0.640_dp, -0.625_dp,
     & -0.691_dp, -0.628_dp, -0.588_dp, -0.577_dp, -0.599_dp,
     & -0.511_dp, -0.482_dp, -0.484_dp, -0.514_dp, -0.568_dp,
     & -0.265_dp, -0.318_dp, -0.382_dp, -0.455_dp, -0.531_dp,
     & -0.013_dp, -0.167_dp, -0.292_dp, -0.394_dp, -0.478_dp,
     &  0.166_dp, -0.056_dp, -0.216_dp, -0.332_dp, -0.415_dp,
     &  0.251_dp,  0.035_dp, -0.122_dp, -0.237_dp, -0.320_dp,
     &  0.221_dp,  0.059_dp, -0.068_dp, -0.168_dp, -0.247_dp,
     &  0.160_dp,  0.055_dp, -0.037_dp, -0.118_dp, -0.189_dp,
     &  0.110_dp,  0.043_dp, -0.022_dp, -0.085_dp, -0.147_dp /
      DATA C1 /-18.396_dp, 84.674_dp, -96.273_dp, 3.927_dp,  55.191_dp /
      DATA D1 / 11.801_dp, 9.079_dp, -0.651_dp, -11.071_dp, -26.545_dp /
C
C  Lyman beta
C
      DATA PROB2/
     & -0.242_dp,  0.060_dp,  0.379_dp,  0.671_dp,  0.894_dp,
     &  0.022_dp,  0.314_dp,  0.569_dp,  0.746_dp,  0.818_dp,
     &  0.273_dp,  0.473_dp,  0.605_dp,  0.651_dp,  0.607_dp,
     &  0.432_dp,  0.484_dp,  0.489_dp,  0.442_dp,  0.343_dp,
     &  0.434_dp,  0.366_dp,  0.294_dp,  0.204_dp,  0.091_dp,
     &  0.304_dp,  0.184_dp,  0.079_dp, -0.025_dp, -0.135_dp,
     &  0.167_dp,  0.035_dp, -0.082_dp, -0.189_dp, -0.290_dp,
     &  0.085_dp, -0.061_dp, -0.183_dp, -0.287_dp, -0.374_dp,
     &  0.032_dp, -0.127_dp, -0.249_dp, -0.344_dp, -0.418_dp,
     & -0.024_dp, -0.167_dp, -0.275_dp, -0.357_dp, -0.420_dp,
     & -0.061_dp, -0.170_dp, -0.257_dp, -0.327_dp, -0.384_dp,
     & -0.047_dp, -0.124_dp, -0.192_dp, -0.252_dp, -0.306_dp,
     & -0.043_dp, -0.092_dp, -0.142_dp, -0.190_dp, -0.238_dp,
     & -0.038_dp, -0.070_dp, -0.107_dp, -0.146_dp, -0.187_dp,
     & -0.030_dp, -0.049_dp, -0.075_dp, -0.106_dp, -0.140_dp /
      DATA C2 / 95.740_dp,  18.489_dp, 14.902_dp, 24.466_dp, 42.456_dp /
      DATA D2 / -6.665_dp, -7.136_dp, -10.605_dp,-15.882_dp,-23.632_dp /
C
C  Balmer alpha
C
      DATA PROB3/
     & -0.484_dp, -0.336_dp, -0.206_dp, -0.111_dp, -0.058_dp,
     & -0.364_dp, -0.264_dp, -0.192_dp, -0.154_dp, -0.144_dp,
     & -0.299_dp, -0.268_dp, -0.250_dp, -0.244_dp, -0.246_dp,
     & -0.319_dp, -0.333_dp, -0.337_dp, -0.336_dp, -0.337_dp,
     & -0.397_dp, -0.414_dp, -0.415_dp, -0.413_dp, -0.420_dp,
     & -0.456_dp, -0.455_dp, -0.451_dp, -0.456_dp, -0.478_dp,
     & -0.446_dp, -0.441_dp, -0.446_dp, -0.469_dp, -0.512_dp,
     & -0.358_dp, -0.381_dp, -0.415_dp, -0.463_dp, -0.522_dp,
     & -0.214_dp, -0.288_dp, -0.360_dp, -0.432_dp, -0.503_dp,
     & -0.063_dp, -0.196_dp, -0.304_dp, -0.394_dp, -0.468_dp,
     &  0.063_dp, -0.108_dp, -0.237_dp, -0.334_dp, -0.409_dp,
     &  0.151_dp, -0.019_dp, -0.148_dp, -0.245_dp, -0.319_dp,
     &  0.149_dp,  0.016_dp, -0.091_dp, -0.177_dp, -0.246_dp,
     &  0.115_dp,  0.023_dp, -0.056_dp, -0.126_dp, -0.189_dp,
     &  0.078_dp,  0.021_dp, -0.036_dp, -0.091_dp, -0.145_dp /
      DATA C3 / -25.088_dp, 145.882_dp,-50.165_dp, 7.902_dp, 51.003_dp /
      DATA D3 / 7.872_dp,  5.592_dp, -2.716_dp, -12.180_dp, -25.661_dp /
C
C  Balmer beta
C
      DATA PROB4/
     & -0.082_dp,  0.163_dp,  0.417_dp,  0.649_dp,  0.829_dp,
     &  0.096_dp,  0.316_dp,  0.515_dp,  0.660_dp,  0.729_dp,
     &  0.242_dp,  0.393_dp,  0.505_dp,  0.556_dp,  0.534_dp,
     &  0.320_dp,  0.373_dp,  0.394_dp,  0.369_dp,  0.290_dp,
     &  0.308_dp,  0.274_dp,  0.226_dp,  0.152_dp,  0.048_dp,
     &  0.232_dp,  0.141_dp,  0.052_dp, -0.046_dp, -0.154_dp,
     &  0.148_dp,  0.020_dp, -0.094_dp, -0.200_dp, -0.299_dp,
     &  0.083_dp, -0.070_dp, -0.195_dp, -0.299_dp, -0.385_dp,
     &  0.031_dp, -0.130_dp, -0.253_dp, -0.348_dp, -0.422_dp,
     & -0.023_dp, -0.167_dp, -0.276_dp, -0.359_dp, -0.423_dp,
     & -0.053_dp, -0.165_dp, -0.254_dp, -0.326_dp, -0.384_dp,
     & -0.038_dp, -0.119_dp, -0.190_dp, -0.251_dp, -0.306_dp,
     & -0.034_dp, -0.088_dp, -0.140_dp, -0.190_dp, -0.239_dp,
     & -0.032_dp, -0.066_dp, -0.103_dp, -0.144_dp, -0.186_dp,
     & -0.027_dp, -0.048_dp, -0.075_dp, -0.106_dp, -0.142_dp /
      DATA C4 / 93.783_dp, 10.066_dp,  9.224_dp, 20.685_dp, 40.136_dp /
      DATA D4 / -5.918_dp, -6.501_dp,-10.130_dp,-15.588_dp,-23.570_dp /
C
C  Paschen alpha
C
      DATA PROB5/
     & -0.819_dp, -0.759_dp, -0.689_dp, -0.612_dp, -0.529_dp,
     & -0.770_dp, -0.707_dp, -0.638_dp, -0.567_dp, -0.498_dp,
     & -0.721_dp, -0.659_dp, -0.595_dp, -0.537_dp, -0.488_dp,
     & -0.671_dp, -0.617_dp, -0.566_dp, -0.524_dp, -0.497_dp,
     & -0.622_dp, -0.582_dp, -0.547_dp, -0.523_dp, -0.516_dp,
     & -0.570_dp, -0.545_dp, -0.526_dp, -0.521_dp, -0.537_dp,
     & -0.503_dp, -0.495_dp, -0.496_dp, -0.514_dp, -0.551_dp,
     & -0.397_dp, -0.418_dp, -0.448_dp, -0.492_dp, -0.547_dp,
     & -0.246_dp, -0.315_dp, -0.384_dp, -0.453_dp, -0.522_dp,
     & -0.080_dp, -0.210_dp, -0.316_dp, -0.406_dp, -0.481_dp,
     &  0.068_dp, -0.107_dp, -0.239_dp, -0.340_dp, -0.418_dp, 
     &  0.177_dp, -0.006_dp, -0.143_dp, -0.246_dp, -0.324_dp,
     &  0.184_dp,  0.035_dp, -0.082_dp, -0.174_dp, -0.249_dp,
     &  0.146_dp,  0.042_dp, -0.046_dp, -0.123_dp, -0.190_dp,
     &  0.103_dp,  0.036_dp, -0.027_dp, -0.088_dp, -0.146_dp /
      DATA C5 / -19.819_dp, 94.981_dp, -79.606_dp, 3.159_dp, 52.106_dp /
      DATA D5 / 10.938_dp, 8.028_dp, -1.267_dp, -11.375_dp, -26.047_dp /
C
C  Paschen beta
C
      DATA PROB6/
     & -0.073_dp,  0.169_dp,  0.415_dp,  0.636_dp,  0.809_dp,  0.102_dp,
     &  0.311_dp,  0.499_dp,  0.639_dp,  0.710_dp,  0.232_dp,  0.372_dp,
     &  0.479_dp,  0.531_dp,  0.514_dp,  0.294_dp,  0.349_dp,  0.374_dp,
     &  0.354_dp,  0.279_dp,  0.278_dp,  0.253_dp,  0.212_dp,  0.142_dp,
     &  0.040_dp,  0.215_dp,  0.130_dp,  0.044_dp, -0.051_dp, -0.158_dp,
     &  0.141_dp,  0.015_dp, -0.097_dp, -0.202_dp, -0.300_dp,  0.080_dp,
     & -0.072_dp, -0.196_dp, -0.299_dp, -0.385_dp,  0.029_dp, -0.130_dp,
     & -0.252_dp, -0.347_dp, -0.421_dp, -0.022_dp, -0.166_dp, -0.275_dp,
     & -0.359_dp, -0.423_dp, -0.050_dp, -0.164_dp, -0.253_dp, -0.325_dp,
     & -0.384_dp, -0.035_dp, -0.118_dp, -0.189_dp, -0.252_dp, -0.306_dp,
     & -0.032_dp, -0.087_dp, -0.139_dp, -0.190_dp, -0.240_dp, -0.029_dp,
     & -0.064_dp, -0.102_dp, -0.143_dp, -0.185_dp, -0.025_dp, -0.046_dp,
     & -0.074_dp, -0.106_dp, -0.142_dp /
      DATA C6 / 111.107_dp, 11.910_dp, 9.857_dp, 21.371_dp, 41.006_dp /
      DATA D6 / -5.899_dp, -6.381_dp,-10.044_dp,-15.574_dp,-23.644_dp /
C
C  Balmer 18
C
      DATA PROB7/
     &  0.005_dp,  0.128_dp,  0.260_dp,  0.389_dp,  0.504_dp,
     &  0.004_dp,  0.109_dp,  0.220_dp,  0.318_dp,  0.389_dp,
     & -0.007_dp,  0.079_dp,  0.162_dp,  0.222_dp,  0.244_dp,
     & -0.018_dp,  0.041_dp,  0.089_dp,  0.106_dp,  0.080_dp,
     & -0.026_dp, -0.003_dp,  0.003_dp, -0.023_dp, -0.086_dp,
     & -0.025_dp, -0.048_dp, -0.087_dp, -0.148_dp, -0.234_dp,
     & -0.008_dp, -0.085_dp, -0.165_dp, -0.251_dp, -0.343_dp,
     &  0.018_dp, -0.111_dp, -0.223_dp, -0.321_dp, -0.407_dp,
     &  0.032_dp, -0.130_dp, -0.255_dp, -0.354_dp, -0.431_dp,
     &  0.014_dp, -0.148_dp, -0.269_dp, -0.359_dp, -0.427_dp,
     & -0.005_dp, -0.140_dp, -0.243_dp, -0.323_dp, -0.386_dp,
     &  0.005_dp, -0.095_dp, -0.178_dp, -0.248_dp, -0.307_dp,
     & -0.002_dp, -0.068_dp, -0.129_dp, -0.187_dp, -0.241_dp,
     & -0.007_dp, -0.049_dp, -0.094_dp, -0.139_dp, -0.186_dp,
     & -0.010_dp, -0.036_dp, -0.067_dp, -0.103_dp, -0.143_dp /
      DATA C7 / 511.318_dp, 1.532_dp,  4.044_dp,  19.266_dp, 41.812_dp /
      DATA D7 / -6.070_dp, -4.528_dp, -8.759_dp, -14.984_dp,-23.956_dp /
      DATA PP / 0.0_dp, 0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp /
      DATA BETA /
     & 1.000_dp, 1.259_dp, 1.585_dp, 1.995_dp, 2.512_dp, 3.162_dp,
     & 3.981_dp, 5.012_dp, 6.310_dp, 7.943_dp, 10.00_dp, 12.59_dp,
     & 15.85_dp, 19.95_dp, 25.12_dp /
C
      IF(B.GT.500.) THEN
C
C  Very large B
C
        B2=B*B
        SOFBET=(1.5_dp/SQRT(B)+27.0_dp/B2)/B2
        RETURN
      END IF
C
C  Other cases
C
      CORR=1.
      B2=B*B
      SB=SQRT(B)
      INDX=7
      MMN=M-N
      IF(N.LE.3.AND.MMN.LE.2) INDX=2*(N-1)+MMN
C
C  Determine relevant Debye range
C
      IM=MIN(INT(5.0_dp*P)+1,4)
      IP=IM+1
      WTPP=50.0_dp*(P-PP(IM))
      WTPM=1.0_dp-WTPP
      IF(B.LE.25.12_dp) THEN
C
        JP=2
   1    IF(B.GT.BETA(JP) .AND. JP.LT.15) THEN
          JP=JP+1
          GO TO 1
        END IF
        JM=JP-1
C
        WTBP=(B-BETA(JM))/(BETA(JP)-BETA(JM))
        WTBM=1.0_dp-WTBP
        CBP=PROPBM(IP,JP,INDX)*WTPP+PROPBM(IM,JP,INDX)*WTPM
        CBM=PROPBM(IP,JM,INDX)*WTPP+PROPBM(IM,JM,INDX)*WTPM
        CORR=1.0_dp+CBP*WTBP+CBM*WTBM
C
C  Get approximate profile for the inner part
C
        PR1=0.0_dp
        PR2=0.0_dp
        WT=dMAX1(MIN(0.5_dp*(10.0_dp-B),1.0_dp),0.0_dp)
        IF(B.LE.10.0_dp) PR1=8.0_dp/(83.0_dp+(2.0_dp+0.95_dp*B2)*B)
        IF(B.GE.8.0_dp)  PR2=(1.5_dp/SB+27.0_dp/B2)/B2
        SOFBET=(PR1*WT+PR2*(1.0_dp-WT))*CORR
      ELSE
C
C  Asymptotic part for medium B's
C
        CC=C(IP,INDX)*WTPP+C(IM,INDX)*WTPM
        DD=D(IP,INDX)*WTPP+D(IM,INDX)*WTPM
        CORR=1.0_dp+DD/(CC+B*SB)
        SOFBET=(1.5_dp/SB+27.0_dp/B2)/B2*CORR
      END IF
C
      RETURN
      END


C***********************************************************************
      REAL*8 FUNCTION AIRVAC(W)
C
C  W is air wavelength in Angstroms. WAVEN is air wavenumber which is 
C  usually good enough, must iterate for exact solution
C
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: sp = selected_real_kind(6, 37)
      REAL(dp) :: W, WAVEN, WNEW
      WAVEN = 1.0D8 / W
      WNEW = W * (1.0000834213D0 + 2406030.0D0 / (1.30D10 - WAVEN**2) +
     &       15997.0D0 / (3.89D9 - WAVEN**2))
      WAVEN = 1.0D8 / WNEW
      WNEW = W * (1.0000834213D0 + 2406030.0D0 / (1.30D10 - WAVEN**2) +
     &       15997.0D0 / (3.89D9 - WAVEN**2))
      WAVEN = 1.0D8 / WNEW
      AIRVAC = W*(1.0000834213D0 + 2406030.0D0 / (1.30D10 - WAVEN**2) +
     &         15997.0D0 / (3.89D9 - WAVEN**2))
C
      RETURN
      END


C***********************************************************************
      REAL*8 FUNCTION VACAIR(W)
C
C  W is vacuum wavelength in Angstroms
C
      implicit none
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: sp = selected_real_kind(6, 37)
      REAL(dp) :: W, WAVEN
      WAVEN = 1.0D8 / W
      VACAIR = W/(1.0000834213D0 + 2406030.0D0 / (1.30D10 - WAVEN**2) +
     &       15997.0D0 / (3.89D9 - WAVEN**2))
C
      RETURN
      END


C***********************************************************************
*
