      program tapet
*
* print the wallpaper for a sos or pos marcs model
* v2. Includes output formatted for multi  BPz 01/02-2019
*
      character*256 file,modnam
      character*256 multimodfile,multitaufile

10    print*, 'what model do you want to print?'
      read(5,'(a)') file
      open (1,file=file,status='old',form='unformatted',
     & err=10,convert='big_endian',recl=1848001)

* this modnam is not the modname from inside the model file.
* it is just the file name

      modnam=file

      print*,' output file for listing?'
      read(5,'(a)') file
      open(20,file=file,status='unknown')
***************
* tau file contains:
*  a depth-scale is read from file dscale.
*  the depth scale can be given either as log column mass or
*  as log tau5000. the unit is given by the first letter of
*  the depth scale type string
*  and should of course be the same in the atmos and the dscale
*  files. the depth scale can either be given explicitly or as
*  the lg depth value in the first and last depth point. this is
*  indicated by the sign of ndep in the dscale file.
*  the file dscale contains :
*
*  depth scale identification line                        format a
*  depth scale type, first non-blank character t or m,
*  t for tau5000 scale, m for mass scale                  format a
*  ndep, dpcon1                                           free format
*
*  the rest of the file depends on the sign of ndep.
*  ndep.lt.0:
*
*  lg depth(1), lg depth(ndep)                            free format
*  the depth points are distributed evenly in lg depth
*
*  ndep.gt.0:
*
*  (lg depth(k), k=1,ndep)                                free format
*
*  dptype is set equal to the first non-blank letter of
*  the depth scale type string.
*  the depth scale unit is determined from the value of dptype
*
* dptype
*  m       lg depth is lg column mass, dpcon1= lg tau5000(1)
*  t       lg depth is lg tau5000,     dpcon1= lg column mass(1)
*  h       lg depth is height (km)     dpcon1= lg tau5000(1)
*
*  if dpcon1 is 0 a start value is calculated in  routine  dpconv
*******************
* mod file contains:
*  atmospheric identification line, (72 characters)       format a
*  depth scale type, starting with t or m,
*  t for tau5000 scale, m for mass scale                  format a
*  log gravity                                            free format
*  number of depth points (ndep),                         free format
*  (log depth(k), t(k), ne(k), v(k), vturb(k),k=1,ndep),  free format
*  ((nhi(i,k),i=1,5),nhii(k),k=1,ndep)                    free format
*
*  v and vturb in km per second, else cgs units
*  if there are no hydrogen populations in the input file
*  they are set to the lte values

      i=index(modnam,'/',back=.true.)
      print*,modnam(i+1:256)
      modnam=modnam(i+1:256)
      write(multimodfile,101) 'atmos.',trim(modnam)
101   format(a,a)
      write(multitaufile,101) 'dscale.',trim(modnam)
      open(30,file=multimodfile,status='unknown')
      open(40,file=multitaufile,status='unknown')

      i=index(modnam,'t')
      print*,'index',i
      print*,modnam(i+1:i+2)
      read(modnam(i+1:i+2),'(i2)',err=99) ixit
      goto 98
99    read(modnam(i+1:i+1),'(i)') ixit
98    continue
      xit=float(ixit)
      print*,'Microturbulence = ',xit,' km/s'
      write(30,30) trim(modnam)
      write(40,30) trim(modnam)
30    format('MARCS ',a)
      write(30,31) 
      write(40,31) 
31    format('TAU(5000) SCALE')

      call listmo(1,1,xit)
      end



      SUBROUTINE LISTMO(MO,IARCH,xit)
C        THIS ROUTINE PRINTS A NUMBER (MO) OF MODELS, WHICH ARE STORED ON
C        FORTRAN UNIT IARCH.
C
      include 'parameter.inc'
      include 'lambda.par'
C
      logical varimacro
      character*4 idrab1,idrab2,idrab3,idrab4,idrab5,idrab6
      character*128 modnam
      DIMENSION ABUND(16),TKORRM(NDP),FCORR(NDP),TAU(NDP),TAUS(NDP),
     *T(NDP),PE(NDP),PG(NDP),PRAD(NDP),PTURB(NDP),XKAPR(NDP),RO(NDP),
     *CP(NDP),CV(NDP),AGRAD(NDP),Q(NDP),U(NDP),V(NDP),ANCONV(NDP),
     *PRESMO(30,NDP),FCONV(NDP),RR(NDP),Z(NDP),EMU(NDP),HNIC(NDP),
     *NJ(16),XLR(20+1),IEL(16),ANJON(16,5),PART(16,5),PROV(mkomp,20+1),
     *ABSKA(mkomp),SPRIDA(mkomp),XLB(lpoint),FLUXME(lpoint),
     *  FLUMAG(lpoint),fluxmec(lpoint),
     &  PEP(16),ABNAME(mkomp),SOURCE(mkomp),DUMMY(11*NDP+2)
      DIMENSION W(lpoint),UW(12),BW(21),VW(25),RW(43),XIW(51),XJW(59),
     &  XKW(39)
      doubleprecision drr(ndp)
      CHARACTER*10 DAG,NAME,NAMEP,KLOCK
      CHARACTER*8 ABNAME,SOURCE
      DIMENSION WAVFLX(10),ptio(ndp)
      COMMON /UTPUT/IREAD,IWRIT
      COMMON /MASSE/RELM,RELLUM,RELATM,RELRAD,BOLMAG
      real    abSc,abTi,abV,abMn,abCo
      real inputabund(92)
      character*2 inputaname(92)
      real macrovel(ndp), macrobeta
      integer modtype
*
      common/specnew9/nlp,xlr,taus
* added for osplot *************
      real app(11*ndp),ptot(ndp)
      common /binstruc/bPPR(NDP),bPPT(NDP),bPP(NDP),bGG(NDP),
     & bZZ(NDP),bDD(NDP),
     & bVV(NDP),bFFC(NDP),bPPE(NDP),bTT(NDP),
     & bTAULN(NDP),NbTAU,IbTER,erad(ndp)
      common /struct/ teff,tau,t,z,ptot,prad,pg,pturb,pe,ro,rr,
     &                xkapr
      common /spectr/ nlb,xlb,w,fluxme
      common /spectrc/fluxmec
      common /pressure/ presmo,ptio
* *************
      DATA UW/0.145,0.436,0.910,1.385,1.843,2.126,2.305,2.241,1.270,
     *0.360,0.128,0.028/
      data BW/0.003,0.026,0.179,0.612,1.903,2.615,2.912,
     *3.005,2.990,2.876,2.681,2.388,2.058,1.725,1.416,1.135,0.840,0.568,
     *0.318,0.126,0.019/
      data VW/0.006,0.077,0.434,1.455,2.207,2.703,2.872,
     *2.738,2.505,2.219,1.890,1.567,1.233,0.918,0.680,0.474,0.312,0.200,
     *0.132,0.096,0.069,0.053,0.037,0.022,0.012/
      DATA  RW/0.03,0.06,0.17,0.28,0.39,0.50,0.60,0.69,0.74,0.79,
     &         0.84,0.88,0.91,0.94,0.96,0.98,0.99,1.00,0.97,0.94,
     &         0.90,0.85,0.79,0.73,0.65,0.57,0.50,0.42,0.36,0.31,
     &         0.24,0.17,0.14,0.11,0.08,0.06,0.05,0.04,0.03,0.02,
     &         0.02,0.01,0.01/
      DATA XIW/0.01,0.01,0.09,0.17,0.26,0.36,0.46,0.56,0.66,0.76,
     &         0.86,0.96,0.97,0.98,0.98,0.99,0.99,1.00,0.99,0.98,
     &         0.96,0.93,0.88,0.84,0.78,0.71,0.64,0.58,0.52,0.47,
     &         0.42,0.36,0.32,0.28,0.24,0.20,0.18,0.15,0.12,0.10,
     &         0.09,0.08,0.06,0.05,0.04,0.03,0.03,0.02,0.02,0.01,
     &         0.01/
      DATA XJW/0.01,0.02,0.02,0.03,0.04,0.06,0.11,0.16,0.26,0.35,
     &         0.49,0.62,0.78,0.93,0.89,0.85,0.82,0.78,0.78,0.78,
     &         0.79,0.80,0.82,0.85,0.89,0.93,0.84,0.75,0.70,0.64,
     &         0.64,0.63,0.63,0.63,0.64,0.66,0.67,0.68,0.69,0.70,
     &         0.70,0.70,0.68,0.66,0.63,0.60,0.53,0.46,0.37,0.27,
     &         0.20,0.14,0.12,0.09,0.08,0.06,0.04,0.02,0.01/
      DATA XKW/0.01,0.035,0.17,0.30,0.46,0.65,0.90,0.97,0.995,
     &         1.00,0.995,0.985,0.975,0.965,0.96,0.955,0.95,0.955,
     &         0.965,0.97,0.97,0.965,0.955,0.945,0.94,0.945,0.95,
     &         0.95,0.955,0.95,0.93,0.88,0.785,0.62,0.46,0.27,0.13,
     &         0.05,0.01/
* Ref HL Johnson ApJ 141,923  1965.
*
      DATA NAME/'LOCAL'/,NAMEP/'PARSONS'/
      DATA A,B/.34785485,.65214515/

*      save a,b,name,namep

      IREAD=5
      IWRIT=20
C
      REWIND IARCH
***** added for osplot ***************

        do i=ndp,1,-1
cc          print*,'oldsta; attempting to read binary model with i=',i
          rewind(iarch)
          read(iarch,err=499)
     &        bppr(1:i),bppt(1:i),bpp(1:i),bgg(1:i),bzz(1:i),bdd(1:i),
     &        bvv(1:i),bffc(1:i),bppe(1:i),btt(1:i),btauln(1:i),nbtau,
     &        ibter,erad(1:i)

          if (nbtau.gt.0.and.nbtau.le.i) then
* we should have found the right dimension
* But we have lost the compatibility with old marcs files that
* had no erad variable! BPz 05/07-2000
            print*,'readmo; binary model had ndp= ',i
            print*,'readmo;       we use now ndp= ',ndp
            print*,"readmo: prev model had ",nbtau," depthpoints"
            goto 497
          endif
499       continue
          if (i.eq.1) then
            print*,'ERROR! We have tried to read the input binary model'
            print*,'ERROR! with all possible values of ndp between 1 '
            print*,'ERROR! and ',ndp,'. And nothing works!'
            print*,'ERROR! Length of common block in stored model seems'
            print*,'ERROR! too long compared to NDP.'
            stop
          endif
        enddo
497     continue
*
c      do k=1,nbtau
c        print*,k,bppr(k),bppt(k),bpp(k),bgg(k),bzz(k)
c        print*,k,bdd(k),bvv(k),bffc(k),bppe(k),btt(k)
c        print*,k,btauln(k),erad(k)
c      enddo
      rewind(iarch)
***********************


      DO 1 IMO=1,MO
* first extract the last few records to get the full chemical composition
*  and the long name. BPz 11/01-2008
      do i=1,2
        read(iarch)
      enddo
      read(iarch) 
     &     teff,flux,g,palfa,pny,py,pbeta,iline,istral,mihal,
     &            idrab1,idrab2,idrab3,idrab4,idrab5,idrab6,
     &            itmax,nel
      do i=1,2
        read(iarch) jtau
      enddo
      ntpo=0
      do k=1,jtau
        read(iarch) kr,tau(k)
        tauk=alog10(tau(k))+10.01
        ktau=tauk
        if(abs(tauk-ktau).gt.0.02) go to 131
        if(ktau.eq.10) k0=k
        ntpo=ntpo+1
  131   continue
      enddo

      read(iarch)(nj(i),i=1,nel),nlp
      do k=1,ntpo
        do j=1,nlp
          if (j.eq.1) then
            do i=1,nel
              read(iarch)
            enddo
          endif
          read(iarch)
        enddo
      enddo
*****************
c  there are 
c 1) old models with nothing after flux record (very old)
c 2) models with longname after flux record (older BPz)
c 3) models with weight record after flux record, followed by long name
c    (Bpz on bohor)
c 4) the v7.3 models with flux record, fluxmec record, abundances record,
c    longname record (BPz v7.3 new models)
c 5) KE models with flux record, fluxmec record, abundance record
c
c flux record
      modtype=0
      read(iarch,err=192) nlb,(xlb(j),fluxme(j),j=1,nlb),(w(j),j=1,nlb)
      read(iarch,end=190,err=191) (fluxmec(j),j=1,nlb)
      read(iarch,err=193,end=193) inputabund,inputaname
      read(iarch,end=199) modnam
      modtype=4
      print*,'Model type 4'
      goto 199
190   modtype=1
      print*,'Model type 1'
      goto 199
191   backspace(iarch)
      read(iarch) modnam
      modtype=2
      print*,'Model type 2'
      goto 199
192   backspace(iarch)
      read(iarch) modnam
      modtype=3
      print*,'Model type 3'
      goto 199
193   modtype=5
      print*,'Model type 5'
      backspace(iarch) 
      read(iarch) inputabund
      do j=1,92
        inputaname(j)=' '
      enddo
      goto 199
* we end up here with the right value for modtype
199   continue
      if (modtype.eq.0) then
        print*,'************'
        print*,' could not find model type!'
        print*,'************'
        stop 
      endif
      rewind iarch
* Now do the full reading:
      READ(IARCH) 
      print*,'******** just before inord record ***********'
      READ(IARCH) INORD,DAG,KLOCK
      print*,'******** just before abund record ***********'
      READ(IARCH,err=999) 
     &     TEFF,FLUX,G,PALFA,PNY,PY,PBETA,ILINE,ISTRAL,MIHAL,
     &            IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &            ITMAX,NEL,(ABUND(I),I=1,NEL),abSc,abTi,abV,abMn,abCo,
     &            macrobeta,macrovel  
      goto 998
999   continue
      backspace(iarch)
      print*,'******** just before abund record again (no macro)****'
      READ(IARCH) 
     &     TEFF,FLUX,G,PALFA,PNY,PY,PBETA,ILINE,ISTRAL,MIHAL,
     &            IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &            ITMAX,NEL,(ABUND(I),I=1,NEL),abSc,abTi,abV,abMn,abCo
      macrobeta=0.
      do i=1,ndp
        macrovel(i)=0.0
      enddo
998   continue
      if (modtype.eq.1.or.modtype.eq.5) then
        write(modnam,195) idrab1,idrab2,idrab3,idrab4,idrab5,idrab6
195     format(6a4)
      endif
* maybe useful for some of KE's models
      print*,' trying double precision for rr(1:ndp)'
      READ(IARCH,err=1950,end=1950) 
     &       JTAU,NCORE,DIFLOG,TAUM,RADIUS,(dRR(K),K=1,JTAU)
      do k=1,jtau
        rr(k)=sngl(drr(k))
      enddo
      goto 1951
1950  continue
      backspace(iarch)
      READ(IARCH) 
     &       JTAU,NCORE,DIFLOG,TAUM,RADIUS,(RR(K),K=1,JTAU)
1951  continue
      print*,JTAU,NCORE,DIFLOG,TAUM,RADIUS,(RR(K),K=1,JTAU)
      backspace(iarch)
      write(20,219)
      GLOG=ALOG10(G)
      write(20,300) TEFF,GLOG,modnam(1:(index(modnam,' ')))
      write(20,283)
      FNORD=0.1*INORD
c  fnord=3.0 = apollo version; =4.0 = Sun version
cccc      fnord=3.0
      relrad=log10(radius/6.9598E10)
      if (radius.le.2.) then 
        write(20,281) FNORD,DAG,KLOCK
        relm=0.
        relatm=0.
        rellum=0.
        bolMAG=99.
      else
        write(20,282) FNORD,DAG,KLOCK
**************
        print*,'radius, rr(2) ',radius,rr(2)
        relm=10**(glog-4.44+2.*relrad)
CCCCC      RELATM=-3.+ALOG10(TSUN/TEFF/RELM)+0.5*RELLUM
        relatm=log10((rr(2)-radius)/radius)
        RELLUM=4.*ALOG10(TEFF/5780.)+2.*RELRAD
        BOLMAG=4.74-2.5*RELLUM
**************
      endif
      write(20,284)
      write(20,200)
C        CONVERT TO 'PHYSICAL FLUX'
      FLUX=3.14159*FLUX
      write(20,201) TEFF,FLUX,G
      if(radius.le.2.) then
        write(20,2001)
      else
        write(20,2002) 10.**relrad,radius,relm,10.**relatm,10.**rellum,
     ,                 bolmag
      endif
      write(20,2003)PALFA,PNY,PY
      IF(PBETA.LE.0. .and. macrobeta.le.0.) write(20,256)
      IF(PBETA.GT.0.0) write(20,257) PBETA
      IF(macrobeta.GT.0.0) then
        varimacro=.false.
        do i=1,jtau
          if (macrovel(i).ne.macrovel(1)) varimacro=.true.
        enddo
        if (varimacro) then
          write(20,259) macrobeta,macrovel(1)/1.e5
        else
          write(20,258) macrobeta,macrovel(1)/1.e5
        endif
      endif
      IF(ISTRAL.LT.1) write(20,250)
      IF(ISTRAL.GE.1) write(20,251)
      IF(ILINE.LT.1) write(20,252)
      IF(ILINE.GE.1) write(20,253)
      write(20,254) MIHAL,modnam(1:index(modnam,' '))
      DO I=1,NEL
        ABUND(I)=ALOG10(ABUND(I))+12.
      enddo
      if (modtype.eq.1.or.modtype.eq.2.or.modtype.eq.3) then
        write(20,202)
        write(20,203) (ABUND(I),I=1,NEL)
        write(20,2021)
        abSc=alog10(abSc)+12.
        abTi=alog10(abTi)+12.
        abV=alog10(abV)+12.
        abMn=alog10(abMn)+12.
        abCo=alog10(abCo)+12.
        write(20,203) absc,abti,abv,abmn,abco
      else
        j=1
        write(20,1202) (inputaname(i),i=(j-1)*20+1,min(92,(j-1)*20+20))
        write(20,1203) (inputabund(i),i=(j-1)*20+1,min(92,(j-1)*20+20))
        do j=2,5
          write(20,1204) (inputaname(i),i=(j-1)*20+1,
     &                            min(92,(j-1)*20+20))
          write(20,1203) (inputabund(i),i=(j-1)*20+1,
     &                            min(92,(j-1)*20+20))
        enddo
      endif
      write(20,255) ITMAX
      READ(IARCH)JTAU,NCORE,DIFLOG,TAUM,RADIUS,(RR(K),K=1,JTAU)
      READ(IARCH)JTAU,(TKORRM(I),I=1,JTAU),(FCORR(K),K=1,JTAU)
      NTPO=0
      DO 3 K=1,JTAU
      READ(IARCH) KR,TAU(K),TAUS(K),Z(K),T(K),PE(K),PG(K),PRAD(K),
     &            PTURB(K),XKAPR(K),RO(K),EMU(K),CP(K),CV(K),
     &            AGRAD(K),Q(K),U(K),V(K),ANCONV(K),HNIC(K),NMOL,
     &            (PRESMO(J,K),J=1,NMOL),PTIO(K)
      TAUK=ALOG10(TAU(K))+10.01
      KTAU=TAUK
      IF(ABS(TAUK-KTAU).GT.0.02) GO TO 31
      IF(KTAU.EQ.10) K0=K
      NTPO=NTPO+1
   31 CONTINUE
    3 CONTINUE
      write(20,204)
      write(20,300) TEFF,GLOG,modnam(1:(index(modnam,' ')))
      write(20,205)
      DO 4 I=1,JTAU
      FCONV(I)=ANCONV(I)*FLUX
      write(20,206) I,TAU(I),T(I),TKORRM(I),FCONV(I),FCORR(I),I
    4 CONTINUE
***
      READ(IARCH)(NJ(I),I=1,NEL),NLP,(XLR(I),I=1,NLP)
     & ,NPROV,NPROVA,NPROVS,(ABNAME(KP),SOURCE(KP),KP=1,NPROV)
***
      write(20,207)
      write(20,300) TEFF,GLOG,modnam(1:(index(modnam,' ')))
      write(20,208) XLR(NLP)
      Z0=Z(K0)
      DO 5 I=1,JTAU
      Z(I)=Z(I)-Z0
* added for osplot*************
        i1=min(i+1,jtau)
        PTOT(I)=PG(I)+PRAD(I)+.5*(pturb(i)+pturb(i1))
      write(20,209) I,TAU(I),TAUS(I),Z(I),T(I),PE(I),PG(I),PRAD(I),
     &             PTURB(I),XKAPR(I),I
    5 CONTINUE
      write(20,210)
      write(20,300) TEFF,GLOG,modnam(1:(index(modnam,' ')))
      write(20,211)
      DO 6 I=1,JTAU
      write(20,212) I,TAU(I),RO(I),EMU(I),CP(I),CV(I),AGRAD(I),Q(I),
     &     U(I),V(I),ANCONV(I),I
    6 CONTINUE
      write(20,213)
      write(20,300) TEFF,GLOG,modnam(1:(index(modnam,' ')))
      write(20,214)
      DO 7 I=1,JTAU
      HNIC(I)=ALOG10(HNIC(I))
      DO 77 J=1,NMOL
      PRESMO(J,I)=ALOG10(PRESMO(J,I))
   77 CONTINUE
       PTIO(I)=AMAX1(PTIO(I),1.E-30)
       PTIOLG=ALOG10(PTIO(I))
      write(20,215) I,HNIC(I),(PRESMO(JJ,I),JJ=1,13),PTIOLG,I
    7 CONTINUE
      write(20,213)
      write(20,300) TEFF,GLOG,modnam(1:(index(modnam,' ')))
      write(20,216)
      DO I=1,JTAU
        write(20,217) I,(PRESMO(J,I),J=14,16),(PRESMO(J,I),J=18,30),I
      enddo
   10 CONTINUE
cc      write(20,260)
ccc      READ(IARCH)(NJ(I),I=1,NEL),NLP,(XLR(I),I=1,NLP)
ccc     & ,NPROV,NPROVA,NPROVS,(ABNAME(KP),SOURCE(KP),KP=1,NPROV)
      HCK=143922240.



** output for multi
      if (xlr(nlp).ne.5000.) then
        print*, 'tau-scale is for lambda = ',xlr(nlp)
        stop 'tau-scale should be for 500nm'
      endif 
      write(30,*) glog
      write(30,*) jtau
      write(40,*) jtau,0.0
      do k=1,jtau
        write(40,30)  log10(taus(k))
        write(30,30)  log10(taus(k)),t(k),pe(k)/t(k)/1.38064852e-16,
     &                0.0,xit
30      format(f11.7,1x,f9.2,1x,1pe12.5,1x0p,f9.5,1x,f9.5)
      enddo
******



      DO 22 KTAU=1,NTPO
      DO 20 IE=1,NEL
      NJP=NJ(IE)
      READ(IARCH) KR,TAUI,TI,PEI,IEL(IE),ABUND(IE),
     &            (ANJON(IE,JJ),JJ=1,NJP),(PART(IE,JJ),JJ=1,NJP)
   20 CONTINUE
*
************************************************************************
* try to improve computation of electron contribution:
* BPz 31/03-2004
*
* first, check that Pg=P(HI)+P(HII)+P(H2)+P(HeI)+Pe (approx, this is neglecting
* any metals, and supposing He is neutral)
* P(He)= fictpressure(H)*abund(He)
*
         Pcheck=10**hnic(kr)*(1.+anjon(1,2)/anjon(1,1))*(1.+abund(2))+
     &          Pei+
     &          10**presmo(2,kr)*(1+2.*abund(2))
         fictPH=10**hnic(kr)*(1.+anjon(1,2)/anjon(1,1))+
     &          10**presmo(2,kr)*2.
         print*,'(Pgas-pcheck)/pgas = ',(pg(kr)-pcheck)/pg(kr),kr
************************************************************************
*
      DO 21 KLAM=1,NLP
      READ(IARCH) KR,TAUIL,(PROV(J,KLAM),J=1,NPROV),
     &            ABSKA(KLAM),SPRIDA(KLAM)
   21 CONTINUE
      IF(KTAU.LE.1) write(20,218)
      IF(KTAU.GT.1) write(20,219)
      write(20,300) TEFF,GLOG,modnam(1:(index(modnam,' ')))
      write(20,220) TAUI
      write(20,221) TI,PEI
      write(20,222)
      PESUM=0.
* treat hydrogen separately, to account for H2 formation
      pep(1)=10**hnic(kr)*anjon(1,2)/anjon(1,1)
      pesum=pesum+pep(1)
      DO IE=2,NEL
        PEP(IE)=0.
        NJP=NJ(IE)
        if (NJP.ge.2) then
          do JJ=2,NJP
            PEP(IE)=PEP(IE)+fictPH*ABUND(IE)*ANJON(IE,JJ)*(JJ-1)
          enddo
          PESUM=PESUM+PEP(IE)
        endif
      enddo
*
*BPz
      print*,' check  (Pe-Pesum)/pe = ', (pe(kr)-pesum)/pei
*
      DO 19 IE=1,NEL
      NJP=NJ(IE)
CCCC      PEP(IE)=PEP(IE)/PESUM
************************************************************************
* BPz:  renormalize to actual Pe. Missing electron donors will then 
*       be accounted for.
*
      PEP(ie)=pep(ie)/pei
*
      if(njp.eq.2) then
        write(20,2232)IEL(IE),ABUND(IE),PEP(IE),(ANJON(IE,JJ),JJ=1,NJP),
     .                (part(ie,jj),jj=1,njp)
      elseif(njp.eq.3) then
        write(20,2233)IEL(IE),ABUND(IE),PEP(IE),(ANJON(IE,JJ),JJ=1,NJP),
     .                (part(ie,jj),jj=1,njp)
      elseif(njp.eq.4) then
        write(20,2234)IEL(IE),ABUND(IE),PEP(IE),(ANJON(IE,JJ),JJ=1,NJP),
     .                (part(ie,jj),jj=1,njp)
      endif
ccccc      write(20,224) (PART(IE,JJ),JJ=1,NJP)
   19 CONTINUE
      write(20,2250)
*
      write(20,225) (ABNAME(KP),KP=1,20)
      HCKT=HCK/TI
      DO 18 KLAM=1,NLP
      ABKLA=ABSKA(KLAM)
      STIM=1.-EXP(-HCKT/XLR(KLAM))
      DO J=1,NPROVA
        PROV(J,KLAM)=PROV(J,KLAM)/ABKLA*STIM
      enddo
      DO J=1,NPROVS
        PROV(NPROVA+J,KLAM)=PROV(NPROVA+J,KLAM)/SPRIDA(KLAM)
      enddo
      write(20,226) XLR(KLAM),ABSKA(KLAM),SPRIDA(KLAM),
     &             (PROV(J,KLAM),J=1,20)
   18 CONTINUE
c  
      write(20,2251) (ABNAME(KP),KP=21,NPROV)
      DO 188 KLAM=1,NLP
        write(20,2261) XLR(KLAM),(PROV(J,KLAM),J=21,NPROV)
  188 CONTINUE 

   22 CONTINUE
*****************
c  there are 
c 1) old models with nothing after flux record (very old)
c 2) models with longname after flux record (older BPz)
c 3) models with weight record after flux record, followed by long name
c    (Bpz on bohor)
c 4) the v7.3 models with flux record, fluxmec record, abundances record, 
c    longname record (BPz v7.3 new models)
c 5) KE models with flux record, fluxmec record, abundance record
c
      if (modtype.eq.3) then
        READ(IARCH) NLB,(XLB(J),FLUXME(J),J=1,NLB)
        read(iarch) (W(J),J=1,NLB)
      else
        READ(IARCH) NLB,(XLB(J),FLUXME(J),J=1,NLB),(W(J),J=1,NLB)
      endif
      if (modtype.eq.4.or.modtype.eq.5) then
        read(iarch,end=79,err=79) (fluxmec(j),j=1,nlb)
        print*,(fluxmec(j),j=10000,10010)
      endif
      goto 78
79    print*,'********************************************************'
      print*,'Error while reading fluxmec in model; skipping that part'
      print*,'********************************************************'
78    continue
C        CONVERT TO 'PHYSICAL' FLUXES
      DO J=1,NLB
        FLUXME(J)=3.14159*FLUXME(J)
        fluxmec(j)=3.14159*fluxmec(j)
      enddo
      DO J=1,NLB
        JNORM=J
        IF(XLB(J).GT.5000.) GOTO 26
      enddo
   26 STMAGN=-2.5*ALOG10(FLUXME(JNORM))
      DO I=1,NLB
        FLUMAG(I)=-2.5*ALOG10(FLUXME(I))-STMAGN
      enddo
C
C        write(20,FLUXES)
C
**** here we jump over fluxes print out
      goto 1111
*
      write(20,521)
      N=0
      DO 52 I=1,NLB,4
      WAVFLX(N+1)=(XLB(I)+XLB(I+1)+XLB(I+2)+XLB(I+3))/4.
      WAVFLX(N+2)=100.*(A*(FLUXME(I)+FLUXME(I+3))+B*(FLUXME(I+1)+
     &      FLUXME(I+2)))
      N=N+2
      IF(N.LT.10) GO TO 52
      write(20,520) WAVFLX
      N=0
   52 CONTINUE
      IF(N.NE.0) write(20,520)(WAVFLX(I),I=1,N)
  520 FORMAT(1X,5(F12.0,E13.5))
  521 FORMAT('1  BELL''S FLUX              (CENTER WAVE & FLUX)'//)
C
C
      write(20,219)
      write(20,300) TEFF,GLOG,modnam(1:(index(modnam,' ')))
      write(20,227)
      LLB=NLB
   60 IF(LLB.LT.NLB) write(20,219)
      IF(LLB.LT.NLB) write(20,300) TEFF,GLOG,
     &       modnam(1:(index(modnam,' ')))
      IF(LLB.LT.NLB) write(20,227)
      write(20,228)
      IRAD = 0
      IF(LLB.LT.NLB) GOTO 65
      I=0
      I2=50
      I3=100
      I4=150
      GOTO 650
   65 CONTINUE
      I = I + 149
      I2 = I2 + 149
      I3 = I3 + 149
      I4 = I4 + 149
  650 CONTINUE
      ILI=2
      IF(LLB.GE.51) ILI=3
      IF(LLB.GE.101)ILI=4
      IF(LLB.GE.151)ILI=5
   66 I4=I4+1
   67 I3=I3+1
   68 I2=I2+1
   69 I=I+1
      IRAD = IRAD + 1
      IF(IRAD.GT.50) GOTO 70
      GOTO(70,61,62,63,64),ILI
   64 write(20,229) XLB(I),FLUXME(I),FLUMAG(I),XLB(I2),FLUXME(I2),
     &             FLUMAG(I2),XLB(I3),FLUXME(I3),FLUMAG(I3),
     &             XLB(I4),FLUXME(I4),FLUMAG(I4)
      IF(I4.GE.NLB) ILI=4
      GOTO 66
   63 write(20,229) XLB(I),FLUXME(I),FLUMAG(I),XLB(I2),FLUXME(I2),
     &             FLUMAG(I2),XLB(I3),FLUXME(I3),FLUMAG(I3)
      IF(I3.GE.NLB) ILI=3
      GOTO 67
   62 write(20,229) XLB(I),FLUXME(I),FLUMAG(I),XLB(I2),FLUXME(I2),
     &             FLUMAG(I2)
      IF(I2.GE.NLB) ILI=2
      GOTO 68
   61 write(20,229) XLB(I),FLUXME(I),FLUMAG(I)
      IF(I.GE.NLB) GOTO 70
      GOTO 69
   70 LLB = LLB - 200
      IF(LLB.GT.0) GOTO 60
*
C
C        COMPUTE U, B, V, R, I AND COLOURS
      SUMCOL=0.
      SUMNOR=0.
      write(20,270)
      DO 40 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.3000..OR.XXLB.GE.4200.) GOTO 40
      IINDEX=INT(XXLB/100.)-29
      VIKT=W(I)*UW(IINDEX)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   40 CONTINUE
      UFLUX=SUMCOL/SUMNOR
      SUMCOL=0.
      SUMNOR=0.
      DO 41 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.3500..OR.XXLB.GE.5600.) GOTO 41
      IINDEX=INT(XXLB/100.)-34
      VIKT=W(I)*BW(IINDEX)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   41 CONTINUE
      BFLUX=SUMCOL/SUMNOR
      SUMCOL=0.
      SUMNOR=0.
      DO 42 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.4700..OR.XXLB.GE.7200.) GOTO 42
      IINDEX=INT(XXLB/100.)-46
      VIKT=W(I)*VW(IINDEX)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   42 CONTINUE
      VFLUX=SUMCOL/SUMNOR
* new colours BPz 290890
      SUMCOL=0.
      SUMNOR=0.
      DO 43 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.5300..OR.XXLB.GE.9500.) GOTO 43
      IINDEX=INT(XXLB/100.)-52
      VIKT=W(I)*RW(IINDEX)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   43 CONTINUE
      RFLUX=SUMCOL/SUMNOR
      SUMCOL=0.
      SUMNOR=0.
      DO 44 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.6900..OR.XXLB.GE.12000.) GOTO 44
      IINDEX=INT(XXLB/100.)-68
      VIKT=W(I)*XIW(IINDEX)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   44 CONTINUE
      XIFLUX=SUMCOL/SUMNOR
      SUMCOL=0.
      SUMNOR=0.
      DO 45 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.9700..OR.XXLB.GE.15600.) GOTO 45
      IINDEX=INT(XXLB/100.)-96
      VIKT=W(I)*XJW(IINDEX)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   45 CONTINUE
      XJFLUX=SUMCOL/SUMNOR
      SUMCOL=0.
      SUMNOR=0.
      DO 46 I=1,NLB
      XXLB=XLB(I)
      IF(XXLB.LT.18100..OR.XXLB.GT.25900.) GOTO 46
      IINDEX=INT((XXLB/200.)-89.5)
      VIKT=W(I)*XKW(IINDEX)
      SUMNOR=SUMNOR+VIKT
      SUMCOL=SUMCOL+VIKT*FLUXME(I)
   46 CONTINUE
      XKFLUX=SUMCOL/SUMNOR
      UMAG=-2.5*ALOG10(UFLUX)-STMAGN
      BMAG=-2.5*ALOG10(BFLUX)-STMAGN
      VMAG=-2.5*ALOG10(VFLUX)-STMAGN
      RMAG=-2.5*ALOG10(RFLUX)-STMAGN
      XIMAG=-2.5*ALOG10(XIFLUX)-STMAGN
      XJMAG=-2.5*ALOG10(XJFLUX)-STMAGN
      XKMAG=-2.5*ALOG10(XKFLUX)-STMAGN
      UB=UMAG-BMAG
      BV=BMAG-VMAG
      UV=UMAG-VMAG
      RI=RMAG-XIMAG
      VR=VMAG-RMAG
      VI=VMAG-XIMAG
      VJ=VMAG-XJMAG
      VK=VMAG-XKMAG
      write(20,300) TEFF,GLOG,modnam(1:(index(modnam,' ')))
      write(20,2698)
      write(20,271) UB,BV,UV
      write(20,272) XLB(JNORM),RI,VR,VI,VJ,VK
      write(20,273) UMAG,BMAG,VMAG,RMAG,XIMAG,xjmag
      write(20,274) XKMAG
***** here we come if we don't print the fluxes nor the colours
1111  continue
    1 CONTINUE
****
  200 FORMAT('M o d e l   P a r a m e t e r s   e t c.'/)
  281 FORMAT(' * The following model was computed by  MARCS35 PP  '
     & ,F5.1,2X,A10,2X,A8,'  *')
  282 FORMAT(' * The following model was computed by  MARCS35 SPH '
     & ,F5.1,2X,A10,2X,A8,'  *')
  283 FORMAT(1X,81('*'))
  284 FORMAT(1X,81('*')////)
  201 FORMAT('Effective temperature',6x,F12.0,'  K'/' Total ',
     +'flux (= sigma*Teff**4)',1PE11.3,'  erg/s/cm**2'/' Acceleration',
     +' of gravity',0PF16.3,'  cm/s**2')
 2001 format(' Plane-parallel model')
 2002 format(' Spherical model with  Radius',1p,e11.3,'  Solar radii',
     ,  8x,'or ',3x,e11.3,' cm',/23x,
     , 'Mass  ',0p,f11.2,'  Solar masses'/23x,'Atm/Radius',2p,f7.2,'  %',
     , /23x,'Lumin.',1p,e11.2,'  Solar luminosities',' or M(bol)',
     , 0p,f8.3)
 2003 format(/' Convection Parameters'/' PALF',
     *'A (L/Hp)=',F5.2,',  PNY (ny)=',F5.2,',  PY (y)=',F5.3/)
  256 FORMAT(' Turbulence pressure is neglected'/)
  257 FORMAT(' Turbulence pressure from convection is included and ',
     ,'PBETA=',F4.2/) 
  258 FORMAT(' Macroturbulence pressure ',
     & 'is included and BETA=',F4.2,' Vmacro= ',f5.2,
     &   ' km/s   constant with depth'/)
  259 FORMAT(' Macroturbulence pressure ',
     & 'is included and BETA=',F4.2,' Vmacro(1)= ',f5.2,
     &   ' km/s   variable with depth'/)
  250 FORMAT(' Convection has been included in this model')
  251 FORMAT(' Convection has  NOT  been included in this model')
  253 FORMAT(' Line blanketing has been included in this model')
  252 FORMAT(' Line blanketing has  NOT  been included in this model')
  254 FORMAT(' The Stroemgren equation has been used for the uppermost',
     *I4,' points'/' Identification    ',A)
  255 FORMAT(//' The following model was obtained after',I3,' itera',
     *'tion(s)')
c...v....1....v....2....v....3....v....4....v....5....v....6....v....7....v....8
  202 FORMAT(//'  Log. Abundances used in model calculations'/2X,'H'
     *,5X,'He',4X,'C',5X,'N',5X,'O',5X,'Ne',4X,'Na',4X,'Mg',4X,'Al',4X,
     *'Si',4X,'S',5X,'K',5X,'Ca',4X,'Cr',4X,'Fe',4X,'Ni')
 2021 format(/2x,'Sc',4x,'Ti',4x,'V',5x,'Mn',4x,'Co')
  203 FORMAT (16F6.2)
 1202 format(//'  Log. Abundances used in model calculations'/
     &         20(2x,a2,3x))
 1204 format(20(2x,a2,3x))
 1203 FORMAT (20F7.3)
c...v....1....v....2....v....3....v....4....v....5....v....6....v....7....v....8
  204 FORMAT(////,'C o r r e c t i o n s   i n   t h e   l a s t   i t e
     * r a t i o n')
  205 FORMAT(/3X,'K',4X,'TauRoss',6X,'T',7X,'Delta(T)',5X,'Fconv',
     *4X,'Delta(Fconv)',2X,'K')
  206 FORMAT(I4,F13.7,F9.2,F11.2,1PE12.3,E12.2,I5)
  207 FORMAT(////,'M o d e l   A t m o s p h e r e     (cgs units)')
  208 FORMAT(/4X,'K',5X,'TauRoss',5X,'Tau(',f6.0,')',1X,
     & 'Geom. Depth',
     & 5X,
     &'T',8X,'Pe',10X,'Pg',9X,'Prad',8X,'Pturb',4X,'KappaRoss',5X,'K')
  209 FORMAT(I5,F12.7,F13.7,1PE13.3,0PF9.0,5(1PE12.3),I6)
c...v....1....v....2....v....3....v....4....v....5....v....6....v....7....v....8
  210 FORMAT(////'T h e r m o d y n a m i c a l   q u a n t i t i e s  
     *and   C o n v e c t i o n   (cgs units)')
  211 FORMAT(/3X,'K',4X,'TauRoss',4X,'Density',5X,'Mu',9X,'Cp',10X,
     *'Cv',7X,'AdGrad',9X,'Q',7X,'Sound vel.',2X,'Conv. vel.',1X,'Fconv/
     *F',3X,'K')
  212 FORMAT(I4,F12.5,1PE11.3,0PF8.3,6(1PE12.3),0PF9.5,I4)
  213 FORMAT(////'L o g a r i t h m i c   M o l e c u l a r   p a r t i 
     *a l   p r e s s u r e s    (cgs units)')
  214 FORMAT(/'  K',3X,'P(H)',2X,'P(H-)',2X,'P(H2)',2X,'P(H2+)',1X,
     *'P(H2O)',1X,'P(OH)',2X,'P(CH)',2X,'P(CO)',2X,'P(CN)',2X,'P(C2)',2X
     *,'P(N2)',2X,'P(O2)',2X,'P(NO)',2X,'P(NH)',1X,'P(TiO)',4X,'K')
  215 FORMAT(I3,15F7.2,2X,I3)
  216 FORMAT(//'  K',3X,'C2H2',3X,'HCN',4X,'C2H',4X,'HS',5X,'SiH',4X
     * ,'C3H',4X,'C3',5X,'CS',5X,'SiC',4X,'SiC2',3X,'NS',5X,'SiN',4X
     * ,'SiO',4X,'SO',5X,'S2',5X,'SiS',5X,'K')
  217 FORMAT(I3,16F7.2,2X,I3)
  218 FORMAT(////'I o n i z a t i o n   c o n d i t i o n s   and   A b 
     *s o r p t i o n  m e c h a n i s m s')
  219 FORMAT(////)
  220 FORMAT(/' *******************'/' * Tau=',F11.7,' *'/' ************
     ********')
  221 FORMAT(/' T=',F7.0,'  Pe=',1PE9.2)
  222 FORMAT(/'  Element  Abundance'
     *                      ,' ElCntb Ionization fractions',17X,'Partiti
     *on functions'/30X,'I',7X,'II',6X,'III',5X,'IV',12X,'I',9X,'II',8X,
     *'III',7X,'IV')
 2232 FORMAT(6X,A2,1PE12.3,0PF6.3,2(F8.4),16x,4x,1p,2e10.2)
 2233 FORMAT(6X,A2,1PE12.3,0PF6.3,3(F8.4),8x,4x,1p,3e10.2)
 2234 FORMAT(6X,A2,1PE12.3,0PF6.3,4(F8.4),4x,1p,4e10.2)
  223 FORMAT(6X,A2,1PE12.3,0PF6.3,4(F8.4))  
  224 FORMAT(62X,4(1PE10.2))
c...v....1....v....2....v....3....v....4....v....5....v....6....v....7....v....8
 2250 FORMAT(//16X,'P e r c e n t a g e s   o f   C o n t i n u o u s'
     * ,'   A b s o r p t i o n   a n d   S c a t t e r i n g')
  225 FORMAT(/'  Wavel     Abs     Scat   ',20(A5))
 2251 format(/'  Wavel     ',24(a5))
  226 FORMAT(F7.0,2(1pE9.2),2p,20F5.1)
 2261 format(f7.0,3x,2p,24f5.1)
  227 FORMAT(' F l u x e s    (Physical fluxes in ergs/s/cm**2/Angstrom)
     *'/)
  228 FORMAT(4(6X,'LAMBDA',5X,'FLUX',7X,'MAGN'))
  229 FORMAT(4(F12.0,1PE11.2,0PF9.3))
  260 FORMAT(///' Molecules other than H2 and H2+ have not been consider
     *ed in the  Temperature - Electron pressure - Pressure  balance')
 2698 FORMAT(//' (UBV Transmission functions after Matthews + Sandage.
     *  Air mass = 1)')
  270 FORMAT(////'G I A N T  L I N E  C O L O U R S')
  271 FORMAT(/8X,'U - B =',F8.3//9X,'B - V =',F8.3//9X,'U - V =',F8.3
     *////)
  272 FORMAT(' TENTATIVE CONTINUUM COLOURS  (R AT',F6.0,'A AND I AT 9000
     *.A)'//9X,'R - I =',F8.3,//9X,'V - R =',F8.3//9X,'V - I =',F8.3
     &//9X,'V - J =',F8.3//9X,'V - K =',F8.3)
  273 FORMAT(////' U =',F10.3,'  B =',F10.3,'  V =',F10.3,'  R =',F10.3,
     *'  I =',F10.3,'  J =',F10.3)
  274 FORMAT(//' K =',F10.3,'  L =',F10.3,'  M =',F10.3,'  N =',F10.3)
  300 FORMAT('Teff=',F6.0,' log g=',F5.2,1X,a)
*
        print*,nlb
        print*,(xlb(j),j=10000,10010)
        print*,(w(j),j=10000,10010)
        print*,(fluxme(j),j=10000,10010)
        print*,(fluxmec(j),j=10000,10010)
        print*,fluxme(107855), lpoint,fluxme(lpoint)
        print*,'*** returning from readmo'
      RETURN
         END
