      program vald3line
* Convert VALD3 "Extract all, element or stellar" long format output to BSYN line data lines
* Bengt Edvardsson April 20, 1999 (for VALD2)
* Adapted to new VALD3 output format October 2011 /BE
* Isotopic identifications and neutral diatomic molecules included August 2013 /BE
*
*   Free format version for first line of each species in VALD3 data
*    BPz 5/02-2019
*
      implicit none
      integer       imax,maxatom
      parameter     (imax=300000)
      parameter     (maxatom=92)
      integer       i,ii,n,nread,nlines,length,ionic,nswap,nskip
      integer       ion(imax),ielem,iel(imax),nspec,nu(imax),irec
      real*8        spec(imax),sss
      real          sunabund(maxatom),ionp1(maxatom),ionp2(maxatom)
      real          w,gflog,chil,chiu,rjupper,gamrad,fdamp,eqw,eqwerr
      real          gam6,dum,gamstark
      character*300 string,filename
      character*500 string4,s4(imax)
      character*200 string1,s1(imax)
      character*92  string2,s2(imax),string3,s3(imax)
      character*91  lowdesig,highdesig
      character*14  full,fullspec(imax)
      character*11  species,nucleus(imax)
      character*8   blipblop,cdum
      character*7   element
      character*2   lele(maxatom),el1(imax),el2(imax),elem1,elem2
      character*2   cion(2),lowcoupling,highcoupling
      character*1   lower,upper
      logical       molecule
*
      data cion/'I ','II'/
      data eqw/0.0/
      data eqwerr/1.0/
      data lele /
     &           'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     &           'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     &           'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     &           'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     &           'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     &           'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     &           'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     &           'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     &           'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     &           'Pa','U '/
*
* Solar abundances ref: Grevesse N., Asplund A., Sauval A.J. 2007,
* Space Science Review 130, 105, DOI 10.1007/s11214-007-9173-7   
      data sunabund /
     & 12.00,10.9275, 1.05,  1.38,  2.70,  8.39,  7.78,  8.66,  4.56,   |  1 -  9
     &  7.84,  6.17,  7.53,  6.37,  7.51,  5.36,  7.14,  5.50,  6.18,   | 10 - 18
     &  5.08,  6.31,  3.17,  4.90,  4.00,  5.64,  5.39,  7.45,  4.92,   | 19 - 27
     &  6.23,  4.21,  4.60,  2.88,  3.58,  2.29,  3.33,  2.56,  3.25,   | 28 - 36
     &  2.60,  2.92,  2.21,  2.58,  1.42,  1.92, -99.0,  1.84,  1.12,   | 37 - 45
     &  1.66,  0.94,  1.77,  1.60,  2.00,  1.00,  2.19,  1.51,  2.24,   | 46 - 54
     &  1.07,  2.17,  1.13,  1.70,  0.58,  1.45, -99.0,  1.00,  0.52,   | 55 - 63
     &  1.11,  0.28,  1.14,  0.51,  0.93,  0.00,  1.08,  0.06,  0.88,   | 64 - 72
     & -0.17,  1.11,  0.23,  1.25,  1.38,  1.64,  1.01,  1.13,  0.90,   | 73 - 81
     &  2.00,  0.65, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0,  0.06,   | 82 - 90
     & -99.0, -0.52 /                                                   | 91 - 92
*
*
* Solar abundances ref: Grevesse Sauval 1998, Space Science Rev 85, 161-174
*     data sunabund /
*    & 12.00, 10.93,  1.10,  1.40,  2.55,  8.52,  7.92,  8.83,  4.56,   |  1 -  9
*    &  8.08,  6.33,  7.58,  6.47,  7.55,  5.45,  7.33,  5.50,  6.40,   | 10 - 18
*    &  5.12,  6.36,  3.17,  5.02,  4.00,  5.67,  5.39,  7.50,  4.92,   | 19 - 27
*    &  6.25,  4.21,  4.60,  2.88,  3.41,  2.37,  3.41,  2.63,  3.31,   | 28 - 36
*    &  2.60,  2.97,  2.24,  2.60,  1.42,  1.92, -99.0,  1.84,  1.12,   | 37 - 45
*    &  1.69,  0.94,  1.77,  1.66,  2.00,  1.00,  2.24,  1.51,  2.17,   | 46 - 54
*    &  1.13,  2.13,  1.17,  1.58,  0.71,  1.50, -99.0,  1.01,  0.51,   | 55 - 63
*    &  1.12, -0.10,  1.14,  0.26,  0.93,  0.00,  1.08,  0.06,  0.88,   | 64 - 72
*    & -0.13,  1.11,  0.28,  1.45,  1.35,  1.80,  1.01,  1.13,  0.90,   | 73 - 81
*    &  1.95,  0.71, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0,  0.09,   | 82 - 90
*    & -99.0, -0.47 /                                                   | 91 - 92
*
* Solar abundances ref: Grevesse Noels Sauval 1996 ASP Conf 99, 117
*     data sunabund /
*    & 12.00, 10.99,  1.16,  1.15,  2.60,  8.55,  7.97,  8.87,  4.56,
*    &  8.08,  6.33,  7.58,  6.47,  7.55,  5.45,  7.33,  5.50,  6.52,
*    &  5.12,  6.36,  3.17,  5.02,  4.00,  5.67,  5.39,  7.50,  4.92,
*    &  6.25,  4.21,  4.60,  2.88,  3.41,  2.37,  3.38,  2.63,  3.23,
*    &  2.60,  2.97,  2.24,  2.60,  1.42,  1.92, -9.00,  1.84,  1.12,
*    &  1.69,  0.94,  1.77,  1.66,  2.00,  1.00,  2.24,  1.51,  2.23,
*    &  1.13,  2.13,  1.17,  1.58,  0.71,  1.50, -9.00,  1.01,  0.51,
*    &  1.12, -0.10,  1.14,  0.26,  0.93,  0.00,  1.08,  0.76,  0.88,
*    & -0.13,  1.11,  0.28,  1.45,  1.35,  1.80,  1.01,  1.17,  0.90,
*    &  1.95,  0.71, -9.00, -9.00, -9.00, -9.00, -9.00, -9.00,  0.09,
*    & -9.00, -0.47 /
*
* Ionization potentials (eV), from file "atomdata"
      data ionp1 /
     & 13.60, 24.59,  5.39,  9.32,  8.30, 11.26, 14.53, 13.62, 17.40,
     & 21.56,  5.14,  7.64,  5.99,  8.15, 10.48, 10.36, 12.97, 15.76,
     &  4.34,  6.11,  6.54,  6.82,  6.74,  6.77,  7.44,  7.87,  7.86,
     &  7.64,  7.73,  9.39,  6.00,  7.88,  9.81,  9.75, 11.81, 14.00,
     &  4.18,  5.69,  6.38,  6.84,  6.88,  7.10,  7.28,  7.36,  7.46,
     &  8.33,  7.57,  8.99,  5.79,  7.34,  8.64,  9.01, 10.45, 12.13,
     &  3.89,  5.21,  5.58,  5.65,  5.42,  5.49,  5.55,  5.63,  5.68,
     &  6.16,  5.85,  5.93,  6.02,  6.10,  6.18,  6.25,  5.43,  7.00,
     &  7.88,  7.98,  7.87,  8.50,  9.10,  9.00,  9.22, 10.44,  6.11,
     &  7.42,  7.29,  8.42,  9.30, 10.75,  4.00,  5.28,  6.90,  6.08,
     &  9.99,  6.00/
*
      data ionp2 /
     &   .00, 54.42, 75.64, 18.21, 25.15, 24.38, 29.60, 35.12, 35.00,
     & 40.96, 47.29, 15.03, 18.83, 16.35, 19.72, 23.33, 23.81, 27.62,
     & 31.63, 11.87, 12.80, 13.58, 14.65, 16.50, 15.64, 16.18, 17.06,
     & 18.17, 20.29, 17.96, 20.51, 15.93, 18.63, 21.19, 21.60, 24.36,
     & 27.50, 11.03, 12.24, 13.13, 14.32, 16.15, 15.26, 16.76, 18.07,
     & 19.42, 21.49, 16.90, 18.87, 14.63, 16.53, 18.60, 19.13, 21.21,
     & 25.10, 10.00, 11.06, 10.85, 10.55, 10.73, 10.90, 11.07, 11.25,
     & 12.10, 11.52, 11.67, 11.80, 11.93, 12.05, 12.17, 13.90, 14.90,
     & 16.20, 17.70, 16.60, 17.50, 20.00, 18.56, 20.50, 18.76, 20.43,
     & 15.03, 16.68, 19.00, 20.00, 21.00, 22.00, 10.14, 12.10, 11.50,
     & 99.99, 12.00/
*
      print*, ' Preliminary version of conversion program'
      print*,' The format for molecular lines is not compatible'
      print*,' with B. Plez''''s version of turbospectrum.'
      print*,' In order to use Barklem et al broadening parameters,'
      print*,' the VALD request should be with extended VdW syntax.'
      print*,' this optin can be selected through Unit selection'
      print*,' in the VALD3 interface.'
      print*
      print *, 'VALD3 ''show line'' (long format) file name?'
      read(*,'(a)') filename
      open(10,file=filename,form='formatted',status='old')
      print *,'output bsyn line file name?'
      read(*,'(a)') filename
      open(11,file=filename,form='formatted',status='unknown')
*
      irec=0
      nread=0
      nlines=0
      nskip=0
      do
          read(10,'(a)',end=99) string
          irec=irec+1
          read(string,'(a8)') blipblop
          if(blipblop(1:1).eq.''''
     &     .and. (blipblop(5:5).eq.'''' 
     &        .or.blipblop(6:6).eq.''''
     &        .or.blipblop(7:7).eq.''''
     &        .or.blipblop(8:8).eq.'''')) then
* This is probably the first record for a long format line
            nread=nread+1
            length=len_trim(string)
            string1=string(1:length)
* string1 contains the wavelength, gf-value etc
* e.g. 'Fe 1',       4080.0614,  -3.664,  4.1777,  5.0,  7.2156,  4.0,  1.420,  1.570,  1.110, 8.280,-4.890,  -7.180, 0.009,
            read(10,'(a)') string
            irec=irec+1
            length=len_trim(string)
            if(length.eq.0) then
              print *,'No lower level data, file record number',irec
              print *,string1
cc            else if(length.ne.92) then
cc              print *,string,' : length=',length
cc              stop 'Wrong length of the second record'
            endif
            string2=string(1:length)
* string2 contains the lower level information
* e.g. '  LS                                                                      3d7.(4F).4p y5F*'
            read(10,'(a)') string
            irec=irec+1
            length=len_trim(string)
            if(length.eq.0) then
              print *,'No upper level data, file record number',irec
              print *,string1
              print *,string2
cc            else if(length.ne.92) then
cc              print *,string,' : length=',length
cc              stop 'Wrong length of the third record'
            endif
            string3=string(1:length)
* string3 contains the upper level information
* e.g. '  LS                                                                     3d6.4s.(6D).7s 7D'
            read(10,'(a)') string4
            irec=irec+1
* string4 contains identification and references
* for extract all:
* e.g. '_                             5 wl:BBSB                          5 gf:BBSB                          5 BBSB                             5 BBSB                             5 BBSB                             5 BBSB                             5 BBSB                             5 BBSB                             5 BBSB                           (12)C(12)C    '
* for extract stellar: :-(
* e.g. '_          Kurucz Fe I 2014Fe               1 wl:K14   1 gf:K14   1 K14   1 K14   1 K14   1 K14   1 K14   1 K14   1 K14'
* e.g. '_                          (12)CH          13 wl:JLIY  13 gf:JLIY  13 JLIY  13 JLIY  13 JLIY  13 JLIY  13 JLIY  13 JLIY  13 JLIY'
* e.g. '_                          (12)C(12)C      23 wl:BBSB  23 gf:BBSB  23 BBSB  23 BBSB  23 BBSB  23 BBSB  23 BBSB  23 BBSB  23 BBSB'
            if (len_trim(string4).gt.300) then
              full=string4(345:358)
            else
              full=string4(29:42)
            endif
            call identify(full,lele,species,ionic,elem1,elem2)
            read(species,'(f11.6)') sss
            ielem=int(sss)
            if(sss.lt.93.) then
              molecule=.false.
              sss=sss+0.0001*dble(ionic)
            else
              molecule=.true.
            endif
              sss=sss+0.0000001*dble(ionic)
* Deselect not wanted lines
            if(sss.lt.3.) then
              print *,'skipping H and He lines: ',string1
*             print *,'input file record number',irec
              nskip=nskip+1
              goto 90
            endif
            if(ionic.ne.1 .and. ionic.ne.2) then
              print *,'skipping ion:         ',ionic,' ',string1
*             print *,'input file record number',irec
              nskip=nskip+1
              goto 90
            endif
            read(string1,*) cdum,dum,dum,chil,dum,chiu
            if(chil.ge.15.) then
              print *,'skipping chil > 15 eV:            ',string1
*             print *,'input file record number',irec
              nskip=nskip+1
              goto 90
            endif
            if(.not.molecule) then
              if(ionic.eq.1) then
                if(chiu.ge.ionp1(ielem)) then
* Remove lines w upper level in the continuum (incl. autoion lines)
                  print *,'skipping bound-free transition:(?)',string1
*                 print *,'input file record number',irec
                  nskip=nskip+1
                  goto 90
                endif
              else
                if(chiu.ge.ionp2(ielem)) then
* Remove lines w upper level in the continuum (incl. autoion lines)
                  print *,'skipping bound-free transition:(?)',string1
*                 print *,'input file record number',irec
                  nskip=nskip+1
                  goto 90
                endif
              endif
            endif
* The rest of the lines should probably be kept:
            nlines=nlines+1
            fullspec(nlines)=full
            nucleus(nlines)=species
            spec(nlines)=sss
            ion(nlines)=ionic
            el1(nlines)=elem1
            el2(nlines)=elem2
            iel(nlines)=ielem
            s1(nlines)=string1
            s2(nlines)=string2
            s3(nlines)=string3
            s4(nlines)=string4
          endif
   90   continue
      enddo
   99 continue
      close(10)
*
* OK data read, bubble sort the lines
  110 continue
      nswap=0
      do i=2,nlines
        if(spec(i-1).gt.spec(i)) then
          sss=spec(i-1)
          full=fullspec(i-1)
          species=nucleus(i-1)
          elem1=el1(i-1)
          elem2=el2(i-1)
          ionic=ion(i-1)
          ielem=iel(i-1)
          string1=s1(i-1)
          string2=s2(i-1)
          string3=s3(i-1)
          string4=s4(i-1)
          spec(i-1)=spec(i)
          fullspec(i-1)=fullspec(i)
          nucleus(i-1)=nucleus(i)
          el1(i-1)=el1(i)
          el2(i-1)=el2(i)
          ion(i-1)=ion(i)
          iel(i-1)=iel(i)
          s1(i-1)=s1(i)
          s2(i-1)=s2(i)
          s3(i-1)=s3(i)
          s4(i-1)=s4(i)
          spec(i)=sss
          fullspec(i)=full
          nucleus(i)=species
          el1(i)=elem1
          el2(i)=elem2
          ion(i)=ionic
          iel(i)=ielem
          s1(i)=string1
          s2(i)=string2
          s3(i)=string3
          s4(i)=string4
          nswap=nswap+1
        endif
      enddo
      if(nswap.ne.0) goto 110
* How many lines of each species are there 'nu()'?
      nspec=1
      nu(1)=1
      full=fullspec(1)
      do i=2,nlines
        if(fullspec(i).eq.full) then
          nu(nspec)=nu(nspec)+1
        else
          nspec=nspec+1
          nu(nspec)=nu(nspec)+1
          full=fullspec(i)
        endif
      enddo
*
*OK, extract data and write
      n=1
      do ii=1,nspec
        molecule=.true.
        if(spec(n).lt.93.) then
          molecule=.false.
        endif
        if(.not.molecule) then
          write(element,'(a2,x,a2,2x)') el1(n),cion(ion(n))
          if(element(2:3).eq.'  ') then
            do i=2,6
              element(i:i)=element(i+1:i+1)
            enddo
            element(7:7)=' '
          endif
          if(element(3:4).eq.'  ') then
            do i=3,6
              element(i:i)=element(i+1:i+1)
            enddo
            element(7:7)=' '
          endif
          write(11,1120) nucleus(n),ion(n),nu(ii),element
 1120     format('''',2x,a11,7x,'''',i5,i10,/,'''',a7,'''')
        else
          write(element,'(a2,a2,x,a2)') el1(n),el2(n),cion(ion(n))
          if(element(2:2).eq.' ') then
            do i=2,6
              element(i:i)=element(i+1:i+1)
            enddo
            element(7:7)=' '
          endif
          if(element(2:3).eq.'  ') then
            do i=2,6
              element(i:i)=element(i+1:i+1)
            enddo
            element(7:7)=' '
          endif
          if(element(3:4).eq.'  ') then
            do i=3,6
              element(i:i)=element(i+1:i+1)
            enddo
            element(7:7)=' '
          endif
          if(element(4:5).eq.'  ') then
            do i=4,6
              element(i:i)=element(i+1:i+1)
            enddo
            element(7:7)=' '
          endif
          write(11,1121) nucleus(n),ion(n),nu(ii),element
 1121     format('''',a11,9x,'''',i5,i10,/,'''',a7,'''')
        endif
        do i=1,nu(ii)
          gam6=0.0
          read(s1(n),*) cdum,w,gflog,chil,dum,chiu,rjupper,
     &         dum,dum,dum,gamrad,gamstark,gam6                ! lande factors (3 dummies); gamstark not used yet
          if(.not.molecule) then
* Atomic species
* get orbital types for atomic lines
            call define_transition(s2(n),s3(n),lower,upper)
* Apply FDAMP, for references See BDP (A&A 275,101) or notes below
            if(ion(n).eq.1) then
              if(iel(n).eq.11) then
* Holweger 1971 A&A 10, 128
                fdamp=2.0
              else if(iel(n).eq.14) then
* Holweger 1973 A&A 26, 275
                fdamp=1.3
              else if(iel(n).eq.20) then
* O'Neill & Smith 1980 A&A 81, 100
                fdamp=1.8
              else if(iel(n).eq.26) then
* Simmons & Blackwell 1982 A&A 112, 209
* Magain & Zhao 1996 A&A 305, 245
                fdamp=1.4
              else
* Mackle et al. 1975 A&A 38, 239
                fdamp=2.5
              endif
            else 
              if(iel(n).eq.20) then
                fdamp=1.4
* from fit of H&K in the HM model to the fts intensity spectrum
              else if(iel(n).eq.38) then
                fdamp=1.8
* from fit of Sr II 4077.724 in the HM model to the fts intensity spectrum
              else if(iel(n).eq.56) then
* Holweger & Muller 1974 Solar Physics 39, 19
                fdamp=3.0
              else
* Mackle et al. 1975 A&A 38, 239
                fdamp=2.5
              endif
            endif
          else
* for molecules:
* Mackle et al. 1975 A&A 38, 239
            fdamp=2.5
            lower='X'
            upper='X'
          endif
          if (gam6.ne.0.) then
* if there is gamma6 data for this line, we store it in fdamp. BPz 05/06-2014
            fdamp=gam6
          endif
          if(gamrad.gt.3.0) then
            gamrad=10.0**gamrad
          else
* place holder:
            gamrad=1.e5
          endif
          call shrinkdesignation(s2(n),lowcoupling,lowdesig)
          call shrinkdesignation(s3(n),highcoupling,highdesig)
          if(n.gt.imax-1) then
            print *,'WARNING More than',imax,' lines of a species',n
            stop 'increase the imax parameter'
          endif
          if(gflog.gt.-10.0 .and. gflog.lt.100.0) then
            write(11,1130) w,chil,gflog,fdamp,
     &      2.*rjupper+1.,gamrad,lower,upper,eqw,eqwerr,
     &      element(1:len_trim(element)),
     &      lowcoupling,':',lowdesig(1:len_trim(lowdesig)),
     &      highcoupling,':',highdesig(1:len_trim(highdesig))
* BPz format changed to accomodate extended VdW information
 1130       format(f10.3,x,f6.3,x,f6.3,x,f8.3,x,f6.1,x,1p,e9.2,0p,
     &             x,'''',a1,'''',x,'''',a1,'''',x,f5.1,x,f6.1,x,'''',
     &             a,x,3a,x,3a,'''')
          else
            write(11,1131) w,chil,gflog,fdamp,
     &      2.*rjupper+1.,gamrad,lower,upper,eqw,eqwerr,
     &      element(1:len_trim(element)),
     &      lowcoupling,':',lowdesig(1:len_trim(lowdesig)),
     &      highcoupling,':',highdesig(1:len_trim(highdesig))
* BPz format changed to accomodate extended VdW information
 1131       format(f10.3,x,f6.3,x,f6.2,x,f8.3,x,f6.1,x,1p,e9.2,0p,
     &             x,'''',a1,'''',x,'''',a1,'''',x,f5.1,x,f6.1,x,'''',
     &             a,x,3a,x,3a,'''')
          endif
          nlines=n
          n=n+1
        enddo
      enddo
      print *,nskip,' lines skipped'
      if(nread.ne.nlines+nskip) then
        print *,'WARNING: lines missed: read=',nread,' written=',nlines
        print *,' skipped=',nskip,'. Where are the rest?'
      endif
      print *,
     &     nlines,' lines of',nspec,' species written to file ',filename
      close(11)
      end
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
      subroutine identify(fullspec,lele,species,ion,lek,lel)
* identify atomic or molecular lines, ionization, and isotopes
*
      implicit none
      integer      length,itype,iso(2),istart1,istart2,iiso2,i,iel(2)
      integer      ichar(14),ion,natom,iatom(2)
      character*14 fullspec,fff
      character*11 species
      character*6  isostring
      character*2  lele(92),lek,lel
      character*1  char(14)
      equivalence (fff,char)
*
      fff=fullspec
      iel(1)=0
      iel(2)=0
      iso(1)=0
      iso(2)=0
      ion=0
      istart1=1
      istart2=0
      length=len_trim(fullspec)
      if(char(1).eq.'(') then
* Cases with isotopes:
        if(char(3).eq.')') then
          read(char(2),'(i1)') iso(1)
          istart1=4
        else if(char(4).eq.')') then
          read(fullspec(2:3),'(i2)') iso(1)
          istart1=5
        else if(char(5).eq.')') then
          read(fullspec(2:4),'(i3)') iso(1)
          istart1=6
        else
          print *,fullspec
          stop 'Unexpected ID, first isotope???'
        endif
        do i=istart1+1,11
          if(char(i).eq.'(') then
            iiso2=i+1
            goto 10
          endif
        enddo
        goto 30
   10   continue
* OK, a second isotope found, i.e. a molecule
        if(char(iiso2+1).eq.')') then
          read(char(iiso2),'(i1)') iso(2)
          istart2=iiso2+2
        else if(char(iiso2+2).eq.')') then
          read(fullspec(iiso2:iiso2+1),'(i2)') iso(2)
          istart2=iiso2+3
        else if(char(iiso2+3).eq.')') then
          read(fullspec(iiso2:iiso2+2),'(i3)') iso(2)
          istart2=iiso2+4
        else
          print *,fullspec
          stop 'Unexpected ID, second isotope???'
        endif
      else
* Is there only a second element isotope?
        do i=2,11
          if(char(i).eq.'(') then
            iiso2=i+1
            goto 20
          endif
        enddo
        goto 30
   20   continue
* OK, only a second isotope found, i.e. a molecule
        if(char(iiso2+1).eq.')') then
          read(char(iiso2),'(i1)') iso(2)
          istart2=iiso2+2
        else if(char(iiso2+2).eq.')') then
          read(fullspec(iiso2:iiso2+1),'(i2)') iso(2)
          istart2=iiso2+3
        else if(char(iiso2+3).eq.')') then
          read(fullspec(iiso2:iiso2+2),'(i3)') iso(2)
          istart2=iiso2+4
        else
          print *,fullspec
          stop 'Unexpected ID, second isotope???'
        endif
      endif
   30 continue
* OK, we now have any isotope information, find 1 or 2 elements:
      natom=0
      lek='  '
      lel='  '
      do i=1,14
        call case(.false.,char(i),ichar(i))
        if(ichar(i).eq.1) then
          natom=natom+1
          iatom(natom)=i
        endif
      enddo
      lek(1:1)=char(iatom(1))
      if(ichar(iatom(1)+1).ne.-1) then
* 1-character atom
        lek(2:2)=' '
        call findatom(lek,lele,iel(1))
        if(char(iatom(1)+1).eq.'+' .and. ichar(iatom(1)+2).eq.3) then
          read(char(iatom(1)+2),'(i1)') ion
          ion=ion+1
        else if(char(iatom(1)+1).eq.'+') then
          ion=2
        else if(char(iatom(1)+1).eq.' ') then
          ion=1
        else if(ichar(iatom(1)+1).eq.1) then
          continue
        else if(char(iatom(1)+1).eq.'(') then
          continue
        else
          print *,fullspec
          stop '1 What is this?'
        endif
*     print *,'full,lek,ion,iel(1-2)',fullspec,lek,ion,iel
      else
* 2-character atom
        lek(2:2)=char(iatom(1)+1)
        call findatom(lek,lele,iel(1))
        if(char(iatom(1)+2).eq.'+' .and. ichar(iatom(1)+3).eq.3) then
          read(char(iatom(1)+3),'(i1)') ion
          ion=ion+1
        else if(char(iatom(1)+2).eq.'+') then
          ion=2
        else if(char(iatom(1)+2).eq.' ') then
          ion=1
        else if(ichar(iatom(1)+2).eq.1) then
          continue
        else if(char(iatom(1)+2).eq.'(') then
          continue
        else
          print *,fullspec
          stop '2 What is this?'
        endif
*     print *,'FULL,lek,ion,iel(1-2)',fullspec,lek,ion,iel
      endif
      if(natom.eq.2) then
        lel(1:1)=char(iatom(2))
        if(ichar(iatom(2)+1).ne.-1) then
* 1-character atom
          lel(2:2)=' '
          call findatom(lel,lele,iel(2))
          if(char(iatom(2)+1).eq.'+' .and. ichar(iatom(2)+2).eq.3) then
            read(char(iatom(2)+2),'(i1)') ion
            ion=ion+1
          else if(char(iatom(2)+1).eq.'+') then
            ion=2
          else if(char(iatom(2)+1).eq.' ') then
            ion=1
          else
            print *,fullspec
            stop '3 What is this?'
          endif
        else
* 2-character atom
          lel(2:2)=char(iatom(2)+1)
          call findatom(lel,lele,iel(2))
          if(iatom(2).lt.13) then
            if(char(iatom(2)+2).eq.'+' .and. ichar(iatom(2)+3).eq.3)then
              read(char(iatom(2)+3),'(i1)') ion
              ion=ion+1
            else if(char(iatom(2)+2).eq.'+') then
              ion=2
            else if(char(iatom(2)+2).eq.' ') then
              ion=1
            else
              print *,fullspec
              stop '4 What is this?'
            endif
          else
            ion=1
          endif
        endif
      endif
* construct species identity number and string
      if(iel(2).eq.0) then
        write(species,'(i2,a1,i3)') iel(1),'.',iso(1)
        if(species(4:4).eq.' ') species(4:4)='0'
        if(species(5:5).eq.' ') species(5:5)='0'
        species(7:11)='     '
      else
        if(iel(1).le.iel(2)) then
          write(species,'(i2,i2,a1,i3,i3)') 
     &                    iel(1),iel(2),'.',iso(1),iso(2)
          if(species(3:3).eq.' ') species(3:3)='0'
          if(species(6:6).eq.' ') species(6:6)='0'
          if(species(7:7).eq.' ') species(7:7)='0'
          if(species(9:9).eq.' ') species(9:9)='0'
          if(species(10:10).eq.' ') species(10:10)='0'
        else
          write(species,'(i2,i2,a1,i3,i3)') 
     &                    iel(2),iel(1),'.',iso(2),iso(1)
          if(species(3:3).eq.' ') species(3:3)='0'
          if(species(6:6).eq.' ') species(6:6)='0'
          if(species(7:7).eq.' ') species(7:7)='0'
          if(species(9:9).eq.' ') species(9:9)='0'
          if(species(10:10).eq.' ') species(10:10)='0'
        endif
      endif
*     print *,'FULLSPEC,LEK,LEL,IEL,ION,ISO,SPECIES',
*    &         fullspec,lek,lel,iel,ion,iso,'|',species,'|'
      return
      end
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
      subroutine case(change,c,ichar)
* changes upper to lower case if  change. and ichar=2
* changes lower to upper case if  change. and ichar=-2
* or
* determines whether "c" is upper case if .not.change.: output ichar=1
* determines whether "c" is lower case if .not.change.: output ichar=-1
* determines whether "c" is a number if .not.change.: output ichar=3
* or
* determines that "c" is not a letter/number if .not.change.: output ichar=0
* 
* c and ichar can be either input or output or both
      integer ichar
      character*1 c
      logical change
*
      if(change) then
        if(ichar.eq.2) then
          if(c.eq.'A') then
                c='a'
          else if(c.eq.'B') then
                     c='b'
          else if(c.eq.'C') then
                     c='c'
          else if(c.eq.'D') then
                     c='d'
          else if(c.eq.'E') then
                     c='e'
          else if(c.eq.'F') then
                     c='f'
          else if(c.eq.'G') then
                     c='g'
          else if(c.eq.'H') then
                     c='h'
          else if(c.eq.'I') then
                     c='i'
          else if(c.eq.'J') then
                     c='j'
          else if(c.eq.'K') then
                     c='k'
          else if(c.eq.'L') then
                     c='l'
          else if(c.eq.'M') then
                     c='m'
          else if(c.eq.'N') then
                     c='n'
          else if(c.eq.'O') then
                     c='o'
          else if(c.eq.'P') then
                     c='p'
          else if(c.eq.'Q') then
                     c='q'
          else if(c.eq.'R') then
                     c='r'
          else if(c.eq.'S') then
                     c='s'
          else if(c.eq.'T') then
                     c='t'
          else if(c.eq.'U') then
                     c='u'
          else if(c.eq.'V') then
                     c='v'
          else if(c.eq.'X') then
                     c='x'
          else if(c.eq.'Y') then
                     c='y'
          else if(c.eq.'Z') then
                     c='z'
          endif
        else if(ichar.eq.-2) then
          if(c.eq.'a') then
                c='A'
          else if(c.eq.'b') then
                     c='B'
          else if(c.eq.'c') then
                     c='C'
          else if(c.eq.'d') then
                     c='D'
          else if(c.eq.'e') then
                     c='E'
          else if(c.eq.'f') then
                     c='F'
          else if(c.eq.'g') then
                     c='G'
          else if(c.eq.'h') then
                     c='H'
          else if(c.eq.'i') then
                     c='I'
          else if(c.eq.'j') then
                     c='J'
          else if(c.eq.'k') then
                     c='K'
          else if(c.eq.'l') then
                     c='L'
          else if(c.eq.'m') then
                     c='M'
          else if(c.eq.'n') then
                     c='N'
          else if(c.eq.'o') then
                     c='O'
          else if(c.eq.'p') then
                     c='P'
          else if(c.eq.'q') then
                     c='Q'
          else if(c.eq.'r') then
                     c='R'
          else if(c.eq.'s') then
                     c='S'
          else if(c.eq.'t') then
                     c='T'
          else if(c.eq.'u') then
                     c='U'
          else if(c.eq.'v') then
                     c='V'
          else if(c.eq.'x') then
                     c='X'
          else if(c.eq.'y') then
                     c='Y'
          else if(c.eq.'z') then
                     c='Z'
          endif
        else
          stop 'subr. change called with wrong data'
        endif
      else
* determine letter case:
        if(c.eq.'a'.or.
     &     c.eq.'b'.or.
     &     c.eq.'c'.or.
     &     c.eq.'d'.or.
     &     c.eq.'e'.or.
     &     c.eq.'f'.or.
     &     c.eq.'g'.or.
     &     c.eq.'h'.or.
     &     c.eq.'i'.or.
     &     c.eq.'j'.or.
     &     c.eq.'k'.or.
     &     c.eq.'l'.or.
     &     c.eq.'m'.or.
     &     c.eq.'n'.or.
     &     c.eq.'o'.or.
     &     c.eq.'p'.or.
     &     c.eq.'q'.or.
     &     c.eq.'r'.or.
     &     c.eq.'s'.or.
     &     c.eq.'t'.or.
     &     c.eq.'u'.or.
     &     c.eq.'v'.or.
     &     c.eq.'w'.or.
     &     c.eq.'x'.or.
     &     c.eq.'y'.or.
     &     c.eq.'z') then
          ichar=-1
        else if(c.eq.'A'.or.
     &     c.eq.'B'.or.
     &     c.eq.'C'.or.
     &     c.eq.'D'.or.
     &     c.eq.'E'.or.
     &     c.eq.'F'.or.
     &     c.eq.'G'.or.
     &     c.eq.'H'.or.
     &     c.eq.'I'.or.
     &     c.eq.'J'.or.
     &     c.eq.'K'.or.
     &     c.eq.'L'.or.
     &     c.eq.'M'.or.
     &     c.eq.'N'.or.
     &     c.eq.'O'.or.
     &     c.eq.'P'.or.
     &     c.eq.'Q'.or.
     &     c.eq.'R'.or.
     &     c.eq.'S'.or.
     &     c.eq.'T'.or.
     &     c.eq.'U'.or.
     &     c.eq.'V'.or.
     &     c.eq.'W'.or.
     &     c.eq.'X'.or.
     &     c.eq.'Y'.or.
     &     c.eq.'Z') then
          ichar=1
        else if(c.eq.'0'.or.
     &     c.eq.'1'.or.
     &     c.eq.'2'.or.
     &     c.eq.'3'.or.
     &     c.eq.'4'.or.
     &     c.eq.'5'.or.
     &     c.eq.'6'.or.
     &     c.eq.'7'.or.
     &     c.eq.'8'.or.
     &     c.eq.'9') then
          ichar=3
        else
          ichar=0
        endif
      endif
      return
      end
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
      subroutine lineid(des,fullspecies,lele,el,iel,ion,isostring)
* Identify species, ionisation state and isotopes
      implicit none
      integer iel(2),imol,ion,isot,iso(2),i,ioffset,maxatom
      parameter (maxatom=92)
      integer ichar2,ichar3,ichar4,ichar5,ichar6
      character des*6,fullspecies*30,isostring*6,onechar*1
      character*2 lele(maxatom),el(2),isoel(2),elem,isoelem
      logical change
      iel(1)=0
      iel(2)=0
      iso(1)=0
      iso(2)=0
      el(1)='  '
      el(2)='  '
      isoel(1)='  '
      isoel(2)='  '
      ion=0
      isostring='000   '
*
* atom or molecule ?
*
      do i=1,maxatom
        if(des(1:2).eq.lele(i)) then
          el(1)=lele(i)
          iel(1)=i
          if(des(2:2).eq.' ') then
*this is an atom/atomic ion
            read(des(3:3),'(i1)') ion
            ioffset=1
            goto 33
          else if(des(3:3).eq.' ') then
*this is an atom/atomic ion
            read(des(4:4),'(i1)') ion
            ioffset=2
            goto 33
          endif
        endif
      enddo
      goto 11
* since this is probably a molecule
   33 continue
* The following for atoms and atomic ions
* Now: isotope known?
      do i=1,26
        if(fullspecies(i:i).eq.'(') then
          if(fullspecies(i+2:i+2).eq.')') then
            read(fullspecies(i+1:i+1),'(i1)') iso(1)
            read(fullspecies(i+3:i+4),'(a2)') isoelem
            if(ioffset.eq.1) isoelem(2:2)=' '
            isoel(1)=isoelem
            goto 22
          else if(fullspecies(i+3:i+3).eq.')') then
            read(fullspecies(i+1:i+2),'(i2)') iso(1)
            read(fullspecies(i+4:i+5),'(a2)') isoelem
            if(ioffset.eq.1) isoelem(2:2)=' '
            isoel(1)=isoelem
            goto 22
          endif
        endif
      enddo
   22 continue
      if(iso(1).ne.0) then
        if(isoel(1).ne.el(1)) then
          print *,'iso,isoel',iso(1),isoel(1)
          stop 'Element confusion in subr. lineid'
        endif
        write(isostring(1:3),'(i3)') iso(1)
        if(isostring(1:1).eq.' ') isostring(1:1)='0'
        if(isostring(2:2).eq.' ') isostring(2:2)='0'
        if(iso(2).ne.0) then
          write(isostring(4:6),'(i3)') iso(2)
          if(isostring(4:4).eq.' ') isostring(4:4)='0'
          if(isostring(5:5).eq.' ') isostring(5:5)='0'
        endif
      endif
      return
*
   11 continue
* Determine (diatomic) molecule species
      onechar=des(2:2)
      call case(.false.,onechar,ichar2)
      onechar=des(3:3)
      call case(.false.,onechar,ichar3)
      onechar=des(4:4)
      call case(.false.,onechar,ichar4)
      onechar=des(5:5)
      call case(.false.,onechar,ichar5)
      onechar=des(6:6)
      call case(.false.,onechar,ichar6)
      if(ichar2.eq.3) then
* this is an isoatomic molecule, probably C2
        if(des(2:2).ne.'2') then
          ion=0
          return
        endif
        if(des(4:4).eq.'1') then
          ion=1
        else if(des(4:4).eq.'2') then 
          ion=2
        else
          ion=0
          return
        endif
        elem(1:1)=des(1:1)
        elem(2:2)=' '
        el(1)=elem
        el(2)=elem
      else if(ichar2.eq.1) then
        if(des(3:3).eq.' ') then
* this is another diatomic, e.g. CH
          if(des(4:4).eq.'1') then
            ion=1
          else if(des(4:4).eq.'2') then 
            ion=2
          else
            ion=0
            return
          endif
          elem(1:1)=des(1:1)
          elem(2:2)=' '
          el(1)=elem
          elem(1:1)=des(2:2)
          elem(2:2)=' '
          el(2)=elem
        endif
      else if(ichar2.eq.-1) then
        el(1)=des(1:2)
        if(des(4:4).eq.' ') then
* possibly MgH, TiO
          el(2)=des(3:4)
        else if(ichar4.eq.-1) then
* possibly NaCl
          if(des(5:5).ne.' ') then
            ion=0
            return
          else
            if(des(6:6).eq.'1') then
              ion=1
            else if(des(6:6).eq.'2') then 
              ion=2
            else
              ion=0
              return
            endif
          endif
        endif
* that exhausts the possibilities that I imagine in VALD 3
      endif
      do i=1,maxatom
        if(el(1).eq.lele(i)) iel(1)=i
        if(el(2).eq.lele(i)) iel(2)=i
      enddo
      end
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
      subroutine findatom(lel,lele,iel)
* get atomic number
      implicit none
      integer     i,iel
      character*2 lel,lele(92)
*
      do i=1,92
        if(lel.eq.lele(i)) then
          iel=i
          goto 12
        endif
      enddo
      print *,'element |',lel,'| does not exist'
      stop 'subr. findatom'
   12 continue
      return
      end
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
      subroutine define_transition(desl,desu,lower,upper)
* adapted to LS coupling for assigning van der Waals damping
* try to understand or guess what kind of transition this is
* desl, desu are level designations e.g. " (1G)sp u3F"
* lower, upper designate "l" quantum design. eg s p d f ..., "X" = unknown
* llower, lupper designate "l" quantum numb. eg 0 1 2 3 ..., -1 = unknown
* GUESS (the lowest possible "l") if one level is undefined
*
      implicit none
      integer llorbit(2),luorbit(2),llower,lupper
      character*91 desl,desu
      character*1 lower,upper
      lower='X'
      upper='X'
      if(desl(4:5).eq.'LS' .or. desl(4:5).eq.'JJ' .or. 
     &   desl(4:5).eq.'JK' .or. desl(4:5).eq.'LK') then
        call identify_level(desl,llorbit)
      else
        llorbit(1)=-1
        llorbit(2)=-1
      endif
      if(desu(4:5).eq.'LS' .or. desu(4:5).eq.'JJ' .or. 
     &   desu(4:5).eq.'JK' .or. desu(4:5).eq.'LK') then
        call identify_level(desu,luorbit)
      else
        luorbit(1)=-1
        luorbit(2)=-1
      endif
      if(luorbit(1).eq.-1 .and. llorbit(1).eq.-1) then
* no understandable orbitals found
        return
      endif
* define transition type if possibly possible
      llower=llorbit(1)
      lupper=luorbit(1)
* these are the defaults, the rest are exceptions
      if(llorbit(1).ge.3 .or. luorbit(1).ge.3) then
* Use the d-f tables (Ref: Barklem, private communication)
        llower=2
        lupper=3
      else if(llorbit(1).ge.0 .and. luorbit(1).ge.0 .and.
     &   (iabs(llorbit(1)-luorbit(1)).ne.1)) then
* if Delta(l) .ne. +-1 then try with the second alternatives
        if(llorbit(2).lt.0) then
* which of the 2 luorbits goes best with llorbit(1)
          if(luorbit(2).ge.0) then
            if((luorbit(1).eq.llorbit(1)) .or.
     &         (iabs(llorbit(1)-luorbit(2)).eq.1)) then
              lupper=luorbit(2)
            endif
          endif
        else
          if(luorbit(2).lt.0) then
* which of the 2 llorbits goes best with luorbit(1)
            if((llorbit(1).eq.luorbit(1)) .or.
     &         (iabs(luorbit(1)-llorbit(2)).eq.1)) then
              llower=llorbit(2)
            endif
          else
* llorbit 1 and 2 and luorbit 1 and 2 are all >=0
* find out which of the 3 remaining combinations make Delta(l)=+-1
* the ordering of these 3 alternatives is arbitrary (I've found no rule)
            if(iabs(llorbit(1)-luorbit(2)).eq.1) then
              lupper=luorbit(2)
            else
              if(iabs(llorbit(2)-luorbit(1)).eq.1) then
                llower=llorbit(2)
              else
                if(iabs(llorbit(2)-luorbit(2)).eq.1) then
* last chance: take both alternatives
                  llower=llorbit(2)
                  lupper=luorbit(2)
                endif
              endif
            endif
          endif
        endif
      else
* GUESS (the lowest possible l) if only one level is undefined:
        if(llower.lt.0 .and. lupper.ge.0) then
          if(lupper.eq.0) then
            llower=1
          else
            llower=lupper-1
          endif
        endif
        if(llower.ge.0 .and. lupper.lt.0) then
          if(llower.eq.0) then
            lupper=1
          else
            lupper=llower-1
          endif
        endif
      endif
* convert to letter labels
      if(llower.lt.0) then
        lower='X'
      else if(llower.eq.0) then
        lower='s'
      else if(llower.eq.1) then
        lower='p'
      else if(llower.eq.2) then
        lower='d'
      else if(llower.eq.3) then
        lower='f'
      else if(llower.eq.4) then
        lower='g'
      else if(llower.eq.5) then
        lower='h'
      else if(llower.eq.6) then
        lower='i'
      else if(llower.eq.7) then
        lower='k'
      else
        print *,'why is llower=',llower,' ???'
        stop 'subr. define_transition'
      endif
      if(lupper.lt.0) then
        upper='X'
      else if(lupper.eq.0) then
        upper='s'
      else if(lupper.eq.1) then
        upper='p'
      else if(lupper.eq.2) then
        upper='d'
      else if(lupper.eq.3) then
        upper='f'
      else if(lupper.eq.4) then
        upper='g'
      else if(lupper.eq.5) then
        upper='h'
      else if(lupper.eq.6) then
        upper='i'
      else if(lupper.eq.7) then
        upper='k'
      else
        print *,'why is lupper=',lupper,' ???'
        stop 'subr. define_transition'
      endif
      return
      end
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
      subroutine identify_level(des,lorb)
* try to tell whether this is an s, p, d, f, g, h or i orbital
* method:
* 1) separate into configuration and term (look for blank), if designation
*    does not consist of 2 or more parts separated by 1 blank: NO RESULT: -1 -1
* 2) throw away term, concentrate on configuration
* 3) find the last lower-case letter, that gives the orbital lorb(1)
*    if any character to the left is also a lower-case letter, that gives the lorb(2)
*    lorb < 0 means no result
*
      implicit none
      character*1 des(91)
      integer i,l,n,lorb(2)
*
      lorb(1)=-1
      lorb(2)=-1
      i=91
   11 continue
      i=i-1
      if(des(i).eq.' '.and.i.gt.4) then
        goto 11
      else if(des(i).eq.' '.and.i.eq.4) then
* no level designation at all given
        return
      endif
* OK: we've found the last character in the term
   12 continue
      i=i-1
      if(des(i).ne.' '.and.i.gt.4) then
        goto 12
      else if(des(i).ne.' '.and.i.eq.4) then
* no configuration data found
        return
      endif
* OK: we've found a blank before the possible term
      if(i.le.4) return
* now look for lower case letters in the configuration:
      n=1
   13 continue
        i=i-1
        if(des(i).eq.'s') then
          l=0
        else if(des(i).eq.'p') then
          l=1
        else if(des(i).eq.'d') then
          l=2
        else if(des(i).eq.'f') then
          l=3
        else if(des(i).eq.'g') then
          l=4
        else if(des(i).eq.'h') then
          l=5
        else if(des(i).eq.'i') then
          l=6
        else if(des(i).eq.'k') then
          l=7
        else
          l=-1
          if(i.gt.5) goto 13
        endif
        lorb(n)=l
      if(n.eq.2) then
        return
      else
        n=n+1
        if(i.gt.5) goto 13
      endif
      return
* lorb(1) and lorb(2) contain the two last l quantum numbers
      end
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
      subroutine shrinkdesignation(string,coupling,desig)
* left-adjust the level designation and remove multiple blanks
      implicit none
      character string*91,coupling*2,desig*91
      integer   i,j,first,last,length
      logical remove
      write(desig,'(91x)')
      coupling=string(4:5)
      if(coupling.eq.'  ') coupling='__'
      last=len_trim(string)
      if(last.gt.6) then
        do j=6,last
          if(string(j:j).ne.' ') then
            first=j
            exit
          endif
        enddo
        write(desig,'(a)') string(first:last)
      endif
* remove multiple blanks
      length=len_trim(desig)
   10 continue
      remove=.false.
      do i=2,length
        if(desig(i:i).eq.' '.and.desig(i+1:i+1).eq.' ') then
          do j=i+2,length
            desig(j-1:j-1)=desig(j:j)
          enddo
          desig(length:length)=' '
          remove=.true.
        endif
      enddo
      if(remove) then
        length=length-1
        goto 10
      endif
      return
      end
*
*
*...v....1....v....2....v....3....v....4....v....5....v....6....v....7..
*
*
