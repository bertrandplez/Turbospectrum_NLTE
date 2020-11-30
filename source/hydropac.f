      subroutine hydropac(lunit,xlboff,nlte,
     &                    nlte_species,maxlevel,modnlevel,b_departure,
     &                    modenergy,modg,modion,modid)
!
! include NLTE capability with departure coefficients as input. BPz 17/11-2020
!
! new version for use with hbop.f  ! H line + continuum bf opacity
! Bpz 02/02-2019
! BPz: inclusion of nlte case. Necessary regardless of case for HI lines, as 
! we need emissivity (stored in source_function), and not only extinction coefficient.

      implicit none

      include 'spectrum.inc'

      character comment*100,species*20
      integer ion,nline,ntau,iter,l,k,i,lstart,lstartp,
     &        nmy,nlbldu,iint,iweak,lunit,maxlam,ll,
     ,        jlcont, nlcont,nf,nfo,nb,nbo,nl,lpos
      doubleprecision ionpot,wave,xl1,xl2,del,xlboff,xlambda
      doubleprecision source_function

! 150 is the max allowed number of H lines
      doubleprecision hlambda(150)
      real xlo(150),xup(150),gf(150),npop(150)
      real cont,total,contrib
      integer nlo(150),nup(150)
      character*9 lname(150)
!
      real fpart(ndp),pe,t,pg,xi,mum,ro,eps,xmyc,
     &     ne(ndp),
     &     nh1(ndp),nhe1(ndp),hnorm,hnormnostim(ndp),stim,
     &     diff,diffp,theta(ndp),hckt(ndp),abso,absos,absocont,
     &     dopple(ndp),xkapr,cross,hlinop,h1bfg,alpha,
     ,     ee(6),d0(ndp),d(lpoint,ndp),absoscont
      logical  lymanalpha, usedam,notfound,contonly,kskip(ndp),lskip
      logical  ldone(lpoint),lineonly
! NLTE
      logical nlte,nlte_species
      real xlsingle,bpl,corr,expcorr,bdl,bdu,ediff
      character*40 modid(*)
      integer modnlevel,modion(*),modnlev,nlevlist,maxlevel
      real b_departure(ndp,0:maxlevel),modenergy(*),modg(*),bd(1000)
      real*8 planck,clight,electron,boltzmann
      real kayser2eV,ev2kayser,hck,kbcgs
*      
      common /atmos/ t(ndp),pe(ndp),pg(ndp),xi(ndp),mum(ndp),ro(ndp),
     &               ntau
      doubleprecision presneutral,presion,presion2,presion3, h1bfgc(30)    
      common /orderedpress/ presneutral(ndp,100),presion(ndp,100),
     &                      presion2(ndp,100),presion3(ndp,100)
      common /pieces/ xl1,xl2,del,eps,nmy,nlbldu,iint,xmyc,iweak
      common/large/ xlambda(lpoint),source_function(ndp,lpoint),
     & maxlam,abso(ndp,lpoint),
     & absos(ndp,lpoint),absocont(ndp,lpoint),absoscont(ndp,lpoint)
      common/rossc/ xkapr(ndp),cross(ndp)
      common/continuum/nlcont,jlcont(lpoint)
      
      character*50 zMC
      integer znlc ,  ii
      real    zxlm,zBP(ndp),zXC(ndp),zS(ndp),zXI(ndp)

      data clight /2.99792458d8/, planck /6.62607015d-34/
      data boltzmann /1.380649d-23/, electron /1.602176634d-19/

!      kayser2ev = 0.0001239854563934909
!      ev2kayser = 8065.4621040899665
      kbcgs=sngl(boltzmann*1.d7)
      kayser2eV = sngl(planck * clight / electron * 1.d2)
      ev2kayser=1./kayser2ev
      hck=sngl(planck * clight / boltzmann * 1.d10)     ! ready for lambda in AA !
! 
      if (nlte_species) then
        do i=1,modnlevel-1
          if (modion(i).ne.1) then
            print*,'hydropac: problem with HI levels'
            stop 'stop in hydropac'
          endif
        enddo
        if (modion(modnlevel).ne.2) then
          print*,'HII level missing in model atom !'
          stop 'stop in hydropac'
        endif
      endif
      modnlev=max(modnlevel-1,0)   ! modnlev=0 if LTE, or =number of HI levels if NLTE

      ldone=.false.
! set populations to 0. (computed as LTE pop in hbop)
      npop=0.0
! compute lines only 
      contonly=.false.
      lineonly=.true.
! read line list
      rewind (lunit)
      read(lunit,*) species,ion,nline
      read(lunit,*) comment
*      print*,species,ion,nline,comment
      if (species(1:6).ne.'01.000'.or.ion.ne.1) then
        print*, 'wrong H line data file!'
        print*,species,ion
        stop 'ERROR!'
      endif
c bsyn uses air wavelengths, and so does hbop.f, except for Lyman series
c      
      nlevlist=0
      do i=1,nline
!       oscillator strengths are computed in hbop.f
        read(lunit,*) hlambda(i),nlo(i),nup(i),xlo(i),xup(i),gf(i),
     &                lname(i)
        nlevlist=max(nlevlist,nup(i))               ! number of levels in line list
! check levels in line list vs levels in model atom !
        if (modnlevel.gt.0) then
          if (nlo(i).le.modnlev) then
            ediff = abs(xlo(i)*ev2kayser - modenergy(nlo(i)))/
     &               modenergy(nlo(i))
            if (ediff.gt.1.e-3) then
              print*,'H level identification problem'
              print*,'model atom ', modenergy(nlo(i)), 
     &               'line list',xlo(i),xlo(i)*ev2kayser
              stop 'hydropac, stopping'
            endif
          endif
          if (nup(i).le.modnlev) then
            ediff = abs(xup(i)*ev2kayser - modenergy(nup(i)))/
     &               modenergy(nup(i))
            if (ediff.gt.1.e-3) then
              print*,'H level identification problem'
              print*,'model atom ', modenergy(nup(i)), 
     &               'line list',xup(i),xup(i)*ev2kayser
              stop 'hydropac, stopping'
            endif
          endif
        endif
!
      enddo
*
      do k=1,ntau
        ne(k)=pe(k)/(t(k)*kbcgs)
        nh1(k)=sngl(presneutral(k,1))/(t(k)*kbcgs)
        nhe1(k)=sngl(presneutral(k,2))/(t(k)*kbcgs)
* hckt=hc/kT ready for lambda in AA.
        hckt(k)=hck/T(k)
        dopple(k)=sqrt( xi(k)**2 * 1. + 
     &                 2.*kbcgs*t(k)/1.6738e-24) /
     &                 sngl(clight*1.d2)
      enddo
! add opacity
      do nl=1,nline
! find out where we have contributing lines
! search if line lies within the interval 
        print*,'H line ',hlambda(nl)
        if (hlambda(nl)-xlambda(1).le.0.) then
! line lies blueward
          lpos=1
        else if (hlambda(nl)-xlambda(maxlam).ge.0.) then
          lpos=maxlam
        else
! locate it by successive splitting in halves
          nb=1
          nbo=nb
          nf=maxlam
          nfo=nf
          do while (nf-nb.gt.1)
            if (hlambda(nl)-xlambda(nb).gt.0.) then
              nbo=nb
              nb=(nb+nf)/2
            else
              nb=nbo
              nf=(nb+nf)/2
            endif
            lpos=nb
          enddo
        endif
!check
        print*,'Hline ',hlambda(nl),lpos,xlambda(lpos),xlambda(lpos+1)
! now compute line profile
        lskip=.false.
        kskip=.false.
        do l=lpos,maxlam
          if (.not.ldone(l).and..not.lskip) then
!  	    print*,' l,xlam ',l,sngl(xlambda(l))
            lskip=.true.
            do k=1,ntau
! departure coefficients for that depth
              if (nlte_species) then
                do i=1,modnlev
                  bd(i)=b_departure(k,i)
                enddo
                do i=max(modnlev+1,1),nlevlist
                  bd(i)=1.0
                enddo
              endif
              if (.not. kskip(k)) then
! include selected H line at that wavelength
!!                call hbop(xlambda(l),nline,nlo,nup,hlambda,
                call hbop(xlambda(l),1,nlo(nl),nup(nl),hlambda(nl),
     &           nh1(k),nhe1(k),ne(k),t(k),dopple(k),npop,
     &           bd,modnlev,
     &           total,cont,contonly,lineonly)
                contrib = (total - cont)/xkapr(k)/ro(k)
                abso(k,l) = abso(k,l) + contrib

! correct source function for the NLTE case. BPz 17/11-2020
! source_function at this stage is the emissivity. It must be divided
! by the opacity to yield the source function.

                if (nlte) then
                  xlsingle=sngl(xlambda(l))
                  if (nlte_species) then
                    expcorr = exp(hckt(k)/xlsingle)
                    corr = (expcorr - 1.)/
     &               ( bd(nlo(nl)) / bd(nup(nl)) * expcorr - 1. )
                    if (k.eq.10) then
                      print*,'check',l,
     &                   expcorr,bd(nlo(nl)),bd(nup(nl)),corr
                    endif
                  else
                    corr=1.0
                  endif
                  source_function(k,l) = source_function(k,l) + 
     &                                 contrib*bpl(T(k),xlsingle)*corr
                endif

! HI bf is already included in babsma.f
!!!!             absocont(k,l) = absocont(k,l) + cont/xkapr(k)/ro(k)
!                print*,xlambda(l),k,contrib,absocont(k,l)
                if(contrib/absocont(k,l) .le. eps) then
                  kskip(k)=.true.
                endif
              endif
! if kskip is .true. at all depths, then lskip becomes .true., and next 
! wavelength is skipped
              lskip=lskip.and.kskip(k)
            enddo
! mark this wavelength as done
!            ldone(l)=.true.
          endif
        enddo
! other side of the profile
        lskip=.false.
        kskip=.false.
        do l=lpos-1,1,-1
          if (.not.ldone(l).and..not.lskip) then
*           print*,' l,xlam ',l,xlambda(l)
            lskip=.true.
            do k=1,ntau
! departure coefficients for that depth
              if (nlte_species) then
                do i=1,modnlev
                  bd(i)=b_departure(k,i)
                enddo
                do i=max(modnlev+1,1),nlevlist
                  bd(i)=1.0
                enddo
              endif
              if (.not. kskip(k)) then
! include selected H line at that wavelength
!!                call hbop(xlambda(l),nline,nlo,nup,hlambda,
                call hbop(xlambda(l),1,nlo(nl),nup(nl),hlambda(nl),
     &           nh1(k),nhe1(k),ne(k),t(k),dopple(k),npop,
     &           b_departure,modnlev,
     &           total,cont,contonly,lineonly)
                contrib = (total - cont)/xkapr(k)/ro(k)
                abso(k,l) = abso(k,l) + contrib

! correct source function for the NLTE case. BPz 17/11-2020
! source_function at this stage is the emissivity. It must be divided
! by the opacity to yield the source function.

                if (nlte) then
                  xlsingle=sngl(xlambda(l))
                  expcorr = exp(hckt(k)/xlsingle)
                  corr = (expcorr - 1.)/
     &               ( bd(nlo(nl)) / bd(nup(nl)) * expcorr - 1. )
                  source_function(k,l) = source_function(k,l) + 
     &                                 contrib*bpl(T(k),xlsingle)*corr
                  if (k.eq.10) then
                    print*,'check',l,
     &                 expcorr,bd(nlo(nl)),bd(nup(nl)),corr
                  endif
                endif

! HI bf is already included in babsma.f
!!!!             absocont(k,l) = absocont(k,l) + cont/xkapr(k)/ro(k)
                if(contrib/absocont(k,l) .le. eps) then
                  kskip(k)=.true.
                endif
              endif
              lskip=lskip.and.kskip(k)
            enddo
! mark this wavelength as done
!            ldone(l)=.true.
          endif
        enddo
! end of loop on H lines
      enddo

      return
      end
