      program interpol_modeles

 
****************************************************************************
*     interpolation of model atmosphere + NLTE departure coefficients
*       
* parameter space of interpolation {Teff,logg,z}                  
* 8 MARCS binary models as input required: 
* ! the order of input models matters !
* ! only standard composition has been tested, no guaranty for peculiar composition ! 
* turbospectrum/babsma format compatible output 

* Interpolation scheme of model atmosphere(@)
* in {Teff,logg,z} stellar parameter space 
*         ^ z
*         |_ _ _ _ _
*         |        /| 
*        /|       / 
*       _ |_ _ _ /  |
*       | |      |  
*         |         |  
*       | |    @ |
*         |--------------> logg
*       |/       | / 
*       /_ _ _ _ ./
*      /          
*     /
*    Teff
*
* Each structure component of each model (T,Pe,Pg,xit,optical depth, kappaross)
* is first resampled on a common optical depth basis tau(5000 or ross) (see resample routine).
*
* Then each component of the atmosphere(@) (T,Pe,Pg,xit,kappa,geom depth +tau5000 and tauross)) 
* is interpolated individually along each dimension between each input model value(*) 
* (physically it is impossible to ensure a simple relationship between models in more than one dimension at the same time, 
* that's why the input models MUST form a "cube" in the stellar parameter space {Teff,logg,z})
* The interpolation is successively done at each optical depth (tau)
* + interpolation point weighted by an empirical value (see TM thesis + manual) 
*
*
*                 ^ T,Pe,Pg,xit or radius 
*                 |
*                 |            *
*                 |           /|
*                 |          @
*                 |         /| | 
*                 |   /    /   
*                 |  /    /  | | 
*         ^       | /    *   
*         |       |/     |   | |
*         |       -----------------------> Teff, logg and z 
*         |      /      low ref up             
*         |     /*      /   / /
*         |    / |\           
*         |   /    \  /   / / 
*         |  /   |  \     
*         | /       /@  / /  
*         |/     |   |\      
*         / _  _ /_ _ /*_ _ _ _ _    
*        /      low ref up            
*       /    
*   tau{5000 or Ross}
*
***************************************************************************
* TM 07/2003  

c  07/2004 resampling of each model on a common optical depth basis                                     
c  06/2006 works for spherical geometry models  
c  09/2007 new calibration of free parameter alpha
c           + modified to read Uppsala ascii models 
c           + kappa interpolation 
c           + rhox calculation 
c           + 2 outputs (babsma and ATLAS/MOOG)
c  10/2007 error estimates
c  10/2011 two non crucial bugs fixed 
c             -unformatted->formatted read for ascii models because there is a couple of trouble makers in the grid  
c             -dimension of taubas was not matching tauR (emo@astro.indiana.edu)          
c  04/2012 unformatted reading of ascii models reinstated (troublemakers hopefully fixed) /BE
****************************************************************************
!  compile with Fortran 90 or 95  
c
c    MB: 2020 added a new module to interpolate in the grid of NLTE departure coefficients
c      
c     use IFPORT
      implicit none
      integer :: file,nfile,k,n,m,ndp,ndepth_ref,out,nlinemod
      parameter (ndp=200)
      parameter (nfile=11)

      integer :: nlte_file, icomp
      
! MB: these TWO PARAMETERS HAVE TO BE SET BY HAND   
c      parameter (nlte_file=14499)
      parameter (nlte_file=14627) !updated by the user. is the length of the number in lines in the auxiliary nlte data file (which is the number of models)
c      parameter (nlte_file=188)
c      parameter (nlte_file=3)
      parameter (icomp=8)
      
      logical :: verif,check,test,extrapol,binary,optimize
      real :: lambda_ref,temp_ref,logg_ref,z_ref,x,y,z,xinf,
     &        xsup,yinf,ysup,zinf,zsup,teffpoint,loggpoint,metpoint
      character*256, dimension (nfile) :: FILE_IN
      character*256 :: FILE_ALT
      real, dimension (:,:), allocatable:: taus,tauR,T,Pe,Pg,xit,rr,
     &        xkapref,rhox,taus_aux,tauR_aux,T_aux,Pe_aux,Pg_aux,
     &        xit_aux,rr_aux,xkapref_aux 
      integer, dimension (nfile) :: ndepth
      real, dimension (nfile) :: xlr,teff,logg,metal
      logical, dimension (nfile) :: sph
      external :: blend_103
      real, external :: inf,sup
      real, dimension(8,3) :: lin_dif,power
******MB      
      character*15, dimension (8) :: coefval
      character*256 :: nlte_binary, nlte_model_list
      integer :: n_dep, n_lev, cnt, cnt1
      character*500, dimension (nlte_file) :: id_model
      character*1000 :: n_comment
      real*8, allocatable :: nlte_tau(:), nlte_data(:,:), 
     &     bvals_aux(:,:,:), bvals(:,:,:), nnbvals(:,:,:),
     &     nntau_aux(:,:), nntau(:,:)
      real, dimension (nlte_file) :: n_logg,n_metal,n_teff,n_abu
      real, dimension (nlte_file) :: nnlogg,nnmetal, nnteff
      integer, dimension (nlte_file) ::  index
      integer, dimension (nlte_file) :: n_dummy
      integer*8, dimension (nlte_file) :: n_pos1, nnpos
      real :: n_result
      integer*4           ::  tmp(1)
      character*65, dimension (nlte_file) :: n_pos
      integer*8 :: int1, stat1, pos1

      real y000
      real y001
      real y010
      real y011
      real y100
      real y101
      real y110
      real y111

******ES/MB: printing the binary files      
      CHARACTER*350 BinaryName, formatBin
      CHARACTER*1000 BinaryHeader
      CHARACTER*10 ABNDstring
      
******MB      
      DATA coefval/'tau500', 'tauross', 'T', 'logpe',  
     &             'logpg', 'xit', 'rr', 'logkapref'/
      
      INTERFACE reec
        subroutine resample(taus,tauR,T,Pe,Pg,xit,rr,xkapref)
        real, dimension (:,:) :: taus,tauR,T,Pe,Pg,xit,rr,xkapref
        end 
      END INTERFACE reec

      write(*,*) '*****************************'
      write(*,*) '* begining of interpolation *'
      write(*,*) '*****************************'  

******* you can choose here to switch of the "optimization" and prefer simple linear interpolation
 
      optimize = .true.

******  read 8 models, put in tables, 
****** check number of layer, reference optical depth, and geometry compatibility ******

      out=9

******MB: loop to 11 : new string added (interpolation coefficients)      
      do file=1,11
         read(*,*) FILE_IN(file)
      end do
      
******MB: read a binary file with NLTE departures
       read(*,*) nlte_binary
       write(*,76) nlte_binary
 76    format('NLTE binary: ', a200)

       read(*,*) nlte_model_list
       write(*,77) nlte_model_list
 77    format('NLTE model list: ', a200)       
       
       read(*,*) temp_ref
       read(*,*) logg_ref
       read(*,*) z_ref
       read(*,*) test

      verif=.true.
      check=.true.
      nlinemod=ndp
      
      allocate(taus_aux(nlinemod,nfile),tauR_aux(nlinemod,nfile),
     & T_aux(nlinemod,nfile),Pe_aux(nlinemod,nfile),
     & Pg_aux(nlinemod,nfile),xit_aux(nlinemod,nfile),
     & rr_aux(nlinemod,nfile),xkapref_aux(nlinemod,nfile))

******MB: read model atmopsheres (only 8 models allowed!)
      read(*,*) binary
      if (binary) then
      do file=1,8
       call extract_bin(FILE_IN(file),teff(file),logg(file),metal(file),
     &  ndepth(file),xlr(file),taus_aux(:,file),tauR_aux(:,file),
     &  T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),xit_aux(:,file),
     &  rr_aux(:,file),sph(file),xkapref_aux(:,file))
        verif=verif.and.(ndepth(file).eq.ndepth(1))
        check=check.and.(xlr(file).eq.xlr(1))
        if (.not.(((sph(1).and.sph(file))).or.
     &    ((.not.(sph(1))).and.(.not.(sph(file)))))) then
         write(*,*) 'geometry compatibility problem with'
         write(*,78) file,teff(file),logg(file),metal(file)     
         stop
        endif  
        write(*,78) file,teff(file),logg(file),metal(file)
      end do
      else
      do file=1,8
        call extract_ascii(FILE_IN(file),teff(file),logg(file),
     &  metal(file),ndepth(file),xlr(file),taus_aux(:,file),
     &  tauR_aux(:,file),T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),
     &  xit_aux(:,file),rr_aux(:,file),sph(file),xkapref_aux(:,file))
        verif=verif.and.(ndepth(file).eq.ndepth(1))
        check=check.and.(xlr(file).eq.xlr(1))
        if (.not.(((sph(1).and.sph(file))).or.
     &    ((.not.(sph(1))).and.(.not.(sph(file)))))) then
        write(*,*) 'geometry compatibility problem with'
        write(*,78) file,teff(file),logg(file),metal(file)     
        stop
       end if   
       write(*,78) file,teff(file),logg(file),metal(file)
      end do
      endif
 78   format('model',i2,'  Teff=',f8.0,'  logg=',f5.2,'  z=',f6.2)

*********MB:read list of model atmospheres for which NLTE departures are available
*********

      open(unit=199, file=nlte_model_list,  form='formatted')
      read(199,*) n_comment
      read(199,*) n_comment
      read(199,*) n_comment

      do cnt=1, nlte_file
         read(199, 79)
     &        id_model(cnt), n_teff(cnt), n_logg(cnt), n_metal(cnt),
     &        n_abu(cnt), n_pos(cnt)
         call str2int(n_pos(cnt),int1,stat1)
         n_pos1(cnt) = int1
c         print*, n_metal(cnt), metal(1)
c        print*, n_pos(cnt), n_pos1(cnt)
      enddo
      
 79   format(a71,2x,f7.1,6x,f5.2,6x,f5.2,4x,f5.2,a65)
      close(199)
      
*********MB: cross-correlate the parameters with those in the input file

      do cnt = 1, 8
        do cnt1 = 1, nlte_file
           if  (teff(cnt).eq.n_teff(cnt1).and.
     &         logg(cnt).le.(n_logg(cnt1)+0.00001).and.
     &         logg(cnt).ge.(n_logg(cnt1)-0.00001).and.
     &         metal(cnt).eq.n_metal(cnt1)) then
               index(cnt) = cnt1
             write(*,*) teff(cnt), logg(cnt), metal(cnt), index(cnt)
            exit
           endif
         enddo
      enddo   

      write(*,*)
     &    'new NLTE array: Teff, Logg, [Fe/H], position / binary'
         
*********MB create a new NLTE array

      do cnt = 1, 8
         nnteff(cnt) = n_teff(index(cnt))
         nnlogg(cnt) = n_logg(index(cnt))
         nnmetal(cnt) = n_metal(index(cnt))
         nnpos(cnt) = n_pos1(index(cnt))
         write(*,81) 'new NLTE array: ',
     &   nnteff(cnt), nnlogg(cnt), nnmetal(cnt), nnpos(cnt)
 81      format(a17, f5.0, f5.2, 2x, f5.2, 2x, i10)
      enddo


*********MB:read NLTE binary
********* do not read the entire binary file, but use pointers set by AUX
********* the pointers were re-ordered to match (T,g,m) of the input MARCS binary models      
      
      open(unit=200,file=nlte_binary,form='unformatted',access="stream")
*     read(200) n_comment
*     rewind(200)

      do cnt=1, 8
        pos1 = nnpos(cnt)

        read(200, pos = pos1) id_model(cnt)
        write(*,fmt='("NLTE: ", i1,1x,A50,1x, i10)')
     &        cnt, trim(adjustl(id_model(cnt))), pos1

        pos1 = pos1 + 500
        read(200, pos = pos1) tmp
        n_dep = tmp(1)
        pos1 = pos1 + icomp
        read(200, pos = pos1) tmp
        n_lev = tmp(1)

        if (cnt.eq.1)  allocate(bvals_aux(n_dep, n_lev, 8))

        if (cnt.eq.1)  allocate(nntau_aux(n_dep,8))
        
        pos1 = pos1 + icomp
        allocate(nlte_tau(n_dep))
        read(200, pos = pos1), nlte_tau
        nntau_aux(:,cnt) = nlte_tau
        deallocate(nlte_tau)

        pos1 = pos1 + n_dep*icomp
        allocate(nlte_data(n_dep, n_lev))
        read(200, pos = pos1), nlte_data
        bvals_aux(:,:,cnt) = nlte_data
        deallocate(nlte_data)
        rewind(200)
      end do
*      stop

*********MB: print the selected departure values into a text file
      
      open(unit=28, file='test_nlte_input.txt')
      do cnt=1,8
        write(28, *) '*'
        write(28, *) trim(adjustl(id_model(cnt)))
        do n=1,n_dep
         write(28,*)  (bvals_aux(n,m,cnt), m=1, n_lev)
       enddo
      enddo
      close(28)
      
      print*, 'NLTE file:  8 models read'  
      close(200)

      allocate(bvals(n_dep, n_lev, icomp))
      bvals = bvals_aux(1:n_dep,:,:)
      
      allocate(nntau(n_dep, icomp))
      nntau = nntau_aux(1:n_dep,:)
      
*********MB:re-order NLTE departure files - NO NEED TO REORDER ANYMORE -
*
      allocate(nnbvals(n_dep, n_lev, out))
      do cnt = 1, 8
*         nnbvals(:,:,cnt) = bvals(:,:,index(cnt))
          nnbvals(:,:,cnt) = bvals(:,:,cnt)
         if (cnt.eq.1) then
c          write(*,*) 'new NLTE array: ', nnbvals(:,1,1)
         endif
      enddo
      
      
*********MB:print checks
      
      write(*,80) temp_ref,logg_ref,z_ref
 80   format('Interpolation point : Teff=',f8.0,'  logg=',f5.2,
     &     '  z=',f6.2)
      
**************check if files are length and depth ref compatible *******
      
      if (.not.(check)) then
         write(*,*) 'All the models do not have the same' 
          write(*,*) 'lambda ref'
          write(*,*) 'no interpolation done'
          stop
       else    
      if (.not.(verif)) then
         write(*,*) 'WARNING : All the models do not have the same' 
         write(*,*) 'number of layers, resampling to',
     &                 ndepth(1),'layers'
      end if
      ndepth_ref=ndepth(1)
      lambda_ref=xlr(1)


*********MB:check if NLTE n_depths are consistent with the input MARCS binary models
*********only  the last entry in the binary file is checked (suboptimal)
      
      if (ndepth_ref.ne.n_dep) then
       write(*,*) 'ERROR: NLTE departures are on a different tau scale!'
      else if (ndepth_ref.eq.n_dep) then
       write(*,*) '*** # depths (NLTE, MARCS) consistent'
      endif
      
********* calculation of the interpolation point(x,y,z) in {Teff,logg,z} space******************
      
       allocate(taus(ndepth_ref,nfile),tauR(ndepth_ref,nfile),
     & T(ndepth_ref,nfile),Pe(ndepth_ref,nfile),Pg(ndepth_ref,nfile),
     & xit(ndepth_ref,nfile),rr(ndepth_ref,nfile),
     &     xkapref(ndepth_ref,nfile))
       
       taus=taus_aux(1:ndepth_ref,:)
       tauR=tauR_aux(1:ndepth_ref,:)
       T=T_aux(1:ndepth_ref,:)
       Pe=Pe_aux(1:ndepth_ref,:)
       Pg=Pg_aux(1:ndepth_ref,:)
       xit=xit_aux(1:ndepth_ref,:)
       rr=rr_aux(1:ndepth_ref,:)
       xkapref=xkapref_aux(1:ndepth_ref,:)

         xinf=inf(teff)
         yinf=inf(logg)
         zinf=inf(metal)
         xsup=sup(teff)
         ysup=sup(logg)
         zsup=sup(metal)
         if (xsup.eq.xinf) then
            teffpoint=0
         else
          teffpoint=(temp_ref-xinf)/(xsup-xinf)
         end if 
         if (ysup.eq.yinf) then
            loggpoint=0
         else      
          loggpoint=(logg_ref-yinf)/(ysup-yinf)
         end if
         if (zsup.eq.zinf) then
            metpoint=0
         else   
          metpoint=(z_ref-zinf)/(zsup-zinf)
         end if
         extrapol=((teffpoint.lt.0).or.(teffpoint.gt.1)
     &         .or.(loggpoint.lt.0).or.(loggpoint.gt.1)
     &         .or.(metpoint.lt.0).or.(metpoint.gt.1))
         if (extrapol) then
         write(*,*) '!!!  WARNING : extrapolation  !!!'
         end if

*     
*******resample each layer of each input model on a common depth basis(tau5000 or tauRoss, see resample routine)*****************
!if you don't want to resample all the model to the same depth scale, just comment the following line  
         call resample(taus,tauR,T,Pe,Pg,xit,rr,xkapref)


****** initialisation of empirical constants for optimized interpolation (see TM thesis)*************
         
         lin_dif(1,1)=0                          !tau5000 vs Teff
         lin_dif(1,2)=0                          ! ...    vs logg
         lin_dif(1,3)=0                          ! ...    vs z    
         lin_dif(2,1)=0                          !tauross vs Teff
         lin_dif(2,2)=0                          ! ...    vs logg 
         lin_dif(2,3)=0                          ! ...    vs z    
         lin_dif(3,1)=0.15                       !T       vs Teff
         lin_dif(3,2)=0.3                        ! ...    vs logg 
         lin_dif(3,3)=1-(temp_ref/4000)**2.0     ! ...    vs z    
         lin_dif(4,1)=0.15                       !logPe   vs Teff
         lin_dif(4,2)=0.06                       ! ...    vs logg 
         lin_dif(4,3)=1-(temp_ref/3500)**2.5     ! ...    vs z    
         lin_dif(5,1)=-0.4                       !logPg   vs Teff
         lin_dif(5,2)=0.06                       ! ...    vs logg  
         lin_dif(5,3)=1-(temp_ref/4100)**4       ! ...    vs z    
         lin_dif(6,1)=0                          !xit     vs Teff
         lin_dif(6,2)=0                          ! ...    vs logg  
         lin_dif(6,3)=0                          ! ...    vs z    
         lin_dif(7,1)=0                          !rr      vs Teff
         lin_dif(7,2)=0                          ! ...    vs logg  
         lin_dif(7,3)=0                          ! ...    vs z    
         lin_dif(8,1)=-0.15                      !logxkapref vs Teff
         lin_dif(8,2)=-0.12                      ! ...    vs logg  
         lin_dif(8,3)=1-(temp_ref/3700)**3.5     ! ...    vs z    
       
         if (optimize) then
          write(*,*) 'optimized interpolation applied for standard compo
     &sition models'
         else
            lin_dif=0.
             write(*,*) 'linear interpolation applied'
         end if   
!these constants are calibrated on a broad range of stellar parameters; scale them now to the present one.
            power(:,1)= 1-(lin_dif(:,1)*(abs(xsup-xinf)/(7000-3800)))
            power(:,2)= 1-(lin_dif(:,2)*(abs(ysup-yinf)/(5-0.0)))
            power(:,3)= 1-(lin_dif(:,3)*(abs(zsup-zinf)/(0-(-4))))
                  
****** interpolation of each component of the atmosphere (taus,teff,Pe,Pg,microt,rr) and at each layer *****************


        do k=1,ndepth_ref
          x=(teffpoint)**power(1,1)
          y=(loggpoint)**power(1,2)
          z=(metpoint)**power(1,3)
          call blend_103(x,y,z,taus(k,1),taus(k,2),
     &     taus(k,3),taus(k,4),taus(k,5),taus(k,6),taus(k,7),taus(k,8)
     &     ,taus(k,out))

          x=(teffpoint)**power(2,1)
          y=(loggpoint)**power(2,2)
          z=(metpoint)**power(2,3)
          call blend_103(x,y,z,tauR(k,1),tauR(k,2),
     &     tauR(k,3),tauR(k,4),tauR(k,5),tauR(k,6),tauR(k,7),tauR(k,8)
     &     ,tauR(k,out))
          
          x=(teffpoint)**power(3,1)
          y=(loggpoint)**power(3,2)
          z=(metpoint)**power(3,3)
          call blend_103(x,y,z,T(k,1),T(k,2),T(k,3),T(k,4)
     &     ,T(k,5),T(k,6),T(k,7),T(k,8),T(k,out))
          
          x=(teffpoint)**power(4,1)
          y=(loggpoint)**power(4,2)
          z=(metpoint)**power(4,3)
          call blend_103(x,y,z,Pe(k,1),Pe(k,2),Pe(k,3),Pe(k,4)
     &     ,Pe(k,5),Pe(k,6),Pe(k,7),Pe(k,8),Pe(k,out))
          
          x=(teffpoint)**power(5,1)
          y=(loggpoint)**power(5,2)
          z=(metpoint)**power(5,3)
          call blend_103(x,y,z,Pg(k,1),Pg(k,2),Pg(k,3),Pg(k,4)
     &     ,Pg(k,5),Pg(k,6),Pg(k,7),Pg(k,8),Pg(k,out))
          
          x=(teffpoint)**power(6,1)
          y=(loggpoint)**power(6,2)
          z=(metpoint)**power(6,3)
          call blend_103(x,y,z,xit(k,1),xit(k,2),
     &     xit(k,3),xit(k,4),xit(k,5),xit(k,6),xit(k,7),xit(k,8)
     &     ,xit(k,out))       
 
          x=(teffpoint)**power(7,1)
          y=(loggpoint)**power(7,2)
          z=(metpoint)**power(7,3)
          call blend_103(x,y,z,rr(k,1),rr(k,2),
     &     rr(k,3),rr(k,4),rr(k,5),rr(k,6),rr(k,7),rr(k,8)
     &     ,rr(k,out))        

          x=(teffpoint)**power(8,1)
          y=(loggpoint)**power(8,2)
          z=(metpoint)**power(8,3)
          call blend_103(x,y,z,xkapref(k,1),xkapref(k,2),
     &     xkapref(k,3),xkapref(k,4),xkapref(k,5),xkapref(k,6),
     &         xkapref(k,7),xkapref(k,8),xkapref(k,out))

c          write(*,fmt="(i2, 9(f10.5,2x))") k,
c     &     xkapref(k,1),xkapref(k,2),
c     &     xkapref(k,3),xkapref(k,4),xkapref(k,5),xkapref(k,6),
c     &     xkapref(k,7),xkapref(k,8),xkapref(k,out)
       end do
       ndepth(out)=ndepth_ref
       xlr(out)=lambda_ref

**********MB: interpolate NLTE departures using the rr (7) coeffs

      do k=1,  ndepth_ref
        do n=1, n_lev
           x=(teffpoint)**power(7,1)
           y=(loggpoint)**power(7,2)
           z= (metpoint)**power(7,3)

           y000 = nnbvals(k,n,1)
           y001 = nnbvals(k,n,2)
           y010 = nnbvals(k,n,3)
           y011 = nnbvals(k,n,4)
           y100 = nnbvals(k,n,5)
           y101 = nnbvals(k,n,6)
           y110 = nnbvals(k,n,7)
           y111 = nnbvals(k,n,8)

           call blend_103(x,y,z,
     &      y000, y001, y010, y011,
     &      y100, y101, y110, y111, n_result)
           
           nnbvals(k,n,out)=n_result
c           if (n.eq.1) then 
c           write(*,fmt="(i2, 9(f10.5,2x))") k,
c     &      nnbvals(k,n,1), nnbvals(k,n,2),nnbvals(k,n,3), 
c     &      nnbvals(k,n,4), nnbvals(k,n,5),nnbvals(k,n,6), 
c     &      nnbvals(k,n,7), nnbvals(k,n,8), n_result
c           endif

           if (n.eq.1) then 
           write(*,fmt="('New b-values: ', i2, f10.5)") k, n_result
           endif
           
        enddo
      end do

**********now calculate rhox*****************
      write(*,*) 'now calculate rhox'
      allocate(rhox(ndepth_ref,nfile))
      do file=1,9
      call calcrhox(tauR(:,file),xkapref(:,file),ndepth_ref,
     &                                                 rhox(:,file))
      enddo

**********calculate estimated error********
      if (optimize) then
         write(*,*) 'now calculate error'
      call calc_error(xinf,xsup,yinf,ysup,zinf,zsup,temp_ref,
     & logg_ref,z_ref)
      endif

********write interpolated model in file nber out (basma compatible format)***********
********MB add print of interpolation coefficients
      
      write(*,*) 'now write result'
c      write(*,*) FILE_IN(out)
c      stop
       open(unit=23,file=FILE_IN(out))
       open(unit=25,file=FILE_IN(out+1))
       open(unit=27,file=FILE_IN(out+2))

       if (sph(1)) then
        write(*,*) 'spherical models'
        write(23,1967) ndepth_ref,xlr(out),logg_ref
 1967   format('''sphINTERPOL''',1x,i3,f8.0,2x,f4.2,1x,'0 0.00')
         do k=1,ndepth_ref
           write(23,1968) taus(k,out),T(k,out),Pe(k,out),
     &                    Pg(k,out),xit(k,out),rr(k,out),taur(k,out)
         enddo  
 1968   format(f8.4,1x,f8.2,3(1x,f8.4),1x,e15.6,1x,f8.4)
        write(25,19671) ndepth_ref,xlr(out)
19671   format('sphINTERPOL',1x,i3,f8.0)
        write(25,19672)
19672   format('    k    log(tau)  T    log(Pe)   log(Pg)    rhox')
         do k=1,ndepth_ref
           write(25,19681) k,taus(k,out),T(k,out),Pe(k,out),Pg(k,out),
     &     rhox(k,out)
19681      format(i5,1x,f8.4,1x,f8.2,2(1x,f8.4),1x,e15.6)
        enddo
        
        else
        write(*,*) 'plane parallel models'
        write(23,1966) ndepth_ref,xlr(out),logg_ref
1966     format('''ppINTERPOL''',1x,i3,f8.0,2x,f4.2,1x,'0 0.00')
         do k=1,ndepth_ref
           write(23,1965) taus(k,out),T(k,out),Pe(k,out),
     &                    Pg(k,out),xit(k,out),rr(k,out),taur(k,out)
1965       format(f8.4,1x,f8.2,3(1x,f8.4),1x,e15.6,1x,f8.4)
          enddo
        write(25,19661) xlr(out)
19661   format('ppINTERPOL',f8.0) 
        write(25,19672)
         do k=1,ndepth_ref
           write(25,19681) k,taus(k,out),T(k,out),Pe(k,out),Pg(k,out),
     &                   rhox(k,out)
         enddo  
      end if
      
******** MB write interpolation coefficients into a TEXT file      
       do k=1,8
         write(27,1969) coefval(k), power(k,:)
       enddo  
 1969  format('# ', a15,3(1x,f10.6))
       write(27,1971) z_ref+7.50
 1971  format(f6.2,1x)
       write(27,1972) ndepth_ref
 1972  format(i3,1x)
       write(27,1973) n_lev
 1973  format(i4,1x)
       do k=1,ndepth_ref
         write(27,1974) taus(k,out)
        enddo
 1974  format(f8.4,1x)
       do n=1,ndepth_ref
         write(27,1975)  (nnbvals(n,m,out), m=1, n_lev)
        enddo
 1975  format(1000(f10.5,2x))
         
       write(23,*) (FILE_IN(file),file=1,8)
       write(25,*) (FILE_IN(file),file=1,8)
       write(27,*) (FILE_IN(file),file=1,8)

       close(23)
       close(25)
       close(27)
c
c
******** MB write interpolation coefficients into a BINARY file  
c
c
c  WRITES DEPARTURE COEFFICENTS AND LOG(TAU500) SCALE INTO
c  UNFORMATTED BINARY FILE
c  THIS BINARY IS INTENDED TO BE USED AS AN INPUT FOR THE TS
c  CREATED BY ES MAY 2020
c
C
C some complicated stuff to have nice formatting of the output binary
C as 'sun_marcs_t5777_4.44_0.00_vmic1_new_A(CA)6.34.bin'
C
c      ABND = 7.0
c      ATOMID = 'Fe'
c      ATMOID = ''
c      
c      write(ABNDstring, '(F10.2)') ABND
c      ABNDstring = TRIM(ADJUSTL(ABNDstring))
c
c      if(len(ADJUSTL(TRIM(ATOMID))).EQ.1) then
c         write(BinaryName, '(A200, A4, A1, A1, A10)')
c     &        TRIM(ATMOID), '_A__',
c     &        ADJUSTL(TRIM(ATOMID)), '_', ABNDstring
c      else
c         write(BinaryName, '(A200, A3, A2, A1, A10)')
c     &        TRIM(ATMOID), '_A_',
c     &        TRIM(ATOMID(1:2)), '_', ABNDstring
c      endif
c      
c      BinaryName = ADJUSTL(TRIM(BinaryName))
c      
c      write(BinaryName, '(A220, A6)') TRIM(BinaryName), '.tsbin'
c      BinaryName = ADJUSTL(TRIM(BinaryName))
c
c      call MOPEN(LBIN4TS, BinaryName, 0,'NEW')
c
c      write(BinaryHeader, '(A1000)')  'This file was created by
c     *Ekaterina Semenova and contains departure coefficients from MULTI'
c
c      write(LBIN4TS) ATMOID
c      write(LBIN4TS) REAL(ABND)
c      write(LBIN4TS) BinaryHeader
c      write(LBIN4TS) ndepth_ref
c      write(LBIN4TS) n_lev
c      write(LBIN4TS) taus(1:ndepth_ref,out)
c      
c      DO I = 1, NK
c          WRITE(LBIN4TS) N(I,1:NDEP)/NSTAR(I,1:NDEP)
c      ENDDO
c      CALL CLOSE(LBIN4TS)
c      END
       
       

********** write a file compatible for sm (used for a control plot) ************
      open (unit=24,file='modele.sm')
      do file=1,8
       do k=1,ndepth_ref
        write(24,*) taus(k,file),T(k,file),Pe(k,file),
     &   Pg(k,file),taur(k,file)   
       end do
      end do
      
!case of a 10th comparison model
      read(*,*) FILE_IN(11)
      if (test) then 
         file = 11
      if (binary) then   
      call extract_bin(FILE_IN(file),teff(file),logg(file),metal(file),
     & ndepth(file),xlr(file),taus_aux(:,file),tauR_aux(:,file),
     & T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),xit_aux(:,file),
     & rr_aux(:,file),sph(file),xkapref_aux(:,file))
      else 
      call extract_ascii(FILE_IN(file),teff(file),logg(file),
     & metal(file),ndepth(file),xlr(file),taus_aux(:,file),
     & tauR_aux(:,file),T_aux(:,file),Pe_aux(:,file),Pg_aux(:,file),
     & xit_aux(:,file),rr_aux(:,file),sph(file),xkapref_aux(:,file))
      endif
      do k=1,ndepth(file)
       write(24,*) taus_aux(k,file),T_aux(k,file),Pe_aux(k,file),
     &               Pg_aux(k,file),tauR_aux(k,file)  
      end do
      end if
      close(24)      

      if (extrapol) then
          write (*,*) 'extrapolation done'
          else
      write (*,*) 'interpolation done'
      end if
      end if

      deallocate(taus,tauR,T,Pe,Pg,xit,taus_aux,tauR_aux,T_aux,Pe_aux,
     & Pg_aux,xit_aux,rr_aux,rr,rhox)

      end


c----------------------------------------------------------------------------------------



c--------------------------------------------------------------------------------
      real function inf(tab)
      implicit none
      integer :: n
      real,dimension(9) :: tab
      inf=tab(1)
      do n=2,8
         if (tab(n).lt.inf) then
            inf=tab(n)
         end if
      end do
      end function

      real function sup(tab)
      implicit none
      integer :: n
      real,dimension(9) :: tab
      sup=tab(1)
      do n=2,8
         if (tab(n).gt.sup) then
            sup=tab(n)
         end if
      end do
      end function



   


c---------------------------------------------------------------------------------
      subroutine blend_103 (r,s,t,x000,x001,x010,x011,x100,x101,x110, 
     & x111, x )
!
!*******************************************************************************
!from http://www.math.iastate.edu/burkardt/f_src/
!
!! BLEND_103 extends scalar point data into a cube.
!
!
!  Diagram:
!
!    011--------------111
!      |               |
!      |               |
!      |               |
!      |               |
!      |               |
!    001--------------101
!
!
!      *---------------*
!      |               |
!      |               |
!      |      rst      |
!      |               |
!      |               |
!      *---------------*
!
!
!    010--------------110
!      |               |
!      |               |
!      |               |
!      |               |
!      |               |
!    000--------------100
!
!
!  Formula:
!
!    Written as a polynomial in R, S and T, the interpolation map has the
!    form:
!
!      X(R,S,T) =
!        1         * ( + x000 )
!      + r         * ( - x000 + x100 )
!      +     s     * ( - x000        + x010 )
!      +         t * ( - x000               + x001 )
!      + r * s     * ( + x000 - x100 - x010                       + x110 )
!      + r     * t * ( + x000 - x100        - x001        + x101 )
!      +     s * t * ( + x000        - x010 - x001 + x011 )
!      + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, T, the coordinates where an interpolated value
!    is desired.
!
!    Input, real X000, X001, X010, X011, X100, X101, X110, X111, the
!    data values at the corners.
!
!    Output, real X, the interpolated data value at (R,S,T).
!
      implicit none
!
      real r
      real s
      real t
      real x
      real x000
      real x001
      real x010
      real x011
      real x100
      real x101
      real x110
      real x111
!
!  Interpolate the interior point.
!
      
      x = 
     & 1.0E+00     * ( + x000 ) 
     & + r         * ( - x000 + x100 ) 
     & +     s     * ( - x000        + x010 ) 
     & +         t * ( - x000               + x001 ) 
     & + r * s     * ( + x000 - x100 - x010                      + x110) 
     & + r     * t * ( + x000 - x100        - x001        + x101 ) 
     & +     s * t * ( + x000        - x010 - x001 + x011 ) 
     & + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110+
     &                                                            x111 )

      return
      end
c--------------------------------------------------------------------------------














c---------------------------------------------------------------------------------
      subroutine extract_bin(FILE,TEFF,grav,metal,ndepth,xlr_ref,tau5,
     &                tauR,temp,prese,presg,xit,rad,sph,kappa)

!extracted from osplot.f 07/2003, to get tau,T,Pe,Pg,microturb from a model
      implicit none
      integer :: ndp,ndepth,k,nlp,nlb
      parameter(ndp=200)
      CHARACTER*117 ADUM
      CHARACTER*256 FILE,file2
      CHARACTER*30 COMMENT
      real metal,grav,xlr_ref,TEFF,GG,radius
      logical :: sph
      real :: tau(ndp),t(ndp),z(ndp),ptot(ndp),prad(ndp),
     &  pg(ndp),pturb(ndp),pe(ndp),ro(ndp),xmass(ndp),xkapr(ndp)
      real :: gradp(ndp),gravity(ndp),pcheck(ndp),rr(ndp),
     &   xit(ndp),geff(ndp),gradptur(ndp),dp(ndp),taus(ndp),xlr(30),
     &   coldens(ndp)
      real :: tau5(ndp),tauR(ndp),temp(ndp),prese(ndp),
     &   presg(ndp),rad(ndp),kappa(ndp)
      real :: xlb(155000),w(155000),fluxme(155000)
      real :: presmo(30,ndp),ptio(ndp)
      real :: bPPR(NDP),bPPT(NDP),bPP(NDP),bGG(NDP),
     & bZZ(NDP),bDD(NDP),
     & bVV(NDP),bFFC(NDP),bPPE(NDP),bTT(NDP),
     & bTAULN(NDP),erad(ndp)
      integer :: NbTAU,IbTER
      common /struct/ tau,t,z,ptot,prad,pg,pturb,pe,ro,rr,taus,xlr,
     &                nlp,xkapr
      common /spectr/ nlb,xlb,w,fluxme
      common /pressure/ presmo,ptio
      common /binstruc/bPPR,bPPT,bPP,bGG,
     & bZZ,bDD,
     & bVV,bFFC,bPPE,bTT,
     & bTAULN,NbTAU,IbTER,erad
      common /radius/ radius
      OPEN(UNIT=10,FILE=FILE,STATUS='OLD',FORM='UNFORMATTED')
c     &     convert='big_endian')
ccc     &     RECL=152600)
*
      CALL READMO(10,NDEPTH,TEFF,GG,metal,sph)
c         open(21,file=file2,status='unknown')
      do k=1,ndepth
        tau5(k)=log10(taus(k))
        tauR(k)=log10(tau(k))
         temp(k)=t(k)
         prese(k)=log10(pe(k))
         presg(k)= log10(pg(k))
         kappa(k)=log10(xkapr(k))
      end do   
         xit=2.0
         xlr_ref=xlr(nlp)
         grav=log10(GG)
*        rad=rr
         if(.not.sph) radius=0.0
         rad=radius-z
c         write(21,1966) ndepth,xlr(nlp),log10(GG)
c1966     format('''INTERPOL''',1x,i3,f8.0,2x,f4.2,1x,'0 0.00')
c         do k=1,ndepth
c           write(21,1965) log10(taus(k)),t(k),log10(pe(k)),
c     &                    log10(pg(k)), xit
c1965       format(f8.4,1x,f8.2,3(x,f8.4))
c         enddo
c      close(21)
      close(10)
      END
C
      SUBROUTINE READMO(IARCH,JTAU,TEFF,G,metal,spherical)
C        THIS ROUTINE READS ONE MODEL, TO GET INFO ON PRAD
C             ( All features taken from listmo )
      PARAMETER (NDP=200)
C
      DIMENSION ABUND(16),TKORRM(NDP),FCORR(NDP),TAU(NDP),TAUS(NDP),
     *T(NDP),PE(NDP),PG(NDP),PRAD(NDP),PTURB(NDP),XKAPR(NDP),RO(NDP),
     *CP(NDP),CV(NDP),AGRAD(NDP),Q(NDP),U(NDP),V(NDP),ANCONV(NDP),
     *PRESMO(30,NDP),FCONV(NDP),RR(NDP),Z(NDP),EMU(NDP),HNIC(NDP)
     *,NJ(16),XLR(30),IEL(16),ANJON(16,5),PART(16,5),PROV(50,20+1),
     *ABSKA(50),SPRIDA(50),XLB(155000),FLUXME(155000),FLUMAG(155000),
     & PEP(16),
     * ABNAME(50),SOURCE(50),PTOT(NDP)
      DIMENSION W(155000),UW(12),BW(21),VW(25)
      CHARACTER*10 DAG,NAME,NAMEP,KLOCK
      CHARACTER*8 ABNAME,SOURCE
      DIMENSION WAVFLX(10)
      dimension PTIO(NDP)
      real*8 dluminosity
      real abSc,abTi,abV,abMn,abCo,metal
      logical spherical
      common /binstruc/ dummy(11*ndp+2),erad(ndp)
      common /struct/ tau,t,z,ptot,prad,pg,pturb,pe,ro,rr,taus,xlr,
     &                nlp,xkapr
      common /spectr/ nlb,xlb,w,fluxme
      common /pressure/ presmo,ptio
      common /radius/ radius
      DATA UW/0.145,0.436,0.910,1.385,1.843,2.126,2.305,2.241,1.270,
     *0.360,0.128,0.028/,BW/0.003,0.026,0.179,0.612,1.903,2.615,2.912,
     *3.005,2.990,2.876,2.681,2.388,2.058,1.725,1.416,1.135,0.840,0.568,
     *0.318,0.126,0.019/,VW/0.006,0.077,0.434,1.455,2.207,2.703,2.872,
     *2.738,2.505,2.219,1.890,1.567,1.233,0.918,0.680,0.474,0.312,0.200,
     *0.132,0.096,0.069,0.053,0.037,0.022,0.012/
      DATA NAME/'LOCAL'/,NAMEP/'PARSONS'/
      DATA A,B/.34785485,.65214515/
      IREAD=5
C
      REWIND IARCH
      READ(IARCH)
      erad=-1.e30
      READ(IARCH) INORD,DAG,KLOCK
      READ(IARCH) TEFF,FLUX,G,PALFA,PNY,PY,PBETA,ILINE,ISTRAL,
     &                MIHAL,IDRAB1,IDRAB2,IDRAB3,IDRAB4,IDRAB5,IDRAB6,
     &            ITMAX,NEL,(ABUND(I),I=1,NEL),abSc,abTi,abV,abMn,abCo
      GLOG=ALOG10(G)
      FNORD=0.1*INORD
C        CONVERT TO 'PHYSICAL FLUX'
      FLUX=3.14159*FLUX
      DO 2 I=1,NEL
    2 ABUND(I)=ALOG10(ABUND(I))+12.
      metal=abund(15)-7.50
      READ(IARCH)JTAU,NCORE,DIFLOG,TAUM,RADIUS,(RR(K),K=1,JTAU)
      if (jtau.gt.ndp) then
        print*, 'ERROR !!! Jtau (number of depths of model) = ',jtau
        print*, ' is larger than ndp!! Increase NDP.'
        stop
      endif
      if (radius.le.2.) then 
         spherical = .false.
      else
         spherical = .true.
         rr=radius-rr
      endif

      READ(IARCH)JTAU,(TKORRM(I),I=1,JTAU),(FCORR(K),K=1,JTAU)
      NTPO=0
      DO 3 K=1,JTAU
        READ(IARCH) KR,TAU(K),TAUS(K),Z(K),T(K),PE(K),PG(K),PRAD(K),
     &              PTURB(K),XKAPR(K),RO(K),EMU(K),CP(K),CV(K),
     &              AGRAD(K),Q(K),U(K),V(K),ANCONV(K),HNIC(K),NMOL,
     &              (PRESMO(J,K),J=1,NMOL),ptio(k)
        TAUK=ALOG10(TAU(K))+10.01
        KTAU=TAUK
        IF(ABS(TAUK-KTAU).GT.0.02) GO TO 31
        IF(KTAU.EQ.10) K0=K
        NTPO=NTPO+1
   31   CONTINUE
    3 CONTINUE
c      Z0=Z(K0)
c      DO 5 I=1,JTAU
c        Z(I)=Z(I)-Z0
c        i1=min(i+1,jtau)
c        PTOT(I)=PG(I)+PRAD(I)+0.5*(pturb(i)+pturb(i1))
    5 CONTINUE
***
      READ(IARCH)(NJ(I),I=1,NEL),NLP,(XLR(I),I=1,NLP)
     & ,NPROV,NPROVA,NPROVS,(ABNAME(KP),SOURCE(KP),KP=1,NPROV)
c      DO 22 KTAU=1,NTPO
c      DO 20 IE=1,NEL
c      NJP=NJ(IE)
c      READ(IARCH) KR,TAUI,TI,PEI,IEL(IE),ABUND(IE),
c     &            (ANJON(IE,JJ),JJ=1,NJP),(PART(IE,JJ),JJ=1,NJP)
c   20 CONTINUE
c      DO 21 KLAM=1,NLP
c      READ(IARCH) KR,TAUIL,(PROV(J,KLAM),J=1,NPROV),
c     &            ABSKA(KLAM),SPRIDA(KLAM)
   21 CONTINUE
   22 continue
c      READ(IARCH) NLB,(XLB(J),FLUXME(J),J=1,NLB),(W(J),J=1,NLB)
C CONVERT TO 'PHYSICAL' FLUXES
c      DO 24 J=1,NLB
c   24 FLUXME(J)=3.14159*FLUXME(J)

c      dluminosity=0.
c      do 25 j=1,nlb
c       dluminosity=dluminosity+fluxme(j)*w(j)
25    continue
c      dluminosity=dluminosity*4.*3.14159*radius**2/3.82d33
c      ddddd=real(dluminosity)
c      print*,'luminosity: ',dluminosity*3.82d33,' erg/s  = ',ddddd,
c     &  ' solar luminosities'
***
      RETURN
         END
*****************************************************************************

c---------------------------------------------------------------------------------
      subroutine extract_ascii(FILE,TEFF,grav,metal,ndepth,xlr_ref,tau5,
     &                tauR,temp,prese,presg,xit,rad,sph,xkapr)
c     adapted from P. DeLaverny   
      implicit none
      integer :: ndp,k
      parameter(ndp=200)
      integer :: imod,idum,ndepth
      CHARACTER*117 ADUM
      CHARACTER*256 FILE,file2
      CHARACTER*30 COMMENT
      CHARACTER*50 MOCODE,blabla      
      real metal,radius,mass,grav,xlr_ref,TEFF,GG,xic
      logical :: sph
      real :: tau(ndp),t(ndp),z(ndp),ptot(ndp),prad(ndp),
     &  pg(ndp),pturb(ndp),pe(ndp),ro(ndp),xmass(ndp),xkapr(ndp)
      real :: gradp(ndp),gravity(ndp),pcheck(ndp),rr(ndp),
     &  xit(ndp),geff(ndp),gradptur(ndp),dp(ndp),taus(ndp),xlr(30),
     &   coldens(ndp)
      real :: tau5(ndp),tauR(ndp),temp(ndp),prese(ndp),
     &   presg(ndp),rad(ndp),emu(ndp),vconv(ndp),fconv(ndp) 
      real :: dimension xlb(155000),w(155000),fluxme(155000)

      

          blabla=''
          imod =10
          OPEN(UNIT=imod,FILE=FILE,STATUS='OLD')
          read(imod,'(a)') mocode
c          print*,mocode,' = mocode'
c          if (mocode(1:1).eq.'p' .or. mocode(1:3).eq.'sun') then
c            print*,' this model is PLANE PARALLEL'
c          else if (mocode(1:1).eq.'s') then
c            print*,' this model is SPHERICALLY SYMMETRIC'
c          else
c            print*,' This model may not be a NewMARCS model!'
c          endif
          sph=(mocode(1:1).eq.'s')
          xlr_ref=5000
          read(imod,*)TEFF
          read(imod,*)
          read(imod,*)grav
          grav=log10(grav)
          read(imod,*)xic
          read(imod,*)mass
          read(imod,*)metal
          read(imod,*)radius
          do while (blabla.ne.'Model structure')
            read(imod,'(a)') blabla
          enddo
          backspace(imod)
          backspace(imod)
          read(imod,*)ndepth
          read(imod,*)
          read(imod,*)          
          do k=1,ndepth
            read(imod,*) idum,tauR(k),tau5(k),rad(k),temp(k),
     &                   Pe(k),Pg(k)
            prese(k)  = log10(Pe(k))
            presg(k)  = log10(Pg(k))
            xit(k) = xic
            rad(k)=radius-rad(k)
          enddo
          read(imod,*)
          do k=1,ndepth
            read(imod,*) idum,tauR(k),xkapr(k),ro(k),emu(k),
     &                   Vconv(k),Fconv(k)
            xkapr(k) = log10(xkapr(k))
          enddo
         close(imod)

      END
C


      subroutine resample(taus,tauR,T,Pe,Pg,xit,rr,xkapref)
      implicit  none 
      integer :: nlinemod,file,nfile,k,i
      real,dimension(:,:) :: taus,tauR,T,Pe,Pg,xit,rr,xkapref
      real,dimension(:,:),allocatable :: taubas
      real,dimension(:),allocatable :: tauresample,taustemp,tauRtemp,
     & Ttemp,Petemp,Pgtemp,xittemp,rrtemp,xkapreftemp

      INTERFACE 
      function SevalSingle(u,x,y)
      REAL,INTENT(IN) :: u  ! abscissa at which the spline is to be evaluated
      REAL,INTENT(IN),DIMENSION(:) :: x ! abscissas of knots
      REAL,INTENT(IN),DIMENSION(:):: y ! ordinates of knots
      real SevalSingle
      end
      END INTERFACE 

      
      nlinemod=size(taus,1)
      nfile=size(taus,2)-3

      allocate(tauresample(nlinemod),taustemp(nlinemod),
     &  tauRtemp(nlinemod),Ttemp(nlinemod)
     & ,Petemp(nlinemod),Pgtemp(nlinemod),xittemp(nlinemod),
     & rrtemp(nlinemod),xkapreftemp(nlinemod),
     &   taubas(nlinemod,size(taus,2)))
      
!!!!choose here the depth basis for interpolation (tau5000,tauRoss)!!!!!
      taubas=tauR
      write(*,*) 'resample models on common depth basis: tauRoss'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call common_depth_basis(tauresample,taubas,nlinemod,nfile)
c now do the resampling with the common tau
       do file=1,nfile
         do k=1,nlinemod

       taustemp(k)=SevalSingle(tauresample(k),taubas(:,file)
     &                                                  ,taus(:,file))
       
       tauRtemp(k)=SevalSingle(tauresample(k),taubas(:,file)
     &                                                  ,tauR(:,file))
       Ttemp(k)=SevalSingle(tauresample(k),taubas(:,file),T(:,file))
       Petemp(k)=SevalSingle(tauresample(k),taubas(:,file),Pe(:,file)) 
       Pgtemp(k)=SevalSingle(tauresample(k),taubas(:,file),Pg(:,file))
       xittemp(k)=SevalSingle(tauresample(k),taubas(:,file),xit(:,file))
       rrtemp(k)=SevalSingle(tauresample(k),taubas(:,file),rr(:,file))
       xkapreftemp(k)=
     &        SevalSingle(tauresample(k),taubas(:,file),xkapref(:,file))
         end do

         taus(:,file)=taustemp
         tauR(:,file)=tauRtemp
         T(:,file)=Ttemp 
         Pe(:,file)=Petemp
         Pg(:,file)=Pgtemp
         xit(:,file)=xittemp
         rr(:,file)=rrtemp
         xkapref(:,file)=xkapreftemp
       end do  

       deallocate(tauresample,taustemp,tauRtemp,Ttemp,Petemp
     &,Pgtemp,xittemp,rrtemp,xkapreftemp,taubas) 

       end       


!*******************************************************************************

      subroutine common_depth_basis(tauresample,tau,nlinemod,nfile)
      implicit none
      integer :: file,nlinemod,nfile
      real,dimension(nlinemod,nfile) :: tau
      real,dimension(nlinemod) :: tauresample

      tauresample=0
c initialize the common tau(5000) with min depth = max of the min depth of the models 
c                                  and max depth = min of the max depth of the models 
c essential for  the resampling with cubic spline
      tauresample(1)=tau(1,1)
      tauresample(nlinemod)=tau(nlinemod,1)
      do file=2,nfile
         if (tauresample(1).lt.tau(1,file)) then
            tauresample(1)=tau(1,file)
          end if
         if (tauresample(nlinemod).gt.tau(nlinemod,file)) then
            tauresample(nlinemod)=tau(nlinemod,file)
         end if
      end do   
      call blend_i_0d1 ( tauresample, nlinemod )
      end




!*******************************************************************************

      subroutine blend_i_0d1 ( x, m )
!
!
!! BLEND_I_0D1 extends indexed scalar data at endpoints along a line.
!http://orion.math.iastate.edu/burkardt/f_src/f_src.html
!
!  Diagram:
!
!    ( X1, ..., ..., ..., ..., ..., XM )
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    15 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X(M).  
!
!    On input, X(1) and X(M) contain scalar values which are to be 
!    interpolated through the entries X(2) through X(M).  It is assumed
!    that the dependence of the data is linear in the vector index I.  
!    
!    On output, X(2) through X(M-1) have been assigned interpolated 
!    values.
!
!    Input, integer M, the number of entries in X.
                               
      implicit none
!
      integer m
!
      integer i
      real r
      real x(m)
!
      do i = 2, m - 1

        r = real ( i - 1 ) / real ( m - 1 )

        call blend_101 ( r, x(1), x(m), x(i) )

      end do

      return
      end

      subroutine blend_101 ( r, x0, x1, x )
!
!*******************************************************************************
!
!! BLEND_101 extends scalar endpoint data to a line.
!http://orion.math.iastate.edu/burkardt/f_src/f_src.html
!
!  Diagram:
!
!    0-----r-----1
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the coordinate where an interpolated value is desired.  
!
!    Input, real X0, X1, the data values at the ends of the line.
!
!    Output, real X, the interpolated data value at (R).
!
      implicit none
!
      real r
      real x
      real x0
      real x1
!
      x = ( 1.0E+00 - r ) * x0 + r * x1

      return
      end
!----------------------------------------------------------------------------

      REAL FUNCTION SevalSingle(u,x,y) 
! ---------------------------------------------------------------------------
!http://www.pdas.com/fmm.htm
!  PURPOSE - Evaluate the cubic spline function
!     Seval=y(i)+b(i)!(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
!           where  x(i) <= u < x(i+1)

!  NOTES- if u<x(1), i=1 is used;if u>x(n), i=n is used

      REAL,INTENT(IN) :: u  ! abscissa at which the spline is to be evaluated
      REAL,INTENT(IN),DIMENSION(:) :: x ! abscissas of knots
      REAL,INTENT(IN),DIMENSION(:):: y ! ordinates of knots
      REAL, DIMENSION(:),allocatable :: b,c,d ! linear,quadratic,cubic coeff

      INTEGER, SAVE :: i=1
      INTEGER :: j, k, n
      REAL:: dx

      INTERFACE 
         subroutine FMMsplineSingle(x, y, b, c, d)
        REAL,DIMENSION(:), INTENT(IN)  :: x ! abscissas of knots
        REAL,DIMENSION(:), INTENT(IN)  :: y ! ordinates of knots
        REAL,DIMENSION(:), INTENT(OUT) :: b ! linear coeff
        REAL,DIMENSION(:), INTENT(OUT) :: c ! quadratic coeff.
        REAL,DIMENSION(:), INTENT(OUT) :: d ! cubic coeff.
         end
      END INTERFACE
!----------------------------------------------------------------------------
      n=SIZE(x)
      allocate(b(n),c(n),d(n))
      call FMMsplineSingle(x, y, b, c, d)
!.....First check if u is in the same interval found on the
!        last call to Seval.............................................
       IF (  (i<1) .OR. (i >= n) ) i=1
       IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
         i=1   ! binary search
         j=n+1
     
         DO
           k=(i+j)/2
           IF (u < x(k)) THEN
             j=k
           ELSE
        i=k
           END IF
           IF (j <= i+1) EXIT
         END DO
       END IF
     
        dx=u-x(i)   ! evaluate the spline
        SevalSingle=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
     
        RETURN
        deallocate(b,c,d)
      END Function SevalSingle  ! -------------------------------------------------------


      SUBROUTINE FMMsplineSingle(x, y, b, c, d)
! ---------------------------------------------------------------------------
!http://www.pdas.com/fmm.htm
! PURPOSE - Compute the coefficients b,c,d for a cubic interpolating spline
!  so that the interpolated value is given by
!    s(x) = y(k) + b(k)*(x-x(k)) + c(k)*(x-x(k))**2 + d(k)*(x-x(k))**3
!      when x(k) <= x <= x(k+1)
!  The end conditions match the third derivatives of the interpolated curve to
!  the third derivatives of the unique polynomials thru the first four and
!  last four points.
!  Use Seval or Seval3 to evaluate the spline.
        REAL,DIMENSION(:), INTENT(IN)  :: x ! abscissas of knots
        REAL,DIMENSION(:), INTENT(IN)  :: y ! ordinates of knots
        REAL,DIMENSION(:), INTENT(OUT) :: b ! linear coeff
        REAL,DIMENSION(:), INTENT(OUT) :: c ! quadratic coeff.
        REAL,DIMENSION(:), INTENT(OUT) :: d ! cubic coeff.
     
        INTEGER:: k,n
        REAL:: t,aux
        REAL,PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
!----------------------------------------------------------------------------
       n=SIZE(x)

       IF (n < 3) THEN   ! Straight line - special case for n < 3
         b(1)=ZERO
         IF (n == 2) b(1)=(y(2)-y(1))/(x(2)-x(1))
         c(1)=ZERO
         d(1)=ZERO
         IF (n < 2) RETURN
         b(2)=b(1)
         c(2)=ZERO
         d(2)=ZERO
         RETURN
       END IF
  
!.....Set up tridiagonal system.........................................
!.    b=diagonal, d=offdiagonal, c=right-hand side
        d(1)=x(2)-x(1)
        c(2)=(y(2)-y(1))/d(1)
       DO k=2,n-1
         d(k)=x(k+1)-x(k)
         b(k)=TWO*(d(k-1)+d(k))
         c(k+1)=(y(k+1)-y(k))/d(k)
         c(k)=c(k+1)-c(k)
       END DO
    
!.....End conditions.  third derivatives at x(1) and x(n) obtained
!.       from divided differences.......................................
       b(1)=-d(1)
       b(n)=-d(n-1)
       c(1)=ZERO
       c(n)=ZERO
       IF (n > 3) THEN
         c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
         c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
         c(1)=c(1)*d(1)*d(1)/(x(4)-x(1))
         c(n)=-c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
       END IF
     
       DO k=2,n    ! forward elimination
         t=d(k-1)/b(k-1)
         b(k)=b(k)-t*d(k-1)
         c(k)=c(k)-t*c(k-1)
       END DO
     
       c(n)=c(n)/b(n)   ! back substitution ( makes c the sigma of text)
       DO k=n-1,1,-1
         c(k)=(c(k)-d(k)*c(k+1))/b(k)
       END DO
     
!.....Compute polynomial coefficients...................................
       b(n)=(y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
       DO k=1,n-1
         b(k)=(y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
         d(k)=(c(k+1)-c(k))/d(k)
         c(k)=THREE*c(k)
       END DO
       c(n)=THREE*c(n)
       d(n)=d(n-1)
     
       RETURN
       END Subroutine FMMsplineSingle ! ---------------------------------------------------

      SUBROUTINE NaturalSplineSingle(x,y,b,c,d)
! ---------------------------------------------------------------------------
! PURPOSE - Construct the natural spline thru a set of points
! NOTES - A natural spline has zero second derivative at both endpoints.
     
       REAL,INTENT(IN),DIMENSION(:):: x,y   ! coordinates of knots
       REAL,INTENT(OUT),DIMENSION(:):: b,c,d  ! cubic coeff.
     
       INTEGER:: k,n
       REAL,PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
!-----------------------------------------------------------------------
       n=SIZE(x)
      
       IF (n < 3) THEN   ! Straight line - special case for n < 3
         b(1)=ZERO
         IF (n == 2) b(1)=(y(2)-y(1))/(x(2)-x(1))
         c(1)=ZERO
         d(1)=ZERO
         b(2)=b(1)
         c(2)=ZERO
         d(2)=ZERO
         RETURN
       END IF

       d(1:n-1) = x(2:n)-x(1:n-1)  ! Put the h-array of the text into array d

!.....Set up the upper triangular system in locations 2 thru n-1 of
!        arrays b and c. B holds the diagonal and c the right hand side.
       b(2)=TWO*(d(1)+d(2))
       c(2)=(y(3)-y(2))/d(2)-(y(2)-y(1))/d(1)
       DO  k=3,n-1
         b(k)=TWO*(d(k-1)+d(k))-d(k-1)*d(k-1)/b(k-1)
       c(k)=(y(k+1)-y(k))/d(k)-(y(k)-y(k-1))/d(k-1)-d(k-1)*c(k-1)/b(k-1)
       END DO
     
       c(n-1)=c(n-1)/b(n-1)   ! Back substitute to get c-array
       DO  k=n-2,2,-1
         c(k)=(c(k)-d(k)*c(k+1))/b(k)
       END DO
       c(1)=ZERO
       c(n)=ZERO   ! c now holds the sigma array of the text 
     
     
!.....Compute polynomial coefficients ..................................
       b(n)=(y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
       DO  k=1,n-1
         b(k)=(y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
         d(k)=(c(k+1)-c(k))/d(k)
         c(k)=THREE*c(k)
       END DO
       c(n)=THREE*c(n)
       d(n)=d(n-1)
       RETURN
  
       END Subroutine NaturalSplineSingle 



!----------------------------------------------------------------
      subroutine calcrhox(tau,kappa,ndepth,rhox)
c     2 ways to calculate rhox : int(ro*dx) or int(1/kappa*dtau)
c      A&A 387, 595-604 (2002)      
      implicit none
      integer :: i,ndepth
      real :: first
      real, dimension(ndepth) :: tau,kappa,rhox,f,x
      real :: tot
      real, external :: rinteg


      f=(1/10**kappa)
      x=10**(tau)
      first = x(1)*f(1)
      tot=rinteg(x,f,rhox,ndepth,first) 
      do i=2,ndepth
           rhox(i) = rhox(i-1) + rhox(i)
      enddo
      end

      real function rinteg(x,f,fint,n,start)
c******************************************************************************
c     This routine is from ATLAS6
c******************************************************************************
      implicit none
      integer :: n,i,n1
      real x(5000), f(5000), fint(5000)
      real a(5000), b(5000), c(5000)
      real :: start

      call parcoe (f,x,a,b,c,n)
      fint(1) = start 
      rinteg = start
      n1 = n - 1
      do 10 i=1,n1
         fint(i+1)= (a(i)+b(i)/2.*(x(i+1)+x(i))+ 
     .     c(i)/3.*((x(i+1)+x(i))*x(i+1)+x(i)*x(i)))*(x(i+1)-x(i))
10    rinteg = rinteg + fint(i+1)

      return
      end 




      subroutine parcoe(f,x,a,b,c,n)

      implicit none
      integer :: n,n1,j,j1
      real f(5000), x(5000), a(5000), b(5000), c(5000)
      real :: d,wt

      c(1)=0.
      b(1)=(f(2)-f(1))/(x(2)-x(1))
      a(1)=f(1)-x(1)*b(1)
      n1=n-1
      c(n)=0.
      b(n)=(f(n)-f(n1))/(x(n)-x(n1))
      a(n)=f(n)-x(n)*b(n) 
      if(n.eq.2)return
      do 1 j=2,n1
      j1=j-1
      d=(f(j)-f(j1))/(x(j)-x(j1)) 
      c(j)=f(j+1)/((x(j+1)-x(j))*(x(j+1)-x(j1)))-f(j)/((x(j)-x(j1))*
     1(x(j+1)-x(j)))+f(j1)/((x(j)-x(j1))*(x(j+1)-x(j1)))
      b(j)=d-(x(j)+x(j1))*c(j)
    1 a(j)=f(j1)-x(j1)*d+x(j)*x(j1)*c(j)
      c(2)=0.
      b(2)=(f(3)-f(2))/(x(3)-x(2))
      a(2)=f(2)-x(2)*b(2)
      c(3)=0.
      b(3)=(f(4)-f(3))/(x(4)-x(3))
      a(3)=f(3)-x(3)*b(3)
      if(c(j).eq.0.)go to 2
      j1=j+1
      wt=abs(c(j1))/(abs(c(j1))+abs(c(j)))
      a(j)=a(j1)+wt*(a(j)-a(j1))
      b(j)=b(j1)+wt*(b(j)-b(j1))
      c(j)=c(j1)+wt*(c(j)-c(j1))
    2 continue
      a(n1)=a(n)
      b(n1)=b(n)
      c(n1)=c(n)
      return
      end 

***********************************************************************************
      subroutine calc_error(xinf,xsup,yinf,ysup,zinf,zsup,teff_ref,
     &   logg_ref,z_ref)
      implicit none
      real :: xinf,xsup,yinf,ysup,zinf,zsup,teff_ref,logg_ref,z_ref
      real ::  error_T,error_Pe,error_Pg,error_kappa
      real :: errorTeffT,errorloggT,errorzT,
     &        errorTeffPe,errorloggPe,errorzPe,
     &         errorTeffPg,errorloggPg,errorzPg,
     &         errorTeffkappa,errorloggkappa,errorzkappa
 
! values read out of the figures of the manual and scaled down o the according step
              errorTeffT=0.055/32
              errorloggT=0.008/5
              errorzT=0.015/4
              errorTeffPe=0.65/32
              errorloggPe=0.4/5
              errorzPe=0.38/4
              errorTeffPg=0.25/32
              errorloggPg=0.23/5
              errorzPg=0.35/4
              errorTeffkappa=0.8/32
              errorloggkappa=0.36/5
              errorzkappa=0.38/4
   

   
      error_T=min((xsup-teff_ref),(teff_ref-xinf))/100*errorTeffT +
     &         min((ysup-logg_ref),(logg_ref-yinf))*errorloggT +
     &         min((zsup-z_ref),(z_ref-zinf))*errorzT 
      write(*,1409) 'estimated max error on T =',error_T*100,'%'

                  
      error_Pe=min((xsup-teff_ref),(teff_ref-xinf))/100*errorTeffPe +
     &         min((ysup-logg_ref),(logg_ref-yinf))*errorloggPe +
     &         min((zsup-z_ref),(z_ref-zinf))*errorzPe 
       write(*,1409) 'estimated max error on Pe =',error_Pe*100,'%'
             

      error_Pg=min((xsup-teff_ref),(teff_ref-xinf))/100*errorTeffPg +
     &         min((ysup-logg_ref),(logg_ref-yinf))*errorloggPg +
     &         min((zsup-z_ref),(z_ref-zinf))*errorzPg 
       write(*,1409) 'estimated max error on Pg =',error_Pg*100,'%'
              

      error_kappa=min((xsup-teff_ref),(teff_ref-xinf))/100
     &                                                *errorTeffkappa +
     &         min((ysup-logg_ref),(logg_ref-yinf))*errorloggkappa +
     &         min((zsup-z_ref),(z_ref-zinf))*errorzkappa 
      write(*,1409) 'estimated max error on kappa =',error_kappa*100,'%'
            
 1409 format(a30,f5.1,a2)
      end
***********************************************************************************
      subroutine str2int(str,int,stat)
      implicit none
      ! Arguments
      character(len=*),intent(in)   :: str
      integer*8,intent(out)         :: int
      integer*8,intent(out)         :: stat
  
      read(str,*,iostat=stat)  int
      end 
