      program read_bin_nlte
! BPz 17/02-2021
! adapted fomr interpolator to read NLTE departure files
! adopted free format for data file. Model name needs to be set between quotes
!
      implicit none
      character fileaux*256, nlte_binary*256, id_model*256, fileout*256
      integer icomp,n_dep,n_lev,k,tmp(1),m,n
      integer*8 pos1
      real n_teff,n_logg,n_metal,n_abu,zref
      real*8, allocatable :: nlte_tau(:), nlte_data(:,:)

      real, dimension(8,3) :: power
      character*15, dimension (8) :: coefval

      DATA coefval/'tau500', 'tauross', 'T', 'logpe',
     &             'logpg', 'xit', 'rr', 'logkapref'/
 

      icomp=8
      print*, 'aux data txt file:'
      read(5,10) fileaux
10    format(a)
      print*, 'nlte data file:'
      read(5,10) nlte_binary
      print*, 'output file for departure coefficients'
      read(5,10) fileout

      open(199,file=fileaux,form='formatted')
      read(199,*)
      read(199, *)
     &        id_model, n_teff, n_logg, n_metal,
     &        n_abu, pos1

      close(199)
      open(unit=200,file=nlte_binary,form='unformatted',access="stream")


        read(200, pos = pos1) id_model
        write(*,fmt='("NLTE: ", i1,1x,A50,1x, i10)')
     &        1, trim(adjustl(id_model)), pos1

        pos1 = pos1 + 500
        read(200, pos = pos1) tmp
        n_dep = tmp(1)
        print*,'ndep',n_dep
        pos1 = pos1 + icomp
        read(200, pos = pos1) tmp
        n_lev = tmp(1)
        print*,'nlev',n_lev
        pos1 = pos1 + icomp
        allocate(nlte_tau(n_dep))
        read(200, pos = pos1), nlte_tau

        pos1 = pos1 + n_dep*icomp
        allocate(nlte_data(n_dep, n_lev))
        read(200, pos = pos1), nlte_data
        rewind(200)
      close(200)

******** MB write interpolation coefficients into a TEXT file
! dummy values

       power=1.0
       zref=0.0
!
       open(unit=27,file=fileout)
       do k=1,8
         write(27,1969) coefval(k), power(k,:)
       enddo 
 1969  format('# ', a15,3(1x,f10.6))
       write(27,1971) zref+7.50
 1971  format(f6.2,1x)
       write(27,1972) n_dep
 1972  format(i3,1x)
       write(27,1973) n_lev
 1973  format(i4,1x)
       do k=1,n_dep
         write(27,1974) log10(nlte_tau(k))
        enddo
 1974  format(1pe12.5,1x)
       do n=1,n_dep
         write(27,1975)  (nlte_data(n,m), m=1, n_lev)
        enddo
 1975  format(1000(1pe12.5,1x))

       write(27,*)

       close(27)
       stop
       end
