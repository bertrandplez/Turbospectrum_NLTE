!
      subroutine read_departure(iunit,departurefile,departbin,maxlevel,
     &                   modnlevel,ndepth,ndepth_read,taumod,
     &              b_departure,abundance_nlte,header_dep1,header_dep2)
!
! read NLTE departure coefficients for one model atom, and one atmosphere model
! created by BPlez 17/04-2020
!
! unit number to open, file name, maximum number of levels, number of levels to read, 
! maximum number of depth points, number of depth points found in file,
! optical depth, departure coefficients
!
      implicit none
      character departurefile*256,oneline*256
      character header_dep1*500,header_dep2*1000
      integer iunit,modnlevel,i,ndepth_read,j,ndepth,maxlevel
      integer modnlevel_read,k
      real taumod(ndepth),b_departure(ndepth,0:maxlevel)
      real abundance_nlte
      logical departbin
      character*15, dimension (8) :: coefval
      real, dimension(8,3) :: power

      if (departbin) then
        open(iunit,file=departurefile,form='unformatted',status='old',
     &     convert='little_endian')
        read(iunit) header_dep1
        read(iunit) abundance_nlte
        read(iunit) header_dep2
        read(iunit) ndepth_read
        print*,taumod
        print*,'read_departure, ndepth ',ndepth_read
        read(iunit) modnlevel_read
        print*,'read_departure, nlevel ',modnlevel_read
        if (ndepth.lt.ndepth_read) then
          print*,'ndepth in departure file ',
     &     ndepth_read,' is too large'
          print*,'increase dimension!'
          stop 'read_departure.f'
        endif
        if (modnlevel.ne.modnlevel_read) then
          print*,'nlevel in atom is',modnlevel,'in departure file ',
     &     modnlevel_read
          print*,'check consistency!'
          stop 'read_departure.f'
        endif
        do i=1,ndepth_read
          read(iunit) taumod(i)
        ! set departure coefficient to 1 for unidentified levels
          b_departure(i,0)=1.0
        enddo
        print*,taumod
        do j=1,modnlevel
          do i=1,ndepth_read
            read(iunit) b_departure(i,j)
          enddo
        enddo
        print*,'read departure coefficients for ',ndepth_read,
     &       ' depths in read_departure'
        close(iunit)
      else
        open(iunit,file=departurefile,status='old')
        do k=1,8
         read(iunit,1969) coefval(k), power(k,:)
        enddo  
 1969  format('# ', a15,3(1x,f10.6))
        read(iunit,1971) abundance_nlte
 1971  format(f6.3,1x)
        read(iunit,1972) ndepth_read
 1972  format(i3,1x)
        print*,'read_departure, ndepth ',ndepth_read
        read(iunit,1973) modnlevel_read
 1973  format(i4,1x)
        print*,'read_departure, nlevel ',modnlevel_read
        if (ndepth.lt.ndepth_read) then
          print*,'ndepth in departure file ',
     &     ndepth_read,' is too large'
          print*,'increase dimension!'
          stop 'read_departure.f'
        endif
        if (modnlevel.ne.modnlevel_read) then
           print*,'nlevel in atom is',modnlevel,'in departure file ',
     &      modnlevel_read
           print*,'check consistency!'
           stop 'read_departure.f'
        endif
        do i=1,ndepth_read
          read(iunit,1974) taumod(i)
          taumod(i)=10.**taumod(i)
! set departure coefficient to 1 for unidentified levels
          b_departure(i,0)=1.0
        enddo
 1974  format(f8.4,1x)
        do i=1,ndepth_read
         read(iunit,1975)  (b_departure(i,j), j=1, modnlevel)
        enddo
 1975  format(1000(f10.5,2x))
        print*,'read departure coefficients for ',ndepth_read,
     &       ' depths in read_departure'

        close(iunit)
      endif

      return
      end
