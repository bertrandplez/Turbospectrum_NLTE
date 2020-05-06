!
      subroutine read_departure(iunit,departurefile,maxlevel,modnlevel,
     &                   ndepth,ndepth_read,taumod,b_departure,
     &                   abundance_nlte,header_dep1,header_dep2)
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
      integer modnlevel_read
      real taumod(ndepth),b_departure(ndepth,maxlevel)
      real abundance_nlte

      open(iunit,file=departurefile,form='unformatted',status='old',
     &     convert='little_endian')
      read(iunit) header_dep1
      read(iunit) abundance_nlte
      read(iunit) header_dep2
      read(iunit) ndepth_read
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
      enddo

      do j=1,modnlevel
        do i=1,ndepth_read
          read(iunit) b_departure(i,j)
        enddo
      enddo

      print*,'read departure coefficients for ',ndepth_read,
     &       ' depths in read_departure'

      close(iunit)

      return
      end
