program readTSbin
implicit none
character*1000 header
character*500 atmosID
integer(4) :: ndep, nk, tmp(1), i,j
integer(8) :: pos1
real(4) :: abu
real(8), allocatable :: tau(:)
real(8), allocatable :: b(:,:)
integer(8) :: positions(8)
open(unit=1, file='1D_NLTE_grid_h20_1_5.bin', form='unformatted',  ACCESS="STREAM")
open(unit=10, file='output.bin',  form='unformatted')
!read(1), header
!print*, header
!rewind(1)

    pos1 = 39697   ! position from the auxilliary file
!    pos1=1
    print*, 'POS=', pos1
    abu = 12.00 ! abundance from the auxilliary file

    ! atmosID
    read(1, pos=pos1) atmosID
    print*,'ATMOS:', trim(adjustl(atmosID))
    write(10) atmosID
    write(10) abu
    write(10) header  ! should be something else...

    ! NDEP
    pos1 = pos1 + 500 !e500 because atmosID is character*500
    ! I have to dump ndep as a 1-element array with Python, so do the following:
!    pos1=pos1+5 ! fix problem for Plez !!??!!
    pos1=pos1+5*5
!    do i=1,40
!    read(1,pos=pos1) ndep
!    print*,pos1,ndep
!    pos1=pos1+1
!    enddo
!    read(1, pos=pos1) tmp(1)
    read(1,pos=pos1),ndep
!    print*, tmp
!    ndep = tmp(1)
    print*, 'NDEP=', ndep
    write(10)  ndep

    ! NK
    pos1 = pos1 + 8
    ! I have to dump nk as a 1-element array with Python, so do the following:
    read(1, pos=pos1), tmp(1)
!    read(1,pos=pos1) nk
    print*,tmp
    nk = tmp(1)
    print*, 'NK=', nk
    write(10) nk
 
    ! DEPTH SCALE, TAU500
    pos1 = pos1 + 8
    allocate(tau(ndep))
    read(1, pos=pos1), tau
    print*, 'TAU:', tau
    do i=1,ndep
      write(10) sngl(tau(i))
    enddo

    ! departure coefficients, an array of a NDEP x NK shape
    pos1 = pos1 +  ndep * 8
    allocate(b(ndep, nk))
    read(1, pos=pos1), b
    do j=1,nk
      do i=1,ndep
        write(10) sngl(b(i,j))
      enddo
    enddo

    deallocate(tau)
    deallocate(b)
close(1)
endprogram
