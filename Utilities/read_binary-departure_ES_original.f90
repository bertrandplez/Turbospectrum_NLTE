program readTSbin
implicit none
character*1000 header
character*500 atmosID
integer(4) :: ndep, nk, tmp(1), i,j
integer(8) :: pos1
real(8), allocatable :: tau(:)
real(8), allocatable :: b(:,:)
integer(8) :: positions(8)
open(unit=1, file='1D_NLTE_grid.bin', form='unformatted',  ACCESS="STREAM")
!read(1), header
!print*, header
!rewind(1)
    pos1 = 1! position from the auxiliary file
    print*, 'POS=', pos1
    ! atmosID
    read(1, pos=pos1) atmosID
    print*,'ATMOS:', atmosID, len(atmosID), sizeof(atmosID)
    ! NDEP
    pos1 = pos1 + 500 !e500 because atmosID is character*500
    ! I have to dump ndep as an 1-element array with Python, so do the following:
    read(1, pos=pos1) tmp
    ndep = tmp(1)
    print*, 'NDEP=', ndep
    ! NK
    pos1 = pos1 + 8
    ! I have to dump nk as an 1-element array with Python, so do the following:
    read(1, pos=pos1), tmp
    nk = tmp(1)
    print*, 'NK=', nk
    ! DEPTH SCALE, TAU500
    pos1 = pos1 + 8
    allocate(tau(ndep))
    read(1, pos=pos1), tau
    print*, 'TAU:', tau
    ! departure coefficients, an array of a NDEP x NK shape
    pos1 = pos1 +  ndep * 8
    allocate(b(ndep, nk))
    read(1, pos=pos1), b
    deallocate(tau)
    deallocate(b)
close(1)
endprogram
