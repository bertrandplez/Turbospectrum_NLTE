MODULE opac_module
!
! created by BPz 2026-08-13
! replaces common block large by allocatable array
!
!      doubleprecision xlambda
!      doubleprecision source_function
!      common/large/ xlambda(lpoint),source_function(ndp,lpoint),
!     & maxlam,ABSO(NDP,lpoint),
!     & absos(ndp,lpoint),absocont(ndp,lpoint),absoscont(ndp,lpoint)

  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: sp = selected_real_kind(6, 37)
  
  integer :: maxlam
  real(dp), allocatable :: source_function(:,:),xlambda(:)
  real(sp), allocatable :: abso(:,:),absos(:,:)
  real(sp), allocatable :: absocont(:,:),absoscont(:,:)

contains

  subroutine allocate_lambdaarray(nlam)
    implicit none
    integer :: nlam

    ! Allocate xlambda array if not already allocated
    if (.not. allocated(xlambda)) then
       allocate(xlambda(nlam))
       maxlam = nlam
    end if

  end subroutine allocate_lambdaarray


  subroutine allocate_opacarrays(ntau,nlam)
    implicit none
    integer, intent(in) :: ntau,nlam

    ! Allocate source_function array if not already allocated
    if (.not. allocated(source_function)) then
       allocate(source_function(ntau,nlam))
    end if

    ! Allocate abso array if not already allocated
    if (.not. allocated(abso)) then
       allocate(abso(ntau,nlam))
    end if

    ! Allocate absos array if not already allocated
    if (.not. allocated(absos)) then
       allocate(absos(ntau,nlam))
    end if

    ! Allocate absocont array if not already allocated
    if (.not. allocated(absocont)) then
       allocate(absocont(ntau,nlam))
    end if

    ! Allocate absoscont array if not already allocated
    if (.not. allocated(absoscont)) then
       allocate(absoscont(ntau,nlam))
    end if

  end subroutine allocate_opacarrays

END MODULE opac_module
