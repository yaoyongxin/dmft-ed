MODULE ED_HAMILTONIAN_SHARED
  USE SF_MISC,    only: assert_shape
  USE SF_CONSTANTS,only:zero
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_SETUP
  implicit none

  !   !> MPI local variables (shared)
  ! #ifdef _MPI
  !   integer          :: MpiComm=MPI_UNDEFINED
  ! #else
  !   integer          :: MpiComm=0
  ! #endif
  !   logical          :: MpiStatus=.false.
  !   integer          :: MpiIerr
  !   integer          :: MpiRank=0
  !   integer          :: MpiSize=1
  !   logical          :: MpiMaster=.true.
  !   integer          :: MpiQ=1
  !   integer          :: MpiR=0

  !   integer          :: MpiIstart
  !   integer          :: MpiIend
  !   integer          :: MpiIshift
  !
  ! integer          :: Dim
  ! integer          :: DimUp
  ! integer          :: DimDw
  !
  integer          :: Hsector=0
  logical          :: Hstatus=.false.
  type(sector_map) :: H
  type(sector_map8) :: H8





end MODULE ED_HAMILTONIAN_SHARED
