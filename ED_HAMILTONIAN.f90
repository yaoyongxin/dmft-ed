MODULE ED_HAMILTONIAN
  USE ED_HAMILTONIAN_SHARED
  USE ED_HAMILTONIAN_STORED_HxV
  USE ED_HAMILTONIAN_DIRECT_HxV
  !
  implicit none
  private

  !>Build sparse hamiltonian of the sector
  public  :: build_Hv_sector
  public  :: delete_Hv_sector
  public  :: vecDim_Hv_sector


  !>Sparse Mat-Vec product using stored sparse matrix
  public  :: spMatVec_cc
#ifdef _MPI
  public  :: spMatVec_MPI_cc
#endif


  !>Sparse Mat-Vec direct on-the-fly product
  public  :: directMatVec_cc
#ifdef _MPI
  public  :: directMatVec_MPI_cc
#endif





contains






  !####################################################################
  !                 MAIN ROUTINES: BUILD/DELETE SECTOR
  !####################################################################
  subroutine build_Hv_sector(isector,Hmat)
    integer                            :: isector,SectorDim
    complex(8),dimension(:,:),optional :: Hmat
    integer                            :: irank
    integer                            :: i,j,Dim
    !
    Hsector=isector
    Hstatus=.true.
    !
    call build_sector(isector,H)
    !
    Dim = getDim(isector)
    !
    !Total Split:
    MpiQ = Dim/MpiSize
    MpiR = 0
    if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)
    !
    MpiIshift = MpiRank*mpiQ! + mpiR
    MpiIstart = MpiRank*mpiQ + 1
    MpiIend   = (MpiRank+1)*mpiQ + mpiR
    !
#ifdef _MPI
    if(MpiStatus.AND.ed_verbose>4)then
       write(LOGfile,*)&
            "         mpiRank,   mpi_Q,   mpi_R,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
       do irank=0,MpiSize-1
          call Barrier_MPI(MpiComm)
          if(MpiRank==irank)then
             write(LOGfile,*)MpiRank,MpiQ,MpiR,MpiIstart,MpiIend,MpiIend-MpiIstart+1
          endif
       enddo
       call Barrier_MPI(MpiComm)
    endif
#endif
    !
    !
    if(present(Hmat))then
       spHtimesV_cc => null()
       call ed_buildh_c(isector,Hmat)
       return
    endif
    !
    select case (ed_sparse_H)
       !
       !
    case (.true.)
       spHtimesV_cc => spMatVec_cc
#ifdef _MPI
       if(MpiStatus)spHtimesV_cc => spMatVec_MPI_cc
#endif
       call ed_buildh_c(isector)
       !
       !
    case (.false.)
       spHtimesV_cc => directMatVec_cc
#ifdef _MPI
       if(MpiStatus)spHtimesV_cc => directMatVec_MPI_cc
#endif
    end select
    !
  end subroutine build_Hv_sector


  subroutine delete_Hv_sector()
    integer :: iud
    call delete_sector(Hsector,H)
    Hstatus=.false.
    !
#ifdef _MPI
    if(MpiStatus)then
       call sp_delete_matrix(MpiComm,spH0)
    else
       call sp_delete_matrix(spH0)
    endif
#else
    call sp_delete_matrix(spH0)
#endif
    !
    spHtimesV_cc => null()
    !
  end subroutine delete_Hv_sector


  function vecDim_Hv_sector(isector) result(vecDim)
    integer :: isector
    integer :: Dim
    integer :: vecDim
    !
    Dim  = getdim(isector)
    !
#ifdef _MPI
    if(MpiStatus)then
       MpiQ = Dim/MpiSize
       MpiR = 0
       if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)
    else
       MpiQ = Dim
       MpiR = 0
    endif
#else
    MpiQ = Dim
    MpiR = 0
#endif
    !
    vecDim=MpiQ + MpiR
    !
  end function vecDim_Hv_sector




end MODULE ED_HAMILTONIAN
