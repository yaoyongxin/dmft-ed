!########################################################################
!PURPOSE  : Build the impurity Hamiltonian
!|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
! |1,2;3...Ns>_UP * |Ns+1,Ns+2;Ns+3,...,2*Ns>_DOWN
!########################################################################
MODULE ED_HAMILTONIAN_MATVEC
  USE SF_MISC,    only: assert_shape
  USE SF_CONSTANTS,only:zero
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_SETUP
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private

  !>Build sparse hamiltonian of the sector
  public  :: vecDim_Hv_sector
  public  :: build_Hv_sector
  public  :: delete_Hv_sector
  !
  !
  !>Sparse Mat-Vec product using stored sparse matrix 
  public  :: spMatVec_cc
#ifdef _MPI
  public  :: spMatVec_MPI_cc
#endif
  !
  !
  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_cc
#ifdef _MPI
  public  :: directMatVec_MPI_cc
#endif
  !
  !
  !> Related auxiliary routines:
  public  :: ed_hamiltonian_matvec_set_MPI
  public  :: ed_hamiltonian_matvec_del_MPI


  !> MPI local variables (shared)
#ifdef _MPI
  integer          :: MpiComm=MPI_UNDEFINED
#else
  integer          :: MpiComm=0
#endif
  logical          :: MpiStatus=.false.
  integer          :: MpiIerr
  integer          :: MpiRank=0
  integer          :: MpiSize=1
  logical          :: MpiMaster=.true.
  integer          :: MpiQ=1
  integer          :: MpiR=0

  integer          :: MpiIstart
  integer          :: MpiIend
  integer          :: MpiIshift
  !
  integer          :: Hsector=0
  logical          :: Hstatus=.false.
  type(sector_map) :: H






contains


  !####################################################################
  !                        AUXILIARY ROUTINES
  !####################################################################
  subroutine ed_hamiltonian_matvec_set_MPI(comm_)
#ifdef _MPI
    integer :: comm_
    MpiComm = comm_
    MpiStatus=.true.
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
    MpiMaster = get_Master_MPI(MpiComm)
#else
    integer,optional :: comm_
#endif
  end subroutine ed_hamiltonian_matvec_set_MPI


  subroutine ed_hamiltonian_matvec_del_MPI()
#ifdef _MPI
    MpiComm = MPI_UNDEFINED
#else
    MpiComm = 0
#endif
    MpiStatus=.false.
    MpiRank=0
    MpiSize=1
    MpiQ=1
    MpiR=0
  end subroutine ed_hamiltonian_matvec_del_MPI






  !####################################################################
  !                 MAIN ROUTINES: BUILD/DELETE SECTOR
  !####################################################################
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
    MpiIshift = MpiRank*mpiQ
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
             write(LOGfile,*)MpiRank,MpiQ,MpiR,MpiIstart,MpiIend,MpiIstart-MpiIend+1
          endif
       enddo
       call Barrier_MPI(MpiComm)
    endif
#endif
    !
    !
    if(present(Hmat))then
       call ed_buildh_c(isector,Hmat)
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)
       call ed_buildh_c(isector)
    case (.false.)
       !nothing to be done: direct matrix-vector product
    end select
    !
  end subroutine build_Hv_sector


  subroutine delete_Hv_sector()
    integer :: iud
    call delete_sector(Hsector,H)
    Hstatus=.false.
    !
    mpiQ = 0
    mpiR = 0
    MpiIshift = 0
    MpiIstart = 0
    MpiIend   = 0
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
  end subroutine delete_Hv_sector










  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine ed_buildH_c(isector,Hmat)
    integer                                :: isector
    complex(8),dimension(:,:),optional     :: Hmat
    integer,dimension(Nlevels)             :: ib
    integer,dimension(Ns)                  :: ibup,ibdw
    integer                                :: Dim
    integer                                :: i,iup,idw
    integer                                :: m,mup,mdw
    integer                                :: ishift,ishift_up,ishift_dw
    integer                                :: j,ms,impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)                :: nup,ndw
    complex(8)                             :: htmp,htmpup,htmpdw
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    integer                                :: first_state_up,last_state_up
    integer                                :: first_state_dw,last_state_dw
    !
    if(.not.Hstatus)stop "ed_buildH_c ERROR: Hsector NOT set"
    isector=Hsector
    !
    Dim = getdim(isector)
    !
    if(present(Hmat))call assert_shape(Hmat,[Dim,Dim],"ed_buildh_main","Hmat")
    !
    !
    !Get diagonal hybridization
    diag_hybr=zero
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             if(bath_type/="replica")then
                diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
             else
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             endif
          enddo
       enddo
    enddo
    !
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0,mpiIstart,mpiIend,mpiIshift)
       call sp_init_matrix(MpiComm,spH0,Dim)
    else
       call sp_init_matrix(spH0,Dim)
    endif
#else
    call sp_init_matrix(spH0,Dim)
#endif


    !-----------------------------------------------!
    states: do i=MpiIstart,MpiIend
       m = H%map(i)
       impi = i-MpiIshift
       ib = bdecomp(m,2*Ns)
       !
       do iorb=1,Norb
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !
       !
       !IMPURITY  HAMILTONIAN
       include "ED_HAMILTONIAN_MATVEC/stored/Himp.f90"
       !
       !LOCAL INTERACTION
       include "ED_HAMILTONIAN_MATVEC/stored/Hint.f90"
       !
       !BATH HAMILTONIAN
       include "ED_HAMILTONIAN_MATVEC/stored/Hbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "ED_HAMILTONIAN_MATVEC/stored/Himp_bath.f90"
       !
       !
    enddo states
    !-----------------------------------------------!
    !
    !
    if(present(Hmat))then
#ifdef _MPI
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0,Hmat)
       else
          call sp_dump_matrix(spH0,Hmat)
       endif
#else
       call sp_dump_matrix(spH0,Hmat)
#endif          
    endif
    !
  end subroutine ed_buildH_c














  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial cmplx(H)*cmplx(V)
  ! - MPI cmplx(H)*cmplx(V)
  !+------------------------------------------------------------------+
  subroutine spMatVec_cc(Nloc,v,Hv)
    integer                         :: Nloc
    complex(8),dimension(Nloc)      :: v
    complex(8),dimension(Nloc)      :: Hv
    integer                         :: i,j
    Hv=zero
    do i=1,Nloc
       matmul: do j=1,spH0%row(i)%Size
          Hv(i) = Hv(i) + spH0%row(i)%vals(j)*v(spH0%row(i)%cols(j))
       end do matmul
    end do
  end subroutine spMatVec_cc


#ifdef _MPI
  subroutine spMatVec_mpi_cc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    integer                             :: i,j
    integer                             :: N,MpiShift
    complex(8),dimension(:),allocatable :: vin
    integer,allocatable,dimension(:)    :: Counts,Offset
    !
    !
    if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_mpi_cc ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_cc ERROR: MpiStatus = F"
    !
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
    !
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    MpiShift = spH0%Ishift
    Hv=0d0
    do i=1,Nloc
       local: do j=1,spH0%loc(i)%Size
          Hv(i) = Hv(i) + spH0%loc(i)%vals(j)*v(spH0%loc(i)%cols(j)-MpiShift)
       end do local
    end do
    !
    allocate(Counts(0:MpiSize-1)) ; Counts(0:)=0
    allocate(Offset(0:MpiSize-1)) ; Offset(0:)=0
    !
    Counts(0:)        = N/MpiSize
    Counts(MpiSize-1) = N/MpiSize+mod(N,MpiSize)
    !
    do i=1,MpiSize-1
       Offset(i) = Counts(i-1) + Offset(i-1)
    enddo
    !
    allocate(vin(N)) ; vin = zero
    call MPI_Allgatherv(&
         v(1:Nloc),Nloc,MPI_Double_Complex,&
         vin      ,Counts,Offset,MPI_Double_Complex,&
         MpiComm,MpiIerr)
    !
    do i=1,Nloc                 !==spH0%Nrow
       matmul: do j=1,spH0%row(i)%Size
          Hv(i) = Hv(i) + spH0%row(i)%vals(j)*vin(spH0%row(i)%cols(j))
       end do matmul
    end do
    !
  end subroutine spMatVec_mpi_cc
#endif












  !####################################################################
  !            SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
  !####################################################################
  subroutine directMatVec_cc(Nloc,vin,Hv)
    integer                                :: Nloc
    complex(8),dimension(Nloc)             :: vin
    complex(8),dimension(Nloc)             :: Hv
    integer                                :: isector
    integer,dimension(Nlevels)             :: ib
    integer,dimension(Ns)                  :: ibup,ibdw
    integer                                :: Dim
    integer                                :: i,iup,idw
    integer                                :: m,mup,mdw
    integer                                :: ishift,ishift_up,ishift_dw
    integer                                :: j,ms,impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)                :: nup,ndw
    complex(8)                             :: htmp,htmpup,htmpdw
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    integer                                :: first_state_up,last_state_up
    integer                                :: first_state_dw,last_state_dw
    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    Dim = getdim(isector)
    !
    if(Nloc/=dim)stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization
    diag_hybr=zero
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             if(bath_type/="replica")then
                diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
             else
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             endif
          enddo
       enddo
    enddo
    !
    Hv=zero
    !-----------------------------------------------!
    states: do i=MpiIstart,MpiIend
       m = H%map(i)
       impi = i-MpiIshift
       ib = bdecomp(m,2*Ns)
       !
       do iorb=1,Norb
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !
       !
       !IMPURITY  HAMILTONIAN
       include "ED_HAMILTONIAN_MATVEC/direct/HxVimp.f90"
       !
       !LOCAL INTERACTION
       include "ED_HAMILTONIAN_MATVEC/direct/HxVint.f90"
       !
       !BATH HAMILTONIAN
       include "ED_HAMILTONIAN_MATVEC/direct/HxVbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "ED_HAMILTONIAN_MATVEC/direct/HxVimp_bath.f90"
    enddo states
    !-----------------------------------------------!
    !
  end subroutine directMatVec_cc



#ifdef _MPI
  subroutine directMatVec_MPI_cc(Nloc,v,Hv)
    integer                                :: Nloc
    complex(8),dimension(Nloc)             :: v
    complex(8),dimension(Nloc)             :: Hv
    integer                                :: N
    complex(8),dimension(:),allocatable    :: vin
    integer,allocatable,dimension(:)       :: Counts,Offset
    integer                                :: isector
    integer,dimension(Nlevels)             :: ib
    integer,dimension(Ns)                  :: ibup,ibdw
    integer                                :: Dim
    integer                                :: i,iup,idw
    integer                                :: m,mup,mdw
    integer                                :: ishift,ishift_up,ishift_dw
    integer                                :: j,ms,impi
    integer                                :: iorb,jorb,ispin,jspin,ibath
    integer                                :: kp,k1,k2,k3,k4
    integer                                :: alfa,beta
    real(8)                                :: sg1,sg2,sg3,sg4
    real(8),dimension(Norb)                :: nup,ndw
    complex(8)                             :: htmp,htmpup,htmpdw
    complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
    logical                                :: Jcondition
    integer                                :: first_state,last_state
    integer                                :: first_state_up,last_state_up
    integer                                :: first_state_dw,last_state_dw
    !
    if(MpiComm==MPI_UNDEFINED)stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    if(.not.Hstatus)stop "directMatVec_MPI_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    Dim = getdim(isector)
    !
    !
    !Get diagonal hybridization
    diag_hybr=zero
    do ibath=1,Nbath
       do ispin=1,Nspin
          do iorb=1,Norb
             if(bath_type/="replica")then
                diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
             else
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             endif
          enddo
       enddo
    enddo
    !
    N=0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    !Reconstruct Vin and get the displacements for AllGatherV call
    allocate(Counts(0:MpiSize-1)) ; Counts(0:)=0
    allocate(Offset(0:MpiSize-1)) ; Offset(0:)=0
    !
    Counts(0:)        = N/MpiSize
    Counts(MpiSize-1) = N/MpiSize+mod(N,MpiSize)
    !
    do i=1,MpiSize-1
       Offset(i) = Counts(i-1) + Offset(i-1)
    enddo
    !
    allocate(vin(N)); vin  = zero
    call MPI_Allgatherv(&
         v(1:Nloc),Nloc,MPI_Double_Complex,&
         vin,Counts,Offset,MPI_Double_Complex,&
         MpiComm,MpiIerr)
    !
    Hv=zero
    !
    !-----------------------------------------------!
    states: do i=MpiIstart,MpiIend
       m = H%map(i)
       impi = i-MpiIshift
       ib = bdecomp(m,2*Ns)
       !
       do iorb=1,Norb
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !
       !
       !IMPURITY  HAMILTONIAN
       include "ED_HAMILTONIAN_MATVEC/direct/HxVimp.f90"
       !
       !LOCAL INTERACTION
       include "ED_HAMILTONIAN_MATVEC/direct/HxVint.f90"
       !
       !BATH HAMILTONIAN
       include "ED_HAMILTONIAN_MATVEC/direct/HxVbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "ED_HAMILTONIAN_MATVEC/direct/HxVimp_bath.f90"
    enddo states
    !-----------------------------------------------!
    !
  end subroutine directMatVec_MPI_cc
#endif








end MODULE ED_HAMILTONIAN_MATVEC







