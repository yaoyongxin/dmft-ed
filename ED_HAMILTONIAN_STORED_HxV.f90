! > BUILD STORED SPARSE HAMILTONIAN of the SECTOR
MODULE ED_HAMILTONIAN_STORED_HxV
  USE ED_HAMILTONIAN_SHARED
  implicit none
  private

  !>Sparse Matric constructors
  public :: ed_buildH_c


  !>Sparse Mat-Vec product using stored sparse matrix
  public  :: spMatVec_cc
#ifdef _MPI
  public  :: spMatVec_MPI_cc
#endif





contains



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
    integer                                :: korb,lorb
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
    integer                                :: isite,jsite,hopndx,ineig
    integer                                :: distance,firstNeig,BosonExp
    real(8)                                :: n_i,n_j
    integer(16)                            :: m_8,k1_8,k2_8
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
                diag_hybr(ispin,iorb,ibath)=dmft_bath%t(ibath)
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
    !
    !-----------------------------------------------!
    !
    !IMPURITY  HAMILTONIAN
    if(plaquette)then
       !
       !STANDARD FOR FERMIONIC (U=inf) MODEL
       BosonExp=1
       firstNeig=1
       !
       if(HardCoreBoson.ne.0)then
          BosonExp=0
          !
          !LOCAL HARD-CORE CONDITION - ALL DIMENSIONS
          if(HardCoreBoson.eq.1)firstNeig=1
          !
          !EXTENDED HARD-CORE CONDITION - 1D
          if((HardCoreBoson.eq.2).and.(Nbath.ne.1))firstNeig=2
          !
          !EXTENDED HARD-CORE CONDITION - 2D
          if((HardCoreBoson.eq.2).and.(Nbath.ne.1))firstNeig=3
          !
       endif
       include "ED_HAMILTONIAN/stored/Hplaquette.f90"
    else
       include "ED_HAMILTONIAN/stored/Himp.f90"
       !
       !LOCAL INTERACTION
       if (Utensor) then
           include "ED_HAMILTONIAN/stored/Hint_tensor.f90"
       else
           include "ED_HAMILTONIAN/stored/Hint.f90"
       endif
       !
       !BATH HAMILTONIAN
       include "ED_HAMILTONIAN/stored/Hbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "ED_HAMILTONIAN/stored/Himp_bath.f90"
    endif

    !
    !
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
    integer                             :: i,j,mpiIerr
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











end MODULE ED_HAMILTONIAN_STORED_HXV
