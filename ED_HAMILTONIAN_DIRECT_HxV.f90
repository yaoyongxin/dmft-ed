! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT
MODULE ED_HAMILTONIAN_DIRECT_HxV
  USE ED_HAMILTONIAN_SHARED
  implicit none
  private


  !>Sparse Mat-Vec direct on-the-fly product
  public  :: directMatVec_cc
#ifdef _MPI
  public  :: directMatVec_MPI_cc
#endif




contains



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
    integer                                :: j,ms
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
                diag_hybr(ispin,iorb,ibath)=dmft_bath%t(ibath)
             endif
          enddo
       enddo
    enddo
    !
    Hv=zero
    !-----------------------------------------------!
    states: do j=MpiIstart,MpiIend
       m    = H%map(j)
       ib   = bdecomp(m,2*Ns)
       !
       do iorb=1,Norb
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !
       !
       !IMPURITY  HAMILTONIAN
       include "ED_HAMILTONIAN/direct/HxVimp.f90"
       !
       !LOCAL INTERACTION
       include "ED_HAMILTONIAN/direct/HxVint.f90"
       !
       !BATH HAMILTONIAN
       include "ED_HAMILTONIAN/direct/HxVbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "ED_HAMILTONIAN/direct/HxVimp_bath.f90"
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
    integer                                :: j,ms
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
    integer :: mpiIerr
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
                diag_hybr(ispin,iorb,ibath)=dmft_bath%t(ibath)
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
    states: do j=MpiIstart,MpiIend
       m  = H%map(j)
       ib = bdecomp(m,2*Ns)
       !
       do iorb=1,Norb
          nup(iorb)=dble(ib(iorb))
          ndw(iorb)=dble(ib(iorb+Ns))
       enddo
       !
       !
       !IMPURITY  HAMILTONIAN
       include "ED_HAMILTONIAN/direct/HxVimp.f90"
       !
       !LOCAL INTERACTION
       include "ED_HAMILTONIAN/direct/HxVint.f90"
       !
       !BATH HAMILTONIAN
       include "ED_HAMILTONIAN/direct/HxVbath.f90"
       !
       !IMPURITY- BATH HYBRIDIZATION
       include "ED_HAMILTONIAN/direct/HxVimp_bath.f90"
    enddo states
    !-----------------------------------------------!
    !
  end subroutine directMatVec_MPI_cc
#endif





end MODULE ED_HAMILTONIAN_DIRECT_HXV
