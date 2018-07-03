!########################################################################
!PURPOSE  : Build the impurity Hamiltonian
!|ImpUP,(2ImpUP),BathUP;,ImpDW,(2ImpDW),BathDW >
! |1,2;3...Ns>_UP * |Ns+1,Ns+2;Ns+3,...,2*Ns>_DOWN
!########################################################################
MODULE ED_HAMILTONIAN_MATVEC
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

  !>build sparse hamiltonian of the sector
  public  :: setup_Hv_sector
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
  !


  !> MPI local variables (shared)
#ifdef _MPI
  integer                            :: MpiComm=MPI_UNDEFINED
#else
  integer                            :: MpiComm=0
#endif
  logical                            :: MpiStatus=.false.
  integer                            :: MpiIerr
  integer                            :: MpiRank=0
  integer                            :: MpiSize=1
  integer                            :: mpiQ=1
  integer                            :: mpiR=0
  !
  integer                            :: Hsector=0
  logical                            :: Hstatus=.false.,Pbool=.true.
  type(sector_map)                   :: H,Hup,Hdw
  !
  integer,dimension(:,:),allocatable :: columns_range





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


  subroutine ed_hamiltonian_matvec_set_ColumnsRange()
#ifdef _MPI
    if(MpiStatus)then
       if(allocated(columns_range))deallocate(columns_range)
       allocate(columns_range(0:MpiSize-1,2))
       call sp_columns_range_matrix(MpiComm,spH0,columns_range(0:,:))
    end if
#endif
  end subroutine ed_hamiltonian_matvec_set_ColumnsRange


  subroutine ed_hamiltonian_matvec_del_ColumnsRange
    if(allocated(columns_range))deallocate(columns_range)
  end subroutine ed_hamiltonian_matvec_del_ColumnsRange







  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine setup_Hv_sector(isector,Hmat)
    integer                            :: isector,SectorDim,irank
    complex(8),dimension(:,:),optional :: Hmat
    !
    Hsector = isector
    Hstatus = .true.
    !
    call build_sector(isector,H)
    !
    SectorDim=getdim(isector)
    if(present(Hmat))then
       if(any( shape(Hmat) /= [SectorDim,SectorDim]))&
            stop "setup_Hv_sector ERROR: size(Hmat) != SectorDim**2"
       call ed_buildH_c(isector,Hmat)
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)               !H is to be stored 
       call ed_buildH_c(isector)
       ! call ed_hamiltonian_matvec_set_ColumnsRange
    case (.false.)
       !nothing to be done: if MPI: get columns range from a dryrun
       ! if(MpiStatus)call ed_buildH_c(isector,dryrun=.true.) !if MPI: get columns range anyway, dryrun
       ! if(MpiStatus)call ed_hamiltonian_matvec_set_ColumnsRange
    end select
    !
  end subroutine setup_Hv_sector



  subroutine delete_Hv_sector(isector)
    integer :: isector
    !
    call delete_sector(isector,H)
    !
    Hsector = 0
    Hstatus = .false.
    !
    if(spH0%status)then
       if(MpiStatus)then          
          call sp_delete_matrix(MpiComm,spH0)
       else
          call sp_delete_matrix(spH0)
       endif
    endif
    !
    ! call ed_hamiltonian_matvec_del_ColumnsRange
    !
  end subroutine delete_Hv_sector








  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine ed_buildH_c(isector,Hmat,dryrun)
    integer,intent(in)                     :: isector
    complex(8),dimension(:,:),optional     :: Hmat
    logical,optional                       :: dryrun
    complex(8),dimension(:,:),allocatable  :: Hredux
    logical                                :: dryrun_
    integer,dimension(Nlevels)             :: ib
    integer,dimension(Ns)                  :: ibup,ibdw
    integer                                :: dim,dimUp,dimDw
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
    if(.not.Hstatus)stop "ed_buildH_c ERROR: HVsector NOT set"
    dryrun_=.false.;if(present(dryrun))dryrun_=dryrun
    !
    dim = getdim(isector)
    !
    if(MpiStatus)then
       call sp_init_matrix(MpiComm,spH0,dim,dryrun=dryrun_)
    else
       call sp_init_matrix(spH0,dim) !serial mode: dryrun is just not needed
    endif
    !
    mpiQ = dim/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
    !
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
    !Get diagonal hybridization
    diag_hybr=zero
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
             enddo
          enddo
       enddo
    else
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             enddo
          enddo
       enddo
    endif
    !
    !-----------------------------------------------!
    include "ED_HAMILTONIAN_MATVEC/build_h.f90"
    !-----------------------------------------------!
    !
    if(present(Hmat))then
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0,Hmat)
          call sp_delete_matrix(MpiComm,spH0)
       else
          call sp_dump_matrix(spH0,Hmat)
          call sp_delete_matrix(spH0)
       endif
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
    integer                      :: Nloc
    complex(8),dimension(Nloc)   :: v
    complex(8),dimension(Nloc)   :: Hv
    integer                      :: i
    type(sparse_element),pointer :: c
    Hv=zero
    do i=1,Nloc
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%cval*v(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine spMatVec_cc


#ifdef _MPI
  subroutine spMatVec_mpi_cc(Nloc,v,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: v,Hv
    !
    integer                             :: N
    complex(8),dimension(:),allocatable :: vin
    integer                             :: i,j,irank,Imin,Imax
    !MPI
    integer                             :: MpiIerr
    integer                             :: MpiRank
    integer                             :: MpiSize
    logical                             :: MpiMaster
    integer                             :: MpiQ
    integer                             :: MpiR,MpiShift
    integer,allocatable,dimension(:)    :: RCounts,SCounts,Offset
    type(sparse_element),pointer        :: c
    !
    !
    if(MpiComm==MPI_UNDEFINED)stop "spMatVec_mpi_cc ERRROR: MpiComm = MPI_UNDEFINED"
    ! if(.not.allocated(columns_range))stop "spMatVec_mpi_cc ERRROR: columns_range NOT allocated"
    MpiRank   = get_rank_MPI(MpiComm)
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    !Get N_total by summing Nloc over all procs 
    N=0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    !Get the chunks again (must equate Nloc for each proc)
    MpiQ = N/MpiSize
    MpiR = 0
    if(MpiRank==MpiSize-1)MpiR=mod(N,MpiSize)
    MpiShift = MpiRank*MpiQ
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    Hv=zero
    do i=1,Nloc
       c => spH0%loc(i)%root%next
       local: do while(associated(c))
          Hv(i) = Hv(i) + c%cval*v(c%col-MpiShift)
          c => c%next
       end do local
    end do
    !
    !Evaluate the non-local contribution.
    !Retrieve the component in the range Min_Column : Max_Column only by looking at ColumnsRange array
    allocate(RCounts(0:MpiSize-1)) ; RCounts(0:)=0
    allocate(SCounts(0:MpiSize-1)) ; SCounts(0:)=0
    allocate(Offset(0:MpiSize-1))  ; Offset(0:)=0
    !
    forall(i=0:MpiSize-1)Offset(i) = i*MpiQ
    ! !
    ! Imin = columns_range(MpiRank,1)
    ! Imax = columns_range(MpiRank,2)
    !
    RCounts(0:)=MpiQ ; RCounts(MpiSize-1)=MpiQ+mod(N,MpiSize)
    RCounts(MpiRank)  = 0
    !
    allocate(vin(N)) ; vin = zero
    !
    do irank=0,MpiSize-1
       call MPI_Gatherv(v,Nloc,MPI_Double_Complex,Vin,RCounts,Offset,MPI_Double_Complex,Irank,MpiComm,MpiIerr)
    enddo
    !
    do i=1,Nloc                 !==spH0%Nrow
       c => spH0%row(i)%root%next
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%cval*vin(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
    deallocate(SCounts,RCounts,Offset,Vin)
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
    integer                                :: dim,dimUp,dimDw
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
    dim=getdim(isector)
    if(Nloc/=dim)stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    mpiQ = dim/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
    !
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
    !Get diagonal hybridization
    diag_hybr=zero
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
             enddo
          enddo
       enddo
    else
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             enddo
          enddo
       enddo
    endif
    !
    Hv=zero
    !-----------------------------------------------!
    include "ED_HAMILTONIAN_MATVEC/build_hxv.f90"
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
    integer,allocatable,dimension(:)       :: SendCounts,Displs
    integer                                :: isector
    integer,dimension(Nlevels)             :: ib
    integer,dimension(Ns)                  :: ibup,ibdw
    integer                                :: dim,dimUp,dimDw
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
    integer                                :: Irank,Imax,Imin
    integer,allocatable,dimension(:)       :: Counts,Offset
    !
    if(.not.Hstatus)stop "directMatVec_MPI_cc ERROR: Hsector NOT set"
    if(MpiComm==MPI_UNDEFINED)stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    !
    isector=Hsector
    !
    dim=getdim(isector)
    !
    !Get diagonal hybridization
    diag_hybr=zero
    if(bath_type/="replica")then
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
             enddo
          enddo
       enddo
    else
       do ibath=1,Nbath
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
             enddo
          enddo
       enddo
    endif
    !
    !Get N_total by summing Nloc over all procs 
    N=0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(N/=dim)stop "directMatVec_MPI_cc ERROR: N != dim(isector)"
    !    
    mpiQ = N/MpiSize
    mpiR = 0
    if(MpiRank==(MpiSize-1))mpiR=mod(N,MpiSize)
    !
    ishift      = MpiRank*mpiQ
    first_state = MpiRank*mpiQ + 1
    last_state  = (MpiRank+1)*mpiQ + mpiR
    !
    !
    ! allocate(vin(N))
    ! allocate(SendCounts(0:MpiSize-1),displs(0:MpiSize-1))
    ! vin                   = zero
    ! SendCounts(0:)        = mpiQ
    ! SendCounts(MpiSize-1) = mpiQ+mod(N,MpiSize)
    ! forall(i=0:MpiSize-1)Displs(i)=i*mpiQ
    ! call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,vin,SendCounts,Displs,MPI_Double_Complex,MpiComm,MpiIerr)
    allocate(vin(N))
    vin = zero
    !
    !Each thread copy its own copy at the right chunk of the input vector.
    !The remaining components are obtained thru successive gatherv
    !Local contribution:
    do i=first_state,last_state
       vin(i) = v(i-Ishift)
    enddo
    !
    !Get the non-local contribution
    allocate(Counts(0:MpiSize-1)) ; Counts(0:)=0
    allocate(Offset(0:MpiSize-1)) ; Offset(0:)=0
    !
    forall(i=0:MpiSize-1)Offset(i) = i*MpiQ
    !
    ! Imin = columns_range(MpiRank,1)
    ! Imax = columns_range(MpiRank,2)
    ! !
    Counts(0:) = MpiQ ; Counts(MpiSize-1) = MpiQ+mod(N,MpiSize)
    Counts(MpiRank)  = 0
    !
    do irank=0,MpiSize-1
       call MPI_Gatherv(v,Nloc,MPI_Double_Complex,Vin,Counts,Offset,MPI_Double_Complex,Irank,MpiComm,MpiIerr)
    enddo
    !
    Hv=zero
    !
    !-----------------------------------------------!
    include "ED_HAMILTONIAN_MATVEC/build_hxv.f90"
    !-----------------------------------------------!
    !
  end subroutine directMatVec_MPI_cc
#endif








end MODULE ED_HAMILTONIAN_MATVEC





























! subroutine spMatVec_mpi_cc(Nloc,v,Hv)
!   integer                             :: Nloc
!   complex(8),dimension(Nloc)          :: v,Hv
!   !
!   integer                             :: N
!   complex(8),dimension(:),allocatable :: vin
!   integer                             :: i,j,irank,Imin,Imax
!   !MPI
!   integer                             :: MpiIerr
!   integer                             :: MpiRank
!   integer                             :: MpiSize
!   logical                             :: MpiMaster
!   integer                             :: MpiQ
!   integer                             :: MpiR,MpiShift,ierr
!   integer,allocatable,dimension(:)    :: RCounts,SCounts,Offset
!   type(sparse_element),pointer        :: c
!   !
!   character(len=MPI_MAX_ERROR_STRING ) :: string
!   !
!   if(MpiComm==MPI_UNDEFINED)stop "spMatVec_mpi_cc ERRROR: MpiComm = MPI_UNDEFINED"
!   if(.not.allocated(columns_range))stop "spMatVec_mpi_cc ERRROR: columns_range NOT allocated"
!   MpiRank   = get_rank_MPI(MpiComm)
!   MpiSize   = get_size_MPI(MpiComm)
!   MpiMaster = get_master_MPI(MpiComm)
!   !
!   !Get N_total by summing Nloc over all procs 
!   N=0
!   call AllReduce_MPI(MpiComm,Nloc,N)
!   !
!   !Get the chunks again (must equate Nloc for each proc)
!   MpiQ = N/MpiSize
!   MpiR = 0
!   if(MpiRank==MpiSize-1)MpiR=mod(N,MpiSize)
!   MpiShift = MpiRank*MpiQ
!   !



!   !Evaluate the local contribution: Hv_loc = Hloc*v
!   Hv=zero
!   do i=1,Nloc
!      c => spH0%loc(i)%root%next
!      local: do while(associated(c))
!         Hv(i) = Hv(i) + c%cval*v(c%col-MpiShift)
!         c => c%next
!      end do local
!   end do
!   !
!   !Evaluate the non-local contribution.
!   !Retrieve the component in the range Min_Column : Max_Column only by looking at ColumnsRange array
!   allocate(RCounts(0:MpiSize-1)) ; RCounts(0:)=0
!   allocate(SCounts(0:MpiSize-1)) ; SCounts(0:)=0
!   allocate(Offset(0:MpiSize-1))  ; Offset(0:)=0
!   !
!   forall(i=0:MpiSize-1)Offset(i) = i*MpiQ
!   !
!   ! Counts(0:)=MpiQ
!   ! Counts(MpiSize-1)=MpiQ+mod(N,MpiSize)
!   ! Counts(MpiRank)  = 0
!   !
!   Imin = columns_range(MpiRank,1)
!   Imax = columns_range(MpiRank,2)
!   !
!   RCounts(Imin:Imax) = MpiQ ; if(Imax==MpiSize-1)RCounts(MpiSize-1) = MpiQ+mod(N,MpiSize)
!   RCounts(MpiRank)   = 0
!   !
!   Scounts(Imin:Imax) = MpiQ ; if(Imax==MpiSize-1)Scounts(MpiSize-1) = MpiQ + mod(N,MpiSize)
!   SCounts(MpiRank)   = 0
!   !
!   ! forall(i=Imin:Imax)Offset(i) = i*MpiQ   
!   !
!   !>DEBUG
!   if(Pbool)then
!      Pbool=.false.

!      if(MpiMaster)write(*,"(10A10)")"mpiRank","mpi_Q","mpi_R","mpi_CHunk","Nloc","N"
!      do irank=0,MpiSize-1
!         if(MpiRank==irank)then
!            write(*,"(10I10)")mpirank,mpiQ,mpiR,mpiQ+mpiR,Nloc,N
!            call MPI_Barrier(MpiComm,mpiierr)
!         endif
!         call MPI_Barrier(MpiComm,mpiierr)
!      enddo
!      call MPI_Barrier(MpiComm,mpiierr)

!      do irank=0,MpiSize-1
!         if(MpiRank==irank)then
!            print*,"rank,SCounts=",irank,">",SCounts,"|",columns_range(IRank,1),columns_range(IRank,2)
!            print*,"rank,RCounts=",irank,">",RCounts,"|",columns_range(IRank,1),columns_range(IRank,2)             
!            call MPI_Barrier(MpiComm,MpiIerr)
!         endif
!         call MPI_Barrier(MpiComm,MpiIerr)
!      enddo
!      if(MpiMaster)print*,""
!      call MPI_Barrier(MpiComm,MpiIerr)
!   endif
!   !<DEBUG
!   allocate(vin(N)) ; vin = zero

!   do irank=0,MpiSize-1
!      call MPI_Gatherv(v,Nloc,MPI_Double_Complex,Vin,RCounts,Offset,MPI_Double_Complex,Irank,MpiComm,MpiIerr)
!      call MPI_ERROR_STRING(MpiIerr,string,i,ierr)
!   enddo
!   !
!   do i=1,Nloc                 !==spH0%Nrow
!      c => spH0%row(i)%root%next
!      matmul: do while(associated(c))
!         Hv(i) = Hv(i) + c%cval*vin(c%col)
!         c => c%next
!      end do matmul
!   end do
!   nullify(c)
!   deallocate(SCounts,RCounts,Offset,Vin)
! end subroutine spMatVec_mpi_cc

! subroutine spMatVec_mpi_cc(Nloc,v,Hv)
!   integer                             :: Nloc
!   complex(8),dimension(Nloc)          :: v
!   complex(8),dimension(Nloc)          :: Hv
!   integer                             :: i
!   integer                             :: N
!   complex(8),dimension(:),allocatable :: vin
!   integer,allocatable,dimension(:)    :: SendCounts,Displs
!   type(sparse_element),pointer        :: c
!   N=0
!   if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_cc ERRROR: MpiComm = MPI_UNDEFINED"
!   call MPI_AllReduce(Nloc,N,1,MPI_Integer,MPI_Sum,MpiComm,MpiIerr)
!   MpiSize = get_Size_MPI(MpiComm)
!   mpiQ = get_Q_MPI(MpiComm,N)
!   mpiR = get_R_MPI(MpiComm,N)
!   allocate(vin(N))
!   allocate(SendCounts(0:MpiSize-1),displs(0:MpiSize-1))
!   vin                   = zero
!   SendCounts(0:)        = mpiQ
!   SendCounts(MpiSize-1) = mpiQ+mod(N,MpiSize)
!   forall(i=0:MpiSize-1)Displs(i)=i*mpiQ
!   call MPI_Allgatherv(v(1:Nloc),Nloc,MPI_Double_Complex,vin,SendCounts,Displs,MPI_Double_Complex,MpiComm,MpiIerr)
!   call MPI_Bcast(vin,N,MPI_Double_Complex,0,MpiComm,MpiIerr)
!   Hv=zero
!   do i=1,Nloc                 !==spH0%Nrow
!      c => spH0%row(i)%root%next       
!      matmul: do while(associated(c))
!         Hv(i) = Hv(i) + c%cval*vin(c%col)
!         c => c%next
!      end do matmul
!   end do
!   nullify(c)
! end subroutine spMatVec_mpi_cc










!   !+------------------------------------------------------------------+
!   !>NORMAL CASE
!   !+------------------------------------------------------------------+
!   subroutine build_H_normal_c(isector,Hmat)
!     complex(8),dimension(:,:),optional     :: Hmat
!     complex(8),dimension(:,:),allocatable  :: Hredux
!     integer                                :: isector
!     type(sector_map)                       :: H,Hup,Hdw
!     integer,dimension(Nlevels)             :: ib
!     integer,dimension(Ns)                  :: ibup,ibdw
!     integer                                :: dim,dimUp,dimDw
!     integer                                :: i,iup,idw
!     integer                                :: m,mup,mdw
!     integer                                :: ishift,ishift_up,ishift_dw
!     integer                                :: j,ms,impi
!     integer                                :: iorb,jorb,ispin,jspin,ibath
!     integer                                :: kp,k1,k2,k3,k4
!     integer                                :: alfa,beta
!     real(8)                                :: sg1,sg2,sg3,sg4
!     real(8),dimension(Norb)                :: nup,ndw
!     complex(8)                             :: htmp,htmpup,htmpdw
!     complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!     logical                                :: Jcondition
!     integer                                :: first_state,last_state
!     integer                                :: first_state_up,last_state_up
!     integer                                :: first_state_dw,last_state_dw
!     !
!     !
!     call build_sector(isector,H)
!     !
!     if(spH0%status)call sp_delete_matrix(spH0) 
!     !
!     dim=getdim(isector)
!     mpiQ = dim/MpiSize
!     mpiR = 0
!     if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
!     call sp_init_matrix(spH0,mpiQ + mpiR)
!     ishift      = MpiRank*mpiQ
!     first_state = MpiRank*mpiQ + 1
!     last_state  = (MpiRank+1)*mpiQ + mpiR
!     !
!     !
!     !Get diagonal hybridization
!     diag_hybr=zero
!     if(bath_type/="replica")then
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),00d0)
!              enddo
!           enddo
!        enddo
!     else
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
!              enddo
!           enddo
!        enddo
!     endif
!     !
!     !-----------------------------------------------!
!     include "ED_HAMILTONIAN/build_h_normal.f90"
!     !-----------------------------------------------!
!     !
!     deallocate(H%map)
!     !
!     if(present(Hmat))then
!        if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_normal_c ERROR: size(Hmat) != dim**2"
!        if(MpiStatus)then
!           allocate(Hredux(dim,dim));Hredux=zero
!           call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
! #ifdef _MPI
!           call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
! #endif
!        else
!           call sp_dump_matrix(spH0,Hmat)
!        endif
!     endif
!     !
!   end subroutine build_H_normal_c


!   !+------------------------------------------------------------------+
!   !>SUPERC CASE
!   !+------------------------------------------------------------------+
!   !DOUBLE COMPLEX
!   subroutine build_H_superc_c(isector,Hmat)
!     complex(8),dimension(:,:),optional     :: Hmat
!     complex(8),dimension(:,:),allocatable  :: Hredux
!     integer                                :: isector
!     type(sector_map)                       :: H,Hup,Hdw
!     integer,dimension(Nlevels)             :: ib
!     integer,dimension(Ns)                  :: ibup,ibdw
!     integer                                :: dim,dimUp,dimDw
!     integer                                :: i,iup,idw
!     integer                                :: m,mup,mdw
!     integer                                :: ishift,ishift_up,ishift_dw
!     integer                                :: j,ms,impi
!     integer                                :: iorb,jorb,ispin,jspin,ibath
!     integer                                :: kp,k1,k2,k3,k4
!     integer                                :: alfa,beta
!     real(8)                                :: sg1,sg2,sg3,sg4
!     real(8),dimension(Norb)                :: nup,ndw
!     complex(8)                             :: htmp,htmpup,htmpdw
!     complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!     logical                                :: Jcondition
!     integer                                :: first_state,last_state
!     integer                                :: first_state_up,last_state_up
!     integer                                :: first_state_dw,last_state_dw
!     !
!     !
!     call build_sector(isector,H)
!     !
!     if(spH0%status)call sp_delete_matrix(spH0) 
!     !
!     dim=getdim(isector)
!     mpiQ = dim/MpiSize
!     mpiR = 0
!     if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
!     call sp_init_matrix(spH0,mpiQ + mpiR)
!     ishift      = MpiRank*mpiQ
!     first_state = MpiRank*mpiQ + 1
!     last_state  = (MpiRank+1)*mpiQ + mpiR
!     !
!     !
!     !Get diagonal hybridization
!     diag_hybr=zero
!     if(bath_type/="replica")then
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),0d0)
!              enddo
!           enddo
!        enddo
!     else
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
!              enddo
!           enddo
!        enddo
!     endif
!     !
!     !-----------------------------------------------!
!     include "ED_HAMILTONIAN/build_h_superc.f90"
!     !-----------------------------------------------!
!     !
!     deallocate(H%map)
!     !
!     if(present(Hmat))then
!        if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_superc_c ERROR: size(Hmat) != dim**2"
!        if(MpiStatus)then
!           allocate(Hredux(dim,dim));Hredux=zero
!           call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
! #ifdef _MPI
!           call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
! #endif
!        else
!           call sp_dump_matrix(spH0,Hmat)
!        endif
!     endif
!     !
!   end subroutine build_H_superc_c


!   !+------------------------------------------------------------------+
!   !>NONSU2 CASE
!   !+------------------------------------------------------------------+
!   subroutine build_H_nonsu2_c(isector,Hmat)
!     complex(8),dimension(:,:),optional     :: Hmat
!     complex(8),dimension(:,:),allocatable  :: Hredux
!     integer                                :: isector
!     type(sector_map)                       :: H,Hup,Hdw
!     integer,dimension(Nlevels)             :: ib
!     integer,dimension(Ns)                  :: ibup,ibdw
!     integer                                :: dim,dimUp,dimDw
!     integer                                :: i,iup,idw
!     integer                                :: m,mup,mdw
!     integer                                :: ishift,ishift_up,ishift_dw
!     integer                                :: j,ms,impi
!     integer                                :: iorb,jorb,ispin,jspin,ibath
!     integer                                :: kp,k1,k2,k3,k4
!     integer                                :: alfa,beta
!     real(8)                                :: sg1,sg2,sg3,sg4
!     real(8),dimension(Norb)                :: nup,ndw
!     complex(8)                             :: htmp,htmpup,htmpdw
!     complex(8),dimension(Nspin,Norb,Nbath) :: diag_hybr
!     logical                                :: Jcondition
!     integer                                :: first_state,last_state
!     integer                                :: first_state_up,last_state_up
!     integer                                :: first_state_dw,last_state_dw    
!     !
!     !
!     call build_sector(isector,H)
!     !
!     if(spH0%status)call sp_delete_matrix(spH0) 
!     !
!     dim=getdim(isector)
!     mpiQ = dim/MpiSize
!     mpiR = 0
!     if(MpiRank==(MpiSize-1))mpiR=mod(dim,MpiSize)
!     call sp_init_matrix(spH0,mpiQ + mpiR)
!     ishift      = MpiRank*mpiQ
!     first_state = MpiRank*mpiQ + 1
!     last_state  = (MpiRank+1)*mpiQ + mpiR
!     !
!     !
!     !Get diagonal hybridization
!     diag_hybr=zero
!     if(bath_type/="replica")then
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dcmplx(dmft_bath%v(ispin,iorb,ibath),0d0)
!              enddo
!           enddo
!        enddo
!     else
!        do ibath=1,Nbath
!           do ispin=1,Nspin
!              do iorb=1,Norb
!                 diag_hybr(ispin,iorb,ibath)=dmft_bath%vr(ibath)
!              enddo
!           enddo
!        enddo
!     endif
!     !
!     !-----------------------------------------------!
!     include "ED_HAMILTONIAN/build_h_nonsu2.f90"
!     !-----------------------------------------------!
!     !
!     deallocate(H%map)
!     !
!     if(present(Hmat))then
!        if(size(Hmat,1)/=dim.OR.size(Hmat,2)/=dim)stop "build_H_nonsu2_c ERROR: size(Hmat) != dim**2"
!        if(MpiStatus)then
!           allocate(Hredux(dim,dim));Hredux=zero
!           call sp_dump_matrix(spH0,Hredux(first_state:last_state,:))
! #ifdef _MPI
!           call MPI_AllReduce(Hredux,Hmat,dim*dim,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
! #endif
!        else
!           call sp_dump_matrix(spH0,Hmat)
!        endif
!     endif
!     !
!   end subroutine build_H_nonsu2_c
!   !+------------------------------------------------------------------+
!   !+------------------------------------------------------------------+
!   !+------------------------------------------------------------------+
