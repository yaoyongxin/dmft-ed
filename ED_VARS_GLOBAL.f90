MODULE ED_VARS_GLOBAL
  USE SF_CONSTANTS
  USE ED_SPARSE_MATRIX
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none


  !-------------------- EFFECTIVE BATH STRUCTURE ----------------------!
  type effective_bath
     real(8),dimension(:,:,:),allocatable        :: e     !local energies [Nspin][Norb][Nbath]/[Nspin][1][Nbath]
     real(8),dimension(:,:,:),allocatable        :: d     !SC amplitues   [Nspin][Norb][Nbath]/[Nspin][1][Nbath]
     real(8),dimension(:,:,:),allocatable        :: v     !spin-keep hyb. [Nspin][Norb][Nbath]
     real(8),dimension(:,:,:),allocatable        :: u     !spin-flip hyb. [Nspin][Norb][Nbath]
     real(8),dimension(:),allocatable            :: t     !diagonal hyb.  [Nbath]
     real(8),dimension(:,:),allocatable          :: o     !OpMat coupling [Nop][Nbath]
     complex(8),dimension(:,:,:,:,:),allocatable :: h     !Replica hamilt [Nspin][Nspin][Norb][Norb][Nbath]
     logical(8),dimension(:,:,:,:,:),allocatable :: mask  !impHloc mask   [Nspin][Nspin][Norb][Norb][Re,Im]
     logical                                     :: status=.false.
  end type effective_bath

  !operatorial representation of the replica bath structure
  !INTERNAL USE (accessed thru functions)
  !=========================================================
  complex(8),dimension(:,:,:,:,:),allocatable    :: OpMat !operators      [Nspin][Nspin][Norb][Norb][Nop]




  !---------------- SECTOR-TO-FOCK SPACE STRUCTURE -------------------!
  type sector_map
     integer,dimension(:),allocatable :: map
     logical                          :: status=.false.
  end type sector_map

  interface map_allocate
     module procedure :: map_allocate_scalar
     module procedure :: map_allocate_vector
  end interface map_allocate

  interface map_deallocate
     module procedure :: map_deallocate_scalar
     module procedure :: map_deallocate_vector
  end interface map_deallocate


  type sector_map8
     integer(8),dimension(:),allocatable :: map
     logical                             :: status=.false.
  end type sector_map8

  interface map_allocate8
     module procedure :: map_allocate_scalar8
     module procedure :: map_allocate_vector8
  end interface map_allocate8

  interface map_deallocate8
     module procedure :: map_deallocate_scalar8
     module procedure :: map_deallocate_vector8
  end interface map_deallocate8




  !------------------ ABTRACT INTERFACES PROCEDURES ------------------!
  !SPARSE MATRIX-VECTOR PRODUCTS USED IN ED_MATVEC
  !cmplxMat*cmplxVec
  abstract interface
     subroutine cc_sparse_HxV(Nloc,v,Hv)
       integer                    :: Nloc
       complex(8),dimension(Nloc) :: v
       complex(8),dimension(Nloc) :: Hv
     end subroutine cc_sparse_HxV
  end interface




  !-------------------------- ED  VARIABLES --------------------------!
  !SIZE OF THE PROBLEM
  !Ns       =              Number of levels per spin
  !Nlevels  = 2*Ns       = Total Number  of levels
  !Nsectors =              Number of sectors
  !=========================================================
  integer                                            :: Ns
  integer                                            :: Nlevels
  integer                                            :: Nsectors
  integer                                            :: Nhel

  !Rank-4 Tensor containing the matrix element of the interaction
  !deafult is Kanamori
  !INTERNAL USE (accessed thru functions)
  !=========================================================
  real(8),allocatable,dimension(:,:,:,:,:)           :: Umat

  !neighboring stride indexes for quartett calculations
  !INTERNAL USE (accessed thru functions)
  !=========================================================
  integer(8),allocatable,dimension(:,:,:)            :: Vstride
  integer(8),allocatable,dimension(:)                :: Neigh
  integer(8),allocatable,dimension(:,:)              :: vec2lat
  integer(8),allocatable,dimension(:,:)              :: lat2vec
  real(8),allocatable,dimension(:,:)                 :: Radius


  !local part of the Hamiltonian
  !INTERNAL USE (accessed thru functions)
  !=========================================================
  complex(8),dimension(:,:,:,:),allocatable          :: impHloc           !local hamiltonian

  !Some maps between sectors and full Hilbert space (pointers)
  !PRIVATE:
  !=========================================================
  integer,allocatable,dimension(:,:)                 :: getsector
  integer,allocatable,dimension(:,:)                 :: getsector_Jz
  integer,allocatable,dimension(:,:)                 :: getCsector
  integer,allocatable,dimension(:,:)                 :: getCDGsector
  integer,allocatable,dimension(:,:,:)               :: getCsector_Jz
  integer,allocatable,dimension(:,:,:)               :: getCDGsector_Jz
  integer,allocatable,dimension(:,:)                 :: getBathStride
  integer,allocatable,dimension(:,:)                 :: impIndex
  integer,allocatable,dimension(:)                   :: getDim,getDimUp,getDimDw
  integer,allocatable,dimension(:)                   :: getNup,getNdw
  integer,allocatable,dimension(:)                   :: getSz
  integer,allocatable,dimension(:)                   :: getN
  integer,allocatable,dimension(:)                   :: gettwoJz
  integer,allocatable,dimension(:)                   :: getmaxtwoJz
  logical,allocatable,dimension(:)                   :: twin_mask
  logical,allocatable,dimension(:)                   :: sectors_mask

  !Effective Bath used in the ED code (this is opaque to user)
  !PRIVATE
  !=========================================================
  type(effective_bath)                               :: dmft_bath


  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================
  type(sparse_matrix_csr)                            :: spH0
  procedure(cc_sparse_HxV),pointer                   :: spHtimesV_cc=>null()


  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================
  integer,allocatable,dimension(:)                   :: neigen_sector
  !--------------- LATTICE WRAP VARIABLES -----------------!
  integer,allocatable,dimension(:,:)                 :: neigen_sectorii
  integer,allocatable,dimension(:)                   :: neigen_totalii
  logical                                            :: trim_state_list=.false.

  !Partition function
  !PRIVATE
  !=========================================================
  real(8)                                            :: zeta_function



  !Impurity Green's function and Self-Energies: (Nspin,Nspin,Norb,Norb,:)
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impSmats,impSAmats
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impSreal,impSAreal
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impGmats,impFmats
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impGreal,impFreal
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impG0mats,impF0mats
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impG0real,impF0real

  !--------------- LATTICE WRAP VARIABLES -----------------!
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Smatsii,Srealii          ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: SAmatsii,SArealii        ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Gmatsii,Grealii          ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Fmatsii,Frealii          ![Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),dimension(:,:,:,:,:)  ,allocatable,save :: imp_density_matrix_ii    ![Nlat][Nspin][Nspin][Norb][Norb]
  !complex(8),dimension(:,:,:,:,:,:),allocatable,save :: bth_density_matrix_ii    ![Nlat][Nspin][Nspin][Norb][Norb][Nbath]

  !Spin Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:)                 :: spinChi_tau
  complex(8),allocatable,dimension(:,:)              :: spinChi_w
  complex(8),allocatable,dimension(:,:)              :: spinChi_iv


  !Diagonal/Off-diagonal charge-charge Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:,:)               :: densChi_tau
  complex(8),allocatable,dimension(:,:,:)            :: densChi_w
  complex(8),allocatable,dimension(:,:,:)            :: densChi_iv

  !Mixed inter-orbital charge-charge Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:,:)               :: densChi_mix_tau
  complex(8),allocatable,dimension(:,:,:)            :: densChi_mix_w
  complex(8),allocatable,dimension(:,:,:)            :: densChi_mix_iv

  !Total (orbital-sum) Density-density Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:)                   :: densChi_tot_tau
  complex(8),allocatable,dimension(:)                :: densChi_tot_w
  complex(8),allocatable,dimension(:)                :: densChi_tot_iv

  !Pair-Pair Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:)                 :: pairChi_tau
  complex(8),allocatable,dimension(:,:)              :: pairChi_w
  complex(8),allocatable,dimension(:,:)              :: pairChi_iv



  !Density and double occupancy
  !PRIVATE (now public but accessible thru routines)
  !=========================================================
  real(8),dimension(:),allocatable                   ::  ed_dens
  real(8),dimension(:),allocatable                   ::  ed_dens_up,ed_dens_dw
  real(8),dimension(:),allocatable                   ::  ed_docc
  real(8),dimension(:),allocatable                   ::  ed_phisc

  !--------------- LATTICE WRAP VARIABLES -----------------!
  real(8),dimension(:,:),allocatable,save            :: nii,dii,mii,pii


  !Local energies and generalized double occupancies
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  real(8)                                            :: ed_Ekin
  real(8)                                            :: ed_Epot
  real(8)                                            :: ed_Eint
  real(8)                                            :: ed_Ehartree
  real(8)                                            :: ed_Eknot
  real(8)                                            :: ed_Dust,ed_Dund,ed_Dse,ed_Dph
  !--------------- LATTICE WRAP VARIABLES -----------------!
  real(8),dimension(:,:),allocatable,save            :: ddii,eii


  !Impurity operators
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:)          :: imp_density_matrix
  complex(8),allocatable,dimension(:,:,:,:,:)        :: bth_density_matrix

  integer,parameter,dimension(3)                     :: Lzdiag = [-1,+1,0]
  integer,parameter,dimension(2)                     :: Szdiag = [1,-1]


  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable                   :: wm,tau,wr,vm


  !File suffixes for printing fine tuning.
  !=========================================================
  character(len=32)                                  :: ed_file_suffix=""       !suffix string attached to the output files.
  character(len=10)                                  :: ineq_site_suffix="_ineq"
  integer                                            :: site_indx_padding=4
  logical                                            :: Jhflag              !spin-exchange and pair-hopping flag.
  logical                                            :: offdiag_gf_flag=.false.



  !This is the internal Mpi Communicator and variables.
  !=========================================================
#ifdef _MPI
  integer                                            :: MpiComm_Global=MPI_UNDEFINED
  integer                                            :: MpiComm=MPI_UNDEFINED
#endif
  integer                                            :: MpiGroup_Global=MPI_GROUP_NULL
  integer                                            :: MpiGroup=MPI_GROUP_NULL
  logical                                            :: MpiStatus=.false.
  logical                                            :: MpiMaster=.true.
  integer                                            :: MpiRank=0
  integer                                            :: MpiSize=1
  integer,allocatable,dimension(:)                   :: MpiMembers
  integer                                            :: mpiQup=0
  integer                                            :: mpiRup=0
  integer                                            :: mpiQdw=0
  integer                                            :: mpiRdw=0
  integer                                            :: mpiQ=0
  integer                                            :: mpiR=0
  integer                                            :: mpiIstart
  integer                                            :: mpiIend
  integer                                            :: mpiIshift
  logical                                            :: mpiAllThreads=.true.



contains


  !=========================================================
  subroutine map_allocate_scalar(H,N)
    type(sector_map) :: H
    integer          :: N
    allocate(H%map(N))
    H%status=.true.
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N)
    type(sector_map),dimension(:) :: H
    integer,dimension(size(H))    :: N
    integer :: i
    do i=1,size(H)
       allocate(H(i)%map(N(i)))
    enddo
  end subroutine map_allocate_vector

  subroutine map_allocate_scalar8(H,N)
    type(sector_map8) :: H
    integer          :: N
    allocate(H%map(N))
    H%status=.true.
  end subroutine map_allocate_scalar8
  !
  subroutine map_allocate_vector8(H,N)
    type(sector_map8),dimension(:) :: H
    integer,dimension(size(H))    :: N
    integer :: i
    do i=1,size(H)
       allocate(H(i)%map(N(i)))
    enddo
  end subroutine map_allocate_vector8



  !=========================================================
  subroutine map_deallocate_scalar(H)
    type(sector_map) :: H
    if(.not.H%status)then
       write(*,*) "WARNING map_deallocate_scalar: H is not allocated"
       return
    endif
    if(allocated(H%map))deallocate(H%map)
    H%status=.false.
  end subroutine map_deallocate_scalar
  !
  subroutine map_deallocate_vector(H)
    type(sector_map),dimension(:) :: H
    integer :: i
    do i=1,size(H)
       deallocate(H(i)%map)
    enddo
  end subroutine map_deallocate_vector

  subroutine map_deallocate_scalar8(H)
    type(sector_map8) :: H
    if(.not.H%status)then
       write(*,*) "WARNING map_deallocate_scalar: H is not allocated"
       return
    endif
    if(allocated(H%map))deallocate(H%map)
    H%status=.false.
  end subroutine map_deallocate_scalar8
  !
  subroutine map_deallocate_vector8(H)
    type(sector_map8),dimension(:) :: H
    integer :: i
    do i=1,size(H)
       deallocate(H(i)%map)
    enddo
  end subroutine map_deallocate_vector8


  !=========================================================
  subroutine ed_set_MpiComm(comm)
#ifdef _MPI
    integer :: comm,ierr
    MpiComm_Global = comm
    MpiComm        = MpiComm_Global
    MpiStatus      = .true.
    MpiSize        = get_Size_MPI(MpiComm_Global)
    MpiRank        = get_Rank_MPI(MpiComm_Global)
    MpiMaster      = get_Master_MPI(MpiComm_Global)
    call Mpi_Comm_group(MpiComm_Global,MpiGroup_Global,ierr)
#else
    integer,optional :: comm
#endif
  end subroutine ed_set_MpiComm

  subroutine ed_del_MpiComm()
#ifdef _MPI
    MpiComm_Global = MPI_UNDEFINED
    MpiComm        = MPI_UNDEFINED
    MpiStatus      = .false.
    MpiSize        = 1
    MpiRank        = 0
    MpiMaster      = .true.
#endif
  end subroutine ed_del_MpiComm


END MODULE ED_VARS_GLOBAL
