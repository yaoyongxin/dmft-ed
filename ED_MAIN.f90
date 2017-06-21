module ED_MAIN
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE, only: state_list,es_delete_espace
  USE ED_AUX_FUNX
  USE ED_IO,only: ed_print_impSigma
  USE ED_SETUP
  USE ED_BATH
  USE ED_HAMILTONIAN_MATVEC
  USE ED_GREENS_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_DIAG
  USE SF_LINALG
  USE SF_IOTOOLS, only: str
  USE SF_TIMER,only: start_timer,stop_timer
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private

  !
  !>INIT ED SOLVER
  !
  interface ed_init_solver
     module procedure :: ed_init_solver_single
     module procedure :: ed_init_solver_lattice
#ifdef _MPI
     module procedure :: ed_init_solver_single_mpi
     module procedure :: ed_init_solver_lattice_mpi
#endif
  end interface ed_init_solver
  !>
  public :: ed_init_solver


  !
  !> ED SOLVER
  !
  interface ed_solve
     module procedure :: ed_solve_single
     module procedure :: ed_solve_lattice
#ifdef _MPI
     module procedure :: ed_solve_single_mpi
     module procedure :: ed_solve_lattice_mpi
#endif
  end interface ed_solve
  !>
  public :: ed_solve



  real(8),dimension(:),allocatable                   :: wr,wm
  character(len=64)                                  :: suffix



contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_single(bath,Hloc)
    real(8),dimension(:),intent(inout) :: bath
    complex(8),intent(in)              :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                            :: check 
    logical,save                       :: isetup=.true.
    integer                            :: i
    logical                            :: MPI_MASTER=.true.
    integer                            :: MPI_ERR
    !
    if(MPI_MASTER)write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure()
    !
    !Init bath:
    call set_Hloc(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "init_ed_solver_single error: wrong bath dimensions"
    !
    bath = 0d0
    !
    call allocate_dmft_bath(dmft_bath)
    if(bath_type=="replica")call init_dmft_bath_mask(dmft_bath)
    call init_dmft_bath(dmft_bath)
    call get_dmft_bath(dmft_bath,bath)
    !
    if(isetup)then
       select case(ed_mode)
       case default
          call setup_pointers_normal
       case ("superc")
          call setup_pointers_superc
       case ("nonsu2")
          call setup_pointers_nonsu2
       end select
    endif
    call deallocate_dmft_bath(dmft_bath)
    isetup=.false.
    !
  end subroutine ed_init_solver_single

#ifdef _MPI
  subroutine ed_init_solver_single_mpi(MpiComm,bath,Hloc)
    integer                            :: MpiComm
    real(8),dimension(:),intent(inout) :: bath
    complex(8),intent(in)              :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                            :: check 
    logical,save                       :: isetup=.true.
    integer                            :: i
    logical                            :: MPI_MASTER=.true.
    integer                            :: MPI_ERR
    !
    MPI_MASTER = get_Master_MPI(MpiComm)
    !
    if(MPI_MASTER)write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure(MpiComm)
    !
    !Init bath:
    call set_hloc(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "init_ed_solver_single error: wrong bath dimensions"
    !
    bath = 0d0
    !
    call allocate_dmft_bath(dmft_bath)
    if(bath_type=="replica")call init_dmft_bath_mask(dmft_bath)
    call init_dmft_bath(dmft_bath)
    call get_dmft_bath(dmft_bath,bath)
    !
    if(isetup)then
       select case(ed_mode)
       case default
          call setup_pointers_normal
       case ("superc")
          call setup_pointers_superc
       case ("nonsu2")
          call setup_pointers_nonsu2
       end select
    endif
    call deallocate_dmft_bath(dmft_bath)
    isetup=.false.
    !
    call MPI_Barrier(MpiComm,MPI_ERR)
    !
  end subroutine ed_init_solver_single_mpi
#endif








  !+-----------------------------------------------------------------------------+!
  !                           INEQUVALENT SITES                                   !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_lattice(bath,Hloc)
    real(8),dimension(:,:)         :: bath ![Nlat][:]
    complex(8),intent(in)          :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer                        :: ilat,Nineq,Nsect
    logical                        :: check_dim
    character(len=5)               :: tmp_suffix
    integer                        :: MPI_ERR
    !
    !
    Nineq = size(bath,1)
    Nsect = get_Nsectors() !< get # sectors to allocate the following array
    if(allocated(neigen_sectorii))deallocate(neigen_sectorii) ; allocate(neigen_sectorii(Nineq,Nsect))
    if(allocated(neigen_totalii))deallocate(neigen_totalii) ; allocate(neigen_totalii(Nineq))
    do ilat=1,Nineq             !all nodes check the bath, u never know...
       !
       ! if(present(Hloc))then
       !    check_dim = check_bath_dimension(bath(ilat,:),Hloc(ilat,:,:,:,:))
       ! else
       !    check_dim = check_bath_dimension(bath(ilat,:))
       ! endif
       ! if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
       !
       ed_file_suffix="_site"//str(ilat,Npad=4)
       !
       ! if(present(Hloc))then
       call ed_init_solver_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
       ! else
       !    call ed_init_solver_single(bath(ilat,:))
       ! endif
       neigen_sectorii(ilat,:) = neigen_sector(:)
       neigen_totalii(ilat)    = lanc_nstates_total
    end do
    !
    ed_file_suffix=""
    !
  end subroutine ed_init_solver_lattice

#ifdef _MPI
  subroutine ed_init_solver_lattice_mpi(MpiComm,bath,Hloc)
    integer                        :: MpiComm
    real(8),dimension(:,:)         :: bath ![Nlat][:]
    complex(8),intent(in)          :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer                        :: ilat,Nineq,Nsect
    logical                        :: check_dim
    character(len=5)               :: tmp_suffix
    integer                        :: MPI_ERR
    !
    !
    Nineq = size(bath,1)
    Nsect = get_Nsectors() !< get # sectors to allocate the following array
    if(allocated(neigen_sectorii))deallocate(neigen_sectorii) ; allocate(neigen_sectorii(Nineq,Nsect))
    if(allocated(neigen_totalii))deallocate(neigen_totalii) ; allocate(neigen_totalii(Nineq))
    do ilat=1,Nineq             !all nodes check the bath, u never know...
       ! if(present(Hloc))then
       !    check_dim = check_bath_dimension(bath(ilat,:),Hloc(ilat,:,:,:,:))
       ! else
       !    check_dim = check_bath_dimension(bath(ilat,:))
       ! endif
       ! if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
       !
       ed_file_suffix="_site"//str(ilat,Npad=4)
       !
       ! if(present(Hloc))then
       call ed_init_solver_single_mpi(MpiComm,bath(ilat,:),Hloc(ilat,:,:,:,:))
       ! else
       !    call ed_init_solver_single_mpi(MpiComm,bath(ilat,:))
       ! endif
       neigen_sectorii(ilat,:) = neigen_sector(:)
       neigen_totalii(ilat)    = lanc_nstates_total
    end do
    !
    call MPI_Barrier(MpiComm,MPI_ERR)
    !
    ed_file_suffix=""
    !
  end subroutine ed_init_solver_lattice_mpi
#endif
















  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+-----------------------------------------------------------------------------+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_single(bath,Hloc)
    integer                         :: MpiComm
    real(8),dimension(:),intent(in) :: bath
    complex(8),optional,intent(in)  :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                         :: check
    logical                         :: MPI_MASTER=.true.
    !
    if(present(Hloc))call set_Hloc(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "ED_SOLVE_SINGLE Error: wrong bath dimensions"
    !
    call allocate_dmft_bath(dmft_bath)
    if(bath_type=="replica")call init_dmft_bath_mask(dmft_bath)
    call set_dmft_bath(bath,dmft_bath)
    if(MPI_MASTER)then
       call write_dmft_bath(dmft_bath,LOGfile)
       call save_dmft_bath(dmft_bath,used=.true.)
    endif
    !
    !ASSOCIATE THE GLOBAL PROCEDURES
    ed_buildh_c  => build_H_c
    !
    select case(ed_sparse_H)
    case (.true.)
       spHtimesV_cc => spMatVec_cc
    case (.false.)
       spHtimesV_cc => directMatVec_cc
    case default
       stop "ed_solve_single ERROR: ed_sparse_H undefined"
    end select
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()         !find target states by digonalization of Hamiltonian
    call observables_impurity()         !obtain impurity observables as thermal averages.  
    call buildgf_impurity()             !build the one-particle impurity Green's functions  & Self-energy
    if(chiflag)call buildchi_impurity() !build the local susceptibilities (spin [todo charge])
    call local_energy_impurity()        !obtain the local energy of the effective impurity problem.
    !
    call deallocate_dmft_bath(dmft_bath)   
    call es_delete_espace(state_list)
    !
    nullify(ed_buildh_c)
    nullify(spHtimesV_cc)
  end subroutine ed_solve_single

#ifdef _MPI
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_single_mpi(MpiComm,bath,Hloc)
    integer                         :: MpiComm
    real(8),dimension(:),intent(in) :: bath
    complex(8),optional,intent(in)  :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                         :: check
    logical                         :: MPI_MASTER=.true.
    !
    MPI_MASTER = get_Master_MPI(MpiComm)
    !
    if(present(Hloc))call set_Hloc(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "ED_SOLVE_SINGLE Error: wrong bath dimensions"
    !
    call allocate_dmft_bath(dmft_bath)
    if(bath_type=="replica")call init_dmft_bath_mask(dmft_bath)
    call set_dmft_bath(bath,dmft_bath)
    if(MPI_MASTER)then
       call write_dmft_bath(dmft_bath,LOGfile)
       call save_dmft_bath(dmft_bath,used=.true.)
    endif
    !
    !ASSOCIATE THE GLOBAL PROCEDURES
    ed_buildh_c  => build_H_c
    !
    select case(ed_sparse_H)
    case (.true.)
       spHtimesV_cc => spMatVec_MPI_cc
    case (.false.)
       spHtimesV_cc => directMatVec_MPI_cc
    case default
       stop "ed_solve_single_mpi ERROR: ed_sparse_H undefined"
    end select
    !
    !SET THE LOCAL COMMUNICATORS IN ALL THE RELEVANT PARTS OF THE CODE:
    call ed_hamiltonian_matvec_set_MPI(MpiComm)
    call ed_diag_set_MPI(MpiComm)
    call ed_observables_set_MPI(MpiComm)
    call ed_greens_functions_set_MPI(MpiComm)
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()         !find target states by digonalization of Hamiltonian
    call observables_impurity()         !obtain impurity observables as thermal averages.  
    call buildgf_impurity()             !build the one-particle impurity Green's functions  & Self-energy
    if(chiflag)call buildchi_impurity() !build the local susceptibilities (spin [todo charge])
    call local_energy_impurity()        !obtain the local energy of the effective impurity problem.
    !
    call deallocate_dmft_bath(dmft_bath)   
    call es_delete_espace(state_list)
    !
    !DELETE THE LOCAL COMMUNICATORS IN ALL THE RELEVANT PARTS OF THE CODE:
    call ed_hamiltonian_matvec_del_MPI()
    call ed_diag_del_MPI()
    call ed_observables_del_MPI()
    call ed_greens_functions_del_MPI()    
    nullify(ed_buildh_c)
    nullify(spHtimesV_cc)
  end subroutine ed_solve_single_mpi
#endif







  





  !+-----------------------------------------------------------------------------+!
  !                          INEQUIVALENT SITES                                   !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_lattice(bath,Hloc,Uloc_ii,Ust_ii,Jh_ii)
    !inputs
    real(8)          :: bath(:,:) ![Nlat][Nb]
    complex(8)       :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    real(8),optional :: Uloc_ii(size(bath,1),Norb)
    real(8),optional :: Ust_ii(size(bath,1))
    real(8),optional :: Jh_ii(size(bath,1))
    ! 
    integer          :: i,j,ilat,iorb,jorb,ispin,jspin
    integer          :: Nsites
    logical          :: check_dim
    character(len=5) :: tmp_suffix
    !
    logical          :: MPI_MASTER=.true.
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    !Allocate the local static observarbles global to the module
    !One can retrieve these values from suitable routines later on
    if(allocated(nii))deallocate(nii)
    if(allocated(dii))deallocate(dii)
    if(allocated(mii))deallocate(mii)
    if(allocated(pii))deallocate(pii)
    if(allocated(eii))deallocate(eii)
    if(allocated(ddii))deallocate(ddii)
    allocate(nii(Nsites,Norb))
    allocate(dii(Nsites,Norb))
    allocate(mii(Nsites,Norb))
    allocate(pii(Nsites,Norb))
    allocate(eii(Nsites,4))
    allocate(ddii(Nsites,4))
    !
    !Allocate the self-energies global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Smatsii))deallocate(Smatsii)
    if(allocated(Srealii))deallocate(Srealii)
    if(allocated(SAmatsii))deallocate(SAmatsii)
    if(allocated(SArealii))deallocate(SArealii)
    allocate(Smatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Srealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SArealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    !Allocate the imp GF global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Gmatsii))deallocate(Gmatsii)
    if(allocated(Grealii))deallocate(Grealii)
    if(allocated(Fmatsii))deallocate(Fmatsii)
    if(allocated(Frealii))deallocate(Frealii)
    allocate(Gmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Grealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Fmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Frealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    if(size(neigen_sectorii,1)<Nsites)stop "ed_solve_lattice error: size(neigen_sectorii,1)<Nsites"
    if(size(neigen_totalii)<Nsites)stop "ed_solve_lattice error: size(neigen_totalii,1)<Nsites"
    !
    !Check the dimensions of the bath are ok:
    do ilat=1,Nsites
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smatsii  = zero 
    Srealii  = zero 
    SAmatsii = zero 
    SArealii = zero 
    Gmatsii  = zero 
    Grealii  = zero 
    Fmatsii  = zero 
    Frealii  = zero 
    nii      = 0d0  
    dii      = 0d0  
    mii      = 0d0  
    pii      = 0d0  
    eii      = 0d0  
    ddii     = 0d0  
    !
    call start_timer
    !
    do ilat = 1, Nsites
       write(*,*)" solves site: "//str(ilat,Npad=4)
       !
       ed_file_suffix="_site"//str(ilat,Npad=4)
       !
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust          = Ust_ii(ilat) 
       if(present(Jh_ii))  Jh           = Jh_ii(ilat) 
       !
       !Solve the impurity problem for the ilat-th site
       neigen_sector(:)   = neigen_sectorii(ilat,:)
       lanc_nstates_total = neigen_totalii(ilat)
       !
       !Call ed_solve in SERIAL MODE!! This is parallel on the ineq. sites
       call ed_solve_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
       !
       neigen_sectorii(ilat,:)   = neigen_sector(:)
       neigen_totalii(ilat)      = lanc_nstates_total
       Smatsii(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
       Srealii(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
       SAmatsii(ilat,:,:,:,:,:) = impSAmats(:,:,:,:,:)
       SArealii(ilat,:,:,:,:,:) = impSAreal(:,:,:,:,:)
       Gmatsii(ilat,:,:,:,:,:)  = impGmats(:,:,:,:,:)
       Grealii(ilat,:,:,:,:,:)  = impGreal(:,:,:,:,:)
       Fmatsii(ilat,:,:,:,:,:)  = impFmats(:,:,:,:,:)
       Frealii(ilat,:,:,:,:,:)  = impFreal(:,:,:,:,:)
       nii(ilat,1:Norb)       = ed_dens(1:Norb)
       dii(ilat,1:Norb)       = ed_docc(1:Norb)
       mii(ilat,1:Norb)       = ed_dens_up(1:Norb)-ed_dens_dw(1:Norb)
       pii(ilat,1:Norb)       = ed_phisc(1:Norb)
       eii(ilat,:)            = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
       ddii(ilat,:)           = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
    enddo
    !
    call stop_timer
    !
    ed_file_suffix=""
    !
  end subroutine ed_solve_lattice


  !FALL BACK: DO A VERSION THAT DOES THE SITES IN PARALLEL USING SERIAL ED CODE
#ifdef _MPI
  subroutine ed_solve_lattice_mpi(MpiComm,bath,Hloc,Uloc_ii,Ust_ii,Jh_ii)
    integer          :: MpiComm
    !inputs
    real(8)          :: bath(:,:) ![Nlat][Nb]
    complex(8)       :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    real(8),optional :: Uloc_ii(size(bath,1),Norb)
    real(8),optional :: Ust_ii(size(bath,1))
    real(8),optional :: Jh_ii(size(bath,1))
    !MPI  auxiliary vars
    complex(8)       :: Smats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Sreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: SAmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: SAreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: Gmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Greal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: Fmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Freal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    real(8)          :: nii_tmp(size(bath,1),Norb)
    real(8)          :: dii_tmp(size(bath,1),Norb)
    real(8)          :: mii_tmp(size(bath,1),Norb)
    real(8)          :: pii_tmp(size(bath,1),Norb)
    real(8)          :: eii_tmp(size(bath,1),4)
    real(8)          :: ddii_tmp(size(bath,1),4)
    !
    integer          :: neigen_sectortmp(size(bath,1),Nsectors)
    integer          :: neigen_totaltmp(size(bath,1))
    ! 
    integer          :: i,j,ilat,iorb,jorb,ispin,jspin
    integer          :: Nsites
    logical          :: check_dim
    character(len=5) :: tmp_suffix
    !
    integer          :: MPI_ID=0
    integer          :: MPI_SIZE=1
    logical          :: MPI_MASTER=.true.
    !
    integer          :: mpi_err 
    !
    MPI_ID     = get_Rank_MPI(MpiComm)
    MPI_SIZE   = get_Size_MPI(MpiComm)
    MPI_MASTER = get_Master_MPI(MpiComm)
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    !Allocate the local static observarbles global to the module
    !One can retrieve these values from suitable routines later on
    if(allocated(nii))deallocate(nii)
    if(allocated(dii))deallocate(dii)
    if(allocated(mii))deallocate(mii)
    if(allocated(pii))deallocate(pii)
    if(allocated(eii))deallocate(eii)
    if(allocated(ddii))deallocate(ddii)
    allocate(nii(Nsites,Norb))
    allocate(dii(Nsites,Norb))
    allocate(mii(Nsites,Norb))
    allocate(pii(Nsites,Norb))
    allocate(eii(Nsites,4))
    allocate(ddii(Nsites,4))
    !
    !Allocate the self-energies global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Smatsii))deallocate(Smatsii)
    if(allocated(Srealii))deallocate(Srealii)
    if(allocated(SAmatsii))deallocate(SAmatsii)
    if(allocated(SArealii))deallocate(SArealii)
    allocate(Smatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Srealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(SAmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(SArealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    !Allocate the imp GF global to the module
    !Once can retrieve these functinos from suitable routines later on
    if(allocated(Gmatsii))deallocate(Gmatsii)
    if(allocated(Grealii))deallocate(Grealii)
    if(allocated(Fmatsii))deallocate(Fmatsii)
    if(allocated(Frealii))deallocate(Frealii)
    allocate(Gmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Grealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Fmatsii(Nsites,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Frealii(Nsites,Nspin,Nspin,Norb,Norb,Lreal))
    !
    if(size(neigen_sectorii,1)<Nsites)stop "ed_solve_lattice error: size(neigen_sectorii,1)<Nsites"
    if(size(neigen_totalii)<Nsites)stop "ed_solve_lattice error: size(neigen_totalii,1)<Nsites"
    neigen_sectortmp = 0
    neigen_totaltmp  = 0
    !
    !Check the dimensions of the bath are ok:
    do ilat=1+MPI_ID,Nsites,MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smats_tmp  = zero
    Sreal_tmp  = zero
    SAmats_tmp = zero
    SAreal_tmp = zero
    Gmats_tmp  = zero
    Greal_tmp  = zero
    Fmats_tmp  = zero
    Freal_tmp  = zero
    nii_tmp    = 0d0
    dii_tmp    = 0d0
    mii_tmp    = 0d0
    pii_tmp    = 0d0
    eii_tmp    = 0d0
    ddii_tmp   = 0d0
    !
    if(MPI_MASTER)call start_timer
    if(.not.MPI_MASTER)LOGfile = 800+MPI_ID
    !
    do ilat = 1 + MPI_ID, Nsites, MPI_SIZE
       write(*,*)str(MPI_ID)//" solves site: "//str(ilat,Npad=4)
       !
       ed_file_suffix="_site"//str(ilat,Npad=4)
       !
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust          = Ust_ii(ilat) 
       if(present(Jh_ii))  Jh           = Jh_ii(ilat) 
       !
       !Solve the impurity problem for the ilat-th site
       neigen_sector(:)   = neigen_sectorii(ilat,:)
       lanc_nstates_total = neigen_totalii(ilat)
       !
       !Call ed_solve in SERIAL MODE!! This is parallel on the ineq. sites
       call ed_solve_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
       !
       neigen_sectortmp(ilat,:)   = neigen_sector(:)
       neigen_totaltmp(ilat)      = lanc_nstates_total
       Smats_tmp(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
       Sreal_tmp(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
       SAmats_tmp(ilat,:,:,:,:,:) = impSAmats(:,:,:,:,:)
       SAreal_tmp(ilat,:,:,:,:,:) = impSAreal(:,:,:,:,:)
       Gmats_tmp(ilat,:,:,:,:,:)  = impGmats(:,:,:,:,:)
       Greal_tmp(ilat,:,:,:,:,:)  = impGreal(:,:,:,:,:)
       Fmats_tmp(ilat,:,:,:,:,:)  = impFmats(:,:,:,:,:)
       Freal_tmp(ilat,:,:,:,:,:)  = impFreal(:,:,:,:,:)
       nii_tmp(ilat,1:Norb)       = ed_dens(1:Norb)
       dii_tmp(ilat,1:Norb)       = ed_docc(1:Norb)
       mii_tmp(ilat,1:Norb)       = ed_dens_up(1:Norb)-ed_dens_dw(1:Norb)
       pii_tmp(ilat,1:Norb)       = ed_phisc(1:Norb)
       eii_tmp(ilat,:)            = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
       ddii_tmp(ilat,:)           = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
    enddo
    !
    call MPI_Barrier(MpiComm,MPI_ERR)
    !
    if(MPI_MASTER)call stop_timer
    !
    ed_file_suffix=""
    !
    neigen_sectorii=0
    neigen_totalii =0
    !
    Smatsii  = zero 
    Srealii  = zero 
    SAmatsii = zero 
    SArealii = zero 
    Gmatsii  = zero 
    Grealii  = zero 
    Fmatsii  = zero 
    Frealii  = zero 
    nii      = 0d0  
    dii      = 0d0  
    mii      = 0d0  
    pii      = 0d0  
    eii      = 0d0  
    ddii     = 0d0  
    call MPI_ALLREDUCE(neigen_sectortmp,neigen_sectorii,Nsites*Nsectors,MPI_INTEGER,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(neigen_totaltmp,neigen_totalii,Nsites,MPI_INTEGER,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(Smats_tmp,Smatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(Sreal_tmp,Srealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(SAmats_tmp,SAmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(SAreal_tmp,SArealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(Gmats_tmp,Gmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(Greal_tmp,Grealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(Fmats_tmp,Fmatsii,Nsites*Nspin*Nspin*Norb*Norb*Lmats,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(Freal_tmp,Frealii,Nsites*Nspin*Nspin*Norb*Norb*Lreal,MPI_DOUBLE_COMPLEX,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(nii_tmp,nii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(dii_tmp,dii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(mii_tmp,mii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(pii_tmp,pii,Nsites*Norb,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(eii_tmp,eii,Nsites*4,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_err)
    call MPI_ALLREDUCE(ddii_tmp,ddii,Nsites*4,MPI_DOUBLE_PRECISION,MPI_SUM,MpiComm,mpi_err)
    !
  end subroutine ed_solve_lattice_mpi
#endif




end module ED_MAIN
