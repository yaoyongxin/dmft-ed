MODULE ED_SETUP
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,file_length
  USE DMFT_INTERACTION
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private

  public :: init_ed_structure
  !
  public :: setup_pointers_normal
  public :: setup_pointers_superc
  public :: setup_pointers_nonsu2
  !
  public :: build_sector
  public :: build_sector_2
  public :: delete_sector
  public :: delete_sector8
  !
  public :: bdecomp, bdecomp8
  public :: bjoin
  !
  public :: c,cdg
  public :: c8,cdg8
  !
  public :: binary_search
  public :: binary_search8
  !
  public :: flip_state
  public :: twin_sector_order
  public :: get_twin_sector
  !
#ifdef _MPI
  public :: scatter_vector_MPI
  public :: scatter_basis_MPI
  public :: gather_vector_MPI
  public :: allgather_vector_MPI
#endif

  ! public :: get_Nsectors
  ! !



contains



  subroutine ed_checks_global
    !
    !CHECKS:
    if(Lfit>Lmats)Lfit=Lmats
    if(Nspin>2)stop "ED ERROR: Nspin > 2 is currently not supported"
    if((.not.plaquette) .and. Norb>3)stop "ED ERROR: Norb > 3 is currently not supported"
    !
    if(ed_mode=="superc")then
       if(Nspin>1)stop "ED ERROR: SC + AFM is currently not supported ."
    endif
    if(ed_mode=="nonsu2")then
       if(Nspin/=2)then
          write(LOGfile,"(A)")"ED msg: ed_mode=nonSU2 with Nspin!=2 is not allowed."
          write(LOGfile,"(A)")"        to enforce spin symmetry up-dw set ed_para=T."
          stop
       endif
    endif
    !
    if(Nspin>1.AND.ed_twin.eqv..true.)then
       write(LOGfile,"(A)")"WARNING: using twin_sector with Nspin>1"
       call sleep(1)
    end if
    !
    if(lanc_method=="lanczos")then
       if(lanc_nstates_total>1)stop "ED ERROR: lanc_method==lanczos available only for lanc_nstates_total==1, T=0"
       if(lanc_nstates_sector>1)stop "ED ERROR: lanc_method==lanczos available only for lanc_nstates_sector==1, T=0"
    endif
    !
    !if(ed_sectors.AND.ed_mode/="normal")then
    !   stop "ED_ERROR: using ed_sectors with ed_mode=[superc,nonsu2] NOT TESTED! Uncomment this line in ED_SETUP if u want to take the risk.."
    !endif
  end subroutine ed_checks_global






  !+------------------------------------------------------------------+
  !PURPOSE  : Setup Dimensions of the problem
  ! Norb    = # of impurity orbitals
  ! Nbath   = # of bath levels (depending on bath_type)
  ! Ns      = # of levels (per spin)
  ! Nlevels = 2*Ns = Total # of levels (counting spin degeneracy 2)
  !+------------------------------------------------------------------+
  subroutine ed_setup_dimensions()
    integer                                           :: maxtwoJz,inJz,dimJz
    integer                                           :: isector,in,shift
    select case(bath_type)
    case default
       Ns = (Nbath+1)*Norb
    case ('hybrid')
       Ns = Nbath+Norb
    case ('replica')
       Ns = Norb*(Nbath+1)
    end select
    Nlevels  = 2*Ns
    !
    if(plaquette)then
       Ns=Norb*Nbath
       Nlevels=Ns
    endif
    !
    select case(ed_mode)
    case default
       Nsectors = (Ns+1)*(Ns+1) !nup=0:Ns;ndw=0:Ns
       Nhel     = 1
       if(plaquette)Nsectors=Ns+1
    case ("superc")
       Nsectors = Nlevels+1     !sz=-Ns:Ns=2*Ns+1=Nlevels+1
       Nhel     = 1
    case("nonsu2")
       if(Jz_basis)then
          isector=0
          do in=0,Nlevels
             !algorithm to find the maximum Jz given the density
             if(in==0.or.in==2*Ns)then
                maxtwoJz=0
             else
                shift=0
                if(in<=Nbath+1)shift=Nbath-in+1
                if(in>=2*Ns-Nbath)shift=Nbath-2*Ns+in+1
                maxtwoJz = 5 + 5*Nbath - abs(in-Ns) - 2*shift
             endif
             !number of available Jz given the maximum value
             dimJz = maxtwoJz + 1
             !count of all the new Jz sectors
             do inJz=1,dimJz
                isector=isector+1
             enddo
          enddo
          Nsectors=isector
          Nhel     = 2
       else
          Nsectors = Nlevels+1     !n=0:2*Ns=2*Ns+1=Nlevels+1
          Nhel     = 2
       endif
    end select
  end subroutine ed_setup_dimensions



  !+------------------------------------------------------------------+
  !PURPOSE  : Init ED structure and calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure(MpiComm)
    integer,optional                                  :: MpiComm
    logical                                           :: control
    real(8),dimension(Nspin,Nspin,Norb,Norb)          :: reHloc         !local hamiltonian, real part
    real(8),dimension(Nspin,Nspin,Norb,Norb)          :: imHloc         !local hamiltonian, imag part
    integer                                           :: i,dim_sector_max(2),iorb,jorb,ispin,jspin
    logical                                           :: MPI_MASTER=.true.
    !
#ifdef _MPI
    if(present(MpiComm))MPI_MASTER=get_Master_MPI(MpiComm)
#endif
    !
    call ed_checks_global
    !
    call ed_setup_dimensions
    !
    dim_sector_max=0
    select case(ed_mode)
    case default
       dim_sector_max(1)=get_normal_sector_dimension(Ns/2)
       dim_sector_max(2)=get_normal_sector_dimension(Ns-Ns/2)
    case ("superc")
       dim_sector_max(1)=get_superc_sector_dimension(0)
    case("nonsu2")
       dim_sector_max(1)=get_nonsu2_sector_dimension(Ns)
       if(plaquette)dim_sector_max(1)=get_normal_sector_dimension(Ns/2)
    end select
    !
    write(LOGfile,"(A)")"Summary:"
    write(LOGfile,"(A)")"--------------------------------------------"
    write(LOGfile,"(A,I15)")'# of levels/spin      = ',Ns
    write(LOGfile,"(A,I15)")'Total size            = ',Nlevels
    write(LOGfile,"(A,I15)")'# of impurities       = ',Norb
    write(LOGfile,"(A,I15)")'# of bath/impurity    = ',Nbath
    write(LOGfile,"(A,I15)")'# of Bath levels/spin = ',Ns-Norb
    write(LOGfile,"(A,2I15)")'Largest Sector(s)    = ',dim_sector_max
    write(LOGfile,"(A,I15)")'Number of sectors     = ',Nsectors
    write(LOGfile,"(A)")"--------------------------------------------"
    call sleep(1)
    !
    !
    !Allocate indexing impHloc
    allocate(impHloc(Nspin,Nspin,Norb,Norb))
    impHloc = zero
    reHloc = 0d0 ; imHloc = 0d0
    !
    !Search and read impHloc
    inquire(file=trim(HLOCfile),exist=control)
    if(control)then
       write(LOGfile,*)"Reading impHloc from file: "//reg(HLOCfile)
       open(50,file=trim(HLOCfile),status='old')
       do ispin=1,Nspin
          do iorb=1,Norb
             read(50,*)((reHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       do ispin=1,Nspin
          do iorb=1,Norb
             read(50,*)((imHloc(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       close(50)
    else
       write(LOGfile,*)"impHloc file not found."
       write(LOGfile,*)"impHloc should be defined elsewhere..."
       call sleep(2)
    endif
    impHloc = dcmplx(reHloc,imHloc)
    if(control)then
       write(LOGfile,"(A)")"H_local:"
       call print_Hloc(impHloc)
    endif
    !
    !
    !Search and read interaction matrix
    if(Utensor)then
       !
       write(LOGfile,*)
       write(LOGfile,*)"Tensor representation fo the interaction asked. Beware, HFMODE is uneffective."
       allocate(Umat(Norb,Norb,Norb,Norb,2))
       Umat = 0d0
       !
       inquire(file=trim(UTENSfile),exist=control)
       if(control)then
          call dmft_interaction_read(Umat,trim(UTENSfile))
       else
          write(LOGfile,*)"Interaction matrix file not found, setting to deafault Kanamori."
          if(present(MpiComm))then
             call dmft_interaction_setKanamori(MpiComm,Umat,Uloc,Ust,Jh,Jx,Jp)
          else
             call dmft_interaction_setKanamori(Umat,Uloc,Ust,Jh,Jx,Jp)
          endif
          call sleep(1)
       endif
       !
       if(MPI_MASTER)call dmft_interaction_print(Umat,"Utensor.in")
       !
    endif
    !
    !
    !Setup the stride for plaquette calculations
    if(plaquette)then
       if(.not.allocated(Vstride))stop "Vstride not allocated."
    endif
    !
    !
    !Allocate indexing arrays
    allocate(impIndex(Norb,2));impIndex=0
    allocate(getDim(Nsectors));getDim=0
    !
    allocate(getDimUp(Nsectors),getDimDw(Nsectors));getDimUp=0;getDimDw=0
    allocate(getNup(Nsectors),getNdw(Nsectors));getNup=0;getNdw=0
    !
    allocate(getSz(Nsectors));getSz=0
    !
    allocate(getN(Nsectors));getN=0
    !
    allocate(gettwoJz(Nsectors));gettwoJz=0
    allocate(getmaxtwoJz(0:Nlevels));getmaxtwoJz=0
    !
    select case(ed_mode)
    case default
       allocate(getSector(0:Ns,0:Ns))
    case ("superc")
       allocate(getSector(-Ns:Ns,1))
    case ("nonsu2")
       if(Jz_basis)then
          allocate(getSector(0:Nlevels,-Nlevels:Nlevels));getSector=0
       else
          allocate(getSector(0:Nlevels,1));getSector=0
       endif
    end select
    getSector=0
    !
    allocate(getCsector(2,Nsectors));getCsector=0
    allocate(getCDGsector(2,Nsectors));getCDGsector=0
    !
    allocate(getCsector_Jz(Norb,Nspin,Nsectors));  getCsector_Jz=-1
    allocate(getCDGsector_Jz(Norb,Nspin,Nsectors));getCDGsector_Jz=-1
    !
    allocate(getBathStride(Norb,Nbath));getBathStride=0
    allocate(twin_mask(Nsectors));
    allocate(sectors_mask(Nsectors))
    allocate(neigen_sector(Nsectors))
    !
    !
    !check finiteT
    finiteT=.true.              !assume doing finite T per default
    if(lanc_nstates_total==1)then     !is you only want to keep 1 state
       finiteT=.false.          !set to do zero temperature calculations
       write(LOGfile,"(A)")"Required Lanc_nstates_total=1 => set T=0 calculation"
    endif
    !
    !
    !check whether lanc_nstates_sector and lanc_states are even (we do want to keep doublet among states)
    if(finiteT)then
       if(mod(lanc_nstates_sector,2)/=0)then
          lanc_nstates_sector=lanc_nstates_sector+1
          write(LOGfile,"(A,I10)")"Increased Lanc_nstates_sector:",lanc_nstates_sector
       endif
       if(mod(lanc_nstates_total,2)/=0)then
          lanc_nstates_total=lanc_nstates_total+1
          write(LOGfile,"(A,I10)")"Increased Lanc_nstates_total:",lanc_nstates_total
       endif
       write(LOGfile,"(A)")"Lanczos FINITE temperature calculation:"
       !
       write(LOGfile,"(A,I3)")"Nstates x Sector = ", lanc_nstates_sector
       write(LOGfile,"(A,I3)")"Nstates   Total  = ", lanc_nstates_total
       call sleep(1)
    else
       write(LOGfile,"(A)")"Lanczos ZERO temperature calculation:"
       call sleep(1)
    endif
    !
    Jhflag=.FALSE.
    if(Norb>1.AND.(Jx/=0d0.OR.Jp/=0d0))Jhflag=.TRUE.
    !
    !
    offdiag_gf_flag=ed_solve_offdiag_gf
    if(bath_type/="normal")offdiag_gf_flag=.true.
    !
    !
    if(nread/=0.d0)then
       i=abs(floor(log10(abs(nerr)))) !modulus of the order of magnitude of nerror
       niter=nloop/3
    endif
    !
    !
    !allocate functions
    allocate(impSmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impSAmats(Nspin,Nspin,Norb,Norb,Lmats)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    allocate(impSAreal(Nspin,Nspin,Norb,Norb,Lreal)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    impSmats=zero
    impSreal=zero
    impSAmats=zero
    impSAreal=zero
    !
    allocate(impGmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impGreal(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impFmats(Nspin,Nspin,Norb,Norb,Lmats)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    allocate(impFreal(Nspin,Nspin,Norb,Norb,Lreal)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    impGmats=zero
    impGreal=zero
    impFmats=zero
    impFreal=zero
    !
    allocate(impG0mats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impG0real(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(impF0mats(Nspin,Nspin,Norb,Norb,Lmats)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    allocate(impF0real(Nspin,Nspin,Norb,Norb,Lreal)) !THIS SHOULD NOT DEPEND ON SPIN: NSPIN=>1
    impG0mats=zero
    impG0real=zero
    impF0mats=zero
    impF0real=zero
    !
    !
    !allocate observables
    allocate(ed_dens(Norb),ed_docc(Norb),ed_phisc(Norb),ed_dens_up(Norb),ed_dens_dw(Norb))
    ed_dens=0d0
    ed_docc=0d0
    ed_phisc=0d0
    ed_dens_up=0d0
    ed_dens_dw=0d0
    !
    allocate(spinChi_tau(Norb+1,0:Ltau))
    allocate(spinChi_w(Norb+1,Lreal))
    allocate(spinChi_iv(Norb+1,0:Lmats))
    !
    allocate(densChi_tau(Norb,Norb,0:Ltau))
    allocate(densChi_w(Norb,Norb,Lreal))
    allocate(densChi_iv(Norb,Norb,0:Lmats))
    allocate(densChi_mix_tau(Norb,Norb,0:Ltau))
    allocate(densChi_mix_w(Norb,Norb,Lreal))
    allocate(densChi_mix_iv(Norb,Norb,0:Lmats))
    allocate(densChi_tot_tau(0:Ltau))
    allocate(densChi_tot_w(Lreal))
    allocate(densChi_tot_iv(0:Lmats))
    !
    allocate(pairChi_tau(Norb,0:Ltau))
    allocate(pairChi_w(Norb,Lreal))
    allocate(pairChi_iv(Norb,0:Lmats))
    !
  end subroutine init_ed_structure





  !+------------------------------------------------------------------+
  !PURPOSE: SETUP THE GLOBAL POINTERS FOR THE ED CALCULAIONS.
  !+------------------------------------------------------------------+
  !
  !NORMAL CASE
  !
  subroutine setup_pointers_normal
    integer                                           :: i,in,dim,isector,jsector
    integer                                           :: nup,ndw,jup,jdw,iorb
    integer                                           :: unit,status,istate
    logical                                           :: IOfile
    integer                                           :: anint
    real(8)                                           :: adouble
    integer                                           :: list_len
    integer,dimension(:),allocatable                  :: list_sector
    isector=0
    do nup=0,Ns
       do ndw=0,Ns
          isector=isector+1
          getSector(nup,ndw)=isector
          getNup(isector)=nup
          getNdw(isector)=ndw
          dim = get_normal_sector_dimension(nup,ndw)
          getDim(isector)=dim
          getDimUp(isector)=get_normal_sector_dimension(nup)
          getDimDw(isector)=get_normal_sector_dimension(ndw)
       enddo
    enddo
    !
    !
    inquire(file="state_list"//reg(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")"Restarting from a state_list file:"
       list_len=file_length("state_list"//reg(ed_file_suffix)//".restart")
       allocate(list_sector(list_len))
       !
       open(free_unit(unit),file="state_list"//reg(ed_file_suffix)//".restart",status="old")
       ! read(unit,*)!read comment line
       ! status=0
       ! do while(status>=0)
       !    read(unit,*,iostat=status)istate,nup,ndw,isector
       !    list_sector(istate)=isector
       !    if(nup/=getnup(isector).OR.ndw/=getndw(isector))&
       !         stop "setup_pointers_normal error: nup!=getnup(isector).OR.ndw!=getndw(isector) "
       ! enddo
       read(unit,*)!read comment line
       status=0
       do while(status>=0)
          !read(unit,"(i6,f18.12,2x,ES19.12,1x,2i3,3x,i3,i10)",iostat=status)istate,adouble,adouble,nup,ndw,isector,anint
          read(unit,*,iostat=status)istate,adouble,adouble,nup,ndw,isector,anint
          list_sector(istate)=isector
          if(nup/=getnup(isector).OR.ndw/=getndw(isector))&
               stop "setup_pointers_normal error: nup!=getnup(isector).OR.ndw!=getndw(isector) "
       enddo
       close(unit)
       !
       lanc_nstates_total = list_len
       do isector=1,Nsectors
          neigen_sector(isector) = max(1,count(list_sector==isector))
       enddo
    else
       do isector=1,Nsectors
          neigen_sector(isector) = min(getdim(isector),lanc_nstates_sector)   !init every sector to required eigenstates
       enddo
    endif
    !
    twin_mask=.true.
    if(ed_twin)then
       ! stop "WARNING: In this updated version with Nup-Ndw factorization the twin-sectors have not been tested!!"
       do isector=1,Nsectors
          nup=getnup(isector)
          ndw=getndw(isector)
          if(nup<ndw)twin_mask(isector)=.false.
       enddo
       write(LOGfile,"(A,I4,A,I4)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
    endif
    !
    do in=1,Norb
       impIndex(in,1)=in
       impIndex(in,2)=in+Ns
    enddo
    !
    select case(bath_type)
    case default
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (iorb-1)*Nbath + i
          enddo
       enddo
    case ('hybrid')
       do i=1,Nbath
          getBathStride(:,i)       = Norb + i
       enddo
    case ('replica')
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = iorb + i*Norb !Norb + (i-1)*Norb + iorb
          enddo
       enddo
    end select
    !
    getCsector=0
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup-1;jdw=ndw;if(jup < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw-1;if(jdw < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(2,isector)=jsector
    enddo
    !
    getCDGsector=0
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup+1;jdw=ndw;if(jup > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw+1;if(jdw > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(2,isector)=jsector
    enddo
  end subroutine setup_pointers_normal

  !
  !SUPERCONDUCTING
  !
  subroutine setup_pointers_superc
    integer                                           :: i,isz,in,dim,isector,jsector
    integer                                           :: sz,iorb,jsz
    integer                                           :: unit,status,istate
    logical                                           :: IOfile
    integer                                           :: anint
    real(8)                                           :: adouble
    integer                                           :: list_len
    integer,dimension(:),allocatable                  :: list_sector
    isector=0
    do isz=-Ns,Ns
       sz=abs(isz)
       isector=isector+1
       getSector(isz,1)=isector
       getSz(isector)=isz
       dim = get_superc_sector_dimension(isz)
       getDim(isector)=dim
    enddo
    !
    !
    inquire(file="state_list"//reg(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       list_len=file_length("state_list"//reg(ed_file_suffix)//".restart")
       allocate(list_sector(list_len))
       !
       open(free_unit(unit),file="state_list"//reg(ed_file_suffix)//".restart",status="old")
       read(unit,*)!read comment line
       status=0
       do while(status>=0)
          !read(unit,"(i6,f18.12,2x,ES19.12,1x,i3,3x,i3,i10)",iostat=status) istate,adouble,adouble,sz,isector,anint
          read(unit,*,iostat=status) istate,adouble,adouble,sz,isector,anint
          list_sector(istate)=isector
          if(sz/=getsz(isector))stop "setup_pointers_superc error: sz!=getsz(isector)."
       enddo
       ! status=0
       ! do while(status>=0)
       !    read(unit,*,iostat=status)istate,sz,isector
       !    list_sector(istate)=isector
       !    if(sz/=getsz(isector))stop "setup_pointers_superc error: sz!=getsz(isector)."
       ! enddo
       close(unit)

       !
       lanc_nstates_total = list_len
       do isector=1,Nsectors
          neigen_sector(isector) = max(1,count(list_sector==isector))
       enddo
    else
       do isector=1,Nsectors
          neigen_sector(isector) = min(getdim(isector),lanc_nstates_sector)   !init every sector to required eigenstates
       enddo
    endif
    twin_mask=.true.
    if(ed_twin)then
       write(LOGfile,*)"USE WITH CAUTION: TWIN STATES IN SC CHANNEL!!";call sleep(1)
       do isector=1,Nsectors
          sz=getsz(isector)
          if(sz>0)twin_mask(isector)=.false.
       enddo
       write(LOGfile,"(A,I4,A,I4)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
    endif
    !
    do in=1,Norb
       impIndex(in,1)=in
       impIndex(in,2)=in+Ns
    enddo
    !
    select case(bath_type)
    case default
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (iorb-1)*Nbath + i
          enddo
       enddo
    case ('hybrid')
       do i=1,Nbath
          getBathStride(:,i)      = Norb + i
       enddo
    case ('replica')
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (i-1)*Norb + iorb
          enddo
       enddo
    end select
    !
    getCsector=0
    !c_up
    do isector=1,Nsectors
       isz=getsz(isector);if(isz==-Ns)cycle
       jsz=isz-1
       jsector=getsector(jsz,1)
       getCsector(1,isector)=jsector
    enddo
    !c_dw
    do isector=1,Nsectors
       isz=getsz(isector);if(isz==Ns)cycle
       jsz=isz+1
       jsector=getsector(jsz,1)
       getCsector(2,isector)=jsector
    enddo
    !
    getCDGsector=0
    !cdg_up
    do isector=1,Nsectors
       isz=getsz(isector);if(isz==Ns)cycle
       jsz=isz+1
       jsector=getsector(jsz,1)
       getCDGsector(1,isector)=jsector
    enddo
    !cdg_dw
    do isector=1,Nsectors
       isz=getsz(isector);if(isz==-Ns)cycle
       jsz=isz-1
       jsector=getsector(jsz,1)
       getCDGsector(2,isector)=jsector
    enddo
  end subroutine setup_pointers_superc

  !
  !NON SU(2) SYMMETRIC
  !
  subroutine setup_pointers_nonsu2
    integer                                           :: i,dim,isector,jsector
    integer                                           :: in,jn,iorb,ispin
    integer                                           :: unit,status,istate
    logical                                           :: IOfile
    integer                                           :: anint
    real(8)                                           :: adouble
    integer                                           :: list_len
    integer,dimension(:),allocatable                  :: list_sector
    integer                                           :: maxtwoJz,twoJz
    integer                                           :: dimJz,inJz,shift
    integer                                           :: twoJz_add,twoJz_del,twoJz_trgt
    isector=0
    if(Jz_basis)then
       !pointers definition
       do in=0,Nlevels
          !
          !algorithm to find the maximum Jz given the density
          if(in==0.or.in==2*Ns)then
             maxtwoJz=0
          else
             shift=0
             if(in<=Nbath+1)shift=Nbath-in+1
             if(in>=2*Ns-Nbath)shift=Nbath-2*Ns+in+1
             maxtwoJz = 5 + 5*Nbath - abs(in-Ns) - 2*shift
          endif
          !
          !number of available Jz given the maximum value
          dimJz = maxtwoJz + 1
          !
          do inJz=1,dimJz
             if(in==0.or.in==2*Ns)then
                twoJz=0
             else
                twoJz = - maxtwoJz + 2*(inJz-1)
             endif
             isector=isector+1
             getN(isector)=in
             gettwoJz(isector)=twoJz
             getmaxtwoJz(in)=maxtwoJz
             getSector(in,twoJz)=isector
             dim = get_nonsu2_sector_dimension_Jz(in,twoJz)
             getDim(isector)=dim
             neigen_sector(isector) = min(dim,lanc_nstates_sector)
          enddo
       enddo
    elseif(plaquette)then
       do in=0,Nlevels
          isector=isector+1
          getSector(in,1)=isector
          getN(isector)=in
          dim = get_normal_sector_dimension(in,0)
          getDim(isector)=dim
          neigen_sector(isector) = min(dim,lanc_nstates_sector)
       enddo
    else
       do in=0,Nlevels
          isector=isector+1
          getSector(in,1)=isector
          getN(isector)=in
          dim = get_nonsu2_sector_dimension(in)
          getDim(isector)=dim
          neigen_sector(isector) = min(dim,lanc_nstates_sector)
       enddo
    endif
    !
    !
    inquire(file="state_list"//reg(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       list_len=file_length("state_list"//reg(ed_file_suffix)//".restart")
       allocate(list_sector(list_len))
       !
       open(free_unit(unit),file="state_list"//reg(ed_file_suffix)//".restart",status="old")
       ! status=0
       ! do while(status>=0)
       !    read(unit,*,iostat=status) istate,in,isector
       !    list_sector(istate)=isector
       !    if(in/=getn(isector))stop "setup_pointers_superc error: n!=getn(isector)."
       ! enddo
       read(unit,*)!read comment line
       status=0
       do while(status>=0)
          !read(unit,"(i6,f18.12,2x,ES19.12,1x,i3,3x,i3,i10)",iostat=status) istate,adouble,adouble,in,isector,anint
          read(unit,*,iostat=status) istate,adouble,adouble,in,isector,anint
          list_sector(istate)=isector
          if(in/=getn(isector))stop "setup_pointers_superc error: n!=getn(isector)."
       enddo
       close(unit)
       !
       lanc_nstates_total = list_len
       do isector=1,Nsectors
          neigen_sector(isector) = max(1,count(list_sector==isector))
       enddo
    else
       do isector=1,Nsectors
          neigen_sector(isector) = min(getdim(isector),lanc_nstates_sector)   !init every sector to required eigenstates
       enddo
    endif
    !
    twin_mask=.true.
    if(ed_twin)then
       write(LOGfile,*)"TWIN STATES IN nonSU2 CHANNEL: NOT TESTED!!"
       call sleep(3)
       do isector=1,Nsectors
          in=getn(isector)
          if(in>Ns)twin_mask(isector)=.false.
          print*,twin_mask(isector),in
       enddo
       write(LOGfile,"(A,I4,A,I4)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
    endif
    !
    do in=1,Norb
       impIndex(in,1)=in
       impIndex(in,2)=in+Ns
    enddo
    !
    select case(bath_type)
    case default
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (iorb-1)*Nbath + i
          enddo
       enddo
    case ('hybrid')
       do i=1,Nbath
          getBathStride(:,i)      = Norb + i
       enddo
    case ('replica')
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (i-1)*Norb + iorb
          enddo
       enddo
    end select
    !
    getCsector=0
    !c_{up,dw}
    do isector=1,Nsectors
       in=getn(isector);if(in==0)cycle
       jn=in-1
       jsector=getsector(jn,1)
       getCsector(1,isector)=jsector
       getCsector(2,isector)=jsector
    enddo
    !
    getCDGsector=0
    !cdg_{up,dw}
    do isector=1,Nsectors
       in=getn(isector);if(in==Nlevels)cycle
       jn=in+1
       jsector=getsector(jn,1)
       getCDGsector(1,isector)=jsector
       getCDGsector(2,isector)=jsector
    enddo

    if(Jz_basis)then
       !
       getCsector_Jz=-1
       !c_{Lz,Sz}
       do isector=1,Nsectors
          in=getn(isector);if(in==0)cycle
          jn=in-1
          !
          twoJz=gettwoJz(isector)
          do iorb=1,Norb
             do ispin=1,Nspin
                twoJz_del  = 2 * Lzdiag(iorb) + Szdiag(ispin)
                twoJz_trgt = twoJz - twoJz_del
                if(abs(twoJz_trgt) > getmaxtwoJz(jn)) cycle
                jsector=getSector(jn,twoJz_trgt)
                getCsector_Jz(iorb,ispin,isector)=jsector
             enddo
          enddo
       enddo
       !
       getCDGsector_Jz=-1
       !cdg_{Lz,Sz}
       do isector=1,Nsectors
          in=getn(isector);if(in==Nlevels)cycle
          jn=in+1
          !
          twoJz=gettwoJz(isector)
          do iorb=1,Norb
             do ispin=1,Nspin
                twoJz_add  = 2 * Lzdiag(iorb) + Szdiag(ispin)
                twoJz_trgt = twoJz + twoJz_add
                if(abs(twoJz_trgt) > getmaxtwoJz(jn)) cycle
                jsector=getSector(jn,twoJz_trgt)
                getCDGsector_Jz(iorb,ispin,isector)=jsector
             enddo
          enddo
       enddo
       !
    endif
  end subroutine setup_pointers_nonsu2





  !+------------------------------------------------------------------+
  !PURPOSE  : return the dimension of a sector
  !+------------------------------------------------------------------+
  !NORMAL
  function get_normal_sector_dimension(nup,ndw) result(dim)
    integer :: nup
    integer,optional :: ndw
    integer :: dim,dimup,dimdw
    if(present(ndw))then
       dimup = binomial(Ns,nup)    !this ensures better evaluation of the dimension
       dimdw = binomial(Ns,ndw)    !as it avoids large numbers
    else
       dimup = binomial(Ns,nup)
       dimdw = 1
    endif
    dim=dimup*dimdw
  end function get_normal_sector_dimension
  !SUPERC
  function get_superc_sector_dimension(mz) result(dim)
    integer :: mz
    integer :: i,dim,Nb
    dim=0
    Nb=Ns-mz
    do i=0,Nb/2
       dim=dim + 2**(Nb-2*i)*binomial(ns,Nb-2*i)*binomial(ns-Nb+2*i,i)
    enddo
  end function get_superc_sector_dimension
  !NONSU2
  function get_nonsu2_sector_dimension(n) result(dim)
    integer :: n
    integer :: dim
    dim=binomial(2*Ns,n)
  end function get_nonsu2_sector_dimension
  !NONSU2 - Jz conserving
  function get_nonsu2_sector_dimension_Jz(n,twoJz) result(dim)
    integer :: n
    integer :: twoJz
    integer :: dim
    integer :: ivec(Ns),jvec(Ns)
    integer :: iup,idw,ibath,iorb
    integer :: nt,twoLz,twoSz
    !
    dim=0
    do idw=0,2**Ns-1
       jvec = bdecomp(idw,Ns)
       do iup=0,2**Ns-1
          ivec = bdecomp(iup,Ns)
          nt   = sum(ivec) + sum(jvec)
          twoLz=0;twoSz=0
          do ibath=0,Nbath
             do iorb=1,Norb
                twoLz = twoLz + 2 * Lzdiag(iorb) * ivec(iorb+Norb*ibath)  &
                     + 2 * Lzdiag(iorb) * jvec(iorb+Norb*ibath)
             enddo
          enddo
          twoSz = (sum(ivec) - sum(jvec))
          !
          if(nt == n .and. twoJz==(twoSz+twoLz) )then
             dim=dim+1
          endif
       enddo
    enddo
  end function get_nonsu2_sector_dimension_Jz




  !+------------------------------------------------------------------+
  !PURPOSE  : constructs the sectors by storing the map to the
  !states i\in Hilbert_space from the states count in H_sector.
  !|ImpUP,BathUP>|ImpDW,BathDW >
  !+------------------------------------------------------------------+
  subroutine build_sector(isector,Hup,Hup8)
    integer                                      :: isector
    type(sector_map)                             :: Hup
    integer                                      :: nup,ndw,sz,nt,twoJz
    integer                                      :: nup_,ndw_,sz_,nt_
    integer                                      :: twoSz_,twoLz_
    integer                                      :: i,ibath,iorb
    integer                                      :: iup,idw
    integer                                      :: p1,p2,p3,p4
    integer                                      :: p5,p6,p7,p8
    integer                                      :: p9,p10,p11,p12
    integer                                      :: p13,p14,p15,p16
    integer                                      :: dim,idim
    integer                                      :: ivec(Ns),jvec(Ns)
    !
    integer                                      :: Nnn,isite,ineig,dimUsed
    integer(16),parameter                        :: zr=0
    integer(16),allocatable                      :: Order(:)
    type(sector_map8),optional                   :: Hup8
    type(sector_map8)                            :: Hup8Full
    select case(ed_mode)
       !
       !
    case default
       nup = getNup(isector)
       ndw = getNdw(isector)
       dim = getDim(isector)
       call map_allocate(Hup,dim)
       dim=0
       do idw=0,2**Ns-1
          jvec  = bdecomp(idw,Ns)
          ndw_  = sum(jvec)
          if(ndw_ /= ndw)cycle
          do iup=0,2**Ns-1
             ivec  = bdecomp(iup,Ns)
             nup_  = sum(ivec)
             if(nup_ /= nup)cycle
             dim      = dim+1
             Hup%map(dim) = iup + idw*2**Ns
          enddo
       enddo
       !
       !
    case ("superc")
       sz  = getSz(isector)
       dim = getDim(isector)
       call map_allocate(Hup,dim)
       dim=0
       do idw=0,2**Ns-1
          jvec = bdecomp(idw,Ns)
          ndw_ = sum(jvec)
          do iup=0,2**Ns-1
             ivec = bdecomp(iup,Ns)
             nup_ = sum(ivec)
             sz_  = nup_ - ndw_
             if(sz_ == sz)then
                dim=dim+1
                Hup%map(dim)=iup + idw*2**Ns
             endif
          enddo
       enddo
       !
       !
    case ("nonsu2")
       if(Jz_basis)then
          nt  = getN(isector)
          dim = getDim(isector)
          twoJz = gettwoJz(isector)
          call map_allocate(Hup,dim)
          dim=0
          do idw=0,2**Ns-1
             jvec = bdecomp(idw,Ns)
             do iup=0,2**Ns-1
                ivec = bdecomp(iup,Ns)
                nt_  = sum(ivec) + sum(jvec)
                twoLz_=0;twoSz_=0
                do ibath=0,Nbath
                   do iorb=1,Norb
                      twoLz_ = twoLz_ + 2 * Lzdiag(iorb) * ivec(iorb+Norb*ibath)  &
                           + 2 * Lzdiag(iorb) * jvec(iorb+Norb*ibath)
                   enddo
                enddo
                twoSz_ = (sum(ivec) - sum(jvec))
                !
                if(nt_ == nt .and. twoJz==(twoSz_+twoLz_) )then
                   dim=dim+1
                   Hup%map(dim)=iup + idw*2**Ns
                endif
             enddo
          enddo
       elseif(plaquette)then
          !
          nt  = getN(isector)
          if((nt.gt.Ns).or.(nt.gt.16))stop "filling too big"
          dim = get_normal_sector_dimension(nt,0)
          call map_allocate8(Hup8Full,dim)
          !
          if(filling.eq.2)then
             write(LOGfile,'(A)')"Building N=2 sector"
             dim=0
             do p1=0,Ns-2
                do p2=p1+1,Ns-1
                   dim=dim+1
                   Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 )
                   if(Hup8Full%map(dim).lt.0) stop
                enddo
             enddo
          elseif(filling.eq.3)then
             write(LOGfile,'(A)')"Building N=3 sector"
             dim=0
             do p1=0,Ns-2
                do p2=p1+1,Ns-1
                   do p3=p2+1,Ns-1
                      dim=dim+1
                      Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 )
                      if(Hup8Full%map(dim).lt.0) stop
                   enddo
                enddo
             enddo
          elseif(filling.eq.4)then
             write(LOGfile,'(A)')"Building N=4 sector"
             dim=0
             do p1=0,Ns-2
                do p2=p1+1,Ns-1
                   do p3=p2+1,Ns-1
                      do p4=p3+1,Ns-1
                         dim=dim+1
                         Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 )
                         if(Hup8Full%map(dim).lt.0) stop
                      enddo
                   enddo
                enddo
             enddo
          elseif(filling.eq.5)then
             write(LOGfile,'(A)')"Building N=5 sector"
             dim=0
             do p1=0,Ns-2
                do p2=p1+1,Ns-1
                   do p3=p2+1,Ns-1
                      do p4=p3+1,Ns-1
                         do p5=p4+1,Ns-1
                            dim=dim+1
                            Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 ) &
                                              + ibset( zr, p5 )
                            if(Hup8Full%map(dim).lt.0) stop
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          elseif(filling.eq.6)then
             write(LOGfile,'(A)')"Building N=6 sector"
             dim=0
             do p1=0,Ns-2
                do p2=p1+1,Ns-1
                   do p3=p2+1,Ns-1
                      do p4=p3+1,Ns-1
                         do p5=p4+1,Ns-1
                            do p6=p5+1,Ns-1
                               dim=dim+1
                               Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 ) &
                                                 + ibset( zr, p5 ) + ibset( zr, p6 )
                               if(Hup8Full%map(dim).lt.0) stop
                            enddo
                         enddo
                     enddo
                  enddo
               enddo
            enddo
         elseif(filling.eq.7)then
             write(LOGfile,'(A)')"Building N=7 sector"
             dim=0
             do p1=0,Ns-2
                do p2=p1+1,Ns-1
                   do p3=p2+1,Ns-1
                      do p4=p3+1,Ns-1
                         do p5=p4+1,Ns-1
                            do p6=p5+1,Ns-1
                               do p7=p6+1,Ns-1
                                  dim=dim+1
                                  Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 ) &
                                                    + ibset( zr, p5 ) + ibset( zr, p6 ) + ibset( zr, p7 )
                                  if(Hup8Full%map(dim).lt.0) stop
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          elseif(filling.eq.8)then
             write(LOGfile,'(A)')"Building N=8 sector"
             dim=0
             do p1=0,Ns-2
                do p2=p1+1,Ns-1
                   do p3=p2+1,Ns-1
                      do p4=p3+1,Ns-1
                         do p5=p4+1,Ns-1
                            do p6=p5+1,Ns-1
                               do p7=p6+1,Ns-1
                                  do p8=p7+1,Ns-1
                                     dim=dim+1
                                     Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 ) &
                                                       + ibset( zr, p5 ) + ibset( zr, p6 ) + ibset( zr, p7 ) + ibset( zr, p8 )
                                     if(Hup8Full%map(dim).lt.0) stop
                                  enddo
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          elseif(filling.eq.9)then
             write(LOGfile,'(A)')"Building N=9 sector"
             dim=0
             do p1=0,Ns-2
                do p2=p1+1,Ns-1
                   do p3=p2+1,Ns-1
                      do p4=p3+1,Ns-1
                         do p5=p4+1,Ns-1
                            do p6=p5+1,Ns-1
                               do p7=p6+1,Ns-1
                                  do p8=p7+1,Ns-1
                                     do p9=p8+1,Ns-1
                                        dim=dim+1
                                        Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 ) &
                                                          + ibset( zr, p5 ) + ibset( zr, p6 ) + ibset( zr, p7 ) + ibset( zr, p8 ) &
                                                          + ibset( zr, p9 )
                                        if(Hup8Full%map(dim).lt.0) stop
                                     enddo
                                  enddo
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          elseif(filling.eq.10)then
             write(LOGfile,'(A)')"Building N=10 sector"
             dim=0
             do p1=0,Ns-2
                do p2=p1+1,Ns-1
                   do p3=p2+1,Ns-1
                      do p4=p3+1,Ns-1
                         do p5=p4+1,Ns-1
                            do p6=p5+1,Ns-1
                               do p7=p6+1,Ns-1
                                  do p8=p7+1,Ns-1
                                     do p9=p8+1,Ns-1
                                        do p10=p9+1,Ns-1
                                           dim=dim+1
                                           Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 ) &
                                                             + ibset( zr, p5 ) + ibset( zr, p6 ) + ibset( zr, p7 ) + ibset( zr, p8 ) &
                                                             + ibset( zr, p9 ) + ibset( zr, p10)
                                           if(Hup8Full%map(dim).lt.0) stop
                                        enddo
                                     enddo
                                  enddo
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          elseif(filling.eq.11)then
             write(LOGfile,'(A)')"Building N=11 sector"
             dim=0
             do p1=0,Ns-2
                do p2=p1+1,Ns-1
                   do p3=p2+1,Ns-1
                      do p4=p3+1,Ns-1
                         do p5=p4+1,Ns-1
                            do p6=p5+1,Ns-1
                               do p7=p6+1,Ns-1
                                  do p8=p7+1,Ns-1
                                     do p9=p8+1,Ns-1
                                        do p10=p9+1,Ns-1
                                           do p11=p10+1,Ns-1
                                              dim=dim+1
                                              Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 ) &
                                                                + ibset( zr, p5 ) + ibset( zr, p6 ) + ibset( zr, p7 ) + ibset( zr, p8 ) &
                                                                + ibset( zr, p9 ) + ibset( zr, p10) + ibset( zr, p11)
                                              if(Hup8Full%map(dim).lt.0) stop
                                           enddo
                                        enddo
                                     enddo
                                  enddo
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          elseif(filling.eq.12)then
             write(LOGfile,'(A)')"Building N=12 sector"
             dim=0
             do p1=0,Ns-2
                do p2=p1+1,Ns-1
                   do p3=p2+1,Ns-1
                      do p4=p3+1,Ns-1
                         do p5=p4+1,Ns-1
                            do p6=p5+1,Ns-1
                               do p7=p6+1,Ns-1
                                  do p8=p7+1,Ns-1
                                     do p9=p8+1,Ns-1
                                        do p10=p9+1,Ns-1
                                           do p11=p10+1,Ns-1
                                              do p12=p11+1,Ns-1
                                                 dim=dim+1
                                                 Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 ) &
                                                                   + ibset( zr, p5 ) + ibset( zr, p6 ) + ibset( zr, p7 ) + ibset( zr, p8 ) &
                                                                   + ibset( zr, p9 ) + ibset( zr, p10) + ibset( zr, p11) + ibset( zr, p12)
                                                 if(Hup8Full%map(dim).lt.0) stop
                                              enddo
                                           enddo
                                        enddo
                                     enddo
                                  enddo
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          elseif(filling.eq.13)then
            write(LOGfile,'(A)')"Building N=13 sector"
            dim=0
            do p1=0,Ns-2
            do p2=p1+1,Ns-1
            do p3=p2+1,Ns-1
            do p4=p3+1,Ns-1
            do p5=p4+1,Ns-1
            do p6=p5+1,Ns-1
            do p7=p6+1,Ns-1
            do p8=p7+1,Ns-1
            do p9=p8+1,Ns-1
            do p10=p9+1,Ns-1
            do p11=p10+1,Ns-1
            do p12=p11+1,Ns-1
            do p13=p12+1,Ns-1
            do p14=p13+1,Ns-1
                dim=dim+1
                Hup8Full%map(dim) =  ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 ) &
                                   + ibset( zr, p5 ) + ibset( zr, p6 ) + ibset( zr, p7 ) + ibset( zr, p8 ) &
                                   + ibset( zr, p9 ) + ibset( zr, p10) + ibset( zr, p11) + ibset( zr, p12) &
                                   + ibset( zr, p13)
                                   if(Hup8Full%map(dim).lt.0) stop
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
          elseif(filling.eq.14)then
            write(LOGfile,'(A)')"Building N=14 sector"
            dim=0
            do p1=0,Ns-2
            do p2=p1+1,Ns-1
            do p3=p2+1,Ns-1
            do p4=p3+1,Ns-1
            do p5=p4+1,Ns-1
            do p6=p5+1,Ns-1
            do p7=p6+1,Ns-1
            do p8=p7+1,Ns-1
            do p9=p8+1,Ns-1
            do p10=p9+1,Ns-1
            do p11=p10+1,Ns-1
            do p12=p11+1,Ns-1
            do p13=p12+1,Ns-1
            do p14=p13+1,Ns-1
               dim=dim+1
               Hup8Full%map(dim) =  ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 ) &
                                  + ibset( zr, p5 ) + ibset( zr, p6 ) + ibset( zr, p7 ) + ibset( zr, p8 ) &
                                  + ibset( zr, p9 ) + ibset( zr, p10) + ibset( zr, p11) + ibset( zr, p12) &
                                  + ibset( zr, p13) + ibset( zr, p14)
                                  if(Hup8Full%map(dim).lt.0) stop
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
          elseif(filling.eq.15)then
             write(LOGfile,'(A)')"Building N=15 sector"
             dim=0
             do p1=0,Ns-2
             do p2=p1+1,Ns-1
             do p3=p2+1,Ns-1
             do p4=p3+1,Ns-1
             do p5=p4+1,Ns-1
             do p6=p5+1,Ns-1
             do p7=p6+1,Ns-1
             do p8=p7+1,Ns-1
             do p9=p8+1,Ns-1
             do p10=p9+1,Ns-1
             do p11=p10+1,Ns-1
             do p12=p11+1,Ns-1
             do p13=p12+1,Ns-1
             do p14=p13+1,Ns-1
             do p15=p14+1,Ns-1
                dim=dim+1
                Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 ) &
                                  + ibset( zr, p5 ) + ibset( zr, p6 ) + ibset( zr, p7 ) + ibset( zr, p8 ) &
                                  + ibset( zr, p9 ) + ibset( zr, p10) + ibset( zr, p11) + ibset( zr, p12) &
                                  + ibset( zr, p13) + ibset( zr, p14) + ibset( zr, p15)
                                  if(Hup8Full%map(dim).lt.0) stop
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
          elseif(filling.eq.16)then
             write(LOGfile,'(A)')"Building N=16 sector"
             dim=0
             do p1=0,Ns-2
             do p2=p1+1,Ns-1
             do p3=p2+1,Ns-1
             do p4=p3+1,Ns-1
             do p5=p4+1,Ns-1
             do p6=p5+1,Ns-1
             do p7=p6+1,Ns-1
             do p8=p7+1,Ns-1
             do p9=p8+1,Ns-1
             do p10=p9+1,Ns-1
             do p11=p10+1,Ns-1
             do p12=p11+1,Ns-1
             do p13=p12+1,Ns-1
             do p14=p13+1,Ns-1
             do p15=p14+1,Ns-1
             do p16=p15+1,Ns-1
                dim=dim+1
                Hup8Full%map(dim) = ibset( zr, p1 ) + ibset( zr, p2 ) + ibset( zr, p3 ) + ibset( zr, p4 ) &
                                  + ibset( zr, p5 ) + ibset( zr, p6 ) + ibset( zr, p7 ) + ibset( zr, p8 ) &
                                  + ibset( zr, p9 ) + ibset( zr, p10) + ibset( zr, p11) + ibset( zr, p12) &
                                  + ibset( zr, p13) + ibset( zr, p14) + ibset( zr, p15) + ibset( zr, p16)
                                  if(Hup8Full%map(dim).lt.0) stop
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo
          endif
          !
          dim = get_normal_sector_dimension(nt,0)
          !
          if(HardCoreBoson.gt.1)then
             !
             allocate(Order(dim));Order=0
             dimUsed=0
             do idim=1,dim
                !
                ivec = bdecomp8(Hup8Full%map(idim),Ns)
                !
                Nnn=0
                do isite=1,Norb*Nbath
                   if(ivec(isite).eq.1)then
                      !Nearest neighbors
                      do ineig=1,size(Vstride(isite,1,:))
                         Nnn = Nnn + ivec(Vstride(isite,1,ineig))
                      enddo
                      !Next-Nearest neighbors
                      if(Nbath.ne.1)then
                        do ineig=1,size(Vstride(isite,2,:))
                           Nnn = Nnn + ivec(Vstride(isite,2,ineig))
                        enddo
                      endif
                   endif
                enddo
                !
                if(Nnn.eq.0)then
                   dimUsed = dimUsed+1
                   Order(dimUsed) = Hup8Full%map(idim)
                endif
             enddo
             write(*,'(A,I10)')"Sector Full dimension: ",dim
             write(*,'(A,I10)')"Sector dimension for HardCore quartet: ",dimUsed
             getDim(isector)=dimUsed
             !
             call map_allocate8(Hup8,dimUsed)
             do idim=1,dimUsed
                Hup8%map(idim) = Order(idim)
             enddo
             write(*,'(A,I10)')"size(Hup8%map): ",size(Hup8%map)
             deallocate(Order)
             call map_deallocate8(Hup8Full)
             dim = dimUsed
             !
          else
             !
             write(*,'(A,I10)')"Sector Full dimension: ",dim
             call map_allocate8(Hup8,dim)
             Hup8%map = Hup8Full%map
             call map_deallocate8(Hup8Full)
             !
          endif
          !
          allocate(Order(dim))
          call sort_array8(Hup8%map,Order)
          !dim = get_normal_sector_dimension(nt,0)
          do iup=1,dim
             jvec = bdecomp8(Hup8%map(iup),Ns)
             if (sum(jvec).ne.filling) then
                write(*,'(3I15,1000I3)')iup,dim,Hup8%map(iup),Ns,jvec
                stop
             endif
          enddo
          !
       else
          nt  = getN(isector)
          dim = getDim(isector)
          call map_allocate(Hup,dim)
          dim=0
          do idw=0,2**Ns-1
             jvec = bdecomp(idw,Ns)
             do iup=0,2**Ns-1
                ivec = bdecomp(iup,Ns)
                nt_  = sum(ivec) + sum(jvec)
                if(nt_ == nt)then
                   dim=dim+1
                   Hup%map(dim)=iup + idw*2**Ns
                endif
             enddo
          enddo
       endif
    end select
  end subroutine build_sector






  subroutine build_sector_2(isector)
    integer                   :: isector
    integer                   :: nup,ndw,sz,nt
    integer                   :: nup_,ndw_,sz_,nt_
    integer                   :: i,ibath
    integer                   :: iup,idw
    integer                   :: dim
    real(8)                   :: Sz_tot,Lz_tot,shift,stride,jzv
    real(8),allocatable       :: Jz(:)
    integer                   :: ivec(Ns),jvec(Ns)
    !
    !
    stride=0.d0
    nt  = getN(isector)
    dim = getDim(isector)

    write(123,*)
    write(123,'(A20,I5)') "sector N",nt
    write(123,'(A20,I5)') "sector dim(N)",dim
    write(123,*)

    write(124,*)
    write(124,'(A20,I5)') "sector N",nt
    write(124,'(A20,I5)') "sector dim(N)",dim

    if(allocated(Jz))deallocate(Jz);allocate(Jz(dim));Jz=0.d0

    dim=0
    !mi guardo tutti i possibili valori di nup,ndw
    !anche quelli con la densità diversa da quella del settore
    do idw=0,2**Ns-1
       jvec = bdecomp(idw,Ns)
       do iup=0,2**Ns-1
          ivec = bdecomp(iup,Ns)
          !ivec e jvec sono la decomposizione in vettore
          !qui controllo la densità che ho ottenuto con lo specifico vettore
          nt_  = sum(ivec) + sum(jvec)
          Lz_tot=0.d0;Sz_tot=0.d0
          do ibath=0,Nbath
             Lz_tot = Lz_tot + 1.d0 * ivec(1+Norb*ibath) + 1.d0 * jvec(1+Norb*ibath)
             Lz_tot = Lz_tot - 1.d0 * ivec(2+Norb*ibath) - 1.d0 * jvec(2+Norb*ibath)
          enddo
          Sz_tot = 0.5*(sum(ivec) - sum(jvec))

          !se la densità è quella del settore metto lo stato in rappresentazione decimale dentro la mappa
          if(nt_ == nt)then
             dim=dim+1
             Jz(dim)=(Lz_tot+Sz_tot)
             !write(123,'(3I3,4X,3I3,4X,3(1F5.2,4X))')ivec,jvec,Lz_tot,Sz_tot,Jz(dim)
             !write(123,'(6I3,4X,6I3,4X,3(1F5.2,4X))')ivec,jvec,Lz_tot,Sz_tot,Jz(dim)
             write(123,'(9I3,4X,9I3,4X,3(1F5.2,4X),9I5)')ivec,jvec,Lz_tot,Sz_tot,Jz(dim),dim
          endif
       enddo
    enddo
    !
    ! algorithm
    if(nt==0.or.nt==2*Ns)then
       jzv=0
    else
       shift=0.
       if(nt<=Nbath+1)shift=Nbath-nt+1
       if(nt>=2*Ns-Nbath)shift=Nbath-2*Ns+nt+1
       jzv = 5/2. +(5/2.*Nbath) -(1/2.)*abs(nt-Ns)-shift
    endif
    !
    write(124,*)
    write(124,*)"-----------------------------"
    write(124,'(2(A20,3X,1F5.2),30F5.2)')"maxval(Jz)",maxval(Jz),"algorithm",Jzv
    write(124,*)
    write(124,'(2(A20,3X,1F5.2))')"minval(Jz)",minval(Jz)
    write(124,*)
    write(124,'(2(A20,3X,1I5))')"degeneracy",int(2.d0*maxval(Jz))+1
    write(124,*)"-----------------------------"
    write(124,*)
  end subroutine build_sector_2

  subroutine delete_sector(isector,Hup)!,Hdw)
    integer                   :: isector
    type(sector_map)          :: Hup
    ! type(sector_map),optional :: Hdw
    call map_deallocate(Hup)
  end subroutine delete_sector

  subroutine delete_sector8(isector,Hup)!,Hdw)
    integer                   :: isector
    type(sector_map8)          :: Hup
    ! type(sector_map),optional :: Hdw
    call map_deallocate8(Hup)
  end subroutine delete_sector8



  !+-------------------------------------------------------------------+
  !PURPOSE: input state |in> of the basis and calculates
  !   |out>=C_pos|in>  OR  |out>=C^+_pos|in> ;
  !   the sign of |out> has the phase convention, pos labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn
    integer               :: l
    if(.not.btest(in,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in,pos-1)
  end subroutine c

  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn
    integer               :: l
    if(btest(in,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibset(in,pos-1)
  end subroutine cdg

  subroutine c8(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer(16),intent(in)    :: in
    integer(16),intent(inout) :: out
    real(8),intent(inout) :: fsgn
    integer               :: l
    if(.not.btest(in,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in,pos-1)
  end subroutine c8

  subroutine cdg8(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer(16),intent(in)    :: in
    integer(16),intent(inout) :: out
    real(8),intent(inout) :: fsgn
    integer               :: l
    if(btest(in,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibset(in,pos-1)
  end subroutine cdg8






  !##################################################################
  !##################################################################
  !TWIN SECTORS ROUTINES:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : Build the re-ordering map to go from sector A(nup,ndw)
  ! to its twin sector B(ndw,nup), with nup!=ndw.
  !+------------------------------------------------------------------+
  subroutine twin_sector_order(isector,order)
    integer                          :: isector
    integer,dimension(:)             :: order
    type(sector_map)                 :: H,Hup,Hdw
    integer                          :: i,dim
    dim = getdim(isector)
    if(size(Order)/=dim)stop "twin_sector_order error: wrong dimensions of *order* array"
    !- build the map from the A-sector to \HHH
    !- get the list of states in \HHH corresponding to sector B twin of A
    !- return the ordering of B-states in \HHH with respect to those of A
    select case(ed_mode)
    case default
       call build_sector(isector,H)
       do i=1,dim
          Order(i)=flip_state(H%map(i))
       enddo
       call sort_array(Order)
       deallocate(H%map)
       !
    case('normal')
       call build_sector(isector,H)
       do i=1,dim
          Order(i)=flip_state(H%map(i))
       enddo
       call sort_array(Order)
       deallocate(H%map)
       !
    end select
  end subroutine twin_sector_order



  !+------------------------------------------------------------------+
  !PURPOSE  : Flip an Hilbert space state m=|{up}>|{dw}> into:
  !
  ! normal: j=|{dw}>|{up}>  , nup --> ndw
  ! superc: j=|{dw}>|{up}>  , sz  --> -sz
  ! nonsu2: j=|{!up}>|{!dw}>, n   --> 2*Ns-n
  !+------------------------------------------------------------------+
  function flip_state(m,n) result(j)
    integer          :: m
    integer,optional :: n
    integer          :: j
    integer          :: ivec(2*Ns),foo(2*Ns),ivup(Ns),ivdw(Ns)
    select case(ed_mode)
    case default    !Exchange UP-config |{up}> with DW-config |{dw}>
       ! Ivup = bdecomp(m,Ns)
       ! Ivdw = bdecomp(n,Ns)
       Ivec = bdecomp(m,2*Ns)
       foo(1:Ns)     =Ivec(Ns+1:2*Ns)!Ivdw
       foo(Ns+1:2*Ns)=Ivec(1:Ns)     !Ivup
    case("superc")  !Invert the overall spin sign: |{up}> <---> |{dw}>
       Ivec = bdecomp(m,2*Ns)
       foo(1:Ns)     =Ivec(Ns+1:2*Ns)
       foo(Ns+1:2*Ns)=Ivec(1:Ns)
    case ("nonsu2") !Exchange Occupied sites (1) with Empty sites (0)
       Ivec = bdecomp(m,2*Ns)
       where(Ivec==1)foo=0
       where(Ivec==0)foo=1
    end select
    !
    j = bjoin(foo,2*Ns)
    !
  end function flip_state


  !+------------------------------------------------------------------+
  !PURPOSE  : get the twin of a given sector (the one with opposite
  ! quantum numbers):
  ! nup,ndw ==> ndw,nup (spin-exchange)
  ! sz      ==> -sz     (total spin flip)
  ! n       ==> 2*Ns-n  (particle hole)
  !+------------------------------------------------------------------+
  function get_twin_sector(isector) result(jsector)
    integer,intent(in) :: isector
    integer :: jsector
    integer :: iup,idw,in,isz
    select case(ed_mode)
    case default
       iup=getnup(isector)
       idw=getndw(isector)
       jsector=getsector(idw,iup)
    case ("superc")
       isz=getsz(isector)
       jsector=getsector(-isz,1)
    case("nonsu2")
       in=getn(isector)
       jsector=getsector(Nlevels-in,1)
    end select
  end function get_twin_sector










  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: Ntot,ivec(Ntot),l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp

  function bdecomp8(i,Ntot) result(ivec)
    integer(16) :: i
    integer :: Ntot,ivec(Ntot),l
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp8



  !+------------------------------------------------------------------+
  !PURPOSE  : input a vector ib(Nlevels) with the binary sequence
  ! and output the corresponding state |i>
  !(corresponds to the recomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bjoin(ib,Ntot) result(i)
    integer                 :: Ntot
    integer,dimension(Ntot) :: ib
    integer                 :: i,j
    i=0
    do j=0,Ntot-1
       i=i+ib(j+1)*2**j
    enddo
  end function bjoin



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial of an integer N!=1.2.3...(N-1).N
  !+------------------------------------------------------------------+
  recursive function factorial(n) result(f)
    integer            :: f
    integer,intent(in) :: n
    if(n<=0)then
       f=1
    else
       f=n*factorial(n-1)
    end if
  end function factorial



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  function binomial(n1,n2) result(nchoos)
    real(8) :: xh
    integer :: n1,n2,i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial



  !+------------------------------------------------------------------+
  !PURPOSE : binary search of a value in an array
  !+------------------------------------------------------------------+
  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       !stop "binary_search error: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search

  recursive function binary_search8(a,value) result(bsresult)
    integer(16),intent(in) :: a(:), value
    integer(16)            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       !stop "binary_search error: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search8(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search8(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search8



  !+------------------------------------------------------------------+
  !PURPOSE : sort array of integer using random algorithm
  !+------------------------------------------------------------------+
  subroutine sort_array(array)
    integer,dimension(:),intent(inout)      :: array
    integer,dimension(size(array))          :: order
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort( array, order, 1, size(array) )
    array=order
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:)                 :: array
      integer, dimension(:)                 :: order
      integer                               :: left
      integer                               :: right
      integer                               :: i
      integer                               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:)                 :: order
      integer                               :: first, second
      integer                               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    function qsort_rand( lower, upper )
      implicit none
      integer                               :: lower, upper
      real(8)                               :: r
      integer                               :: qsort_rand
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    function compare(f,g)
      integer                               :: f,g
      integer                               :: compare
      compare=1
      if(f<g)compare=-1
    end function compare
  end subroutine sort_array

  subroutine sort_array8(array,order)
    implicit none
    integer(16),dimension(:)                    :: array
    integer(16),dimension(size(array))          :: order
    integer(16),dimension(size(array))          :: backup
    integer(16)                                 :: i
    integer(16)                                 :: lf,rg
    lf=1
    rg=size(array)
    forall(i=1:size(array))order(i)=i
    call qsort_sort(array, order,lf, rg)
    do i=1,size(array)
       backup(i)=array(order(i))
    enddo
    array=backup
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer(16), dimension(:) :: array
      integer(16), dimension(:) :: order
      integer(16)               :: left
      integer(16)               :: right
      integer(16)               :: i
      integer(16)               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer(16), dimension(:) :: order
      integer(16)               :: first, second
      integer(16)               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    integer(16) function qsort_rand( lower, upper )
      integer(16)            :: lower, upper
      real(8)               :: r
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    !---------------------------------------------!
    function compare(f,g)
      implicit none
      integer(16)               :: f,g
      integer(16)               :: compare
      if(f<g) then
         compare=-1
      else
         compare=1
      endif
    end function compare
  end subroutine sort_array8



#ifdef _MPI
  !! Scatter V into the arrays Vloc on each thread: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine scatter_vector_MPI(MpiComm,v,vloc)
    integer                          :: MpiComm
    complex(8),dimension(:)          :: v    !size[N]
    complex(8),dimension(:)          :: vloc !size[Nloc]
    integer                          :: i,irank,Nloc,N
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
    if( MpiComm == MPI_UNDEFINED ) stop "scatter_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "scatter_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    Vloc=0
    call MPI_Scatterv(V,Counts,Offset,MPI_DOUBLE_COMPLEX,Vloc,Nloc,MPI_DOUBLE_COMPLEX,0,MpiComm,MpiIerr)
    !
    return
  end subroutine scatter_vector_MPI


  subroutine scatter_basis_MPI(MpiComm,v,vloc)
    integer                   :: MpiComm
    complex(8),dimension(:,:) :: v    !size[N,N]
    complex(8),dimension(:,:) :: vloc !size[Nloc,Neigen]
    integer                   :: N,Nloc,Neigen,i
    N      = size(v,1)
    Nloc   = size(vloc,1)
    Neigen = size(vloc,2)
    if( size(v,2) < Neigen ) stop "error scatter_basis_MPI: size(v,2) < Neigen"
    !
    do i=1,Neigen
       call scatter_vector_MPI(MpiComm,v(:,i),vloc(:,i))
    end do
    !
    return
  end subroutine scatter_basis_MPI


  !! AllGather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine gather_vector_MPI(MpiComm,vloc,v)
    integer                          :: MpiComm
    complex(8),dimension(:)          :: vloc !size[Nloc]
    complex(8),dimension(:)          :: v    !size[N]
    integer                          :: i,irank,Nloc,N
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
    if( MpiComm == MPI_UNDEFINED ) stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "gather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    call MPI_Gatherv(Vloc,Nloc,MPI_DOUBLE_COMPLEX,V,Counts,Offset,MPI_DOUBLE_COMPLEX,0,MpiComm,MpiIerr)
    !
    return
  end subroutine gather_vector_MPI


  !! AllGather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine allgather_vector_MPI(MpiComm,vloc,v)
    integer                          :: MpiComm
    complex(8),dimension(:)          :: vloc !size[Nloc]
    complex(8),dimension(:)          :: v    !size[N]
    integer                          :: i,irank,Nloc,N
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
    if( MpiComm == MPI_UNDEFINED ) stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "gather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    call MPI_AllGatherv(Vloc,Nloc,MPI_DOUBLE_COMPLEX,V,Counts,Offset,MPI_DOUBLE_COMPLEX,MpiComm,MpiIerr)
    !
    return
  end subroutine Allgather_vector_MPI
#endif





end MODULE ED_SETUP






! !+------------------------------------------------------------------+
! !PURPOSE  : return the Number of sectors using local variables
! ! This is similar to setup Dimensions but does only return Nsectors
! ! It is used in initialization to allocate arrays which require
! ! to know Nsectors before really initializing ED
! !+------------------------------------------------------------------+
! function get_Nsectors() result(Nsectors)
!   integer :: Nsectors
!   integer :: Ns
!   !
!   select case(bath_type)
!   case default
!      Ns = (Nbath+1)*Norb
!   case ('hybrid')
!      Ns = Nbath+Norb
!   case ('replica')
!      Ns = Norb*(Nbath+1)
!   end select
!   !
!   select case(ed_mode)
!   case default
!      Nsectors = (Ns+1)*(Ns+1) !nup=0:Ns;ndw=0:Ns
!   case ("superc")
!      Nsectors = 2*Ns+1     !sz=-Ns:Ns=2*Ns+1=Nlevels+1
!   case("nonsu2")
!      Nsectors = 2*Ns+1     !n=0:2*Ns=2*Ns+1=Nlevels+1
!   end select
! end function get_Nsectors
