program ed_SOC_bethe
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none


  !#########   VARIABLES DECLARATION   #########
  integer                                        :: iloop,unit
  integer                                        :: Nlat,ilat
  integer                                        :: Nso,io,jo
  integer                                        :: iorb,jorb,ispin,jspin
  logical                                        :: converged,converged_n
  real(8)                                        :: wmixing
  character(len=60)                              :: finput
  character(len=32)                              :: hkfile
  !Mpi:
  integer                                        :: comm,rank,ier
  logical                                        :: master
  !Bath:
  integer                                        :: Nb
  real(8),dimension(:),allocatable               :: Bath
  !Local functions:
  complex(8),dimension(:,:,:,:,:),allocatable    :: Smats,Gmats
  complex(8),dimension(:,:,:,:,:),allocatable    :: Sreal,Greal
  complex(8),dimension(:),allocatable            :: Ftest
  !Weiss&Hybridization functions:
  complex(8),dimension(:,:,:,:,:),allocatable    :: Field,Field_old
  !Frequency meshes:
  integer                                        :: iw
  real(8)                                        :: dw
  real(8),dimension(:),allocatable               :: wr,wm
  !Energy meshes and Bethe DoS:
  integer                                        :: ie,Le
  real(8),dimension(4)                           :: semiW          !semi-bandwidth from inputfile
  real(8),dimension(:),allocatable               :: Dband          !semi-bandwidth used
  real(8),dimension(:),allocatable               :: de
  real(8),dimension(:,:),allocatable             :: DoS,er
  !K-dependent Hamiltonian:
  integer                                        :: ik,Nk,Nkpath,Lk
  real(8),dimension(:),allocatable               :: Wtk
  complex(8),dimension(:,:,:),allocatable        :: Hk
  !Local Hamiltonian:
  complex(8),dimension(:,:),allocatable          :: Hloc_so
  complex(8),dimension(:,:,:,:),allocatable      :: Hloc_nn
  !Variables for the model:
  real(8)                                        :: lambda_soc
  real(8),dimension(4)                           :: Elevels        !Local energies from inputfile
  real(8),dimension(:),allocatable               :: Eloc           !Local energies used
  character(len=32)                              :: lattice
  !Observables
  real(8),dimension(:),allocatable               :: dens
  complex(8),dimension(:,:),allocatable          :: rho
  complex(8),dimension(:,:),allocatable          :: Op
  complex(8),dimension(:,:),allocatable          :: Gtmp
  !Miscellaneous:
  logical                                        :: socsym,calcG0
  logical                                        :: symtest
  !Dummy variables:
  real(8)                                        :: dum,cum
  complex(8),dimension(:,:),allocatable          :: Odum


  !#########   MPI INITIALIZATION   #########
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  !rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  !#########    VARIABLE PARSING    #########
  call parse_cmd_variable(finput,      "FINPUT",  default='inputED_S.in')
  !
  call parse_input_variable(nk,        "NK",       finput,default=10)
  call parse_input_variable(nkpath,    "NKPATH",   finput,default=500)
  call parse_input_variable(Le,        "LE",       finput,default=1000)
  call parse_input_variable(semiW,     "SEMIW",    finput,default=[1d0,1d0,1d0,1d0,1d0])
  call parse_input_variable(Elevels,   "ELEVELS",  finput,default=[0d0,0d0,0d0,0d0,0d0])
  call parse_input_variable(wmixing,   "WMIXING",  finput,default=0.5d0)
  call parse_input_variable(lambda_soc,"SOC",      finput,default=0.0d0)
  call parse_input_variable(socsym,    "SOCSYM",   finput,default=.false.)
  call parse_input_variable(symtest,   "SYMTEST",  finput,default=.true.)
  call parse_input_variable(calcG0,    "CALCG0LOC",finput,default=.true.)
  call parse_input_variable(lattice,   "LATTICE",  finput,default="Bethe")
  call parse_input_variable(Utensor,   "UTENSOR",  finput,default=.false.)
  !
  call ed_read_input(trim(finput),comm)
  lattice=reg(lattice)
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,     "NORB")
  call add_ctrl_var(Nspin,    "NSPIN")
  call add_ctrl_var(beta,     "BETA")
  call add_ctrl_var(xmu,      "XMU")
  call add_ctrl_var(wini,     "WINI")
  call add_ctrl_var(wfin,     "WFIN")
  call add_ctrl_var(eps,      "eps")
  call add_ctrl_var(ed_para,  "ED_PARA")
  call add_ctrl_var(bath_type,"BATH_TYPE")
  call add_ctrl_var(Jz_basis, "JZ_BASIS")
  call add_ctrl_var(replica_operators, "REPLICA_OPERATORS")


  !#########      CONTROL FLAGS     #########
  dum=0d0
  do iorb=1,Norb
     dum=dum+abs(Elevels(iorb))
  enddo
  if ((dum.ne.0d0).and.(lambda_soc.ne.0d0)) stop "Crystal field and SOC are not compatible"


  !#########   FIELDS ALLOCATION    #########
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats)    ); Smats=zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats)    ); Gmats=zero
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal)    ); Sreal=zero
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal)    ); Greal=zero
  allocate(Field(Nspin,Nspin,Norb,Norb,Lmats)    ); Field=zero
  !
  allocate(Field_old(Nspin,Nspin,Norb,Norb,Lmats)); Field_old=zero
  allocate(Ftest(Lmats)); Ftest=zero
  allocate(dens(Norb)); dens=0d0
  allocate(Op(Nspin*Norb,Nspin*Norb)); Op=zero
  allocate(rho(Nspin*Norb,Nspin*Norb)); rho=zero
  allocate(Gtmp(Nspin*Norb,Nspin*Norb)); Gtmp=zero
  allocate(wr(Lreal)); wr=linspace(wini,wfin,Lreal,mesh=dw)
  allocate(wm(Lmats)); wm = pi/beta*dble(2*arange(1,Lmats)-1)
  !
  allocate(Odum(Norb,Norb)); Odum=zero


  !######### LATTICE INITIALIZATION #########
  Nso=Nspin*Norb
  !
  allocate(Eloc(Nso))
  Eloc(:Norb) = Elevels(:Norb)
  Eloc(Norb+1:) = Elevels(:Norb)
  allocate(Dband(Nso))
  Dband(:Norb)= semiW(:Norb)
  Dband(Norb+1:)= semiW(:Norb)
  !
  allocate(Hloc_so(Nso,Nso)              ); Hloc_so=zero
  allocate(Hloc_nn(Nspin,Nspin,Norb,Norb)); Hloc_nn=zero
  !
  if (lattice.eq."Bethe") then
     !
     allocate(er(Nso,Le)); er=0d0
     allocate(DoS(Nso,Le)); DoS=0d0
     allocate(de(Nso)); de=0d0
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb
           er(io,:) = linspace(-Dband(io),Dband(io),Le,mesh=de(io))! + Eloc(io)
           DoS(io,:) = dens_bethe(er(io,:),Dband(io))*de(io)
        enddo
     enddo
     !
     Hloc_so = diag(Eloc) + lambda_soc * atomic_SOC()
     if (lambda_soc.ne.0d0) call Cbasis_to_Ybasis(Hloc_so,"Hloc_Bethe")
     !
  else
     !
     call build_hk()
     !
  endif
  !
  Hloc_nn = so2nn_reshape(Hloc_so,Nspin,Norb)
  if (master) then
     call TB_write_Hloc(Hloc_so,"Hloc")
     if(lambda_soc.ne.0d0) then
        call Ybasis_to_Jbasis(Hloc_so) ; call TB_write_Hloc(Hloc_so,"Hloc_J")
        call Jbasis_to_Cbasis(Hloc_so) ; call TB_write_Hloc(Hloc_so,"Hloc_C")
        call Cbasis_to_Ybasis(Hloc_so) ; call TB_write_Hloc(Hloc_so,"Hloc_Y")
     endif
  endif
  !
  if (calcG0) then
     !
     if (lattice.eq."Bethe") then
        call dmft_gloc_matsubara(Comm,er,DoS,Hloc_so,Gmats,Smats)               !provides only diagonal one
        Gmats(1,2,1,3,:)=-(Gmats(1,1,2,2,:)-Gmats(1,1,1,1,:))/sqrt(2.d0) ; Gmats(2,1,3,1,:)=Gmats(1,2,1,3,:)
        Gmats(1,2,3,2,:)=+(Gmats(1,1,2,2,:)-Gmats(1,1,1,1,:))/sqrt(2.d0) ; Gmats(2,1,2,3,:)=Gmats(1,2,3,2,:)
        call dmft_gloc_realaxis(Comm,er,DoS,Hloc_so,Greal,Sreal)                !provides only diagonal one
        Greal(1,2,1,3,:)=-(Greal(1,1,2,2,:)-Greal(1,1,1,1,:))/sqrt(2.d0) ; Greal(2,1,3,1,:)=Greal(1,2,1,3,:)
        Greal(1,2,3,2,:)=+(Greal(1,1,2,2,:)-Greal(1,1,1,1,:))/sqrt(2.d0) ; Greal(2,1,2,3,:)=Greal(1,2,3,2,:)
     else
        call dmft_gloc_matsubara(Comm,Hk,Wtk,Gmats,Smats)
        call dmft_gloc_realaxis(Comm,Hk,Wtk,Greal,Sreal)
     endif
     !
     if(lambda_soc.ne.0d0) then
        call print_G("0Y")
        do iw=1,Lmats
           Gtmp=nn2so_reshape(Gmats(:,:,:,:,iw),Nspin,Norb)
           call Ybasis_to_Jbasis(Gtmp,"Gmats")
           Gmats(:,:,:,:,iw)=so2nn_reshape(Gtmp,Nspin,Norb)
        enddo
        do iw=1,Lreal
           Gtmp=nn2so_reshape(Greal(:,:,:,:,iw),Nspin,Norb)
           call Ybasis_to_Jbasis(Gtmp,"Greal")
           Greal(:,:,:,:,iw)=so2nn_reshape(Gtmp,Nspin,Norb)
        enddo
        call print_G("0J")
     else
        call print_G("0J")
     endif
     !
  endif
  !



  !######### SOLVER INITIALIZATION  #########
  if (lambda_soc.eq.0d0) then
     Jz_basis = .false.
     bath_type = "normal"
     Nb=get_bath_dimension()
  else
     Jz_basis = .true.
     bath_type = "replica"
     ed_mode = "nonsu2"
     Nb=get_bath_dimension(Hloc_nn)
     !
     if (replica_operators) then
        Op = eye(Nso)
        call set_replica_operators(Op,1)
        Op=atomic_SOC()
        call Cbasis_to_Ybasis(Op,"SOC")
        call set_replica_operators(Op,2)
        call set_replica_operators(-1.d0*Op,3)
     endif
     !
  endif
  if(master)write(LOGfile,*)"Bath_size:",Nb
  allocate(Bath(Nb));Bath=0.0d0
  !
  call ed_init_solver(Comm,Bath,Hloc_nn)
  !
  if (lambda_soc.ne.0d0.and.Utensor) call ed_rotate_interaction(Comm,orbital_Lz_rotation_Norb())
  !
  Odum=zero
  dum=cos(pi/5.d0)
  cum=sin(pi/5.d0)
  Odum(1,1)=cmplx(+dum,0d0)
  Odum(1,2)=cmplx(-cum,0d0)
  Odum(2,1)=cmplx(+cum,0d0)
  Odum(2,2)=cmplx(+dum,0d0)
  Odum(3,3)=cmplx(1.d0,0d0)
  !call ed_rotate_interaction(Comm,orbital_Lz_rotation_Norb())
  !

  !#########       DMFT CYCLE       #########
  iloop=0 ; converged=.false. ; converged_n=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if (master) call start_loop(iloop,nloop,"DMFT-loop")


     !solve impurity
     call ed_solve(comm,Bath)


     !get sigmas
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     if ((lambda_soc.ne.0d0).and.socsym) then
        call SOC_symmetrize(Smats,"Y",ed_para)
        call SOC_symmetrize(Sreal,"Y",ed_para)
     endif


     !get Gloc
     if (lattice.eq."Bethe") then
        !these two are providing only diagonal entries but that's only for check.
        call dmft_gloc_matsubara(Comm,er,DoS,Hloc_so,Gmats,Smats)
        call dmft_gloc_realaxis(Comm,er,DoS,Hloc_so,Greal,Sreal)
     else
        call dmft_gloc_matsubara(Comm,Hk,Wtk,Gmats,Smats)
        call dmft_gloc_realaxis(Comm,Hk,Wtk,Greal,Sreal)
     endif


     !print everything to check for the correct symmetry in {Y} basis
     call print_G("Y")
     call print_Sigma("Y")
     if(lambda_soc.ne.0d0) then
        do iw=1,Lreal
           Gtmp=nn2so_reshape(Greal(:,:,:,:,iw),Nspin,Norb)
           call Ybasis_to_Jbasis(Gtmp)
           Greal(:,:,:,:,iw)=so2nn_reshape(Gtmp,Nspin,Norb)
        enddo
        call print_G("J")
     endif


     !Get the Weiss field/Delta function to be fitted
     if (lattice.eq."Bethe") then
        call ed_get_gimp_matsubara(Gmats)
        call dmft_self_consistency(Comm,Gmats,Field,Hloc_nn,SCtype=cg_scheme,wbands=Dband)
     else
        call dmft_self_consistency(Comm,Gmats,Smats,Field,Hloc_nn,SCtype=cg_scheme)
     endif
     if (master) then
        if (symtest) then
           call dmft_print_gf_matsubara(Field,cg_scheme,iprint=3)
        else
           call dmft_print_gf_matsubara(Field,cg_scheme,iprint=1)
        endif
     endif


     ! Mixing
     if (iloop>1) Field = wmixing*Field + (1.d0-Field)*Field_old
     Field_old=Field


     !Perform the SELF-CONSISTENCY by fitting the new bath
     if (bath_type=="normal") then
        if (ed_para) then
           call ed_chi2_fitgf(Comm,Field,bath,ispin=1)
           call spin_symmetrize_bath(bath,save=.true.)
           call orb_symmetrize_bath(bath,save=.true.)
        else
           call ed_chi2_fitgf(Comm,Field,bath)
        endif
     else
        call ed_chi2_fitgf(Comm,Field,bath)
     endif


     !Get observables
     if (master) then
        call ed_get_density_matrix(rho)
        if (lambda_soc.ne.0d0) then
           call Ybasis_to_Cbasis(rho)
           call ed_print_density_matrix(rho,suff="C")
           call Cbasis_to_Jbasis(rho)
           call ed_print_density_matrix(rho,suff="J")
        endif
        !
       !call ed_get_quantum_SOC_operators_single()
    endif


     !Check convergence (if required change chemical potential)
     if (master) then
        if (bath_type=="normal") then
           do ispin=1,Nspin
              do iorb=1,Norb
                 Ftest=Ftest+Field(ispin,ispin,iorb,iorb,:)/Nso
              enddo
           enddo
        else
           Ftest=(Field(1,1,3,3,:)+Field(1,2,3,2,:))/2.d0
        endif
        converged = check_convergence(Ftest,dmft_error,nsuccess,nloop,reset=.false.)
        if (nread/=0d0) then
           call ed_get_dens(dens)
           !call ed_search_variable(xmu,sum(dens),converged)
           call search_chempot(xmu,sum(dens),converged_n)
        endif
        converged = converged .and. converged_n
     endif
     call Bcast_MPI(comm,converged)
     call Bcast_MPI(comm,xmu)
     !
     if(master)call end_loop
     !
  enddo




contains



   !+-------------------------------------------------------------------------+!
   !PURPOSE:
   !+-------------------------------------------------------------------------+!
   subroutine build_hk()
     implicit none
     !internal
     real(8),dimension(3)                        :: bk_x,bk_y,bk_z
     !
     bk_x = [1.d0,0.d0,0.d0]*2*pi
     bk_y = [0.d0,1.d0,0.d0]*2*pi
     bk_z = [0.d0,0.d0,1.d0]*2*pi
     !
     if (lattice.eq."Square") then
        call TB_set_bk(bk_x,bk_y)
        Lk=Nk*Nk
        allocate(Hk(Nso,Nso,Lk));Hk=zero
        allocate(Wtk(Lk));Wtk=1.d0/Lk
        call TB_build_model(Hk,hk_t2g_Sz,Nso,[Nk,Nk])
        if(master)call TB_write_hk(Hk,"Hk.dat",Nso,Norb,1,1,[Nk,Nk])
     elseif (lattice.eq."Cubic") then
        call TB_set_bk(bk_x,bk_y,bk_z)
        Lk=Nk*Nk*Nk
        allocate(Hk(Nso,Nso,Lk));Hk=zero
        allocate(Wtk(Lk));Wtk=1.d0/Lk
        call TB_build_model(Hk,hk_t2g_Sz,Nso,[Nk,Nk,Nk])
        if(master)call TB_write_hk(Hk,"Hk.dat",Nso,Norb,1,1,[Nk,Nk,Nk])
     endif
     !
     Hloc_so = sum(Hk(:,:,:),dim=3)/Lk
     !
     if (lambda_soc.ne.0d0) then
        do ik=1,Lk
           call Cbasis_to_Ybasis(Hk(:,:,ik),"Hk_"//str(ik))
        enddo
        call Cbasis_to_Ybasis(Hloc_so,"Hloc_"//lattice)
     endif
     !
   end subroutine build_hk



   !+-------------------------------------------------------------------------+!
   !PURPOSE:
   !+-------------------------------------------------------------------------+!
   function hk_t2g_Sz(kvec,Nso_) result(H_k)
     implicit none
     integer                                     :: Nso_
     real(8),dimension(:)                        :: kvec
     complex(8),dimension(Nso_,Nso_)             :: H_k
     !internal
     real(8)                                     :: kx,ky,kz
     real(8)                                     :: t, fact
     !
     if (lattice.eq."Square") then
        kx=kvec(1)
        ky=kvec(2)
        kz=pi/2.d0
        if (cos(kz).gt.1e-8) stop "square cos(kz) not vanishing"
        fact=4.d0
     elseif (lattice.eq."Cubic") then
        kx=kvec(1)
        ky=kvec(2)
        kz=kvec(3)
        fact=6.d0
     endif
     !
     H_k=zero
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb
           !
           t = Dband(io)/fact
           H_k(io,io) = -2 * t * (cos(kx) + cos(ky) + cos(kz)) + Eloc(io)
           !
        enddo
     enddo
     !
     H_k = H_k + lambda_soc * atomic_SOC()
     !
   end function hk_t2g_Sz



   !+-------------------------------------------------------------------------+!
   !PURPOSE:
   !+-------------------------------------------------------------------------+!
   subroutine Cbasis_to_Ybasis(A,varname)
     implicit none
     complex(8),dimension(:,:),intent(inout)     :: A
     character(len=*),intent(in),optional        :: varname
     !iternal
     complex(8),dimension(size(A,1),size(A,2))   :: B
     complex(8),dimension(6,6)                   :: U,Udag
     character(len=32)                           :: assertname
     !
     assertname="--"
     if(present(varname))assertname=varname
     call assert_shape(A,[Nspin*Norb,Nspin*Norb],"Cbasis_to_Ybasis",reg(assertname))
     !
     B=A
     A=zero
     !
     U=orbital_Lz_rotation_NorbNspin()
     Udag=transpose(conjg(U))
     !
     A=matmul(Udag,matmul(B,U))
     !
     do io=1,Nso
        do jo=1,Nso
           if (abs(A(io,jo))<1e-4) A(io,jo)=zero
        enddo
     enddo
     !
   end subroutine Cbasis_to_Ybasis
   !
   subroutine Ybasis_to_Cbasis(A,varname)
     implicit none
     complex(8),dimension(:,:),intent(inout)     :: A
     character(len=*),intent(in),optional        :: varname
     !iternal
     complex(8),dimension(size(A,1),size(A,2))   :: B
     complex(8),dimension(6,6)                   :: U,Udag
     character(len=32)                           :: assertname
     !
     assertname="--"
     if(present(varname))assertname=varname
     call assert_shape(A,[Nspin*Norb,Nspin*Norb],"Ybasis_to_Cbasis",reg(assertname))
     !
     B=A
     A=zero
     !
     U=transpose(conjg(orbital_Lz_rotation_NorbNspin()))
     Udag=transpose(conjg(U))
     !
     A=matmul(Udag,matmul(B,U))
     !
     do io=1,Nso
        do jo=1,Nso
           if (abs(A(io,jo))<1e-4) A(io,jo)=zero
        enddo
     enddo
     !
  end subroutine Ybasis_to_Cbasis
  !
  subroutine Cbasis_to_Jbasis(A,varname)
    implicit none
    complex(8),dimension(:,:),intent(inout)     :: A
    character(len=*),intent(in),optional        :: varname
    !iternal
    complex(8),dimension(size(A,1),size(A,2))   :: B
    complex(8),dimension(6,6)                   :: U,Udag
    character(len=32)                           :: assertname
    !
    assertname="--"
    if(present(varname))assertname=varname
    call assert_shape(A,[Nspin*Norb,Nspin*Norb],"Cbasis_to_Jbasis",reg(assertname))
    !
    B=A
    A=zero
    !
    U=atomic_SOC_rotation()
    Udag=transpose(conjg(U))
    !
    A=matmul(Udag,matmul(B,U))
    !
    do io=1,Nso
      do jo=1,Nso
          if (abs(A(io,jo))<1e-4) A(io,jo)=zero
      enddo
    enddo
    !
  end subroutine Cbasis_to_Jbasis
  !
  subroutine Jbasis_to_Cbasis(A,varname)
    implicit none
    complex(8),dimension(:,:),intent(inout)     :: A
    character(len=*),intent(in),optional        :: varname
    !iternal
    complex(8),dimension(size(A,1),size(A,2))   :: B
    complex(8),dimension(6,6)                   :: U,Udag
    character(len=32)                           :: assertname
    !
    assertname="--"
    if(present(varname))assertname=varname
    call assert_shape(A,[Nspin*Norb,Nspin*Norb],"Jbasis_to_Cbasis",reg(assertname))
    !
    B=A
    A=zero
    !
    U=transpose(conjg(atomic_SOC_rotation()))
    Udag=transpose(conjg(U))
    !
    A=matmul(Udag,matmul(B,U))
    !
    do io=1,Nso
      do jo=1,Nso
          if (abs(A(io,jo))<1e-4) A(io,jo)=zero
      enddo
    enddo
    !
  end subroutine Jbasis_to_Cbasis
  !
  subroutine Ybasis_to_Jbasis(A,varname)
    implicit none
    complex(8),dimension(:,:),intent(inout)     :: A
    character(len=*),intent(in),optional        :: varname
    !iternal
    complex(8),dimension(size(A,1),size(A,2))   :: B
    complex(8),dimension(6,6)                   :: U,Udag
    character(len=32)                           :: assertname
    !
    assertname="--"
    if(present(varname))assertname=varname
    call assert_shape(A,[Nspin*Norb,Nspin*Norb],"Cbasis_to_Jbasis",reg(assertname))
    !
    B=A
    A=zero
    !
    U=matmul(transpose(conjg(orbital_Lz_rotation_NorbNspin())),atomic_SOC_rotation())
    Udag=transpose(conjg(U))
    !
    A=matmul(Udag,matmul(B,U))
    !
    do io=1,Nso
      do jo=1,Nso
          if (abs(A(io,jo))<1e-4) A(io,jo)=zero
      enddo
    enddo
    !
  end subroutine Ybasis_to_Jbasis
  !
  subroutine Jbasis_to_Ybasis(A,varname)
    implicit none
    complex(8),dimension(:,:),intent(inout)     :: A
    character(len=*),intent(in),optional        :: varname
    !iternal
    complex(8),dimension(size(A,1),size(A,2))   :: B
    complex(8),dimension(6,6)                   :: U,Udag
    character(len=32)                           :: assertname
    !
    assertname="--"
    if(present(varname))assertname=varname
    call assert_shape(A,[Nspin*Norb,Nspin*Norb],"Cbasis_to_Jbasis",reg(assertname))
    !
    B=A
    A=zero
    !
    U=matmul(transpose(conjg(atomic_SOC_rotation())),orbital_Lz_rotation_NorbNspin())
    Udag=transpose(conjg(U))
    !
    A=matmul(Udag,matmul(B,U))
    !
    do io=1,Nso
      do jo=1,Nso
          if (abs(A(io,jo))<1e-4) A(io,jo)=zero
      enddo
    enddo
    !
  end subroutine Jbasis_to_Ybasis



   !+-------------------------------------------------------------------------+!
   !PURPOSE:
   !+-------------------------------------------------------------------------+!
   subroutine print_G(suff)
     implicit none
     character(len=*),optional,intent(in)        :: suff
     character(len=32)                           :: Gm,Gr
     !
     if (present(suff)) then
        Gm=reg("G"//suff//"mats")
        Gr=reg("G"//suff//"real")
     else
        Gm=reg("Gmats")
        Gr=reg("Greal")
     endif
     !
     if (master) then
        if (symtest) then
           call dmft_print_gf_matsubara(Gmats,Gm,iprint=3)
           call dmft_print_gf_realaxis(Greal,Gr,iprint=3)
        else
           call dmft_print_gf_matsubara(Gmats,Gm,iprint=1)
           call dmft_print_gf_realaxis(Greal,Gr,iprint=1)
        endif
     endif
     !
   end subroutine print_G

   subroutine print_Sigma(suff)
     implicit none
     character(len=*),optional,intent(in)        :: suff
     character(len=32)                           :: Sm,Sr
     !
     if (present(suff)) then
        Sm=reg("S"//suff//"mats")
        Sr=reg("S"//suff//"real")
     else
        Sm=reg("Smats")
        Sr=reg("Sreal")
     endif
     !
     if (master) then
        if (symtest) then
           call dmft_print_gf_matsubara(Smats,Sm,iprint=3)
           call dmft_print_gf_realaxis(Sreal,Sr,iprint=3)
        else
           call dmft_print_gf_matsubara(Smats,Sm,iprint=1)
           call dmft_print_gf_realaxis(Sreal,Sr,iprint=1)
        endif
     endif
     !
   end subroutine print_Sigma



   !+-------------------------------------------------------------------------+!
   !PURPOSE: Inversion test
   !+-------------------------------------------------------------------------+!
   subroutine inversion_test(A,B,tol)
     implicit none
     complex (8), intent(in)                     ::   A(Nspin*Norb,Nspin*Norb)
     complex (8), intent(in)                     ::   B(Nspin*Norb,Nspin*Norb)
     real    (8), intent(in)                     ::   tol
     integer                                     ::   dime
     if (size(A).ne.size(B)) then
        write(LOGfile,*) "Matrices not equal cannot perform inversion test"
        stop
     endif
     dime=maxval(shape(A))
     if (abs(float(dime)-real(sum(matmul(A,B)))).gt.tol) then
        write(LOGfile,'(A30)') "inversion test fail"
     endif
   end subroutine inversion_test



end program ed_SOC_bethe
