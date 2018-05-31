!The twisted bilayer graphene flat bands effective model can be found in arXiv:1805.06819
! this driver needs the file 'eff_hopping.dat' to work
program ed_effective_tBLG
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  !
  real(8),parameter                             :: eV_to_meV=1000d0
  integer                                       :: i,rx,ry,iloop,Lk,Nso,Nvso,Nvalley
  logical                                       :: converged
  integer                                       :: ispin,ilat,Nhopping_wannier
  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath,Bath_prev
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Greal
  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk
  complex(8),allocatable,dimension(:,:)         :: graphHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  real(8),allocatable,dimension(:)              :: Wtk
  integer,allocatable,dimension(:)              :: ik2ix,ik2iy
  real(8),dimension(2)                          :: L1,L2,G1,G2
  real(8),dimension(2)                          :: bk1,bk2,bklen
  real(8),dimension(2)                          :: pointM,pointK1,pointK2
  !variables for the model:
  integer                                       :: Nk,Nkpath,m,n,p,q,unit
  real(8)                                       :: Lm,wmixing,re,im
  character(len=32)                             :: finput,fwannier
  !
  real(8),dimension(2)                          :: Eout
  real(8),allocatable,dimension(:)              :: dens
  !
  complex(8),dimension(:),allocatable           :: wannier_hopping ![Nhopping_wannier]
  integer,dimension(:,:),allocatable            :: wannier_orbital_index ![Nhopping_wannier,Norb]
  real(8),dimension(:,:),allocatable            :: wannier_Rgrid ![Nhopping_wannier,2=dim]
  !
  !Parse additional variables && read Input && read H(k)^2x2
  call parse_cmd_variable(finput,"FINPUT",default='input_tBLG_wannier.conf')
  call parse_input_variable(fwannier,"FWANNIER",finput,default='eff_hopping.data')
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  !
  call ed_read_input(trim(finput))
  !
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")


  if(Norb/=2.OR.Nspin/=1)stop "Wrong setup from input file: Norb=2,Nspin=1"
  Nvalley=2                        !decoupled Vallyes (+,-)
  Nso=Nspin*Norb
  Nvso=Nvalley*Nso
  !
  Lm=13.422225299852479
  L1=[sqrt(3.d0)/2.d0,1.d0/2.d0]*Lm
  L2=[0.d0,1.d0]*Lm
  G1=((2.d0*pi)/(L2(2)*L1(1)-L1(2)*L2(1)))*[L2(2),-L2(1)]
  G2=((2.d0*pi)/(L2(2)*L1(1)-L1(2)*L2(1)))*[-L1(2),L1(1)]

  !RECIPROCAL LATTICE VECTORS:
  bk1=G1
  bk2=G2

  call TB_set_bk(G1,G2)
  pointM =(G1+G2)/2.d0 
  pointK1=(G1+2.d0*G2)/3.d0 
  pointK2=(2.d0*G1+G2)/3.d0

  !read the hopping constants from file'eff_hopping.dat'!!!!!
  Nhopping_wannier= file_length(str(fwannier))
  open(free_unit(unit),file=str(fwannier))
  allocate(wannier_Rgrid(Nhopping_wannier,2))
  allocate(wannier_orbital_index(Nhopping_wannier,Norb))
  allocate(wannier_hopping(Nhopping_wannier))
  !
  do i=1,Nhopping_wannier !1162
     read(unit,*) rx,ry,p,q,re,im
     if(p>2.OR.q>2)stop "ed_TBG_effective_wannier ERROR: p OR q > 2 in wannier file"
     wannier_Rgrid(i,:)         = rx*L1 + ry*L2
     wannier_orbital_index(i,:) = [p,q]
     wannier_hopping(i)         = eV_to_meV*(re+xi*im) !transforms to meV
  end do
  close(unit)


  !Allocate Weiss Field:
  allocate(Weiss(Nvalley,Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  allocate(Smats(Nvalley,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Gmats(Nvalley,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Sreal(Nvalley,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Greal(Nvalley,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
  allocate(Hloc(Nvalley,Nspin,Nspin,Norb,Norb));Hloc=zero

  !Build the Hamiltonian on a grid or on  path
  call build_hk()
  Hloc = lso2nnn_reshape(graphHloc,Nvalley,Nspin,Norb)

  
  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath(Nvalley,Nb))
  allocate(Bath_prev(Nvalley,Nb))
  call ed_init_solver(Bath,Hloc)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(Bath,Hloc)
     call ed_get_sigma_matsubara(Smats,Nvalley)

     ! compute the local gf:
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats)

     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats,Smats,Weiss,Hloc)
     else
        call dmft_delta(Gmats,Smats,Weiss,Hloc)
     endif

     !Fit the new bath, starting from the old bath + the supplied Weiss
     call ed_chi2_fitgf(Bath,Weiss,Hloc,ispin=1)

     !MIXING:
     if(iloop>1)Bath=wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     converged = check_convergence(Weiss(:,1,1,1,1,:),dmft_error,nsuccess,nloop)

     call end_loop
  enddo

  call dmft_print_gf_matsubara(Gmats,"Gmats",iprint=1)

  ! extract and print retarded self-energy and Green's function 
  call ed_get_sigma_real(Sreal,Nvalley)
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)


  call dmft_kinetic_energy(Hk,Wtk,Smats)



contains



  !--------------------------------------------------------------------!
  ! Twisted Bilayer Graphene (theta=1.05) effective  HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_TBGeffective_model(kpoint,Nvso) result(hk)
    real(8),dimension(:)                :: kpoint
    integer                             :: Nvso,i,m,n,p,q
    complex(8)                          :: uni
    real(8)                             :: k_dot_R
    complex(8),dimension(Nvso,Nvso)     :: hk       
    real(8),dimension(2)                :: RR
    !
    hk=zero
    do i=1,Nhopping_wannier     !1162
       p=wannier_orbital_index(i,1)
       q=wannier_orbital_index(i,2)
       !
       k_dot_R = dot_product(kpoint,wannier_Rgrid(i,:))!RR)        ! RR = m*L1 + n*L2
       !
       hk(p,q)    = hk(p,q)     + exp( xi*k_dot_R)*wannier_hopping(i)
       hk(p+2,q+2)= hk(p+2,q+2) + exp(-xi*k_dot_R)*wannier_hopping(i)
       !
    end do
    !
  end function hk_TBGeffective_model






  !---------------------------------------------------------------------
  !PURPOSE: Get BL Graphene Model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional                              :: file
    integer                                                :: i,j,ik
    integer                                                :: ix,iy
    real(8)                                                :: kx,ky  
    integer                                                :: iorb,jorb
    integer                                                :: isporb,jsporb
    integer                                                :: ispin,jspin
    integer                                                :: unit
    complex(8),dimension(Nvalley,Nspin,Nspin,Norb,Norb,Lmats) :: Gmats
    real(8),dimension(2)                                   :: kvec
    real(8)                                                :: blen,area_hex,area_rect,points_in,points_tot
    real(8),allocatable,dimension(:)                       :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable                     :: KPath

    Lk= Nk*Nk

    write(LOGfile,*)"Build effective H(k) of twisted Bilayer Graphene (theta=1.05):",Lk
    write(LOGfile,*)"# of SO-bands     :",Nvso

    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nvso,Nvso,Lk));Hk=zero
    allocate(wtk(Lk));Wtk=0d0
    call TB_build_model(Hk,hk_TBGeffective_model,Nvso,[Nk,Nk])
    Wtk = 1d0/Lk


    if(present(file))then
       call TB_write_hk(Hk,"Hkr_tBLG_wannier.data",&
            No=Nvso,&
            Nd=Norb,&
            Np=0,&
            Nineq=1,&
            Nkvec=[Nk,Nk])
    endif
    !
    allocate(graphHloc(Nvso,Nvso))
    graphHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(graphHloc))<1.d-4)graphHloc=0d0
    call TB_write_Hloc(graphHloc)
    call TB_write_Hloc(graphHloc,'Hloc_tBLG_wannier.dat')
    !
    !
    allocate(Kpath(4,2))
    KPath(1,:)=pointK1
    KPath(2,:)=[0,0]
    KPath(3,:)=pointM
    KPath(4,:)=pointK2
    call TB_Solve_model(hk_TBGeffective_model,Nvso,KPath,Nkpath,&
         colors_name=[red1,blue1,green1,purple1],&
         points_name=[character(len=10) :: "K1","G","M","K2"],&
         file="Eigenbands_tBLG_wannier.nint")

    !Build the local GF:
    call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nvalley,Nspin,Nspin,Norb,Norb,Lmats))
    call dmft_print_gf_matsubara(Gmats,"LG0",iprint=1)
    !
  end subroutine build_hk








end program ed_effective_tBLG



