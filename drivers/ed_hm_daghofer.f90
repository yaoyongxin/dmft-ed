program ed_hm_3bands
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                                       :: iloop,Lk,Nso,Nlso,Nlat
  logical                                       :: converged
  integer                                       :: ispin,ilat!,i,j

  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath,Bath_prev

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Greal

  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk
  complex(8),allocatable,dimension(:,:)         :: modelHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  real(8),allocatable,dimension(:)              :: Wtk

  integer,allocatable,dimension(:)              :: ik2ix,ik2iy
  real(8),dimension(2)                          :: bk1,bk2

  !variables for the model:
  integer                                       :: Nk,Nkpath
  real(8)                                       :: ts,wmixing
  character(len=32)                             :: finput
  character(len=32)                             :: hkfile
  logical                                       :: spinsym,bathsym
  !
  real(8),dimension(2)                          :: Eout
  real(8),allocatable,dimension(:)              :: dens
  !
  real(8),dimension(:,:),allocatable            :: Zmats
  complex(8),dimension(:,:,:),allocatable       :: Zfoo
  complex(8),allocatable,dimension(:,:,:,:,:)   :: S0

  !modify daghofer hamiltonian
  real(8)                                       :: alpha,theta

  !parse input file
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')

  call parse_input_variable(alpha,"ALPHA",finput,default=1.d0)
  call parse_input_variable(theta,"THETA",finput,default=0.d0)

  !Parse additional variables && read Input && read H(k)^2x2
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS","inputED.conf",default=1d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  !
  call parse_input_variable(bathsym,"BATHSYM",finput,default=.false.)
  !
  call ed_read_input(trim(finput))

  if(Norb/=3)stop "Wrong setup from input file: Norb=3"
  Nlat=1
  Nso=Nspin*Norb
  Nlso=Nlat*Nso



  !Allocate Weiss Field:
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb));Hloc=zero
  allocate(S0(Nlat,Nspin,Nspin,Norb,Norb));S0=zero
  allocate(Zmats(Nlso,Nlso));Zmats=eye(Nlso)
  allocate(Zfoo(Nlat,Nso,Nso));Zfoo=0d0
  allocate(dens(Norb));dens=0d0

  !Build the Hamiltonian on a grid or on a path
  call build_hk(trim(hkfile))
  Hloc = lso2nnn_reshape(modelHloc,Nlat,Nspin,Norb)

  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath(Nlat,Nb))
  allocate(Bath_prev(Nlat,Nb))
  call ed_init_solver(Bath,Hloc)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(Bath,Hloc,iprint=0)

     call ed_get_sigma_matsubara(Smats,Nlat)
     call ed_get_dens(dens)

     ! compute the local gf:
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats,iprint=4)

     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats,Smats,Weiss,Hloc,iprint=0)
     else
        call dmft_delta(Gmats,Smats,Weiss,Hloc,iprint=0)
     endif

     !Fit the new bath, starting from the old bath + the supplied Weiss
     call ed_chi2_fitgf(Bath,Weiss,Hloc,ispin=1)
     
     !copy iorb component in jorb to enforce symmetry
     if(bathsym) call copy_component_bath(Bath(1,:),1,1,Bath(1,:),1,2)

     !MIXING:
     if(iloop>1)Bath=wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath

     converged = check_convergence(Weiss(:,1,1,1,1,:),dmft_error,nsuccess,nloop)

     if(NREAD/=0.d0) call search_chemical_potential(xmu,sum(dens),converged)


     call end_loop
  enddo

  call ed_get_sigma_real(Sreal,Nlat)
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal,iprint=4)


  ! save self-energy on disk
  call save_sigma_mats(Smats)
  call save_sigma_real(Sreal)

contains




  !--------------------------------------------------------------------!
  !Lattice Hamitonian:
  !--------------------------------------------------------------------!
  function hk_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: kx,ky
    real(8)                         :: t1,t2,t3,t4,t5,t6,t7,t8,dxy,xmu_tb
    !
    kx=kpoint(1)
    ky=kpoint(2)
    !
    hk(:,:) = zero
    !
    !daghofer model 
    t1  =  0.02d0
    t2  =  0.06d0
    t3  =  0.03d0
    t4  = -0.01d0
    t5  =  0.2d0*alpha
    t6  =  0.3d0*alpha
    t7  = -0.2d0*alpha
    t8  = -t7/2.d0
    dxy =  0.4d0-theta
    !
    xmu_tb = 0.212d0
    !
    hk(1,1) = 2.d0*t2*cos(kx) + 2.d0*t1*cos(ky) + 4.d0*t3*cos(kx)*cos(ky) - xmu_tb
    hk(2,2) = 2.d0*t1*cos(kx) + 2.d0*t2*cos(ky) + 4.d0*t3*cos(kx)*cos(ky) - xmu_tb
    hk(3,3) = 2.d0*t5*(cos(kx)+cos(ky)) + 4.d0*t6*cos(kx)*cos(ky) + dxy   - xmu_tb
    hk(1,2) = 4.d0*t4*sin(kx)*sin(ky)
    hk(1,3) = 2.d0*t7*sin(kx)*xi + 4.d0*t8*sin(kx)*cos(ky)*xi
    hk(2,3) = 2.d0*t7*sin(ky)*xi + 4.d0*t8*sin(ky)*cos(kx)*xi 
    !
    hk(2,1) = hk(1,2)
    hk(3,1) = dconjg(hk(1,3))
    hk(3,2) = dconjg(hk(2,3))
    !
  end function hk_model






  !---------------------------------------------------------------------
  !PURPOSE: get model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional             :: file
    integer                               :: i,j,ik
    integer                               :: ix,iy
    real(8)                               :: kx,ky  
    integer                               :: iorb,jorb
    integer                               :: isporb,jsporb
    integer                               :: ispin,jspin
    integer                               :: unit
    complex(8),dimension(Nlso,Nlso,Lmats) :: Gmats,fooSmats
    complex(8),dimension(Nlso,Nlso,Lreal) :: Greal,fooSreal
    real(8),dimension(2)                  :: kvec
    real(8)                               :: blen,area_hex,area_rect,points_in,points_tot
    real(8),allocatable,dimension(:)      :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable    :: kpath

    Lk= Nk*Nk

    write(LOGfile,*)"Build H(k)    :",Lk
    write(LOGfile,*)"# of SO-bands :",Nlso

    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nlso,Nlso,Lk));Hk=zero
    allocate(wtk(Lk));Wtk=0d0

    call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])

    call TB_build_model(Hk, hk_model, Nlso, [Nk,Nk])

    Wtk = 1d0/Lk

    if(present(file))&
         call TB_write_hk(Hk, trim(file), &
         No = Nlso,Nd = Norb,Np = 0,Nineq = 1,&
         Nkvec=[Nk,Nk])
    !
    allocate(modelHloc(Nlso,Nlso))
    modelHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(modelHloc))<1.d-4)modelHloc=0d0

    !path: G X M G
    allocate(kpath(4,3))
    kpath(1,:)=[0d0,0d0,0d0]
    kpath(2,:)=[ pi,0d0,0d0]
    kpath(3,:)=[ pi, pi,0d0]
    kpath(4,:)=[0d0,0d0,0d0]
    call TB_solve_model(hk_model,Nlso,kpath,Nkpath,&
         colors_name=[red1,green1,blue1],&
         points_name=[character(len=10) :: "G","X","M", "G"],&
         file="Eigenbands.nint")



    !Build the local GF:
    Gmats=zero
    Greal=zero
    fooSmats =zero
    fooSreal =zero
    call add_ctrl_var(beta,"BETA")
    call add_ctrl_var(xmu,"xmu")
    call add_ctrl_var(wini,"wini")
    call add_ctrl_var(wfin,"wfin")
    call add_ctrl_var(eps,"eps")
    call dmft_gloc_matsubara(Hk,Wtk,Gmats,fooSmats,iprint=1)
    call dmft_gloc_realaxis(Hk,Wtk,Greal,fooSreal,iprint=1)
    !
  end subroutine build_hk



  !----------------------------------------------------------------------------------------!
  ! purpose: save the matsubare local self-energy on disk
  !----------------------------------------------------------------------------------------!
  subroutine save_sigma_mats(Smats)
    complex(8),intent(inout)         :: Smats(:,:,:,:,:,:)
    character(len=30)                :: suffix
    integer                          :: ilat,ispin,iorb
    real(8),dimension(:),allocatable :: wm

    if(size(Smats,2)/=Nspin) stop "save_sigma: error in dim 2. Nspin"
    if(size(Smats,3)/=Nspin) stop "save_sigma: error in dim 3. Nspin"
    if(size(Smats,4)/=Norb) stop "save_sigma: error in dim 4. Norb"
    if(size(Smats,5)/=Norb) stop "save_sigma: error in dim 5. Norb"

    allocate(wm(Lmats))

    wm = pi/beta*(2*arange(1,Lmats)-1)
    write(LOGfile,*)"write spin-orbital diagonal elements:"
    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
          call store_data("LSigma"//trim(suffix),Smats(:,ispin,ispin,iorb,iorb,:),wm)
       enddo
    enddo

  end subroutine save_sigma_mats


  !----------------------------------------------------------------------------------------!
  ! purpose: save the real local self-energy on disk
  !----------------------------------------------------------------------------------------!
  subroutine save_sigma_real(Sreal)
    complex(8),intent(inout)         :: Sreal(:,:,:,:,:,:)
    character(len=30)                :: suffix
    integer                          :: ilat,ispin,iorb
    real(8),dimension(:),allocatable :: wm,wr

    if(size(Sreal,2)/=Nspin) stop "save_sigma: error in dim 2. Nspin"
    if(size(Sreal,3)/=Nspin) stop "save_sigma: error in dim 3. Nspin"
    if(size(Sreal,4)/=Norb) stop "save_sigma: error in dim 4. Norb"
    if(size(Sreal,5)/=Norb) stop "save_sigma: error in dim 5. Norb"

    allocate(wr(Lreal))

    wr = linspace(wini,wfin,Lreal)
    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          call store_data("LSigma"//trim(suffix),Sreal(:,ispin,ispin,iorb,iorb,:),wr)
       enddo
    enddo

  end subroutine save_sigma_real




end program ed_hm_3bands



