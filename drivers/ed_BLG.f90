program ed_BLG
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  abstract interface
     function hk_model_template(kpoint,Nlso) result(hk)
       real(8),dimension(:)            :: kpoint
       integer                         :: Nlso
       complex(8),dimension(Nlso,Nlso) :: hk
     end function hk_model_template
  end interface


  procedure(hk_model_template),pointer          :: hk_graphene_model
  character(len=32)                             :: stacking   

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
  complex(8),allocatable,dimension(:,:)         :: graphHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  real(8),allocatable,dimension(:)              :: Wtk

  integer,allocatable,dimension(:)              :: ik2ix,ik2iy
  real(8),dimension(2)                          :: d1,d2,d3
  real(8),dimension(2)                          :: a1,a2,a3
  real(8),dimension(2)                          :: bk1,bk2,bklen
  real(8),dimension(2)                          :: pointM,pointK,pointKp

  !variables for the model:
  integer                                       :: Nk,Nkpath
  real(8)                                       :: ts,t0,t3,tsp,phi,delta,Mh,wmixing
  character(len=32)                             :: finput
  character(len=32)                             :: hkfile
  logical                                       :: spinsym
  !
  real(8),dimension(2)                          :: Eout
  real(8),allocatable,dimension(:)              :: dens
  !
  real(8),dimension(:,:),allocatable            :: Zmats
  complex(8),dimension(:,:,:),allocatable       :: Zfoo
  complex(8),allocatable,dimension(:,:,:,:,:)   :: S0


  !Parse additional variables && read Input && read H(k)^2x2
  call parse_cmd_variable(finput,"FINPUT",default='inputGRAPHENE.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS","inputGRAPHENE.conf",default=1d0)
  call parse_input_variable(t0,"T0","inputGRAPHENE.conf",default=0.142d0*ts)
  call parse_input_variable(t3,"T3","inputGRAPHENE.conf",default=0.107d0*ts)
  call parse_input_variable(mh,"MH","inputGRAPHENE.conf",default=0d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  !
  call parse_input_variable(stacking,"STACKING",finput,default='AB')
  !
  select case(trim(stacking))
  case ('AA')
     hk_graphene_model => hk_BLG_AA_model
  case ('AB')
     hk_graphene_model => hk_BLG_AB_model
  end select

  call ed_read_input(trim(finput))

  if(Norb/=1)stop "Wrong setup from input file: Norb=1"
  Nlat=4
  Nso=Nspin*Norb
  Nlso=Nlat*Nso


  !LATTICE BASIS:
  ! nearest neighbor: A-->B, B-->A
  d1= [  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= [  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= [ -1d0     , 0d0           ]
  !
  !
  !next nearest-neighbor displacements: A-->A, B-->B, cell basis
  a1 = d2-d3                    !a*sqrt(3)[sqrt(3)/2,-1/2]
  a2 = d3-d1                    !a*sqrt(3)[-sqrt(3)/2,-1/2]
  a3 = d1-d2                    !a*sqrt(3)[0, 1]
  !
  !
  !RECIPROCAL LATTICE VECTORS:
  bklen=4d0*pi/sqrt(3d0)
  bk1=bklen*[ sqrt(3d0)/2d0 ,  1d0/2d0 ]
  bk2=bklen*[ sqrt(3d0)/2d0 , -1d0/2d0 ]

  call TB_set_bk([1d0,0d0]*bklen,[0d0,1d0]*bklen)



  pointM = [2*pi/3, 0d0             ]
  pointK = [2*pi/3, 2*pi/3/sqrt(3d0)]
  pointKp= [2*pi/3,-2*pi/3/sqrt(3d0)]


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
     call ed_solve(Bath,Hloc)

     call ed_get_sigma_matsubara(Smats,Nlat)

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

  call dmft_print_gf_matsubara(Smats,"Smats",iprint=1)

  ! extract and print retarded self-energy and Green's function 
  call ed_get_sigma_real(Sreal,Nlat)
  call dmft_print_gf_realaxis(Sreal,"Sreal",iprint=1)
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)


contains




  !--------------------------------------------------------------------!
  ! Bilayer Graphene HAMILTONIAN:
  !--------------------------------------------------------------------!
  !
  ! AA-stacking
  function hk_BLG_AA_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)                :: kpoint
    integer                             :: Nlso,ii,jj
    complex(8),dimension(Nlso,Nlso)     :: hk       
    complex(8),dimension(Nlso/2,Nlso/2) :: hk0       
    real(8)                             :: h0,hx,hy,hz
    real(8)                             :: kdotd(3),kdota(3)
    !(k.d_j)
    kdotd(1) = dot_product(kpoint,d1)
    kdotd(2) = dot_product(kpoint,d2)
    kdotd(3) = dot_product(kpoint,d3)
    !(k.a_j)
    kdota(1) = dot_product(kpoint,a1)
    kdota(2) = dot_product(kpoint,a2)
    kdota(3) = dot_product(kpoint,a3)
    !
    h0 = 2*tsp*cos(phi)*sum( cos(kdota(:)) )
    hx =-ts*sum( cos(kdotd(:)) )
    hy =-ts*sum( sin(kdotd(:)) )
    hz = 2*tsp*sin(phi)*sum( sin(kdota(:)) ) + Mh 
    hk0 = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z
    !
    ! non-interacting graphene layers
    do ii=1,Nlso
       do jj=1,Nlso
          hk(ii,jj)=(0.d0,0.d0)
       end do
    end do
    !
    do ii=1,2
       do jj=1,2
          hk(ii,jj)=hk0(ii,jj)
          hk(ii+2,jj+2)=hk0(ii,jj)
       end do
    end do
    ! 
    ! interaction between planes
    hk(3,1)=t0
    hk(1,3)=t0    
    !
    hk(2,4)=t0
    hk(4,2)=t0    
    !
  end function hk_BLG_AA_model


  !
  ! AB-stacking
  function hk_BLG_AB_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)                :: kpoint
    integer                             :: Nlso,ii,jj
    complex(8),dimension(Nlso,Nlso)     :: hk
    complex(8),dimension(Nlso/2,Nlso/2) :: hk0
    real(8)                             :: h0,hx,hy,hz
    real(8)                             :: kdotd(3),kdota(3)
    !(k.d_j)
    kdotd(1) = dot_product(kpoint,d1)
    kdotd(2) = dot_product(kpoint,d2)
    kdotd(3) = dot_product(kpoint,d3)
    !(k.a_j)
    kdota(1) = dot_product(kpoint,a1)
    kdota(2) = dot_product(kpoint,a2)
    kdota(3) = dot_product(kpoint,a3)
    !
    h0 = 2*tsp*cos(phi)*sum( cos(kdota(:)) )
    hx =-ts*sum( cos(kdotd(:)) )
    hy =-ts*sum( sin(kdotd(:)) )
    hz = 2*tsp*sin(phi)*sum( sin(kdota(:)) ) + Mh
    hk0 = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z
    !
    ! non-interacting graphene layers
    do ii=1,Nlso
       do jj=1,Nlso
          hk(ii,jj)=(0.d0,0.d0)
       end do
    end do
    !
    do ii=1,2
       do jj=1,2
          hk(ii,jj)=hk0(ii,jj)
          hk(ii+2,jj+2)=hk0(ii,jj)
       end do
    end do
    !
    ! interaction between planes
    hk(3,2)=t0
    hk(2,3)=t0
    !
    ! <<<<<<< HEAD
    ! hk(1,4)=t3*hk0(2,1)
    ! hk(4,1)=t3*hk0(1,2)
    ! =======
    hk(1,4)=t3*hk0(1,2)
    hk(4,1)=t3*hk0(2,1)
    ! >>>>>>> 17a0ecd2fcea8e738ea0d35d160286970e3ee752
    !
  end function hk_BLG_AB_model



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
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: Gmats,fooSmats
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: Greal,fooSreal
    real(8),dimension(2)                                   :: kvec
    real(8)                                                :: blen,area_hex,area_rect,points_in,points_tot
    real(8),allocatable,dimension(:)                       :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable                     :: KPath

    Lk= Nk*Nk

    write(LOGfile,*)"Build H(k) of AA Bilayer Graphene:",Lk
    write(LOGfile,*)"# of SO-bands     :",Nlso

    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nlso,Nlso,Lk));Hk=zero
    allocate(wtk(Lk));Wtk=0d0
    call TB_build_model(Hk,hk_graphene_model,Nlso,[Nk,Nk])
    Wtk = 1d0/Lk


    if(present(file))then
       call TB_write_hk(Hk,"Hkrfile_BLG_AA.data",&
            No=Nlso,&
            Nd=Norb,&
            Np=0,&
            Nineq=1,&
            Nkvec=[Nk,Nk])
    endif
    !
    allocate(graphHloc(Nlso,Nlso))
    graphHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(graphHloc))<1.d-4)graphHloc=0d0
    call TB_write_Hloc(graphHloc)
    call TB_write_Hloc(graphHloc,'Hloc.txt')
    !
    !


    allocate(Kpath(4,2))
    KPath(1,:)=pointK
    KPath(2,:)=[0,0]
    KPath(3,:)=pointM
    KPath(4,:)=pointK
    call TB_Solve_model(hk_graphene_model,Nlso,KPath,Nkpath,&
         colors_name=[red1,blue1,green1,purple1],&
         points_name=[character(len=10) :: "K","G","M","K"],&
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
    call dmft_gloc_matsubara(Hk,Wtk,Gmats,fooSmats)
    call dmft_print_gf_matsubara(Gmats,"LG0",iprint=1)
    call dmft_gloc_realaxis(Hk,Wtk,Greal,fooSreal)
    call dmft_print_gf_realaxis(Greal,"LG0",iprint=1)
    !
  end subroutine build_hk








end program ed_BLG



