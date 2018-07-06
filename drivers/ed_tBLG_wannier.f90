program ed_effective_tBLG
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS

!!!!!!!!:::::::: The twisted bilayer graphene flat bands effective model can be found in arXiv:1805.06819:::::!!!!!!!!
  ! this driver needs the file 'eff_hopping.dat' to work


  implicit none

  real(8),parameter                             :: eV_to_meV=1000d0
  integer                                       :: i,iloop,Lk,Nso,Nlso,Nlat
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
  real(8),dimension(2)                          :: L1,L2,G1,G2
  real(8),dimension(2)                          :: bk1,bk2,bklen
  real(8),dimension(2)                          :: pointM,pointK1,pointK2

  !variables for the model:
  integer                                       :: Nk,Nkpath,m,n,p,q
  real(8)                                       :: Lm,Mh,wmixing,re,im
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
  complex(8),dimension(1162)                    :: hopping   !1162 is the number of hopping constants
  integer,dimension(1162,4)                     :: integer_input


  !Parse additional variables && read Input && read H(k)^2x2
  call parse_cmd_variable(finput,"FINPUT",default='inputGRAPHENE.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(mh,"MH","inputGRAPHENE.conf",default=0d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  !
  call ed_read_input(trim(finput))

  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")


  if(Norb/=1)stop "Wrong setup from input file: Norb=1"
  Nlat=4
  Nso=Nspin*Norb
  Nlso=Nlat*Nso
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

!!!! all the hopping constants are contained in the file 'eff_hopping.dat'!!!!!

  open(10,file='eff_hopping.dat')
  do i=1,1162
     read(10,*) integer_input(i,1),integer_input(i,2),integer_input(i,3),integer_input(i,4),re,im
     hopping(i)=eV_to_meV*(re+xi*im) !transforms to meV
  end do
  close(10)

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

  !Build the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))
  Hloc = lso2nnn_reshape(graphHloc,Nlat,Nspin,Norb)

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

  call dmft_print_gf_matsubara(Gmats,"Gmats",iprint=1)
  call dmft_print_gf_matsubara(Smats,"Smats",iprint=1)

  ! extract and print retarded self-energy and Green's function 
  call ed_get_sigma_real(Sreal,Nlat)
  call dmft_print_gf_realaxis(Sreal,"Sreal",iprint=1)
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)


contains



  !--------------------------------------------------------------------!
  ! Twisted Bilayer Graphene (theta=1.05) effective  HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_TBGeffective_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)                :: kpoint
    integer                             :: Nlso,i,j,ii,jj,m,n,p,q
    complex(8)                          :: uni
    real(8)                             :: re,im
    complex(8),dimension(Nlso,Nlso)     :: hk       
    complex(8),dimension(Nlso/2,Nlso/2) :: hk_plus,hk_minus
    real(8),dimension(2)                :: pos_1,pos_2,r_12  
    real(8),dimension(2)                :: r_AB,r_BA,RR

!    r_BA=[1.d0/2.d0,sqrt(3.d0)/2.d0]*(Lm/sqrt(3.d0))
!    r_AB=[-1.d0/2.d0,sqrt(3.d0)/2.d0]*(Lm/sqrt(3.d0))

    hk_plus=zero
    hk_minus=zero

    do i=1,1162

       m=integer_input(i,1)
       n=integer_input(i,2)    
       p=integer_input(i,3)
       q=integer_input(i,4)

       RR=dfloat(m)*L1+dfloat(n)*L2 

!       if (p.eq.1) then  
!          pos_1=RR+r_BA
!       else
!          pos_1=RR+r_AB
!       end if

!       if (q.eq.1) then  
!          pos_2=r_BA
!       else
!          pos_2=r_AB
!       end if

!       r_12=pos_1-pos_2

       hk_plus(p,q)=hk_plus(p,q)+(exp(xi*dot_product(kpoint,RR)))*hopping(i)
       hk_minus(p,q)=hk_minus(p,q)+(exp(-xi*dot_product(kpoint,RR)))*hopping(i)

    end do

    hk=zero

    do i=1,2
       do j=1,2
          hk(i,j)=hk_plus(i,j)
       end do
    end do

    do i=3,4
       do j=3,4
          hk(i,j)=hk_minus(i-2,j-2)
       end do
    end do

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
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: Gmats
    real(8),dimension(2)                                   :: kvec
    real(8)                                                :: blen,area_hex,area_rect,points_in,points_tot
    real(8),allocatable,dimension(:)                       :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable                     :: KPath

    Lk= Nk*Nk

    write(LOGfile,*)"Build effective H(k) of twisted Bilayer Graphene (theta=1.05):",Lk
    write(LOGfile,*)"# of SO-bands     :",Nlso

    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nlso,Nlso,Lk));Hk=zero
    allocate(wtk(Lk));Wtk=0d0
    ! allocate(kxgrid(Nk),kygrid(Nk))
    ! ik=0
    ! do iy=1,Nk
    !    ky = dble(iy-1)/Nk
    !    do ix=1,Nk
    !       ik=ik+1
    !       kx=dble(ix-1)/Nk
    !       kvec = kx*bk1 + ky*bk2
    !       kxgrid(ix) = kvec(1)
    !       kygrid(iy) = kvec(2)
    !       Hk(:,:,ik) = hk_TBGeffective_model(kvec,Nlso)
    !    enddo
    ! enddo
    call TB_build_model(Hk,hk_TBGeffective_model,Nlso,[Nk,Nk])
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
    KPath(1,:)=pointK1
    KPath(2,:)=[0,0]
    KPath(3,:)=pointM
    KPath(4,:)=pointK2
    call TB_Solve_model(hk_TBGeffective_model,Nlso,KPath,Nkpath,&
         colors_name=[red1,blue1,green1,purple1],&
         points_name=[character(len=10) :: "K1","G","M","K2"],&
         file="Eigenbands.nint")

    !Build the local GF:
    call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    call dmft_print_gf_matsubara(Gmats,"LG0",iprint=1)
    !
  end subroutine build_hk








end program ed_effective_tBLG



