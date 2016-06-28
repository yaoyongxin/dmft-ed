program ed_bhz_2d_edge
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI_INEQ
  USE MPI
#endif
  implicit none

  integer                                       :: iloop
  integer                                       :: Nlso
  integer                                       :: Nso
  integer                                       :: Nineq
  integer                                       :: ilat,iy,iorb,ispin,ineq,i
  logical                                       :: converged
  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath_ineq
  real(8),allocatable,dimension(:,:)            :: Bath_prev
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal,Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal
  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hkr
  complex(8),allocatable,dimension(:,:)         :: bhzHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc,Hloc_ineq,S0

  !gamma matrices:
  complex(8),allocatable,dimension(:,:)         :: gamma1
  complex(8),allocatable,dimension(:,:)         :: gamma2
  complex(8),allocatable,dimension(:,:)         :: gamma5
  real(8),allocatable,dimension(:)              :: Wtk
  real(8),allocatable,dimension(:)              :: kxgrid
  real(8),dimension(:,:),allocatable            :: kpath
  integer                                       :: Nk,Ly,Nkpath
  real(8)                                       :: e0,mh,lambda,wmixing
  logical                                       :: spinsym,tridiag,lrsym,rebuild_sigma
  character(len=60)                             :: finput
  character(len=32)                             :: hkfile
  real(8),dimension(:,:),allocatable            :: Zmats
  complex(8),dimension(:,:,:),allocatable       :: Zfoo


#ifdef _MPI_INEQ
  ! START MPI !
  call MPI_INIT(mpiERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpiID,mpiERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpiSIZE,mpiERR)
  write(*,"(A,I4,A,I4,A)")'Processor ',mpiID,' of ',mpiSIZE,' is alive'
  call MPI_BARRIER(MPI_COMM_WORLD,mpiERR)
#endif

  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ_EDGE.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(Ly,"Ly",finput,default=20)
  call parse_input_variable(Nkpath,"NKPATH",finput,default=501)
  call parse_input_variable(tridiag,"TRIDIAG",finput,default=.true.)
  call parse_input_variable(mh,"MH",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(e0,"e0",finput,default=1d0)
  call parse_input_variable(lrsym,"LRSYM",finput,default=.true.)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(rebuild_sigma,"REBUILD_SIGMA",finput,default=.false.)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  !
  call ed_read_input(trim(finput))

  !set the global number of lattice sites equal to the number of layers along the y-axis
  Nlat = Ly
  Nineq= Ly
  if(lrsym)then
     if(mod(Ly,2)/=0)stop "Wrong setup from input file: Ly%2 > 0 (odd number of sites)"
     Nineq=Ly/2
     print*,"Using L-R Symmetry. Solve",Nineq," of",Nlat," sites."
     call sleep(2)
  endif

  !set the local number of total spin-orbitals (4)
  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso  = Nspin*Norb

  !set the total lattice-spin-orbit dimension:
  Nlso=Nlat*Nspin*Norb



  !Allocate Functions:
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb));Hloc=zero
  allocate(S0(Nlat,Nspin,Nspin,Norb,Norb));S0=zero
  allocate(Zmats(Nlso,Nlso));Zmats=eye(Nlso)
  allocate(Zfoo(Nlat,Nso,Nso));Zfoo=0d0
  allocate(Weiss(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Smats_ineq=zero
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal));Sreal_ineq=zero
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats));Gmats_ineq=zero
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb));Hloc_ineq=zero


  !Buil the Hamiltonian on a grid or on  path
  call build_hkr(trim(hkfile))
  Hloc = lso2nnn_reshape(bhzHloc,Nlat,Nspin,Norb)
  do ineq=1,Nineq
     ilat = ineq2ilat(ineq)
     Hloc_ineq(ineq,:,:,:,:) = Hloc(ilat,:,:,:,:)
  enddo


  !Setup solver
  Nb=get_bath_size()
  allocate(Bath_ineq(Nineq,Nb) )
  allocate(Bath_prev(Nineq,Nb) )
  call ed_init_solver_lattice(Bath_ineq)


  if(rebuild_sigma)then
     call ed_rebuild_sigma_lattice(Bath_ineq,Hloc_ineq,iprint=1)
     call ed_get_sigma_matsubara_lattice(Smats_ineq,Nineq)
     do ilat=1,Nlat
        ineq = ilat2ineq(ilat)
        S0(ilat,:,:,:,:)      = Smats_ineq(ineq,:,:,:,:,1)
     enddo
     do ilat=1,Nlat
        Zfoo(ilat,:,:)        = select_block(ilat,S0)
        do iorb=1,Nso
           i = iorb + (ilat-1)*Nso
           Zmats(i,i)  = 1.d0/( 1.d0 + abs( dimag(Zfoo(ilat,iorb,iorb))/(pi/beta) ))
        enddo
     enddo
     if(mpiID==0)call build_eigenbands()
     stop
  endif





  !DMFT loop:
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(mpiID==0) call start_loop(iloop,nloop,"DMFT-loop")   
     ! solve the impurities on each inequivalent y-layer
     call ed_solve_lattice(Bath_ineq,Hloc_ineq,iprint=1)
     ! retrieve the self-energies
     ! store the 1st Matsubara freq. into S0, used to get H_topological = Hk + S0
     call ed_get_sigma_matsubara_lattice(Smats_ineq,Nineq)
     call ed_get_sigma_real_lattice(Sreal_ineq,Nineq)
     do ilat=1,Nlat
        ineq = ilat2ineq(ilat)
        Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
        Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
        S0(ilat,:,:,:,:)      = Smats_ineq(ineq,:,:,:,:,1)
     enddo
     do ilat=1,Nlat
        Zfoo(ilat,:,:)        = select_block(ilat,S0)
        do iorb=1,Nso
           i = iorb + (ilat-1)*Nso
           Zmats(i,i)  = 1.d0/( 1.d0 + abs( dimag(Zfoo(ilat,iorb,iorb))/(pi/beta) ))
        enddo
     enddo
     ! compute the local gf:
     call ed_get_gloc_lattice(Hkr,Wtk,Gmats,Greal,Smats,Sreal,iprint=1,tridiag=tridiag)
     do ineq=1,Nineq
        ilat = ineq2ilat(ineq)
        Gmats_ineq(ineq,:,:,:,:,:) = Gmats(ilat,:,:,:,:,:)
     enddo
     ! compute the Weiss field (only the Nineq ones)
     call ed_get_weiss_lattice(Gmats_ineq,Smats_ineq,Weiss,Hloc_ineq,iprint=1)
     ! fit baths and mix result with old baths
     call ed_chi2_fitgf_lattice(Bath_ineq,Weiss,Hloc_ineq,ispin=1)
     if(spinsym)then
        call spin_symmetrize_bath(Bath_ineq)
     else
        call ed_chi2_fitgf_lattice(Bath_ineq,Weiss,Hloc_ineq,ispin=2)
     endif
     Bath_ineq=wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath_ineq
     if(mpiID==0)converged = check_convergence(Weiss(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
#ifdef _MPI_INEQ
     call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpiERR)
#endif
     if(mpiID==0)call end_loop
  enddo

  if(mpiID==0)call build_eigenbands()

#ifdef _MPI_INEQ
  call MPI_FINALIZE(mpiERR)
#endif

contains



  !+-----------------------------------------------------------------------------+!
  !PURPOSE: build the BHZ Hamiltonian H(k_x,R_y) on the STRIPE along Y
  !+-----------------------------------------------------------------------------+!
  subroutine build_hkr(file)
    character(len=*),optional          :: file
    integer :: i,ik
    !
    !SETUP THE GAMMA MATRICES:
    allocate(gamma1(Nso,Nso),gamma2(Nso,Nso),gamma5(Nso,Nso))
    gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x )
    gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y )
    gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z )
    !
    !SETUP THE H(kx,Ry):
    if(mpiID==0)then
       write(LOGfile,*)"Build H(kx,y) for BHZ-stripe:"
       write(*,*)"# of kx-points     :",Nk
       write(*,*)"# of y-layers      :",Nlat
    endif
    !
    if(allocated(Kxgrid))deallocate(Kxgrid)
    allocate(Kxgrid(Nk))
    if(allocated(Hkr))deallocate(Hkr)
    allocate(Hkr(Nlso,Nlso,Nk))
    kxgrid = kgrid(Nk)
    Hkr    = TB_build_model(bhz_edge_model,Ly,Nso,kxgrid,[0d0],[0d0],pbc=.false.)
    if(mpiID==0)call write_hk_w90("Hkrfile.in",&
         No=Nlso,&
         Nd=Norb,&
         Np=0,&
         Nineq=Ly,&
         Hk=Hkr,&
         kxgrid=kxgrid,kygrid=[0d0],kzgrid=[0d0])
    if(allocated(Wtk))deallocate(Wtk)
    allocate(Wtk(Nk))
    Wtk = 1d0/Nk
    !
    !SETUP THE LOCAL PART Hloc(Ry)
    allocate(bhzHloc(Nlso,Nlso))
    bhzHloc = extract_Hloc(Hkr,Nlat,Nspin,Norb)
  end subroutine build_hkr





  !----------------------------------------------------------------------------------------!
  ! purpose: read the local self-energy from disk
  !----------------------------------------------------------------------------------------!
  subroutine read_sigma_matsubara(Self)
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Self
    character(len=30)                             :: suffix
    integer                                       :: ilat,ispin,iorb
    real(8),dimension(:),allocatable              :: wm
    call assert_shape(Self,[Nineq,Nspin,Nspin,Norb,Norb,Lmats],"read_sigma_matsubara","Self_ineq")
    allocate(wm(Lmats))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    if(mpiID==0)then
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
             call read_data("LSigma"//trim(suffix),Self(:,ispin,ispin,iorb,iorb,:),wm)
          enddo
       enddo
    endif
  end subroutine read_sigma_matsubara

  subroutine read_sigma_real(Self)
    complex(8),allocatable,dimension(:,:,:,:,:,:) :: Self
    character(len=30)                             :: suffix
    integer                                       :: ilat,ispin,iorb
    real(8),dimension(:),allocatable              :: wr
    call assert_shape(Self,[Nineq,Nspin,Nspin,Norb,Norb,Lreal],"read_sigma_real","Self_ineq")
    allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    if(mpiID==0)then
       do ispin=1,Nspin
          do iorb=1,Norb
             suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
             call read_data("LSigma"//trim(suffix),Self(:,ispin,ispin,iorb,iorb,:),wr)
          enddo
       enddo
    endif
  end subroutine read_sigma_real






  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve H_BHZ(k_x,R_y) along the 1d -pi:pi path in the BZ.
  !+-----------------------------------------------------------------------------+!  
  subroutine build_eigenbands(kpath_)
    real(8),dimension(:,:),optional    :: kpath_
    real(8),dimension(:,:),allocatable :: kpath
    integer                            :: Npts
    character(len=64)                  :: file
    if(present(kpath_))then
       if(mpiID==0)write(LOGfile,*)"Solve H(kx,y) along a given path:"
       Npts = size(kpath_,1)
       allocate(kpath(Npts,size(kpath_,2)))
       kpath=kpath_
       file="Eigenbands_path.nint"
    else
       !PRINT H(kx,Ry) ALONG A -pi:pi PATH
       if(mpiID==0)write(LOGfile,*)"Solve H(kx,y) along [-pi:pi]:"
       Npts=3
       allocate(Kpath(Npts,1))
       kpath(1,:)=[-1]*pi
       kpath(2,:)=[ 0]*pi
       kpath(3,:)=[ 1]*pi
       file="Eigenbands.nint"
    endif
    call TB_solve_path(bhz_edge_model,Ly,Nso,kpath,Nkpath,&
         colors_name=[red1,gray88,blue1,gray88,blue1,gray88,red1,gray88],&
         points_name=[character(len=10) :: "-pi","0","pi"],&
         file="Eigenbands.nint",pbc=.false.)
    ! call solve_HkR_along_BZpath(bhz_edge_model,Ly,Nso,kpath,Nkpath,reg(file),pbc=.false.)
  end subroutine build_eigenbands







  !+-----------------------------------------------------------------------------+!
  !PURPOSE: the BHZ-edge model hamiltonian
  !+-----------------------------------------------------------------------------+!
  !BHZ on a stripe geometry;
  function bhz_edge_model(kpoint,Nlat,N,pbc) result(Hrk)
    real(8),dimension(:)                :: kpoint
    real(8)                             :: kx
    integer                             :: Nlat,N
    complex(8),dimension(N,N)           :: Hmat,Tmat,TmatH
    complex(8),dimension(Nlat*N,Nlat*N) :: Hrk
    integer                             :: i,Idmin,Idmax,Itmin,Itmax
    logical                             :: pbc
    kx=kpoint(1)
    Hrk=zero
    Hmat=h0_rk_bhz(kx,N)
    Tmat=t0_rk_bhz(N)
    TmatH=conjg(transpose(Tmat))
    do i=1,Nlat
       Idmin=1+(i-1)*N
       Idmax=      i*N
       Hrk(Idmin:Idmax,Idmin:Idmax)=Hmat + dreal(select_block(i,S0)) !< H(k) + Re(Sigma_iy(:Nso,:Nso;omega=0))
    enddo
    do i=1,Nlat-1
       Idmin=1 + (i-1)*N
       Idmax=        i*N
       Itmin=1 +     i*N
       Itmax=    (i+1)*N
       Hrk(Idmin:Idmax,Itmin:Itmax)=Tmat
       Hrk(Itmin:Itmax,Idmin:Idmax)=TmatH
    enddo
    if(pbc)then
       Itmin=1+(Nlat-1)*N
       Itmax=0+Nlat*N
       Hrk(1:N,Itmin:Itmax)=TmatH
       Hrk(Itmin:Itmax,1:N)=Tmat
    endif
    Hrk = matmul(Zmats,Hrk)
  end function bhz_edge_model

  function h0_rk_bhz(kx,N) result(H)
    real(8)                    :: kx
    integer                    :: N
    complex(8),dimension(N,N)  :: H
    H = (mh-e0*cos(kx))*gamma5 + lambda*sin(kx)*gamma1
  end function h0_rk_bhz

  function t0_rk_bhz(N) result(H)
    integer                    :: N
    complex(8),dimension(N,N)  :: H
    H = -0.5d0*e0*gamma5 + xi*0.5d0*lambda*gamma2
  end function T0_rk_bhz





  function ilat2ineq(ilat) result(ineq)
    integer,intent(in) :: ilat
    integer            :: ineq
    ineq=ilat
    if( lrsym .AND. (ilat>Nineq) )ineq=Nlat-ilat+1
  end function ilat2ineq

  function ineq2ilat(ineq) result(ilat)
    integer,intent(in) :: ineq
    integer            :: ilat
    ilat=ineq
    if(ineq>Nineq)stop "ineq2ilat error: called with ineq > Nineq"
  end function ineq2ilat


end program ed_bhz_2d_edge
