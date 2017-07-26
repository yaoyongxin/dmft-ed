program ed_wsm_3d
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                       :: iloop,Lk,Nso,test=0,uno,due,tre
  logical                       :: converged
  !Bath:
  integer                       :: Nb
  real(8),allocatable           :: Bath(:),Bath_Prev(:)
  !The local hybridization function:
  complex(8),allocatable        :: Weiss(:,:,:,:,:)
  complex(8),allocatable        :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable        :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable        :: Hk(:,:,:),wsmHloc(:,:),sigmaWSM(:,:)
  real(8),allocatable           :: Wtk(:)
  integer,allocatable           :: ik2ix(:),ik2iy(:),ik2iz(:)
  !variables for the model:
  integer                       :: Nk,Nkpath,CHERN_INDEX,ICOMP,JCOMP
  real(8)                       :: e0,mh,lambda,bx,by,bz,BIA,uno_r,due_r,tre_r
  real(8)                       :: wmixing
  real(8),dimension(5)          :: PLANE_INDEX
  character(len=16)             :: finput
  character(len=32)             :: hkfile
  logical                       :: spinsym,getpoles
  complex(8),dimension(4,4)     :: Gamma1,Gamma2,Gamma3,Gamma5
  real(8),dimension(4,3)        :: TEST2
  real(8),dimension(4)          :: TEST1,TEST3 
  real(8),dimension(3)        :: kpoint
  call parse_cmd_variable(finput,"FINPUT",default='inputED_WSM.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=30)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=1d0)
  call parse_input_variable(e0,"E0",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.5d0)
  call parse_input_variable(bx,"BX",finput,default=0.3d0)
  call parse_input_variable(by,"BY",finput,default=0.d0)
  call parse_input_variable(bz,"BZ",finput,default=0.d0)
  call parse_input_variable(BIA,"BIA",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  !
  call ed_read_input(trim(finput))

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")


  !SETUP THE GAMMA MATRICES:
  gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gamma3=kron_pauli( pauli_tau_x,-pauli_sigma_x)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)


  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb
  !
  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(SigmaWSM(Nso,Nso))
  call set_sigmaWSM()           !this set sigma_WSM(0) to zero
  !
  !
  !
  !
  !
  !Buil the Hamiltonian on a grid or on path
  call build_hk(trim(hkfile))!
  !
  !
  !Setup solver
  Nb=get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_Prev(Nb))
  call ed_init_solver(bath,j2so(wsmHloc))
  !
  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     !
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     !
     !
     ! retrieve the self-energies
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     call build_z2_indices( so2j(Smats(:,:,:,:,1),Nso) ) !so2j transforms a Nspin:Nspin:Norb:Norb into a Nso:Nso matrix
     !
     !
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats,iprint=3)
     !
     !
     !
     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats,Smats,Weiss,Hloc=j2so(wsmHloc),iprint=3)
     else
        call dmft_delta(Gmats,Smats,Weiss,Hloc=j2so(wsmHloc),iprint=3)
     endif
     !
     !
     !
     !Fit the new bath, starting from the old bath + the supplied Weiss/Delta
     if(ed_mode=="normal")then
        call ed_chi2_fitgf(Weiss,bath,ispin=1)
        call ed_chi2_fitgf(Weiss,bath,ispin=2)
     else
        call ed_chi2_fitgf(Weiss,bath)
     endif
     !
     !
     !
     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_Prev
     Bath_Prev=Bath
     !
     !
     !
     converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop)
     !
     call end_loop
     !
  enddo

  
  ! compute the local gf:
  !call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal,iprint=3)
  
  !Get kinetic energy:
  !call dmft_kinetic_energy(Hk,Wtk,Smats)


  !Get 3d Bands from Top. Hamiltonian
  call solve_hk_topological( so2j(Smats(:,:,:,:,1),Nso) )

  !Find out if it is a semimetal
  call is_weyl_brutal( Nso, so2j(Smats(:,:,:,:,1),Nso) )

  call chern_retriever([0.0d0,0.0d0,0.0d0])
  write(*,*) "ADESSO DEI TEST"
  !Find out weyl point chirality
  do uno=0,10
    uno_r=-pi+(pi/5)*uno
    do due=0,10
      due_r=-pi+(pi/5)*due
      do tre=0,10
       tre_r=-pi+(pi/5)*tre
       write(*,*) uno_r,due_r,tre_r
       call chern_retriever([uno_r,due_r,tre_r])
      enddo
    enddo
  enddo
contains



  !---------------------------------------------------------------------
  !PURPOSE: GET WSM HAMILTONIAN IN THE FULL BZ
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,k,ik,iorb,jorb,io,ispin
    integer                             :: ix,iy,iz
    real(8)                             :: kx,ky,kz
    real(8)                             :: foo
    integer                             :: unit
    !
    !get H(k) and solve the non-interacting problem along a path in 3d:
    call build_hk_path()
    !
    !Get H(k) in the BZ:    
    Lk=Nk**3
    !
    write(LOGfile,*)"Build H(k) for WSM:"
    write(*,*)"# of k-points     :",Lk
    !
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk))
    allocate(Wtk(Lk))
    call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])
    call TB_build_model(Hk,hk_weyl,Nso,[Nk,Nk,Nk])
    Wtk = 1d0/Lk
    if(present(file))call TB_write_hk(Hk,trim(file),Nso,&
         Nd=Norb,Np=1,Nineq=1,&
         Nkvec=[Nk,Nk,Nk])
    !
    !
    !   
    !GET LOCAL PART OF THE HAMILTONIAN
    if(allocated(wsmHloc))deallocate(wsmHloc)
    allocate(wsmHloc(Nso,Nso))
    wsmHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(wsmHloc))<1.d-9)wsmHloc=0d0
    call TB_write_Hloc(wsmHloc)
    !
  end subroutine build_hk




  !---------------------------------------------------------------------
  !PURPOSE: GET THE WSM HAMILTONIAN ALONG A 3D PATH IN THE BZ
  !---------------------------------------------------------------------
  subroutine build_hk_path(kpath_)
    integer                                :: i,j
    integer                                :: Npts
    real(8),dimension(:,:),optional        :: kpath_
    real(8),dimension(:,:),allocatable     :: kpath
    character(len=64)                      :: file
    !
    !This routine build the H(k) along the GXMG path in BZ, Hk(k) is constructed along this path.
    !
    sigmaWSM=zero
    if(present(kpath_))then
       write(LOGfile,*)"Build H(k) WSM along a given path:"
       Npts = size(kpath_,1)
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,size(kpath_,2)))      
       kpath=kpath_
       file="Eig_path.nint"
    else
       write(LOGfile,*)"Build H(k) WSM along the high-symmetry path:"
       Npts = 8
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_m1
       kpath(2,:)=kpoint_x2
       kpath(3,:)=kpoint_gamma
       kpath(4,:)=kpoint_x1
       kpath(5,:)=kpoint_m2
       kpath(6,:)=kpoint_r
       kpath(7,:)=kpoint_x3
       kpath(8,:)=kpoint_gamma
       file="Eigenbands.nint"
    endif
    !
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk))
    allocate(Wtk(Lk))
    !
    call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])
    call TB_build_model(Hk,hk_weyl,Nso,kpath,Nkpath)
    Wtk = 1d0/Lk
    !
    call TB_solve_model(hk_weyl,Nso,kpath,Nkpath,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=10) :: "M","X","G","X1","A","R","Z","G"],&
         file="Eigenband.nint")
  end subroutine build_hk_path


  subroutine solve_hk_topological(sigma)
    integer                                :: i,j
    integer                                :: Npts
    complex(8),dimension(Nso,Nso)          :: sigma(Nso,Nso)
    real(8),dimension(:,:),allocatable     :: kpath
    !
    !This routine build the H(k) along the GXMG path in BZ, Hk(k) is constructed along this path.
    write(LOGfile,*)"Build H_TOP(k) WSM along the X-G-M-R-Z-A-G-Z path:"
    !
    Npts = 8
    Lk=(Npts-1)*Nkpath
    allocate(kpath(Npts,3))
    kpath(1,:)=kpoint_m1
    kpath(2,:)=kpoint_x2
    kpath(3,:)=kpoint_gamma
    kpath(4,:)=kpoint_x1
    kpath(5,:)=kpoint_m2
    kpath(6,:)=kpoint_r
    kpath(7,:)=kpoint_x3
    kpath(8,:)=kpoint_gamma
    call set_sigmaWSM(sigma)
    call TB_solve_model(hk_weyl,Nso,kpath,Nkpath,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=10) :: "M","X","G","X1","A","R","Z","G"],&
         file="Eig_Htop.ed")
  end subroutine solve_hk_topological




  !--------------------------------------------------------------------!
  !WSM HAMILTONIAN:
  !--------------------------------------------------------------------!
  subroutine set_SigmaWSM(sigma)
    complex(8),dimension(Nso,Nso),optional :: sigma(Nso,Nso)
    sigmaWSM = zero;if(present(sigma))sigmaWSM=sigma
  end subroutine set_SigmaWSM

  function hk_weyl(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky,kz
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_weyl: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    kz=kpoint(3)
    !
    Hk          = zero
    Hk(1:2,1:2) = &
         (Mh - e0*(cos(kx) + cos(ky) + cos(kz)) )*pauli_tau_z +&
         lambda*sin(kx)*pauli_tau_x + lambda*sin(ky)*pauli_tau_y +&
         by*pauli_tau_y + bz*pauli_tau_z
    Hk(3:4,3:4) = conjg( &
         (Mh - e0*(cos(-kx) + cos(-ky) + cos(-kz)) )*pauli_tau_z +&
         lambda*sin(-kx)*pauli_tau_x + lambda*sin(-ky)*pauli_tau_y +&
         by*pauli_tau_y - bz*pauli_tau_z)
    Hk(1:2,3:4) = lambda*sin(kz)*pauli_tau_x -xi*BIA*pauli_tau_y + bx*pauli_tau_z
    Hk(3:4,1:2) = lambda*sin(kz)*pauli_tau_x +xi*BIA*pauli_tau_y + bx*pauli_tau_z
    !
    !add the SigmaWSM term to get Topologial Hamiltonian if required:
    Hk = Hk + dreal(SigmaWSM)
    !
  end function hk_weyl

  !--------------------------------------------------------------------!
  !WSM REGION FINDER:
  !--------------------------------------------------------------------!

  subroutine is_weyl_brutal(N,sigma)
    integer                                 :: weyl_index, unit,test_brutal=0,i,j,k,mash_thickness
    real                                    :: step
    complex(8),dimension(Nso,Nso)           :: sigma(Nso,Nso)
    integer                                 :: N
    real(8),dimension(:),allocatable        :: kpoint
    real(8),dimension(N)                    :: eval
    complex(8),dimension(N,N)               :: ham

    if(N/=4)stop "hk_weyl: error in N dimensions"
    !
    !
    weyl_index=0
    call set_sigmaWSM(sigma)
    !
    mash_thickness=10
    step=2*pi/(mash_thickness*Nk)
    !
    write(*,*) "Starting Weyl point search"
    !
    kpoint=[-pi,-pi,-pi]
    xloop: do i=0,mash_thickness+Nk-1
      kpoint(1)=-pi+step*i
      yloop: do j=0,mash_thickness*Nk-1
        kpoint(2)=-pi+step*j
        zloop: do k=0,mash_thickness*Nk-1
          kpoint(3)=-pi+step*k
          test_brutal=test_brutal+1
          ham=hk_weyl(kpoint,N)
          call eigh(ham,Eval)
          if (abs(Eval(3))+abs(Eval(2)) .lt. 0.1) then
            weyl_index=1
            write(*,*) "Brutal iterations: found Weyl point at ",kpoint/pi, "value is ", abs(Eval(3))
            call chern_retriever(kpoint)
            !exit xloop
          end if
        end do zloop
      end do yloop
    end do xloop
    unit=free_unit()
    open(unit,file="is_weyl.ed")
    write(unit,*)weyl_index
    close(unit)
  end subroutine is_weyl_brutal


  subroutine retrieve_third_band(n,x,fvec,iflag)
    integer                                 :: iflag,n, INFO
    real(8),dimension(N)                    :: x
    real(8),dimension(N)                    :: fvec
    complex(8),dimension(Nso,Nso)           :: ham
    complex(8),dimension(Nso,Nso)           :: sigma
    real(8),dimension(Nso)                  :: eval
    !
    sigma=so2j(Smats(:,:,:,:,1),Nso)
    call set_sigmaWSM(sigma)
    ham=hk_weyl(x,Nso)
    call eigh(ham,Eval)
    fvec=[Eval(3),0.0d0,0.0d0]
    test=test+1
  end subroutine retrieve_third_band


  !--------------------------------------------------------------------!
  !CHERN NUMBER EVALUATION  (NEEDS FIXING AND OPTIMIZATION)
  !--------------------------------------------------------------------!

  subroutine chern_retriever(kpoint)   !integrates Berry flux on a cubic surface around the point kpoint
    real(8)                         :: z2
    real(8),dimension(3)            :: kpoint
    real(8)                         :: e,a,b
    integer                         :: unit
    e=pi/10
    write(*,*) "Testing chirality"
    write(*,*) "Integrating on a cube around ",kpoint
    z2=0
    !
    a=kpoint(1)
    b=kpoint(2)
    !
    PLANE_INDEX=[kpoint(3)+e,1.0d0,1.0d0,2.0d0,3.0d0]
    !
    z2=z2-simps2d(chern_flux_integrable,[a-e,a+e],[b-e,b+e],N0=200,iterative=.false.)   
    write(*,*) "Done face 1, Berry flux now is ",z2 
    !
    PLANE_INDEX=[kpoint(3)-e,-1.0d0,1.0d0,2.0d0,3.0d0]
    !
    z2=z2-simps2d(chern_flux_integrable,[a-e,a+e],[b-e,b+e],N0=200,iterative=.false.)   
    write(*,*) "Done face 2, Berry flux now is ",z2
    !
    a=kpoint(3)
    b=kpoint(1)
    !
    PLANE_INDEX=[kpoint(2)+e,1.0d0,3.0d0,1.0d0,2.0d0]
    !
    z2=z2-simps2d(chern_flux_integrable,[a-e,a+e],[b-e,b+e],N0=200,iterative=.false.)   
    write(*,*) "Done face 3, Berry flux now is ",z2 
    !
    PLANE_INDEX=[kpoint(2)-e,-1.0d0,3.0d0,1.0d0,2.0d0]
    !
    z2=z2-simps2d(chern_flux_integrable,[a-e,a+e],[b-e,b+e],N0=200,iterative=.false.)   
    write(*,*) "Done face 4, Berry flux now is ",z2 
    !
    a=kpoint(2)
    b=kpoint(3)
    !
    PLANE_INDEX=[kpoint(1)+e,1.0d0,2.0d0,3.0d0,1.0d0]
    !
    z2=z2-simps2d(chern_flux_integrable,[a-e,a+e],[b-e,b+e],N0=200,iterative=.false.)   
    write(*,*) "Done face 5, Berry flux now is ",z2    
    !
    PLANE_INDEX=[kpoint(1)-e,-1.0d0,2.0d0,3.0d0,1.0d0]
    !
    z2=z2-simps2d(chern_flux_integrable,[a-e,a+e],[b-e,b+e],N0=200,iterative=.false.)   
    write(*,*) "Done face 6, Berry flux now is ",z2 
    !
    write(*,*) "Chirailty is ",z2/(2*pi)
    unit=free_unit()
    open(unit,file="Chirality.ed")
    write(unit,*) z2/(2*pi)
    close(unit)
  end subroutine chern_retriever


  !subroutine chern_retriever(kpoint)   !integrates Berry connection on a path round kpoint (Stokes theorem) 
    !real(8)                         :: z2
    !real(8),dimension(3)            :: kpoint
    !real(8)                         :: e,a,b
    !integer                         :: unit
    !e=pi/50
    !write(*,*) "Testing chirality"
    !write(*,*) "Integrating on a square around ",kpoint
    !z2=0
    !!
    !a=kpoint(1)
    !!
    !PLANE_INDEX=[kpoint(2)-e,kpoint(3),1.0d0,1.0d0,2.0d0]   ![fixed component on the square, fixed component overall, versor,
                                                            !!varying varying component index, fixed component index]
    !!
    !z2=z2-simps(chern_flux_integrable, a-e, a+e, 200)   
    !write(*,*) "Done side 1, Berry flux now is ",z2 
    !!
    !a=kpoint(1)
    !!
    !PLANE_INDEX=[kpoint(2)+e,kpoint(3),-1.0d0,1.0d0,2.0d0]
    !!
    !z2=z2-simps(chern_flux_integrable, a-e, a+e, 200)   
    !write(*,*) "Done side 2, Berry flux now is ",z2 
    !!
    !a=kpoint(2)
    !!
    !PLANE_INDEX=[kpoint(1)-e,kpoint(3),-1.0d0,2.0d0,1.0d0]
    !!
    !z2=z2-simps(chern_flux_integrable, a-e, a+e, 200)   
    !write(*,*) "Done side 3, Berry flux now is ",z2 
    !!
    !a=kpoint(2)
    !!
    !PLANE_INDEX=[kpoint(1)+e,kpoint(3),1.0d0,2.0d0,1.0d0]
    !!
    !z2=z2-simps(chern_flux_integrable, a-e, a+e, 200)   
    !write(*,*) "Done side 4, Berry flux now is ",z2 
    !!
    !!
    !z2=z2/(2*pi)
    !if (z2<0.01d0) then
      !z2=0.0d0
    !end if
    !write(*,*) "Chirailty is ",z2
    !unit=free_unit()
    !open(unit,file="Chirality.ed")
    !write(unit,*) z2
    !close(unit)
  !end subroutine chern_retriever

  function chern_flux_integrable(kpoint_) result(flux) !gives scalar product between Berry flux and normal vector; takes 2d input,
    real(8),dimension(3)                     :: kpoint
    real(8),dimension(:)                     :: kpoint_
    real(8)                                  :: flux
    real(8),dimension(3)                     :: curvature 
    real(8),dimension(3)                     :: normal
    integer                                  :: i,j,k
    !
    normal=[0.0d0,0.0d0,0.0d0]
    kpoint=[0.0d0,0.0d0,0.0d0]
    i=IDINT(PLANE_INDEX(3))
    j=IDINT(PLANE_INDEX(4))
    k=IDINT(PLANE_INDEX(5))
    normal(k)=PLANE_INDEX(2)
    kpoint(i)=kpoint_(1)
    kpoint(j)=kpoint_(2)
    kpoint(k)=PLANE_INDEX(1)
    call get_Berry_Curvature(kpoint,curvature)
    flux=dot_product(curvature,normal)  
  end function chern_flux_integrable




  !function chern_flux_integrable(kpoint_) result(flux) !gives scalar product between Berry flux and normal vector; takes 2d input,
    !real(8),dimension(3)                     :: kpoint
    !real(8)                                  :: kpoint_
    !real(8)                                  :: flux
    !real(8),dimension(3)                     :: connection 
    !real(8),dimension(3)                     :: normal
    !integer                                  :: versor,varying_comp,fixed_comp
    !!
    !normal=[0.0d0,0.0d0,0.0d0]
    !kpoint=[0.0d0,0.0d0,0.0d0]
    !versor=IDINT(PLANE_INDEX(3))
    !varying_comp=IDINT(PLANE_INDEX(4))
    !fixed_comp=IDINT(PLANE_INDEX(5))
    !normal(varying_comp)=versor
    !kpoint(varying_comp)=kpoint_
    !kpoint(fixed_comp)=PLANE_INDEX(1)
    !kpoint(3)=PLANE_INDEX(2)
    !connection=get_Berry_Connection(kpoint)
    !flux=dot_product(connection,normal)  
  !end function chern_flux_integrable






  function get_occupied_state(kpoint,M) result(BlochStates) !takes nth occupied state of Hamiltonian in kpoint (n is passed by
    real(8),dimension(:),intent(in)        :: kpoint
    !
    integer                                :: M
    complex (8),dimension(M,M)             :: Eigvec ![Nlso][Nlso]
    real(8),dimension(M)                   :: Eigval ![Nlso]
    complex(8),dimension(M)                :: BlochStates ![Nlso]
    !
    !
    Eigvec=hk_weyl(kpoint,M)
    call eigh(Eigvec,Eigval)
    BlochStates = Eigvec(:,ICOMP)
    if (REALPART(BlochStates(1))<0.0d0) then
      BlochStates=-1.0d0*BlochStates
    end if
  end function get_occupied_state




  function get_Berry_connection(kpoint) result(BerryConnection) !gives Berry Connection in kpoint by its definition (eg Turner 2013)
    real(8),dimension(:),intent(in)                    :: kpoint
    !
    integer                                            :: i
    real(8),dimension(size(kpoint))                    :: BerryConnection ![N_dimensions]    
    complex(8),dimension(Nso)                          :: BlochStates
    complex(8),dimension(Nso,size(kpoint))             :: TmpKmat
    !
    !
    BerryConnection=[0.0d0,0.0d0,0.0d0]
    do ICOMP=1,Norb
      BlochStates=get_occupied_state(kpoint,Nso)
      call djacobian_complex(kpoint,TmpKMat,Nso)
      do i = 1,size(kpoint)
        BerryConnection(i) = BerryConnection(i) - IMAGPART(dot_product(BlochStates, TmpKmat(:,i)))
      enddo
    enddo
  end function get_Berry_connection




  subroutine get_Berry_Curvature(kpoint,BerryCurvature) !Does curl to get Berry Curvature
    real(8),dimension(:)                      :: kpoint
    !
    real(8),dimension(:,:),allocatable        :: TmpKmat ![N_dimension][N_dimension]
    real(8),dimension(3)                      :: BerryCurvature ![N_dimensions]    
    !
    allocate(TmpKMat(size(kpoint),size(kpoint))) 
    !
    !Do the rotor to get Berry Curvature
    call djacobian(get_Berry_connection,kpoint,TmpKMat)
    BerryCurvature(1)=TmpKMat(3,2)-TmpKMat(2,3)
    BerryCurvature(2)=TmpKMat(1,3)-TmpKMat(3,1)
    BerryCurvature(3)=TmpKMat(2,1)-TmpKMat(1,2)
    !
    deallocate(TmpKMat)
  end subroutine get_Berry_Curvature





  subroutine djacobian_complex(kpoint,jacobian,M)
    real(8),dimension(:)                        :: kpoint
    real(8),dimension(size(kpoint))             :: kpoint1,kpoint2,delta
    complex(8),dimension(M)                     :: eigenstate1,eigenstate2
    complex(8),dimension(M,size(kpoint))        :: jacobian
    integer                                     :: i,j,M
    !
    do j=1,size(kpoint)
      delta=zero
      delta(j)=0.001d0
      kpoint1=kpoint+delta
      kpoint2=kpoint-delta
      eigenstate1=get_occupied_state(kpoint1,M)
      eigenstate2=get_occupied_state(kpoint2,M)
      jacobian(:,j)=(eigenstate1-eigenstate2)/(2*delta(j))
    enddo
  end subroutine djacobian_complex

  



  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!
  subroutine build_z2_indices(sigma0)
    integer                                :: unit
    complex(8),dimension(Nso,Nso),optional :: sigma0(Nso,Nso)
    complex(8),dimension(Nso,Nso)          :: sigma0_
    !
    integer,dimension(4)                   :: z2
    !
    sigma0_=zero;if(present(sigma0))sigma0_=sigma0
    !Evaluate the Z2 index:
    !STRONG TI
    z2(1) = z2_number(reshape( [ [0,0,0] , [1,0,0] , [1,1,0] , [0,1,0] , [0,1,1] , [0,0,1] , [1,0,1] , [1,1,1] ] , [3,8] )*pi,sigma0_)
    !WEAK TI
    !K=1: n_1=1, n_2,3=0,1
    z2(2) = z2_number(reshape( [ [1,0,0] , [1,1,0] , [1,1,1] , [1,0,1] ] , [3,4])*pi,sigma0_)
    !K=2: n_2=1, n_1,2=0,1
    z2(3) = z2_number(reshape( [ [0,1,0] , [0,1,1] , [1,1,1] , [1,1,0] ] , [3,4])*pi,sigma0_)
    !k=3: n_3=1, n_1,2=0,1
    z2(4) = z2_number(reshape( [ [0,0,1] , [0,1,1] , [1,1,1] , [1,0,1] ] , [3,4])*pi,sigma0_)
    unit=free_unit()
    open(unit,file="z2_invariant.ed")
    write(unit,*)z2
    close(unit)
  end subroutine build_z2_indices

  function z2_number(ktrims,sigma0) result(z2)
    real(8),dimension(:,:),intent(in)            :: ktrims
    complex(8),dimension(Nso,Nso)                :: sigma0
    complex(8),dimension(Nso,Nso,size(Ktrims,2)) :: Htrims
    complex(8),dimension(size(Ktrims,2))         :: Delta
    integer                                      :: z2
    integer                                      :: i,j,Ntrim,itrim,Nocc
    !
    Ntrim=size(Ktrims,2)
    !
    do itrim=1,Ntrim
       Htrims(:,:,itrim) = hk_weyl(Ktrims(:,itrim),Nso) + sigma0(:,:)
       Delta(itrim)=-sign(1d0,dreal(Htrims(1,1,itrim)))
    enddo
    !
    z2=product(Delta(:))
    if(z2>0)then
       z2=0
    else
       z2=1
    end if
    !
  end function z2_number



  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index


  function so2j(fg,Nso) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nso,Nso)               :: g
    integer                                     :: Nso,i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function so2j

  function j2so(fg) result(g)
    complex(8),dimension(Nso,Nso)               :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2so


end program ed_wsm_3d
