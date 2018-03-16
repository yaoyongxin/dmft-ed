program ed_wsm_3d
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer                       :: iloop,Lk,Nso
  logical                       :: converged,small_error,sigma_symmetric
  !Bath:
  integer                       :: Nb
  real(8),allocatable           :: Bath(:),Bath_Prev(:)
  !The local hybridization function:
  complex(8),allocatable        :: Weiss(:,:,:,:,:)
  complex(8),allocatable        :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable        :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable        :: Hk(:,:,:),wsmHloc(:,:),sigmaWSM(:,:),Zmats(:,:)
  real(8),allocatable           :: Wtk(:)
  integer,allocatable           :: ik2ix(:),ik2iy(:),ik2iz(:)
  !variables for the model:
  integer                       :: Nk,Nkpath
  real(8)                       :: e0,mh,sigma_difference,sigma_symmetry,lambda,bx,by,bz,BIA
  real(8)                       :: wmixing
  character(len=16)             :: finput
  character(len=32)             :: hkfile
  logical                       :: orbsym,getpoles,usez
  complex(8),dimension(4,4)     :: emat,soxmat,soymat,sozmat,bxmat,bymat,bzmat,BIAmat
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED_WSM.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=30)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=1d0)
  call parse_input_variable(sigma_symmetry,"SIGMA_SYMMETRY",finput,default=10.0d0)
  call parse_input_variable(e0,"E0",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.5d0)
  call parse_input_variable(bx,"BX",finput,default=0.3d0)
  call parse_input_variable(by,"BY",finput,default=0.d0)
  call parse_input_variable(bz,"BZ",finput,default=0.d0)
  call parse_input_variable(BIA,"BIA",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(orbsym,"ORBSYM",finput,default=.false.)
  call parse_input_variable(usez,"USEZ",finput,default=.false.)
  !
  !READ INPUT FILE
  !
  call ed_read_input(trim(finput))
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  !
  !SETUP THE GAMMA MATRICES:
  emat = kron_pauli( pauli_sigma_0, pauli_tau_z)
  soxmat=kron_pauli( pauli_sigma_z, pauli_tau_x)
  soymat=kron_pauli( pauli_sigma_0, pauli_tau_y)
  sozmat=kron_pauli( pauli_sigma_x, pauli_tau_x)
  bxmat =kron_pauli( pauli_sigma_x, pauli_tau_z)
  bymat =kron_pauli( pauli_sigma_y, pauli_tau_0)
  bzmat =kron_pauli( pauli_sigma_z, pauli_tau_z)
  BIAmat=kron_pauli( pauli_sigma_y, pauli_tau_y)
  !
  !PRELIMINARY CHECKS
  !
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
  allocate(Zmats(Nso,Nso))
  !
  !SET sigma_WSM(0) to zero
  !
  call set_sigmaWSM()
  !
  !Buil the Hamiltonian on a grid or on path
  !
  call build_hk(trim(hkfile))
  !
  !Setup solver
  !
  Nb=get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_Prev(Nb))
  call ed_init_solver(bath,j2so(wsmHloc))
  !
  !DMFT loop
  !
  small_error=.false.                 ! flag for dmft error
  sigma_symmetric=.false.             ! flag for sigma symmetry breaking - usually huge aka disabled
  iloop=0;converged=.false.           ! loop end flags
  !
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath)
     !
     ! retrieve the self-energies
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)
     call build_z2_indices( so2j(Smats(:,:,:,:,1),Nso) ) !so2j transforms a Nspin:Nspin:Norb:Norb into a Nso:Nso matrix
     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=3)
     !
     ! compute the Weiss field (only the Nineq ones)
     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats,Smats,Weiss,Hloc=j2so(wsmHloc))
     else
        call dmft_delta(Gmats,Smats,Weiss,Hloc=j2so(wsmHloc))
     endif
     !
     !Fit the new bath, starting from the old bath + the supplied Weiss/Delta
     if(ed_mode=="normal")then
        call ed_chi2_fitgf(Weiss,bath,ispin=1)
        call ed_chi2_fitgf(Weiss,bath,ispin=2)
     else
        call ed_chi2_fitgf(Weiss,bath)
     endif
     !
     !if flag is set, symmetrize the bath
     !
     if(orbsym)then
     call copy_component_bath(-bath,1,1,bath,1,2,1)
     call copy_component_bath(-bath,2,1,bath,2,2,1)
     !
     call copy_component_bath(bath,1,1,bath,1,2,2)
     call copy_component_bath(bath,2,1,bath,2,2,2)
     endif
     !
     !MIXING the current bath with the previous:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_Prev
     Bath_Prev=Bath
     !
     !First check: DMFT error
     !
     small_error = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop)
     !
     !Second check: Sigma symmetry
     !
     sigma_difference=MAXVAL(Abs(Smats(1,1,1,1,:)+CONJG(Smats(1,1,2,2,:))))
     !
     if (sigma_difference .lt. sigma_symmetry) then
        sigma_symmetric=.true.
        write(*,"(A,ES15.7)")bold_green("Sigma symmetry is wrong by "),sigma_difference
     else
        sigma_symmetric=.false.
        write(*,"(A,ES15.7)")bold_red("Sigma symmetry is wrong by "),sigma_difference
     endif
     !
     !converge if all checks are positive
     !
     converged=small_error.AND.sigma_symmetric
     !
     call end_loop
     !
  enddo
  !
  ! compute the local gf:
  call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=3)
  !
  !Get kinetic energy:
  call dmft_kinetic_energy(Hk,Wtk,Smats)
  !
  !Get 3d Bands from Top. Hamiltonian
  call solve_hk_topological( so2j(Smats(:,:,:,:,1),Nso) )
  !
  !Find out if it is a semimetal
  call is_weyl( Nso, so2j(Smats(:,:,:,:,1),Nso) )
 
 


contains



  !--------------------------------------------------------------------!
  !PURPOSE: Define the WSM Hamiltonian
  !--------------------------------------------------------------------!
  !
  function hk_weyl(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N,i
    real(8)                   :: kx,ky,kz
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_weyl: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    kz=kpoint(3)
    !
    Hk          = zero
    !
    Hk=(Mh - e0*(cos(kx) + cos(ky) + cos(kz)))*emat+&
        lambda*(sin(kx)*soxmat - sin(ky)*soymat + sin(kz)*sozmat)+&
        BIA*BIAmat + bx*bxmat + by*bymat + bz*bzmat
    !
    !
    !add the SigmaWSM term to get Topologial Hamiltonian if required:
    Hk = Hk + dreal(SigmaWSM)
    !
    if (usez) then
      Zmats=zero
      do i=1,Nso
        Zmats(i,i)  = 1.d0/abs( 1.d0 -  dimag(sigmaWSM(i,i))/(pi/beta) )
      end do
      Hk = matmul(Zmats,Hk)
    endif 
    !
  end function hk_weyl


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
    !GET LOCAL PART OF THE HAMILTONIAN
    if(allocated(wsmHloc))deallocate(wsmHloc)
    allocate(wsmHloc(Nso,Nso))
    wsmHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(wsmHloc))<1.d-9)wsmHloc=0d0
    call TB_write_Hloc(wsmHloc)
    !
  end subroutine build_hk


  !--------------------------------------------------------------------!
  !PURPOSE: Set the Self-Energy
  !--------------------------------------------------------------------!
  !
  subroutine set_SigmaWSM(sigma)
    complex(8),dimension(Nso,Nso),optional :: sigma(Nso,Nso)
    sigmaWSM = zero;if(present(sigma))sigmaWSM=sigma
  end subroutine set_SigmaWSM
  


  !--------------------------------------------------------------------!
  !PURPOSE: Solve the topological Hamiltonian
  !--------------------------------------------------------------------!
  
  
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
         points_name=[charactClaudio er(len=10) :: "G","X1","G","X","G","Z","G","M","G","A1","G","A","G","R"],&
         file="Eig_Htop.ed")
    if (usez) then
      write(*,*) "Z11=",Zmats(1,1)
      write(*,*) "Z22=",Zmats(2,2)
      write(*,*) "Z33=",Zmats(3,3)
      write(*,*) "Z44=",Zmats(4,4)
    endif
  end subroutine solve_hk_topological




  !-------------------------------------------------------------------------------!
  !PURPOSE: Adaptive routine for Weyl point finder
  !-------------------------------------------------------------------------------!
  !
  subroutine is_weyl(N,sigma)
    integer                                 :: i,j,k,n_found
    real(8)                                 :: halfwidth,tmp,e_low
    complex(8),dimension(Nso,Nso)           :: sigma(Nso,Nso)
    real(8),dimension(4)                    :: diff_array
    integer                                 :: N,unit
    real(8),dimension(100,4)                :: possible_weyls
!preliminary checks
    if(N/=4)stop "hk_weyl: error in N dimensions"
    if(.not. usez) then
      write(*,*) "WARNING: It is better to renormalize Htop (USEZ=true), setting it automatically."
      usez=.true.
    endif
!
    call set_sigmaWSM(sigma)
!
    possible_weyls=zero
!initialize a fake Weyl point in Gamma
    possible_weyls(1,:)=[0.0d0,0.0d0,0.0d0,1.0d0]
!first step cube is all the BZ
    halfwidth=pi
!initial energy threshold
    e_low=0.1d0
!Main loop,10 steps should be enough
    write(*,*) "Starting Weyl point search"
    convloop: do i=1,10
      write(*,*) " "
      write(*,*) "step ", i
      write(*,*) "cube size is ",2*halfwidth
      write(*,*) "Energy threshold is ",e_low
!check cube around each possible weyl point
      call grid_checker(possible_weyls,halfwidth,e_low,n_found)
!the new cube half-size is half the minimum distance of the found WPs
      do j=1,n_found-1
        do k=j+1,n_found
          diff_array=abs(possible_weyls(j,:)-possible_weyls(k,:))
          tmp=minval(diff_array, MASK=(diff_array > 0.000001d0))
          if (tmp < halfwidth) then
            halfwidth=0.5*tmp   !!!!!should I multiply by 0.4?
          endif
        enddo
      enddo
!if something is strange default to a decreasing function of the step
      if (halfwidth > pi/(i**3)) then
        halfwidth=(pi/i**3)
      endif
!exit from main loop if decent convergences is achieved or no points are left
      if (e_low < 0.000001 .or. n_found .eq. 0) then
        exit convloop
      endif
    enddo convloop 
    write(*,*) "Ended Weyl point search"
    !
!print nonzero chirality points
    unit=free_unit()
    open(unit,file="Chirality.ed",position="append")
    do i=1,100
      if(abs(possible_weyls(i,4)) > 0.0d0) then
        write(unit,'(4F16.6)') possible_weyls(i,:)
      endif
    enddo
    close(unit)
  end subroutine is_weyl


!--------------------------------------------------------------------!
!Purpose: build the cubes in which to evaluate chiralities
!--------------------------------------------------------------------!

  subroutine grid_checker(weyl_centers,cubehalf_size,E_thresh,w)
    real(8)                           :: cubehalf_size,step,E_thresh,E_thresh_tmp
    real(8),dimension(3)              :: kpoint
    real(8),dimension(100,4)           :: weyl_centers,weyl_centers_tmp
    integer                           :: i,j,k,p,w,mash_thickness
    complex(8),dimension(4,4)         :: ham
    real(8),dimension(4)              :: eval
!Discretize the cube    
    mash_thickness=2
    step=2*cubehalf_size/(mash_thickness*Nk)
    write(*,*) "kstep is ",step
    weyl_centers_tmp=zero
!Number of new Weyls
    w=1
!Check cubes around each input Weyl point
    weylloop: do p=1,100
!If there are too many weyl points stop
      if (w .gt. 100) then
        write(*,*) "too many weyls"
        exit weylloop
      endif
!Only do cubes around Weyls
      if (abs(weyl_centers(p,4)) > 0.6d0) then
        write(*,*)"doing around",weyl_centers(p,:)
!Dynamic threshold for the current cube
        E_thresh_tmp=E_thresh
!K-mash in the cube
        xloop: do i=0,mash_thickness*Nk-1
          kpoint(1)=weyl_centers(p,1) - cubehalf_size + step*i
          yloop: do j=0,mash_thickness*Nk-1
            kpoint(2)=weyl_centers(p,2) - cubehalf_size + step*j
            zloop: do k=0,mash_thickness*Nk-1
              kpoint(3)=weyl_centers(p,3) - cubehalf_size + step*k
!Eigenvalues and gap
              ham=hk_weyl(kpoint,4)
              call eigh(ham,Eval)
!If gap is smaller than the dynamic threshold
              if (abs(Eval(3) - Eval(2)) .lt. E_thresh_tmp) then 
!Decrease threshold (leaving some room upwards)
                E_thresh_tmp=MINVAL([1.1*(abs(Eval(3) - Eval(2))),E_thresh])
!Calculate chirality
                call chern_retriever(kpoint,step,weyl_centers_tmp(w,4))
!If it's Weyl add it and increase count
                if (abs(weyl_centers_tmp(w,4))>0.6) then
                  write(*,*)"found point ",w
                  weyl_centers_tmp(w,1)=kpoint(1)
                  weyl_centers_tmp(w,2)=kpoint(2)
                  weyl_centers_tmp(w,3)=kpoint(3)
                  w=w+1
                end if
              end if
            end do zloop
          end do yloop
        end do xloop
      endif
    end do weylloop
!Update quantities, w is actually w-1, E_thresh is indicative
    w=w-1
    E_thresh=E_thresh_tmp
!Write stuff
    write(*,*)"Number of points found: ",w
    weyl_centers=weyl_centers_tmp
    do i=1,w+1
      write(*,*) weyl_centers(i,:)
    end do
    write(*,*) "New energy threshold=",E_thresh
  end subroutine grid_checker

  !--------------------------------------------------------------------!
  !PURPOSE: Evaluation of the Chirality of the candidate Weyl point
  !--------------------------------------------------------------------!
  !
  subroutine chern_retriever(kpoint,e,z2)   !integrates Berry flux on a cubic surface around the point kpoint
    real(8)                         :: z2
    complex(8)                      :: phase
    real(8),dimension(3)            :: kpoint,kpoint_
    real(8)                         :: e,cubesize
    integer                         :: perm,j,q,N,face_indx,face_sign,side_indx,side_sign,path_indx,border,varying_indx,run_direction
    integer,dimension(6,6)          :: permutations
    integer,dimension(4)            :: versor
    complex(8),dimension(4,4)       :: BlochOld,BlochNew
    complex(8),dimension(2,2)       :: OverlapMatrix
    real(8),dimension(4)            :: Eigval
    !write(*,*) "Cubesize is",2*e
    N=500
    !write(*,*) "Testing chirality"
    !write(*,*) "Integrating on a cube around ",kpoint
    z2=0.0d0
    !
    permutations(1,:)=[ 3, 1, 2,-1,-2, 1] ![index fixing the face, 4*index fixing the path side, coordinate varying on side labelled by #5]
    permutations(2,:)=[-3,-2,-1, 2, 1, 2] 
    permutations(3,:)=[ 1, 2, 3,-2,-3, 2]
    permutations(4,:)=[-1,-3,-2, 3, 2, 3]
    permutations(5,:)=[ 2, 3, 1,-3,-1, 3]
    permutations(6,:)=[-2,-1,-3, 1, 3, 1]
    versor=[1,-1,-1,1]
    !
    do perm=1,6
      face_indx=ABS(permutations(perm,1))   !which plane am I parallel to?
      face_sign=SIGN(1,permutations(perm,1)) !top-botton, left-right face, now I have selected one
      phase=1.0d0   !do path around this face
      !
      kpoint_=zero  
      BlochOld=zero
      BlochNew=zero
      OverlapMatrix=zero
      !
      kpoint_(face_indx)=kpoint(face_indx)+e*face_sign !fix one coordinate of the running point
      !
      do border=2,5 !run along border
        !
        side_indx=ABS(permutations(perm,border))      !which direction am I parallel to?
        side_sign=SIGN(1,permutations(perm,border))   !top-bottom, left-right side, now I've selected one
        varying_indx=ABS(permutations(perm,border+1)) !the point runs along the other coordinate
        run_direction=face_sign*versor(border-1)                !direction to run along to have a closed path
        !
        kpoint_(side_indx)=kpoint(side_indx)+e*side_sign !fix another coordinate of the running point
        !
        do path_indx=1,2*N !run on one side of the border
          !
          if (all(BlochOld(:,1) .eq. 0.0d0)) then !if I have yet to start, do one cycle
            kpoint_(varying_indx)=kpoint(varying_indx)-e*run_direction
            BlochOld=hk_weyl(kpoint_,4)
            call eigh(BlochOld,Eigval)  !store old eigenvectors
          endif
          kpoint_(varying_indx)=kpoint(varying_indx)-e*run_direction+(e*path_indx/N)*run_direction
          BlochNew=hk_weyl(kpoint_,4)
          call eigh(BlochNew,Eigval)  !store new eigenvectors
          do j=1,2
            do q=1,2
              OverlapMatrix(j,q)=dot_product(BlochOld(:,j),BlochNew(:,q))
            enddo
          enddo
          BlochOld=BlochNew !update eigenvectors
          phase=phase*det(OverlapMatrix)
        enddo !end run on one side of the border
      enddo !end run on the border
      z2=z2+IMAG(log(phase))!  !sum the phase to the face
    enddo !end run on face
    if (abs(z2)>0.6d0)then  !if it's a point, append to Chirality file
      z2=SIGN(1.0d0,z2/(2*pi))
    else
      z2=0.0d0
    endif
  end subroutine chern_retriever

  !--------------------------------------------------------------------!
  !PURPOSE: TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
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

