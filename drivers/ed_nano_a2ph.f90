program ed_nano_adiabatic
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  abstract  interface
     subroutine drive_template(Vijt,Uijt,Hij0,Hijt,time_step,N)
        integer                   :: N
        real(8)                   :: time_step
        complex(8),dimension(N,N) :: Vijt !drive
        complex(8),dimension(N,N) :: Uijt !drive time derivative
        complex(8),dimension(N,N) :: Hij0 !time-independent Hamiltonian
        complex(8),dimension(N,N) :: Hijt !time-  dependent Hamiltonian
     end subroutine drive_template
  end interface  



  integer                                         :: iloop
  logical                                         :: converged
  integer                                         :: ilat,ineq,ispin,iorb
  !bath:
  integer                                         :: Nb
  real(8),allocatable                             :: Bath_prev(:,:),Bath_ineq(:,:)
  !local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Sreal,Sreal_ineq ![Nlat*(Nspin*Norb)**2*Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Greal,Greal_ineq
  real(8), allocatable,dimension(:)               :: dens,dens_ineq
  real(8), allocatable,dimension(:)               :: docc,docc_ineq
  !hamiltonian input:
  complex(8),allocatable                          :: Hij(:,:,:) ![Nlso][Nlso][Nk]
  complex(8),allocatable                          :: Hij_static(:,:,:) ![Nlso][Nlso][Nk]
  complex(8),allocatable                          :: nanoHloc(:,:),Hloc(:,:,:,:,:),Hloc_ineq(:,:,:,:,:)
  integer                                         :: Nk,Nlso,Nlo,Nineq,Nlat
  integer,dimension(:),allocatable                :: lat2ineq,ineq2lat
  integer,dimension(:),allocatable                :: sb_field_sign
  !
  real(8)                                         :: wmixing,Eout(2)
  !
  !input files:
  character(len=32)                               :: finput
  character(len=32)                               :: nfile,hijfile
  !
  logical                                         :: phsym
  logical                                         :: leads
  !
  !non-local Green's function:
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gijreal
  !
  !hybridization function to environment
  complex(8),dimension(:,:,:),allocatable         :: Hyb_mats ![Nlso][Nlso][Lmats]
  complex(8),dimension(:,:,:),allocatable         :: Hyb_real ![Nlso][Nlso][Lreal]
  !
  !drive time/frequency variables
  integer                                         :: itime,Ltime
  real(8),dimension(:),allocatable                :: time ![Ltime]
  complex(8),dimension(:,:),allocatable           :: Vij,Uij ![Nlso][Nlso]
  !
  procedure(drive_template),pointer               :: drive_model
  character(len=32)                               :: drive
  !
  integer                                         :: unit
  logical                                         :: exist,noint,igetgf


  call parse_cmd_variable(finput,"FINPUT",default='inputED_NANO.conf')
  call parse_input_variable(nfile,"NFILE",finput,default="nano.in")
  call parse_input_variable(hijfile,"HIJFILE",finput,default="hij.in")
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(phsym,"phsym",finput,default=.false.)
  ! parse environment & transport flags
  call parse_input_variable(leads,"leads",finput,default=.false.)
  ! parse time variable
  call parse_input_variable(Ltime,"LTIME",finput,default=100)
  !
  call parse_input_variable(noint,"NOINT",finput,default=.true.)
  call parse_input_variable(igetgf,"IGETGF",finput,default=.false.)

  ! parse drive model
  call parse_input_variable(drive,"DRIVE",finput,default='default')
   

  select case(trim(drive))
     case ('angle')
        drive_model => angle
     case ('flux')
        drive_model => flux
     case ('pump')
        drive_model => pump
     case default
        write(*,*) "error: select a drive"
        stop
  end select






  ! read input
  call ed_read_input(trim(finput))


  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  ! set input structure hamiltonian
  call build_Hij([nfile,hijfile])
  ! store static hamiltonian
  Hij_static=Hij

  ! allocate weiss field:
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  ! allocate self-energy
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  ! allocate Green's function
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gijreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  ! allocate Hloc
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb))
  ! allocate density
  allocate(dens(Nlat))
  allocate(dens_ineq(Nineq))
  ! allocate double occupations
  allocate(docc(Nlat))
  allocate(docc_ineq(Nineq))

  !Hloc = reshape_Hloc(nanoHloc,Nlat,Nspin,Norb)
  Hloc = lso2nnn_reshape(nanoHloc,Nlat,Nspin,Norb)


  ! change suffix for prints
  call set_gf_suffix(".ed")


  ! allocate hybridization matrix
  if(leads)then
     call set_hyb()
  endif


  ! allocate effective time [0:pi2] 
  allocate(time(Ltime))
  if(Ltime==1)then
     time=0.d0
  else
     time = linspace(0.d0,pi2,Ltime)
  endif
  ! write interval
  !do itime=1,Ltime
  !   write(*,*)itime,time(itime)
  !enddo


  ! allocate drive
  allocate(Vij(Nlso,Nlso),Uij(Nlso,Nlso))


  ! setup solver
  Nb=get_bath_dimension()
  allocate(Bath_ineq(Nineq,Nb))
  allocate(Bath_prev(Nineq,Nb))



  Sreal=zero
  !call dmft_gij_realaxis(Hij,[1d0],Gijreal,Sreal)
  !call dmft_print_gij_realaxis(Gijreal,"Gij_static",iprint=2)


  !start time loop
  do itime=1,Ltime

     ! set drive & its derivative, overwrite Hij to include drive
     call set_drive(drive_model,Vij,Uij,Hij_static,Hij,time(itime))

     
     ! extract nanoHloc and reshape to get Hloc at the actual time
     nanoHloc = extract_Hloc(Hij,Nlat,Nspin,Norb)
     Hloc = lso2nnn_reshape(nanoHloc,Nlat,Nspin,Norb)
     !
     ! set Hloc for each inequivalent site
     do ineq=1,Nineq
        ilat = ineq2lat(ineq)
        Hloc_ineq(ineq,:,:,:,:) = Hloc(ilat,:,:,:,:)
     enddo
  
     ! ----------------------------------------------------------------- 
     ! non-interacting case
     ! ----------------------------------------------------------------- 
     if(noint)then
        Sreal=zero
        !if(itime==1)then 
        !  call dmft_gij_realaxis(Hij,[1d0],Gijreal,Sreal)
        !  call dmft_print_gij_realaxis(Gijreal,"Gij_drive",iprint=2)
        !endif
     else

     ! ----------------------------------------------------------------- 
     ! perform L\'anczos for interacting system
     ! ----------------------------------------------------------------- 
     !
     ! init solver
     call ed_init_solver(Bath_ineq,Hloc_ineq)
     
     ! break SU(2) symmetry for magnetic solutions
     do ineq=1,Nineq
        ilat = ineq2lat(ineq)
        if(Nspin>1) call break_symmetry_bath(Bath_ineq(ineq,:),sb_field,dble(sb_field_sign(ineq)))
     enddo
   
     ! dmft loop
     iloop=0 ; converged=.false.
     do while(.not.converged.AND.iloop<nloop) 
        iloop=iloop+1
        call start_loop(iloop,nloop,"DMFT-loop")   
        bath_prev=bath_ineq
   
        ! solve impurities on each inequivalent site:
        call ed_solve(bath_ineq,Hloc_ineq)
   
        ! retrieve self-energies and occupations(Nineq,Norb=1)
        call ed_get_sigma_matsubara(Smats_ineq,Nineq)
        call ed_get_sigma_real(Sreal_ineq,Nineq)
        call ed_get_dens(dens_ineq,Nineq,iorb=1)
        call ed_get_docc(docc_ineq,Nineq,iorb=1)
   
        ! spread self-energies and occupation to all lattice sites
        do ilat=1,Nlat
           ineq = lat2ineq(ilat)
           dens(ilat) = dens_ineq(ineq)
           docc(ilat) = docc_ineq(ineq)
           Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
           Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
        enddo
   
        ! compute the local gf:
        !
        if(leads)then
           call dmft_set_Gamma_matsubara(hyb_mats)
        endif
        call dmft_gloc_matsubara(Hij,[1d0],Gmats,Smats)
        call dmft_print_gf_matsubara(Gmats,"LG",iprint=1)
        do ineq=1,Nineq
           ilat = ineq2lat(ineq)
           Gmats_ineq(ineq,:,:,:,:,:) = Gmats(ilat,:,:,:,:,:)
        enddo
        !
        if(leads)then
           call dmft_set_Gamma_realaxis(hyb_real)
        endif
        call dmft_gloc_realaxis(Hij,[1d0],Greal,Sreal)
        call dmft_print_gf_realaxis(Greal,"LG",iprint=1)
        do ineq=1,Nineq
           ilat = ineq2lat(ineq)
           Greal_ineq(ineq,:,:,:,:,:) = Greal(ilat,:,:,:,:,:)
        enddo
   
        ! compute the Weiss field
        if(cg_scheme=="weiss")then
           call dmft_weiss(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq)
        else
           call dmft_delta(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq)
        endif
   
        ! fit baths and mix result with old baths
        do ispin=1,Nspin
           call ed_chi2_fitgf(bath_ineq,Weiss_ineq,Hloc_ineq,ispin)
        enddo
   
        if(phsym)then
           do ineq=1,Nineq
              call ph_symmetrize_bath(bath_ineq(ineq,:),save=.true.)
           enddo
        endif
        Bath_ineq=wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
   
        converged = check_convergence(Weiss_ineq(1,1,1,1,1,:),dmft_error,nsuccess,nloop)
        ! alternative convergency criteria
        !converged = check_convergence_local(docc_ineq,dmft_error,nsuccess,nloop)
        if(NREAD/=0.d0) call search_chemical_potential(xmu,sum(dens)/Nlat,converged)
   
        call end_loop()
     end do
   
     ! save self-energy on disk
     call dmft_print_gf_matsubara(Smats_ineq,"LSigma",iprint=1)
     call dmft_print_gf_realaxis(Sreal_ineq,"LSigma",iprint=1)
   
     ! compute kinetic energy at convergence
     !call dmft_kinetic_energy(Hij,[1d0],Smats)
   

     ! get occupations at convergence and write as a function of time
     call ed_get_dens(dens_ineq,Nineq,iorb=1)
     unit = free_unit()
     inquire(file="observables_time.ed",exist=exist)
     if(exist)then
       open(unit,file="observables_time.ed",status="old",position="append")
     else
       open(unit,file="observables_time.ed")
     endif
     write(unit,'(1f16.9)',advance='no')time(itime)
     do ineq=1,Nineq
        write(unit,'(1f16.9)',advance='no')dens_ineq(ineq)
     enddo
     write(unit,*) !newline
     close(unit)

 
     endif
     ! ----------------------------------------------------------------- 


     ! ----------------------------------------------------------------- 
     ! set embedding matrix
     ! ----------------------------------------------------------------- 
     if(leads)then
        call dmft_set_Gamma_realaxis(hyb_real) !needed for dmft_gij_realaxis
     endif

     ! ----------------------------------------------------------------- 
     ! evalueate and print Green's function
     ! ----------------------------------------------------------------- 
     call dmft_gloc_realaxis(Hij,[1d0],Greal,Sreal)
     if(igetgf)then
        call dmft_print_gf_realaxis(Greal,"LG",iprint=1)
     endif
   
     ! ----------------------------------------------------------------- 
     ! evaluate DC/AC transport properties
     ! ----------------------------------------------------------------- 
     !
     ! evaluates the non-local Green's function
     if(.not.allocated(Gijreal))then
        allocate(Gijreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
     endif
     call dmft_gij_realaxis(Hij,[1d0],Gijreal,Sreal)
     !
     ! evaluate transmission function, DC & AC current
     call ed_transport_acdc(Gijreal,time(itime))
     !
     ! deallocate non-local Green's function
     deallocate(Gijreal)

     !stop

  enddo ! end of time loop
   


contains



  !----------------------------------------------------------------------------------------!
  ! purpose: build real-space Hamiltonian for a nanostructure of size [Nlat*Nspin*Norb]**2
  !----------------------------------------------------------------------------------------!
  subroutine build_Hij(file)
    character(len=*)     :: file(2)
    integer              :: ilat,jlat,iorb,jorb,is,js,ispin,ie
    integer              :: i,isite,iineq,iineq0,isign
    integer              :: EOF
    character, parameter :: tab = achar ( 9 )
    integer              :: unit,ineq_count
    integer              :: Ns,Ne,Nb,Nk         ! #atoms, #inequivalent, #bands
    real(8)              :: ret,imt
    logical              :: blank_at_right
    character(len=1)     :: next,prev
    character(len=6)     :: site,sign
    write(LOGfile,*)"Build H(R_i,R_j) for a NANO object:"
    ! readin generic input
    ! allocate & fill inequivalent list
    unit = free_unit()
    open(unit,file=trim(file(1)),status='old')
    read(unit,*)Ns,Ne,Nb
    !Checks:
    if(Nb/=Norb)stop "build_Hij error: Nb read from file != Norb in input.conf"
    Nk   = 1
    Nb   = Norb
    Nlat = Ns
    Nineq= Ne
    Nlso = Nlat*Nspin*Norb
    allocate(lat2ineq(Nlat),ineq2lat(Nineq))
    read(unit,"(A1)",advance='no',IOSTAT=EOF)next
    site  = next
    isite = 0
    i     = 0
    do 
       prev=next
       read(unit,"(A1)",advance='no',IOSTAT=EOF)next
       blank_at_right = ((prev/=' '.AND.prev/=tab).AND.(next==' '.OR.next==tab))
       if(.not.blank_at_right)then
          site=trim(site)//next
       else
          read(site,"(I6)")isite
          site=""
          i=i+1
          if(i>Nlat)stop "build_Hij error: lattice index > Nlat read from file"
          lat2ineq(i)=isite+1
       endif
       if(EOF<0)exit
    enddo
    if(i<Nlat)stop "build_Hij error: lattice index < Nlat read from file"
    write(*,*)"# of sites      :",Nlat
    write(*,*)"# of ineq sites :",Nineq
    write(*,*)"# of bands      :",Norb
    !
    ineq_count=1
    iineq=lat2ineq(Nlat)
    do i=Nlat,2,-1
       iineq0=lat2ineq(i-1)!iineq
       iineq =lat2ineq(i)
       if(iineq/=iineq0)then
          ineq2lat(iineq)=i
          ineq_count=ineq_count+1
       endif
       !if(ineq_count==Nineq)exit
    enddo
    iineq=lat2ineq(1)
    ineq2lat(1)=iineq
    !close(unit) ! do not close unit if readin info below
    !
    ! allocate & fill sign list of symmetry-breaking field
    allocate(sb_field_sign(Nineq))
    sign  = next
    isign = 0
    i     = 0
    do 
       prev=next
       read(unit,"(A1)",advance='no',IOSTAT=EOF)next
       blank_at_right = ((prev/=' '.AND.prev/=tab).AND.(next==' '.OR.next==tab))
       if(.not.blank_at_right)then
          sign=trim(sign)//next
       else
          read(sign,"(I6)")isign
          sign=""
          i=i+1
          if(i>Nineq)stop "build_Hij error: lattice index > Nineq read from file"
          sb_field_sign(i)=isign
       endif
       if(EOF<0)exit
    enddo
    close(unit)
    !
    ! allocate and initialize H(r_i,r_j)
    allocate(Hij(Nlso,Nlso,Nk))
    Hij = zero 
    unit = free_unit()
    open(unit,file=trim(file(2)),status='old')
    do !while(EOF>=0)
       read(unit,*,IOSTAT=EOF)ilat,iorb,jlat,jorb,ret,imt
       ilat=ilat+1
       iorb=iorb+1
       jlat=jlat+1
       jorb=jorb+1
       if(EOF<0)exit
       do ispin=1,Nspin
          is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
          js = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb
          ! symmetric hopping
          Hij(is,js,1)=dcmplx(ret, imt) 
          Hij(js,is,1)=dcmplx(ret,-imt) ! symmetrize hopping
       enddo
    enddo
    close(unit)
    !
    ! basis vectors must be defined
    call TB_set_bk([1d0,0d0,0d0],[0d0,1d0,0d0],[0d0,0d0,1d0])
    call TB_write_hk(Hk=Hij,file="Hij_nano.data",&
         No=Nlso,&
         Nd=Norb,&
         Np=0,&
         Nineq=Nineq,&
         Nkvec=[1,1,1])
    !
    allocate(nanoHloc(Nlso,Nlso))
    nanoHloc = extract_Hloc(Hij,Nlat,Nspin,Norb)
    !
    !save lat2ineq,ineq2lat arrays
    unit=free_unit()
    open(unit,file="lat2ineq.ed")
    do ilat=1,Nlat
       write(unit,*)ilat,lat2ineq(ilat)
    enddo
    close(unit)
    unit=free_unit()
    open(unit,file="ineq2lat.ed")
    do i=1,Nineq
       write(unit,*)i,ineq2lat(i)
    enddo
    close(unit)
  end subroutine build_Hij



  !----------------------------------------------------------------------------------------!
  ! purpose: read the matsubara local self-energy from disk
  !----------------------------------------------------------------------------------------!
  subroutine read_sigma_mats(Smats)
    complex(8),intent(inout)         :: Smats(:,:,:,:,:,:)
    character(len=30)                :: suffix
    integer                          :: ilat,ispin,iorb
    real(8),dimension(:),allocatable :: wm,wr

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
          call sread("LSigma"//trim(suffix),wm,Smats(:,ispin,ispin,iorb,iorb,:))
       enddo
    enddo

  end subroutine read_sigma_mats


  !----------------------------------------------------------------------------------------!
  ! purpose: read the real local self-energy from disk
  !----------------------------------------------------------------------------------------!
  subroutine read_sigma_real(Sreal)
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
    write(LOGfile,*)"write spin-orbital diagonal elements:"
    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          call sread("LSigma"//trim(suffix),wr,Sreal(:,ispin,ispin,iorb,iorb,:))
       enddo
    enddo

  end subroutine read_sigma_real


  !----------------------------------------------------------------------------------------!
  ! purpose: read the local self-energy from disk
  !----------------------------------------------------------------------------------------!
  subroutine read_sigma(Smats,Sreal)
    complex(8),intent(inout)         :: Smats(:,:,:,:,:,:)
    complex(8),intent(inout)         :: Sreal(:,:,:,:,:,:)
    character(len=30)                :: suffix
    integer                          :: ilat,ispin,iorb
    real(8),dimension(:),allocatable :: wm,wr

    if(size(Smats,2)/=Nspin) stop "save_sigma: error in dim 2. Nspin"
    if(size(Smats,3)/=Nspin) stop "save_sigma: error in dim 3. Nspin"
    if(size(Smats,4)/=Norb) stop "save_sigma: error in dim 4. Norb"
    if(size(Smats,5)/=Norb) stop "save_sigma: error in dim 5. Norb"

    if(size(Sreal,2)/=Nspin) stop "save_sigma: error in dim 2. Nspin"
    if(size(Sreal,3)/=Nspin) stop "save_sigma: error in dim 3. Nspin"
    if(size(Sreal,4)/=Norb) stop "save_sigma: error in dim 4. Norb"
    if(size(Sreal,5)/=Norb) stop "save_sigma: error in dim 5. Norb"

    allocate(wm(Lmats))
    allocate(wr(Lreal))

    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)
    write(LOGfile,*)"write spin-orbital diagonal elements:"
    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
          call sread("LSigma"//trim(suffix),wm,Smats(:,ispin,ispin,iorb,iorb,:))
       enddo
    enddo
    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          call sread("LSigma"//trim(suffix),wr,Sreal(:,ispin,ispin,iorb,iorb,:))
       enddo
    enddo

  end subroutine read_sigma



  !----------------------------------------------------------------------------------------!
  ! purpose: evaluate 
  !  - conductance (without vertex corrections) 
  !  - bias-driven DC & AC current
  ! for a nanostructure on the real axis, given the non-local Green's function 
  ! and the L/R hybridization matrix
  !----------------------------------------------------------------------------------------!
  subroutine ed_transport_acdc(Gret,time_step)
    complex(8),intent(inout)              :: Gret(:,:,:,:,:,:,:)  ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    real(8),intent(in)                    :: time_step
    ! auxiliary variables for matmul        
    complex(8),dimension(:,:),allocatable :: GR,HR,GA,HL,Re,Le,Te ![Nlo][Nlo]
    complex(8),dimension(:,:),allocatable :: Vt                   ![Nlo][Nlo]
    !
    complex(8),dimension(:,:),allocatable :: transe_dc            ![Nspin][Lreal]
    complex(8),dimension(:,:),allocatable :: transl_ac,transr_ac  ![Nspin][Lreal]
    !
    real(8),dimension(:),allocatable      :: g                    ![Nspin]
    !
    integer,dimension(:,:),allocatable    :: rmask,lmask          ![Nlat][Nlat]
    !
    real(8),dimension(:),allocatable      :: wr
    real(8)                               :: wmesh
    !
    real(8),dimension(:),allocatable      :: jcurr_dc,jcurr_ac    ![Nspin]
    real(8)                               :: lbias,rbias
    !
    integer                               :: Nlso,Nlo
    integer                               :: ilat,jlat,ispin,jspin,iorb,jorb
    integer                               :: io,jo,is,js,i,ixmu
    integer                               :: unit,unit_in,unit_out,eof,lfile
    logical                               :: exist
    character(len=30)                     :: suffix
    !
    Nlso = Nlat*Nspin*Norb
    Nlo  = Nlat*Norb
    !
    allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal,mesh=wmesh)
    !
    !find index corresponding to chemical potential
    do i=1,Lreal
       !write(*,*) xmu, wr(i), i
       !if((wr(i)-xmu)<wmesh/2.d0) ixmu=i
       if((wr(i)-0.d0)<wmesh/2.d0) ixmu=i
    enddo
    !write(*,*) 'index selected: ',ixmu," corresponding to: ",wr(ixmu)
    !write(*,*) 'also see wr(ixmu-1: ',wr(ixmu-1)," and wr(ixmu+1): ",wr(ixmu+1)
       

    ! allocate variables for matrix-matrix multiplication
    allocate(GR(Nlo,Nlo));GR=zero
    allocate(HR(Nlo,Nlo));HR=zero
    allocate(GA(Nlo,Nlo));GA=zero
    allocate(HL(Nlo,Nlo));HL=zero
    allocate(Re(Nlo,Nlo));Re=zero
    allocate(Le(Nlo,Nlo));Le=zero
    allocate(Te(Nlo,Nlo));Te=zero
    allocate(Vt(Nlo,Nlo));Vt=zero

    ! set masks
    allocate(lmask(Nlat,Nlat),rmask(Nlat,Nlat))
    lmask(:,:)=0
    rmask(:,:)=0
    lfile = file_length("lmask.in")
    unit = free_unit()
    open(unit,file='lmask.in',status='old')
    do i=1,lfile
       read(unit,*) ilat, jlat
       ilat=ilat+1
       jlat=jlat+1
       lmask(ilat,jlat)=1
       write(6,*) ilat,jlat,lmask(ilat,jlat)
    enddo
    close(unit)
    lfile = file_length("rmask.in")
    unit = free_unit()
    open(unit,file='rmask.in',status='old')
    do i=1,lfile
       read(unit,*) ilat, jlat
       ilat=ilat+1
       jlat=jlat+1
       rmask(ilat,jlat)=1
       write(6,*) ilat,jlat,rmask(ilat,jlat)
    enddo
    close(unit)

    ! allocate spin-resolved transmission coefficient
    allocate(transe_dc(Nspin,Lreal),transl_ac(Nspin,Lreal),transr_ac(Nspin,Lreal))

    do ispin=1,Nspin
       do i=1,Lreal
          ! fill auxiliary matrix [Nlso]**2
          do ilat=1,Nlat
             do jlat=1,Nlat
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb +  (ilat-1)*Norb
                      jo = jorb +  (jlat-1)*Norb
                      is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb !== ilat
                      js = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb !== jlat
                      !
                      ! retarded Green's function
                      GR(io,jo)=Gret(ilat,jlat,ispin,ispin,iorb,jorb,i)
                      !
                      ! set \Gamma matrix for L/R according to masks to select L-subset OR R-subset
                      ! R-subset
                      HR(io,jo)=zero
                      if(rmask(ilat,jlat)==1) HR(io,jo) = cmplx(2.d0*dimag(Hyb_real(is,js,i)),0d0)
                      ! L-subset
                      HL(io,jo)=zero
                      if(lmask(ilat,jlat)==1) HL(io,jo) = cmplx(2.d0*dimag(Hyb_real(is,js,i)),0d0)
                      !
                      ! time derivative of the drive, spin-resolved
                      Vt(io,jo)=dconjg(Uij(is,js))
                      !
                   enddo
                enddo
             enddo
          enddo
          ! advanced Green's function
          GA=conjg(transpose(GR))
          !
          ! get DC transmission function as T(ispin,i)=Tr[Gadvc*Hybl*Gret*Hybr]
          Re = matmul(GR,HR)
          Le = matmul(GA,HL)
          Te = matmul(Le,Re)
          transe_dc(ispin,i) = trace_matrix(Te,Nlo)
          !
          ! get left AC transmission function as T(ispin,i)=Tr[(dV/dt)*Gadv*Hybl*Gret]
          Re = matmul(GA,matmul(HL,GR))
          Te = matmul(Vt,Re)
          transl_ac(ispin,i) = trace_matrix(Te,Nlo)
          !
          ! get right AC transmission function as T(ispin,i)=Tr[(dV/dt)*Gadv*Hybr*Gret]
          Re = matmul(GA,matmul(HR,GR))
          Te = matmul(Vt,Re)
          transr_ac(ispin,i) = trace_matrix(Te,Nlo)
 
       enddo
    enddo
    !
    if(ed_verbose==3)then
       ! write transport coefficients of disk
       do ispin=1,Nspin
          !suffix="_s"//reg(txtfy(ispin))//"_realw.ed"
          suffix="_t"//reg(txtfy(itime))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          call splot("Te_dc"//trim(suffix),wr,transe_dc(ispin,:))
       enddo
       do ispin=1,Nspin
          suffix="_s"//reg(txtfy(ispin))//"_realw.ed"
          call splot("Tl_ac"//trim(suffix),wr,transl_ac(ispin,:))
       enddo
       do ispin=1,Nspin
          suffix="_s"//reg(txtfy(ispin))//"_realw.ed"
          call splot("Tr_ac"//trim(suffix),wr,transr_ac(ispin,:))
       enddo
    endif
    
    deallocate(GR,HR,GA,HL)
    deallocate(rmask,lmask)
    deallocate(Re,Le)
    deallocate(Te)



    ! evaluate conductance from transmission coefficient
    ! g(ispin)=T(ispin,0)
    !
    ! allocate spin-resolved conductance
    allocate(g(Nspin));g(:)=0.d0
    !
    unit_out= free_unit()
    inquire(file="gcond_dc.ed",exist=exist)
    if (exist) then
      open(unit_out,file="gcond_dc.ed",status="old",position="append")
    else
      open(unit_out,file="gcond_dc.ed")
    endif
    ! 
    ! write effective time [0:pi2]
    write(unit_out,'(3f16.9)',advance='no')time_step
    do ispin=1,Nspin
       ! extract conductance at the chemical potential
       g(ispin) = real(transe_dc(ispin,ixmu))
       ! write spin-resolved conductance on disk
       write(unit_out,'(1f16.9)',advance='no')g(ispin)
    enddo
    close(unit_out)
    !
    deallocate(g)


    ! evaluate spin-resolved DC current as:
    ! J = \int_{-\infty}^{\infty} de T_dc(e) (f_L(e)-f_R(e))
    !
    allocate(jcurr_dc(Nspin));jcurr_dc=0.d0
    !
    unit_in = free_unit()
    open(unit_in,file='jbias.in',status='old')
    !
    unit_out= free_unit()
    inquire(file="jbias_dc.ed",exist=exist)
    if(exist)then
      open(unit_out,file="jbias_dc.ed",status="old",position="append")
    else
      open(unit_out,file="jbias_dc.ed")
    endif
    do
       read(unit_in,*,IOSTAT=EOF)lbias,rbias
       if(EOF<0)exit
       ! write L/R bias voltages
       write(unit_out,'(3f16.9)',advance='no')time_step,lbias,rbias
       jcurr_dc=0.d0
       do ispin=1,Nspin
           do i=1,Lreal
              ! compute current by integration
              jcurr_dc(ispin) = jcurr_dc(ispin) + transe_dc(ispin,ixmu)* &
                                (fermi(wr(i)-xmu-lbias,beta)-fermi(wr(i)-xmu-rbias,beta))* &
                                abs(wfin-wini)/Lreal
           enddo
           ! write spin-resolved current on disk
           write(unit_out,'(1e16.9)',advance='no')jcurr_dc(ispin)
       enddo
       write(unit_out,*) ! newline
    enddo
    close(unit_in)
    close(unit_out)
    !
    deallocate(jcurr_dc)
    !
    deallocate(transe_dc)



    ! evaluate spin-resolved AC current (at T=0) as:
    ! J = Re{T_ac(xmu)}/(2\pi)
    !
    allocate(jcurr_ac(Nspin));jcurr_ac=0.d0
    !
    unit_out= free_unit()
    inquire(file="jbias_ac.ed",exist=exist)
    if(exist)then
      open(unit_out,file="jbias_ac.ed",status="old",position="append")
    else
      open(unit_out,file="jbias_ac.ed")
    endif
    write(unit_out,'(1f16.9)',advance='no')time_step
    jcurr_ac=0.d0
    do ispin=1,Nspin
        ! exctract left current at the chemical potential 
        jcurr_ac(ispin) = -real(transl_ac(ispin,ixmu))/pi2
        ! write spin-resolved current on disk
        write(unit_out,'(1e16.91x)',advance='no')jcurr_ac(ispin)
    enddo
    do ispin=1,Nspin
        ! exctract right current at the chemical potential 
        jcurr_ac(ispin) = -real(transr_ac(ispin,ixmu))/pi2
        ! write spin-resolved current on disk
        write(unit_out,'(1e16.91x)',advance='no')jcurr_ac(ispin)
    enddo
    write(unit_out,*) ! newline
    close(unit_out)
    !
    deallocate(jcurr_ac)
    !
    deallocate(transl_ac,transr_ac)


  end subroutine ed_transport_acdc



  !----------------------------------------------------------------------------------------!
  ! purpose: define the hybridization matrix of size [Nlat][Nlat][Nspin][Norb][Norb][Lreal] 
  ! reading the parameters from an input file
  !----------------------------------------------------------------------------------------!
  subroutine set_hyb()
    integer                                 :: ilat,jlat,ispin,jspin,iorb,jorb,io,jo,i,Nlso
    integer                                 :: k,kmax
    integer                                 :: unit,l,lfile
    ! leads
    integer                                 :: ikind,ilead,Nlead
    real(8)                                 :: D,mu,V,epsk
    complex(8)                              :: ksum
    complex(8),dimension(:,:,:),allocatable :: lead_real,lead_mats ![Nlead][Nspin][Lreal/Lmats]
    real(8),dimension(:),allocatable        :: wr,wm
    character(50)                           :: suffix
    !
    Nlso = Nlat*Nspin*Norb
    !
    kmax=10000
    !
    allocate(wm(Lmats),wr(Lreal))
    wm = pi/beta*(2*arange(1,Lmats)-1)
    wr = linspace(wini,wfin,Lreal)

    ! initialize embedding hybridization function
    allocate(Hyb_mats(Nlso,Nlso,Lmats))
    allocate(Hyb_real(Nlso,Nlso,Lreal))
    Hyb_mats=zero
    Hyb_real=zero

    ! determine Nleads & allocate lead matrix
    lfile = file_length("lead.in")
    unit = free_unit()
    open(unit,file='lead.in',status='old')
    read(unit,*)Nlead
    allocate(lead_real(Nlead,Nspin,Lreal))
    allocate(lead_mats(Nlead,Nspin,Lmats))
    lead_real(:,:,:)=zero
    ! lead file setup lead by kind, half-bandwitdh (D) and chemical potential (mu)
    do l=1,lfile-1 ! because Nlead was read separately above
       read(unit,*) ilead, ispin, D, mu, ikind
       ilead=ilead+1
       ispin=ispin+1
       if(ilead>Nlead)stop "set_hyb error: in input file 'lead.in' ilead > Nlead"
       if(ispin>Nspin)stop "set_hyb error: non-spin degenerate leads for Nspin=1 calculation"
       suffix="_ilead"//reg(txtfy(ilead))//"_s"//reg(txtfy(ispin))
       !
       ! set the lead's Green's function, depending on ikind
       if(ikind==0)then
          ! flat DOS (analytic)
          write(*,*) "flat DOS (analytic)"
          lead_real(ilead,ispin,:)=dcmplx( log(abs((D+wr(:)+mu)/(D-wr(:)-mu))) , -pi*heaviside(D-abs(wr(:)+mu)) )/(2d0*D)
       elseif(ikind==1)then
          ! flat DOS (k-sum)
          write(*,*) "flat DOS (k-sum)"
          do i=1,Lreal
             ksum=zero
             do k=1,kmax
                epsk = -D + 2*D/kmax*(k-1)
                ksum = ksum + 1d0/( wr(i)+xi*eps+mu - epsk)
             enddo
             lead_real(ilead,ispin,i)=ksum/kmax
          enddo
       elseif(ikind==2)then
          ! broad-band limit
          write(*,*) "broad-band limit (analytic)" 
          lead_real(ilead,ispin,:)=dcmplx(0d0,-1.d0*pi) ! to ensure DOS normalization
       elseif(ikind==3)then
          ! semicircular DOS (k-sum) 
          write(*,*) "semicircular DOS (k-sum)"
          do i=1,Lreal
             ksum=zero
             do k=1,kmax
                epsk = -D + 2*D/kmax*(k-1)
                ksum = ksum + (4d0/(pi*kmax))*sqrt(1d0-(epsk/D)**2)/( wr(i)+xi*eps+mu - epsk)
             enddo
             lead_real(ilead,ispin,i)=ksum
          enddo
       elseif(ikind==4)then
          ! readin hk DOS
          write(*,*) "readin hk DOS to be implemented and benchmarked w/ w2dynamics"
          stop
       else
          write(*,*) "set_hyb error: in input file 'lead.in' invalid ikind"
          stop
       endif
       ! store lead(s) DOS on disk
       suffix="_ilead"//reg(txtfy(ilead))//"_s"//reg(txtfy(ispin))//"_realw.ed"
       call splot("lead"//trim(suffix),wr,lead_real(ilead,ispin,:))
       call get_matsubara_gf_from_dos(wr,lead_real(ilead,ispin,:),lead_mats(ilead,ispin,:),beta)
       suffix="_ilead"//reg(txtfy(ilead))//"_s"//reg(txtfy(ispin))//"_iw.ed"
       call splot("lead"//trim(suffix),wm,lead_mats(ilead,ispin,:))
    enddo
    close(unit)
    !
    ! hybridization file determine lead-site connections 
    lfile = file_length("vij.in")
    unit = free_unit()
    open(unit,file='vij.in',status='old')
    do i=1,lfile
       read(unit,*) ilat, iorb, jlat, jorb, ilead, V
       ilat=ilat+1
       iorb=iorb+1
       jlat=jlat+1
       jorb=jorb+1
       ilead=ilead+1
       if((iorb>Norb).or.(jorb>Norb))stop "set_hyb error: in input file 'vij.in' i/jorb > Norb"
       if((ilat>Nlat).or.(jlat>Nlat))stop "set_hyb error: in input file 'vij.in' i/jlat > Nlat"
       if(ilead>Nlead)stop "set_hyb error: in input file 'vij.in' ilead > Nlead"
       do ispin=1,Nspin
          ! get stride and set matrix element: no symmetrization
          io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb !== ilat
          jo = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb !== jlat
          Hyb_real(io,jo,:)=Hyb_real(io,jo,:)+lead_real(ilead,ispin,:)*V**2
          Hyb_mats(io,jo,:)=Hyb_mats(io,jo,:)+lead_mats(ilead,ispin,:)*V**2
          suffix="_i"//reg(txtfy(ilat))//"_j"//reg(txtfy(jlat))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          call splot("Hyb"//trim(suffix),wr,Hyb_real(io,jo,:))
       enddo
    enddo
    close(unit)
    deallocate(lead_real,lead_mats,wr,wm)
    
  end subroutine set_hyb



  !----------------------------------------------------------------------------------------!
  ! purpose: define the adiabatic non-equilibrium drive 
  !----------------------------------------------------------------------------------------!
  subroutine set_drive(drive_model,Vijt,Uijt,Hij0,Hijt,time_step)
    complex(8),dimension(:,:),intent(inout)   :: Vijt ![Nlso][Nlso]
    complex(8),dimension(:,:),intent(inout)   :: Uijt ![Nlso][Nlso]
    complex(8),dimension(:,:,:),intent(inout) :: Hij0 ![Nlso][Nlso][Nk]
    complex(8),dimension(:,:,:),intent(inout) :: Hijt ![Nlso][Nlso][Nk]
    real(8),intent(in)                        :: time_step
    integer                                   :: i,j
    integer                                   :: unit
    !
    ! set size of the drive from time-independent Hamiltonian
    Nlso=size(Hij0,1)
    !
    ! reset drive
    Vijt(:,:) = zero
    Uijt(:,:) = zero
    !
    ! set drive and time-dependent Hamiltonian
    call drive_model(Vijt,Uijt,Hij0(:,:,1),Hijt(:,:,1),time_step,Nlso)
    !
    ! write drive on HDD
    if(ed_verbose==3)then
       unit = free_unit()
       inquire(file="drive.out",exist=exist)
       if(exist)then
         open(unit,file="drive.out",status="old",position="append")
       else
         open(unit,file="drive.out")
       endif
       do i=1,Nlso
          do j=1,Nlso
             write(unit,'(1f16.9,2i6,8f16.9)')time_step,i,j,&
                dreal(Vijt(i,j)),dimag(Vijt(i,j)),dreal(Uijt(i,j)),dimag(Uijt(i,j)),&
                dreal(Hij0(i,j,1)),dimag(Hij0(i,j,1)),dreal(Hijt(i,j,1)),dimag(Hijt(i,j,1))
          enddo
       enddo
       close(unit)
    endif
    !
  end subroutine set_drive





  !----------------------------------------------------------------------------------------!
  ! drive: time-dependent angle between phenyl rings
  !----------------------------------------------------------------------------------------!
  subroutine angle(Vijt,Uijt,Hij0,Hijt,time_step,N)
    integer                   :: N
    real(8)                   :: time_step
    complex(8),dimension(N,N) :: Vijt !drive
    complex(8),dimension(N,N) :: Uijt !drive time derivative
    complex(8),dimension(N,N) :: Hij0 !time-independent Hamiltonian
    complex(8),dimension(N,N) :: Hijt !time-  dependent Hamiltonian
    !
    integer                   :: i,ilat,jlat,iorb,jorb,is,js
    integer                   :: lfile
    real(8)                   :: Vg,delta,omega

    ! ----------------------------------------------------------
    ! explicit form of time-dependent Hamiltonian & drive:
    !
    ! H_{ij}(t) = t_{ph-ph} = t cos (\omega t + delta)
    ! V_{ij}(t) = t cos(\omega t + \delta)
    ! U_{ij}(t) = - t \omega sin(\omega t + \delta)
    !
    ! time_step: effective variable time (i.e., includes \omega)
    ! t_{ph-ph} : hopping between phenyl ring
    ! t         : nearest neighbor C-C hopping 
    ! delta     : phase of the drive (in units of pi)
    ! omega     : frequency of the drive
    ! ----------------------------------------------------------

    ! set time-dependent Hamiltonian to the static one
    Hijt(:,:)=Hij0(:,:)

    ! readin drive parameters from input (see format in the read below)
    lfile = file_length("drive.in")
    unit = free_unit()
    open(unit,file='drive.in',status='old')
    do i=1,lfile
       read(unit,*) ilat, iorb, jlat, jorb, ispin, delta, omega
       ilat=ilat+1
       iorb=iorb+1
       jlat=jlat+1
       jorb=jorb+1
       ispin=ispin+1
       !
       if((ilat>Nlat).or.(jlat>Nlat))stop "set_drive error: in input file 'drive.in' i/jlat > Nlat"
       if((iorb>Norb).or.(jorb>Norb))stop "set_drive error: in input file 'drive.in' i/jorb > Norb"
       if((ispin>Nspin))stop "set_drive error: in input file 'drive.in' ispin > Nspin"
       !
       ! get stride and set matrix element
       is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
       js = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb
       !
       ! set drive: 
       Vijt(is,js)=cos(delta*pi+time_step)
       Vijt(js,is)=dconjg(Vijt(is,js)) ! hermitian drive
       !
       ! set drive time derivative:
       Uijt(is,js)=-omega*sin(delta*pi+time_step)
       Uijt(js,is)=dconjg(Uijt(is,js)) ! hermitian drive
       !
       ! set time-dependent Hamiltonian
       Hijt(is,js)=Hij0(is,js)*Vijt(is,js)
       Hijt(js,is)=Hij0(js,is)*Vijt(js,is) ! hermitian Hamiltonian
       !
    enddo
    close(unit)
    !
  end subroutine angle




  !----------------------------------------------------------------------------------------!
  ! drive: sinusoidal flux 
  !----------------------------------------------------------------------------------------!
  subroutine flux(Vijt,Uijt,Hij0,Hijt,time_step,N)
    integer                   :: N
    real(8)                   :: time_step
    complex(8),dimension(N,N) :: Vijt !drive
    complex(8),dimension(N,N) :: Uijt !drive time derivative
    complex(8),dimension(N,N) :: Hij0 !time-independent Hamiltonian
    complex(8),dimension(N,N) :: Hijt !time-  dependent Hamiltonian
    !
    integer                   :: i,ilat,jlat,iorb,jorb,is,js
    integer                   :: lfile
    real(8)                   :: phidc,phiac,delta,omega

    ! ----------------------------------------------------------
    ! explicit form of time-dependent Hamiltonian & drive:
    !
    ! H_{ij}(t) = t V_{ij}(t)
    ! V_{ij}(t) = \exp(\imath (\phi_{DC} + \phi_{AC} cos(\omega t + \delta) ))
    ! U_{ij}(t) = -\imath \phi_{AC} \omega sin(\omega t + \delta) * V_{ij}(t)
    !
    ! time_step: effective variable time (i.e., includes \omega)
    ! phidc    : DC part of the drive
    ! phiac    : AC part of the drive
    ! delta    : phase of the drive (in units of pi)
    ! omega    : frequency of the drive
    ! ----------------------------------------------------------

    ! set time-dependent Hamiltonian to the static one
    Hijt(:,:)=Hij0(:,:)

    ! readin drive parameters from input (see format in the read below)
    lfile = file_length("drive.in")
    unit = free_unit()
    open(unit,file='drive.in',status='old')
    do i=1,lfile
       read(unit,*) ilat, iorb, jlat, jorb, ispin, phidc, phiac, delta, omega
       ilat=ilat+1
       iorb=iorb+1
       jlat=jlat+1
       jorb=jorb+1
       ispin=ispin+1
       !
       if((ilat>Nlat).or.(jlat>Nlat))stop "set_drive error: in input file 'drive.in' i/jlat > Nlat"
       if((iorb>Norb).or.(jorb>Norb))stop "set_drive error: in input file 'drive.in' i/jorb > Norb"
       if((ispin>Nspin))stop "set_drive error: in input file 'drive.in' ispin > Nspin"
       !
       ! get stride and set matrix element
       is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
       js = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb
       !
       Vijt(is,js)=exp(xi*(phidc+phiac*cos(delta*pi+time_step)))
       Vijt(js,is)=dconjg(Vijt(is,js)) ! hermitian drive
       !
       ! set drive time derivative:
       Uijt(is,js)=-xi*phiac*omega*sin(delta*pi+time_step)*exp(xi*(phidc+phiac*cos(delta*pi+time_step)))
       Uijt(js,is)=dconjg(Uijt(is,js)) ! hermitian drive
       !
       ! set time-dependent Hamiltonian
       Hijt(is,js)=Hij0(is,js)*Vijt(is,js)
       Hijt(js,is)=Hij0(js,is)*Vijt(js,is) ! hermitian Hamiltonian
       !
    enddo
    close(unit)
    !
  end subroutine flux





  !----------------------------------------------------------------------------------------!
  ! drive: both ph-ph angle and sinusoidal flux 
  !----------------------------------------------------------------------------------------!
  subroutine pump(Vijt,Uijt,Hij0,Hijt,time_step,N)
    integer                   :: N
    real(8)                   :: time_step
    complex(8),dimension(N,N) :: Vijt !drive
    complex(8),dimension(N,N) :: Uijt !drive time derivative
    complex(8),dimension(N,N) :: Hij0 !time-independent Hamiltonian
    complex(8),dimension(N,N) :: Hijt !time-  dependent Hamiltonian
    !
    integer                   :: i,ilat,jlat,iorb,jorb,is,js
    integer                   :: ientry
    real(8)                   :: phidc,phiac,psi,omega,angle,delta


    ! ----------------------------------------------------------
    ! explicit form of time-dependent Hamiltonian & drive:
    !
    ! angle:
    ! ----------------------------------------------------------
    ! H_{ij}(t) = t_{ph-ph} = t cos (\omega t + delta)
    ! V_{ij}(t) = t cos(\omega t + \delta)
    ! U_{ij}(t) = - t \omega sin(\omega t + \delta)
    !
    ! flux:
    ! ----------------------------------------------------------
    ! H_{ij}(t) = t V_{ij}(t)
    ! V_{ij}(t) = \exp(\imath (\phi_{DC} + \phi_{AC} cos(\omega t + \delta) ))
    ! U_{ij}(t) = -\imath \phi_{AC} \omega sin(\omega t + \delta) * V_{ij}(t)
    !
    ! time_step: effective variable time (i.e., includes \omega)
    ! t_{ph-ph} : hopping between phenyl ring
    ! t         : nearest neighbor C-C hopping 
    ! delta     : phase of the angle (in units of pi)
    ! phidc     : DC part of the drive
    ! phiac     : AC part of the drive
    ! pis       : phase of the flux (in units of pi)
    ! omega     : frequency of the drive
    ! ----------------------------------------------------------


    ! set time-dependent Hamiltonian to the static one
    Hijt(:,:)=Hij0(:,:)

    ! readin drive parameters from input (see format in the read below)
    unit = free_unit()
    open(unit,file='drive.in',status='old')

    ! -------------------------------------------------------------
    ! *** angle
    ! -------------------------------------------------------------
    ! read # entries for this drive
    read(unit,*) ientry
    do i=1,ientry
       ! read angle
       read(unit,*) ilat, iorb, jlat, jorb, ispin, delta, omega
       ilat=ilat+1
       iorb=iorb+1
       jlat=jlat+1
       jorb=jorb+1
       ispin=ispin+1
       !
       if((ilat>Nlat).or.(jlat>Nlat))stop "set_drive error: in input file 'drive.in' i/jlat > Nlat"
       if((iorb>Norb).or.(jorb>Norb))stop "set_drive error: in input file 'drive.in' i/jorb > Norb"
       if((ispin>Nspin))stop "set_drive error: in input file 'drive.in' ispin > Nspin"
       !
       ! get stride and set matrix element
       is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
       js = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb
       !
       ! set drive: 
       Vijt(is,js)=cos(delta*pi+time_step)
       Vijt(js,is)=dconjg(Vijt(is,js)) ! hermitian drive
       !
       ! set drive time derivative:
       Uijt(is,js)=-omega*sin(delta*pi+time_step)
       Uijt(js,is)=dconjg(Uijt(is,js)) ! hermitian drive
       !
       ! set time-dependent Hamiltonian
       Hijt(is,js)=Hij0(is,js)*Vijt(is,js)
       Hijt(js,is)=Hij0(js,is)*Vijt(js,is) ! hermitian Hamiltonian
    enddo
 
    ! -------------------------------------------------------------
    ! *** flux left ring
    ! -------------------------------------------------------------
    read(unit,*) ientry
    do i=1,ientry
       ! read flux
       read(unit,*) ilat, iorb, jlat, jorb, ispin, phidc, phiac, psi, omega
       ilat=ilat+1
       iorb=iorb+1
       jlat=jlat+1
       jorb=jorb+1
       ispin=ispin+1
       !
       if((ilat>Nlat).or.(jlat>Nlat))stop "set_drive error: in input file 'drive.in' i/jlat > Nlat"
       if((iorb>Norb).or.(jorb>Norb))stop "set_drive error: in input file 'drive.in' i/jorb > Norb"
       if((ispin>Nspin))stop "set_drive error: in input file 'drive.in' ispin > Nspin"
       !
       ! get stride and set matrix element
       is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
       js = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb
       !
       Vijt(is,js)=exp(xi*(phidc+phiac*cos(psi*pi+time_step)))
       Vijt(js,is)=dconjg(Vijt(is,js)) ! hermitian drive
       !
       ! set drive time derivative:
       Uijt(is,js)=-xi*phiac*omega*sin(psi*pi+time_step)*exp(xi*(phidc+phiac*cos(psi*pi+time_step)))
       Uijt(js,is)=dconjg(Uijt(is,js)) ! hermitian drive
       !
       ! set time-dependent Hamiltonian
       Hijt(is,js)=Hij0(is,js)*Vijt(is,js)
       Hijt(js,is)=Hij0(js,is)*Vijt(js,is) ! hermitian Hamiltonian
    enddo

    ! -------------------------------------------------------------
    ! *** flux right ring (rotating)
    ! -------------------------------------------------------------
    read(unit,*) ientry
    do i=1,ientry
       ! read flux
       read(unit,*) ilat, iorb, jlat, jorb, ispin, phidc, phiac, psi, omega
       ilat=ilat+1
       iorb=iorb+1
       jlat=jlat+1
       jorb=jorb+1
       ispin=ispin+1
       !
       if((ilat>Nlat).or.(jlat>Nlat))stop "set_drive error: in input file 'drive.in' i/jlat > Nlat"
       if((iorb>Norb).or.(jorb>Norb))stop "set_drive error: in input file 'drive.in' i/jorb > Norb"
       if((ispin>Nspin))stop "set_drive error: in input file 'drive.in' ispin > Nspin"
       !
       ! get stride and set matrix element
       is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
       js = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb
       !
       Vijt(is,js)=exp(xi*(phidc*cos(delta*pi+time_step)+phiac*cos(delta*pi+time_step)*cos(psi*pi+time_step)))
       Vijt(js,is)=dconjg(Vijt(is,js)) ! hermitian drive
       !
       ! set drive time derivative:
       Uijt(is,js)=-xi*phidc*omega*sin(delta*pi+time_step)&
                   -xi*phiac*omega&
                   *(sin(delta*pi+time_step)*cos(psi*pi+time_step)+cos(delta*pi+time_step)*sin(psi*pi+time_step))&
                   *exp(xi*(phidc*cos(delta*pi+time_step)+phiac*cos(delta*pi+time_step)*cos(psi*pi+time_step)))
       Uijt(js,is)=dconjg(Uijt(is,js)) ! hermitian drive
       !
       ! set time-dependent Hamiltonian
       Hijt(is,js)=Hij0(is,js)*Vijt(is,js)
       Hijt(js,is)=Hij0(js,is)*Vijt(js,is) ! hermitian Hamiltonian
    enddo

    close(unit)
    !
  end subroutine pump





  function trace_matrix(M,dim) result(tr)
    integer                       :: dim
    complex(8),dimension(dim,dim) :: M
    complex(8) :: tr
    integer                       :: i
    tr=dcmplx(0d0,0d0)
    do i=1,dim
       tr=tr+M(i,i)
    enddo
  end function trace_matrix


  function extract_Hloc(Hk,Nlat,Nspin,Norb) result(Hloc)
    complex(8),dimension(:,:,:)                 :: Hk
    integer                                     :: Nlat,Nspin,Norb
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Hloc
    !
    integer                                     :: iorb,ispin,ilat,is
    integer                                     :: jorb,jspin,js
    Hloc = zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hloc(is,js) = sum(Hk(is,js,:))/size(Hk,3)
                enddo
             enddo
          enddo
       enddo
    enddo
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
  end function extract_Hloc




end program ed_nano_adiabatic
