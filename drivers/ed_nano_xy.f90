program ed_nano_isoc
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none

  integer                                       :: iloop!,Nx
  logical                                       :: converged
  integer                                       :: ilat,ineq,ispin,iorb,jspin
  !bath:
  integer                                       :: Nb
  real(8),allocatable                           :: Bath_prev(:,:),Bath_ineq(:,:)
  !local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal,Sreal_ineq ![Nlat*(Nspin*Norb)**2*Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal,Greal_ineq
  real(8), allocatable,dimension(:)             :: dens,dens_ineq
  real(8), allocatable,dimension(:)             :: Sz,Sz_ineq
  real(8), allocatable,dimension(:)             :: docc,docc_ineq
  !hamiltonian input:
  complex(8),allocatable                        :: Hij(:,:,:) ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk==1]
  complex(8),allocatable                        :: nanoHloc(:,:),Hloc(:,:,:,:,:),Hloc_ineq(:,:,:,:,:)
  integer                                       :: Nk,Nlso,Nineq,Nlat,Nk_chi0,Nq_chi0
  integer,dimension(:),allocatable              :: lat2ineq,ineq2lat
  integer,dimension(:),allocatable              :: sb_field_sign
  !
  real(8)                                       :: wmixing,Eout(2),Nread_bkp
  !input files:
  character(len=32)                             :: finput
  character(len=32)                             :: nfile,hijfile,hisocfile
  !
  !hybridization function to environment
  real(8)                                       :: e1(2),e2(2),nss(2,2)
  integer                                       :: unitHIJ,jlat
  complex(8),allocatable                        :: Hij_test(:,:,:)
  integer                                       :: is,js,jorb,top_index,top_jndex
  complex(8),allocatable                        :: Hij_top(:,:),Smats_rep(:,:),xmu_mat(:,:)
  logical                                       :: Ioptimize
  real(8),dimension(2)                          :: aparams
  real(8),dimension(:,:),allocatable            :: Sig

  integer :: comm
  logical :: master


  call Init_MPI(comm)
  master = get_master_MPI()

  call parse_cmd_variable(finput,"FINPUT",default='inputED_NANO.conf')
  call parse_input_variable(nfile,"NFILE",finput,default="nano.in")
  call parse_input_variable(hijfile,"HIJFILE",finput,default="hij.in")
  call parse_input_variable(hisocfile,"HISOCFILE",finput,default="hisoc.in")
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call ed_read_input(trim(finput),comm)


  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")
  call add_ctrl_var(ULOC,"uloc")

  ! set input structure hamiltonian
  call build_Hij([nfile,hijfile,hisocfile])


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
  ! allocate Hloc
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb))
  ! allocate density
  allocate(dens(Nlat),Sz(Nlat))
  allocate(dens_ineq(Nineq),Sz_ineq(Nineq))
  ! allocate double occupations
  allocate(docc(Nlat))
  allocate(docc_ineq(Nineq))

  !Hloc = reshape_Hloc(nanoHloc,Nlat,Nspin,Norb)
  Hloc = lso2nnn_reshape(nanoHloc,Nlat,Nspin,Norb)

  allocate(Sig(Nineq,3))

  ! setup solver
  Nb=get_bath_dimension()

  allocate(Bath_ineq(Nineq,Nb))
  allocate(Bath_prev(Nineq,Nb))
  call ed_init_solver(Bath_ineq, Hloc_ineq)

  do ineq=1,Nineq
     ilat = ineq2lat(ineq)
     ! break SU(2) symmetry for magnetic solutions
     if(Nspin>1) call break_symmetry_bath(Bath_ineq(ineq,:),sb_field,dble(sb_field_sign(ineq)))
     Hloc_ineq(ineq,:,:,:,:) = Hloc(ilat,:,:,:,:)
  enddo


  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")   
     bath_prev=bath_ineq

     ! solve impurities on each inequivalent site:
     call ed_solve(comm,bath_ineq,Hloc_ineq)

     ! retrieve self-energies and occupations(Nineq,Norb=1)
     call ed_get_sigma_matsubara(Smats_ineq,Nineq)
     call ed_get_sigma_real(Sreal_ineq,Nineq)
     call ed_get_dens(dens_ineq,Nineq,iorb=1)
     call ed_get_docc(docc_ineq,Nineq,iorb=1)
     call ed_get_mag(Sz_ineq,Nineq,iorb=1)


     ! spread self-energies and occupation to all lattice sites
     do ilat=1,Nlat
        ineq = lat2ineq(ilat)
        dens(ilat) = dens_ineq(ineq)
        docc(ilat) = docc_ineq(ineq)
        Sz(ilat) = Sz_ineq(ineq)
        Smats(ilat,:,:,:,:,:) = Smats_ineq(ineq,:,:,:,:,:)
        Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
     enddo

     ! compute the local gf:
     !
     call dmft_gloc_matsubara(comm,Hij,[1d0],Gmats,Smats)
     do ineq=1,Nineq
        ilat = ineq2lat(ineq)
        Gmats_ineq(ineq,:,:,:,:,:) = Gmats(ilat,:,:,:,:,:)
     enddo
     !

     if(master)open(345,file="magXY.ed")
     do ineq=1,Nineq
        Sig(ineq,1) = 0.5d0*(fft_get_density(Gmats_Ineq(ineq,1,2,1,1,:),beta)+fft_get_density(Gmats_Ineq(ineq,2,1,1,1,:),beta))
        Sig(ineq,2) = 0.5d0*(fft_get_density(xi*Gmats_Ineq(ineq,1,2,1,1,:),beta)+fft_get_density(xi*Gmats_Ineq(ineq,2,1,1,1,:),beta))
        Sig(ineq,3) = 0.5d0*(fft_get_density(Gmats_Ineq(ineq,1,1,1,1,:),beta)-fft_get_density(Gmats_Ineq(ineq,2,2,1,1,:),beta))
        if(master)write(345,*)ineq,Sig(ineq,:)
     enddo
     if(master)close(345)



     ! compute the Weiss field
     if(cg_scheme=="weiss")then
        call dmft_weiss(comm,Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq)
     else
        call dmft_delta(comm,Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq)
     endif

     ! fit baths and mix result with old baths
     do ispin=1,Nspin
        call ed_chi2_fitgf(comm,bath_ineq,Weiss_ineq,Hloc_ineq)
     enddo

     Bath_ineq=wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev

     if(master)then
        converged = check_convergence(Weiss_ineq(1,1,1,1,1,:),dmft_error,nsuccess,nloop)
        ! alternative convergency criteria
        !converged = check_convergence_local(docc_ineq,dmft_error,nsuccess,nloop)

        !if(NREAD/=0.d0) call search_chemical_potential(xmu,sum(dens)/Nlat,converged)
        if(NREAD/=0.d0)call ed_search_variable(xmu,sum(dens)/Nlat,converged)
     endif
     call bcast_MPI(comm,converged)
     call bcast_MPI(comm,xmu)


     if(master)call end_loop()
  end do


  if(master)call dmft_print_gf_matsubara(Gmats,"Gmats",iprint=6)
  call dmft_gloc_realaxis(comm,Hij,[1d0],Greal,Sreal)
  do ineq=1,Nineq
     ilat = ineq2lat(ineq)
     Greal_ineq(ineq,:,:,:,:,:) = Greal(ilat,:,:,:,:,:)
  enddo
  if(master)call dmft_print_gf_realaxis(Greal,"Greal",iprint=6)



  ! compute kinetic energy at convergence
  ! call dmft_kinetic_energy(comm,Hij,[1d0],Smats)


  call Finalize_MPI()

contains





  !----------------------------------------------------------------------------------------!
  ! purpose: build real-space Hamiltonian for a nanostructure of size [Nlat*Nspin*Norb]**2
  !----------------------------------------------------------------------------------------!
  subroutine build_Hij(file)
    character(len=*)     :: file(3)
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
    !
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
          ! hermitian hopping
          Hij(is,js,1)=dcmplx(ret, imt)
          Hij(js,is,1)=dcmplx(ret,-imt)
       enddo
    enddo
    close(unit)
    !


    ! if Nspin!=2 raise error
    if(Nspin/=2)stop "build_Hij error: cannot implement intrinsic SOC with Nspin/=2"
    unit = free_unit()
    open(unit,file=trim(file(3)),status='old')
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
          ! hermitian (imaginary) hopping w/ spin-antisymmetric intrinsic SOC
          Hij(is,js,1)=dcmplx(ret, imt)*(1-2*(ispin-1))
          Hij(js,is,1)=dcmplx(ret,-imt)*(1-2*(ispin-1))
       enddo
    enddo
    close(unit)
    !
    !! basis vectors must be defined
    call TB_set_bk([1d0,0d0,0d0],[0d0,1d0,0d0],[0d0,0d0,1d0])
    call TB_write_hk(Hk=Hij,file="Hij_nano.data",&
         Nlat=Nlat,Nspin=Nspin,Norb=Norb,Nkvec=[1,1,1])
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


end program ed_nano_isoc
