program ed_LVO_hetero
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  !
  !#############################################
  !#                                           #
  !#       THIS CODE IS MPI COMPILED ONLY      #
  !#                                           #
  !#############################################
  !
  !
  !#########   VARIABLEs DECLARATION   #########
  !
  integer                                             :: iloop
  integer                                             :: i,j,io,jo,ndx
  integer                                             :: iorb,jorb
  integer                                             :: ispin,jspin
  integer                                             :: ilat,jlat,Nlat
  integer                                             :: ifreq,Lfreq
  integer                                             :: ilayer,Nlayer
  integer                                             :: NlNsNo,NsNo
  logical                                             :: converged
  real(8)                                             :: wmixing
  character(len=60)                                   :: finput
  character(len=32)                                   :: geometry
  character(len=32)                                   :: z_symmetry
  !Mpi:
  integer                                             :: comm,rank,ier
  logical                                             :: master
  !Bath:
  integer                                             :: Nb
  real(8)   ,allocatable,dimension(:,:)               :: Bath
  real(8)   ,allocatable,dimension(:)                 :: Bath_single
  !Hamiltoninas: 
  integer                                             :: ik,Nk,Lk,Nkpath
  complex(8),allocatable,dimension(:,:,:)             :: Hk
  complex(8),allocatable,dimension(:,:)               :: Hloc_lso
  complex(8),allocatable,dimension(:,:,:,:,:)         :: Hloc_nnn
  real(8)   ,allocatable,dimension(:)                 :: Wtk
  !local dmft fields:
  complex(8),allocatable,dimension(:,:,:,:,:,:)       :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:)       :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:,:,:,:)       :: field,field_old
  !Irreducible dmft fields:
  complex(8),allocatable,dimension(:,:,:,:,:,:)       :: Smats_hetero,Sreal_hetero
  complex(8),allocatable,dimension(:,:,:,:,:)         :: Smats_single,Sreal_single
  !meshes:
  real(8)                                             :: dw
  real(8)   ,allocatable,dimension(:)                 :: wr,wm
  !convergence test:
  integer                                             :: memory
  complex(8),allocatable,dimension(:)                 :: conv_funct
  !custom variables for chempot search:
  logical                                             :: converged_n
  integer                                             :: conv_n_loop=0
  real(8)                                             :: sumdens,xmu_old
  real(8)   ,allocatable,dimension(:,:)               :: orb_dens_lat,orb_mag_lat
  real(8)   ,allocatable,dimension(:)                 :: orb_dens_single,orb_mag_single
  logical                                             :: look4n=.true.
  !custom variables misc:
  logical                                             :: computeG0loc
  logical                                             :: lattice_flag=.true.
  logical                                             :: bulk_magsym
  complex(8),allocatable,dimension(:,:)               :: U,Udag
  complex(8),allocatable,dimension(:,:)               :: zeta
  complex(8),allocatable,dimension(:,:,:)             :: Gloc
  !Ek calculation:
  logical                                             :: computeEk
  real(8)                                             :: Ek
  real(8)   ,allocatable,dimension(:)                 :: Ekm
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)     :: Gkmats

  !#########   MPI INITIALIZATION   #########
  !
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  master = get_Master_MPI(comm)
  !
  !
  !#########    VARIABLE PARSING    #########
  !
  !
  call parse_cmd_variable(finput,           "FINPUT",              default='inputED_LVO.in')
  call parse_input_variable(nk,             "NK",finput,           default=10)
  call parse_input_variable(NLAT,           "NLAT",finput,         default=4)
  call parse_input_variable(nkpath,         "NKPATH",finput,       default=20)
  call parse_input_variable(wmixing,        "WMIXING",finput,      default=0.5d0)
  call parse_input_variable(computeG0loc,   "COMPUTEG0loc",finput, default=.false.)
  call parse_input_variable(geometry,       "GEOMETRY",finput,     default="bulk")
  call parse_input_variable(z_symmetry,     "ZSYMMETRY",finput,    default="FERRO")
  call parse_input_variable(bulk_magsym,    "BULKMAGSYM",finput,   default=.false.)
  call parse_input_variable(memory,         "MEMORY",finput,       default=3)
  call parse_input_variable(computeEk,      "COMPUTEEK",finput,    default=.false.)
  !
  call ed_read_input(trim(finput),comm)
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(Nlat,"nlat")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(Lfit,"Lfit")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  call add_ctrl_var(ed_para,"ed_para")
  !
  geometry=reg(geometry)
  z_symmetry=reg(z_symmetry)
  if (geometry=="bulk".and.ed_para)    lattice_flag=.false.
  if (geometry=="bulk".and.bulk_magsym)lattice_flag=.false.
  if (geometry=="bulk".and.Nlat/=4) stop
  Nlayer=Nlat/2
  NlNsNo=Nlat*Nspin*Norb
  NsNo=Nspin*Norb
  !
  !
  !##################       ALLOCATION       ##################
  !
  !
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));                 Smats=zero
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));                 Gmats=zero
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));                 Sreal=zero
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));                 Greal=zero
  allocate(field(Nlat,Nspin,Nspin,Norb,Norb,Lmats));                 field=zero
  !
  allocate(Hk(NlNsNo,NlNsNo,Nk*Nk*Nk));                              Hk=zero
  allocate(Hloc_lso(NlNsNo,NlNsNo));                                 Hloc_lso=zero
  allocate(Hloc_nnn(Nlat,Nspin,Nspin,Norb,Norb));                    Hloc_nnn=zero
  allocate(Wtk(Nk*Nk*Nk));                                           Wtk=1.d0/(Nk*Nk*Nk)
  allocate(U(NlNsNo,NlNsNo));                                        U=zero
  allocate(Udag(NlNsNo,NlNsNo));                                     Udag=zero
  !
  allocate(wr(Lreal));wr=0.0d0;                                      wr=linspace(wini,wfin,Lreal,mesh=dw)
  allocate(wm(Lmats));wm=0.0d0;                                      wm = pi/beta*real(2*arange(1,Lmats)-1,8)
  !
  Lfit=min(int((Uloc(1)+3.)*(beta/pi))+100,Lmats)
  if(master)write(LOGfile,'(a12,I6,2(a12,F10.3))')"Lfit:",Lfit,"iwmax:",(pi/beta)*(2*Lfit-1),"U+2D:",Uloc(1)+3.
  !
  !
  !##################        BUILD Hk        ##################
  !
  !
  call read_myhk("LVO_hr.dat","Hk.dat","Hloc","Kpoints.dat",.true.)
  !from here Hloc_lso and Hloc_nnn are in the diagonal basis
  !
  !##################          BATH          ##################
  !
  !
  if (bath_type/="replica") then
     Nb=get_bath_dimension()
  else
     Nb=get_bath_dimension(Hloc_nnn(1,:,:,:,:))
  endif
  if(master)write(LOGfile,*)"   Bath_size: ",Nb," layers: ",Nlayer
  allocate(Bath(Nlayer,Nb));    Bath=0.0d0
  allocate(Bath_single(Nb));    Bath_single=0.0d0
  !
  !
  !##################      INIT SOLVER       ##################
  !
  !
  if (lattice_flag)then
     call ed_init_solver(Comm,Bath,Hloc_nnn(1:Nlat:2,:,:,:,:))
  else
     call ed_init_solver(Comm,Bath_single,Hloc_nnn(1,:,:,:,:))
  endif
  !
  !
  !##################          DMFT          ##################
  !
  !
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !--------  solve impurity (diagonal basis)  --------
     if (lattice_flag)then
        call ed_solve(comm,Bath,Hloc_nnn(1:Nlat:2,:,:,:,:))
     else
        call ed_solve(comm,Bath_single,Hloc_nnn(1,:,:,:,:))
     endif
     !--------    get sigmas   (diagonal basis)  --------
     if (lattice_flag)then
        allocate(Smats_hetero(Nlayer,Nspin,Nspin,Norb,Norb,Lmats));Smats_hetero=zero
        allocate(Sreal_hetero(Nlayer,Nspin,Nspin,Norb,Norb,Lreal));Sreal_hetero=zero
        call ed_get_sigma_matsubara(Smats_hetero)
        call ed_get_sigma_real(Sreal_hetero)
        call sigma_symmetrization()
        deallocate(Smats_hetero,Sreal_hetero)
     else
        allocate(Smats_single(Nspin,Nspin,Norb,Norb,Lmats));Smats_single=zero
        allocate(Sreal_single(Nspin,Nspin,Norb,Norb,Lreal));Sreal_single=zero
        call ed_get_sigma_matsubara(Smats_single)
        call ed_get_sigma_real(Sreal_single)
        call sigma_symmetrization()
        deallocate(Smats_single,Sreal_single)
     endif
     !
     !------  get local Gf's  (Wannier90 basis)  ------
     call dmft_gloc_matsubara(Comm,Hk,Wtk,Gmats,Smats,mpi_split='k')
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=6)
     !
     !------    get field     (Wannier90 basis)  --------
     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats,Smats,field,Hloc_nnn)
        call dmft_print_gf_matsubara(field,"Weiss",iprint=6)
     elseif(cg_scheme=='delta')then
        call dmft_delta(Gmats,Smats,field,Hloc_nnn)
        call dmft_print_gf_matsubara(field,"Delta",iprint=6)
     endif
     !
     !------   rotate field   (diagonal basis)  --------
     call rotate_local_funct(field,U)
     !call dmft_print_gf_matsubara(field,"Weiss_rot",iprint=6)
     !
     !------    mix field     ------
     if(iloop>1)then
        field = wmixing*field + (1.d0-wmixing)*field_old
     endif
     field_old=field
     !
     !------    fit field     ------
     if (lattice_flag)then
        if(ed_para)then
           call ed_chi2_fitgf(Comm,Bath,field(1:Nlat:2,:,:,:,:,:),Hloc_nnn(1:Nlat:2,:,:,:,:),ispin=1)
           call spin_symmetrize_bath(Bath,save=.true.)
        else
           call ed_chi2_fitgf(Comm,Bath,field(1:Nlat:2,:,:,:,:,:),Hloc_nnn(1:Nlat:2,:,:,:,:))
        endif
     else
        call set_Hloc(Hloc_nnn(1,:,:,:,:))
        if(ed_para)then
           call ed_chi2_fitgf(Comm,field(1,:,:,:,:,:),Bath_single,ispin=1)
           call spin_symmetrize_bath(Bath_single,save=.true.)
        else
           call ed_chi2_fitgf(Comm,field(1,:,:,:,:,:),Bath_single)
        endif
     endif
     !
     !---- each loop operations -----
     if(master)then
        !
        !chemical potential find
        converged_n=.true.
        sumdens=0d0
        xmu_old=xmu
        if(lattice_flag)then
           allocate(orb_dens_lat(Nlayer,Norb));orb_dens_lat=0.d0
           allocate(orb_mag_lat(Nlayer,Norb));orb_mag_lat=0.d0
           call ed_get_dens(orb_dens_lat,Nlayer)
           call ed_get_mag(orb_mag_lat,Nlayer)
           do ilayer=1,Nlayer
              sumdens=sumdens+sum(orb_dens_lat(ilayer,:))/float(Nlayer)
              write(LOGfile,*)
              write(LOGfile,'(A7,I3,100F10.4)')"  Nlat: ",ilayer,orb_dens_lat(ilayer,:),orb_mag_lat(ilayer,:)
              write(LOGfile,*)
           enddo
           deallocate(orb_dens_lat,orb_mag_lat)
        else
           allocate(orb_dens_single(Norb));orb_dens_single=0.d0
           allocate(orb_mag_single(Norb));orb_mag_single=0.d0
           call ed_get_dens(orb_dens_single)
           write(LOGfile,*)
           write(LOGfile,'(A7,I3,100F10.4)')"  Nlat: ",1,orb_dens_single(:),orb_mag_single(:)
           write(LOGfile,*)
           sumdens=sum(orb_dens_single)
           deallocate(orb_dens_single,orb_mag_single)
        endif
        write(LOGfile,*)"  n avrg:",sumdens
        !
        if(nread/=0.d0.and.look4n)then
           converged_n=.false.
           if(iloop>=2)call search_chempot(xmu,sumdens,converged_n)
        endif
        if(converged_n)then
           conv_n_loop=conv_n_loop+1
        else
           conv_n_loop=0
        endif
        !
        !convergence
        write(LOGfile,*)
        write(LOGfile,*) "   ------------------- convergence --------------------"
        allocate(conv_funct(Lmats));conv_funct=zero
        if (lattice_flag)then
           do i=1,Lmats
              conv_funct(i)=sum(nnn2lso_reshape(field(1:Nlat:2,:,:,:,:,i),Nlayer,Nspin,Norb))
           enddo
        else
           do i=1,Lmats
              conv_funct(i)=sum(nn2so_reshape(field(1,:,:,:,:,i),Nspin,Norb))
           enddo
        endif
        if(converged_n)converged = check_convergence(conv_funct,dmft_error,nsuccess,nloop)
        write(LOGfile,'(a35,L3)') "sigma converged",converged
        write(LOGfile,'(a35,L3)') "dens converged",converged_n
        converged = converged .and. converged_n
        write(LOGfile,'(a35,L3)') "total converged",converged
        write(LOGfile,'(a35,I3)') "global iloop",iloop
        write(LOGfile,'(a35,I3)') "times dens is ok",conv_n_loop
        write(LOGfile,*) "   ----------------------------------------------------"
        write(LOGfile,*)
        deallocate(conv_funct)
     endif
     call Bcast_MPI(Comm,xmu)
     call Bcast_MPI(Comm,converged)
     call MPI_Barrier(Comm,ier)
     !
     !kinetic energy
     if(computeEk)then
        allocate(Ekm(Lmats));                                        Ekm=0d0
        allocate(Gkmats(Nk*Nk*Nk,Nlat,Nspin,Nspin,Norb,Norb,Lmats)); Gkmats=zero
        do ik=1,Nk*Nk*Nk
           call dmft_gk_matsubara(Comm,Hk(:,:,ik),1.d0/(Nk*Nk*Nk),Gkmats(ik,:,:,:,:,:,:),Smats)
        enddo
        do i=1,Lmats
          do ik=1,Nk*Nk*Nk
              Ekm(i)=Ekm(i)+trace(matmul(Hk(:,:,ik),nnn2lso_reshape(Gkmats(ik,:,:,:,:,:,i),Nlat,Nspin,Norb)))/(Nk*Nk*Nk)
           enddo
        enddo
        Ek=sum(Ekm)/beta
        write(LOGfile,*) "  Ekin:",Ek
        if(master)then
           open(unit=106,file="Ekin_all.dat",status="unknown",action="write",position="append")
           write(106,'(I5,F15.7)')iloop, Ek
           close(106)
        endif
        deallocate(Ekm,Gkmats)
     endif
     !
     if(master)call end_loop
     !
  enddo
  !
  !
  !##################    POST-PROCESSING     ##################
  !
  !
  !------   compute Ekin   ------
  if(master.and.computeEk)then
     open(unit=106,file="Ekin_last.dat",status="unknown",action="write",position="rewind")
     write(106,'(F15.7)')Ek
     close(106)
  endif
  !
  !------ compute Gloc(wr) ------
  if(nread==0.d0)then
     call dmft_gloc_realaxis(Comm,Hk,Wtk,Greal,Sreal,mpi_split='k')
  else
     if(master) then
        allocate(Gloc(NlNsNo,NlNsNo,Lreal));Gloc=zero
        allocate(zeta(NlNsNo,NlNsNo))      ;zeta=zero
        do i=1,Lreal
           zeta=zero
           do ik=1,Nk*Nk*Nk
              zeta=Hk(:,:,ik)+nnn2lso_reshape(Sreal(:,:,:,:,:,i),Nlat,Nspin,Norb)
              Gloc(:,:,i)=Gloc(:,:,i) + inverse_g0k(dcmplx(wr(i),eps),zeta,Nlat,xmu)/(Nk*Nk*Nk)
           enddo
          Greal(:,:,:,:,:,i)=lso2nnn_reshape(Gloc(:,:,i),Nlat,Nspin,Norb)
        enddo
     endif
  endif
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=6)
  call rotate_local_funct(Greal,U)
  call dmft_print_gf_realaxis(Greal,"Gloc_rot",iprint=6)
  !
  !------   compute Bands  ------
  if(master)call build_eigenbands("LVO_hr.dat","Bands.dat","Hk_path.dat","Kpoints_path.dat",Sreal)
  !
  call finalize_MPI()
  !
  !
  !
contains
  !
  !
  !
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Perform the symmetry operations on the Sigma in the diagonal basis
  !         coming from the solver
  !+------------------------------------------------------------------------------------------+!
  subroutine sigma_symmetrization()
    implicit none
    complex(8),allocatable,dimension(:,:)             ::   Sigma_in
    complex(8),allocatable,dimension(:,:)             ::   Sigma_out
    !
    write(LOGfile,*) "Self-energy symmetrization"
    !
    allocate(Sigma_in(NlNsNo,NlNsNo));Sigma_in=zero
    allocate(Sigma_out(NlNsNo,NlNsNo));Sigma_out=zero
    !
    if (lattice_flag)then
       if(ed_para)then
          !
          !I'm plugging impS in the neighboring site same plane no spin-flip
          do ilayer=1,Nlayer
             write(LOGfile,'(2(A,I5))') " plugghing impS nr.",ilayer," into ilat nr. ",2*ilayer-1,"&",2*ilayer
             !site 1 in plane - same spin - same layer
             Smats(2*ilayer-1,:,:,:,:,:)=Smats_hetero(ilayer,:,:,:,:,:)
             Sreal(2*ilayer-1,:,:,:,:,:)=Sreal_hetero(ilayer,:,:,:,:,:)
             !site 2 in plane - same spin - same layer
             Smats(2*ilayer,:,:,:,:,:)  =Smats_hetero(ilayer,:,:,:,:,:)
             Sreal(2*ilayer,:,:,:,:,:)  =Sreal_hetero(ilayer,:,:,:,:,:)
          enddo
          !
       elseif(.not.ed_para)then
          !
          do ilayer=1,Nlayer
             write(LOGfile,'(2(A,I5))') " plugghing impS nr.",ilayer," into ilat nr. ",2*ilayer-1
             !site 1 in plane - same spin - same layer
             Smats(2*ilayer-1,:,:,:,:,:)=Smats_hetero(ilayer,:,:,:,:,:)
             Sreal(2*ilayer-1,:,:,:,:,:)=Sreal_hetero(ilayer,:,:,:,:,:)
             write(LOGfile,'(2(A,I5))') " plugghing spin-flipped impS nr.",ilayer," into ilat nr. ",2*ilayer
             !site 2 in plane - flip spin - same layer
             Smats(2*ilayer,1,1,:,:,:)  =Smats_hetero(ilayer,2,2,:,:,:)
             Smats(2*ilayer,2,2,:,:,:)  =Smats_hetero(ilayer,1,1,:,:,:)
             Sreal(2*ilayer,1,1,:,:,:)  =Sreal_hetero(ilayer,2,2,:,:,:)
             Sreal(2*ilayer,2,2,:,:,:)  =Sreal_hetero(ilayer,1,1,:,:,:)
          enddo
          !
       endif
    else
       if(ed_para)then
          !
          do ilat=1,Nlat
             Smats(ilat,:,:,:,:,:)=Smats_single
             Sreal(ilat,:,:,:,:,:)=Sreal_single
          enddo
          !
       elseif(.not.ed_para)then
          if(z_symmetry=="ANTIFERRO")then
             write(LOGfile,'(A,I5)') " AFM in all directions"
             !
             Smats(1,:,:,:,:,:)=Smats_single
             Smats(4,:,:,:,:,:)=Smats(1,:,:,:,:,:)
             Smats(2,1,1,:,:,:)=Smats_single(2,2,:,:,:)
             Smats(2,2,2,:,:,:)=Smats_single(1,1,:,:,:)
             Smats(3,:,:,:,:,:)=Smats(2,:,:,:,:,:)
             !
          elseif(z_symmetry=="FERRO")then
             write(LOGfile,'(A,I5)') " AFM in plane - ferro between planes"
             !
             Smats(1,:,:,:,:,:)=Smats_single
             Smats(3,:,:,:,:,:)=Smats(1,:,:,:,:,:)
             Smats(2,1,1,:,:,:)=Smats_single(2,2,:,:,:)
             Smats(2,2,2,:,:,:)=Smats_single(1,1,:,:,:)
             Smats(4,:,:,:,:,:)=Smats(2,:,:,:,:,:)
             !
          endif
       endif
    endif
    !
    deallocate(Sigma_in,Sigma_out)
    call rotate_local_funct(Smats,Udag)
    call rotate_local_funct(Sreal,Udag)
    !
  end subroutine sigma_symmetrization
  !
  !
  !
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Read the Non interacting Hamiltonian from  file
  !         Also this is just for testing the correct interface with the 
  !         translator of the W90 output
  !         The re-ordering part can be used or not, depending on what the user of W90 did.
  !+------------------------------------------------------------------------------------------+!
  subroutine read_myhk(fileHR,fileHk,fileHloc,fileKpoints,local_diagonal_basis)
    implicit none
    character(len=*)            ,intent(in)           ::   fileHR
    character(len=*)            ,intent(in)           ::   fileHk
    character(len=*)            ,intent(in)           ::   fileHloc
    character(len=*)            ,intent(in)           ::   fileKpoints
    logical                     ,intent(in)           ::   local_diagonal_basis
    logical                                           ::   IOfile
    integer                                           ::   P(NlNsNo,NlNsNo)
    real(8)                                           ::   mu,bk_x(3),bk_y(3),bk_z(3)
    real(8)         ,allocatable                      ::   Aw(:,:)
    integer         ,allocatable                      ::   Nkvec(:)
    real(8)         ,allocatable                      ::   Kvec(:,:)
    complex(8)      ,allocatable                      ::   Hloc(:,:)
    complex(8)      ,allocatable                      ::   Gmats(:,:,:),Greal(:,:,:)
    complex(8)      ,allocatable                      ::   Gso(:,:,:,:,:,:)
    !
    !
    bk_x = [1.d0,0.d0,0.d0]*2*pi
    bk_y = [0.d0,1.d0,0.d0]*2*pi
    bk_z = [0.d0,0.d0,1.d0]*2*pi
    call TB_set_bk(bk_x,bk_y,bk_z)
    !
    call Hk_order(P)
    !
    Lk=Nk*Nk*Nk
    if(master)write(LOGfile,*)" Bulk tot k-points:",Lk
    allocate(Hloc(NlNsNo,NlNsNo));Hloc=zero
    allocate(Kvec(Lk,3));Kvec=0d0
    allocate(Nkvec(3));Nkvec=0
    Nkvec=[Nk,Nk,Nk]
    !
    inquire(file=fileHk,exist=IOfile)
    !
    if(IOfile)then
       write(LOGfile,*) " Reading existing Hk from:  ",fileHk
       call TB_read_hk(Hk,fileHk,Nspin*Norb*Nlat,1,1,Nlat,Nkvec,Kvec)
    else
       write(LOGfile,*) " Transforming HR from:  ",fileHR
       call TB_hr_to_hk(Hk,fileHR,Nspin,Norb,Nlat,Nkvec,P,Kvec,fileHk,fileKpoints)
    endif
    !
    Hloc=zero
    Hloc=sum(Hk,dim=3)/Lk
    call TB_write_Hloc(Hloc,reg(fileHloc//".w90"))
    if(local_diagonal_basis)then
       call build_rotations(Hloc,Hloc_lso,U)
       Udag=transpose(conjg(U))
    else
       Hloc_lso=Hloc
    endif
    Hloc_nnn=lso2nnn_reshape(Hloc_lso,Nlat,Nspin,Norb)
    call TB_write_Hloc(Hloc_lso,reg(fileHloc//".used"))
    write(LOGfile,*) " H(k) and Hloc linked"
    !
    !
    !-----  Build the local GF in the spin-orbital Basis   -----
    if(computeG0loc)then
       mu=15.429
       !
       !matsu freq
       allocate(Gmats(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lmats));Gmats=zero
       do ik=1,Lk
          do i=1,Lmats
             Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k(xi*wm(i),Hk(:,:,ik),Nlat,mu)/Lk
          enddo
       enddo
       !
       allocate(Gso(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gso=zero
       do i=1,Lmats
          Gso(:,:,:,:,:,i)=lso2nnn_reshape(Gmats(:,:,i),Nlat,Nspin,Norb)
       enddo
       !
       if(master) call dmft_print_gf_matsubara(Gso,"G0loc",iprint=6)
       deallocate(Gso,Gmats)
       !
       !real freq
       allocate(Greal(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal));Greal=zero
       do ik=1,Lk
          do i=1,Lreal
             Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps),Hk(:,:,ik),Nlat,mu)/Lk
          enddo
       enddo
       !
       allocate(Gso(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gso=zero
       do i=1,Lreal
          Gso(:,:,:,:,:,i)=lso2nnn_reshape(Greal(:,:,i),Nlat,Nspin,Norb)
       enddo
       !
       inquire(file="Aw0.dat",exist=IOfile)
       if(.not.IOfile)then
          open(unit=106,file="Aw0.dat",status="unknown",action="write",position="rewind")
          allocate(Aw(Norb,Lreal));Aw=0d0
          do i=1,Lreal
             do iorb=1,Norb
                do ilat=1,Nlat
                   io = iorb + (ilat-1)*Norb*Nspin
                   Aw(iorb,i)=Aw(iorb,i)-aimag(Greal(io,io,i))/pi
                enddo
             enddo
             write(106,'(100F15.7)') wr(i), Aw(:,i), sum(Aw(:,i))
          enddo
          deallocate(Aw)
          close(106)
       endif
       !
       if(master) call dmft_print_gf_realaxis(Gso,"G0loc",iprint=6)
       deallocate(Gso,Greal)
    endif
    !
    deallocate(Hloc,Kvec,Nkvec)
    !
  end subroutine read_myhk
  !
  !
  !
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Rotations that diagonalizes the lattice local hamiltonians
  !+------------------------------------------------------------------------------------------+!
  subroutine build_rotations(Hloc_lso_in,Hloc_lso_out,rot_lso_,Hloc_nnn_out_)
    implicit none
    complex(8),allocatable,intent(in)                 ::   Hloc_lso_in(:,:)
    complex(8),allocatable,intent(out)                ::   Hloc_lso_out(:,:)
    complex(8),allocatable,intent(out),optional       ::   rot_lso_(:,:)
    complex(8),allocatable,intent(out),optional       ::   Hloc_nnn_out_(:,:,:,:,:)
    complex(8),allocatable                            ::   rot_lso(:,:)
    complex(8),allocatable                            ::   Hloc_nnn_out(:,:,:,:,:)
    complex(8),allocatable                            ::   Hloc_nnn_in(:,:,:,:,:)
    complex(8),allocatable                            ::   rot_nnn(:,:,:,:,:)
    complex(8),allocatable                            ::   arg(:,:)
    real(8)   ,allocatable                            ::   eig(:)
    !
    allocate(Hloc_nnn_in(Nlat,Nspin,Nspin,Norb,Norb));Hloc_nnn_in=zero
    allocate(rot_nnn(Nlat,Nspin,Nspin,Norb,Norb));rot_nnn=zero
    allocate(Hloc_lso_out(NlNsNo,NlNsNo));Hloc_lso_out=zero
    allocate(arg(Norb,Norb));arg=zero
    allocate(eig(Norb));eig=0.d0
    allocate(rot_lso(NlNsNo,NlNsNo));rot_lso=zero
    allocate(Hloc_nnn_out(Nlat,Nspin,Nspin,Norb,Norb));Hloc_nnn_out=zero
    !
    Hloc_nnn_in=lso2nnn_reshape(Hloc_lso_in,Nlat,Nspin,Norb)
    !
    do ilat=1,Nlat
       do ispin=1,Nspin
          arg=zero;eig=0.0d0
          arg=Hloc_nnn_in(ilat,ispin,ispin,:,:)
          call eigh(arg,eig,'V','U')
          rot_nnn(ilat,ispin,ispin,:,:)=arg
       enddo
    enddo
    !
    !outputs
    rot_lso=nnn2lso_reshape(rot_nnn,Nlat,Nspin,Norb)
    Hloc_lso_out=matmul(transpose(conjg(rot_lso)),matmul(Hloc_lso_in,rot_lso))
    Hloc_nnn_out=lso2nnn_reshape(Hloc_lso_out,Nlat,Nspin,Norb)
    if(present(rot_lso_))rot_lso_=rot_lso
    if(present(Hloc_nnn_out_))Hloc_nnn_out_=Hloc_nnn_out
    !
    deallocate(Hloc_nnn_in,rot_nnn,arg,eig,rot_lso,Hloc_nnn_out)
    !
  end subroutine build_rotations
  !
  !
  !
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Rotate function with specific rotation
  !+------------------------------------------------------------------------------------------+!
  subroutine rotate_local_funct(funct,rotation)
    implicit none
    complex(8),allocatable,intent(inout)              ::   funct(:,:,:,:,:,:)
    complex(8),allocatable,intent(in)                 ::   rotation(:,:)
    complex(8),allocatable                            ::   funct_in(:,:)
    complex(8),allocatable                            ::   funct_out(:,:)
    !
    allocate(funct_in(NlNsNo,NlNsNo));funct_in=zero
    allocate(funct_out(NlNsNo,NlNsNo));funct_out=zero
    Lfreq=size(funct,6)
    !
    do ifreq=1,Lfreq
       funct_in=zero;funct_out=zero
       funct_in=nnn2lso_reshape(funct(:,:,:,:,:,ifreq),Nlat,Nspin,Norb);funct(:,:,:,:,:,ifreq)=zero
       funct_out=matmul(transpose(conjg(rotation)),matmul(funct_in,rotation))
       funct(:,:,:,:,:,ifreq)=lso2nnn_reshape(funct_out,Nlat,Nspin,Norb)
    enddo
    deallocate(funct_in,funct_out)
    !
  end subroutine rotate_local_funct
  !
  !
  !
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: solve H(k) along path in the BZ.
  !+------------------------------------------------------------------------------------------+!
  subroutine build_eigenbands(fileHR,fileband,fileHk_path,fileKpoints_path,Sreal_)
    implicit none
    character(len=*)            ,intent(in)           ::   fileHR
    character(len=*)            ,intent(in),optional  ::   fileband,fileHk_path,fileKpoints_path
    complex(8)      ,allocatable,intent(in),optional  ::   Sreal_(:,:,:,:,:,:)
    integer                                           ::   Npts,Nkpathread
    integer                                           ::   P(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    real(8)         ,allocatable                      ::   kpath(:,:),kgrid(:,:),Scorr(:,:)
    complex(8)      ,allocatable                      ::   Hkpath(:,:,:)
    complex(8)      ,allocatable                      ::   Gkreal(:,:,:,:,:,:,:)
    complex(8)      ,allocatable                      ::   Gkr(:,:,:,:)
    type(rgb_color) ,allocatable                      ::   colors(:),colors_orb(:)
    logical                                           ::   IOfile
    !
    call Hk_order(P)
    !
    allocate(colors(NlNsNo))
    allocate(colors_orb(Norb))
    colors_orb=[red1,green1,blue1]
    do i=1,Nspin*Nlat
       colors(1+(i-1)*Norb:Norb+(i-1)*Norb)=colors_orb
    enddo
    !
    write(LOGfile,*)
    write(LOGfile,*)"Build bulk H(k) along the path M-R-G-M-X-G-X"
    write(LOGfile,*)
    Npts = 7
    Lk=(Npts-1)*Nkpath
    allocate(kpath(Npts,3))
    kpath(1,:)=kpoint_M1
    kpath(2,:)=kpoint_R
    kpath(3,:)=kpoint_Gamma
    kpath(4,:)=kpoint_M1
    kpath(5,:)=kpoint_X1 
    kpath(6,:)=kpoint_Gamma
    kpath(7,:)=kpoint_X1
    !
    allocate(kgrid(Lk,3))  ;kgrid=0d0
    allocate(Hkpath(NlNsNo,NlNsNo,Lk));Hkpath=zero
    allocate(Gkreal(Lk,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
    allocate(Scorr(NlNsNo,NlNsNo));Scorr=0d0
    if(present(Sreal_))Scorr=real(nnn2lso_reshape(Sreal_(:,:,:,:,:,1),Nlat,Nspin,Norb),8)
    !
    inquire(file=fileHk_path,exist=IOfile)
    !
    if(IOfile)then
       write(LOGfile,*) "   Reading existing Hkpath on: ",fileHk_path

       call TB_read_hk(Hkpath,fileHk_path,Nspin*Norb*Nlat,Nkpathread,kpath,kgrid)
       if(Nkpathread.ne.Nkpath) stop "Eigenbands wrong Nkpath readed"
    else
       write(LOGfile,*) "   Solving model on path"
       call TB_solve_model(   fileHR,Nspin,Norb,Nlat,kpath,Nkpath,colors               &
                          ,   [character(len=20) ::'M', 'R', 'G', 'M', 'X', 'G', 'X']  &
                          ,   P                                                        &
                          ,   fileband                                                 &
                          ,   fileHk_path                                              &
                          ,   fileKpoints_path                                         &
                          ,   Scorr                                                    &
                          ,   Hkpath                                                   &
                          ,   kgrid                                                    )
       !
    endif
    !
    allocate(Gkr(Lk,Nspin*Nlat*Norb,Nspin*Nlat*Norb,Lreal));Gkr=zero
    do ik=1,Lk
       call dmft_gk_realaxis(Hkpath(:,:,ik),1.d0/Lk,Gkreal(ik,:,:,:,:,:,:),Sreal)
       !faccio questa cosa qui sotto per separare per bene i due blocchi di spin
       do ifreq=1,Lreal
          Gkr(ik,:,:,ifreq)=matmul(P,matmul(nnn2lso_reshape(Gkreal(ik,:,:,:,:,:,ifreq),Nlat,Nspin,Norb),transpose(P)))
       enddo
    enddo
    !
    open(unit=106,file='Akw_s1.dat',status='unknown',action='write',position='rewind')
    open(unit=107,file='Akw_s2.dat',status='unknown',action='write',position='rewind')
    do ifreq=1,Lreal
       write(106,'(9000F18.12)')wr(ifreq),(30.*trace(-aimag(Gkr(ik,1:Nlat*Norb,1:Nlat*Norb,ifreq))/pi)+10.*ik,ik=1,Lk)
       write(107,'(9000F18.12)')wr(ifreq),(30.*trace(-aimag(Gkr(ik,1+Nlat*Norb:Nlat*Nspin*Norb,1+Nlat*Norb:Nlat*Nspin*Norb,ifreq))/pi)+10.*ik,ik=1,Lk)
    enddo
    close(106)
    close(107)
    write(LOGfile,*)"Im done on the path"
    !
    deallocate(kgrid,Hkpath,Gkreal,Scorr,Gkr,colors,colors_orb)
    !
  end subroutine build_eigenbands
  !
  !
  !
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Put in proper order the input coming from w90. Depends on the w90 users.
  !+------------------------------------------------------------------------------------------+!
  subroutine Hk_order(Porder)
    implicit none
    integer   ,intent(out)                       :: Porder(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    integer   ,allocatable,dimension(:,:)        :: shift2!,shift1,shift3
!   integer                                      :: P1(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    integer                                      :: P2(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
!   integer                                      :: P3(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    integer                                      :: io1,jo1,io2,jo2
    !
    !-----  Ordering 1: same orbital position for each Nlat block   -----
!    allocate(shift1(2,2));shift1=0
!    shift1(1,:)=[1,2]
!    shift1(2,:)=[7,8]
!    P1=0;P1=int(eye(Nlat*Nspin*Norb))
!    do i=1,size(shift1,1)
!       do j=1,2
!          P1(shift1(i,j),shift1(i,j))=0
!       enddo
!       P1(shift1(i,1),shift1(i,2))=1
!       P1(shift1(i,2),shift1(i,1))=1
!    enddo
!    P1(1+Nlat*Norb:Nlat*Norb*Nspin,1+Nlat*Norb:Nlat*Norb*Nspin)=P1(1:Nlat*Norb,1:Nlat*Norb)
    !
    !-----  Ordering 2: as used in the code [[[Norb],Nspin],Nlat]   -----
    ndx=0
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             !input
             io1 = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
             !output
             io2 = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             !
             if(io1.ne.io2) ndx=ndx+1
          enddo
       enddo
    enddo
    allocate(shift2(ndx,2));shift2=0
    ndx=0
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             !input
             io1 = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
             !output
             io2 = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin
             !
             if(io1.ne.io2) then
                ndx=ndx+1
                shift2(ndx,:)=[io1,io2]
             endif
          enddo
       enddo
    enddo
    P2=0;P2=int(eye(Nlat*Nspin*Norb))
    do i=1,size(shift2,1)
       do j=1,2
          P2(shift2(i,j),shift2(i,j))=0
       enddo
       P2(shift2(i,1),shift2(i,2))=1
    enddo
    !
    !-------------  Ordering 3: swapping site 2 with site 3   -----------
!    allocate(shift3(6,2));shift3=0
!    do i=1,Nspin*Norb
!       shift3(i,:)=[Nspin*Norb+i,2*Nspin*Norb+i]
!    enddo
!    ndx=0
!    P3=0;P3=int(eye(Nlat*Nspin*Norb))
!    do i=1,size(shift3,1)
!       do j=1,2
!          P3(shift3(i,j),shift3(i,j))=0
!       enddo
!       P3(shift3(i,1),shift3(i,2))=1
!       P3(shift3(i,2),shift3(i,1))=1
!    enddo
    !
    !--------------------  Global reordering   -------------------------
    !Porder=matmul(P1,P2)
    Porder=P2
    !Porder=matmul(P1,matmul(P2,P3))
    !
  end subroutine Hk_order
  !
  !
  !
  !____________________________________________________________________________________________!
  !                                       Gfs
  !____________________________________________________________________________________________!
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: G0_loc functions
  !+------------------------------------------------------------------------------------------+!
  function inverse_g0k(iw,hk_,Nlat,mu_) result(g0k)
    implicit none
    complex(8),intent(in)                                  :: iw
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)  :: hk_
    real(8),intent(in),optional                            :: mu_
    real(8)                                                :: mu
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)  :: g0k,g0k_tmp
    integer                                                :: i,ndx,Nlat
    integer (kind=4), dimension(6)                         :: ipiv
    integer (kind=1)                                       :: ok
    integer (kind=4), parameter                            :: lwork=2000
    complex (kind=8), dimension(lwork)                     :: work
    real    (kind=8), dimension(lwork)                     :: rwork
    !
    mu=0.d0
    if(present(mu_))mu=mu_
    g0k=zero;g0k_tmp=zero
    !
    g0k=(iw+mu)*eye(Nlat*Nspin*Norb)-hk_
    g0k_tmp=g0k
    !
    call inv(g0k)
    call inversion_test(g0k,g0k_tmp,1.e-6,Nlat)
  end function inverse_g0k
  !
  !
  !
  !____________________________________________________________________________________________!
  !                                     utilities
  !____________________________________________________________________________________________!
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Inversion test
  !+------------------------------------------------------------------------------------------+!
  subroutine inversion_test(A,B,tol,Nlat)
    implicit none
    integer (kind=4), intent(in)   ::   Nlat
    complex (kind=8), intent(in)   ::   A(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    complex (kind=8), intent(in)   ::   B(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    real    (kind=4), intent(in)   ::   tol
    real    (kind=4)               ::   error
    integer (kind=2)               ::   dime

    if (size(A).ne.size(B)) then
       write(LOGfile,*) "Matrices not equal cannot perform inversion test"
       stop
    endif
    dime=maxval(shape(A))
    error=abs(float(dime)-real(sum(matmul(A,B))))
    if (error.gt.tol) write(LOGfile,*) "inversion test fail",error
  end subroutine inversion_test
  !
  !
  !
end program ed_LVO_hetero








!    kpath( 1,:)=kpoint_Gamma
!    kpath( 2,:)=kpoint_X1
!    kpath( 3,:)=kpoint_M1
!    kpath( 4,:)=kpoint_X2
!    kpath( 5,:)=kpoint_Gamma
!    kpath( 6,:)=kpoint_X3
!    kpath( 7,:)=kpoint_M3
!    kpath( 8,:)=kpoint_R
!    kpath( 9,:)=kpoint_M2
!    kpath(10,:)=kpoint_X1
!    kpath(11,:)=kpoint_M1
!    kpath(12,:)=kpoint_X2
!    kpath(13,:)=kpoint_Gamma
!    kpath(14,:)=kpoint_X3
!    kpath(15,:)=kpoint_M3
!    kpath(16,:)=kpoint_R
!    kpath(17,:)=kpoint_M2
!    kpath(18,:)=kpoint_X3
