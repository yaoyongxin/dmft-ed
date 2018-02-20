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
  integer                                        :: iloop,i,j,ndx
  integer                                        :: Nlat,ilat
  integer                                        :: ilayer,Nlayer
  integer                                        :: io,jo,ik
  integer                                        :: iorb,jorb,ispin,jspin
  logical                                        :: converged
  real(8)                                        :: wmixing
  character(len=60)                              :: finput
  character(len=32)                              :: hkfile
  character(len=32)                              :: geometry
  character(len=32)                              :: z_symmetry
  !Mpi:
  integer                                        :: comm,rank,ier
  logical                                        :: master
  !Bath:
  integer                                        :: Nb,unit
  real(8)   ,allocatable,dimension(:,:)          :: Bath
  real(8)   ,allocatable,dimension(:)            :: Bath_single
  !Local lattice functions:
  complex(8),allocatable,dimension(:,:,:,:,:,:)  :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:)  :: Gmats,Greal
  !Weiss&Hybridization functions
  complex(8),allocatable,dimension(:,:,:,:,:,:)  :: field
  complex(8),allocatable,dimension(:,:,:,:,:,:)  :: field_red,field_red_old
  complex(8),allocatable,dimension(:,:,:,:,:,:,:):: field_red_mem
  complex(8),allocatable,dimension(:,:,:,:,:)    :: field_single,field_single_old
  complex(8),allocatable,dimension(:,:,:,:,:,:)  :: field_single_mem
  !reduced functions
  complex(8),allocatable,dimension(:,:,:,:,:,:)  :: Smats_red,Sreal_red
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Smats_single,Sreal_single
  !Hmiltonian input:
  integer                                        :: Nk,Nkpath
  real(8)   ,allocatable,dimension(:)            :: Wtk
  complex(8),allocatable,dimension(:,:,:)        :: Hk
  complex(8),allocatable,dimension(:,:)          :: Hloc_nso
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Hloc_nnn
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Hloc_red_nnn
  complex(8),allocatable,dimension(:,:,:,:)      :: Hloc_single_nn
  !custom variables for convergence test:
  complex(8),allocatable,dimension(:)            :: conv_funct
  complex(8),allocatable,dimension(:,:)          :: conv_funct_mem
  !custom variables for chempot search:
  character(len=32)                              :: ed_file_suffix
  logical                                        :: converged_n,upprshft
  integer                                        :: conv_n_loop=0,Nlat_max
  real(8)                                        :: xmu_start
  real(8)                                        :: dw,sumdens,xmu_old
  real(8)   ,allocatable,dimension(:)            :: wr,wm
  real(8)   ,allocatable,dimension(:,:)          :: orb_dens
  real(8)   ,allocatable,dimension(:)            :: orb_dens_red
  logical                                        :: look4n=.true.
  !custom variables misc:
  logical                                        :: computeG0loc
  logical                                        :: lattice_flag=.true.
  logical                                        :: bulk_magsym
  complex(8),allocatable,dimension(:,:,:)        :: Gloc
  complex(8),allocatable,dimension(:,:)          :: zeta
  !
  !#########   MPI INITIALIZATION   #########
  !
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  master = get_Master_MPI(comm)
  !
  !#########    VARIABLE PARSING    #########
  !
  call parse_cmd_variable(finput,           "FINPUT",              default='inputED_LVO.in')
  call parse_input_variable(hkfile,         "HKFILE",finput,       default="Hk.dat")
  call parse_input_variable(nk,             "NK",finput,           default=10)
  call parse_input_variable(NLAT,           "NLAT",finput,         default=4)
  call parse_input_variable(nkpath,         "NKPATH",finput,       default=20)
  call parse_input_variable(wmixing,        "WMIXING",finput,      default=0.5d0)
  call parse_input_variable(computeG0loc,   "COMPUTEG0loc",finput, default=.false.)
  call parse_input_variable(geometry,       "GEOMETRY",finput,     default="bulk")
  call parse_input_variable(z_symmetry,     "ZSYMMETRY",finput,    default="FERRO")
  call parse_input_variable(bulk_magsym,    "BULKMAGSYM",finput,   default=.false.)
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
  call add_ctrl_var(ed_file_suffix,"ed_file_suffix")
  !
  geometry=reg(geometry)
  z_symmetry=reg(z_symmetry)
  if (geometry=="bulk".and.ed_para)    lattice_flag=.false.
  if (geometry=="bulk".and.bulk_magsym)lattice_flag=.false.
  if (geometry=="bulk".and.Nlat/=4) stop
  xmu_start=xmu
  Nlayer=Nlat/2
  !
  !##################       ALLOCATION       !##################
  !
  !1)Global qualntities
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));            Smats=zero
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));            Gmats=zero
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));            Sreal=zero
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));            Greal=zero
  allocate(field(Nlat,Nspin,Nspin,Norb,Norb,Lmats));            field=zero
  !
  allocate(Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Nk*Nk*Nk));       Hk=zero
  allocate(Hloc_nso(Nlat*Nspin*Norb,Nlat*Nspin*Norb));          Hloc_nso=zero
  allocate(Hloc_nnn(Nlat,Nspin,Nspin,Norb,Norb));               Hloc_nnn=zero
  !
  !2)single site quantities (used for bulk+para)
  allocate(Hloc_single_nn(Nspin,Nspin,Norb,Norb));              Hloc_single_nn=zero
  allocate(Smats_single(Nspin,Nspin,Norb,Norb,Lmats));          Smats_single=zero
  allocate(Sreal_single(Nspin,Nspin,Norb,Norb,Lreal));          Sreal_single=zero
  allocate(field_single(Nspin,Nspin,Norb,Norb,Lmats));          field_single=zero
  allocate(field_single_old(Nspin,Nspin,Norb,Norb,Lmats));      field_single_old=zero
  allocate(field_single_mem(Nspin,Nspin,Norb,Norb,Lmats,3));    field_single_mem=zero
  !
  !3)Layer quantities (used for bulk+afm,hetero+para,hetero+afm)
  allocate(Hloc_red_nnn(Nlayer,Nspin,Nspin,Norb,Norb));         Hloc_red_nnn=zero
  allocate(Smats_red(Nlayer,Nspin,Nspin,Norb,Norb,Lmats));      Smats_red=zero
  allocate(Sreal_red(Nlayer,Nspin,Nspin,Norb,Norb,Lreal));      Sreal_red=zero
  allocate(field_red(Nlayer,Nspin,Nspin,Norb,Norb,Lmats));      field_red=zero
  allocate(field_red_old(Nlayer,Nspin,Nspin,Norb,Norb,Lmats));  field_red_old=zero
  allocate(field_red_mem(Nlayer,Nspin,Nspin,Norb,Norb,Lmats,3));field_red_mem=zero
  !
  allocate(conv_funct(Lmats));                                  conv_funct=zero
  allocate(wr(Lreal));wr=0.0d0;                                 wr=linspace(wini,wfin,Lreal,mesh=dw)
  allocate(wm(Lmats));wm=0.0d0;                                 wm = pi/beta*real(2*arange(1,Lmats)-1,8)
  allocate(Wtk(Nk*Nk*Nk));                                      Wtk=1.d0/(Nk*Nk*Nk)
  !
  Lfit=min(int((Uloc(1)+3.)*(beta/pi))+100,Lmats)
  if(master)write(LOGfile,'(a12,I6,2(a12,F10.3))')"Lfit:",Lfit,"iwmax:",(pi/beta)*(2*Lfit-1),"U+2D:",Uloc(1)+3.
  !
  !##################        BUILD Hk        ##################
  !
  call read_myhk("LVO_hr.dat","Hk.dat","Hloc.dat","Kpoints.dat")
  Hloc_single_nn=Hloc_nnn(1,:,:,:,:)
  do ilayer=1,Nlayer
     Hloc_red_nnn(ilayer,:,:,:,:)=Hloc_nnn(2*ilayer-1,:,:,:,:)
  enddo
  !
  !call ed_read_impSigma_lattice(Nlat)
  !call ed_get_Sreal(Sreal,Nlat)
  !if(master)call build_eigenbands("LVO_hr.dat","Bands","Hk_path.dat","Kpoints_path.dat",Sreal)
  !stop
  !
  !##################          BATH          ##################
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
  !##################      INIT SOLVER       ##################
  !
  if (lattice_flag)then
     call ed_init_solver(Comm,Bath,Hloc_red_nnn)
  else
     call ed_init_solver(Comm,Bath_single,Hloc_single_nn)
  endif
  !
  !##################          DMFT          ##################
  !
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !------  solve impurity  ------
     if (lattice_flag)then
        call ed_solve(comm,Bath,Hloc_red_nnn)
     else
        call ed_solve(comm,Bath_single,Hloc_single_nn)
     endif
     !
     !------    get sigmas    ------
     if (lattice_flag)then
        call ed_get_sigma_matsubara(Smats_red)
        call ed_get_sigma_real(Sreal_red)
        if(ed_para)then
           !I'm plugging impS in the neighboring site same plane no spin-flip
           do ilayer=1,Nlayer
              write(LOGfile,'(2(A,I5))') " plugghing impS nr.",ilayer," into ilat nr. ",2*ilayer-1,"&",2*ilayer
              !site 1 in plane - same spin - same layer
              Smats(2*ilayer-1,:,:,:,:,:)=Smats_red(ilayer,:,:,:,:,:)
              Sreal(2*ilayer-1,:,:,:,:,:)=Sreal_red(ilayer,:,:,:,:,:)
              !site 2 in plane - same spin - same layer
              Smats(2*ilayer,:,:,:,:,:)  =Smats_red(ilayer,:,:,:,:,:)
              Sreal(2*ilayer,:,:,:,:,:)  =Sreal_red(ilayer,:,:,:,:,:)
           enddo
        elseif(.not.ed_para)then
           !I'm plugging impS in the neighboring site same plane flipping the spin
           do ilayer=1,Nlayer
              write(LOGfile,'(2(A,I5))') " plugghing impS nr.",ilayer," into ilat nr. ",2*ilayer-1
              !site 1 in plane - same spin - same layer
              Smats(2*ilayer-1,:,:,:,:,:)=Smats_red(ilayer,:,:,:,:,:)
              Sreal(2*ilayer-1,:,:,:,:,:)=Sreal_red(ilayer,:,:,:,:,:)
              write(LOGfile,'(2(A,I5))') " plugghing spin-flipped impS nr.",ilayer," into ilat nr. ",2*ilayer
              !site 2 in plane - flip spin - same layer
              Smats(2*ilayer,1,1,:,:,:)  =Smats_red(ilayer,2,2,:,:,:)
              Smats(2*ilayer,2,2,:,:,:)  =Smats_red(ilayer,1,1,:,:,:)
              Sreal(2*ilayer,1,1,:,:,:)  =Sreal_red(ilayer,2,2,:,:,:)
              Sreal(2*ilayer,2,2,:,:,:)  =Sreal_red(ilayer,1,1,:,:,:)
           enddo
        endif
     else
        call ed_get_sigma_matsubara(Smats_single)
        call ed_get_sigma_real(Sreal_single)
        if(ed_para)then
           !I'm plugging same impS in all the sites no spin-flip
           do ilat=1,Nlat
              write(LOGfile,'(A,I5)') " plugghing impS into ilat nr. ",ilat
              !site 1 in plane - same spin - same sinlge sigma
              Smats(ilat,:,:,:,:,:)=Smats_single
              Sreal(ilat,:,:,:,:,:)=Sreal_single
           enddo
        elseif(.not.ed_para)then
           if(z_symmetry=="ANTIFERRO")then
              write(LOGfile,'(A,I5)') " AFM in all directions"
              !site 1 == reference
              Smats(1,:,:,:,:,:)=Smats_single
              Sreal(1,:,:,:,:,:)=Sreal_single
              !site 2 flip-reference
              Smats(2,1,1,:,:,:)=Smats_single(2,2,:,:,:)
              Smats(2,2,2,:,:,:)=Smats_single(1,1,:,:,:)
              Sreal(2,1,1,:,:,:)=Sreal_single(2,2,:,:,:)
              Sreal(2,2,2,:,:,:)=Sreal_single(1,1,:,:,:)
              !site 3 flip-reference
              Smats(3,1,1,:,:,:)=Smats_single(2,2,:,:,:)
              Smats(3,2,2,:,:,:)=Smats_single(1,1,:,:,:)
              Sreal(3,1,1,:,:,:)=Sreal_single(2,2,:,:,:)
              Sreal(3,2,2,:,:,:)=Sreal_single(1,1,:,:,:)
              !site 4 == reference
              Smats(4,:,:,:,:,:)=Smats_single
              Sreal(4,:,:,:,:,:)=Sreal_single
           elseif(z_symmetry=="FERRO")then
              write(LOGfile,'(A,I5)') " AFM in plane - ferro between planes"
              !I'm plugging same impS in all the sites flipping the spin
              do ilayer=1,Nlayer
                 write(LOGfile,'(A,I5)') " plugghing impS into ilat nr. ",2*ilayer-1
                 !site 1 in plane - same spin - same sinlge sigma
                 Smats(2*ilayer-1,:,:,:,:,:)=Smats_single
                 Sreal(2*ilayer-1,:,:,:,:,:)=Sreal_single
                 write(LOGfile,'(A,I5)') " plugghing spin-flipped impS into ilat nr. ",2*ilayer
                 !site 2 in plane - flip spin - same sinlge sigma
                 Smats(2*ilayer,1,1,:,:,:)  =Smats_single(2,2,:,:,:)
                 Smats(2*ilayer,2,2,:,:,:)  =Smats_single(1,1,:,:,:)
                 Sreal(2*ilayer,1,1,:,:,:)  =Sreal_single(2,2,:,:,:)
                 Sreal(2*ilayer,2,2,:,:,:)  =Sreal_single(1,1,:,:,:)
              enddo
           endif
        endif
     endif
     !
     !------  get local Gf's  ------
     call dmft_gloc_matsubara(Comm,Hk,Wtk,Gmats,Smats,mpi_split='k')
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=6)
     !
     !------    get field     ------
     if(cg_scheme=='weiss')then
        call dmft_weiss(Gmats,Smats,field,Hloc_nnn)
        call dmft_print_gf_matsubara(field,"Weiss",iprint=6)
     elseif(cg_scheme=='delta')then
        call dmft_delta(Gmats,Smats,field,Hloc_nnn)
        call dmft_print_gf_matsubara(field,"Delta",iprint=6)
     endif
     !
     !------  reduce field    ------
     field_single=field(1,:,:,:,:,:)
     do ilayer=1,Nlayer
        field_red(ilayer,:,:,:,:,:)=field(2*ilayer-1,:,:,:,:,:)
     enddo
     !
     !------    mix field     ------
     if(iloop>1)then
        field_red    = wmixing*field_red    + (1.d0-wmixing)*field_red_old
        field_single = wmixing*field_single + (1.d0-wmixing)*field_single_old
     endif
     field_single_old=field_single
     field_red_old=field_red
     !
     !------    mem field     ------
     if(lattice_flag)then
        do ndx=2,3
           field_red_mem(:,:,:,:,:,:,ndx-1)=field_red_mem(:,:,:,:,:,:,ndx)
        enddo
        field_red_mem(:,:,:,:,:,:,3)=field_red_old
     else
        do ndx=2,3
           field_single_mem(:,:,:,:,:,ndx-1)=field_single_mem(:,:,:,:,:,ndx)
        enddo
        field_single_mem(:,:,:,:,:,3)=field_single_old
     endif
     !
     !------    fit field     ------
     if (lattice_flag)then
        call ed_chi2_fitgf(Comm,Bath,field_red,Hloc_red_nnn)
     else
        call set_Hloc(Hloc_single_nn)
        call ed_chi2_fitgf(Comm,field_single,bath_single)
     endif
     !
     !each loop operations
     if(master)then
        !
        !e - chemical potential find
        converged_n=.true.
        sumdens=0d0
        xmu_old=xmu
        if(lattice_flag)then
           allocate(orb_dens(Nlayer,Norb));orb_dens=0.d0
           call ed_get_dens(orb_dens,Nlayer)
           do ilayer=1,Nlayer
              sumdens=sumdens+sum(orb_dens(ilayer,:))/float(Nlayer)
              write(LOGfile,*)"  Nlat:",ilayer,orb_dens(ilayer,:)
           enddo
           deallocate(orb_dens)
        else
           allocate(orb_dens_red(Norb));orb_dens_red=0.d0
           call ed_get_dens(orb_dens_red)
           write(LOGfile,*)"  Nlat:",1,orb_dens_red(:)
           sumdens=sum(orb_dens_red)
           deallocate(orb_dens_red)
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
        !f - convergence
        write(LOGfile,*)
        write(LOGfile,*) "   ------------------- convergence --------------------"
        if (lattice_flag)then
           if(ed_para)then
              do i=1,Lmats
                 conv_funct(i)=sum(nnn2lso_reshape(field_red(:,:,:,:,:,i),Nlayer,Nspin,Norb))
              enddo
           else
              do ndx=1,3
                 do i=1,Lmats
                    conv_funct_mem(i,ndx)=sum(nnn2lso_reshape(field_red_mem(:,:,:,:,:,i,ndx),Nlayer,Nspin,Norb))
                 enddo
              enddo
           endif
        else
           if(ed_para)then
              do i=1,Lmats
                 conv_funct(i)=sum(nn2so_reshape(field_single(:,:,:,:,i),Nspin,Norb))
              enddo
           else
              do ndx=1,3
                 do i=1,Lmats
                    conv_funct_mem(i,ndx)=sum(nn2so_reshape(field_single_mem(:,:,:,:,i,ndx),Nspin,Norb))
                 enddo
              enddo
           endif
        endif
        if(converged_n)then
           if(ed_para)then
              converged = check_convergence(conv_funct,dmft_error,nsuccess,nloop)
           else
              do ndx=1,3
                 converged = check_convergence(conv_funct_mem(:,ndx),dmft_error,nsuccess,nloop,file="error_mem"//str(ndx)//".err",index=ndx)
              enddo
           endif
        endif
        write(LOGfile,'(a35,L3)') "sigma converged",converged
        write(LOGfile,'(a35,L3)') "dens converged",converged_n
        converged = converged .and. converged_n
        write(LOGfile,'(a35,L3)') "total converged",converged
        write(LOGfile,'(a35,I3)') "global iloop",iloop
        write(LOGfile,'(a35,I3)') "times dens is ok",conv_n_loop
        write(LOGfile,*) "   ----------------------------------------------------"
        write(LOGfile,*)
     endif
     !
     call Bcast_MPI(Comm,xmu)
     call Bcast_MPI(Comm,converged)
     call MPI_Barrier(Comm,ier)
     !
     if(master)call end_loop
     !
  enddo
  !
  !
  !#########      compute Ekin      #########
  !
  call dmft_kinetic_energy(Comm,Hk,Wtk,Smats)
  !
  !#########     BUILD Gloc(wr)     #########
  !
  if(nread==0.d0)then
     call dmft_gloc_realaxis(Comm,Hk,Wtk,Greal,Sreal,mpi_split='k')
  else
     if(master) then
        allocate(Gloc(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal));Gloc=zero
        allocate(zeta(Nlat*Nspin*Norb,Nlat*Nspin*Norb))      ;zeta=zero
        do i=1,Lreal
           zeta=zero
           do ik=1,Nk*Nk*Nk
              zeta=Hk(:,:,ik)+nnn2lso_reshape(Sreal(:,:,:,:,:,i),Nlat,Nspin,Norb)
              Gloc(:,:,i)=Gloc(:,:,i) + inverse_g0k(dcmplx(wr(i),eps),zeta,Nlat,xmu)/(Nk*Nk*Nk)
           enddo
          Greal(:,:,:,:,:,i)=lso2nnn_reshape(Gloc(:,:,i),Nlat,Nspin,Norb)
        enddo
        call dmft_print_gf_realaxis(Greal,"Gloc",iprint=6)
     endif
  endif
  !
  !#########    BUILD Hk ON PATH    #########
  !
  if(master)call build_eigenbands("LVO_hr.dat","Bands.dat","Hk_path.dat","Kpoints_path.dat",Sreal)
  !
  call finalize_MPI()
  !
  !
contains


  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Read the Non interacting Hamiltonian from  file
  !         Also this is just for testing the correct interface with the 
  !         translator of the W90 output
  !         The re-ordering part can be used or not, depending on what the user of W90 did.
  !+------------------------------------------------------------------------------------------+!
  subroutine read_myhk(fileHR,fileHk,fileHloc,fileKpoints)
    implicit none
    character(len=*),intent(in)                  :: fileHR
    character(len=*),intent(in)                  :: fileHk
    character(len=*),intent(in)                  :: fileHloc
    character(len=*),intent(in)                  :: fileKpoints
    integer                                      :: ispin,iorb,jspin,jorb,ilat,jlat
    integer                                      :: ik,Lk,i,j
    integer                                      :: io1,jo1,io2,jo2,ndx
    real(8)                                      :: mu
    integer   ,allocatable,dimension(:)          :: Nkvec
    real(8)   ,allocatable,dimension(:,:)        :: Kvec
    complex(8),allocatable,dimension(:,:)        :: Hloc
    complex(8),allocatable,dimension(:,:,:)      :: ham_k,ham_k_aux
    integer                                      :: P(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    real(8),dimension(3)                         :: bk_x,bk_y,bk_z
    logical                                      :: IOfile
    complex(8)                                   :: Gmats(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lmats)
    real(8)                                      :: Aw(Norb,Lreal)
    complex(8)                                   :: Greal(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal)
    complex(8),allocatable                       :: Gso(:,:,:,:,:,:)
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
    allocate(ham_k(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lk))     ;ham_k=zero
    allocate(ham_k_aux(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lk)) ;ham_k_aux=zero
    allocate(Hloc(Nlat*Nspin*Norb,Nlat*Nspin*Norb))         ;Hloc=zero
    allocate(Kvec(Lk,3));Kvec=0d0
    allocate(Nkvec(3));Nkvec=0
    Nkvec=[Nk,Nk,Nk]
    !
    inquire(file=fileHk,exist=IOfile)
    !
    if(IOfile)then
       write(LOGfile,*) " Reading existing Hk"
       call TB_read_hk(ham_k_aux,fileHk,Nspin*Norb*Nlat,1,1,Nlat,Nkvec,Kvec)
    else
       write(LOGfile,*) " Transforming HR from:  ",fileHR
       call TB_hr_to_hk(ham_k_aux,fileHR,Nspin,Norb,Nlat,Nkvec,P,Kvec,fileHk,fileKpoints)
    endif
    ham_k=zero;Hloc=zero
    ham_k=ham_k_aux
    Hloc=sum(ham_k,dim=3)/Lk
    call TB_write_Hloc(Hloc,fileHloc)
    ham_k_aux=zero
    !
    !-------------- Linking to external variables --------------
    Hk=ham_k
    Hloc_nso=Hloc
    Hloc_nnn=lso2nnn_reshape(Hloc,Nlat,Nspin,Norb)
    write(LOGfile,*) " H(k) and Hloc linked"
    !
    !-----  Build the local GF in the spin-orbital Basis   -----
    if(computeG0loc)then
       mu=15.429
       !
       !matsu freq
       Gmats=zero
       do ik=1,Lk
          do i=1,Lmats
             Gmats(:,:,i)=Gmats(:,:,i) + inverse_g0k(xi*wm(i),ham_k(:,:,ik),Nlat,mu)/Lk
          enddo
       enddo
       !
       allocate(Gso(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gso=zero
       do i=1,Lmats
          Gso(:,:,:,:,:,i)=lso2nnn_reshape(Gmats(:,:,i),Nlat,Nspin,Norb)
       enddo
       !
       if(master)      call dmft_print_gf_matsubara(Gso,"G0loc",iprint=6)
       deallocate(Gso)
       !
       !real freq
       Greal=zero
       do ik=1,Lk
          do i=1,Lreal
             Greal(:,:,i)=Greal(:,:,i) + inverse_g0k(dcmplx(wr(i),eps),ham_k(:,:,ik),Nlat,mu)/Lk
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
          Aw=0d0
          do i=1,Lreal
             do iorb=1,Norb
                do ilat=1,Nlat
                   io = iorb + (ilat-1)*Norb*Nspin
                   Aw(iorb,i)=Aw(iorb,i)-aimag(Greal(io,io,i))/pi
                enddo
             enddo
             write(106,'(100F15.7)') wr(i), Aw(:,i), sum(Aw(:,i))
          enddo
          close(106)
       endif
       !
       if(master)      call dmft_print_gf_realaxis(Gso,"G0loc",iprint=6)
       deallocate(Gso)
    endif
    !
  end subroutine read_myhk



  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Put in proper order the input coming from w90. Depends on the w90 users.
  !+------------------------------------------------------------------------------------------+!
  subroutine Hk_order(Porder)
    implicit none
    integer   ,intent(out)                       :: Porder(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    integer   ,allocatable,dimension(:,:)        :: shift1,shift2,shift3
    integer                                      :: P1(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    integer                                      :: P2(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    integer                                      :: P3(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    integer                                      :: ispin,iorb,jspin,jorb,ilat,jlat
    integer                                      :: io1,jo1,io2,jo2,ndx,i,j
    !
    !-----  Ordering 1: same orbital position for each Nlat block   -----
    allocate(shift1(2,2));shift1=0
    shift1(1,:)=[1,2]
    shift1(2,:)=[7,8]
    P1=0;P1=int(eye(Nlat*Nspin*Norb))
    do i=1,size(shift1,1)
       do j=1,2
          P1(shift1(i,j),shift1(i,j))=0
       enddo
       P1(shift1(i,1),shift1(i,2))=1
       P1(shift1(i,2),shift1(i,1))=1
    enddo
    P1(1+Nlat*Norb:Nlat*Norb*Nspin,1+Nlat*Norb:Nlat*Norb*Nspin)=P1(1:Nlat*Norb,1:Nlat*Norb)
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
    allocate(shift3(6,2));shift3=0
    do i=1,Nspin*Norb
       shift3(i,:)=[Nspin*Norb+i,2*Nspin*Norb+i]
    enddo
    ndx=0
    P3=0;P3=int(eye(Nlat*Nspin*Norb))
    do i=1,size(shift3,1)
       do j=1,2
          P3(shift3(i,j),shift3(i,j))=0
       enddo
       P3(shift3(i,1),shift3(i,2))=1
       P3(shift3(i,2),shift3(i,1))=1
    enddo
    !
    !--------------------  Global reordering   -------------------------
    Porder=matmul(P1,P2)
    !Porder=matmul(P1,matmul(P2,P3))
    !
  end subroutine Hk_order



  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: solve H(k) along path in the BZ.
  !+------------------------------------------------------------------------------------------+!
  subroutine build_eigenbands(fileHR,fileband,fileHk_path,fileKpoints_path,Sreal_)
    implicit none
    character(len=*),intent(in)                     :: fileHR
    character(len=*),intent(in),optional            :: fileband,fileHk_path,fileKpoints_path
    complex(8)     ,allocatable,intent(in),optional :: Sreal_(:,:,:,:,:,:)
    integer                                         :: ispin,iorb,jspin,jorb,ilat,jlat
    integer                                         :: io1,jo1,io2,jo2,ndx,ik,ifreq
    integer                                         :: Npts,Lk,Nkpathread,Mpier
    integer                                         :: P(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    real(8)        ,allocatable,dimension(:,:)      :: kpath,kgrid,Scorr
    complex(8)     ,allocatable,dimension(:,:,:)    :: Hkpath
    complex(8)     ,allocatable                     :: Gkreal(:,:,:,:,:,:,:)
    complex(8)     ,allocatable                     :: Gkr(:,:,:,:)
    type(rgb_color),allocatable,dimension(:)        :: colors,colors_orb
    logical                                         :: IOfile
    !
    call Hk_order(P)
    !
    allocate(colors(Nspin*Norb*Nlat))
    allocate(colors_orb(Norb))
    colors_orb=[red1,green1,blue1]
    do i=1,Nspin*Nlat
       colors(1+(i-1)*Norb:Norb+(i-1)*Norb)=colors_orb
    enddo
    !
    write(LOGfile,*)"Build bulk H(k) along the path M-R-G-M-X-G-X"
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
    allocate(Hkpath(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lk));Hkpath=zero
    allocate(Gkreal(Lk,Nlat,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
    allocate(Scorr(Nlat*Nspin*Norb,Nlat*Nspin*Norb));Scorr=0d0
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
       call TB_solve_model(   fileHR,Nspin,Norb,Nlat,kpath,Nkpath                      &
                          ,   colors                                                   &
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
    do ik=1,Lk
       call dmft_gk_realaxis(Hkpath(:,:,ik),1.d0/Lk,Gkreal(ik,:,:,:,:,:,:),Sreal)
    enddo
    !
    allocate(Gkr(Lk,Nspin*Nlat*Norb,Nspin*Nlat*Norb,Lreal));Gkr=zero
    do ik=1,Lk
       do ifreq=1,Lreal
          Gkr(ik,:,:,ifreq)=matmul(P,matmul(nnn2lso_reshape(Gkreal(ik,:,:,:,:,:,ifreq),Nlat,Nspin,Norb),transpose(P)))
       enddo
    enddo
    !
    open(unit=106,file='Akw_s1.dat',status='unknown',action='write',position='rewind')
    open(unit=107,file='Akw_s2.dat',status='unknown',action='write',position='rewind')
    do ifreq=1,Lreal
       write(106,'(9000F18.12)')wr(ifreq),(10.*trace(-aimag(Gkr(ik,1:Nlat*Norb,1:Nlat*Norb,ifreq))/pi)+10.*ik,ik=1,Lk)
       write(107,'(9000F18.12)')wr(ifreq),(10.*trace(-aimag(Gkr(ik,1+Nlat*Norb:Nlat*Nspin*Norb,1+Nlat*Norb:Nlat*Nspin*Norb,ifreq))/pi)+10.*ik,ik=1,Lk)
    enddo
    close(106)
    close(107)
    write(LOGfile,*)"Im done on the path"
    !
  end subroutine build_eigenbands



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
    call inversion_test(g0k,g0k_tmp,1.e-5,Nlat)
  end function inverse_g0k



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
