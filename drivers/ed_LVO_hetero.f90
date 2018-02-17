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
  integer                                        :: iloop,i,j
  integer                                        :: Nlat,ilat
  integer                                        :: io,jo
  integer                                        :: iorb,jorb,ispin,jspin
  logical                                        :: converged
  real(8)                                        :: wmixing
  character(len=60)                              :: finput
  character(len=32)                              :: hkfile
  !Mpi:
  integer                                        :: comm,rank,ier
  logical                                        :: master
  !Bath:
  integer                                        :: Nb,unit
  real(8),allocatable                            :: Bath(:,:),Bath_red(:,:)
  !Local functions:
  complex(8),allocatable,dimension(:,:,:,:,:,:)  :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:)  :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:,:)  :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:)  :: Greal
  !Weiss&Hybridization functions
  complex(8),allocatable,dimension(:,:,:,:,:,:)  :: Weiss,Weiss_old,Weiss_red
  complex(8),allocatable,dimension(:,:,:,:,:,:)  :: Delta,Delta_old,Delta_red
  !Hmiltonian input:
  integer                                        :: Nk,Nkpath
  real(8)   ,allocatable,dimension(:)            :: Wtk
  complex(8),allocatable,dimension(:,:,:)        :: Hk
  complex(8),allocatable,dimension(:,:)          :: Hloc_nso
  complex(8),allocatable,dimension(:,:,:,:,:)    :: Hloc_nnn,Hloc_red_nnn
  !custom variables for convergence test:
  complex(8),allocatable,dimension(:)            :: conv_funct
  !custom variables for chempot search:
  character(len=32)                              :: ed_file_suffix
  logical                                        :: converged_n,upprshft
  integer                                        :: conv_n_loop=0,Nlat_max
  real(8)   ,allocatable,dimension(:)            :: bottom,top,shift
  real(8)                                        :: Alvl=0.d0
  real(8)                                        :: dw,sumdens,xmu_old
  real(8),allocatable,dimension(:)               :: wr,wm
  real(8),allocatable,dimension(:,:)             :: orb_dens
  logical                                        :: look4n=.true.
  !custom variables misc:
  logical                                        :: computeG0loc
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
  call parse_cmd_variable(finput,           "FINPUT",              default='inputED_LVO_hetero.in')
  call parse_input_variable(hkfile,         "HKFILE",finput,       default="Hk.dat")
  call parse_input_variable(nk,             "NK",finput,           default=10)
  call parse_input_variable(NLAT,           "NLAT",finput,         default=4)
  call parse_input_variable(nkpath,         "NKPATH",finput,       default=20)
  call parse_input_variable(wmixing,        "WMIXING",finput,      default=0.5d0)
  call parse_input_variable(computeG0loc,   "COMPUTEG0loc",finput, default=.false.)
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
  !#########       ALLOCATION       #########
  !
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));            Smats=zero
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));            Gmats=zero
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));            Sreal=zero
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));            Greal=zero
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats));            Weiss=zero
  allocate(Delta(Nlat,Nspin,Nspin,Norb,Norb,Lmats));            Delta=zero
  allocate(weiss_old(Nlat,Nspin,Nspin,Norb,Norb,Lmats));        weiss_old=zero
  allocate(delta_old(Nlat,Nspin,Nspin,Norb,Norb,Lmats));        delta_old=zero
  allocate(weiss_red(2,Nspin,Nspin,Norb,Norb,Lmats));           weiss_red=zero
  allocate(delta_red(2,Nspin,Nspin,Norb,Norb,Lmats));           delta_red=zero
  !
  allocate(Hk(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Nk*Nk*Nk));       Hk=zero
  allocate(Hloc_nso(Nlat*Nspin*Norb,Nlat*Nspin*Norb));          Hloc_nso=zero
  allocate(Hloc_nnn(Nlat,Nspin,Nspin,Norb,Norb));               Hloc_nnn=zero
  allocate(Hloc_red_nnn(2,Nspin,Nspin,Norb,Norb));              Hloc_red_nnn=zero
  !
  allocate(conv_funct(Lmats));                                  conv_funct=zero
  allocate(wr(Lreal));wr=0.0d0;wr=linspace(wini,wfin,Lreal,mesh=dw)
  allocate(wm(Lmats));wm=0.0d0;wm = pi/beta*real(2*arange(1,Lmats)-1,8)
  allocate(Wtk(Nk*Nk*Nk));Wtk=1.d0/(Nk*Nk*Nk)
  !
  Lfit=min(int((Uloc(1)+3.)*(beta/pi))+100,Lmats)
  if(master)write(LOGfile,'(a12,I6,2(a12,F10.3))')"Lfit:",Lfit,"iwmax:",(pi/beta)*(2*Lfit-1),"U+2D:",Uloc(1)+3.
  !
  !#########        BUILD Hk        #########
  !
  call read_myhk("LVO_hr.dat","Hk.dat","Hloc.dat","Kpoints.dat")
  Hloc_red_nnn=Hloc_nnn(1:2,:,:,:,:)
  !call ed_read_impSigma_lattice(Nlat)
  !call ed_get_Sreal(Sreal,Nlat)
  !if(master)call build_eigenbands("LVO_hr.dat","Bands.dat","Hk_path.dat","Kpoints_path.dat",Sreal)
  !stop
  !
  !#########          BATH          #########
  !
  if (bath_type/="replica") then
     Nb=get_bath_dimension()
  else
     Nb=get_bath_dimension(Hloc_nnn(1,:,:,:,:))
  endif
  if(master)write(LOGfile,*)"   Bath_size:",Nb
  allocate(Bath(Nlat,Nb));      Bath=0.0d0
  allocate(Bath_red(2,Nb));     Bath_red=0.0d0
  !
  !#########      INIT SOLVER       #########
  !
  if (ed_para)then
     call ed_init_solver(Comm,Bath_red,Hloc_red_nnn)
  else
     call ed_init_solver(Comm,Bath,Hloc_nnn)
  endif
  !
  !#########          DMFT          #########
  !
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !------  solve impurity  ------
     if (ed_para)then
        call ed_solve(comm,Bath_red,Hloc_red_nnn)
     else
        call ed_solve(comm,Bath,Hloc_nnn)
     endif
     !
     !------    get sigmas    ------
     if (ed_para)then
        call ed_get_sigma_matsubara(Smats(1:2,:,:,:,:,:))  ;Smats(3:4,:,:,:,:,:)=Smats(1:2,:,:,:,:,:)
        call ed_get_sigma_real(Sreal(1:2,:,:,:,:,:))       ;Sreal(3:4,:,:,:,:,:)=Sreal(1:2,:,:,:,:,:)
     else
        call ed_get_sigma_matsubara(Smats)
        call ed_get_sigma_real(Sreal)
     endif

     !
     !------  get local Gf's  ------
     call dmft_gloc_matsubara(Comm,Hk,Wtk,Gmats,Smats,mpi_split='k')
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=6)
     !
     !-  operations on Weiss/Delta -
     if(cg_scheme=='weiss')then
        !get Weiss
        call dmft_weiss(Gmats,Smats,Weiss,Hloc_nnn)
        call dmft_print_gf_matsubara(Weiss,"WeissG0",iprint=6)
        !mix Weiss
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_old
        !old Weiss
        Weiss_old=Weiss
        !fit Weiss
        if(ed_para) then
           Weiss_red=Weiss(1:2,:,:,:,:,:)
           if (bath_type/="replica") then
              if(ed_para)then
                 call ed_chi2_fitgf(Comm,bath_red,Weiss_red,Hloc_red_nnn,ispin=1)
                 call spin_symmetrize_bath(bath_red,.true.)
              else
                 call ed_chi2_fitgf(Comm,bath_red,Weiss_red,Hloc_red_nnn)
              endif
           else
              call ed_chi2_fitgf(Comm,bath_red,Weiss_red,Hloc_red_nnn)
           endif
           Weiss(1:2,:,:,:,:,:)=Weiss_red;Weiss(3:4,:,:,:,:,:)=Weiss_red
        else
           if (bath_type/="replica") then
              if(ed_para)then
                 call ed_chi2_fitgf(Comm,bath,Weiss,Hloc_nnn,ispin=1)
                 call spin_symmetrize_bath(bath,.true.)
              else
                 call ed_chi2_fitgf(Comm,bath,Weiss,Hloc_nnn)
              endif
           else
              call ed_chi2_fitgf(Comm,bath,Weiss,Hloc_nnn)
           endif
        endif
     elseif(cg_scheme=='delta')then
        !get Delta
        call dmft_delta(Gmats,Smats,Delta,Hloc_nnn)
        call dmft_print_gf_matsubara(Delta,"Delta",iprint=6)
        !mix Delta
        if(iloop>1)Delta = wmixing*Delta + (1.d0-wmixing)*Delta_old
        !old Delta
        Delta_old=Delta
        !fit Delta
        if(ed_para) then
           Delta_red=Delta(1:2,:,:,:,:,:)
           if (bath_type/="replica") then
              if(ed_para)then
                 call ed_chi2_fitgf(Comm,bath_red,Delta_red,Hloc_red_nnn,ispin=1)
                 call spin_symmetrize_bath(bath_red,.true.)
              else
                 call ed_chi2_fitgf(Comm,bath_red,Delta_red,Hloc_red_nnn)
              endif
           else
              call ed_chi2_fitgf(Comm,bath_red,Delta_red,Hloc_red_nnn)
           endif
           Delta(1:2,:,:,:,:,:)=Delta_red;Delta(3:4,:,:,:,:,:)=Delta_red
        else
           if (bath_type/="replica") then
              if(ed_para)then
                 call ed_chi2_fitgf(Comm,bath,Delta,Hloc_nnn,ispin=1)
                 call spin_symmetrize_bath(bath,.true.)
              else
                 call ed_chi2_fitgf(Comm,bath,Delta,Hloc_nnn)
              endif
           else
              call ed_chi2_fitgf(Comm,bath,Delta,Hloc_nnn)
           endif
        endif
     endif
     !
     !each loop operations
     if(master)then
        !
        !e - chemical potential find
        converged_n=.true.
        xmu_old=xmu
        if(ed_para)then
           Nlat_max=2
        else
           Nlat_max=4
        endif
        allocate(orb_dens(Nlat_max,Norb));orb_dens=0.d0
        call ed_get_dens(orb_dens,Nlat_max)
        sumdens=0d0
        do ilat=1,Nlat_max
           sumdens=sumdens+sum(orb_dens(ilat,:))/float(Nlat_max)
           write(LOGfile,*)"  Nlat:",ilat,orb_dens(1,:)
        enddo
        write(LOGfile,*)"  n avrg:",sumdens
        deallocate(orb_dens)
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
        if(cg_scheme=='weiss')then
           do i=1,Lmats
              conv_funct(i)=sum(nnn2lso_reshape(weiss(:,:,:,:,:,i),Nlat,Nspin,Norb))
           enddo
        else
           do i=1,Lmats
              conv_funct(i)=sum(nnn2lso_reshape(delta(:,:,:,:,:,i),Nlat,Nspin,Norb))
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
  !#########     BUILD Gloc(wr)     #########
  !
  call dmft_gloc_realaxis(Comm,Hk,Wtk,Greal,Sreal,mpi_split='k')
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=6)
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
    Porder=matmul(P1,matmul(P2,P3))
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
    open(unit=106,file='Akw.dat',status='unknown',action='write',position='rewind')
    do ifreq=1,Lreal
       write(106,'(9000F18.12)')wr(ifreq),(2.*trace(-aimag(nnn2lso_reshape(Gkreal(ik,:,:,:,:,:,ifreq),Nlat,Nspin,Norb))/pi)+10.*ik,ik=1,Lk)
    enddo
    close(106)
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
