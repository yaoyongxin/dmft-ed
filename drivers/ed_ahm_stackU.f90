!########################################################
!PURPOSE  :solve the attractive (A) Hubbard model
!          on a stripe geometry with periodic Hubbard moduation
!AUTHOR   :G Mazza & A Amaricci
!########################################################
program ed_stripe
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  !
  real(8),dimension(:),allocatable   :: wm,wr  
  !
  complex(8),allocatable             :: Hk(:,:,:)              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
  complex(8),allocatable             :: Hk_stripe(:,:,:,:)              ![2][Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
  complex(8),allocatable             :: Hloc(:,:,:,:,:)       ![Nlat][Nspin][Nspin][Norb][Norb]
  complex(8),allocatable             :: Hloc_(:,:,:,:,:)      ![Nindep][Nspin][Nspin][Norb][Norb]
  complex(8),allocatable             :: Hijk(:,:)              ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk]
  !
  complex(8),allocatable             :: Smats(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable             :: Sreal(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable             :: Gmats(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable             :: Greal(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable             :: Delta(:,:,:,:,:,:,:)  ![2][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  !
  complex(8),allocatable             :: Smats_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable             :: Sreal_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable             :: Gmats_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable             :: Greal_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable             :: Delta_(:,:,:,:,:,:,:) ![2][Nindep][Nspin][Nspin][Norb][Norb][Lmats]
  !
  real(8),allocatable,dimension(:,:) :: bath,bath_old
  real(8),allocatable,dimension(:,:) :: bath_,bath_old_
  integer                            :: Nb
  !
  integer                            :: Nindep,Nsymm
  !
  real(8),allocatable,dimension(:)   :: nii,dii,pii,eii,dpii
  real(8),allocatable,dimension(:)   :: nii_,dii_,pii_,eii_,dpii_
  !  
  real(8)                            :: ts,t_perp,alpha_join
  logical                            :: slab_perp
  real(8),allocatable,dimension(:,:) :: Usite,Usite_
  real(8),allocatable,dimension(:,:) :: Uij,Unodes
  

  real(8)                            :: Uamplitude
  real(8)                            :: delta_hop,u_max,u_min,t_min,t_max
  real(8)                            :: wmixing
  real(8),allocatable,dimension(:)   :: wt
  logical,allocatable,dimension(:)   :: hk_symm,hk_symm_
  integer                            :: Uperiod,Nperiod
  integer                            :: i,iloop,ik
  integer                            :: Lk,Lk_

  integer                            :: Nrow,Ncol
  integer                            :: row,col,ilat,i_ind

  integer                            :: unit,unit_lanc
  logical                            :: pbc_row,pbc_col
  logical                            :: symmetry_flag,rdmft_phsym,tb_dos
  integer                            :: symmetry_type
  !
  integer                            :: Xperiod,Yperiod
  integer                            :: N_Xperiod,N_Yperiod
  integer                            :: N_Umax,N_Umin
  integer                            :: Ni_Umax,Ni_Umin
  logical                            :: Xpbc,Ypbc
  logical                            :: print_GF,conv_weiss
  integer                            :: Nx,Ny,ix,iy,iz,ik_l,iw
  integer,dimension(:),allocatable   :: n_lanc_,n_lanc
!  real(8),dimension(:),allocatable   :: n_lanc
  integer                            :: i_lanc,delta_n

  real(8)                            :: lanc_error

  logical                            :: conv_dmft
  logical                            :: converged
  logical                            :: conv_nlanc
  !
  real(8)                            :: xr,dk,kx,ky,peso
  real(8),dimension(:),allocatable   :: epsik_layer,t_layer,wreal,dos
  real(8),dimension(:),allocatable   :: epsik_slab,t_slab,epsik
  !
  !--------------- LATTICE WRAP VARIABLES -----------------!
  !Lattice size:
  !PUBLIC USE: (should be set by routine)
  !=========================================================
  integer                                            :: Nlat
  !Symmetry operations
  !=========================================================
  integer,allocatable,dimension(:)                   :: indep_list
  integer,dimension(:),allocatable                   :: map_lat2ind
  integer,dimension(:,:),allocatable                 :: map_ind2lat
  !=========================================================
  integer,dimension(:),allocatable                   :: icol,irow
  integer,dimension(:,:),allocatable                 :: ij2site
  ! real(8),dimension(:,:),allocatable               :: H0
  integer                                         :: inho_type


  !MPI VARIABLES:
  integer                                       :: comm,rank,ierr
  logical                                       :: master
  real(8)                                       :: temp

  real(8),allocatable :: Etmp(:)
  complex(8) :: tmpG

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  !+------------------------+!
  !+- READ INPUT VARIABLES +-!
  !+------------------------+!
  call parse_input_variable(ts,"TS","inputRDMFT.in",default=0.5d0)
  call parse_input_variable(t_perp,"T_perp","inputRDMFT.in",default=0.5d0)
  call parse_input_variable(alpha_join,"alpha_join","inputRDMFT.in",default=1.d0)
  call parse_input_variable(wmixing,"WMIXING","inputRDMFT.in",default=0.5d0)
  call parse_input_variable(N_Umax,"N_Umax","inputRDMFT.in",default=4)
  call parse_input_variable(N_Umin,"N_Umin","inputRDMFT.in",default=4)
  call parse_input_variable(Nx,"Nx","inputRDMFT.in",default=10)
  call parse_input_variable(u_max,"u_max","inputRDMFT.in",default=0.d0)
  call parse_input_variable(u_min,"u_min","inputRDMFT.in",default=0.d0)
  call parse_input_variable(print_GF,"print_GF","inputRDMFT.in",default=.false.)
  call parse_input_variable(rdmft_phsym,"RDMFT_PHSYM","inputRDMFT.in",default=.true.)
  call parse_input_variable(inho_type,"INHO_TYPE","inputRDMFT.in",default=1)
  call parse_input_variable(temp,"TEMP","inputRDMFT.in",default=0.002d0)
  call parse_input_variable(conv_weiss,"CONV_WEISS","inputRDMFT.in",default=.true.)
  call parse_input_variable(tb_dos,"tb_dos","inputRDMFT.in",default=.false.)
  call parse_input_variable(slab_perp,"slab_perp","inputRDMFT.in",default=.false.)
  !
  call ed_read_input("inputRDMFT.in",comm)
  !call set_store_size(1024)  
  beta=1.d0/temp
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  !write(*,*) order_of_magnitude(43.2d0)

  
  lanc_error=1.d-10


  !+-----------------------------+!
  !+- BUILD LATTICE HAMILTONIAN -+!
  !+-----------------------------+!
  
  !+- 2d mesh -+!
  Lk_=(Nx/2+1)*(Nx/2+2)/2
  allocate(epsik_layer(Lk_),epsik_slab(Nx))
  Lk = Lk_*Nx
  allocate(epsik(Lk),wt(Lk))
  !+- linear disperison -+!
  do ix=1,Nx     
     kx=dble(ix)/dble(Nx)*2.d0*pi-pi
     epsik_slab(ix) = -2.d0*cos(kx)
  end do
  !
  ik=0
  ik_l=0
  do ix=0,Nx/2
     do iy=0,ix
        !+- Layer density of states -+!
        ik_l=ik_l+1
        kx=dble(ix)/dble(Nx)*2.d0*pi-pi
        ky=dble(iy)/dble(Nx)*2.d0*pi-pi
        epsik_layer(ik_l) = -2*( dcos(kx) + dcos(ky) )
        if(ix==0) then
           peso=1.d0
        elseif(ix==Nx/2) then
           if(iy==0) then
              peso=2.d0
           elseif(iy==Nx/2) then
              peso=1.d0
           else
              peso=4.d0
           end if
        else
           if(iy==ix) then
              peso=4.d0
           elseif(iy==0) then
              peso=4.d0
           else
              peso=8.d0
           end if
        end if
        !+- slab density of states -+!
        do iz=1,Nx
           ik=ik+1
           wt(ik) = peso/dble(Nx**2)/Nx
           epsik(ik) = epsik_layer(ik_l) + epsik_slab(iz)
        end do
     end do
  end do

  !
  !+-> hoppings <-+!
  !
  Xperiod=N_Umax+N_Umin
  Nlat=Xperiod  
  allocate(t_layer(Nlat),t_slab(Nlat))
  t_layer=ts
  !
  t_slab=t_perp
  t_slab(1:N_Umin)=t_perp*alpha_join
  t_slab(Nlat)=t_perp*alpha_join
  !
  
  !
  allocate(Hk(Nlat,Nlat,Lk))
  Hk=0.d0
  ik=0
  ik_l=0
  do ix=0,Nx/2
     do iy=0,ix
        ik_l=ik_l+1
        do iz=1,Nx
           ik=ik+1           
           !
           do ilat=1,Nlat
              Hk(ilat,ilat,ik) = Hk(ilat,ilat,ik)+epsik_layer(ik_l)*t_layer(ilat)
              if(ilat.lt.Nlat) then
                 Hk(ilat,ilat+1,ik) = Hk(ilat,ilat+1,ik)-t_slab(ilat)
                 Hk(ilat+1,ilat,ik) = Hk(ilat+1,ilat,ik)-t_slab(ilat)
              end if
           end do
           !
           kx=dble(iz)/dble(Nx)*2.d0*pi-pi
           Hk(1,Nlat,ik) = Hk(1,Nlat,ik) -t_slab(Nlat)*exp(+xi*kx)
           Hk(Nlat,1,ik) = Hk(Nlat,1,ik) -t_slab(Nlat)*exp(-xi*kx)
           !
        end do
     end do
  end do
  if(master) call splot("thop_layer.data",(/(dble(i),i=1,Nlat)/),t_layer)
  if(master) call splot("thop_slab.data",(/(dble(i),i=1,Nlat)/),t_slab)
  !
  call MPI_Barrier(Comm,ierr)
  if(tb_dos) then
     allocate(wreal(Lreal),dos(Lreal),Hijk(Nlat,Nlat),Etmp(Nlat))
     wreal = linspace(wini,wfin,Lreal)
     
     dos=0.d0
     do ik=1,Lk
        Hijk=Hk(:,:,ik)
        Etmp=eigvals(Hijk)
        do ilat=1,Nlat
           do iw=1,Lreal
              tmpG=1.d0/(wreal(iw)+xi*eps-Etmp(ilat))
              dos(iw)=dos(iw) - 1.d0/pi*dimag(tmpG)*wt(ik)/dble(Nlat)
           end do
        end do
     end do
     if(master) call splot("bare_dos.data",wreal,dos)
  end if

  call MPI_Barrier(Comm,ierr)  
  !if(tb_dos) stop
  ! 
  !

  !
  Nlat=size(Hk,1)
  allocate(Hk_stripe(2,Nlat,Nlat,Lk))
  Hk_stripe=zero
  Hk_stripe(1,:,:,:) = Hk
  Hk_stripe(2,:,:,:) = -Hk
  !

  !+-------------------------+!
  !+- BUILD LATTICE DETAILS -+!
  !+-------------------------+!
  allocate(Usite(Nlat,Norb),Uij(Xperiod,Yperiod))
  Usite(1:N_Umin,:) = U_min
  Usite(N_Umin+1:Nlat,:) = U_max
  if(master) call splot("Usite.data",(/(dble(i),i=1,Nlat)/),Usite(:,1))
  


  !+- independent sites -+!
  Nindep = 0
  if(mod(N_Umin,2).eq.0) then
     Ni_Umin=N_Umin/2
  else
     Ni_Umin=(N_Umin+1)/2
  end if
  !
  if(mod(N_Umax,2).eq.0) then
     Ni_Umax=N_Umax/2
  else
     Ni_Umax=(N_Umax+1)/2
  end if
  Nindep=Ni_Umax+Ni_Umin
  !maps independnt sites!
  allocate(indep_list(Nindep),map_lat2ind(Nlat))
  do i_ind=1,Ni_Umin
     indep_list(i_ind) = i_ind
  end do
  do i_ind=Ni_Umin+1,Nindep
     indep_list(i_ind) = i_ind+(N_Umin-Ni_Umin)
  end do
  if(master) call splot("indep_sites.data",(/(dble(i),i=1,Nindep)/),dble(indep_list))  
  map_lat2ind=0
  do i_ind=1,Ni_Umin
     ilat=indep_list(i_ind)
     map_lat2ind(ilat)=i_ind
     ilat=N_Umin+1-indep_list(i_ind)
     map_lat2ind(ilat)=i_ind
  end do
  do i_ind=Ni_Umin+1,Nindep
     ilat=indep_list(i_ind)
     map_lat2ind(ilat)=i_ind
     ilat=Nlat+1-i_ind+Ni_Umin
     map_lat2ind(ilat)=i_ind
  end do
  if(master) call splot("map_lat2ind.data",(/(dble(i),i=1,Nlat)/),dble(map_lat2ind))

  !


  !



  !+----------------------------------+!  
  !+- ALLOCATE GF & INITIALIZE BATHS -+!
  !+----------------------------------+!  
  ! Matsubara and Real freq
  allocate(wm(Lmats),wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  wm(:)  = pi/beta*real(2*arange(1,Lmats)-1,8)
  
  ! Local Hamiltonian
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb));Hloc=0.d0
  ! Observables
  allocate(nii(Nlat))
  allocate(dii(Nlat))
  allocate(pii(Nlat))
  allocate(eii(Nlat))
  allocate(dpii(Nlat))  

  ! Self energies
  allocate(Smats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  ! Green function
  allocate(Gmats(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(2,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  ! Impurity-bath hybritizations
  allocate(Delta(2,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  !
  Nb=get_bath_dimension()
  allocate(bath(Nlat,Nb))
  allocate(bath_old(Nlat,Nb))
  !
  !



  ! Independent sites
  allocate(Hloc_(Nindep,Nspin,Nspin,Norb,Norb))
  ! Observables
  allocate(nii_(Nindep))
  allocate(dii_(Nindep))
  allocate(pii_(Nindep))
  allocate(eii_(Nindep))
  allocate(dpii_(Nindep))

  allocate(n_lanc_(Nindep))  
  allocate(n_lanc(Nindep))  
  ! Self energies
  allocate(Smats_(2,Nindep,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal_(2,Nindep,Nspin,Nspin,Norb,Norb,Lreal))
  ! Green function
  allocate(Gmats_(2,Nindep,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal_(2,Nindep,Nspin,Nspin,Norb,Norb,Lreal))
  ! Impurity-bath hybritizations
  allocate(Delta_(2,Nindep,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Usite_(Nindep,Norb))
  !
  Nb=get_bath_dimension()
  allocate(bath_(Nindep,Nb))
  allocate(bath_old_(Nindep,Nb))
  
  do i_ind=1,Nindep
     !bath_(i_ind,:) = bath(indep_list(i_ind),:)
     Hloc_(i_ind,:,:,:,:) = Hloc(indep_list(i_ind),:,:,:,:)
     Usite_(i_ind,:) = Usite(indep_list(i_ind),:)
  end do


  call ed_init_solver(comm,bath_,Hloc_)
  call MPI_Barrier(Comm,ierr)
  !
  if(rdmft_phsym)then
     do i=1,Nindep
        call ph_symmetrize_bath(bath_(i,:))
     enddo
  endif
  call MPI_Barrier(Comm,ierr)



  ! !+-------------+!
  ! !+- DMFT LOOP -+!
  ! !+-------------+!
  !
  !
  call ed_get_neigen_total(n_lanc,Nindep)
  !
  !
  unit=free_unit()
  if(master) open(unit,file='average_observables_loop.dat',position="append")  
  unit_lanc=free_unit()
  if(master) open(unit_lanc,file='nlanc_error.err',position="append")  

  iloop=0 ; converged=.false.; conv_nlanc=.false.; conv_dmft=.false.
  i_lanc=0
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     !
     if(master) call start_loop(iloop,nloop,"DMFT-loop")
     call MPI_Barrier(Comm,ierr)     
     !
     bath_old_=bath_
     !
     call ed_solve(comm,bath_,Hloc_,Uloc_ii=Usite_)
     !
     call ed_get_dens(nii_,Nindep,1)
     call ed_get_docc(dii_,Nindep,1)
     call ed_get_epot(eii_,Nindep)
     call ed_get_phisc(pii_,Nindep,1)     
     call ed_get_neigen_total(n_lanc_,Nindep)
     forall(i=1:Nindep) dpii_(i)=dii_(i)**2.0+pii_(i)**2.d0

     call ed_get_sigma_matsubara(Smats_(1,:,:,:,:,:,:),Nindep)
     call ed_get_sigma_real(Sreal_(1,:,:,:,:,:,:),Nindep)     !
     call ed_get_self_matsubara(Smats_(2,:,:,:,:,:,:),Nindep)
     call ed_get_self_real(Sreal_(2,:,:,:,:,:,:),Nindep)
     !
     do ilat=1,Nlat
        i_ind=map_lat2ind(ilat)
        Smats(:,ilat,:,:,:,:,:) = Smats_(:,i_ind,:,:,:,:,:)
        Sreal(:,ilat,:,:,:,:,:) = Sreal_(:,i_ind,:,:,:,:,:)
        nii(ilat) = nii_(i_ind)
        dii(ilat) = dii_(i_ind)
        pii(ilat) = pii_(i_ind)
        eii(ilat) = eii_(i_ind)
        dpii(ilat) = dpii_(i_ind)
     end do
     if(master) write(unit,*) dble(iloop),sum(nii)/dble(Nlat),sum(dii)/dble(Nlat),sum(pii)/dble(Nlat),sum(eii)/dble(Nlat),sum(dpii)/dble(Nlat) 
     !
     call dmft_gloc_matsubara(Comm,Hk_stripe,wt,Gmats,Smats,hk_symm=hk_symm)
     call dmft_gloc_realaxis(Comm,Hk_stripe,wt,Greal,Sreal,hk_symm=hk_symm)
     !
     do i_ind=1,Nindep
        ilat=indep_list(i_ind)
        Gmats_(:,i_ind,:,:,:,:,:)=Gmats(:,ilat,:,:,:,:,:)
        Greal_(:,i_ind,:,:,:,:,:)=Greal(:,ilat,:,:,:,:,:)
     end do
     !
     if(cg_scheme=='weiss') then
        call dmft_weiss(Comm,Gmats_(1,:,:,:,:,:,:),Gmats_(2,:,:,:,:,:,:),&
             Smats_(1,:,:,:,:,:,:),Smats_(2,:,:,:,:,:,:),Delta_,Hloc_)
     else
        call dmft_delta(Comm,Gmats_(1,:,:,:,:,:,:),Gmats_(2,:,:,:,:,:,:),&
             Smats_(1,:,:,:,:,:,:),Smats_(2,:,:,:,:,:,:),Delta_,Hloc_)
     end if
     call MPI_Barrier(Comm,ierr)
     !  
     call ed_chi2_fitgf(Comm,bath_,Delta_,Hloc_)
     call MPI_Barrier(Comm,ierr)
     !  
     bath_=wmixing*bath_ + (1.d0-wmixing)*bath_old_
     if(rdmft_phsym)then
        do i=1,Nindep
           call ph_symmetrize_bath(bath_(i,:))
        enddo
     endif
     call MPI_Barrier(Comm,ierr)
     !
     !
     !
     
     if(master) then
        
        delta_n=0
        do i=1,Nindep
           delta_n = delta_n+abs(n_lanc(i)-n_lanc_(i))
        end do
        write(unit_lanc,*) iloop,delta_n
        if(delta_n.eq.0) then
           i_lanc=i_lanc+1
        else
           i_lanc=0
        end if
        if(i_lanc.ge.2) then
           conv_nlanc=.true.
        else
           conv_nlanc=.false.
        end if
        
        !conv_nlanc = check_convergence_global(n_lanc,lanc_error,Nsuccess_lanc,nloop,id=0,file="lanc_error.err")
        if(conv_weiss) then
           conv_dmft = check_convergence(Delta(1,:,1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
        else
           conv_dmft = check_convergence_global(dpii,dmft_error,Nsuccess,nloop,reset=.false.,id=0,file="error.err")
        end if
        converged=conv_nlanc.AND.conv_dmft
     end if
     call bcast_MPI(comm,converged)

     n_lanc=n_lanc_


     !    !+-> print 
     if(master) call print_sc_out(converged)
     if(master) call end_loop
     !    !
  enddo
  converged=.true.
  if(master) call print_sc_out(converged)
  if(master) close(unit)
  if(master) close(unit_lanc)

  call stripe_energy
  call finalize_MPI
  
  !+******************************************************************+!
  !+******************************************************************+!
  !+******************************************************************+!
CONTAINS
  !+----------------------+!
  !+- AUXILIARY ROUTINES -+!
  !+----------------------+!
  subroutine print_sc_out(converged)
    implicit none
    integer                              :: i,j,is,row,col,ilat,jlat
    real(8)                              :: nimp,phi,ccdw,docc
    real(8),dimension(Nlat)              :: cdwii,rii,sii,zii
    real(8),dimension(Nrow,Ncol)         :: dij,nij,cij,pij
    real(8),dimension(Nrow)              :: grid_x
    real(8),dimension(Ncol)              :: grid_y
    real(8)                              :: mean,sdev,var,skew,kurt
    real(8),dimension(2,Nlat)            :: data_covariance
    real(8),dimension(2,2)               :: covariance_nd
    real(8),dimension(2)                 :: data_mean,data_sdev
    logical                              :: converged
    complex(8),dimension(2,Lmats)        :: aGmats,aSmats
    complex(8),dimension(2,Lreal)        :: aGreal,aSreal
    character(len=4)                     :: loop
    integer,dimension(6)                 :: units,N_min
    real(8)                              :: Eint,Ekin,Epot,Eout(2)
    complex(8),dimension(1,1,Nlat,Lmats) :: Sigma_tmp
    complex(8),dimension(1,1,Nlat,Lmats) :: SigmaA_tmp

    !Get CDW "order parameter"
    do is = 1,Nlat
       ilat = mod(is-1,Nlat)+1
       jlat = (is-1)/Nlat+1
       cdwii(is) = (-1.d0)**(ilat+jlat)*(nii(is)-1.d0)
       write(*,'(3(I3),F18.10,I3)') is,ilat,jlat,cdwii(is),Nlat
    end do
    nimp = sum(nii)/dble(Nlat)
    phi  = sum(pii)/dble(Nlat)
    docc = sum(dii)/dble(Nlat)
    ccdw = sum(cdwii)/dble(Nlat)
    !+- use splot -+!
    call splot("nVSisite.data",dble(arange(1,Nlat)),nii)
    call splot("phiVSisite.data",dble(arange(1,Nlat)),pii)
    call splot("doccVSisite.data",dble(arange(1,Nlat)),dii)
    !
    ! Plots at convergence
    if(converged)then
       ! GF & Sigma
       if(print_GF) then
          call dmft_print_gf_matsubara(Gmats(1,:,:,:,:,:,:),"LG_iw.data",4)
          call dmft_print_gf_matsubara(Gmats(2,:,:,:,:,:,:),"LF_iw.data",4)
          !
          call dmft_print_gf_realaxis(Greal(1,:,:,:,:,:,:),"LG_realw.data",4)
          call dmft_print_gf_realaxis(Greal(2,:,:,:,:,:,:),"LF_realw.data",4)
          !
          call dmft_print_gf_matsubara(Smats(1,:,:,:,:,:,:),"LSigma_iw.data",4)
          call dmft_print_gf_matsubara(Smats(2,:,:,:,:,:,:),"LSelf_iw.data",4)
          !
          call dmft_print_gf_realaxis(Sreal(1,:,:,:,:,:,:),"LSigma_realw.data",4)
          call dmft_print_gf_realaxis(Sreal(2,:,:,:,:,:,:),"LSelf_realw.data",4)
       end if

       ! Observables
       do is=1,Nlat
          cdwii(is) = (-1.d0)**(is)*(nii(is)-1.d0)
          sii(is)   = dimag(Smats(1,is,1,1,1,1,1))-&
               wm(1)*(dimag(Smats(1,is,1,1,1,1,2))-dimag(Smats(1,is,1,1,1,1,1)))/(wm(2)-wm(1))
          rii(is)   = dimag(Gmats(1,is,1,1,1,1,1))-&
               wm(1)*(dimag(Gmats(1,is,1,1,1,1,2))-dimag(Gmats(1,is,1,1,1,1,1)))/(wm(2)-wm(1))
          zii(is)   = 1.d0/( 1.d0 + abs( dimag(Smats(1,is,1,1,1,1,1))/wm(1) ))
       enddo
       rii=abs(rii)
       sii=abs(sii)
       zii=abs(zii)

       ! write(*,*)"Print col, row"
       ! units = free_units(6)
       ! open(units(1),file='n_col.data')
       ! open(units(2),file='docc_col.data')
       ! open(units(3),file='phi_col.data')
       ! open(units(4),file='n_row.data')
       ! open(units(5),file='docc_row.data')
       ! open(units(6),file='phi_row.data')
       ! !
       ! do col=1,Ncol
       !    grid_y(col)=col
       !    do row=1,Nrow
       !       grid_x(row)  = row
       !       i            = ij2site(row,col)
       !       nij(row,col) = nii(i)
       !       dij(row,col) = dii(i)
       !       pij(row,col) = pii(i)
       !    enddo
       ! enddo
       ! !
       ! do row=1,Nrow
       !    write(units(1),'(100(f18.10))') dble(row),nij(row,:)
       !    write(units(2),'(100(f18.10))') dble(row),dij(row,:)
       !    write(units(3),'(100(f18.10))') dble(row),pij(row,:)
       ! end do
       ! !
       ! do col=1,Ncol
       !    write(units(4),'(100(f18.10))') dble(col),nij(:,col)
       !    write(units(5),'(100(f18.10))') dble(col),dij(:,col)
       !    write(units(6),'(100(f18.10))') dble(col),pij(:,col)
       ! end do
       ! !
       ! do i=1,6
       !    close(units(i))
       ! end do
       !
    end if
    !
    !
    write(*,*)"Exit print_sc_out"
  end subroutine print_sc_out


  subroutine stripe_energy
    integer :: ispin,jspin,iorb,jorb,is,js
    integer,dimension(6)                 :: units,N_min
    real(8)                              :: Eint,Ekin,Epot,Eout(2)
    real(8),dimension(Nlat)              :: Ekin_
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb,Lmats) :: Sigma_tmp
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb,Lmats) :: SigmaA_tmp
    !
    Eint=0.d0
    Ekin=0.d0
    Epot=0.d0
    Epot=sum(eii)/dble(Nlat)
    write(*,*)"DMFT_kinetic_energy:"
    call dmft_kinetic_energy(Comm,Hk,Wt,Smats(1,:,:,:,:,:,:),Smats(2,:,:,:,:,:,:),Ekin=Ekin_)
    write(*,*)"done..."
    !
    Ekin = sum(Ekin_)/dble(Nlat)
    Eint=Ekin+Epot
    unit=free_unit()
    if(master) then
       open(unit,file='internal_energy.data')
       write(unit,'(10(F18.10))') Eint,Ekin,Epot
       close(unit)
    end if
    !
  end subroutine stripe_energy



  !   subroutine  get_Hk_2dsquare(Nx,Xpbc,Ny,Ypbc,Lx,Ly,Hk_lat)
  !     implicit none
  !     integer                  :: Nx,Ny
  !     integer                  :: Lx,Ly
  !     logical                  :: Xpbc,Ypbc
  !     real(8),allocatable      :: kxgrid(:),kygrid(:),kxgrid_(:),kygrid_(:)    
  !     real(8),allocatable      :: epsik(:)
  !     real(8),allocatable      :: Hlat(:,:)
  !     integer                 ::ik,Nx_,Ny_,Nlat_,Nlat_sub,Lk_,Lk
  !     integer                  :: ilat,jlat,ilat_,jlat_
  !     integer                  :: ix_,iy_,ix_sub,iy_sub,ilat_sub,jlat_sub

  !     integer,allocatable      :: map_lattice_new2old_(:,:,:)
  !     integer,allocatable      :: map_lattice_new2old(:,:)
  !     integer,allocatable      :: map_lattice_old2new_sites(:),map_lattice_old2new_sublat(:)
  !     real(8),allocatable      :: thop_lattice(:,:,:,:)
  !     type(vect2D),allocatable :: Rlat_(:),kVect_(:)
  !     complex(8),allocatable   :: Hk_lat(:,:,:)
  !     complex(8)               :: psi_ki,psi_kj
  !     !
  !     if(mod(Nx,Lx)/=0) stop "wrong X-periodicity"
  !     if(mod(Ny,Ly)/=0) stop "wrong Y-periodicity"


  !     !+- STEP-1: build the lattice hamiltonian using FT of the energy-momentum dispersion with the right boundary conditions.
  !     allocate(kxgrid(Nx),kygrid(Ny))    
  !     if(Xpbc) then
  !        kxgrid = linspace(0.d0,2.d0*pi,Nx,istart=.true.,iend=.false.)
  !     else
  !        kxgrid = linspace(0.d0,pi,Nx,istart=.false.,iend=.false.)
  !     end if
  !     if(Ypbc) then
  !        kygrid = linspace(0.d0,2.d0*pi,Ny,istart=.true.,iend=.false.)
  !     else
  !        kygrid = linspace(0.d0,pi,Ny,istart=.false.,iend=.false.)  
  !     end if
  !     Nlat=Nx*Ny; Lk = Nx*Ny; allocate(epsik(Lk))
  !     ik=0
  !     do ix=1,Nx
  !        do iy=1,Ny
  !           ik = ik + 1
  !           epsik(ik) = -2*ts*( dcos(kxgrid(ix)) + dcos(kygrid(iy)))
  !        end do
  !     end do
  !     allocate(Hlat(Nlat,Nlat))
  !     call k2latticeFT(kxgrid,Xpbc,kygrid,Ypbc,epsik,Hlat)        
  !     !+- STEP-2: build the hopping matrices for a lattice of different periodicity -+!
  !     !           1 <= Lx =< Nx  ; 1 <= Ly =< Ny            !
  !     Nx_=Nx/Lx;Ny_=Ny/Ly    
  !     Nlat_=Nx_*Ny_;Nlat_sub=Lx*Ly    
  !     allocate(Rlat_(Nlat_))
  !     allocate(map_lattice_new2old(Nlat_,Nlat_sub),map_lattice_new2old_(Nlat_,Lx,Ly))
  !     allocate(map_lattice_old2new_sites(Nlat),map_lattice_old2new_sublat(Nlat))
  !     allocate(thop_lattice(Nlat_,Nlat_,Nlat_sub,Nlat_sub))    
  !     do ix_=1,Nx_
  !        do iy_=1,Ny_
  !           ilat_=square_lattice_stride(ix_,iy_,Nx_)        
  !           Rlat_(ilat_)%x= dble(ix_)*dble(Lx)
  !           Rlat_(ilat_)%y= dble(iy_)*dble(Ly)
  !           do ix_sub=1,Lx
  !              do iy_sub=1,Ly
  !                 ilat_sub = square_lattice_stride(ix_sub,iy_sub,Lx)              
  !                 ix = (ix_-1)*Lx+ix_sub
  !                 iy = (iy_-1)*Ly+iy_sub
  !                 ilat = square_lattice_stride(ix,iy,Nx)
  !                 map_lattice_new2old_(ilat_,ix_sub,iy_sub) = ilat
  !                 map_lattice_new2old(ilat_,ilat_sub) = ilat
  !                 map_lattice_old2new_sites(ilat) = ilat_
  !                 map_lattice_old2new_sublat(ilat) = ilat_sub
  !              end do
  !           end do
  !        end do
  !     end do
  !     thop_lattice=0.d0
  !     do ilat_=1,Nlat_
  !        do jlat_=1,Nlat_
  !           do ilat_sub=1,Nlat_sub
  !              do jlat_sub=1,Nlat_sub
  !                 ilat=map_lattice_new2old(ilat_,ilat_sub)
  !                 jlat=map_lattice_new2old(jlat_,jlat_sub)
  !                 thop_lattice(ilat_,jlat_,ilat_sub,jlat_sub) = Hlat(ilat,jlat)
  !              end do
  !           end do
  !        end do
  !     end do
  !     !
  !     !+- STEP 3 -+! 
  !     !   build-up the new momentum grid corresponding to the periodicity (Lx,Ly)
  !     Lk_ =Nx_*Ny_;allocate(kxgrid_(Nx_),kygrid_(Ny_),kVect_(Lk_))
  !     if(Xpbc) then
  !        kxgrid_ = linspace(0.d0,2.d0*pi/dble(Lx),Nx_,istart=.true.,iend=.false.)
  !     else
  !        kxgrid_ = linspace(0.d0,pi/dble(Lx),Nx_,istart=.false.,iend=.false.)
  !     end if
  !     if(Ypbc) then
  !        kygrid_ = linspace(0.d0,2.d0*pi/dble(Ly),Ny_,istart=.true.,iend=.false.)
  !     else
  !        kygrid_ = linspace(0.d0,pi/dble(Ly),Ny_,istart=.false.,iend=.false.)
  !     end if
  !     ik=0
  !     do ix_=1,Nx_
  !        do iy_=1,Ny_
  !           ik=ik+1
  !           kVect_(ik)%x=kxgrid_(ix_)
  !           kVect_(ik)%y=kygrid_(iy_)
  !        end do
  !     end do
  !     !+- STEP 4 -+! 
  !     !   For each (ilat_sub,jlat_sub) transform back thop_lattice(:,:,ilat_sub,jlat_sub) to momentum space
  !     if(allocated(Hk_lat)) deallocate(Hk_lat)
  !     allocate(Hk_lat(Nlat_sub,Nlat_sub,Lk_))
  !     do ik=1,Lk_
  !        Hk_lat(:,:,ik) = 0.d0
  !        do ilat_=1,Nlat_
  !           do jlat_=1,Nlat_
  !              !
  !              psi_ki=1.d0
  !              psi_kj=1.d0             
  !              if(Xpbc) then
  !                 psi_ki=psi_ki*sqrt(1.d0/dble(Nx_))*exp(-xi*kVect_(ik)%x*Rlat_(ilat_)%x)
  !                 psi_kj=psi_kj*sqrt(1.d0/dble(Nx_))*exp(-xi*kVect_(ik)%x*Rlat_(jlat_)%x)
  !              else
  !                 psi_ki=psi_ki*sqrt(2.d0/dble(Nx_+1))*sin(kVect_(ik)%x*Rlat_(ilat_)%x)
  !                 psi_kj=psi_kj*sqrt(2.d0/dble(Nx_+1))*sin(kVect_(ik)%x*Rlat_(jlat_)%x)
  !              end if
  !              if(Ypbc) then
  !                 psi_ki=psi_ki*sqrt(1.d0/dble(Ny_))*exp(-xi*kVect_(ik)%y*Rlat_(ilat_)%y)
  !                 psi_kj=psi_kj*sqrt(1.d0/dble(Ny_))*exp(-xi*kVect_(ik)%y*Rlat_(jlat_)%y)
  !              else
  !                 psi_ki=psi_ki*sqrt(2.d0/dble(Ny_+1))*sin(kVect_(ik)%y*Rlat_(ilat_)%y)
  !                 psi_kj=psi_kj*sqrt(2.d0/dble(Ny_+1))*sin(kVect_(ik)%y*Rlat_(jlat_)%y)
  !              end if
  !              Hk_lat(:,:,ik) = Hk_lat(:,:,ik) + thop_lattice(ilat_,jlat_,:,:)*psi_kj*conjg(psi_ki)
  !           end do
  !        end do
  !     end do

  !     unit=free_unit()
  !     open(unit,file="hk_symmetric_points")
  !     Hk_symm=.true.
  !     do ik=1,Lk_
  !        do ilat_sub=1,Nlat_sub
  !           do jlat_sub=1,Nlat_sub
  !              if(Hk(ilat_sub,jlat_sub,ik) /= Hk(jlat_sub,ilat_sub,ik)) Hk_symm(ik)=.false.
  !           end do
  !        end do
  !        if(Hk_symm(ik)) write(unit,'(4(F18.10))') kVect_(ik)%x,kVect_(ik)%y
  !     end do
  !   end subroutine get_Hk_2dsquare



  !   subroutine k2latticeFT(kx_grid,Xpbc,ky_grid,Ypbc,epsik,Hlat) 
  !     real(8)                     :: kx_grid(:),ky_grid(:)
  !     logical                     :: Xpbc,Ypbc
  !     real(8)                     :: epsik(size(kx_grid)*size(ky_grid))
  !     real(8)                     :: Hlat(size(kx_grid)*size(ky_grid),size(kx_grid)*size(ky_grid))
  !     !
  !     type(vect2D),allocatable   :: kgrid(:),Rlat(:)
  !     integer                     :: Nx,Ny
  !     !    complex(8),dimension(Nx,Ny) :: Hk
  !     !type(vect2D),dimension(:),allocatable   :: 

  !     real(8) :: arg,test
  !     complex(8) :: psi_ki,psi_kj
  !     integer :: i,j,ilat,jlat,ill,ill_
  !     integer :: Lk,Nlat    
  !     !
  !     integer,allocatable :: stride_ill2ij(:)
  !     real(8),allocatable :: Hlat_tmp(:,:)
  !     !
  !     Nx=size(kx_grid);Ny=size(ky_grid);Nlat = Nx*Ny;Lk=Nlat
  !     allocate(Rlat(Nlat),kgrid(Lk))
  !     !+- build lattice vectors -+!
  !     ik=0
  !     do i=1,Nx
  !        do j=1,Ny
  !           ik=ik+1
  !           kgrid(ik)%x=kx_grid(i)
  !           kgrid(ik)%y=ky_grid(j)          
  !           ilat=(j-1)*Nx+i
  !           Rlat(ilat)%x = dble(i)
  !           Rlat(ilat)%y = dble(j)
  !        end do
  !     end do
  !     !+-------------------------+!
  !     allocate(Hlat_tmp(Nlat,Nlat))
  !     write(LOGfile,*) "Building full Real-space hamiltonian from Fourier transform."
  !     if(master)call start_timer
  !     do ill=1+mpiID,Nlat*Nlat,mpiSIZE
  !        ilat=mod(ill-1,Nlat)+1
  !        jlat=(ill-1)/Nlat+1
  !        Hlat_tmp(ilat,jlat) = 0.d0
  !        do ik=1,Lk
  !           !
  !           psi_ki=1.d0
  !           psi_kj=1.d0             
  !           if(Xpbc) then
  !              psi_ki=psi_ki*sqrt(1.d0/dble(Nx))*exp(xi*kgrid(ik)%x*Rlat(ilat)%x)
  !              psi_kj=psi_kj*sqrt(1.d0/dble(Nx))*exp(xi*kgrid(ik)%x*Rlat(jlat)%x)
  !           else
  !              psi_ki=psi_ki*sqrt(2.d0/dble(Nx+1))*sin(kgrid(ik)%x*Rlat(ilat)%x)
  !              psi_kj=psi_kj*sqrt(2.d0/dble(Nx+1))*sin(kgrid(ik)%x*Rlat(jlat)%x)
  !           end if
  !           if(Ypbc) then
  !              psi_ki=psi_ki*sqrt(1.d0/dble(Ny))*exp(xi*kgrid(ik)%y*Rlat(ilat)%y)
  !              psi_kj=psi_kj*sqrt(1.d0/dble(Ny))*exp(xi*kgrid(ik)%y*Rlat(jlat)%y)
  !           else
  !              psi_ki=psi_ki*sqrt(2.d0/dble(Ny+1))*sin(kgrid(ik)%y*Rlat(ilat)%y)
  !              psi_kj=psi_kj*sqrt(2.d0/dble(Ny+1))*sin(kgrid(ik)%y*Rlat(jlat)%y)
  !           end if
  !           !
  !           Hlat_tmp(ilat,jlat) = Hlat_tmp(ilat,jlat) + epsik(ik)*psi_ki*conjg(psi_kj)
  !        end do
  !        !if(master) write(*,*) ill,Nlat*Nlat
  !        if(master)call eta(ill,Nlat*Nlat,unit=LOGfile)
  !     end do
  !     if(mpiId==0) call stop_timer
  !     call MPI_ALLREDUCE(Hlat_tmp,Hlat,Nlat*Nlat,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,MPIerr)
  !   end subroutine k2latticeFT


  !   function square_lattice_stride(ix,iy,Ncol) result(ilat) 
  !     integer :: ix,iy,Ncol 
  !     integer :: ilat
  !     ilat = (iy-1)*Ncol + ix
  !     return
  !   end function square_lattice_stride



  !   !+-----------------------------------------------------------------------------+!
  !   !PURPOSE:  
  !   ! find all inequivalent sites with respect the user defined symmetry operations 
  !   ! Build and check maps from the full(independent) lattice to the independent 
  !   ! (full) lattice                        !
  !   !+-----------------------------------------------------------------------------+!
  !   subroutine get_independent_sites(symmetry_operations,Nsymm,Nindep)
  !     integer,intent(in)                 :: Nsymm
  !     integer,intent(inout)              :: Nindep
  !     integer                            :: i,unit,isymm
  !     integer,dimension(Nlat)            :: tmp_search
  !     integer,dimension(Nlat,Nsymm)      :: tmp_map
  !     integer                            :: i_ind,check_maps
  !     character(len=5)                   :: tmp_suffix
  !     ! integer,dimension(:),allocatable   :: map_lat2ind
  !     ! integer,dimension(:,:),allocatable :: map_ind2lat
  !     !integer,allocatable,dimension(:)   :: indep_list
  !     interface
  !        function symmetry_operations(site_in) result(sites_out)
  !          implicit none
  !          integer                  :: site_in
  !          integer,allocatable      :: sites_out(:)
  !        end function symmetry_operations
  !     end interface
  !     !+- search and get number of independent sites -+!
  !     tmp_search=0
  !     i_ind=0
  !     do i=1,Nlat
  !        tmp_map(i,:)=symmetry_operations(i)
  !        if(tmp_search(i).ge.0) then
  !           i_ind=i_ind+1
  !           tmp_search(i)=i
  !           do isymm=1,Nsymm
  !              tmp_search(tmp_map(i,isymm))=-1
  !           end do
  !        end if
  !     end do
  !     write(*,*) Nlat
  !     !
  !     Nindep=i_ind
  !     ! (remember: each site is connected with Nsymm sites (+ 1 = identity)) !
  !     allocate(indep_list(Nindep),map_lat2ind(Nlat),map_ind2lat(Nindep,Nsymm+1))  
  !     !
  !     !+- get list of independent sites -+!
  !     i_ind=0
  !     unit=free_unit()    
  !     if(master) open(unit,file='independent_sites.lattice')
  !     do i=1,Nlat
  !        if(tmp_search(i).ge.0) then
  !           i_ind=i_ind+1
  !           indep_list(i_ind) = tmp_search(i)
  !           if(master) write(unit,*) dble(icol(indep_list(i_ind))),dble(irow(indep_list(i_ind))),dble(indep_list(i_ind)) 
  !        end if
  !     end do
  !     if(master) close(unit)
  !     !+-  build maps -+!
  !     !
  !     write(*,*) "NINDEP",Nindep
  !     write(*,*) indep_list
  !     do i_ind=1,Nindep
  !        map_lat2ind(indep_list(i_ind))=i_ind
  !        do isymm=1,Nsymm
  !           map_lat2ind(tmp_map(indep_list(i_ind),isymm))=i_ind
  !        end do
  !     end do
  !     ! 
  !     do i_ind=1,Nindep
  !        unit=free_unit()
  !        !write(tmp_suffix,'(I4.4)') i_ind
  !        ed_file_suffix="_site"//reg(txtfy(i_ind,Npad=4))!trim(tmp_suffix)
  !        if(master) open(unit,file='equivalents'//trim(tmp_suffix)//'.lattice')
  !        map_ind2lat(i_ind,1) = indep_list(i_ind)
  !        if(master) write(unit,*) icol(indep_list(i_ind)),irow(indep_list(i_ind))
  !        do isymm=1,Nsymm
  !           map_ind2lat(i_ind,isymm+1) = tmp_map(indep_list(i_ind),isymm)
  !           if(master) write(unit,*) icol(tmp_map(indep_list(i_ind),isymm)),irow(tmp_map(indep_list(i_ind),isymm))
  !        end do
  !        if(master) close(unit)
  !     end do
  !     !+- check maps +-!
  !     do i_ind=1,Nindep
  !        do isymm=1,Nsymm+1
  !           check_maps=map_ind2lat(i_ind,isymm)
  !           if(i_ind /= map_lat2ind(check_maps)) stop "WRONG MAPS"
  !        end do
  !     end do
  !     ed_file_suffix=""
  !     ! allocate(Ineq_sites_list(Nindep))
  !     ! Ineq_sites_list = indep_list
  !   end subroutine get_independent_sites


  ! subroutine get_k_hamiltonian_stripe_fast(Nrow,Ncol,Lk_row,Lk_col,tb_dos_)
  !   USE DMFT_VECTORS
  !   implicit none
  !   integer               :: Nrow,Lk_row
  !   integer               :: Ncol,Lk_col
  !   integer               :: Lx,Ly
  !   logical,optional      :: tb_dos_
  !   logical               :: tb_dos
  !   integer               :: Nsquare
  !   logical               :: pbc_row,pbc_col
  !   logical               :: symm
  !   logical               :: k_connect
  !   integer               :: i,jj,j,k,row,col,link(4),Lk,ik
  !   integer               :: unit,unitk,ix,iy
  !   real(8),allocatable      :: kx_grid(:),ky_grid(:)
  !   type(vect2D),allocatable   :: BZ_kgrid(:)
  !   type(vect2D) :: RnnU,RnnD,RnnL,RnnR
  !   real(8),allocatable   :: htmp(:,:),wk(:)    
  !   real(8),dimension(:,:),allocatable  :: tL,tR,tU,tD,H0
  !   real(8) :: tmp_exp
  !   real(8),dimension(:),allocatable :: tcol,trow
  !   real(8),dimension(:),allocatable :: wreal,dos

  !   real(8),dimension(:,:),allocatable :: Hijk
  !   real(8),dimension(:),allocatable :: Etmp
  !   complex(8) :: tmpG
  !   integer :: iw
  !   !

  !   tb_dos=.false.
  !   if(present(tb_dos_)) tb_dos=tb_dos_

  !   allocate(wreal(Lreal),dos(Lreal))
  !   wreal = linspace(wini,wfin,Lreal)
  !   !
  !   Lx=Lk_col
  !   Ly=Lk_row
  !   Nlat=Ncol*Nrow
  !   Lk=Lx*Ly
  !   allocate(Hijk(Nlat,Nlat),Etmp(Nlat))


  !   !
  !   allocate(tcol(Ncol),trow(Nrow))
  !   do row=0,Nrow-1
  !      trow(row+1) = -ts
  !   end do
  !   do col=1,Ncol!-1
  !      tcol(col) = -1.d0*(ts+delta_hop*ts*dsin(2.d0*pi*dble(col)/dble(Ncol))**dble(inho_type))
  !      if(master) write(320,*) dble(col),tcol(col)
  !   end do
  !   ! call MPI_Barrier(Comm,ierr)
  !   ! stop


  !   allocate(kx_grid(Lx),ky_grid(Ly),BZ_Kgrid(Lk))
  !   do ix=1,Lx
  !      kx_grid(ix) = 2*pi/dble(Lx)*dble(ix-1)
  !   end do
  !   do iy=1,Ly
  !      ky_grid(iy) = 2*pi/dble(Ly)*dble(iy-1)
  !   end do
  !   ik=0
  !   do ix=1,Lx
  !      do iy=1,Ly
  !         !
  !         ik=ik+1
  !         BZ_Kgrid(ik)%x=kx_grid(ix)
  !         BZ_Kgrid(ik)%y=ky_grid(iy)          
  !         !
  !      end do
  !   end do
  !   !
  !   RnnU%x= 0.d0
  !   RnnU%y= 1.d0
  !   !
  !   RnnD%x= 0.d0
  !   RnnD%y=-1.d0
  !   !
  !   RnnL%x=-1.d0
  !   RnnL%y= 0.d0
  !   !
  !   RnnR%x= 1.d0
  !   RnnR%y= 0.d0
  !   !
  !   ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk] 
  !   allocate(H0(Nlat,Nlat))
  !   allocate(Hk(Nlat,Nlat,Lk),tL(Nlat,Nlat),tR(Nlat,Nlat),tU(Nlat,Nlat),tD(Nlat,Nlat))
  !   !+-----------------------------------------------------+!
  !   !+- GET real-space block hamiltonian H0 (aka H-IRRED) -+!
  !   !+-----------------------------------------------------+!
  !   H0=0.d0
  !   unit=free_unit()
  !   if(master) open(unit,file='rdmft_sites.lattice')
  !   allocate(icol(Nlat),irow(Nlat))
  !   allocate(ij2site(Nrow,Ncol))
  !   do row=0,Nrow-1
  !      do col=0,Ncol-1
  !         i=col+ 1 + row*Ncol
  !         !
  !         irow(i)=row+1
  !         icol(i)=col+1
  !         ij2site(row+1,col+1)=i
  !         !
  !         if(master) write(unit,*) dble(col+1),dble(row+1)

  !         ! right hop
  !         link(1)= i + 1     
  !         if((col+1)==Ncol) then
  !            link(1)=0  
  !         end if

  !         ! left  hop
  !         link(3)= i - 1    
  !         if((col-1)<0) then
  !            link(3)=0  
  !         end if

  !         ! up hop
  !         link(2)= i + Ncol 
  !         if((row+1)==Nrow) then
  !            link(2)=0  
  !         end if

  !         ! down  hop
  !         link(4)= i - Ncol 
  !         if((row-1)<0) then
  !            link(4)=0  
  !         end if
  !         !          
  !         do jj=1,4
  !            if(link(jj)>0) then
  !               select case(jj)
  !               case(1) !+- right hop
  !                  H0(i,link(jj))=tcol(col+1)
  !               case(2) !+- up hop
  !                  H0(i,link(jj))=trow(row+1)
  !               case(3) !+- left hop
  !                  H0(i,link(jj))=tcol(col)
  !               case(4) !+- down hop
  !                  H0(i,link(jj))=trow(row)
  !               end select
  !            end if
  !         enddo
  !         !
  !      enddo
  !   enddo

  !   if(master) then
  !      do ilat=1,Nlat
  !         write(330,'(60F10.4)') H0(ilat,:)
  !      end do
  !   end if
  !   if(master) close(unit)

  !   !+------------------------------------------------------------------------+!
  !   !+- GET real-space connection between blocks hamiltonian (aka H-CONNECT) -+!
  !   !+------------------------------------------------------------------------+
  !   Hk = 0.d0
  !   tL=0.d0;tR=0.d0;tU=0.d0;tD=0.d0
  !   !  tLeft 
  !   do row=0,Nrow-1
  !      col=0
  !      i=col+ 1 + row*Ncol
  !      link=0
  !      link(3)=Ncol+row*Ncol
  !      do jj=1,4
  !         if(link(jj)>0) then
  !            tL(i,link(jj))=-ts 
  !         end if
  !      enddo
  !   end do
  !   !+- build up tRight -+!
  !   do row=0,Nrow-1
  !      col=Ncol-1
  !      i=col+ 1 + row*Ncol
  !      link=0
  !      link(1)=1+row*Ncol  
  !      do jj=1,4
  !         if(link(jj)>0) then
  !            tR(i,link(jj))=-ts
  !         end if
  !      enddo
  !   end do
  !   !+- build up tDown -+!    
  !   do col=0,Ncol-1       
  !      row=0
  !      i = ij2site(row+1,col+1)
  !      link=0
  !      row = Nrow-1      
  !      link(4) = ij2site(row+1,col+1)
  !      do jj=1,4
  !         if(link(jj)>0) then
  !            tD(i,link(jj))=-ts
  !         end if
  !      enddo
  !   end do
  !   !+- build up tUp -+!    
  !   do col=0,Ncol-1       
  !      row=Nrow-1
  !      i = ij2site(row+1,col+1)
  !      link=0
  !      row = 0
  !      link(1) = ij2site(row+1,col+1)
  !      do jj=1,4
  !         if(link(jj)>0) then
  !            tU(i,link(jj))=-ts
  !         end if
  !      enddo
  !   end do

  !   !+- local fourier transform -+!

  !   do ik=1,Lk
  !      Hk(:,:,ik) = H0
  !      !+- UP -+!
  !      tmp_exp = BZ_Kgrid(ik).dot.RnnU
  !      Hk(:,:,ik) = HK(:,:,ik) + tU(:,:)*exp(xi*tmp_exp)
  !      !+- DOWN -+!
  !      tmp_exp = BZ_Kgrid(ik).dot.RnnD
  !      Hk(:,:,ik) = HK(:,:,ik) + tD(:,:)*exp(xi*tmp_exp)
  !      !+- LEFT -+!
  !      tmp_exp = BZ_Kgrid(ik).dot.RnnL
  !      Hk(:,:,ik) = HK(:,:,ik) + tL(:,:)*exp(xi*tmp_exp)
  !      !+- RIGHT -+!
  !      tmp_exp = BZ_Kgrid(ik).dot.RnnR
  !      Hk(:,:,ik) = HK(:,:,ik) + tR(:,:)*exp(xi*tmp_exp)
  !   end do
  !   !
  !   unit=free_unit()
  !   open(unit,file="hk_symmetric_points")
  !   Hk_symm=.true.
  !   do ik=1,Lk
  !      do i=1,Nlat
  !         do j=1,Nlat
  !            if(Hk(i,j,ik) /= Hk(j,i,ik)) Hk_symm(ik)=.false.
  !         end do
  !      end do
  !      if(Hk_symm(ik)) write(unit,'(4(F18.10))') BZ_Kgrid(ik)%x,BZ_Kgrid(ik)%y
  !   end do
  !   !

  !   !
  ! end subroutine get_k_hamiltonian_stripe_fast



  subroutine get_k_hamiltonian_stripe(Nrow,Ncol,Lk_row,Lk_col,tb_dos_)
    USE DMFT_VECTORS
    implicit none
    integer               :: Nrow,Lk_row
    integer               :: Ncol,Lk_col
    logical,optional      :: tb_dos_
    logical               :: tb_dos
    integer               :: Lx,Ly
    integer               :: Nsquare
    logical               :: pbc_row,pbc_col
    logical               :: symm
    logical               :: k_connect
    integer               :: i,jj,j,k,row,col,link(4),Lk,ik
    integer               :: unit,unitk,ix,iy
    real(8),allocatable      :: kx_grid(:),ky_grid(:)
    type(vect2D),allocatable   :: BZ_kgrid(:)
    type(vect2D) :: RnnU,RnnD,RnnL,RnnR
    real(8),allocatable   :: htmp(:,:),wk(:)    
    real(8),dimension(:,:),allocatable  :: tL,tR,tU,tD,H0
    real(8) :: tmp_exp,tmpR
    real(8),dimension(:),allocatable :: tcol,trow,tlink_col,tlink_row

    real(8),dimension(:),allocatable :: wreal,dos
    real(8),dimension(:,:),allocatable :: Hijk
    real(8),dimension(:),allocatable :: Etmp
    complex(8) :: tmpG
    integer :: iw
    real(8) :: t_max,t_min,x_half
    !

    tb_dos=.false.
    if(present(tb_dos_)) tb_dos=tb_dos_
    
    allocate(wreal(Lreal),dos(Lreal))
    wreal = linspace(wini,wfin,Lreal)

    !
    Lx=Lk_col
    Ly=Lk_row
    Nlat=Ncol*Nrow
    Lk=Lx*Ly
    !
    allocate(Hijk(Nlat,Nlat),Etmp(Nlat))
    
    u_min=u_min/ts
    u_max=u_max/ts
    !
    t_max = abs(Uloc(1)/u_min)
    t_min = abs(Uloc(1)/u_max)
    ts=(t_max+t_min)*0.5d0
    ! u_max=Uloc(1)-u_max
    ! u_min=Uloc(1)-u_min 
    if(master) write(176,*) u_min,u_max,t_max,t_min

    allocate(tcol(Ncol),trow(Nrow))
    do row=1,Nrow
       trow(row) = -ts
    end do
    if(master) open(unit,file='thop_irr.out')

    if(mod(Ncol,2).eq.0) then
       x_half=dble(Ncol)/2.d0
    else
       x_half=dble(Ncol+1)/2.d0
    end if
    
    tcol=0.d0
    do col=1,Ncol!-1
       tmpR=0.5d0*dble(col-1)+0.5d0*dble(col)
       tcol(col)=-ts
       !tcol(col) = tcol(col) - (t_max-ts)*smooth_theta(tmpR,0.d0)*smooth_theta(-tmpR,-dble(Ncol)*0.5d0)
       tcol(col) = tcol(col) - (t_max-ts)*smooth_theta(tmpR,0.d0)*smooth_theta(-tmpR,-x_half)
       tcol(col) = tcol(col) - (t_min-ts)*smooth_theta(tmpR,dble(Ncol)*0.5)*smooth_theta(-tmpR,-dble(Ncol))
       if(master) write(unit,*) tmpR,tcol(col)
    end do
    if(master) close(unit)

    allocate(tlink_col(Nrow),tlink_row(Ncol))
    do row=1,Nrow
       tlink_col(row)=tcol(Ncol) 
    end do
    if(master) open(unit,file='thop_connect.out')
    do col=1,Ncol
       tmpR=dble(col)-1.d0
       tlink_row(col)=-ts
       !tlink_row(col)=tlink_row(col)-(t_max-ts)*smooth_theta(tmpR,0.d0)*smooth_theta(-tmpR,-dble(Ncol)*0.5d0)
       tlink_row(col)=tlink_row(col)-(t_max-ts)*smooth_theta(tmpR,0.d0)*smooth_theta(-tmpR,-x_half)
       tlink_row(col)=tlink_row(col)-(t_min-ts)*smooth_theta(tmpR,dble(Ncol)*0.5)*smooth_theta(-tmpR,-dble(Ncol))
       if(master) write(unit,*) tmpR,tlink_row(col)
    end do
    if(master) close(unit)

    allocate(kx_grid(Lx),ky_grid(Ly),BZ_Kgrid(Lk))
    do ix=1,Lx
       kx_grid(ix) = 2*pi/dble(Lx)*dble(ix-1)
    end do
    do iy=1,Ly
       ky_grid(iy) = 2*pi/dble(Ly)*dble(iy-1)
    end do
    ik=0
    do ix=1,Lx
       do iy=1,Ly
          !
          ik=ik+1
          BZ_Kgrid(ik)%x=kx_grid(ix)
          BZ_Kgrid(ik)%y=ky_grid(iy)          
          !
       end do
    end do
    !
    RnnU%x= 0.d0
    RnnU%y= 1.d0
    !
    RnnD%x= 0.d0
    RnnD%y=-1.d0
    !
    RnnL%x=-1.d0
    RnnL%y= 0.d0
    !
    RnnR%x= 1.d0
    RnnR%y= 0.d0
    !
    ![Nlat*Norb*Nspin][Nlat*Norb*Nspin][Nk] 
    allocate(H0(Nlat,Nlat))
    allocate(Hk(Nlat,Nlat,Lk),tL(Nlat,Nlat),tR(Nlat,Nlat),tU(Nlat,Nlat),tD(Nlat,Nlat))
    !+-----------------------------------------------------+!
    !+- GET real-space block hamiltonian H0 (aka H-IRRED) -+!
    !+-----------------------------------------------------+!
    H0=0.d0
    unit=free_unit()
    if(master) open(unit,file='rdmft_sites.lattice')
    allocate(icol(Nlat),irow(Nlat))
    allocate(ij2site(Nrow,Ncol))
    do row=0,Nrow-1
       do col=0,Ncol-1
          i=col+ 1 + row*Ncol
          !
          irow(i)=row+1
          icol(i)=col+1
          ij2site(row+1,col+1)=i
          !
          if(master) write(unit,*) dble(col+1),dble(row+1)

          ! right hop
          link(1)= i + 1     
          if((col+1)==Ncol) then
             link(1)=0  
          end if

          ! left  hop
          link(3)= i - 1    
          if((col-1)<0) then
             link(3)=0  
          end if

          ! up hop
          link(2)= i + Ncol 
          if((row+1)==Nrow) then
             link(2)=0  
          end if

          ! down  hop
          link(4)= i - Ncol 
          if((row-1)<0) then
             link(4)=0  
          end if
          !          
          do jj=1,4
             if(link(jj)>0) then
                select case(jj)
                case(1) !+- right hop
                   H0(i,link(jj))=tcol(col+1)
                case(2) !+- up hop
                   H0(i,link(jj))=trow(row+1)
                case(3) !+- left hop
                   H0(i,link(jj))=tcol(col)
                case(4) !+- down hop
                   H0(i,link(jj))=trow(row)
                end select
             end if
          enddo
          !
       enddo
    enddo

    if(master) then
       do ilat=1,Nlat
          write(330,'(60F10.4)') H0(ilat,:)
       end do
    end if
    if(master) close(unit)

    !+------------------------------------------------------------------------+!
    !+- GET real-space connection between blocks hamiltonian (aka H-CONNECT) -+!
    !+------------------------------------------------------------------------+
    Hk = 0.d0
    tL=0.d0;tR=0.d0;tU=0.d0;tD=0.d0
    !  tLeft 
    do row=0,Nrow-1
       col=0
       i=col+ 1 + row*Ncol
       link=0
       link(3)=Ncol+row*Ncol !+-> col=Ncol-1
       do jj=1,4
          if(link(jj)>0) then
             tL(i,link(jj))=tlink_col(row+1)!-ts
          end if
       enddo
    end do
    !+- build up tRight -+!
    do row=0,Nrow-1
       col=Ncol-1
       i=col+ 1 + row*Ncol
       link=0
       link(1)=1+row*Ncol  !+-> col=0
       do jj=1,4
          if(link(jj)>0) then
             tR(i,link(jj))=tlink_col(row+1)!-ts
          end if
       enddo
    end do
    !+- build up tDown -+!    
    do col=0,Ncol-1       
       row=0
       i = ij2site(row+1,col+1)
       link=0
       row = Nrow-1      
       link(4) = ij2site(row+1,col+1) !+->row=Nrow-1
       do jj=1,4
          if(link(jj)>0) then
             tD(i,link(jj))=tlink_row(col+1)!-ts
          end if
       enddo
    end do
    !+- build up tUp -+!    
    do col=0,Ncol-1       
       row=Nrow-1
       i = ij2site(row+1,col+1)
       link=0
       row = 0
       link(1) = ij2site(row+1,col+1)  !+-> row=0
       do jj=1,4
          if(link(jj)>0) then
             tU(i,link(jj))=tlink_row(col+1)!-ts
          end if
       enddo
    end do

    !+- local fourier transform -+!
    do ik=1,Lk
       Hk(:,:,ik) = H0
       !+- UP -+!
       tmp_exp = BZ_Kgrid(ik).dot.RnnU
       Hk(:,:,ik) = HK(:,:,ik) + tU(:,:)*exp(xi*tmp_exp)
       !+- DOWN -+!
       tmp_exp = BZ_Kgrid(ik).dot.RnnD
       Hk(:,:,ik) = HK(:,:,ik) + tD(:,:)*exp(xi*tmp_exp)
       !+- LEFT -+!
       tmp_exp = BZ_Kgrid(ik).dot.RnnL
       Hk(:,:,ik) = HK(:,:,ik) + tL(:,:)*exp(xi*tmp_exp)
       !+- RIGHT -+!
       tmp_exp = BZ_Kgrid(ik).dot.RnnR
       Hk(:,:,ik) = HK(:,:,ik) + tR(:,:)*exp(xi*tmp_exp)
    end do
    !
    unit=free_unit()
    if(master) open(unit,file="hk_symmetric_points")
    Hk_symm=.true.
    do ik=1,Lk
       do i=1,Nlat
          do j=1,Nlat
             if(Hk(i,j,ik) /= Hk(j,i,ik)) Hk_symm(ik)=.false.
          end do
       end do
       if(master) then
          if(Hk_symm(ik)) write(unit,'(4(F18.10))') BZ_Kgrid(ik)%x,BZ_Kgrid(ik)%y
       end if
    end do
    if(master) close(unit)
    !


    if(tb_dos) then
       dos=0.d0
       do ik=1,Lk
          Hijk=Hk(:,:,ik)
          Etmp=eigvals(Hijk)
          do ilat=1,Nlat
             do iw=1,Lreal
                tmpG=1.d0/(wreal(iw)+xi*eps-Etmp(ilat))
                dos(iw)=dos(iw) - 1.d0/pi*dimag(tmpG)*wt(ik)/dble(Nlat)
             end do
          end do
       end do
       if(master) call splot("bare_dos.data",wreal,dos)
    end if

  end subroutine get_k_hamiltonian_stripe


  function smooth_theta(x,x0) result(f)
    real(8) :: x,x0,f
    ! real(8) :: xr    
    ! xr=0.5
    if(x.lt.x0) then
       f=0.d0
    else
       if(x.gt.x0+xr) then
          f=1.00
       else
          f=(1.d0 - 1.5d0*cos(pi*(x-x0)/xr) + 0.5d0*(cos(pi*(x-x0)/xr))**3)*0.5d0
       end if
    end if
  end function smooth_theta




end program ed_stripe




