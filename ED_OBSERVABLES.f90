!########################################################################
!PURPOSE  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_OBSERVABLES
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_BATH
  USE ED_AUX_FUNX
  !
  implicit none
  private
  !
  public :: observables_impurity
  public :: local_energy_impurity


  logical,save                          :: iolegend=.true.
  real(8),dimension(:),allocatable      :: dens,dens_up,dens_dw
  real(8),dimension(:),allocatable      :: docc
  real(8),dimension(:),allocatable      :: magZ,magX,magY
  real(8),dimension(:),allocatable      :: phisc
  real(8),dimension(:,:),allocatable    :: sz2,n2
  real(8),dimension(:,:),allocatable    :: exct_s0
  real(8),dimension(:,:),allocatable    :: exct_tz
  complex(8),dimension(:,:),allocatable :: exct_tx
  complex(8),dimension(:,:),allocatable :: exct_ty
  real(8),dimensioN(:,:),allocatable    :: zimp,simp
  real(8)                               :: s2tot
  real(8)                               :: Egs
  !


contains 



  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine observables_impurity()
    integer,dimension(Nlevels)      :: ib
    integer                         :: i,j
    integer                         :: istate
    integer                         :: isector,jsector
    integer                         :: idim,jdim
    integer                         :: isz,jsz
    integer                         :: iorb,jorb,ispin,jspin,isite,jsite,ibath
    integer                         :: numstates
    integer                         :: r,m,k
    real(8)                         :: sgn,sgn1,sgn2
    real(8)                         :: gs_weight
    real(8)                         :: Ei
    real(8)                         :: peso
    real(8)                         :: norm
    real(8),dimension(Norb)         :: nup,ndw,Sz,nt
    complex(8),dimension(:),pointer :: gscvec
    type(sector_map)                :: H,HJ
    complex(8),allocatable          :: vvinit(:)
    real(8),dimension(Norb,Norb)    :: theta_upup,theta_dwdw,theta_updw,theta_dwup
    !
    !LOCAL OBSERVABLES:
    ! density, 
    ! double occupancy, 
    ! magnetization, 
    ! orbital//spin correlations  
    ! superconducting order parameter, etc..
    allocate(dens(Norb),dens_up(Norb),dens_dw(Norb))
    allocate(docc(Norb))
    allocate(phisc(Norb),magX(Norb),magY(Norb))
    allocate(exct_S0(Norb,Norb),exct_Tz(Norb,Norb))
    allocate(exct_Tx(Norb,Norb),exct_Ty(Norb,Norb))
    allocate(magZ(Norb),sz2(Norb,Norb),n2(Norb,Norb))
    allocate(simp(Norb,Nspin),zimp(Norb,Nspin))
    !
    Egs     = state_list%emin
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    phisc   = 0.d0    
    magz    = 0.d0
    magx    = 0.d0
    magy    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
    exct_s0 = 0d0
    exct_tz = 0d0
    exct_tx = zero
    exct_ty = zero
    theta_upup=0d0
    theta_updw=0d0
    theta_dwup=0d0
    theta_dwdw=0d0
    !
    numstates=state_list%size
    do istate=1,numstates
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       !
#ifdef _MPI
       if(MpiStatus)then
          gscvec => es_return_cvector(MpiComm,state_list,istate)
       else
          gscvec => es_return_cvector(state_list,istate)
       endif
#else
       gscvec => es_return_cvector(state_list,istate)
#endif
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       idim    = getdim(isector)
       !
       if(Mpimaster)then
          call build_sector(isector,H)
          !
          do i=1,idim
             m=H%map(i)
             ib = bdecomp(m,2*Ns)
             !
             gs_weight=peso*abs(gscvec(i))**2
             !
             !Get operators:
             do iorb=1,Norb
                nup(iorb)= dble(ib(iorb))
                ndw(iorb)= dble(ib(iorb+Ns))
                sz(iorb) = (nup(iorb) - ndw(iorb))/2.d0
                nt(iorb) =  nup(iorb) + ndw(iorb)
             enddo
             !
             !Evaluate averages of observables:
             do iorb=1,Norb
                dens(iorb)     = dens(iorb)      +  nt(iorb)*gs_weight
                dens_up(iorb)  = dens_up(iorb)   +  nup(iorb)*gs_weight
                dens_dw(iorb)  = dens_dw(iorb)   +  ndw(iorb)*gs_weight
                docc(iorb)     = docc(iorb)      +  nup(iorb)*ndw(iorb)*gs_weight
                magz(iorb)     = magz(iorb)      +  (nup(iorb)-ndw(iorb))*gs_weight
                sz2(iorb,iorb) = sz2(iorb,iorb)  +  (sz(iorb)*sz(iorb))*gs_weight
                n2(iorb,iorb)  = n2(iorb,iorb)   +  (nt(iorb)*nt(iorb))*gs_weight
                do jorb=iorb+1,Norb
                   sz2(iorb,jorb) = sz2(iorb,jorb)  +  (sz(iorb)*sz(jorb))*gs_weight
                   sz2(jorb,iorb) = sz2(jorb,iorb)  +  (sz(jorb)*sz(iorb))*gs_weight
                   n2(iorb,jorb)  = n2(iorb,jorb)   +  (nt(iorb)*nt(jorb))*gs_weight
                   n2(jorb,iorb)  = n2(jorb,iorb)   +  (nt(jorb)*nt(iorb))*gs_weight
                enddo
             enddo
             s2tot = s2tot  + (sum(sz))**2*gs_weight
          enddo
          if(associated(gscvec))nullify(gscvec)
          call delete_sector(isector,H)
       endif
    enddo
    !

    !SUPERCONDUCTING ORDER PARAMETER
    if(ed_mode=="superc")then
       do ispin=1,Nspin
          do iorb=1,Norb
             do istate=1,state_list%size
                !
                isector = es_return_sector(state_list,istate)
                Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
                if(MpiStatus)then
                   gscvec => es_return_cvector(MpiComm,state_list,istate)
                else
                   gscvec => es_return_cvector(state_list,istate)
                endif
#else
                gscvec => es_return_cvector(state_list,istate)
#endif
                !
                peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
                peso = peso/zeta_function
                !
                idim    = getdim(isector)
                !
                if(Mpimaster)then
                   call build_sector(isector,H)
                   !GET <(C_UP + CDG_DW)(CDG_UP + C_DW)> = 
                   !<C_UP*CDG_UP> + <CDG_DW*C_DW> + <C_UP*C_DW> + <CDG_DW*CDG_UP> = 
                   !<N_UP> + < 1 - N_DW> + 2*<PHI>
                   isz = getsz(isector)
                   if(isz<Ns)then
                      jsz     = isz+1
                      jsector = getsector(jsz,1)
                      jdim    = getdim(jsector)
                      !
                      allocate(vvinit(jdim));vvinit=zero
                      call build_sector(jsector,HJ)
                      do i=1,idim
                         m=H%map(i)
                         ib = bdecomp(m,2*Ns)
                         if(ib(iorb)==0)then
                            call cdg(iorb,m,r,sgn)
                            j=binary_search(HJ%map,r)
                            vvinit(j) = sgn*gscvec(i)
                         endif
                      enddo
                      do i=1,idim
                         m=H%map(i)
                         ib = bdecomp(m,2*Ns)
                         if(ib(iorb+Ns)==1)then
                            call c(iorb+Ns,m,r,sgn)
                            j=binary_search(HJ%map,r)
                            vvinit(j) = vvinit(j) + sgn*gscvec(i)
                         endif
                      enddo
                      call delete_sector(jsector,HJ)
                      phisc(iorb) = phisc(iorb) + dot_product(vvinit,vvinit)*peso
                   endif
                   call delete_sector(isector,H)
                endif
                if(associated(gscvec)) nullify(gscvec)
                if(allocated(vvinit))deallocate(vvinit)
                !
             enddo
             phisc(iorb) = 0.5d0*(phisc(iorb) - dens_up(iorb) - (1.d0-dens_dw(iorb)))
          enddo
       enddo
    end if


    !EVALUATE <SX> AND <SY>
    if(ed_mode=="nonsu2")then
       ispin=1
       jspin=2
       do iorb=1,Norb
          isite=impIndex(iorb,1)
          jsite=impIndex(iorb,2)
          do istate=1,state_list%size
             !
             isector = es_return_sector(state_list,istate)
             Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
             if(MpiStatus)then
                gscvec => es_return_cvector(MpiComm,state_list,istate)
             else
                gscvec => es_return_cvector(state_list,istate)
             endif
#else
             gscvec => es_return_cvector(state_list,istate)
#endif
             !
             peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
             peso = peso/zeta_function
             !
             idim    = getdim(isector)
             !GET <(CDG_UP + CDG_DW)(C_UP + C_DW)> = 
             !<CDG_UP*C_UP> + <CDG_DW*C_DW> + <CDG_UP*C_DW + CDG_DW*C_UP> = 
             !<N_UP> + <N_DW> + 2*<Sx>
             !<Sx> = <CDG_UP*C_DW + CDG_DW*C_UP>
             jsector = getCsector(ispin,isector)
             jdim    = getdim(jsector)
             if(jsector/=0)then
                if(MPImaster)then
                   call build_sector(isector,H)
                   call build_sector(jsector,HJ)
                   allocate(vvinit(jdim));vvinit=zero
                   do i=1,idim
                      m=H%map(i)
                      ib = bdecomp(m,2*Ns)
                      if(ib(isite)==1)then
                         call c(isite,m,r,sgn)
                         j=binary_search(HJ%map,r)
                         vvinit(j) = sgn*gscvec(i)
                      endif
                   enddo
                   do i=1,idim
                      m=H%map(i)
                      ib = bdecomp(m,2*Ns)
                      if(ib(jsite)==1)then
                         call c(jsite,m,r,sgn)
                         j=binary_search(HJ%map,r)
                         vvinit(j) = vvinit(j) + sgn*gscvec(i)
                      endif
                   enddo
                   call delete_sector(isector,H)
                   call delete_sector(jsector,HJ)
                   magx(iorb) = magx(iorb) + dot_product(vvinit,vvinit)*peso
                   if(allocated(vvinit))deallocate(vvinit)
                endif
             endif
             !
             !GET <(-i*CDG_UP + CDG_DW)(i*C_UP + C_DW)> = 
             !<CDG_UP*C_UP> + <CDG_DW*C_DW> - i<CDG_UP*C_DW - CDG_DW*C_UP> = 
             !<N_UP> + <N_DW> + 2*<Sy>
             !<Sy> = -i/2<CDG_UP*C_DW - CDG_DW*C_UP>
             jsector = getCsector(ispin,isector)
             jdim    = getdim(jsector)
             if(jsector/=0)then
                if(MPImaster)then
                   call build_sector(isector,H)
                   call build_sector(jsector,HJ)
                   allocate(vvinit(jdim));vvinit=zero
                   do i=1,idim
                      m=H%map(i)
                      ib = bdecomp(m,2*Ns)
                      if(ib(isite)==1)then
                         call c(isite,m,r,sgn)
                         j=binary_search(HJ%map,r)
                         vvinit(j) = xi*sgn*gscvec(i)
                      endif
                   enddo
                   do i=1,idim
                      m=H%map(i)
                      ib = bdecomp(m,2*Ns)
                      if(ib(jsite)==1)then
                         call c(jsite,m,r,sgn)
                         j=binary_search(HJ%map,r)
                         vvinit(j) = vvinit(j) + sgn*gscvec(i)
                      endif
                   enddo
                   call delete_sector(isector,H)
                   call delete_sector(jsector,HJ)
                   magy(iorb) = magy(iorb) + dot_product(vvinit,vvinit)*peso
                   if(allocated(vvinit))deallocate(vvinit)
                endif
             endif
             !
             if(associated(gscvec)) nullify(gscvec)
             !
          enddo
          magx(iorb) = 0.5d0*(magx(iorb) - dens_up(iorb) - dens_dw(iorb))
          magy(iorb) = 0.5d0*(magy(iorb) - dens_up(iorb) - dens_dw(iorb))
       enddo
    end if


    !
    !EVALUATE EXCITON OP <S_ab> AND <T^x,y,z_ab>
    !<S_ab>  :=   <C^+_{a,up}C_{b,up} + C^+_{a,dw}C_{b,dw}>
    !<T^z_ab>:=   <C^+_{a,up}C_{b,up} - C^+_{a,dw}C_{b,dw}>
    !<T^x_ab>:=   <C^+_{a,up}C_{b,dw} + C^+_{a,dw}C_{b,up}>
    !<T^y_ab>:= -i<C^+_{a,up}C_{b,dw} - C^+_{a,dw}C_{b,up}>
    if(ed_mode=="nonsu2")then
       do istate=1,state_list%size
          isector = es_return_sector(state_list,istate)
          Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
          if(MpiStatus)then
             gscvec => es_return_cvector(MpiComm,state_list,istate)
          else
             gscvec => es_return_cvector(state_list,istate)
          endif
#else
          gscvec => es_return_cvector(state_list,istate)
#endif
          !
          peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
          peso = peso/zeta_function
          !
          idim    = getdim(isector)
          !    
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                !
                !\Theta_upup = <v|v>, |v> = (C_aup + C_bup)|>
                jsector = getCsector(1,isector)
                jdim    = getdim(jsector)
                if(jsector/=0)then
                   if(MPImaster)then
                      call build_sector(isector,H)
                      call build_sector(jsector,HJ)
                      allocate(vvinit(jdim));vvinit=zero
                      isite=impIndex(iorb,1)!+c_a,up
                      jsite=impIndex(jorb,1)! c_b,up
                      do i=1,idim
                         m  = H%map(i)
                         ib = bdecomp(m,2*Ns)
                         if(ib(jsite)/=1)cycle
                         call c(jsite,m,r,sgn)
                         j=binary_search(HJ%map,r)
                         vvinit(j) = sgn*gscvec(i)
                      enddo
                      do i=1,idim
                         m  = H%map(i)
                         ib = bdecomp(m,2*Ns)
                         if(ib(isite)/=1)cycle
                         call c(isite,m,r,sgn)
                         j=binary_search(HJ%map,r)
                         vvinit(j) = vvinit(j) + sgn*gscvec(i)
                      enddo
                      call delete_sector(isector,H)
                      call delete_sector(jsector,HJ)
                      theta_upup(iorb,jorb) = theta_upup(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                      if(allocated(vvinit))deallocate(vvinit)
                   endif
                endif
                !
                !\Theta_dwdw = <v|v>, |v> = (C_adw + C_bdw)|>
                jsector = getCsector(2,isector)
                jdim    = getdim(jsector)
                if(jsector/=0)then
                   if(MPImaster)then
                      call build_sector(isector,H)
                      call build_sector(jsector,HJ)
                      allocate(vvinit(jdim));vvinit=zero
                      isite=impIndex(iorb,2)!+c_a,dw
                      jsite=impIndex(jorb,2)! c_b,dw
                      do i=1,idim
                         m  = H%map(i)
                         ib = bdecomp(m,2*Ns)
                         if(ib(jsite)/=1)cycle
                         call c(jsite,m,r,sgn)
                         j=binary_search(HJ%map,r)
                         vvinit(j) = sgn*gscvec(i)
                      enddo
                      do i=1,idim
                         m  = H%map(i)
                         ib = bdecomp(m,2*Ns)
                         if(ib(isite)/=1)cycle
                         call c(isite,m,r,sgn)
                         j=binary_search(HJ%map,r)
                         vvinit(j) = vvinit(j) + sgn*gscvec(i)
                      enddo
                      call delete_sector(isector,H)
                      call delete_sector(jsector,HJ)
                      theta_dwdw(iorb,jorb) = theta_dwdw(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                      if(allocated(vvinit))deallocate(vvinit)
                   endif
                endif
                !
                !\Theta_updw = <v|v>, |v> = (C_aup + C_bdw)|>
                jsector = getCsector(1,isector)
                jdim    = getdim(jsector)
                if(jsector/=0)then
                   if(MPImaster)then
                      call build_sector(isector,H)
                      call build_sector(jsector,HJ)
                      allocate(vvinit(jdim));vvinit=zero
                      isite=impIndex(iorb,1)!+c_a,up
                      jsite=impIndex(jorb,2)! c_b,dw
                      do i=1,idim
                         m  = H%map(i)
                         ib = bdecomp(m,2*Ns)
                         if(ib(jsite)/=1)cycle
                         call c(jsite,m,r,sgn)
                         j=binary_search(HJ%map,r)
                         vvinit(j) = sgn*gscvec(i)
                      enddo
                      do i=1,idim
                         m  = H%map(i)
                         ib = bdecomp(m,2*Ns)
                         if(ib(isite)/=1)cycle
                         call c(isite,m,r,sgn)
                         j=binary_search(HJ%map,r)
                         vvinit(j) = vvinit(j) + sgn*gscvec(i)
                      enddo
                      call delete_sector(isector,H)
                      call delete_sector(jsector,HJ)
                      theta_updw(iorb,jorb) = theta_updw(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                      if(allocated(vvinit))deallocate(vvinit)
                   endif
                endif
                !
                !\Theta_dwup = <v|v>, |v> = (C_adw + C_bup)|>
                jsector = getCsector(1,isector)
                jdim    = getdim(jsector)
                if(jsector/=0)then
                   if(MPImaster)then
                      call build_sector(isector,H)
                      call build_sector(jsector,HJ)
                      allocate(vvinit(jdim));vvinit=zero
                      isite=impIndex(iorb,2)!+c_a,dw
                      jsite=impIndex(jorb,1)! c_b,up
                      do i=1,idim
                         m  = H%map(i)
                         ib = bdecomp(m,2*Ns)
                         if(ib(jsite)/=1)cycle
                         call c(jsite,m,r,sgn)
                         j=binary_search(HJ%map,r)
                         vvinit(j) = sgn*gscvec(i)
                      enddo
                      do i=1,idim
                         m  = H%map(i)
                         ib = bdecomp(m,2*Ns)
                         if(ib(isite)/=1)cycle
                         call c(isite,m,r,sgn)
                         j=binary_search(HJ%map,r)
                         vvinit(j) = vvinit(j) + sgn*gscvec(i)
                      enddo
                      call delete_sector(isector,H)
                      call delete_sector(jsector,HJ)
                      theta_dwup(iorb,jorb) = theta_dwup(iorb,jorb) + dot_product(vvinit,vvinit)*peso
                      if(allocated(vvinit))deallocate(vvinit)
                   endif
                endif
             enddo
          enddo
          !
          if(associated(gscvec)) nullify(gscvec)
          !
       enddo
       !
       if(MPImaster)then
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                exct_s0(iorb,jorb) = 0.5d0*(theta_upup(iorb,jorb) + theta_dwdw(iorb,jorb) - dens(iorb) - dens(jorb))
                exct_tz(iorb,jorb) = 0.5d0*(theta_upup(iorb,jorb) - theta_dwdw(iorb,jorb) - magZ(iorb) - magZ(jorb))
                exct_tx(iorb,jorb) = 0.5d0*(theta_updw(iorb,jorb) + theta_dwup(iorb,jorb) - dens(iorb) - dens(jorb))
                exct_ty(iorb,jorb) = -xi*0.5d0*(theta_updw(iorb,jorb) - theta_dwup(iorb,jorb) - magZ(iorb) + magZ(jorb))
             enddo
          enddo
       endif

       !
       !IMPURITY DENSITY MATRIX
       if(allocated(imp_density_matrix))deallocate(imp_density_matrix)
       allocate(imp_density_matrix(Nspin,Nspin,Norb,Norb))
       imp_density_matrix=zero
       do istate=1,state_list%size
          !
          isector = es_return_sector(state_list,istate)
          Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
          if(MpiStatus)then
             gscvec => es_return_cvector(MpiComm,state_list,istate)
          else
             gscvec => es_return_cvector(state_list,istate)
          endif
#else
          gscvec => es_return_cvector(state_list,istate)
#endif
          !
          peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
          peso = peso/zeta_function
          !
          idim  = getdim(isector)
          if(Mpimaster)then
             call build_sector(isector,H)
             !
             !Diagonal densities
             do ispin=1,Nspin
                do iorb=1,Norb
                   isite=impIndex(iorb,ispin)
                   do m=1,idim
                      i=H%map(m)
                      ib = bdecomp(i,2*Ns)
                      imp_density_matrix(ispin,ispin,iorb,iorb) = imp_density_matrix(ispin,ispin,iorb,iorb) + &
                           peso*ib(isite)*conjg(gscvec(m))*gscvec(m)
                   enddo
                enddo
             enddo
             !off-diagonal
             do ispin=1,Nspin
                do jspin=1,Nspin
                   do iorb=1,Norb
                      do jorb=1,Norb
                         if((ed_mode=="normal").and.(ispin/=jspin))cycle
                         if((bath_type=="normal").and.(iorb/=jorb))cycle
                         if(bath_type=="replica".and.Jz_basis)then
                            if((.not.dmft_bath%mask(ispin,jspin,iorb,jorb,1)).and.(.not.dmft_bath%mask(ispin,jspin,iorb,jorb,2)))cycle
                         endif
                         isite=impIndex(iorb,ispin)
                         jsite=impIndex(jorb,jspin)
                         do m=1,idim
                            i=H%map(m)
                            ib = bdecomp(i,2*Ns)
                            if((ib(jsite)==1).and.(ib(isite)==0))then
                               call c(jsite,i,r,sgn1)
                               call cdg(isite,r,k,sgn2)
                               j=binary_search(H%map,k)
                               imp_density_matrix(ispin,jspin,iorb,jorb) = imp_density_matrix(ispin,jspin,iorb,jorb) + &
                                    peso*sgn1*gscvec(m)*sgn2*conjg(gscvec(j))
                            endif
                         enddo
                      enddo
                   enddo
                enddo
             enddo
             !
             if(associated(gscvec))nullify(gscvec)
             call delete_sector(isector,H)
          endif
          !
       enddo
    endif
    !
    !
    !BATH DENSITY MATRIX (only if bath_type=="replica")
    if(bath_type=="replica")then
       if(allocated(bth_density_matrix)) deallocate(bth_density_matrix);allocate(bth_density_matrix(Nspin,Nspin,Norb,Norb,Nbath))
       bth_density_matrix=zero
       do istate=1,state_list%size
          !
          isector = es_return_sector(state_list,istate)
          Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
          if(MpiStatus)then
             gscvec => es_return_cvector(MpiComm,state_list,istate)
          else
             gscvec => es_return_cvector(state_list,istate)
          endif
#else
          gscvec => es_return_cvector(state_list,istate)
#endif
          !
          peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
          peso = peso/zeta_function
          !
          idim  = getdim(isector)
          !
          if(Mpimaster)then             
             call build_sector(isector,H)
             !
             !Diagonal densities
             do ibath=1,Nbath
                do ispin=1,Nspin
                   do iorb=1,Norb
                      isite = iorb + ibath*Norb + (ispin-1)*Ns
                      do m=1,idim
                         i=H%map(m)
                         ib = bdecomp(i,2*Ns)
                         bth_density_matrix(ispin,ispin,iorb,iorb,ibath) = bth_density_matrix(ispin,ispin,iorb,iorb,ibath) + &
                              peso*ib(isite)*conjg(gscvec(m))*gscvec(m)
                      enddo
                   enddo
                enddo
                !off-diagonal
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            if((ed_mode=="normal").and.(ispin/=jspin))cycle
                            if((bath_type=="normal").and.(iorb/=jorb))cycle
                            if(Jz_basis.and.(.not.dmft_bath%mask(ispin,jspin,iorb,jorb,1)).and.(.not.dmft_bath%mask(ispin,jspin,iorb,jorb,2)))cycle
                            isite = iorb + ibath*Norb + (ispin-1)*Ns
                            jsite = jorb + ibath*Norb + (jspin-1)*Ns
                            do m=1,idim
                               i=H%map(m)
                               ib = bdecomp(i,2*Ns)
                               if((ib(jsite)==1).and.(ib(isite)==0))then
                                  call c(jsite,i,r,sgn1)
                                  call cdg(isite,r,k,sgn2)
                                  j=binary_search(H%map,k)
                                  bth_density_matrix(ispin,jspin,iorb,jorb,ibath) = bth_density_matrix(ispin,jspin,iorb,jorb,ibath) + &
                                       peso*sgn1*gscvec(m)*sgn2*conjg(gscvec(j))
                               endif
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
             !
             if(associated(gscvec))nullify(gscvec)
             call delete_sector(isector,H)
          endif
       enddo
    endif
    !
    if(MPIMASTER)then
       call get_szr
       if(iolegend)call write_legend
       call write_observables()

       write(LOGfile,"(A,10f18.12,f18.12,A)")"dens"//reg(ed_file_suffix)//"=",(dens(iorb),iorb=1,Norb),sum(dens)
       select case(ed_mode)
       case default
          write(LOGfile,"(A,10f18.12,A)")    "docc"//reg(ed_file_suffix)//"=",(docc(iorb),iorb=1,Norb)
       case("superc")
          write(LOGfile,"(A,20f18.12,A)")    "phi "//reg(ed_file_suffix)//"=",(phisc(iorb),iorb=1,Norb),(abs(uloc(iorb))*phisc(iorb),iorb=1,Norb)
       case("nonsu2")
          write(LOGfile,"(A,10f18.12,A)")    "magX"//reg(ed_file_suffix)//"=",(magX(iorb),iorb=1,Norb)
          write(LOGfile,"(A,10f18.12,A)")    "magY"//reg(ed_file_suffix)//"=",(magY(iorb),iorb=1,Norb)
          if(Nspin==2)then
             write(LOGfile,"(A,10f18.12,A)")    "magZ"//reg(ed_file_suffix)//"=",(magz(iorb),iorb=1,Norb)
          endif
          write(LOGfile,"(A,20f18.12,A)")    "exS0"//reg(ed_file_suffix)//"=",((exct_S0(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
          write(LOGfile,"(A,20f18.12,A)")    "exTz"//reg(ed_file_suffix)//"=",((exct_Tz(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
          write(LOGfile,"(A,20f18.12,A)")    "exTx"//reg(ed_file_suffix)//"=",((exct_Tx(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
          write(LOGfile,"(A,20f18.12,A)")    "exTy"//reg(ed_file_suffix)//"=",((exct_Ty(iorb,jorb),jorb=iorb+1,Norb),iorb=1,Norb)
       end select

    endif
    !
    do iorb=1,Norb
       ed_dens_up(iorb)=dens_up(iorb)
       ed_dens_dw(iorb)=dens_dw(iorb)
       ed_dens(iorb)   =dens(iorb)
       ed_docc(iorb)   =docc(iorb)
       ed_phisc(iorb)  =phisc(iorb)
       ed_mag(1,iorb)  =magX(iorb)
       ed_mag(2,iorb)  =magY(iorb)
       ed_mag(3,iorb)  =magZ(iorb)
    enddo
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_dens_up)
       call Bcast_MPI(MpiComm,ed_dens_dw)
       call Bcast_MPI(MpiComm,ed_dens)
       call Bcast_MPI(MpiComm,ed_docc)
       call Bcast_MPI(MpiComm,ed_phisc)
       call Bcast_MPI(MpiComm,ed_mag)
       if(allocated(imp_density_matrix))call Bcast_MPI(MpiComm,imp_density_matrix)
       if(bath_type=="replica")then
          call Bcast_MPI(MpiComm,bth_density_matrix)
       endif
    endif
#endif
    !
    deallocate(dens,docc,phisc,dens_up,dens_dw,magz,sz2,n2,magx,magy,exct_s0,exct_tx,exct_ty,exct_tz)
    deallocate(simp,zimp)
  end subroutine observables_impurity







  !+-------------------------------------------------------------------+
  !PURPOSE  : Get internal energy from the Impurity problem.
  !+-------------------------------------------------------------------+
  subroutine local_energy_impurity()
    integer,dimension(Nlevels)      :: ib
    integer                         :: i,j
    integer                         :: istate
    integer                         :: isector
    integer                         :: idim
    integer                         :: iorb,jorb,ispin
    integer                         :: numstates
    integer                         :: m,k1,k2,k3,k4
    real(8)                         :: sg1,sg2,sg3,sg4
    real(8)                         :: Egs,gs_weight
    real(8)                         :: Ei
    real(8)                         :: peso
    real(8)                         :: norm
    real(8),dimension(Norb)         :: nup,ndw
    real(8),dimension(Nspin,Norb)   :: eloc
    complex(8),dimension(:),pointer :: gscvec
    type(sector_map)                :: H
    logical                         :: Jcondition
    !
    Egs     = state_list%emin
    ed_Ehartree= 0.d0
    ed_Eknot   = 0.d0
    ed_Epot    = 0.d0
    ed_Dust    = 0.d0
    ed_Dund    = 0.d0
    ed_Dse     = 0.d0
    ed_Dph     = 0.d0
    !
    !Get diagonal part of Hloc
    do ispin=1,Nspin
       do iorb=1,Norb
          eloc(ispin,iorb)=impHloc(ispin,ispin,iorb,iorb)
       enddo
    enddo
    !
    numstates=state_list%size
    do istate=1,numstates
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          gscvec => es_return_cvector(MpiComm,state_list,istate)
       else
          gscvec => es_return_cvector(state_list,istate)
       endif
#else
       gscvec => es_return_cvector(state_list,istate)
#endif
       !
       idim  = getdim(isector)
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(Mpimaster)then
          !
          call build_sector(isector,H)
          !
          do i=1,idim
             m=H%map(i)
             ib = bdecomp(m,2*Ns)
             !
             gs_weight=peso*abs(gscvec(i))**2
             !
             !Get operators:
             do iorb=1,Norb
                nup(iorb)= dble(ib(iorb))
                ndw(iorb)= dble(ib(iorb+Ns))
             enddo
             !
             !start evaluating the Tr(H_loc) to estimate potential energy
             !
             !LOCAL ENERGY
             ed_Eknot = ed_Eknot + dot_product(eloc(1,:),nup)*gs_weight + dot_product(eloc(Nspin,:),ndw)*gs_weight
             !==> HYBRIDIZATION TERMS I: same or different orbitals, same spins.
             do iorb=1,Norb
                do jorb=1,Norb
                   !SPIN UP
                   if((ib(iorb)==0).AND.(ib(jorb)==1))then
                      call c(jorb,m,k1,sg1)
                      call cdg(iorb,k1,k2,sg2)
                      j=binary_search(H%map,k2)
                      if(Jz_basis.and.j==0)cycle
                      !WARNING: note that the previous line, and all the other hereafter, are equivalent to:
                      !if(Jz_basis.and.(.not.dmft_bath%mask(1,1,iorb,jorb,1)).and.(.not.dmft_bath%mask(1,1,iorb,jorb,2)))cycle
                      ed_Eknot = ed_Eknot + impHloc(1,1,iorb,jorb)*sg1*sg2*gscvec(i)*conjg(gscvec(j))
                   endif
                   !SPIN DW
                   if((ib(iorb+Ns)==0).AND.(ib(jorb+Ns)==1))then
                      call c(jorb+Ns,m,k1,sg1)
                      call cdg(iorb+Ns,k1,k2,sg2)
                      j=binary_search(H%map,k2)
                      if(Jz_basis.and.j==0)cycle
                      ed_Eknot = ed_Eknot + impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2*gscvec(i)*conjg(gscvec(j))
                   endif
                enddo
             enddo
             !==> HYBRIDIZATION TERMS II: same or different orbitals, opposite spins.
             if(ed_mode=="nonsu2")then
                do iorb=1,Norb
                   do jorb=1,Norb
                      !UP-DW
                      if((impHloc(1,Nspin,iorb,jorb)/=zero).AND.(ib(iorb)==0).AND.(ib(jorb+Ns)==1))then
                         call c(jorb+Ns,m,k1,sg1)
                         call cdg(iorb,k1,k2,sg2)
                         j=binary_search(H%map,k2)
                         if(Jz_basis.and.j==0)cycle
                         ed_Eknot = ed_Eknot + impHloc(1,Nspin,iorb,jorb)*sg1*sg2*gscvec(i)*conjg(gscvec(j))
                      endif
                      !DW-UP
                      if((impHloc(Nspin,1,iorb,jorb)/=zero).AND.(ib(iorb+Ns)==0).AND.(ib(jorb)==1))then
                         call c(jorb,m,k1,sg1)
                         call cdg(iorb+Ns,k1,k2,sg2)
                         j=binary_search(H%map,k2)
                         if(Jz_basis.and.j==0)cycle
                         ed_Eknot = ed_Eknot + impHloc(Nspin,1,iorb,jorb)*sg1*sg2*gscvec(i)*conjg(gscvec(j))
                      endif
                   enddo
                enddo
             endif
             !
             !DENSITY-DENSITY INTERACTION: SAME ORBITAL, OPPOSITE SPINS
             !Euloc=\sum=i U_i*(n_u*n_d)_i
             !ed_Epot = ed_Epot + dot_product(uloc,nup*ndw)*gs_weight
             do iorb=1,Norb
                ed_Epot = ed_Epot + Uloc(iorb)*nup(iorb)*ndw(iorb)*gs_weight
             enddo
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, OPPOSITE SPINS
             !Eust=\sum_ij Ust*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             !    "="\sum_ij (Uloc - 2*Jh)*(n_up_i*n_dn_j + n_up_j*n_dn_i)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      ed_Epot = ed_Epot + Ust*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                      ed_Dust = ed_Dust + (nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))*gs_weight
                   enddo
                enddo
             endif
             !
             !DENSITY-DENSITY INTERACTION: DIFFERENT ORBITALS, PARALLEL SPINS
             !Eund = \sum_ij Und*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Ust-Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             !    "="\sum_ij (Uloc-3*Jh)*(n_up_i*n_up_j + n_dn_i*n_dn_j)
             if(Norb>1)then
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      ed_Epot = ed_Epot + (Ust-Jh)*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                      ed_Dund = ed_Dund + (nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))*gs_weight
                   enddo
                enddo
             endif
             !
             !SPIN-EXCHANGE (S-E) TERMS
             !S-E: Jh *( c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up )  (i.ne.j) 
             if(Norb>1.AND.Jhflag)then
                do iorb=1,Norb
                   do jorb=1,Norb
                      Jcondition=((iorb/=jorb).AND.&
                           (ib(jorb)==1)      .AND.&
                           (ib(iorb+Ns)==1)   .AND.&
                           (ib(jorb+Ns)==0)   .AND.&
                           (ib(iorb)==0))
                      if(Jcondition)then
                         call c(jorb,m,k1,sg1)
                         call c(iorb+Ns,k1,k2,sg2)
                         call cdg(jorb+Ns,k2,k3,sg3)
                         call cdg(iorb,k3,k4,sg4)
                         j=binary_search(H%map,k4)
                         if(Jz_basis.and.j==0)cycle
                         ed_Epot = ed_Epot + Jx*sg1*sg2*sg3*sg4*gscvec(i)*conjg(gscvec(j))!gs_weight
                         ed_Dse  = ed_Dse  + sg1*sg2*sg3*sg4*gscvec(i)*conjg(gscvec(j))!gs_weight
                      endif
                   enddo
                enddo
             endif
             !
             !PAIR-HOPPING (P-H) TERMS
             !P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
             !P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
             if(Norb>1.AND.Jhflag)then
                do iorb=1,Norb
                   do jorb=1,Norb
                      Jcondition=((iorb/=jorb).AND.&
                           (ib(jorb)==1)      .AND.&
                           (ib(jorb+Ns)==1)   .AND.&
                           (ib(iorb+Ns)==0)   .AND.&
                           (ib(iorb)==0))
                      if(Jcondition)then
                         call c(jorb,m,k1,sg1)
                         call c(jorb+Ns,k1,k2,sg2)
                         call cdg(iorb+Ns,k2,k3,sg3)
                         call cdg(iorb,k3,k4,sg4)
                         j=binary_search(H%map,k4)
                         if(Jz_basis.and.j==0)cycle
                         ed_Epot = ed_Epot + Jp*sg1*sg2*sg3*sg4*gscvec(i)*conjg(gscvec(j))!gs_weight
                         ed_Dph  = ed_Dph  + sg1*sg2*sg3*sg4*gscvec(i)*conjg(gscvec(j))!gs_weight
                      endif
                   enddo
                enddo
             endif
             !
             !HARTREE-TERMS CONTRIBUTION:
             if(hfmode)then
                !ed_Ehartree=ed_Ehartree - 0.5d0*dot_product(uloc,nup+ndw)*gs_weight + 0.25d0*sum(uloc)*gs_weight
                do iorb=1,Norb
                   ed_Ehartree=ed_Ehartree - 0.5d0*uloc(iorb)*(nup(iorb)+ndw(iorb))*gs_weight + 0.25d0*uloc(iorb)*gs_weight
                enddo
                if(Norb>1)then
                   do iorb=1,Norb
                      do jorb=iorb+1,Norb
                         ed_Ehartree=ed_Ehartree - 0.5d0*Ust*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.25d0*Ust*gs_weight
                         ed_Ehartree=ed_Ehartree - 0.5d0*(Ust-Jh)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))*gs_weight + 0.25d0*(Ust-Jh)*gs_weight
                      enddo
                   enddo
                endif
             endif
          enddo
          if(associated(gscvec))nullify(gscvec)
          call delete_sector(isector,H)
       endif
    enddo
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_Epot)
       call Bcast_MPI(MpiComm,ed_Eknot)
       call Bcast_MPI(MpiComm,ed_Ehartree)
       call Bcast_MPI(MpiComm,ed_Dust)
       call Bcast_MPI(MpiComm,ed_Dund)
       call Bcast_MPI(MpiComm,ed_Dse)
       call Bcast_MPI(MpiComm,ed_Dph)
    endif
#endif
    !
    ed_Epot = ed_Epot + ed_Ehartree
    !
    if(ed_verbose==3)then
       write(LOGfile,"(A,10f18.12)")"<Hint>  =",ed_Epot
       write(LOGfile,"(A,10f18.12)")"<V>     =",ed_Epot-ed_Ehartree
       write(LOGfile,"(A,10f18.12)")"<E0>    =",ed_Eknot
       write(LOGfile,"(A,10f18.12)")"<Ehf>   =",ed_Ehartree    
       write(LOGfile,"(A,10f18.12)")"Dust    =",ed_Dust
       write(LOGfile,"(A,10f18.12)")"Dund    =",ed_Dund
       write(LOGfile,"(A,10f18.12)")"Dse     =",ed_Dse
       write(LOGfile,"(A,10f18.12)")"Dph     =",ed_Dph
    endif
    if(MPIMASTER)then
       call write_energy_info()
       call write_energy()
    endif
    !
    !
  end subroutine local_energy_impurity



  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : get scattering rate and renormalization constant Z
  !+-------------------------------------------------------------------+
  subroutine get_szr()
    integer                  :: ispin,iorb
    real(8)                  :: wm1,wm2
    wm1 = pi/beta ; wm2=3d0*pi/beta
    do ispin=1,Nspin
       do iorb=1,Norb
          simp(iorb,ispin) = dimag(impSmats(ispin,ispin,iorb,iorb,1)) - &
               wm1*(dimag(impSmats(ispin,ispin,iorb,iorb,2))-dimag(impSmats(ispin,ispin,iorb,iorb,1)))/(wm2-wm1)
          zimp(iorb,ispin)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ispin,ispin,iorb,iorb,1))/wm1 ))
       enddo
    enddo
  end subroutine get_szr



  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_legend()
    integer :: unit,iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_info.ed")
    select case(ed_mode)
    case default
       write(unit,"(A1,90(A10,6X))")"#",&
            (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(Norb+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(2*Norb+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(3*Norb+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(4*Norb+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
            reg(txtfy(5*Norb+1))//"s2",&
            reg(txtfy(5*Norb+2))//"egs",&
            ((reg(txtfy(5*Norb+2+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
            ((reg(txtfy((5+Norb)*Norb+2+(iorb-1)*Norb+jorb))//"n2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
            ((reg(txtfy((5+2*Norb)*Norb+2+(ispin-1)*Nspin+iorb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
            ((reg(txtfy((6+2*Norb)*Norb+2+Nspin+(ispin-1)*Nspin+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    case ("superc")
       write(unit,"(A1,90(A10,6X))")"#",&
            (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(Norb+iorb))//"phi_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(2*Norb+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(3*Norb+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(4*Norb+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(5*Norb+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
            reg(txtfy(6*Norb+1))//"s2",&
            reg(txtfy(6*Norb+2))//"egs",&
            ((reg(txtfy(6*Norb+2+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
            ((reg(txtfy((6+Norb)*Norb+2+(iorb-1)*Norb+jorb))//"n2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
            ((reg(txtfy((6+2*Norb)*Norb+2+(ispin-1)*Nspin+iorb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
            ((reg(txtfy((7+2*Norb)*Norb+2+Nspin+(ispin-1)*Nspin+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    case("nonsu2")
       write(unit,"(A1,90(A10,6X))")"#",&
            (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(Norb+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(2*Norb+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(3*Norb+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(4*Norb+iorb))//"magX_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(5*Norb+iorb))//"magY_"//reg(txtfy(iorb)),iorb=1,Norb),&
            (reg(txtfy(6*Norb+iorb))//"magZ_"//reg(txtfy(iorb)),iorb=1,Norb),&            
            reg(txtfy(7*Norb+1))//"s2",&
            reg(txtfy(7*Norb+2))//"egs",&
            ((reg(txtfy(7*Norb+2+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
            ((reg(txtfy((7+Norb)*Norb+2+(iorb-1)*Norb+jorb))//"n2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
            ((reg(txtfy((7+2*Norb)*Norb+2+(ispin-1)*Nspin+iorb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
            ((reg(txtfy((8+2*Norb)*Norb+2+Nspin+(ispin-1)*Nspin+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    end select
    close(unit)
    !
    if(ed_mode=="nonsu2")then
       open(unit,file="exciton_info.ed")
       write(unit,"(A1,6(A10,6X))")"#","1S_0","2T_z","3reT_x","4imT_x","5reT_y","6imT_y"
       close(unit)
    endif
    !
    unit = free_unit()
    open(unit,file="parameters_info.ed")
    write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
         (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
    close(unit)
    !
    iolegend=.false.
  end subroutine write_legend

  subroutine write_energy_info()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_info.ed")
    write(unit,"(A1,90(A14,1X))")"#",&
         reg(txtfy(1))//"<Hi>",&
         reg(txtfy(2))//"<V>=<Hi-Ehf>",&
         reg(txtfy(3))//"<Eloc>",&
         reg(txtfy(4))//"<Ehf>",&
         reg(txtfy(5))//"<Dst>",&
         reg(txtfy(6))//"<Dnd>",&
         reg(txtfy(7))//"<Dse>",&
         reg(txtfy(8))//"<Dph>"
    close(unit)
  end subroutine write_energy_info


  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    integer :: unit
    integer :: iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_all"//reg(ed_file_suffix)//".ed",position='append')
    select case(ed_mode)
    case default
       write(unit,"(90(F15.9,1X))")&
            (dens(iorb),iorb=1,Norb),&
            (docc(iorb),iorb=1,Norb),&
            (dens_up(iorb),iorb=1,Norb),&
            (dens_dw(iorb),iorb=1,Norb),&
            (magz(iorb),iorb=1,Norb),&
            s2tot,egs,&
            ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
    case ("superc")
       write(unit,"(90(F15.9,1X))")&
            (dens(iorb),iorb=1,Norb),&
            (phisc(iorb),iorb=1,Norb),&
            (docc(iorb),iorb=1,Norb),&
            (dens_up(iorb),iorb=1,Norb),&
            (dens_dw(iorb),iorb=1,Norb),&
            (magz(iorb),iorb=1,Norb),&
            s2tot,egs,&
            ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
    case ("nonsu2")
       write(unit,"(90(F15.9,1X))")&
            (dens(iorb),iorb=1,Norb),&
            (docc(iorb),iorb=1,Norb),&
            (dens_up(iorb),iorb=1,Norb),&
            (dens_dw(iorb),iorb=1,Norb),&
            (magX(iorb),iorb=1,Norb),&
            (magY(iorb),iorb=1,Norb),&
            (magZ(iorb),iorb=1,Norb),&
            s2tot,egs,&
            ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
       !
       unit = free_unit()
       open(unit,file="exciton_last"//reg(ed_file_suffix)//".ed")
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             write(unit,"(90(F15.9,1X))")&
                  exct_s0(iorb,jorb),exct_tz(iorb,jorb),&
                  dreal(exct_tx(iorb,jorb)),dimag(exct_tx(iorb,jorb)),&
                  dreal(exct_ty(iorb,jorb)),dimag(exct_ty(iorb,jorb))
          enddo
       enddo
       close(unit)
       !
    end select
    close(unit)    
    !
    unit = free_unit()
    open(unit,file="parameters_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp
    close(unit)
    !
    unit = free_unit()
    open(unit,file="observables_last"//reg(ed_file_suffix)//".ed")
    select case(ed_mode)
    case default
       write(unit,"(90(F15.9,1X))")&
            (dens(iorb),iorb=1,Norb),&
            (docc(iorb),iorb=1,Norb),&
            (dens_up(iorb),iorb=1,Norb),&
            (dens_dw(iorb),iorb=1,Norb),&
            (magz(iorb),iorb=1,Norb),&
            s2tot,egs,&
            ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    case ("superc")
       write(unit,"(90(F15.9,1X))")&
            (dens(iorb),iorb=1,Norb),&
            (phisc(iorb),iorb=1,Norb),&
            (docc(iorb),iorb=1,Norb),&
            (dens_up(iorb),iorb=1,Norb),&
            (dens_dw(iorb),iorb=1,Norb),&
            (magz(iorb),iorb=1,Norb),&
            s2tot,egs,&
            ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    case ("nonsu2")
       write(unit,"(90(F15.9,1X))")&
            (dens(iorb),iorb=1,Norb),&
            (docc(iorb),iorb=1,Norb),&
            (dens_up(iorb),iorb=1,Norb),&
            (dens_dw(iorb),iorb=1,Norb),&
            (magX(iorb),iorb=1,Norb),&
            (magY(iorb),iorb=1,Norb),&
            (magZ(iorb),iorb=1,Norb),&
            s2tot,egs,&
            ((sz2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((n2(iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((zimp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
    end select
    close(unit)         
  end subroutine write_observables

  subroutine write_energy()
    integer :: unit
    unit = free_unit()
    open(unit,file="energy_last"//reg(ed_file_suffix)//".ed")
    write(unit,"(90F15.9)")ed_Epot,ed_Epot-ed_Ehartree,ed_Eknot,ed_Ehartree,ed_Dust,ed_Dund,ed_Dse,ed_Dph
    close(unit)
  end subroutine write_energy



end MODULE ED_OBSERVABLES
