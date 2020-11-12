!########################################################################
!PURPOSE  : Obtain some physical quantities and print them out
!########################################################################
MODULE ED_OBSERVABLES
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy,str
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
  public :: observables_plaquette


  logical,save                       :: iolegend=.true.
  real(8),dimension(:),allocatable   :: dens,dens_up,dens_dw
  real(8),dimension(:),allocatable   :: docc
  real(8),dimension(:),allocatable   :: magz
  real(8),dimension(:),allocatable   :: phisc
  real(8),dimension(:,:),allocatable :: sz2,n2
  real(8),dimensioN(:,:),allocatable :: zimp,simp
  real(8)                            :: s2tot
  real(8)                            :: Egs
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
    !
    !LOCAL OBSERVABLES:
    ! density,
    ! double occupancy,
    ! magnetization,
    ! orbital//spin correlations
    ! superconducting order parameter, etc..
    allocate(dens(Norb),dens_up(Norb),dens_dw(Norb))
    allocate(docc(Norb))
    allocate(phisc(Norb))
    allocate(magz(Norb),sz2(Norb,Norb),n2(Norb,Norb))
    allocate(simp(Norb,Nspin),zimp(Norb,Nspin))
    !
    Egs     = state_list%emin
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    phisc   = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    n2      = 0.d0
    s2tot   = 0.d0
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
    endif
    write(LOGfile,"(A,10f18.12,f18.12,A)")"dens"//reg(ed_file_suffix)//"=",(dens(iorb),iorb=1,Norb),sum(dens)
    select case(ed_mode)
    case default
       write(LOGfile,"(A,10f18.12,A)")    "docc"//reg(ed_file_suffix)//"=",(docc(iorb),iorb=1,Norb)
    case("superc")
       write(LOGfile,"(A,20f18.12,A)")    "phi "//reg(ed_file_suffix)//"=",(phisc(iorb),iorb=1,Norb),(abs(uloc(iorb))*phisc(iorb),iorb=1,Norb)
    end select
    if(Nspin==2)then
       write(LOGfile,"(A,10f18.12,A)") "mag "//reg(ed_file_suffix)//"=",(magz(iorb),iorb=1,Norb)
    endif
    !
    do iorb=1,Norb
       ed_dens_up(iorb)=dens_up(iorb)
       ed_dens_dw(iorb)=dens_dw(iorb)
       ed_dens(iorb)   =dens(iorb)
       ed_docc(iorb)   =docc(iorb)
       ed_phisc(iorb)  =phisc(iorb)
    enddo
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,ed_dens_up)
       call Bcast_MPI(MpiComm,ed_dens_dw)
       call Bcast_MPI(MpiComm,ed_dens)
       call Bcast_MPI(MpiComm,ed_docc)
       call Bcast_MPI(MpiComm,ed_phisc)
       if(allocated(imp_density_matrix))call Bcast_MPI(MpiComm,imp_density_matrix)
       if(bath_type=="replica")then
          call Bcast_MPI(MpiComm,bth_density_matrix)
       endif
    endif
#endif
    !
    deallocate(dens,docc,phisc,dens_up,dens_dw,magz,sz2,n2)
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
    end select
    close(unit)
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




  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine observables_plaquette()
    integer,dimension(Nlevels)      :: ib
    integer                         :: i,j
    integer                         :: istate
    integer                         :: isector,jsector
    integer                         :: idim,jdim
    integer                         :: numstates
    integer                         :: c1,c2,c3,c4
    integer                         :: cd1,cd2,cd3,cd4
    integer(16)                     :: r,m,k
    integer(16)                     :: k1,k2,k3
    integer(16)                     :: r1,r2,r3
    real(8)                         :: sgn,sgn1,sgn2
    real(8)                         :: sgn3,sgn4,sgn5
    real(8)                         :: sgn6,sgn7,sgn8
    real(8)                         :: gs_weight
    real(8)                         :: Ei
    real(8)                         :: peso
    complex(8),dimension(:),pointer :: gscvec
    type(sector_map)                :: H
    type(sector_map8)               :: H8
    logical                         :: Jcondition,Jleg
    !
    integer                         :: distance
    integer                         :: isite,jsite,hopndx,ineig
    integer                         :: row,col
    integer                         :: unit_,ier
    integer                         :: BosonExp
    real(8)                         :: n1i,n2i,n3i,n4i
    real(8)                         :: n1j,n2j,n3j,n4j
    !
    real(8)                         :: NormCorr
    real(8)                         :: ntot,Ektot,Eptot,Xop
    real(8)                         :: NNr,GGr,XXr,OOx,OOy,PPr,absOOx,absOOy
    !
    real(8),allocatable             :: denslatt(:,:,:),denslatt_tot(:,:)
    real(8),allocatable             :: Eklatt(:,:,:)  ,Epotlatt(:,:,:,:)
    real(8),allocatable             :: Ncorrfunc(:,:)   !single particle density-density
    real(8),allocatable             :: Gcorrfunc(:,:)   !single particle propagator
    real(8),allocatable             :: Ocorrfunc(:,:,:) !single particle ordpar
    real(8),allocatable             :: Xcorrfunc(:,:)   !quartet density-density
    real(8),allocatable             :: Pcorrfunc(:,:)   !quartet ordpar
    !
    !
    Egs       = state_list%emin
    !
    numstates=state_list%size
    !
    allocate(denslatt(Nbath,Norb,numstates))     ; denslatt  = 0.d0
    allocate(denslatt_tot(Nbath,Norb))           ; denslatt_tot = 0.d0
    allocate(Ncorrfunc(size(Neigh),numstates))   ; Ncorrfunc = 0d0 ; NNr = 0d0
    allocate(Gcorrfunc(size(Neigh),numstates))   ; Gcorrfunc = 0d0 ; GGr = 0d0
    allocate(Ocorrfunc(size(Neigh),numstates,2)) ; Ocorrfunc = 0d0 ; OOx = 0d0 ; OOy = 0d0
    allocate(Xcorrfunc(size(Neigh),numstates))   ; Xcorrfunc = 0d0 ; XXr = 0d0
    allocate(Pcorrfunc(size(Neigh),numstates))   ; Pcorrfunc = 0d0 ; PPr = 0d0
    !
    allocate(Eklatt(Nbath,Norb,numstates))       ; Eklatt    = 0.d0
    allocate(Epotlatt(Nbath,Norb,0:size(Neigh),numstates)) ; Epotlatt = 0.d0
    !
    Xop = 0d0
    !
    BosonExp=1
    if(HardCoreBoson.ne.0)BosonExp=0
    !
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
          !
          !write(*,*)Ei,Egs,zeta_function
          call build_sector(isector,H,H8)
          !do i=1,idim
          !   m=H8%map(i)
          !   ib = bdecomp8(m,Ns)
          !   write(*,"(A,50I5)")"i,idim,m,Ns  ", i,idim,m,Ns,size(H8%map),(ib(isite),isite=1,Norb*Nbath)
          !enddo
          !
          do i=1,idim
             m=H8%map(i)
             ib = bdecomp8(m,Ns)
             !
             gs_weight=peso*abs(gscvec(i))**2
             !
             do isite=1,Norb*Nbath
                !
                row = vec2lat(isite,1)
                col = vec2lat(isite,2)
                !
                if(isite.ne.Xstride(isite,1)) stop "something wrong"
                n1i = dble(ib(Xstride(isite,1)))
                n2i = dble(ib(Xstride(isite,2)))
                n3i = dble(ib(Xstride(isite,3)))
                if(Nbath.ne.1)n4i = dble(ib(Xstride(isite,4)))
                !
                if(isnan(n1i).or.isnan(n2i).or.isnan(n3i).or.isnan(n4i))then
                   write(*,"(4F6.2)") n1i,n2i,n3i,n4i
                   write(*,"(4I6)") distance,ineig,jsite,isite
                   stop "ni"
                endif
                !
                !
                !Local observables on isite
                denslatt(row,col,istate) = denslatt(row,col,istate) + n1i*gs_weight
                denslatt_tot(row,col) = denslatt_tot(row,col) + n1i*gs_weight
                Xop = Xop + (n1i*n2i*n3i*n4i) * gs_weight/(Nbath*Norb)
                Epotlatt(row,col,0,istate) = Epotlatt(row,col,0,istate) + Ust * Radius(row,col) * n1i*gs_weight
                !
                !
                !Non-local observables. Loop over the sites located on the different radius
                do distance=1,size(Neigh)
                   do ineig=1,Neigh(distance)
                      !
                      jsite = Vstride(isite,distance,ineig)
                      !
                      n1j = dble(ib(Xstride(jsite,1)))
                      n2j = dble(ib(Xstride(jsite,2)))
                      n3j = dble(ib(Xstride(jsite,3)))
                      if(Nbath.ne.1)n4j = dble(ib(Xstride(jsite,4)))
                      !
                      if(isnan(n1j).or.isnan(n2j).or.isnan(n3j).or.isnan(n4j))then
                         write(*,"(4F6.2)") n1j,n2j,n3j,n4j
                         write(*,"(4I6)") distance,ineig,jsite,isite
                         stop "nj"
                      endif
                      !
                      !
                      !Single particle correlation function: <N(R)N(0)>
                      Ncorrfunc(distance,istate) = Ncorrfunc(distance,istate) + (n1i*n1j)*gs_weight / Neigh(distance)
                      !
                      !
                      !Quartet correlation function: <X(R)X(0)> - REDO_LADDER
                      !if(.not.((Nbath.eq.2).and.(mod(distance,2).eq.0)))then
                      Xcorrfunc(distance,istate) = Xcorrfunc(distance,istate) &
                                                 + (n1i*n2i*n3i*n4i) * (n1j*n2j*n3j*n4j) * gs_weight / (Nbath*Norb)
                      !
                      !
                      !Potential energy: same as <N(R)N(0)> but with the energy at a given R
                      Epotlatt(row,col,distance,istate) = Epotlatt(row,col,distance,istate) + (n1i*n1j)*gs_weight*Vmat(distance)
                      !
                      !
                      !Kinetic energy: hopping considered only for first neighbors
                      if(distance.eq.1)then
                         hopndx = Vstride(isite,distance,ineig)
                         if((ib(isite)==1).and.(ib(hopndx)==0))then
                            call c8(isite,m,r,sgn1)
                            call cdg8(hopndx,r,k,sgn2)
                            j=binary_search8(H8%map,k)
                            if((j.eq.0).and.(HardCoreBoson.eq.0))write(*,'(A,14I06)')" Ekin ",i,j,m,r,k,isite,jsite,hopndx
                            if(j.ne.0)Eklatt(row,col,istate) = Eklatt(row,col,istate) + peso*(sgn1**BosonExp)*gscvec(i)*(sgn2**BosonExp)*conjg(gscvec(j))*Thopping
                         endif
                      endif
                      !
                      !
                      !Single particle static propagator
                      hopndx = Vstride(isite,distance,ineig)
                      if((ib(isite)==1).and.(ib(hopndx)==0))then
                         call c8(isite,m,r,sgn1)
                         call cdg8(hopndx,r,k,sgn2)
                         j=binary_search8(H8%map,k)
                         if((j.eq.0).and.(HardCoreBoson.eq.0))write(*,'(A,14I06)')" Gprop ",i,j,m,r,k,isite,jsite,hopndx
                         if(j.ne.0)Gcorrfunc(distance,istate) = Gcorrfunc(distance,istate) + peso*(sgn1**BosonExp)*gscvec(i)*(sgn2**BosonExp)*conjg(gscvec(j))/Neigh(distance)
                      endif
                      !
                      !Pair correlations
                      if(Nbath.ne.1)then
                         !
                         c1 = Xstride(isite,1)
                         c2 = Xstride(isite,2)
                         c3 = Xstride(isite,3)
                         c4 = Xstride(isite,4)
                         !
                         cd1 = Xstride(jsite,1)
                         cd2 = Xstride(jsite,2)
                         cd3 = Xstride(jsite,3)
                         cd4 = Xstride(jsite,4)
                         !
                         !for the full quartet - REDO_LADDER
                         Jleg=.true.
                         !if((Nbath.eq.2).and.(mod(distance,2).eq.0))Jleg=.false.
                         Jcondition = (ib(c1) ==1).and.(ib(c2) ==1).and.(ib(c3) ==1).and.(ib(c4) ==1).and. &
                                      (ib(cd1)==0).and.(ib(cd2)==0).and.(ib(cd3)==0).and.(ib(cd4)==0).and.(filling.gt.4).and.Jleg
                         if(Jcondition)then
                            call   c8(c1 ,m ,r1,sgn1)
                            call   c8(c2 ,r1,r2,sgn2)
                            call   c8(c3 ,r2,r3,sgn3)
                            call   c8(c4 ,r3,r ,sgn4)
                            call cdg8(cd1,r ,k1,sgn5)
                            call cdg8(cd2,k1,k2,sgn6)
                            call cdg8(cd3,k2,k3,sgn7)
                            call cdg8(cd4,k3,k ,sgn8)
                            j=binary_search8(H8%map,k)
                            if((j.eq.0).and.(HardCoreBoson.eq.0))write(*,'(A,140I6)')" PairX ",distance,c1,c2,c3,c4,cd1,cd2,cd3,cd4,isite,jsite
                            if(j.ne.0)Pcorrfunc(distance,istate) = Pcorrfunc(distance,istate) + peso*gscvec(i)*conjg(gscvec(j))*(sgn1*sgn2*sgn3*sgn4*sgn5*sgn6*sgn7*sgn8)**BosonExp
                         endif
                         !
                         !pair on the x
                         Jcondition = (ib(c1)==1).and.(ib(c2)==1).and.(ib(cd1)==0).and.(ib(cd2)==0)
                         if(Jcondition)then
                            call   c8(c1 ,m ,r1,sgn1)
                            call   c8(c2 ,r1,r ,sgn2)
                            call cdg8(cd1,r ,k1,sgn3)
                            call cdg8(cd2,k1,k ,sgn4)
                            j=binary_search8(H8%map,k)
                            if((j.eq.0).and.(HardCoreBoson.eq.0))write(*,'(A,140I6)')" PairO ",distance,c1,cd1,isite,jsite
                            if(j.ne.0)Ocorrfunc(distance,istate,1) = Ocorrfunc(distance,istate,1) + peso*gscvec(i)*conjg(gscvec(j))*(sgn1*sgn2*sgn3*sgn4)**BosonExp
                         endif
                         !
                         !pair on the y
                         Jcondition = (ib(c1)==1).and.(ib(c4)==1).and.(ib(cd1)==0).and.(ib(cd4)==0).and.Jleg
                         if(Jcondition)then
                            call   c8(c1 ,m ,r1,sgn1)
                            call   c8(c4 ,r1,r ,sgn2)
                            call cdg8(cd1,r ,k1,sgn3)
                            call cdg8(cd4,k1,k ,sgn4)
                            j=binary_search8(H8%map,k)
                            if((j.eq.0).and.(HardCoreBoson.eq.0))write(*,'(A,140I6)')" PairO ",distance,c1,cd1,isite,jsite
                            if(j.ne.0)Ocorrfunc(distance,istate,2) = Ocorrfunc(distance,istate,2) + peso*gscvec(i)*conjg(gscvec(j))*(sgn1*sgn2*sgn3*sgn4)**BosonExp
                         endif
                         !
                      endif
                      !
                   enddo !ineig
                enddo !distance
             enddo !isite
             !
          enddo !states in the block
          !
          if(associated(gscvec))nullify(gscvec)
          call delete_sector8(isector,H8)
          !
          !
          !Single particle correlation function normalization
          NormCorr=0d0
          do distance=1,size(Neigh)
             NormCorr=NormCorr+Ncorrfunc(distance,istate)
          enddo
          if(NormCorr.gt.0d0)Ncorrfunc(:,istate) = Ncorrfunc(:,istate)/NormCorr
          !
          !
          !Quartet correlation function normalization
          NormCorr=0d0
          do distance=1,size(Neigh) !if((Nbath.eq.2).and.(mod(distance,2).eq.0))cycle
             NormCorr=NormCorr+Xcorrfunc(distance,istate)
          enddo
          if(NormCorr.gt.0d0)Xcorrfunc(:,istate) = Xcorrfunc(:,istate)/NormCorr
          !
          !
          unit_ = free_unit()
          open(unit=unit_,file="n_gs"//reg(str(istate))//".DAT",status='unknown',position='rewind',action='write',form='formatted')
          do row=1,Nbath
             write(unit_,'(30(E22.10,1X))') (denslatt(row,col,istate),col=1,Norb)
          enddo
          close(unit_)
          !
          unit_ = free_unit()
          open(unit=unit_,file="Ek_gs"//reg(str(istate))//".DAT",status='unknown',position='rewind',action='write',form='formatted')
          do row=1,Nbath
             write(unit_,'(30(E22.10,1X))') (Eklatt(row,col,istate),col=1,Norb)
          enddo
          close(unit_)
          !
          do distance=0,size(Neigh)
             unit_ = free_unit()
             open(unit=unit_,file="Epot_r"//reg(str(distance))//"_gs"//reg(str(istate))//".DAT",status='unknown',position='rewind',action='write',form='formatted')
             do row=1,Nbath
                 write(unit_,'(30(E22.10,1X))') (Epotlatt(row,col,distance,istate),col=1,Norb)
             enddo
             close(unit_)
          enddo
          !
          unit_ = free_unit()
          open(unit=unit_,file="DensityCorrelators_gs"//reg(str(istate))//".DAT",status='unknown',position='rewind',action='write',form='formatted')
          do distance=1,size(Neigh)
             write(unit_,'(1I3,30(E22.10,1X))') distance,distprint(distance),Ncorrfunc(distance,istate),Xcorrfunc(distance,istate),Gcorrfunc(distance,istate)
          enddo
          close(unit_)
          !
          unit_ = free_unit()
          open(unit=unit_,file="PairCorrelators_gs"//reg(str(istate))//".DAT",status='unknown',position='rewind',action='write',form='formatted')
          do distance=1,size(Neigh)
             !if((Nbath.eq.2).and.(mod(distance,2).eq.0))cycle
             write(unit_,'(1I3,30(E22.10,1X))') distance,distprint(distance),Pcorrfunc(distance,istate),Ocorrfunc(distance,istate,1),Ocorrfunc(distance,istate,2)
          enddo
          close(unit_)
          !
          !
       endif ! master
       !
    enddo ! end loop on degenrate groundstates considered
    !
    !
    if(Mpimaster)then
      do istate=1,numstates
         ntot=0d0
         do isite=1,Norb*Nbath
            row = vec2lat(isite,1)
            col = vec2lat(isite,2)
            ntot = ntot + denslatt(row,col,istate)
         enddo
         write(*,"(A,I4,A,1F10.5)")"#GS= ",istate, " NTOT= ",ntot
      enddo
      !
      NNr=0d0
      do istate=1,numstates
         do distance=1,size(Neigh)
            NNr = NNr + Ncorrfunc(distance,istate)*distprint(distance)/numstates
         enddo
      enddo
      !
      GGr=0d0
      do istate=1,numstates
         do distance=1,size(Neigh)
            GGr = GGr + Gcorrfunc(distance,istate)*distprint(distance)/numstates
         enddo
      enddo
      !
      XXr=0d0
      do istate=1,numstates
         do distance=1,size(Neigh) !-2;if((Nbath.eq.2).and.(mod(distance,2).eq.0))cycle
            XXr = XXr + Xcorrfunc(distance,istate)*distprint(distance)/numstates
         enddo
      enddo
      !
      PPr=0d0
      do istate=1,numstates
         do distance=1,size(Neigh) !-2;if((Nbath.eq.2).and.(mod(distance,2).eq.0))cycle
            PPr = PPr + Pcorrfunc(distance,istate)/numstates
         enddo
      enddo
      !
      OOx=0d0;absOOx=0d0
      OOy=0d0;absOOy=0d0
      do istate=1,numstates
         do distance=1,size(Neigh)
            OOx = OOx + Ocorrfunc(distance,istate,1)/numstates
            absOOx = absOOx + abs(Ocorrfunc(distance,istate,1))/numstates
            OOy = OOy + Ocorrfunc(distance,istate,2)/numstates
            absOOy = absOOy + abs(Ocorrfunc(distance,istate,2))/numstates
         enddo
      enddo
      !
      Ektot=0d0
      do istate=1,numstates
         do row=1,Nbath
            do col=1,Norb
               Ektot=Ektot+Eklatt(row,col,istate)/(Nbath*Norb)
            enddo
         enddo
      enddo
      !
      Eptot=0d0
      do istate=1,numstates
         do distance=0,size(Neigh)
            do row=1,Nbath
               do col=1,Norb
                  Eptot=Eptot+Epotlatt(row,col,distance,istate)/(Nbath*Norb)
               enddo
            enddo
         enddo
      enddo
      !
      !
      unit_ = free_unit()
      open(unit=unit_,file="nlatt.DAT",status='unknown',position='rewind',action='write',form='formatted')
      do row=1,Nbath
        write(unit_,'(30(E22.10,1X))') (denslatt_tot(row,col),col=1,Norb)
      enddo
      close(unit_)
      !
      !
      unit_ = free_unit()
      open(unit=unit_,file="observables_last.ed",status='unknown',position='rewind',action='write',form='formatted')
      !                                                           3                 5     6         7           8   9   10  11  12  13  14
      write(unit_,'(2I3,1F10.5,1I5,30(E22.10,1X))') Nbath,Norb,Thopping,numstates,Ektot,Eptot,Egs/(Nbath*Norb),Xop,NNr,XXr,GGr,PPr,OOx,OOy!,absOOx,absOOy
      close(unit_)
      unit_ = free_unit()
      open(unit_,file="parameters_last.ed")
      write(unit_,"(90F15.9)")xmu,beta,uloc(1),uloc(2),uloc(3),Ust,Jh,Jx,Jp
      close(unit_)
      !
      !
    endif
    !
    call MPI_Barrier(MpiComm,ier)
    !
    deallocate(denslatt,denslatt_tot)
    deallocate(Eklatt,Epotlatt)
    deallocate(Ncorrfunc)
    !
  end subroutine observables_plaquette



end MODULE ED_OBSERVABLES
