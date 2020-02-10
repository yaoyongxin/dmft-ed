MODULE ED_GF_NONSU2
  USE ED_GF_SHARED
  implicit none
  private


  public :: build_gf_nonsu2
  public :: build_sigma_nonsu2


contains



  !+------------------------------------------------------------------+
  !                            NONSU2
  !+------------------------------------------------------------------+
  subroutine build_gf_nonsu2()
    integer :: iorb,jorb,ispin,jspin,i,io,jo,ifreq
    integer :: isect0,numstates
    !
    if(.not.allocated(impGmats))stop "build_gf_nonsu2: Gmats not allocated"
    if(.not.allocated(impGreal))stop "build_gf_nonsu2: Greal not allocated"
    impGmats=zero
    impGreal=zero
    !
    select case(bath_type)
    case("normal")
       !
       !Here we evaluate the same orbital, same spin GF: G_{aa}^{ss}(z)
       do ispin=1,Nspin
          do iorb=1,Norb
             write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//str(ispin)
             if(MPIMASTER)call start_timer()
             call lanc_build_gf_nonsu2_diagOrb_diagSpin_c(iorb,ispin)
             if(MPIMASTER)call stop_timer(LOGfile)
          enddo
       enddo
       !same orbital, different spin GF: G_{aa}^{ss'}(z)
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.ne.jspin).and.(iorb.eq.jorb)) then
                      write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                      if(MPIMASTER)call start_timer()
                      call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                      if(MPIMASTER)call stop_timer(LOGfile)
                   endif
                enddo
             enddo
          enddo
       enddo
       !Here we put the symmetry manipulation
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.ne.jspin).and.(iorb.eq.jorb)) then
                      !
                      impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                      !
                      impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                      !
                   endif
                enddo
             enddo
          enddo
       enddo
       !
    case("hybrid")
       !
       !Here we evaluate the same orbital, same spin GF: G_{aa}^{ss}(z)
       do ispin=1,Nspin
          do iorb=1,Norb
             write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//str(ispin)
             if(MPIMASTER)call start_timer()
             call lanc_build_gf_nonsu2_diagOrb_diagSpin_c(iorb,ispin)
             if(MPIMASTER)call stop_timer(LOGfile)
          enddo
       enddo
       !
       !same orbital, different spin GF: G_{aa}^{ss'}(z)
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.ne.jspin).and.(iorb.eq.jorb)) then
                      write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                      if(MPIMASTER)call start_timer()
                      call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                      if(MPIMASTER)call stop_timer(LOGfile)
                   endif
                enddo
             enddo
          enddo
       enddo
       !Here we put the symmetry manipulation
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.ne.jspin).and.(iorb.eq.jorb)) then
                      !
                      impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                      !
                      impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                      !
                   endif
                enddo
             enddo
          enddo
       enddo
       !
       !Here we evaluate the different orbital, same spin GF: G_{ab}^{ss}(z)
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.eq.jspin).and.(iorb.ne.jorb)) then
                      write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                      if(MPIMASTER)call start_timer()
                      call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                      if(MPIMASTER)call stop_timer(LOGfile)
                   endif
                enddo
             enddo
          enddo
       enddo
       !Here we put the symmetry manipulation
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.eq.jspin).and.(iorb.ne.jorb)) then
                      !
                      impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                      !
                      impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                      !
                   endif
                enddo
             enddo
          enddo
       enddo
       !
       !Here we evaluate the different orbital, different spin GF: G_{ab}^{ss'}(z)
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.ne.jspin).and.(iorb.ne.jorb)) then
                      write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                      if(MPIMASTER)call start_timer()
                      call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                      if(MPIMASTER)call stop_timer(LOGfile)
                   endif
                enddo
             enddo
          enddo
       enddo
       !Here we put the symmetry manipulation
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.ne.jspin).and.(iorb.ne.jorb)) then
                      !
                      impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                      !
                      impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                      !
                   endif
                enddo
             enddo
          enddo
       enddo
       !
    case("replica")
       !
       !Here we evaluate the same orbital, same spin GF: G_{aa}^{ss}(z)
       do ispin=1,Nspin
          do iorb=1,Norb
             write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//str(ispin)
             if(MPIMASTER)call start_timer()
             call lanc_build_gf_nonsu2_diagOrb_diagSpin_c(iorb,ispin)
             if(MPIMASTER)call stop_timer(LOGfile)
          enddo
       enddo
       !
       !same orbital, different spin GF: G_{aa}^{ss'}(z)
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.ne.jspin).and.(iorb.eq.jorb)) then
                      if((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv. .false.)&
                           .and.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv. .false.))cycle
                      write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                      if(MPIMASTER)call start_timer()
                      call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                      if(MPIMASTER)call stop_timer(LOGfile)
                   endif
                enddo
             enddo
          enddo
       enddo
       !Here we put the symmetry manipulation
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.ne.jspin).and.(iorb.eq.jorb)) then
                      if((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv. .false.)&
                           .and.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv. .false.))cycle
                      !
                      impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                      !
                      impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                      !
                   endif
                enddo
             enddo
          enddo
       enddo
       !
       !Here we evaluate the different orbital, same spin GF: G_{ab}^{ss}(z)
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.eq.jspin).and.(iorb.ne.jorb)) then
                      if((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv. .false.)&
                           .and.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv. .false.))cycle
                      write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                      if(MPIMASTER)call start_timer()
                      call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                      if(MPIMASTER)call stop_timer(LOGfile)
                   endif
                enddo
             enddo
          enddo
       enddo
       !Here we put the symmetry manipulation
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.eq.jspin).and.(iorb.ne.jorb)) then
                      if((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv. .false.)&
                           .and.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv. .false.))cycle
                      !
                      impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                      !
                      impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                      !
                   endif
                enddo
             enddo
          enddo
       enddo
       !
       !Here we evaluate the different orbital, different spin GF: G_{ab}^{ss'}(z)
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.ne.jspin).and.(iorb.ne.jorb)) then
                      if((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv. .false.)&
                           .and.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv. .false.))cycle
                      write(LOGfile,"(A)")"Get G_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//str(jspin)
                      if(MPIMASTER)call start_timer()
                      call lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
                      if(MPIMASTER)call stop_timer(LOGfile)
                   endif
                enddo
             enddo
          enddo
       enddo
       !Here we put the symmetry manipulation
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   if((ispin.ne.jspin).and.(iorb.ne.jorb)) then
                      if((dmft_bath%mask(ispin,jspin,iorb,jorb,1).eqv. .false.)&
                           .and.(dmft_bath%mask(ispin,jspin,iorb,jorb,2).eqv. .false.))cycle
                      !
                      impGmats(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGmats(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGmats(jspin,jspin,jorb,jorb,:))
                      !
                      impGreal(ispin,jspin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,jspin,iorb,jorb,:) &
                           - (one+xi)*impGreal(ispin,ispin,iorb,iorb,:) &
                           - (one+xi)*impGreal(jspin,jspin,jorb,jorb,:))
                      !
                   endif
                enddo
             enddo
          enddo
       enddo
       !
       !if(ed_para)then
       !   call SOC_jz_symmetrize(impGmats,dmft_bath%mask)
       !   call SOC_jz_symmetrize(impGreal,dmft_bath%mask)
       !endif
       !
    end select
  end subroutine build_gf_nonsu2








  !PURPOSE: Evaluate the same orbital IORB, same spin ISPIN impurity GF.
  subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin_c(iorb,ispin)
    complex(8),allocatable           :: vvinit(:),vvloc(:)
    real(8),allocatable              :: alfa_(:),beta_(:)
    integer                          :: iorb,ispin,isite,isector,istate
    integer                          :: idim,jsector
    integer                          :: jdim,vecDim
    integer                          :: ib(Nlevels)
    integer                          :: m,i,j,r
    real(8)                          :: sgn,norm2
    complex(8)                       :: cnorm2
    integer                          :: Nitermax,Nlanc
    type(sector_map) :: HI,HJ
    !
    isite=impIndex(iorb,ispin)
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       idim  = getdim(isector)
       call build_sector(isector,HI)
       !
       !ADD ONE PARTICLE with IORB,ISPIN:
       !
       if(Jz_basis)then
          jsector = getCDGsector_Jz(iorb,ispin,isector)
       else
          jsector = getCDGsector(ispin,isector)
       endif
       !
       if(getN(isector)/=Nlevels.and.jsector>=0)then
          if(Jz_basis)then
             if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(3(A,1F5.1,1X))")&
                  "add Jz:",Lzdiag(iorb)+Szdiag(ispin)/2.,&
                  "  from:",gettwoJz(isector)/2.,"  to:",gettwoJz(jsector)/2.
          else
             if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(A15,I3)")' add particle:',getn(jsector)
          endif
          !
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             allocate(vvinit(jdim)) ; vvinit=zero
             !
             call build_sector(jsector,HJ)
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(isite)==0)then
                   call cdg(isite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
             call delete_sector(jsector,HJ)
             norm2=dot_product(vvinit,vvinit)
             vvinit=vvinit/sqrt(norm2)
          endif
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          !
          call build_Hv_sector(jsector)
#ifdef _MPI
          if(MpiStatus)then
             call Bcast_MPI(MpiComm,norm2)
             vecDim = vecDim_Hv_sector(jsector)
             allocate(vvloc(vecDim))
             call scatter_vector_MPI(MpiComm,vvinit,vvloc)
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvloc,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          cnorm2=one*norm2
          call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,iorb,ispin,ispin)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !REMOVE ONE PARTICLE with ISPIN:
       !
       if(Jz_basis)then
          jsector = getCsector_Jz(iorb,ispin,isector)
       else
          jsector = getCsector(ispin,isector)
       endif
       !
       if(getN(isector)/=0.and.jsector>=0)then
          if(Jz_basis)then
             if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(3(A,1F5.1,1X))")&
                  "del Jz:",Lzdiag(iorb)+Szdiag(ispin)/2.,&
                  "  from:",gettwoJz(isector)/2.,"  to:",gettwoJz(jsector)/2.
          else
             if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(A15,I3)")' del particle:',getn(jsector)
          endif
          !
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             allocate(vvinit(jdim)) ; vvinit=zero
             !
             call build_sector(jsector,HJ)
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(isite)==1)then
                   call c(isite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
             call delete_sector(jsector,HJ)
             norm2=dot_product(vvinit,vvinit)
             vvinit=vvinit/sqrt(norm2)
          endif
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          !
          call build_Hv_sector(jsector)
#ifdef _MPI
          if(MpiStatus)then
             call Bcast_MPI(MpiComm,norm2)
             vecDim = vecDim_Hv_sector(jsector)
             allocate(vvloc(vecDim))
             call scatter_vector_MPI(MpiComm,vvinit,vvloc)
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvloc,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          cnorm2=one*norm2
          call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin,ispin)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       nullify(state_cvec)
       call delete_sector(isector,HI)
       !
    enddo
    return
  end subroutine lanc_build_gf_nonsu2_diagOrb_diagSpin_c




  !PURPOSE: Evaluate the same different orbital IORB,JORB, different spin ISPIN,JSPIN impurity GF.
  subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin_c(iorb,jorb,ispin,jspin)
    integer                          :: iorb,jorb,ispin,jspin,isite,jsite,isector,istate
    integer                          :: idim,jsector,vecDim
    integer                          :: jdim,jdim_old
    integer                          :: ib(Nlevels)
    integer                          :: m,i,j,r
    real(8)                          :: sgn,norm2
    complex(8)                       :: cnorm2
    complex(8),allocatable           :: vvinit(:),vvloc(:)
    real(8),allocatable              :: alfa_(:),beta_(:)
    integer                          :: Nitermax,Nlanc
    type(sector_map)                 :: HI,HJ
    !
    isite=impIndex(iorb,ispin)
    jsite=impIndex(jorb,jspin)
    !
    !
    do istate=1,state_list%size
       isector     =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       idim  = getdim(isector)
       call build_sector(isector,HI)
       !
       !
       !APPLY (c^+_{jorb,jspin} + c^+_{iorb,ispin})|gs>
       if(Jz_basis)then
          jsector = getCDGsector_Jz(iorb,ispin,isector)
       else
          jsector = getCDGsector(ispin,isector)
       endif
       !
       if(getN(isector)/=Nlevels.and.jsector>=0)then
          if(Jz_basis)then
             if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(3(A,1F5.1,1X))")&
                  "  add Jz:",Lzdiag(iorb)+Szdiag(ispin)/2.,&
                  "  from:",gettwoJz(isector)/2.,"  to:",gettwoJz(jsector)/2.
          else
             if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(A15,I3)")' add particle:',getn(jsector)
          endif
          !
          jdim_old  = getdim(jsector)
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             allocate(vvinit(jdim)) ; vvinit=zero
             !
             call build_sector(jsector,HJ)
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(isite)==0)then
                   call cdg(isite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
             if(Jz_basis)then
                jsector = getCDGsector_Jz(jorb,jspin,isector)
                if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(3(A,1F5.1,1X))")&
                     "  add Jz:",Lzdiag(jorb)+Szdiag(jspin)/2.,&
                     "  from:",gettwoJz(isector)/2.,"  to:",gettwoJz(jsector)/2.
                jdim  = getdim(jsector)
                if(jdim/=jdim_old)stop "lanczos builgf dimensional error"
                call delete_sector(jsector,HJ)
                call build_sector(jsector,HJ)
             endif
             !
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(jsite)==0)then
                   call cdg(jsite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = vvinit(j) + sgn*state_cvec(m)
                endif
             enddo
             call delete_sector(jsector,HJ)
             norm2=dot_product(vvinit,vvinit)
             vvinit=vvinit/sqrt(norm2)
          endif
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          !
          call build_Hv_sector(jsector)
#ifdef _MPI
          if(MpiStatus)then
             call Bcast_MPI(MpiComm,norm2)
             vecDim = vecDim_Hv_sector(jsector)
             allocate(vvloc(vecDim))
             call scatter_vector_MPI(MpiComm,vvinit,vvloc)
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvloc,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          cnorm2=one*norm2
          call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,jspin)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !APPLY (c_{jorb,jspin} + c_{iorb,ispin})|gs>
       if(Jz_basis)then
          jsector = getCsector_Jz(iorb,ispin,isector)
       else
          jsector = getCsector(ispin,isector)
       endif
       !
       if(getN(isector)/=0.and.jsector>=0)then
          if(Jz_basis)then
             if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(3(A,1F5.1,1X))")&
                  "  del Jz:",Lzdiag(iorb)+Szdiag(ispin)/2.,&
                  "  from:",gettwoJz(isector)/2.,"  to:",gettwoJz(jsector)/2.
          else
             if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(A15,I3)")' del particle:',getn(jsector)
          endif
          !
          jdim_old  = getdim(jsector)
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             allocate(vvinit(jdim)) ; vvinit=zero
             !
             call build_sector(jsector,HJ)
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(isite)==1)then
                   call c(isite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
             !
             if(Jz_basis)then
                jsector = getCsector_Jz(jorb,jspin,isector)
                if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(3(A,1F5.1,1X))")&
                     "  del Jz:",Lzdiag(jorb)+Szdiag(jspin)/2.,&
                     "  from:",gettwoJz(isector)/2.,"  to:",gettwoJz(jsector)/2.
                jdim  = getdim(jsector)
                if(jdim/=jdim_old)stop "lanczos builgf dimensional error"
                call delete_sector(jsector,HJ)
                call build_sector(jsector,HJ)
             endif
             !
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(jsite)==1)then
                   call c(jsite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = vvinit(j) + sgn*state_cvec(m)
                endif
             enddo
             call delete_sector(jsector,HJ)
             norm2=dot_product(vvinit,vvinit)
             vvinit=vvinit/sqrt(norm2)
          endif
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          !
          call build_Hv_sector(jsector)
#ifdef _MPI
          if(MpiStatus)then
             call Bcast_MPI(MpiComm,norm2)
             vecDim = vecDim_Hv_sector(jsector)
             allocate(vvloc(vecDim))
             call scatter_vector_MPI(MpiComm,vvinit,vvloc)
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvloc,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          cnorm2=one*norm2
          call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,jspin)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !
       !APPLY (+i*c^+_{jorb,jspin} + c^+_{iorb,ispin})|gs>
       if(Jz_basis)then
          jsector = getCDGsector_Jz(iorb,ispin,isector)
       else
          jsector = getCDGsector(ispin,isector)
       endif
       !
       if(getN(isector)/=Nlevels.and.jsector>=0)then
          if(Jz_basis)then
             if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(3(A,1F5.1,1X))")&
                  "  add Jz:",Lzdiag(iorb)+Szdiag(ispin)/2.,&
                  "  from:",gettwoJz(isector)/2.,"  to:",gettwoJz(jsector)/2.
          else
             if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(A15,I3)")' add particle:',getn(jsector)
          endif
          !
          jdim_old  = getdim(jsector)
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             allocate(vvinit(jdim)) ; vvinit=zero
             !
             call build_sector(jsector,HJ)
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(isite)==0)then
                   call cdg(isite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
             !
             if(Jz_basis)then
                jsector = getCDGsector_Jz(jorb,jspin,isector)
                if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(3(A,1F5.1,1X))")&
                     "  add Jz:",Lzdiag(jorb)+Szdiag(jspin)/2.,&
                     "  from:",gettwoJz(isector)/2.,"  to:",gettwoJz(jsector)/2.
                jdim  = getdim(jsector)
                if(jdim/=jdim_old)stop "lanczos builgf dimensional error"
                call delete_sector(jsector,HJ)
                call build_sector(jsector,HJ)
             endif
             !
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(jsite)==0)then
                   call cdg(jsite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = vvinit(j) + xi*sgn*state_cvec(m)
                endif
             enddo
             call delete_sector(jsector,HJ)
             norm2=dot_product(vvinit,vvinit)
             vvinit=vvinit/sqrt(norm2)
          endif
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          !
          call build_Hv_sector(jsector)
#ifdef _MPI
          if(MpiStatus)then
             call Bcast_MPI(MpiComm,norm2)
             vecDim = vecDim_Hv_sector(jsector)
             allocate(vvloc(vecDim))
             call scatter_vector_MPI(MpiComm,vvinit,vvloc)
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvloc,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          cnorm2=1*xi*norm2
          call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,1,iorb,jorb,ispin,jspin)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !APPLY (-xi*c_{jorb,jspin} + c_{iorb,ispin})|gs>
       if(Jz_basis)then
          jsector = getCsector_Jz(iorb,ispin,isector)
       else
          jsector = getCsector(ispin,isector)
       endif
       !
       if(getN(isector)/=0.and.jsector>=0)then
          if(Jz_basis)then
             if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(3(A,1F5.1,1X))")&
                  "  del Jz:",Lzdiag(iorb)+Szdiag(ispin)/2.,&
                  "  from:",gettwoJz(isector)/2.,"  to:",gettwoJz(jsector)/2.
          else
             if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(A15,I3)")' del particle:',getn(jsector)
          endif
          !
          jdim_old  = getdim(jsector)
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             allocate(vvinit(jdim)) ; vvinit=zero
             !
             call build_sector(jsector,HJ)
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(isite)==1)then
                   call c(isite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
             !
             if(Jz_basis)then
                jsector = getCsector_Jz(jorb,jspin,isector)
                if(ed_verbose==3.and.MPIMASTER)write(LOGfile,"(3(A,1F5.1,1X))")&
                     "  del Jz:",Lzdiag(jorb)+Szdiag(jspin)/2.,&
                     "  from:",gettwoJz(isector)/2.,"  to:",gettwoJz(jsector)/2.
                jdim  = getdim(jsector)
                if(jdim/=jdim_old)stop "lanczos builgf dimensional error"
                call delete_sector(jsector,HJ)
                call build_sector(jsector,HJ)
             endif
             !
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(jsite)==1)then
                   call c(jsite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = vvinit(j) - xi*sgn*state_cvec(m)
                endif
             enddo
             call delete_sector(jsector,HJ)
             norm2=dot_product(vvinit,vvinit)
             vvinit=vvinit/sqrt(norm2)
          endif
          !
          nlanc=min(jdim,lanc_nGFiter)
          allocate(alfa_(nlanc),beta_(nlanc))
          !
          call build_Hv_sector(jsector)
#ifdef _MPI
          if(MpiStatus)then
             call Bcast_MPI(MpiComm,norm2)
             vecDim = vecDim_Hv_sector(jsector)
             allocate(vvloc(vecDim))
             call scatter_vector_MPI(MpiComm,vvinit,vvloc)
             call sp_lanc_tridiag(MpiComm,spHtimesV_cc,vvloc,alfa_,beta_)
          else
             call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
          endif
#else
          call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
#endif
          cnorm2=1*xi*norm2
          call add_to_lanczos_gf_nonsu2(cnorm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin,jspin)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       nullify(state_cvec)
       call delete_sector(isector,HI)
       !
    enddo
    return
  end subroutine lanc_build_gf_nonsu2_mixOrb_mixSpin_c




  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################






  subroutine add_to_lanczos_gf_nonsu2(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,ispin,jspin)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc
    integer                                    :: isign,iorb,jorb,ispin,jspin
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    if((finiteT).and.(beta*(Ei-Egs).lt.200))then
       pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    elseif(.not.finiteT)then
       pesoBZ = vnorm2/zeta_function
    else
       pesoBZ=0.d0
    endif
    !
    !pesoBZ = vnorm2/zeta_function
    !if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    !
    ! itype=(3+isign)/2
    diag             = 0.d0
    subdiag          = 0.d0
    Z                = eye(Nlanc)
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       do i=1,Lmats
          iw=xi*wm(i)
          impGmats(ispin,jspin,iorb,jorb,i)=impGmats(ispin,jspin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
       do i=1,Lreal
          iw=dcmplx(wr(i),eps)
          impGreal(ispin,jspin,iorb,jorb,i)=impGreal(ispin,jspin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf_nonsu2








  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################









  !+------------------------------------------------------------------+
  !PURPOSE  : Build the Self-energy functions, NONSU2 case
  !+------------------------------------------------------------------+
  subroutine build_sigma_nonsu2
    integer                                           :: i,j,isign,unit(7),iorb,jorb,ispin,jspin,io,jo
    complex(8)                                        :: fg0
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats,invimpG0mats,impDeltamats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: invG0real,invimpG0real,impDeltareal
    complex(8),dimension(Nspin*Norb,Nspin*Norb)       :: invGimp,Foo
    character(len=20)                                 :: suffix
    !
    !
    ! if(.not.allocated(wm))allocate(wm(Lmats))
    ! if(.not.allocated(wr))allocate(wr(Lreal))
    ! wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    ! wr     = linspace(wini,wfin,Lreal)
    !
    impG0mats=zero
    impG0real=zero
    invG0mats = zero
    invG0real = zero
    !
    !Get G0^-1
    invG0mats(:,:,:,:,:)=invg0_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
    invG0real(:,:,:,:,:)=invg0_bath_real(dcmplx(wr(:),eps),dmft_bath)
    !Get impDelta_anderson
    impDeltamats(:,:,:,:,:)=delta_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
    impDeltareal(:,:,:,:,:)=delta_bath_real(dcmplx(wr(:),eps),dmft_bath)
    !Get inverse functions
    invimpG0mats=invG0mats
    invimpG0real=invG0real
    !
    select case(bath_type)
       !
    case ("normal")
       !
       !Get Gimp^-1 - Matsubara freq.
       do i=1,Lmats
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      if (iorb.eq.jorb) then
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         invGimp(io,jo) = impGmats(ispin,jspin,iorb,jorb,i)
                      endif
                   enddo
                enddo
             enddo
          enddo
          call inv(invGimp)!<--- get [G_{imp}]^-1
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      if (iorb.eq.jorb) then
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         impSmats(ispin,jspin,iorb,jorb,i) = invG0mats(ispin,jspin,iorb,jorb,i) - invGimp(io,jo) !<-- calG0_imp^-1 - Gimp^-1
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
       !Get Gimp^-1 - Real freq.
       do i=1,Lreal
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      if (iorb.eq.jorb) then
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         invGimp(io,jo) = impGreal(ispin,jspin,iorb,jorb,i)
                      endif
                   enddo
                enddo
             enddo
          enddo
          call inv(invGimp)!<--- get [G_{imp}]^-1
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      if (iorb.eq.jorb) then
                         io = iorb + (ispin-1)*Norb
                         jo = jorb + (jspin-1)*Norb
                         impSreal(ispin,jspin,iorb,jorb,i) = invG0real(ispin,jspin,iorb,jorb,i) - invGimp(io,jo) !<-- calG0_imp^-1 - Gimp^-1
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    case ("hybrid","replica")
       !
       !Get Gimp^-1 - Matsubara freq.
       do i=1,Lmats
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      invGimp(io,jo) = impGmats(ispin,jspin,iorb,jorb,i)
                   enddo
                enddo
             enddo
          enddo
          call inv(invGimp)!<--- get [G_{imp}]^-1
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      impSmats(ispin,jspin,iorb,jorb,i) = invG0mats(ispin,jspin,iorb,jorb,i) - invGimp(io,jo) !<-- calG0_imp^-1 - Gimp^-1
                   enddo
                enddo
             enddo
          enddo
       enddo
       !Get Gimp^-1 - Real freq.
       do i=1,Lreal
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      invGimp(io,jo) = impGreal(ispin,jspin,iorb,jorb,i)
                   enddo
                enddo
             enddo
          enddo
          call inv(invGimp)!<--- get [G_{imp}]^-1
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb
                      jo = jorb + (jspin-1)*Norb
                      impSreal(ispin,jspin,iorb,jorb,i) = invG0real(ispin,jspin,iorb,jorb,i) - invGimp(io,jo) !<-- calG0_imp^-1 - Gimp^-1
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    end select
    !
    !Get G0and:
    impG0mats(:,:,:,:,:) = g0and_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
    impG0real(:,:,:,:,:) = g0and_bath_real(dcmplx(wr(:),eps),dmft_bath)
    !
    !
    !
    ! if(allocated(wm))deallocate(wm)
    ! if(allocated(wr))deallocate(wr)
    !
  end subroutine build_sigma_nonsu2





END MODULE ED_GF_NONSU2
