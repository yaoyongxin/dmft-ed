MODULE ED_GF_SUPERC
  USE ED_GF_SHARED
  implicit none
  private

  public :: build_gf_superc
  public :: build_sigma_superc



contains



  !+------------------------------------------------------------------+
  !                        SUPERC
  !+------------------------------------------------------------------+
  subroutine build_gf_superc()
    integer    :: iorb,jorb,ispin,i,isign
    complex(8) :: barGmats(Norb,Lmats),barGreal(Norb,Lreal)
    !
    !
    if(.not.allocated(auxGmats))allocate(auxGmats(3,Lmats))
    if(.not.allocated(auxGreal))allocate(auxGreal(3,Lreal))
    auxgmats=zero
    auxGreal=zero
    barGmats=zero
    barGreal=zero
    !
    ispin=1                       !in this channel Nspin=2 is forbidden. check in ED_AUX_FUNX.
    do iorb=1,Norb
       auxGmats=zero
       auxGreal=zero
       write(LOGfile,"(A)")"Get G&F_l"//str(iorb)//"_s"//str(ispin)
       if(MPIMASTER)call start_timer()
       call lanc_build_gf_superc_c(iorb)
       if(MPIMASTER)call stop_timer(LOGfile)
       !
       impGmats(ispin,ispin,iorb,iorb,:) = auxGmats(1,:) !this is G_{iorb,iorb} = G_{up,up;iorb,iorb}
       impGreal(ispin,ispin,iorb,iorb,:) = auxGreal(1,:)
       barGmats(                 iorb,:) = auxGmats(2,:) !this is \bar{G}_{iorb,iorb} = \bar{G}_{dw,dw;iorb,iorb}
       barGreal(                 iorb,:) = auxGreal(2,:)
       impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(3,:)-auxGmats(1,:)-auxGmats(2,:))
       impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(3,:)-auxGreal(1,:)-auxGreal(2,:))
       !
       ! Comment out this and following lines marked with >anomal to use the more general algorithm
       ! for the evaluation of the anomalous gf
       ! >ANOMAL
       ! impFmats(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGmats(3,:)-(one-xi)*auxGmats(1,:)-(one-xi)*auxGmats(2,:))
       ! impFreal(ispin,ispin,iorb,iorb,:) = 0.5d0*(auxGreal(3,:)-(one-xi)*auxGreal(1,:)-(one-xi)*auxGreal(2,:))
       ! <ANOMAL
    enddo
    !
    !now we add the other mixed/anomalous GF in for the bath_type="hybrid" case
    if(bath_type=='hybrid')then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
             if(MPIMASTER)call start_timer()
             call lanc_build_gf_superc_mix_c(iorb,jorb)
             if(MPIMASTER)call stop_timer()
             impGmats(ispin,ispin,iorb,jorb,:) = auxGmats(3,:)
             impGreal(ispin,ispin,iorb,jorb,:) = auxGreal(3,:)
          enddo
       enddo
       !
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             impFmats(ispin,ispin,iorb,jorb,:) = 0.5d0*( impGmats(ispin,ispin,iorb,jorb,:) - &
                  (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*barGmats(jorb,:) )
             impFreal(ispin,ispin,iorb,jorb,:) = 0.5d0*( impGreal(ispin,ispin,iorb,jorb,:) - &
                  (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*barGreal(jorb,:) )
          enddo
       enddo
    endif
    deallocate(auxGmats,auxGreal)
  end subroutine build_gf_superc










  subroutine lanc_build_gf_superc_c(iorb)
    complex(8),allocatable :: vvinit(:),vvloc(:)
    real(8),allocatable    :: alfa_(:),beta_(:)  
    integer                :: iorb,isector,istate
    integer                :: idim,jsector,vecDim
    integer                :: jdim,isz,jsz
    integer                :: ib(Nlevels)
    integer                :: m,i,j,r,numstates
    real(8)                :: sgn,norm2,norm0
    complex(8)             :: cnorm2
    integer                :: Nitermax,Nlanc
    type(sector_map)       :: HI,HJ
    !
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
       !
       !EVALUATE c^+_{up,iorb}|v> --> Gaux(1) = G_{iorb,iorb}
       jsector = getCDGsector(1,isector)
       if(jsector/=0)then 
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             if(ed_verbose==3)write(LOGfile,"(A23,I3)")'apply c^+_up:',getsz(jsector)
             allocate(vvinit(jdim)) ; vvinit=zero
             !
             call build_sector(jsector,HJ)
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(iorb)==0)then
                   call cdg(iorb,i,r,sgn)
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,1,ichan=1)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !EVALUATE c_{up,iorb}|v> --> Gaux(1) = G_{iorb,iorb}
       jsector = getCsector(1,isector)
       if(jsector/=0)then 
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             if(ed_verbose==3)write(LOGfile,"(A23,I3)")'apply c_up:',getsz(jsector)
             allocate(vvinit(jdim)) ; vvinit=zero
             !
             call build_sector(jsector,HJ)
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(iorb)==1)then
                   call c(iorb,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
             call delete_sector(jsector,HJ)
             !
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,-1,ichan=1)
          !
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !EVALUATE c_{dw,iorb}|v> --> Gaux(2) = barG_{iorb,iorb}
       jsector = getCsector(2,isector)
       if(jsector/=0)then 
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             if(ed_verbose==3)&
                  write(LOGfile,"(A23,I3)")'apply c_dw:',getsz(jsector)
             allocate(vvinit(jdim)) ; vvinit=zero
             call build_sector(jsector,HJ)
             vvinit=zero
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(iorb+Ns)==1)then
                   call c(iorb+Ns,i,r,sgn)
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,1,ichan=2)
          !
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !EVALUATE c^+_{dw,iorb}|v> --> Gaux(2) = barG_{iorb,iorb}
       jsector = getCDGsector(2,isector)
       if(jsector/=0)then 
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             if(ed_verbose==3)&
                  write(LOGfile,"(A23,I3)")'apply c^+_dw:',getsz(jsector)
             allocate(vvinit(jdim)) ; vvinit=zero
             call build_sector(jsector,HJ) !note that here you are doing twice the map building...
             vvinit=zero
             do m=1,idim                     !loop over |gs> components m
                i=HI%map(m)                    !map m to Hilbert space state i
                ib = bdecomp(i,2*Ns)            !i into binary representation
                if(ib(iorb+Ns)==0)then           !if impurity is empty: proceed
                   call cdg(iorb+Ns,i,r,sgn)
                   j=binary_search(HJ%map,r)      !map r back to  jsector
                   vvinit(j) = sgn*state_cvec(m)  !build the cdg_up|gs> state
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,-1,ichan=2)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !EVALUATE [c^+_{up,iorb} + c_{dw,iorb}]|gs> --> A_{iorb,iorb}
       isz = getsz(isector)
       if(isz<Ns)then
          jsz   = isz+1
          jsector = getsector(jsz,1)
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             if(ed_verbose==3)write(LOGfile,"(A23,I3)")'apply c^+_up + c_dw:',getsz(jsector)
             allocate(vvinit(jdim)) ; vvinit=zero
             call build_sector(jsector,HJ)
             vvinit=zero
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(iorb)==0)then
                   call cdg(iorb,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(iorb+Ns)==1)then
                   call c(iorb+Ns,i,r,sgn)
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,1,ichan=3)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !EVALUATE [c_{up,iorb} + c^+_{dw,iorb}]|gs>  --> A_{iorb,iorb}
       isz = getsz(isector)
       if(isz>-Ns)then
          jsz   = isz-1
          jsector = getsector(jsz,1)
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             if(ed_verbose==3)&
                  write(LOGfile,"(A23,I3)")'apply c_up + c^+_dw:',getsz(jsector)
             allocate(vvinit(jdim)) ; vvinit=zero
             call build_sector(jsector,HJ)
             vvinit=zero
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(iorb)==1)then
                   call c(iorb,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(iorb+Ns)==0)then
                   call cdg(iorb+Ns,i,r,sgn)
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,-1,ichan=3)
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
  end subroutine lanc_build_gf_superc_c








  subroutine lanc_build_gf_superc_mix_c(iorb,jorb)
    complex(8),allocatable           :: vvinit(:),vvloc(:)
    real(8),allocatable              :: alfa_(:),beta_(:)  
    integer                          :: iorb,jorb,isector,istate
    integer                          :: idim,jsector,isite
    integer                          :: jdim,isz,jsz,jsite
    integer                          :: ib(Nlevels)
    integer                          :: m,i,j,r,numstates,vecDim
    real(8)                          :: sgn,norm2,norm0
    complex(8)                       :: cnorm2
    integer                          :: Nitermax,Nlanc
    type(sector_map) :: HI,HJ
    !
    isite=impIndex(iorb,1)  !orbital alfa_up
    jsite=impIndex(jorb,2)  !orbital beta_dw
    !
    numstates=state_list%size
    !   
    !
    do istate=1,numstates
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
       !EVALUATE [c^+_{up,iorb} + c_{dw,jorb}]|gs> --> A_{iorb,jorb}
       isz = getsz(isector)
       if(isz<Ns)then
          jsz   = isz+1
          jsector = getsector(jsz,1)
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             if(ed_verbose==3)write(LOGfile,"(A23,I3)")'apply c^+_{up,iorb} + c_{dw,jorb}:',getsz(jsector)
             allocate(vvinit(jdim)) ; vvinit=zero
             call build_sector(jsector,HJ)
             vvinit=zero
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(isite)==0)then
                   call cdg(isite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,1,ichan=3)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !EVALUATE [c_{up,iorb} + c^+_{dw,jorb}]|gs>  --> A_{iorb,jorb}
       isz = getsz(isector)
       if(isz>-Ns)then
          jsz   = isz-1
          jsector = getsector(jsz,1)
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             if(ed_verbose==3)&
                  write(LOGfile,"(A23,I3)")'apply c_{up,iorb} + c^+_{dw,jorb}:',getsz(jsector)
             allocate(vvinit(jdim)) ; vvinit=zero
             call build_sector(jsector,HJ)
             vvinit=zero
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(isite)==1)then
                   call c(isite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
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
          call add_to_lanczos_gf_superc(one*norm2,state_e,alfa_,beta_,-1,ichan=3)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !EVALUATE [c^+_{up,iorb} + xi*c_{dw,jorb}]|gs> --> -xi*B_{iorb,jorb}
       isz = getsz(isector)
       if(isz<Ns)then
          jsz   = isz+1
          jsector = getsector(jsz,1)
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             if(ed_verbose==3)write(LOGfile,"(A23,I3)")'apply c^+_{up,iorb} + xi*c_{dw,horb}:',getsz(jsector)
             allocate(vvinit(jdim)) ; vvinit=zero
             call build_sector(jsector,HJ)
             vvinit=zero
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(isite)==0)then
                   call cdg(isite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(jsite)==1)then
                   call c(jsite,i,r,sgn)
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
          call add_to_lanczos_gf_superc(-xi*norm2,state_e,alfa_,beta_,1,ichan=3)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !EVALUATE [c_{up,iorb} - xi*c^+_{dw,jorb}]|gs> --> -xi*B_{iorb,jorb}
       isz = getsz(isector)
       if(isz>-Ns)then
          jsz   = isz-1
          jsector = getsector(jsz,1)
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             if(ed_verbose==3)&
                  write(LOGfile,"(A23,I3)")'apply c_{up,iorb} - xi*c^+_{dw,jorb}:',getsz(jsector)
             allocate(vvinit(jdim)) ; vvinit=zero
             call build_sector(jsector,HJ)
             vvinit=0.d0
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(isite)==1)then
                   call c(isite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = sgn*state_cvec(m)
                endif
             enddo
             do m=1,idim
                i=HI%map(m)
                ib = bdecomp(i,2*Ns)
                if(ib(jsite)==0)then
                   call cdg(jsite,i,r,sgn)
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
          call add_to_lanczos_gf_superc(-xi*norm2,state_e,alfa_,beta_,-1,ichan=3)
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
    !
    return
  end subroutine lanc_build_gf_superc_mix_c




  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################








  subroutine add_to_lanczos_gf_superc(vnorm2,Ei,alanc,blanc,isign,ichan)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,ichan
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    ! if((finiteT).and.(beta*(Ei-Egs).lt.200))then
    !    pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    ! elseif(.not.finiteT)then
    !   pesoBZ = vnorm2/zeta_function
    !else
    !   pesoBZ=0.d0
    !endif
    !
    pesoBZ = vnorm2/zeta_function
    if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    !
    itype=(3+isign)/2
    diag             = 0.d0
    subdiag          = 0.d0
    Z                = eye(Nlanc)
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call tql2(Nlanc,diag,subdiag,Z,ierr)
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       do i=1,Lmats
          iw=xi*wm(i)
          auxGmats(ichan,i)=auxGmats(ichan,i) + peso/(iw-isign*de)
       enddo
       do i=1,Lreal
          iw=dcmplx(wr(i),eps)
          auxGreal(ichan,i)=auxGreal(ichan,i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf_superc






  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################
  !###################################################################################################







  !+------------------------------------------------------------------+
  !PURPOSE  : Build the Self-energy functions, SUPERC case
  !+------------------------------------------------------------------+
  subroutine build_sigma_superc
    integer                                               :: i,ispin,iorb
    real(8)                                               :: det_mats(Lmats)
    complex(8)                                            :: det_real(Lreal)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)     :: invG0mats,invF0mats,invGmats,invFmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)     :: invG0real,invF0real,invGreal,invFreal
    complex(8),dimension(2*Nspin*Norb,2*Nspin*Norb)       :: invGimp
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    invG0mats = zero
    invF0mats = zero
    invGmats  = zero
    invFmats  = zero
    invG0real = zero
    invF0real = zero
    invGreal  = zero
    invFreal  = zero
    !
    !Get G0^-1,F0^-1
    ispin=1
    invG0mats(ispin,ispin,:,:,:) = invg0_bath_mats(ispin,ispin,dcmplx(0d0,wm(:)),dmft_bath)
    invF0mats(ispin,ispin,:,:,:) = invf0_bath_mats(ispin,ispin,dcmplx(0d0,wm(:)),dmft_bath)
    !
    invG0real(ispin,ispin,:,:,:) = invg0_bath_real(ispin,ispin,dcmplx(wr(:),eps),dmft_bath)
    invF0real(ispin,ispin,:,:,:) = invf0_bath_real(ispin,ispin,dcmplx(wr(:),eps),dmft_bath)
    !
    select case(bath_type)
    case default
       !      
       !Get Gimp^-1
       do iorb=1,Norb
          det_mats  =  abs(impGmats(ispin,ispin,iorb,iorb,:))**2 + (impFmats(ispin,ispin,iorb,iorb,:))**2
          invGmats(ispin,ispin,iorb,iorb,:) = conjg(impGmats(ispin,ispin,iorb,iorb,:))/det_mats
          invFmats(ispin,ispin,iorb,iorb,:) = impFmats(ispin,ispin,iorb,iorb,:)/det_mats
          !
          det_real  = -impGreal(ispin,ispin,iorb,iorb,:)*conjg(impGreal(ispin,ispin,iorb,iorb,Lreal:1:-1)) - impFreal(ispin,ispin,iorb,iorb,:)**2
          invGreal(ispin,ispin,iorb,iorb,:) =  -conjg(impGreal(ispin,ispin,iorb,iorb,Lreal:1:-1))/det_real(:)
          invFreal(ispin,ispin,iorb,iorb,:) =  -impFreal(ispin,ispin,iorb,iorb,:)/det_real(:)
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSAmats=zero
       impSreal=zero
       impSAreal=zero
       do iorb=1,Norb
          impSmats(ispin,ispin,iorb,iorb,:)  = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
          impSAmats(ispin,ispin,iorb,iorb,:) = invF0mats(ispin,ispin,iorb,iorb,:) - invFmats(ispin,ispin,iorb,iorb,:)
          !
          impSreal(ispin,ispin,iorb,iorb,:)  = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
          impSAreal(ispin,ispin,iorb,iorb,:) = invF0real(ispin,ispin,iorb,iorb,:) - invFreal(ispin,ispin,iorb,iorb,:)
       enddo
       !
    case ("hybrid")
       !
       !Get Gimp^-1
       do i=1,Lmats
          invGimp=zero
          invGimp(1:Norb,1:Norb)               = impGmats(ispin,ispin,:,:,i)
          invGimp(1:Norb,Norb+1:2*Norb)        = impFmats(ispin,ispin,:,:,i)
          invGimp(Norb+1:2*Norb,1:Norb)        = impFmats(ispin,ispin,:,:,i)
          invGimp(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(impGmats(ispin,ispin,:,:,i))
          call inv(invGimp)
          invGmats(ispin,ispin,:,:,i) = invGimp(1:Norb,1:Norb)
          invFmats(ispin,ispin,:,:,i) = invGimp(1:Norb,Norb+1:2*Norb)
       enddo
       do i=1,Lreal
          invGimp=zero
          invGimp(1:Norb,1:Norb)               = impGreal(ispin,ispin,:,:,i)
          invGimp(1:Norb,Norb+1:2*Norb)        = impFreal(ispin,ispin,:,:,i)
          invGimp(Norb+1:2*Norb,1:Norb)        = impFreal(ispin,ispin,:,:,i)
          invGimp(Norb+1:2*Norb,Norb+1:2*Norb) =-conjg(impGreal(ispin,ispin,:,:,Lreal-i+1))
          call inv(invGimp)
          invGreal(ispin,ispin,:,:,i) =  invGimp(1:Norb,1:Norb)
          invFreal(ispin,ispin,:,:,i) =  invGimp(1:Norb,Norb+1:2*Norb)
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSAmats=zero
       impSreal=zero
       impSAreal=zero
       !
       impSmats(ispin,ispin,:,:,:)  = invG0mats(ispin,ispin,:,:,:) - invGmats(ispin,ispin,:,:,:)
       impSAmats(ispin,ispin,:,:,:) = invF0mats(ispin,ispin,:,:,:) - invFmats(ispin,ispin,:,:,:)
       !
       impSreal(ispin,ispin,:,:,:)  = invG0real(ispin,ispin,:,:,:) - invGreal(ispin,ispin,:,:,:)
       impSAreal(ispin,ispin,:,:,:) = invF0real(ispin,ispin,:,:,:) - invFreal(ispin,ispin,:,:,:)
       !
    end select
    !
    !Get G0and:
    impG0mats(ispin,ispin,:,:,:) = g0and_bath_mats(ispin,ispin,dcmplx(0d0,wm(:)),dmft_bath)
    impF0mats(ispin,ispin,:,:,:) = f0and_bath_mats(ispin,ispin,dcmplx(0d0,wm(:)),dmft_bath)
    !
    impG0real(ispin,ispin,:,:,:) = g0and_bath_real(ispin,ispin,dcmplx(wr(:),eps),dmft_bath)
    impF0real(ispin,ispin,:,:,:) = f0and_bath_real(ispin,ispin,dcmplx(wr(:),eps),dmft_bath)
    !!
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  end subroutine build_sigma_superc





END MODULE ED_GF_SUPERC
