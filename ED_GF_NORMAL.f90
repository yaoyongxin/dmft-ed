MODULE ED_GF_NORMAL
  USE ED_GF_SHARED
  implicit none
  private


  public :: build_gf_normal
  public :: build_sigma_normal


contains



  !+------------------------------------------------------------------+
  !                        NORMAL
  !+------------------------------------------------------------------+
  subroutine build_gf_normal()
    integer :: iorb,jorb,ispin,i
    logical :: MaskBool
    !
    !
    !NORMAL: (default)
    do ispin=1,Nspin
       do iorb=1,Norb
          write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_s"//str(ispin)
          if(MPIMASTER)call start_timer
          call lanc_build_gf_normal_c(iorb,ispin)
          if(MPIMASTER)call stop_timer(unit=logfile)
       enddo
    enddo
    !
    !HYBRID:
    if(bath_type/="normal")then
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                !if(hybrid)always T; if(replica)T iff following condition is T
                MaskBool=.true.   
                if(bath_type=="replica")MaskBool=&
                     (dmft_bath%mask(ispin,ispin,iorb,jorb,1)).OR.(dmft_bath%mask(ispin,ispin,iorb,jorb,2))
                if(.not.MaskBool)cycle
                !
                write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
                if(MPIMASTER)call start_timer
                call lanc_build_gf_normal_mix_c(iorb,jorb,ispin)
                if(MPIMASTER)call stop_timer(unit=logfile)
             enddo
          enddo
       enddo
       !Put here off-diagonal manipulation by symmetry:
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                !
                MaskBool=.true.   
                if(bath_type=="replica")MaskBool=&
                     (dmft_bath%mask(ispin,ispin,iorb,jorb,1)).OR.(dmft_bath%mask(ispin,ispin,iorb,jorb,2))
                if(.not.MaskBool)cycle
                !
                impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,ispin,iorb,jorb,:) &
                     - (one-xi)*impGmats(ispin,ispin,iorb,iorb,:) - (one-xi)*impGmats(ispin,ispin,jorb,jorb,:))
                impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,ispin,iorb,jorb,:) &
                     - (one-xi)*impGreal(ispin,ispin,iorb,iorb,:) - (one-xi)*impGreal(ispin,ispin,jorb,jorb,:))
                !>>ACTHUNG: this relation might not be true, it depends on the value of the impHloc_ij
                ! if impHloc_ij is REAL then it is true. if CMPLX hermiticity must be ensured
                impGmats(ispin,ispin,jorb,iorb,:) = impGmats(ispin,ispin,iorb,jorb,:)
                impGreal(ispin,ispin,jorb,iorb,:) = impGreal(ispin,ispin,iorb,jorb,:)
             enddo
          enddo
       enddo
    endif
    !
    ! IF SOMEONE WANTS TO GET BACK THIS I AM GONNA GET ANGRY.
    ! I DEBUGGED THIS 2 or 3 TIMES ALREADY.
    ! IF YOU WANT TO TAKE THE RISK, FINE BUT AT LEAST SET ED_PARA=FALSE TO
    ! AVOID MIS-BEHAVIORS TO OTHERS.
    ! REMARK: THE NSPIN=2 CASE IS MEANT FOR SOLVING MAGNETIC SYSTEMS. IF YOU
    ! SET ED_PARA=TRUE YOU ARE SCREWING UP THE CALCUTIONS OF OTHERS WHO WANT
    ! TO DO MAGNETISM. SO IF YOU KNOW YOU NEED NSPIN=2 AND NO-MAGNETISM
    ! THEN YOU WILL ADJUST YOUR ED_PARA.
    ! I HONESTLY DON'T SEE ANY REASON WHY ED_PARA FLAG SHOULD EXIST ANYWAY.
    ! IF YOU NEED TO SUPPRESS MAGNETISM JUST DO AS THE SMART PEOPLE DO: DON'T BREAK
    ! THE SPIN SYMMETRY AND COPY BATH_UP INTO BATH_DW ONCE BATH_UP IS DETERMINED.
    ! if(ed_para.and.Nspin==2)then
    !    impGmats(1,1,:,:,:) = (impGmats(1,1,:,:,:)+impGmats(2,2,:,:,:))/2.d0
    !    impGmats(2,2,:,:,:) =  impGmats(1,1,:,:,:)
    !    impGreal(1,1,:,:,:) = (impGreal(1,1,:,:,:)+impGreal(2,2,:,:,:))/2.d0
    !    impGreal(2,2,:,:,:) =  impGreal(1,1,:,:,:)
    ! endif
    !
  end subroutine build_gf_normal











  !################################################################
  !################################################################
  !################################################################
  !################################################################








  subroutine lanc_build_gf_normal_c(iorb,ispin)
    complex(8),allocatable           :: vvinit(:),vvloc(:)
    real(8),allocatable              :: alfa_(:),beta_(:)
    integer                          :: iorb,ispin,isite,isector,istate
    integer                          :: idim,jsector
    integer                          :: jdim,vecDim
    integer                          :: ib(Nlevels)
    integer                          :: m,i,j,r,numstates
    real(8)                          :: sgn,norm2,norm0
    complex(8)                       :: cnorm2
    integer                          :: Nitermax,Nlanc
    type(sector_map)                 :: HI,HJ
    !
    isite=impIndex(iorb,ispin)
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
       !ADD ONE PARTICLE:
       jsector = getCDGsector(ispin,isector)
       if(jsector/=0)then 
          !
          jdim  = getdim(jsector)
          !
          !The Op|gs> is worked out by the master only:
          if(MpiMaster)then
             if(ed_verbose==3)write(LOGfile,"(A,2I3)")' add particle:',getnup(jsector),getndw(jsector)
             !
             allocate(vvinit(jdim)) ;  vvinit=zero
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
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,iorb,ispin)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !REMOVE ONE PARTICLE:
       jsector = getCsector(ispin,isector)
       if(jsector/=0)then
          !
          jdim  = getdim(jsector)
          !
          if(MpiMaster)then
             !
             if(ed_verbose==3)write(LOGfile,"(A,2I3)")' del particle:',getnup(jsector),getndw(jsector)
             !
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
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin)
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
  end subroutine lanc_build_gf_normal_c







  !################################################################
  !################################################################
  !################################################################
  !################################################################







  subroutine lanc_build_gf_normal_mix_c(iorb,jorb,ispin)
    integer                          :: iorb,jorb,ispin,isite,jsite,isector,istate
    integer                          :: idim,jsector
    integer                          :: jdim,vecDim
    integer                          :: ib(Nlevels)
    integer                          :: m,i,j,r,numstates
    real(8)                          :: sgn,norm2,norm0
    complex(8)                       :: cnorm2
    complex(8),allocatable           :: vvinit(:),vvloc(:)
    real(8),allocatable              :: alfa_(:),beta_(:)
    integer                          :: Nitermax,Nlanc
    type(sector_map)                 :: HI,HJ
    !
    isite=impIndex(iorb,ispin)  !orbital 1
    jsite=impIndex(jorb,ispin)  !orbital 2
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
       !EVALUATE (c^+_iorb + c^+_jorb)|gs>
       jsector = getCDGsector(ispin,isector)
       if(jsector/=0)then 
          jdim  = getdim(jsector)
          if(MpiMaster)then
             if(ed_verbose==3)write(LOGfile,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
             !
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
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !EVALUATE (c_iorb + c_jorb)|gs>
       jsector = getCsector(ispin,isector)
       if(jsector/=0)then
          jdim   = getdim(jsector)
          if(MpiMaster)then
             if(ed_verbose==3)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
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
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !EVALUATE (c^+_iorb + i*c^+_jorb)|gs>
       jsector = getCDGsector(ispin,isector)
       if(jsector/=0)then 
          jdim  = getdim(jsector)
          if(MpiMaster)then
             if(ed_verbose==3)write(LOGfile,"(A,2I3,I15)")' add particle:',getnup(jsector),getndw(jsector),jdim
             allocate(vvinit(jdim)); vvinit=zero
             !
             call build_sector(jsector,HJ)
             !
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
                if(ib(jsite)==0)then
                   call cdg(jsite,i,r,sgn)
                   j=binary_search(HJ%map,r)
                   vvinit(j) = vvinit(j) + xi*sgn*state_cvec(m)
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
          call add_to_lanczos_gf_normal(-xi*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin)
          call delete_Hv_sector()
          !
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
          if(allocated(vvloc))deallocate(vvloc)
       endif
       !
       !EVALUATE (c_iorb - xi*c_jorb)|gs>
       jsector = getCsector(ispin,isector)
       if(jsector/=0)then
          jdim   = getdim(jsector)
          !
          if(MpiMaster)then
             if(ed_verbose==3)write(LOGfile,"(A,2I3,I15)")' del particle:',getnup(jsector),getndw(jsector),jdim
             allocate(vvinit(jdim));vvinit=zero
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
          call add_to_lanczos_gf_normal(-xi*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin)
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
  end subroutine lanc_build_gf_normal_mix_c








  !################################################################
  !################################################################
  !################################################################
  !################################################################














  subroutine add_to_lanczos_gf_normal(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,ispin)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,iorb,jorb,ispin
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
    !THIS IS HARMLESS BUT WE ARE SLOWLY GOING TO CHANGE EVERYWHERE
    ! FROM TQL2 TO LAPACK EIGH.
    ! diag             = 0.d0
    ! subdiag          = 0.d0
    ! Z                = eye(Nlanc)
    ! diag(1:Nlanc)    = alanc(1:Nlanc)
    ! subdiag(2:Nlanc) = blanc(2:Nlanc)
    ! call tql2(Nlanc,diag,subdiag,Z,ierr)
    !
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       do i=1,Lmats
          iw=xi*wm(i)
          impGmats(ispin,ispin,iorb,jorb,i)=impGmats(ispin,ispin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
       do i=1,Lreal
          iw=dcmplx(wr(i),eps)
          impGreal(ispin,ispin,iorb,jorb,i)=impGreal(ispin,ispin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf_normal







  !################################################################
  !################################################################
  !################################################################
  !################################################################












  subroutine build_sigma_normal
    integer                                           :: i,ispin,iorb
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats,invGmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: invG0real,invGreal
    complex(8),dimension(Norb,Norb)                   :: invGimp
    !
    ! if(.not.allocated(wm))allocate(wm(Lmats))
    ! if(.not.allocated(wr))allocate(wr(Lreal))
    ! wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    ! wr     = linspace(wini,wfin,Lreal)
    !
    invG0mats = zero
    invGmats  = zero
    invG0real = zero
    invGreal  = zero
    !
    !Get G0^-1
    invG0mats(:,:,:,:,:) = invg0_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
    invG0real(:,:,:,:,:) = invg0_bath_real(dcmplx(wr(:),eps),dmft_bath)
    !
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       !
       !Get Gimp^-1
       do ispin=1,Nspin
          do iorb=1,Norb
             invGmats(ispin,ispin,iorb,iorb,:) = one/impGmats(ispin,ispin,iorb,iorb,:)
             invGreal(ispin,ispin,iorb,iorb,:) = one/impGreal(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSreal=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             impSmats(ispin,ispin,iorb,iorb,:) = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
             impSreal(ispin,ispin,iorb,iorb,:) = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !
    case ("hybrid","replica")   !Diagonal in spin only. Full Orbital structure
       !
       !Get Gimp^-1
       do ispin=1,Nspin
          do i=1,Lmats
             invGimp = impGmats(ispin,ispin,:,:,i)
             call inv(invGimp)
             invGmats(ispin,ispin,:,:,i)=invGimp
          enddo
          !
          do i=1,Lreal
             invGimp = impGreal(ispin,ispin,:,:,i)
             call inv(invGimp)
             invGreal(ispin,ispin,:,:,i)=invGimp
          enddo
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSreal=zero
       do ispin=1,Nspin
          impSmats(ispin,ispin,:,:,:) = invG0mats(ispin,ispin,:,:,:) - invGmats(ispin,ispin,:,:,:)
          !
          impSreal(ispin,ispin,:,:,:) = invG0real(ispin,ispin,:,:,:) - invGreal(ispin,ispin,:,:,:)
       enddo
       !
    end select
    !
    !Get G0and:
    impG0mats(:,:,:,:,:) = g0and_bath_mats(dcmplx(0d0,wm(:)),dmft_bath)
    impG0real(:,:,:,:,:) = g0and_bath_real(dcmplx(wr(:),eps),dmft_bath)
    !!
    !
    ! if(allocated(wm))deallocate(wm)
    ! if(allocated(wr))deallocate(wr)
    !
  end subroutine build_sigma_normal


END MODULE ED_GF_NORMAL
