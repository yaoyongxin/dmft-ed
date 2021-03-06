!+-------------------------------------------------------------+
!PURPOSE  : Chi^2 interface for Hybrid bath normal phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_hybrid_normal(fg,bath_,ispin)
  complex(8),dimension(:,:,:)          :: fg
  real(8),dimension(:,:),intent(inout) :: bath_
  integer                              :: ispin
  real(8),dimension((1+Norb)*Nbath)    :: a
  integer                              :: iter,stride,ifirst,ilast,i,j,corb,l
  integer                              :: iorb,jorb
  real(8)                              :: chi
  complex(8),dimension(size(fg,3))     ::  g0
  logical                              :: check
  type(effective_bath)                 :: dmft_bath
  complex(8)                           :: fgand
  real(8)                              :: w
  character(len=20)                    :: suffix
  integer                              :: unit
  if(size(fg,1)/=Norb)stop "CHI2FIT: wrong dimension 1 in chi2_input"
  if(size(fg,2)/=Norb)stop "CHI2FIT: wrong dimension 2 in chi2_input"
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_irred: wrong bath dimensions"
  allocate(getIorb(Norb*(Norb+1)/2),getJorb(Norb*(Norb+1)/2))
  corb=0
  do iorb=1,Norb
     do jorb=iorb,Norb
        corb=corb+1
        getIorb(corb)=iorb
        getJorb(corb)=jorb
     enddo
  enddo
  totNorb=corb
  if(totNorb/=(Norb*(Norb+1)/2))stop "CHI2FIT: Error counting the orbitals"
  !
  Ldelta = Lfit
  if(Ldelta>size(fg,3))Ldelta=size(fg,3)
  allocate(Fdelta(totNorb,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  do i=1,totNorb
     Fdelta(i,1:Ldelta) = fg(getIorb(i),getJorb(i),1:Ldelta)
  enddo
  forall(i=1:Ldelta)Xdelta(i)=pi/beta*real(2*i-1,8)
  select case(Cg_weight)
  case default
     Wdelta=dble(Ldelta)
  case(1)
     Wdelta=1.d0
  case(2)
     Wdelta=(/(real(i,8),i=1,Ldelta)/)
  case(3)
     Wdelta=Xdelta
  end select
  !
  call allocate_bath(dmft_bath)
  Spin_indx=ispin
  a(:) = bath_(ispin,:)
  if(cg_method==0)then
     if(cg_scheme=='weiss')then
        call fmin_cg(a,chi2_weiss_hybrd,iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
     else
        call fmin_cg(a,chi2_delta_hybrd,dchi2_delta_hybrd,iter,chi,itmax=cg_niter,ftol=cg_Ftol,istop=cg_stop,eps=cg_eps)
     endif
  elseif(cg_method==1)then
     if(cg_scheme=='weiss')then
        call fmin_cgminimize(a,chi2_weiss_hybrd,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     else
        call fmin_cgminimize(a,chi2_delta_hybrd,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     endif
  elseif(cg_method==2)then
     if(cg_scheme=='weiss')then
        call fmin_cgplus(a,chi2_weiss_hybrd,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     else
        call fmin_cgplus(a,chi2_delta_hybrd,dchi2_delta_hybrd,iter,chi,itmax=cg_niter,ftol=cg_Ftol)
     endif
  else
     stop "ED_CHI2FIT: error cg_method > 2"
  end if
  bath_(ispin,:) = a(:)
  call set_bath(bath_,dmft_bath)
  if(ed_verbose<5)write(LOGfile,"(A,ES18.9,A,I5)") 'chi^2|iter'//reg(ed_file_suffix)//'=',chi," | ",iter,"  <--  all Orbs, Spin"//reg(txtfy(ispin))
  if(ed_verbose<2)then
     suffix="_ALLorb_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
     unit=free_unit()
     open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
     write(unit,"(ES18.9,1x,I5)") chi,iter
     close(unit)
  endif
  if(ed_verbose<2)call write_bath(dmft_bath,LOGfile)
  unit=free_unit()
  open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
  call write_bath(dmft_bath,unit)
  close(unit)
  if(ed_verbose<3)call write_fit_result(ispin)
  call deallocate_bath(dmft_bath)
  deallocate(Fdelta,Xdelta,Wdelta)
  deallocate(getIorb,getJorb)
contains
  subroutine write_fit_result(ispin)
    integer                              :: i,j,l,m,iorb,jorb,ispin,jspin
    real(8)                              :: w
    character(len=20)                    :: suffix
    integer                              :: unit
    complex(8),dimension(Norb,Norb)      :: gwf
    complex(8),dimension(totNorb,Ldelta) :: fgand
    do i=1,Ldelta
       w=Xdelta(i)
       if(cg_scheme=='weiss')then
          do l=1,Norb
             gwf(l,l)=xi*w + xmu - impHloc(ispin,ispin,l,l) - delta_bath_mats(ispin,l,l,xi*w,dmft_bath)
             do m=l+1,Norb
                gwf(l,m) = - impHloc(ispin,ispin,l,m) - delta_bath_mats(ispin,l,m,xi*w,dmft_bath)
                gwf(m,l) = - impHloc(ispin,ispin,m,l) - delta_bath_mats(ispin,m,l,xi*w,dmft_bath)
             enddo
          enddo
          call matrix_inverse(gwf)
       else
          do l=1,Norb
             gwf(l,l)= delta_bath_mats(ispin,l,l,xi*w,dmft_bath)
             do m=l+1,Norb
                gwf(l,m) = delta_bath_mats(ispin,l,m,xi*w,dmft_bath)
                gwf(m,l) = delta_bath_mats(ispin,m,l,xi*w,dmft_bath)
             enddo
          enddo
       endif
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          fgand(l,i)=gwf(iorb,jorb)
       enddo
    enddo
    !
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//reg(ed_file_suffix)
       unit=free_unit()
       open(unit,file="fit_delta"//reg(suffix)//".ed")
       do i=1,Ldelta
          write(unit,"(5F24.15)")Xdelta(i),dimag(Fdelta(l,i)),dimag(fgand(l,i)),&
               dreal(Fdelta(l,i)),dreal(fgand(l,i))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_hybrid_normal




!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
! The gradient is evaluated analytically d\chi^2.
! \Delta_Anderson  and its gradient d\Delta_Anderson are evaluated 
! as separate functions.  
! HYBRID bath. 
! SPIN  DIAGONAL ; ORBITAL NON-DIAGONAL
!+-------------------------------------------------------------+
function chi2_delta_hybrd(a) result(chi2)
  real(8),dimension(:)                 :: a
  real(8),dimension(totNorb)           :: chi_orb
  complex(8),dimension(Ldelta)         :: Delta_orb
  real(8)                              :: chi2
  integer                              :: i,l,iorb,jorb
  chi_orb = 0.d0 
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     do i=1,Ldelta
        Delta_orb(i)= fg_delta_hybrd(xdelta(i),iorb,jorb,a)
     enddo
     chi_orb(l) = sum(abs(Fdelta(l,:)-Delta_orb(:))**2/Wdelta(:))
  enddo
  chi2=sum(chi_orb)
end function chi2_delta_hybrd
! the analytic GRADIENT of \chi^2
!+-------------------------------------------------------------+
function dchi2_delta_hybrd(a) result(dchi2)
  real(8),dimension(:)                 :: a
  real(8),dimension(size(a))           :: dchi2
  real(8),dimension(totNorb,size(a))   :: df
  complex(8),dimension(Ldelta)         :: g0
  complex(8),dimension(Ldelta,size(a)) :: dg0
  integer                              :: i,j,l,iorb,jorb
  df=0.d0
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     do i=1,Ldelta
        g0(i)    = fg_delta_hybrd(xdelta(i),iorb,jorb,a)
        dg0(i,:) = grad_fg_delta_hybrd(xdelta(i),iorb,jorb,a)
     enddo
     do j=1,size(a)
        df(l,j)=&
             sum( dreal(Fdelta(l,:)-g0(:))*dreal(dg0(:,j))/Wdelta(:) ) + &
             sum( dimag(Fdelta(l,:)-g0(:))*dimag(dg0(:,j))/Wdelta(:) )
     enddo
  enddo
  dchi2 = -2.d0*sum(df,1)     !sum over all orbital indices
end function dchi2_delta_hybrd
! the \Delta_Anderson function used in \chi^2 and d\chi^2
! \Delta_ab = \sum_l V^a_l*V^b_l/(iw-e_l)
!+-------------------------------------------------------------+
function fg_delta_hybrd(w,orb1,orb2,a) result(gg)
  real(8)                      :: w
  integer                      :: orb1,orb2
  real(8),dimension(:)         :: a
  real(8),dimension(Nbath)     :: eps
  real(8),dimension(Norb,Nbath):: vps
  complex(8)                   :: gg,gwf(Norb,Norb)
  integer                      :: i,l,m
  eps=a(1:Nbath)
  do l=1,Norb
     vps(l,:)=a((l*Nbath+1):((l+1)*Nbath))
  enddo
  gg=zero
  do i=1,Nbath
     gg=gg + vps(orb1,i)*vps(orb2,i)/(xi*w-eps(i))
  enddo
end function fg_delta_hybrd
! the gradient d\Delta_Anderson function used in d\chi^2
! d\Delta_ab = \grad_{V^c_k,e_k}\sum_l V^a_l*V^c_l/(iw-e_l)
!+-------------------------------------------------------------+
function grad_fg_delta_hybrd(w,orb1,orb2,a) result(dgz)
  real(8)                         :: w,sgn
  integer                         :: orb1,orb2
  real(8),dimension(:)            :: a
  real(8),dimension(Nbath)        :: eps
  real(8),dimension(Norb,Nbath)   :: vps
  complex(8),dimension(size(a))   :: dgz
  integer                         :: i,l
  dgz=zero
  eps=a(1:Nbath)
  do l=1,Norb
     vps(l,:)=a((l*Nbath+1):((l+1)*Nbath))
  enddo
  !
  do i=1,Nbath
     dgz(i)    = vps(orb1,i)*vps(orb2,i)/(xi*w-eps(i))**2
     if(orb1==orb2)then
        dgz(orb1*Nbath+i) = 2.d0*Vps(orb1,i)/(xi*w-eps(i))
     else
        dgz(orb1*Nbath+i) = Vps(orb2,i)/(xi*w-eps(i))
        dgz(orb2*Nbath+i) = Vps(orb1,i)/(xi*w-eps(i))
     endif
  enddo
end function grad_fg_delta_hybrd




!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distanec of G_0  function 
! The Gradient is not evaluated, so the minimization requires 
! a numerical estimate of the gradient. 
! G_0 is evaluated in a function below
! HYBRID bath.
! SPIN DIAGONAL & ORBITAL OFF-DIAGONAL
!+-------------------------------------------------------------+
function chi2_weiss_hybrd(a) result(chi2)
  real(8),dimension(:)                 :: a
  real(8),dimension(totNorb)           :: chi_orb
  complex(8),dimension(Ldelta)         :: Delta_orb
  real(8)                              :: chi2
  integer                              :: i,l,iorb,jorb
  chi_orb = 0.d0 
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     do i=1,Ldelta
        Delta_orb(i)= fg_weiss_hybrd(xdelta(i),iorb,jorb,a)
     enddo
     chi_orb(l) = sum(abs(Fdelta(l,:)-Delta_orb(:))**2/Wdelta(:))
  enddo
  chi2=sum(chi_orb)
end function chi2_weiss_hybrd
!
! the inverse non-interacting GF (~inverse Weiss Field) 
! used in \chi^2(\caG_0 - G_0) 
! {\GG_0}_ab = {[iw_n + xmu - H_loc - \Delta]^-1}_ab
!+-------------------------------------------------------------+
function fg_weiss_hybrd(w,orb1,orb2,a) result(gg)
  real(8)                      :: w
  integer                      :: orb1,orb2
  real(8),dimension(:)         :: a
  complex(8)                   :: gg,gwf(Norb,Norb)
  integer                      :: i,l,m,ispin
  ispin=Spin_indx
  do l=1,Norb
     gwf(l,l) = xi*w + xmu - impHloc(ispin,ispin,l,l) - fg_delta_hybrd(w,l,l,a) !sum(vps(l,:)**2/(xi*w-eps(:)))
     do m=l+1,Norb
        gwf(l,m) = -impHloc(ispin,ispin,l,m) - fg_delta_hybrd(w,l,m,a) !sum(vps(l,:)*vps(m,:)/(xi*w-eps(:)))
        gwf(m,l) = -impHloc(ispin,ispin,m,l) - fg_delta_hybrd(w,m,l,a) !sum(vps(m,:)*vps(l,:)/(xi*w-eps(:)))
     enddo
  enddo
  call matrix_inverse(gwf)
  gg = gwf(orb1,orb2)
end function fg_weiss_hybrd
