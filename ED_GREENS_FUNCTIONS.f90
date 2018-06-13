! include "ED_GF_SHARED.f90"
! include "ED_GF_NORMAL.f90"
! include "ED_GF_SUPERC.f90"
! include "ED_GF_NONSU2.f90"
! ! include "ED_GF_CHISPIN.f90"
! ! include "ED_GF_CHIDENS.f90"
! ! include "ED_GF_CHIPAIR.f90"
MODULE ED_GREENS_FUNCTIONS
  USE ED_GF_SHARED
  USE ED_GF_NORMAL
  USE ED_GF_SUPERC
  USE ED_GF_NONSU2
  implicit none

  public :: buildGf_impurity
  public :: buildChi_impurity

contains



  !+------------------------------------------------------------------+
  ! GREENS FUNCTIONS CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildgf_impurity()
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*dble(2*arange(1,Lmats)-1)
    wr     = linspace(wini,wfin,Lreal)
    !
    impGmats=zero
    impGreal=zero
    impFmats=zero
    impFreal=zero
    !
    impSmats = zero
    impSreal = zero
    impSAmats = zero
    impSAreal = zero
    !
    impG0mats=zero
    impG0real=zero
    impF0mats=zero
    impF0real=zero
    !
    if(.not.allocated(impDeltamats)) allocate(impDeltamats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(invimpG0mats)) allocate(invimpG0mats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(invimpGmats))  allocate( invimpGmats(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(impDeltareal)) allocate(impDeltareal(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(invimpG0real)) allocate(invimpG0real(Nspin,Nspin,Norb,Norb,Lreal))
    if(.not.allocated(invimpGreal))  allocate( invimpGreal(Nspin,Nspin,Norb,Norb,Lreal))
    impDeltamats=zero
    invimpGmats=zero
    invimpG0mats=zero
    impDeltareal=zero
    invimpGreal=zero
    invimpG0real=zero
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    select case(ed_mode)
    case default
       call build_gf_normal()
       call build_sigma_normal()
    case ("superc")
       call build_gf_superc()
       call build_sigma_superc()
    case ("nonsu2")
       call build_gf_nonsu2()
       call build_sigma_nonsu2()
    end select
    !
    if(MPI_MASTER)then
       if(ed_print_Sigma)call ed_print_impSigma()
       if(ed_print_G)call ed_print_impG()
       if(ed_print_G0)call ed_print_impG0()
    endif
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    if(allocated(invimpG0mats))deallocate(invimpG0mats)
    if(allocated(invimpGmats))deallocate(invimpGmats)
    if(allocated(impDeltamats))deallocate(impDeltamats)
    if(allocated(invimpG0real))deallocate(invimpG0real)
    if(allocated(invimpGreal))deallocate(invimpGreal)
    if(allocated(impDeltareal))deallocate(impDeltareal)
  end subroutine buildgf_impurity








  !+------------------------------------------------------------------+
  ! SUSCEPTIBILITY CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildChi_impurity()
    integer :: i,iorb,jorb
    !
    call allocate_grids
    !
    ! !
    ! !BUILD SPIN SUSCEPTIBILITY
    ! if(.not.allocated(spinChi_tau)) stop "buildChi_impurity: spinChi_tau not allocated"
    ! if(.not.allocated(spinChi_w))  stop "buildChi_impurity: spinChi_w not allocated"
    ! if(.not.allocated(spinChi_iv)) stop "buildChi_impurity: spinChi_iv not allocated"
    ! spinChi_tau=zero
    ! spinChi_w=zero
    ! spinChi_iv=zero
    ! write(LOGfile,"(A)")"Get impurity spin Chi:"
    ! do iorb=1,Norb
    !    write(LOGfile,"(A)")"Get Chi_spin_l"//reg(txtfy(iorb))
    !    if(MPI_MASTER)call start_timer()
    !    call lanc_ed_build_spinChi_c(iorb)
    !    if(MPI_MASTER)call stop_timer(LOGfile)
    ! enddo
    ! if(Norb>1)then
    !    write(LOGfile,"(A)")"Get Chi_spin_tot"
    !    if(MPI_MASTER)call start_timer()
    !    call lanc_ed_build_spinChi_tot_c()
    !    if(MPI_MASTER)call stop_timer(LOGfile)
    ! endif
    ! spinChi_tau = SpinChi_tau/zeta_function
    ! spinChi_w   = spinChi_w/zeta_function
    ! spinChi_iv  = spinChi_iv/zeta_function
    ! !
    ! ! !BUILD CHARGE SUSCEPTIBILITY
    ! if(.not.allocated(densChi_tau)) stop "buildChi_impurity: densChi_tau not allocated"
    ! if(.not.allocated(densChi_w))  stop "buildChi_impurity: densChi_w not allocated"
    ! if(.not.allocated(densChi_iv)) stop "buildChi_impurity: densChi_iv not allocated"
    ! if(.not.allocated(densChi_mix_tau))stop "buildChi_impurity: densChi_mix_tau not allocated"
    ! if(.not.allocated(densChi_mix_w))  stop "buildChi_impurity: densChi_mix_w not allocated"
    ! if(.not.allocated(densChi_mix_iv)) stop "buildChi_impurity: densChi_mix_iv not allocated"
    ! if(.not.allocated(densChi_tot_tau))stop "buildChi_impurity: densChi_tot_tau not allocated"
    ! if(.not.allocated(densChi_tot_w))  stop "buildChi_impurity: densChi_tot_w not allocated"
    ! if(.not.allocated(densChi_tot_iv)) stop "buildChi_impurity: densChi_tot_iv not allocated"
    ! densChi_tau=zero
    ! densChi_w=zero
    ! densChi_iv=zero
    ! densChi_mix_tau=zero
    ! densChi_mix_w=zero
    ! densChi_mix_iv=zero
    ! densChi_tot_tau=zero
    ! densChi_tot_w=zero
    ! densChi_tot_iv=zero
    ! !
    ! write(LOGfile,"(A)")"Get impurity dens Chi:"
    ! do iorb=1,Norb
    !    write(LOGfile,"(A)")"Get Chi_dens_diag_l"//reg(txtfy(iorb))
    !    if(MPI_MASTER)call start_timer()
    !    call lanc_ed_build_densChi_diag_c(iorb)
    !    if(MPI_MASTER)call stop_timer(LOGfile)
    ! enddo
    ! !
    ! !
    ! if(Norb>1)then
    !    do iorb=1,Norb
    !       do jorb=iorb+1,Norb
    !          write(LOGfile,"(A)")"Get Chi_dens_offdiag_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
    !          if(MPI_MASTER)call start_timer()
    !          call lanc_ed_build_densChi_offdiag_c(iorb,jorb)
    !          if(MPI_MASTER)call stop_timer(LOGfile)
    !       end do
    !    end do
    !    do iorb=1,Norb
    !       do jorb=iorb+1,Norb
    !          denschi_w(iorb,jorb,:) = 0.5d0*( denschi_w(iorb,jorb,:) - (one+xi)*denschi_w(iorb,iorb,:) - (one+xi)*denschi_w(jorb,jorb,:))
    !       enddo
    !    enddo
    !    !
    !    do iorb=1,Norb
    !       do jorb=1,Norb
    !          write(LOGfile,"(A)")"Get Chi_dens_offdiag_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
    !          if(MPI_MASTER)call start_timer()
    !          call lanc_ed_build_densChi_mix_c(iorb,jorb)
    !          if(MPI_MASTER)call stop_timer(LOGfile)
    !       end do
    !    end do
    !    !
    !    write(LOGfile,"(A)")"Get Chi_dens_tot"
    !    if(MPI_MASTER)call start_timer()
    !    call lanc_ed_build_densChi_tot_c()
    !    if(MPI_MASTER)call stop_timer(LOGfile)
    ! endif
    ! !
    ! denschi_tau = Denschi_tau/zeta_function
    ! denschi_w   = denschi_w/zeta_function
    ! denschi_iv  = denschi_iv/zeta_function
    ! !
    ! !
    ! !BUILD PAIR SUSCEPTIBILITY
    ! if(.not.allocated(pairChi_tau))stop "buildChi_impurity: pairChi_tau not allocated"
    ! if(.not.allocated(pairChi_w))  stop "buildChi_impurity: pairChi_w not allocated"
    ! if(.not.allocated(pairChi_iv)) stop "buildChi_impurity: pairChi_iv not allocated"
    ! pairChi_tau=zero
    ! pairChi_w=zero
    ! pairChi_iv=zero
    ! write(LOGfile,"(A)")"Get impurity pair Chi:"
    ! do iorb=1,Norb
    !    write(LOGfile,"(A)")"Get Chi_pair_l"//reg(txtfy(iorb))
    !    if(MPI_MASTER)call start_timer()
    !    call lanc_ed_build_pairChi_c(iorb)
    !    if(MPI_MASTER)call stop_timer(LOGfile)
    ! enddo
    ! pairChi_tau = PairChi_tau/zeta_function
    ! pairChi_w   = pairChi_w/zeta_function
    ! pairChi_iv  = pairChi_iv/zeta_function
    ! !
    ! !PRINTING:
    ! if(MPI_MASTER)call ed_print_impChi()
    ! !
    !
    call deallocate_grids
  end subroutine buildChi_impurity





END MODULE ED_GREENS_FUNCTIONS
