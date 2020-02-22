!+-------------------------------------------------------------------+
!PURPOSE  : Allocate the ED bath
!+-------------------------------------------------------------------+
subroutine allocate_dmft_bath(dmft_bath_)
  type(effective_bath) :: dmft_bath_
  if(dmft_bath_%status)call deallocate_dmft_bath(dmft_bath_)
  !
  select case(bath_type)
  case default
     !
     select case(ed_mode)
     case default                                         !normal [N,Sz]
        allocate(dmft_bath_%e(Nspin,Norb,Nbath))          !local energies of the bath
        allocate(dmft_bath_%v(Nspin,Norb,Nbath))          !same-spin hybridization
     case ("superc")                                      !superc [Sz]
        allocate(dmft_bath_%e(Nspin,Norb,Nbath))          !local energies of the bath
        allocate(dmft_bath_%d(Nspin,Norb,Nbath))          !local SC order parameters the bath
        allocate(dmft_bath_%v(Nspin,Norb,Nbath))          !same-spin hybridization
     case ("nonsu2")                                      !nonsu2 [N]
        allocate(dmft_bath_%e(Nspin,Norb,Nbath))          !local energies of the bath
        allocate(dmft_bath_%v(Nspin,Norb,Nbath))          !same-spin hybridization
        allocate(dmft_bath_%u(Nspin,Norb,Nbath))          !spin-flip hybridization
     end select
     !
  case('hybrid')
     !
     select case(ed_mode)
     case default                                         !normal  [N,Sz]
        allocate(dmft_bath_%e(Nspin,1,Nbath))             !local energies of the bath
        allocate(dmft_bath_%v(Nspin,Norb,Nbath))          !same-spin hybridization
     case ("superc")                                      !superc  [Sz]
        allocate(dmft_bath_%e(Nspin,1,Nbath))             !local energies of the bath
        allocate(dmft_bath_%d(Nspin,1,Nbath))             !local SC order parameters the bath
        allocate(dmft_bath_%v(Nspin,Norb,Nbath))          !same-spin hybridization
     case ("nonsu2")                                      !nonsu2 case [N] qn
        allocate(dmft_bath_%e(Nspin,1,Nbath))             !local energies of the bath
        allocate(dmft_bath_%v(Nspin,Norb,Nbath))          !same-spin hybridization
        allocate(dmft_bath_%u(Nspin,Norb,Nbath))          !spin-flip hybridization
     end select
     !
  case('replica')
     !
     allocate(dmft_bath_%h(Nspin,Nspin,Norb,Norb,Nbath))  !replica H_k of the bath
     allocate(dmft_bath_%t(Nbath))                        !replica hybridization
     allocate(dmft_bath_%o(Nop,Nbath))                    !replica coupling with fixed operator
     allocate(dmft_bath_%mask(Nspin,Nspin,Norb,Norb,2))   !mask on components
     !
  end select
  dmft_bath_%status=.true.
end subroutine allocate_dmft_bath




!+-------------------------------------------------------------------+
!PURPOSE  : Deallocate the ED bath
!+-------------------------------------------------------------------+
subroutine deallocate_dmft_bath(dmft_bath_)
  type(effective_bath) :: dmft_bath_
  if(allocated(dmft_bath_%e))   deallocate(dmft_bath_%e)
  if(allocated(dmft_bath_%d))   deallocate(dmft_bath_%d)
  if(allocated(dmft_bath_%v))   deallocate(dmft_bath_%v)
  if(allocated(dmft_bath_%u))   deallocate(dmft_bath_%u)
  if(allocated(dmft_bath_%t))   deallocate(dmft_bath_%t)
  if(allocated(dmft_bath_%o))   deallocate(dmft_bath_%o)
  if(allocated(dmft_bath_%h))   deallocate(dmft_bath_%h)
  if(allocated(dmft_bath_%mask))deallocate(dmft_bath_%mask)
  dmft_bath_%status=.false.
end subroutine deallocate_dmft_bath




!+------------------------------------------------------------------+
!PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or
!reading previous (converged) solution
!+------------------------------------------------------------------+
subroutine init_dmft_bath(dmft_bath_)
  type(effective_bath) :: dmft_bath_
  complex(8)           :: hrep_aux(Nspin*Norb,Nspin*Norb)
  real(8)              :: hybr_aux_R,hybr_aux_I
  real(8)              :: hrep_aux_R(Nspin*Norb,Nspin*Norb)
  real(8)              :: hrep_aux_I(Nspin*Norb,Nspin*Norb)
  real(8)              :: re,im
  integer              :: i,unit,flen,Nh
  integer              :: io,jo,iorb,ispin,jorb,jspin,iop
  logical              :: IOfile
  real(8)              :: de,noise_tot
  real(8),allocatable  :: noise_b(:),noise_s(:),noise_o(:),noise_m(:)
  character(len=21)    :: space
  if(.not.dmft_bath_%status)stop "init_dmft_bath error: bath not allocated"
  allocate(noise_b(Nbath));noise_b=0.d0
  allocate(noise_s(Nspin));noise_s=0.d0
  allocate(noise_o(Norb)); noise_o=0.d0
  allocate(noise_m(Nop));  noise_m=0.d0
  call random_number(noise_b)
  call random_number(noise_s)
  call random_number(noise_o)
  call random_number(noise_m)
  noise_b=noise_b*ed_bath_noise_thr
  noise_s=noise_s*ed_bath_noise_thr
  noise_o=noise_o*ed_bath_noise_thr
  noise_m=noise_m*ed_bath_noise_thr

  select case(bath_type)
  case default
     !Get energies:
     dmft_bath_%e(:,:,1)    =-hwband + noise_b(1)
     dmft_bath_%e(:,:,Nbath)= hwband + noise_b(Nbath)
     Nh=Nbath/2
     if(mod(Nbath,2)==0.and.Nbath>=4)then
        de=hwband/max(Nh-1,1)
        dmft_bath_%e(:,:,Nh)  = -1.d-3 + noise_b(Nh)
        dmft_bath_%e(:,:,Nh+1)=  1.d-3 + noise_b(Nh+1)
        do i=2,Nh-1
           dmft_bath_%e(:,:,i)   =-hwband + (i-1)*de  + noise_b(i)
           dmft_bath_%e(:,:,Nbath-i+1)= hwband - (i-1)*de + noise_b(Nbath-i+1)
        enddo
     elseif(mod(Nbath,2)/=0.and.Nbath>=3)then
        de=hwband/Nh
        dmft_bath_%e(:,:,Nh+1)= 0.0d0 + noise_b(Nh+1)
        do i=2,Nh
           dmft_bath_%e(:,:,i)        =-hwband + (i-1)*de + noise_b(i)
           dmft_bath_%e(:,:,Nbath-i+1)= hwband - (i-1)*de + noise_b(Nbath-i+1)
        enddo
     endif
     !Get spin-keep yhbridizations
     do i=1,Nbath
        dmft_bath_%v(:,:,i)=max(0.1d0,1.d0/sqrt(dble(Nbath)))+noise_b(i)
     enddo
     !Get SC amplitudes
     if(ed_mode=="superc")dmft_bath_%d(:,:,:)=deltasc
     !Get spin-flip hybridizations
     if(ed_mode=="nonsu2")then
        do i=1,Nbath
           dmft_bath_%u(:,:,i) = dmft_bath_%v(:,:,i)*ed_vsf_ratio+noise_b(i)
        enddo
     endif
     !
  case('replica')
     !BATH INITIALIZATION
     dmft_bath_%h=zero
     dmft_bath_%o=zero
     if (replica_operators) then
        do i=1,Nbath
           do iop=1,Nop
              dmft_bath_%o(iop,i)=1.d0*((-1)**iop)+noise_m(iop)+noise_b(i)
              dmft_bath_%h(:,:,:,:,i) = dmft_bath_%h(:,:,:,:,i) + dmft_bath_%o(iop,i) * OpMat(:,:,:,:,iop)
           enddo
        enddo
     else
        do i=1,Nbath
           !
           dmft_bath_%h(:,:,:,:,i)=impHloc-(xmu+noise_b(i))*so2nn_reshape(eye(Nspin*Norb),Nspin,Norb)
           !
        enddo
     endif
     !HYBR. INITIALIZATION
     dmft_bath_%t=zero
     do i=1,Nbath
        noise_tot=noise_b(i)
        dmft_bath_%t(i)=0.5d0+noise_b(i)
     enddo
     !
     deallocate(noise_b,noise_s,noise_o)
     !
  end select
  !
  !Read from file if exist:
  !
  inquire(file=trim(Hfile)//trim(ed_file_suffix)//".restart",exist=IOfile)
  if(IOfile)then
     write(LOGfile,"(A)")'Reading bath from file'//trim(Hfile)//trim(ed_file_suffix)//".restart"
     unit = free_unit()
     flen = file_length(trim(Hfile)//trim(ed_file_suffix)//".restart")
     !
     open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
     !
     select case(bath_type)
     case default
        !
        read(unit,*)
        select case(ed_mode)
        case default
           do i=1,min(flen,Nbath)
              read(unit,*)((&
                   dmft_bath_%e(ispin,iorb,i),&
                   dmft_bath_%v(ispin,iorb,i),&
                   iorb=1,Norb),ispin=1,Nspin)
           enddo
        case ("superc")
           do i=1,min(flen,Nbath)
              read(unit,*)((&
                   dmft_bath_%e(ispin,iorb,i),&
                   dmft_bath_%d(ispin,iorb,i),&
                   dmft_bath_%v(ispin,iorb,i),&
                   iorb=1,Norb),ispin=1,Nspin)
           enddo
        case("nonsu2")
           do i=1,min(flen,Nbath)
              read(unit,*)((&
                   dmft_bath_%e(ispin,iorb,i),&
                   dmft_bath_%v(ispin,iorb,i),&
                   dmft_bath_%u(ispin,iorb,i),&
                   iorb=1,Norb),ispin=1,Nspin)
           enddo
        end select
        !
     case ('hybrid')
        read(unit,*)
        !
        select case(ed_mode)
        case default
           do i=1,min(flen,Nbath)
              read(unit,*)(&
                   dmft_bath_%e(ispin,1,i),&
                   (&
                   dmft_bath_%v(ispin,iorb,i),&
                   iorb=1,Norb),&
                   ispin=1,Nspin)
           enddo
        case ("superc")
           do i=1,min(flen,Nbath)
              read(unit,*)(&
                   dmft_bath_%e(ispin,1,i),&
                   dmft_bath_%d(ispin,1,i),&
                   (&
                   dmft_bath_%v(ispin,iorb,i),&
                   iorb=1,Norb),&
                   ispin=1,Nspin)
           enddo
        case ("nonsu2")
           do i=1,min(flen,Nbath)
              read(unit,*)(&
                   dmft_bath_%e(ispin,1,i),&
                   (&
                   dmft_bath_%v(ispin,iorb,i),&
                   dmft_bath_%u(ispin,iorb,i),&
                   iorb=1,Norb),&
                   ispin=1,Nspin)
           enddo
        end select
        !
     case ('replica')
        !
        select case(ed_mode)
        case ("normal","nonsu2")
           do i=1,Nbath
              hrep_aux_R=0.0d0;hrep_aux_I=0.0d0
              hybr_aux_R=0.0d0;hybr_aux_I=0.0d0
              hrep_aux=zero
              do io=1,Nspin*Norb
                 if(io==1)read(unit,"(90(F21.12,1X))")     hybr_aux_R,hybr_aux_I,(hrep_aux_R(io,jo),jo=1,Nspin*Norb),(hrep_aux_I(io,jo),jo=1,Nspin*Norb)
                 if(io/=1)read(unit,"(2a21,90(F21.12,1X))")   space  ,   space  ,(hrep_aux_R(io,jo),jo=1,Nspin*Norb),(hrep_aux_I(io,jo),jo=1,Nspin*Norb)
              enddo
              read(unit,*)
              hrep_aux=cmplx(hrep_aux_R,hrep_aux_I)
              dmft_bath_%h(:,:,:,:,i)=so2nn_reshape(hrep_aux,Nspin,Norb)
              dmft_bath_%t(i)=hybr_aux_R
              if (replica_operators) then
                 read(unit,"(90(F21.12,1X))") (dmft_bath_%o(iop,i),iop=1,Nop)
                 read(unit,*)
                 do iop=2,Nop
                    dmft_bath_%h(:,:,:,:,i) = dmft_bath_%h(:,:,:,:,i) + dmft_bath_%o(iop,i) * OpMat(:,:,:,:,iop)
                 enddo
              endif
           enddo
           !
        case ("superc")
           !
        end select
        !
     end select
     close(unit)
  endif
end subroutine init_dmft_bath

!+-------------------------------------------------------------------+
!PURPOSE  : set the mask based on impHloc in the replica bath topology
!+-------------------------------------------------------------------+
subroutine init_dmft_bath_mask(dmft_bath_)
  type(effective_bath),intent(inout) :: dmft_bath_
  integer                            :: iorb,ispin,jorb,jspin
  integer                            :: ibath,iop
  integer                            :: io,jo
  !
  dmft_bath_%mask=.false.
  !
  if (replica_operators) then
     !
     ! Mask initialization on user given operators (check all bath and all operators)
     if(.not.(allocated(OpMat))) stop "operatos not allocated on mask initialization"
     if(bath_type/="replica") stop "operators bath representation requires bath_type=replica"
     !
     do iop=1,Nop
        do ispin=1,Nspin
           do iorb=1,Norb
              !Re-diagonal elements always present
              dmft_bath_%mask(ispin,ispin,iorb,iorb,1)=.true.
              !Im-diagonal elements checked
              if(abs(aimag(OpMat(ispin,ispin,iorb,iorb,iop))).gt.1e-6)stop "one element of operatorial representation causes non Hermiticity of the bath"
              do jspin=1,Nspin
                 do jorb=1,Norb
                    io = iorb + (ispin-1)*Norb
                    jo = jorb + (jspin-1)*Norb
                    if(io/=jo)then
                       !Re
                       if( abs(real(OpMat(ispin,jspin,iorb,jorb,iop))).gt.1e-6)then
                          dmft_bath_%mask(ispin,jspin,iorb,jorb,1)=.true.
                       endif
                       !Im
                       if(abs(aimag(OpMat(ispin,jspin,iorb,jorb,iop))).gt.1e-6)then
                          dmft_bath_%mask(ispin,jspin,iorb,jorb,2)=.true.
                       endif
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
     !
  else
     !
     ! Mask initialization on Hloc
     if(.not.(allocated(impHloc))) stop "impHloc not allocated on mask initialization"
     !
     do ispin=1,Nspin
        do iorb=1,Norb
           !Re-diagonal elements always present
           dmft_bath_%mask(ispin,ispin,iorb,iorb,1)=.true.
           !Im-diagonal elements checked
           if(abs(aimag(impHloc(ispin,ispin,iorb,iorb))).gt.1e-6)stop "impHloc is not Hermitian"
           !off-diagonal elements
           do jspin=1,Nspin
              do jorb=1,Norb
                 io = iorb + (ispin-1)*Norb
                 jo = jorb + (jspin-1)*Norb
                 if(io/=jo)then
                    !Re
                    if( abs(real(impHloc(ispin,jspin,iorb,jorb))).gt.1e-6)then
                       dmft_bath_%mask(ispin,jspin,iorb,jorb,1)=.true.
                    endif
                    !Im
                    if(abs(aimag(impHloc(ispin,jspin,iorb,jorb))).gt.1e-6)then
                       dmft_bath_%mask(ispin,jspin,iorb,jorb,2)=.true.
                    endif
                 endif
              enddo
           enddo
        enddo
     enddo
     !
  endif

  !
end subroutine init_dmft_bath_mask


!+-------------------------------------------------------------------+
!PURPOSE  : write out the bath to a given unit with
! the following column formatting:
! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
!+-------------------------------------------------------------------+
subroutine write_dmft_bath(dmft_bath_,unit)
  type(effective_bath) :: dmft_bath_
  integer,optional     :: unit
  integer              :: unit_
  integer              :: i
  integer              :: io,jo,iorb,ispin,iop
  complex(8)           :: hybr_aux
  complex(8)           :: hrep_aux(Nspin*Norb,Nspin*Norb)
  unit_=LOGfile;if(present(unit))unit_=unit
  if(.not.dmft_bath_%status)stop "write_dmft_bath error: bath not allocated"
  select case(bath_type)
  case default
     !
     select case(ed_mode)
     case default
        write(unit_,"(90(A21,1X))")&
             ((&
             "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             iorb=1,Norb),ispin=1,Nspin)
        do i=1,Nbath
           write(unit_,"(90(F21.12,1X))")((&
                dmft_bath_%e(ispin,iorb,i),&
                dmft_bath_%v(ispin,iorb,i),&
                iorb=1,Norb),ispin=1,Nspin)
        enddo
     case ("superc")
        write(unit_,"(90(A21,1X))")&
             ((&
             "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             "Dk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)) ,&
             "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             iorb=1,Norb),ispin=1,Nspin)
        do i=1,Nbath
           write(unit_,"(90(F21.12,1X))")((&
                dmft_bath_%e(ispin,iorb,i),&
                dmft_bath_%d(ispin,iorb,i),&
                dmft_bath_%v(ispin,iorb,i),&
                iorb=1,Norb),ispin=1,Nspin)
        enddo
     case ("nonsu2")
        write(unit_,"(90(A21,1X))")&
             ((&
             "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             "Vak_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             "Vbk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             iorb=1,Norb), ispin=1,Nspin)
        do i=1,Nbath
           write(unit,"(90(F21.12,1X))")((&
                dmft_bath_%e(ispin,iorb,i),&
                dmft_bath_%v(ispin,iorb,i),&
                dmft_bath_%u(ispin,iorb,i),&
                iorb=1,Norb),ispin=1,Nspin)
        enddo
     end select
     !
  case('hybrid')
     !
     select case(ed_mode)
     case default
        write(unit_,"(90(A21,1X))")(&
             "#Ek_s"//reg(txtfy(ispin)),&
             ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),&
             ispin=1,Nspin)
        do i=1,Nbath
           write(unit_,"(90(F21.12,1X))")(&
                dmft_bath_%e(ispin,1,i),&
                (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),&
                ispin=1,Nspin)
        enddo
     case ("superc")
        write(unit_,"(90(A21,1X))")(&
             "#Ek_s"//reg(txtfy(ispin)),&
             "Dk_s"//reg(txtfy(ispin)) ,&
             ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),&
             ispin=1,Nspin)
        do i=1,Nbath
           write(unit_,"(90(F21.12,1X))")(&
                dmft_bath_%e(ispin,1,i),&
                dmft_bath_%d(ispin,1,i),&
                (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),&
                ispin=1,Nspin)
        enddo
     case ("nonsu2")
        write(unit_,"(90(A21,1X))")(&
             "#Ek_s"//reg(txtfy(ispin)),&
             (&
             "Vak_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             "Vbk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
             iorb=1,Norb),&
             ispin=1,Nspin)
        do i=1,Nbath
           write(unit_,"(90(F21.12,1X))")(&
                dmft_bath_%e(ispin,1,i),    &
                (dmft_bath_%v(ispin,iorb,i),dmft_bath_%u(ispin,iorb,i),iorb=1,Norb),&
                ispin=1,Nspin)
        enddo
     end select
     !
  case ('replica')
     !
     select case(ed_mode)
     case ("normal","nonsu2")
        do i=1,Nbath
           hrep_aux=zero;hrep_aux=nn2so_reshape(dmft_bath_%h(:,:,:,:,i),Nspin,Norb)
           hybr_aux=dmft_bath_%t(i)
           do io=1,Nspin*Norb
              if(unit_==LOGfile)then
                 if(io==1) write(unit_,"(2F8.3,a5,90(F8.3,1X))") real(hybr_aux),aimag(hybr_aux),"|",( real(hrep_aux(io,jo)),jo=1,Nspin*Norb),&
                     (aimag(hrep_aux(io,jo)),jo=1,Nspin*Norb)
                 if(io/=1) write(unit_,"(2a8,a5,90(F8.3,1X))")        "  "     ,      "  "     ,"|",( real(hrep_aux(io,jo)),jo=1,Nspin*Norb),&
                      (aimag(hrep_aux(io,jo)),jo=1,Nspin*Norb)
              else
                 if(io==1)write(unit_,"(90(F21.12,1X))")         real(hybr_aux),aimag(hybr_aux),(real(hrep_aux(io,jo)),jo=1,Nspin*Norb),(aimag(hrep_aux(io,jo)),jo=1,Nspin*Norb)
                 if(io/=1)write(unit_,"(2a21,90(F21.12,1X))")         "  "     ,     "  "      ,(real(hrep_aux(io,jo)),jo=1,Nspin*Norb),(aimag(hrep_aux(io,jo)),jo=1,Nspin*Norb)
              endif
           enddo
           write(unit_,*)
           if (replica_operators) then
              if(unit_==LOGfile) then
                 write(unit_,"(90(F8.3,1X))") (dmft_bath_%o(iop,i),iop=1,Nop)
              else
                 write(unit_,"(90(F21.12,1X))") (dmft_bath_%o(iop,i),iop=1,Nop)
              endif
              write(unit_,*)
           endif
        enddo
        !
     case ("superc")
        !
     end select
     !
  end select
end subroutine write_dmft_bath






!+-------------------------------------------------------------------+
!PURPOSE  : save the bath to a given file using the write bath
! procedure and formatting:
!+-------------------------------------------------------------------+
subroutine save_dmft_bath(dmft_bath_,file,used)
  type(effective_bath)      :: dmft_bath_
  character(len=*),optional :: file
  character(len=256)        :: file_
  logical,optional          :: used
  logical                   :: used_
  character(len=16)         :: extension
  integer                   :: unit_
  ! if(ED_MPI_ID==0)then
  if(.not.dmft_bath_%status)stop "save_dmft_bath error: bath is not allocated"
  used_=.false.;if(present(used))used_=used
  extension=".restart";if(used_)extension=".used"
  file_=reg(reg(Hfile)//reg(ed_file_suffix)//reg(extension))
  if(present(file))file_=reg(file)
  unit_=free_unit()
  open(unit_,file=reg(file_))
  call write_dmft_bath(dmft_bath_,unit_)
  close(unit_)
  ! endif
end subroutine save_dmft_bath




!+-------------------------------------------------------------------+
!PURPOSE  : set the bath components from a given user provided
! bath-array
!+-------------------------------------------------------------------+
subroutine set_dmft_bath(bath_,dmft_bath_)
  real(8),dimension(:)   :: bath_
  type(effective_bath)   :: dmft_bath_
  integer                :: stride,io,jo,i
  integer                :: iorb,ispin,jorb,jspin,ibath,maxspin,iop
  logical                :: check
  real(8)                :: element_R,element_I
  !
  if(.not.dmft_bath_%status)stop "set_dmft_bath error: bath not allocated"
  check = check_bath_dimension(bath_)
  if(.not.check)stop "set_dmft_bath error: wrong bath dimensions"
  !
  select case(bath_type)
  case default
     !
     select case(ed_mode)
        !
     case default
        stride = 0
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%e(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        stride = Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%v(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        !
     case ("superc")
        stride = 0
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%e(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        stride = Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%d(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        stride = 2*Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%v(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        !
     case("nonsu2")
        stride = 0
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%e(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        stride = Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%v(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        stride = 2*Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 dmft_bath_%u(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
     end select
     !
  case ('hybrid')
     !
     select case(ed_mode)
     case default
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              dmft_bath_%e(ispin,1,i) = bath_(io)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 dmft_bath_%v(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        !
     case ("superc")
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              dmft_bath_%e(ispin,1,i) = bath_(io)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              dmft_bath_%d(ispin,1,i) = bath_(io)
           enddo
        enddo
        stride = 2*Nspin*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 dmft_bath_%v(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        !
     case("nonsu2")
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              dmft_bath_%e(ispin,1,i) = bath_(io)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 dmft_bath_%v(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
        stride = Nspin*Nbath + Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 dmft_bath_%u(ispin,iorb,i) = bath_(io)
              enddo
           enddo
        enddo
     end select
     !
  case ('replica')
     !
     select case(ed_mode)
     case("normal")
        !
        dmft_bath_%h=zero
        dmft_bath_%t=0d0
        dmft_bath_%o=0d0
        i = 0
        !
        if(replica_operators)then
           !
           do ibath=1,Nbath
              do iop=1,Nop
                 ![a, b, c]_k1  [a, b, c]_k2 [a, b, c]_k3 ...
                 i=i+1
                 dmft_bath_%o(iop,ibath)=real(bath_(i))
                 dmft_bath_%h(:,:,:,:,ibath) = dmft_bath_%h(:,:,:,:,ibath) + dmft_bath_%o(iop,ibath) * OpMat(:,:,:,:,iop)
              enddo
           enddo
           !
        else
           !
           if(ed_para)then
              Maxspin=1
           else
              Maxspin=2
           endif
           !
           do ibath=1,Nbath
              !all diagonal non-vanishing terms in imploc - real only - always present
              do ispin=1,Maxspin
                 do iorb=1,Norb
                    i=i+1
                    dmft_bath_%h(ispin,ispin,iorb,iorb,ibath)=cmplx(bath_(i),0d0)
                 enddo
              enddo
              !all off-diagonal non-vanishing terms in imploc - all spin sectors
              do ispin=1,Maxspin
                 do iorb=1,Norb
                    do jorb=1,Norb
                       io = iorb + (ispin-1)*Norb
                       jo = jorb + (ispin-1)*Norb
                       if(io.gt.jo)cycle !only upper triangular saved. Hermiticity imposed here
                       element_R=0.0d0;element_I=0.0d0
                       if(dmft_bath_%mask(ispin,ispin,iorb,jorb,1)) then
                          i=i+1
                          element_R=bath_(i)
                       endif
                       if(dmft_bath_%mask(ispin,ispin,iorb,jorb,2)) then
                          i=i+1
                          element_I=bath_(i)
                       endif
                       dmft_bath_%h(ispin,ispin,iorb,jorb,ibath)=cmplx(element_R,element_I)
                       !hermiticity
                       if(iorb/=jorb)dmft_bath_%h(ispin,ispin,jorb,iorb,ibath)=conjg(dmft_bath_%h(ispin,ispin,iorb,jorb,ibath))
                       !spin-conservation
                       if(Maxspin==1)dmft_bath_%h(2,2,iorb,jorb,ibath)=dmft_bath_%h(1,1,iorb,jorb,ibath)
                       if(Maxspin==1)dmft_bath_%h(2,2,jorb,iorb,ibath)=dmft_bath_%h(1,1,jorb,iorb,ibath)
                    enddo
                 enddo
              enddo
           enddo
           !
        endif
        !
        !all Re[Hybr]
        do ibath=1,Nbath
           i=i+1
           dmft_bath_%t(ibath)=real(bath_(i))
        enddo
        !
     case("nonsu2")
        !
        dmft_bath_%h=zero
        dmft_bath_%t=zero
        i = 0
        !
        if(replica_operators)then
           !
           do ibath=1,Nbath
              do iop=1,Nop
                 ![a, b, c]_k1  [a, b, c]_k2 [a, b, c]_k3 ...
                 i=i+1
                 dmft_bath_%o(iop,ibath)=real(bath_(i))
                 dmft_bath_%h(:,:,:,:,ibath) = dmft_bath_%h(:,:,:,:,ibath) + dmft_bath_%o(iop,ibath) * OpMat(:,:,:,:,iop)
              enddo
           enddo
           !
        else
           !
           do ibath=1,Nbath
              !all diagonal non-vanishing terms in imploc - real only - always present
              do ispin=1,Nspin
                 do iorb=1,Norb
                    i=i+1
                    dmft_bath_%h(ispin,ispin,iorb,iorb,ibath)=cmplx(bath_(i),0d0)
                 enddo
              enddo
              !all off-diagonal non-vanishing terms in imploc - all spin sectors
              do ispin=1,Nspin
                 do iorb=1,Norb
                    do jspin=1,Nspin
                       do jorb=1,Norb
                          io = iorb + (ispin-1)*Norb
                          jo = jorb + (jspin-1)*Norb
                          if(io.gt.jo)cycle !only upper triangular saved. Hermiticity imposed here
                          element_R=0.0d0;element_I=0.0d0
                          if(dmft_bath_%mask(ispin,jspin,iorb,jorb,1)) then
                             i=i+1
                             element_R=bath_(i)
                          endif
                          if(dmft_bath_%mask(ispin,jspin,iorb,jorb,2)) then
                             i=i+1
                             element_I=bath_(i)
                          endif
                          dmft_bath_%h(ispin,jspin,iorb,jorb,ibath)=cmplx(element_R,element_I)
                          !hermiticity
                          if((ispin==jspin).and.(iorb/=jorb))dmft_bath_%h(ispin,ispin,jorb,iorb,ibath)=conjg(dmft_bath_%h(ispin,ispin,iorb,jorb,ibath))
                          if((ispin/=jspin).and.(iorb==jorb))dmft_bath_%h(jspin,ispin,iorb,iorb,ibath)=conjg(dmft_bath_%h(ispin,jspin,iorb,iorb,ibath))
                          if((ispin/=jspin).and.(iorb/=jorb))dmft_bath_%h(jspin,ispin,jorb,iorb,ibath)=conjg(dmft_bath_%h(ispin,jspin,iorb,jorb,ibath))
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           !
        endif
        !
        !all Re[Hybr]
        do ibath=1,Nbath
           i=i+1
           dmft_bath_%t(ibath)=real(bath_(i))
        enddo
        !
     case ("superc")

        !
     end select
     !
  end select
end subroutine set_dmft_bath



!+-------------------------------------------------------------------+
!PURPOSE  : copy the bath components back to a 1-dim array
!+-------------------------------------------------------------------+
subroutine get_dmft_bath(dmft_bath_,bath_)
  type(effective_bath)   :: dmft_bath_
  real(8),dimension(:)   :: bath_
  integer                :: stride,io,jo,i
  integer                :: iorb,ispin,jorb,jspin,ibath,maxspin,iop
  logical                :: check
  if(.not.dmft_bath_%status)stop "get_dmft_bath error: bath not allocated"
  check=check_bath_dimension(bath_)
  if(.not.check)stop "get_dmft_bath error: wrong bath dimensions"
  !
  select case(bath_type)
  case default
     !
     select case(ed_mode)
     case default
        stride = 0
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) = dmft_bath_%e(ispin,iorb,i)
              enddo
           enddo
        enddo
        stride = Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) = dmft_bath_%v(ispin,iorb,i)
              enddo
           enddo
        enddo
        !
     case ("superc")
        stride = 0
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) = dmft_bath_%e(ispin,iorb,i)
              enddo
           enddo
        enddo
        stride = Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) =  dmft_bath_%d(ispin,iorb,i)
              enddo
           enddo
        enddo
        stride = 2*Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) =  dmft_bath_%v(ispin,iorb,i)
              enddo
           enddo
        enddo
        !
     case("nonsu2")
        stride = 0
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) = dmft_bath_%e(ispin,iorb,i)
              enddo
           enddo
        enddo
        stride = Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) = dmft_bath_%v(ispin,iorb,i)
              enddo
           enddo
        enddo
        stride = 2*Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
                 bath_(io) = dmft_bath_%u(ispin,iorb,i)
              enddo
           enddo
        enddo
     end select
     !
  case ('hybrid')
     !
     select case(ed_mode)
     case default
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              bath_(io) =  dmft_bath_%e(ispin,1,i)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 bath_(io) =  dmft_bath_%v(ispin,iorb,i)
              enddo
           enddo
        enddo
        !
     case ("superc")
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              bath_(io) =  dmft_bath_%e(ispin,1,i)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              bath_(io) =  dmft_bath_%d(ispin,1,i)
           enddo
        enddo
        stride = 2*Nspin*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 bath_(io) =  dmft_bath_%v(ispin,iorb,i)
              enddo
           enddo
        enddo
        !
     case("nonsu2")
        stride = 0
        do ispin=1,Nspin
           do i=1,Nbath
              io = stride + i + (ispin-1)*Nbath
              bath_(io) = dmft_bath_%e(ispin,1,i)
           enddo
        enddo
        stride = Nspin*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 bath_(io) = dmft_bath_%v(ispin,iorb,i)
              enddo
           enddo
        enddo
        stride = Nspin*Nbath + Nspin*Norb*Nbath
        do ispin=1,Nspin
           do iorb=1,Norb
              do i=1,Nbath
                 io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
                 bath_(io) = dmft_bath_%u(ispin,iorb,i)
              enddo
           enddo
        enddo
     end select
     !
  case ('replica')
     !
     select case(ed_mode)
     case("normal")
        !
        i = 0
        !
        if(replica_operators)then
           !
           do ibath=1,Nbath
              do iop=1,Nop
                 ![a, b, c]_k1  [a, b, c]_k2 [a, b, c]_k3 ...
                 i=i+1
                 bath_(i)=real(dmft_bath_%o(iop,ibath))
              enddo
           enddo
           !
        else
           !
           if(ed_para)then
              Maxspin=1
           else
              Maxspin=2
           endif
           !
           do ibath=1,Nbath
              !all diagonal non-vanishing terms in imploc - always present
              do ispin=1,Maxspin
                 do iorb=1,Norb
                    i=i+1
                    bath_(i)=real(dmft_bath_%h(ispin,ispin,iorb,iorb,ibath))
                 enddo
              enddo
              !all off-diagonal non-vanishing terms in imploc - diagonal spin sectors
              do ispin=1,Maxspin
                 do iorb=1,Norb
                    do jorb=1,Norb
                       io = iorb + (ispin-1)*Norb
                       jo = jorb + (ispin-1)*Norb
                       if(io.gt.jo)cycle !only upper triangular saved. Hermiticity imposed in set_dmft_bath
                       if(dmft_bath_%mask(ispin,ispin,iorb,jorb,1)) then
                          i=i+1
                          bath_(i)=real(dmft_bath_%h(ispin,ispin,iorb,jorb,ibath))
                       endif
                       if(dmft_bath_%mask(ispin,ispin,iorb,jorb,2)) then
                          i=i+1
                          bath_(i)=aimag(dmft_bath_%h(ispin,ispin,iorb,jorb,ibath))
                       endif
                    enddo
                 enddo
              enddo
           enddo
           !
        endif
        !
        !all Re[Hybr]
        do ibath=1,Nbath
           i=i+1
           bath_(i)=real(dmft_bath_%t(ibath))
        enddo
        !
     case("nonsu2")
        !
        i = 0
        !
        if(replica_operators)then
           !
           do ibath=1,Nbath
              do iop=1,Nop
                 ![a, b, c]_k1  [a, b, c]_k2 [a, b, c]_k3 ...
                 i=i+1
                 bath_(i)=real(dmft_bath_%o(iop,ibath))
              enddo
           enddo
           !
        else
           !
           do ibath=1,Nbath
              !all diagonal non-vanishing terms in imploc - real only - always present
              do ispin=1,Nspin
                 do iorb=1,Norb
                    i=i+1
                    bath_(i)=real(dmft_bath_%h(ispin,ispin,iorb,iorb,ibath))
                 enddo
              enddo
              !all off-diagonal non-vanishing terms in imploc - all spin sectors
              do ispin=1,Nspin
                 do iorb=1,Norb
                    do jspin=1,Nspin
                       do jorb=1,Norb
                          io = iorb + (ispin-1)*Norb
                          jo = jorb + (jspin-1)*Norb
                          if(io.gt.jo)cycle !only upper triangular saved. Hermiticity imposed in set_dmft_bath
                          if(dmft_bath_%mask(ispin,jspin,iorb,jorb,1)) then
                             i=i+1
                             bath_(i)=real(dmft_bath_%h(ispin,jspin,iorb,jorb,ibath))
                          endif
                          if(dmft_bath_%mask(ispin,jspin,iorb,jorb,2)) then
                             i=i+1
                             bath_(i)=aimag(dmft_bath_%h(ispin,jspin,iorb,jorb,ibath))
                          endif
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           !
        endif
        !
        !all Re[Hybr]
        do ibath=1,Nbath
           i=i+1
           bath_(i)=real(dmft_bath_%t(ibath))
        enddo
        !
     case ("superc")
        stride = 0

        !
     end select
     !
  end select
end subroutine get_dmft_bath




!##################################################################
!
!     W_hyb PROCEDURES
!
!##################################################################
!+-----------------------------------------------------------------------------+!
!PURPOSE: build up the all-spin hybridization matrix W_{ss`}
!+-----------------------------------------------------------------------------+!
function get_Whyb_matrix_1orb(v,u) result(w)
  real(8),dimension(Nspin,Nbath)       :: v,u
  real(8),dimension(Nspin,Nspin,Nbath) :: w
  integer                              :: ispin
  !
  if(ed_para)then
     do ispin=1,Nspin
        w(ispin,ispin,:) = v(1,:)
     enddo
     w(1,Nspin,:) = u(1,:)
     w(Nspin,1,:) = u(1,:)
  else
     do ispin=1,Nspin
        w(ispin,ispin,:) = v(ispin,:)
     enddo
     w(1,Nspin,:) = u(1,:)
     w(Nspin,1,:) = u(2,:)
  endif
  !
end function get_Whyb_matrix_1orb

function get_Whyb_matrix_Aorb(v,u) result(w)
  real(8),dimension(Nspin,Norb,Nbath)       :: v,u
  real(8),dimension(Nspin,Nspin,Norb,Nbath) :: w
  integer                                   :: ispin
  !
  if(ed_para)then
     do ispin=1,Nspin
        w(ispin,ispin,:,:) = v(1,:,:)
     enddo
     w(1,Nspin,:,:) = u(1,:,:)
     w(Nspin,1,:,:) = u(1,:,:)
  else
     do ispin=1,Nspin
        w(ispin,ispin,:,:) = v(ispin,:,:)
     enddo
     w(1,Nspin,:,:) = u(1,:,:)
     w(Nspin,1,:,:) = u(2,:,:)
  endif
  !
end function get_Whyb_matrix_Aorb

function get_Whyb_matrix_dmft_bath(dmft_bath_) result(w)
  type(effective_bath)                      :: dmft_bath_
  real(8),dimension(Nspin,Nspin,Norb,Nbath) :: w
  integer                                   :: ispin
  !
  if(ed_para)then
     do ispin=1,Nspin
        w(ispin,ispin,:,:) = dmft_bath_%v(1,:,:)
     enddo
     w(1,Nspin,:,:) = dmft_bath_%u(1,:,:)
     w(Nspin,1,:,:) = dmft_bath_%u(1,:,:)
  else
     do ispin=1,Nspin
        w(ispin,ispin,:,:) = dmft_bath_%v(ispin,:,:)
     enddo
     w(1,Nspin,:,:) = dmft_bath_%u(1,:,:)
     w(Nspin,1,:,:) = dmft_bath_%u(2,:,:)
  endif
  !
end function get_Whyb_matrix_dmft_bath
