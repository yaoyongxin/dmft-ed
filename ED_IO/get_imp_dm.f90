subroutine ed_get_density_matrix_single(dm_,dm_eig_,dm_rot_)
  implicit none
  !passed
  complex(8),allocatable,intent(out)           :: dm_(:,:)
  complex(8),allocatable,intent(out),optional  :: dm_rot_(:,:)
  real(8),allocatable,intent(out),optional     :: dm_eig_(:)
  !internal
  logical                                      :: diag_dm
  !
  if (.not.allocated(imp_density_matrix)) then
     write(LOGfile,"(A)") "imp_density_matrix is not allocated"
     stop
  endif
  !
  diag_dm=present(dm_eig_).and.present(dm_rot_)
  !
  if(allocated(dm_))deallocate(dm_);allocate(dm_(Nspin*Norb,Nspin*Norb));dm_ = zero
  if (diag_dm) then
     if(allocated(dm_eig_))deallocate(dm_eig_);allocate(dm_eig_(Nspin*Norb))           ; dm_eig_ = 0.0d0
     if(allocated(dm_rot_))deallocate(dm_rot_);allocate(dm_rot_(Nspin*Norb,Nspin*Norb)); dm_rot_ = zero
  endif
  !
  dm_ = nn2so_reshape(imp_density_matrix,Nspin,Norb)
  !
  if(bath_type=="normal")then
     !
     if(diag_dm)then
        dm_rot_=eye(Nspin*Norb)
        dm_eig_=diagonal(dm_)
        call ed_print_density_matrix(dm_,dm_rot_,dm_eig_)
     else
        call ed_print_density_matrix(dm_)
     endif
     !
  else
     !
     if(diag_dm)then
        dm_rot_=dm_
        call eigh(dm_rot_,dm_eig_,'V','U')
        call ed_print_density_matrix(dm_,dm_rot_,dm_eig_)
    else
       call ed_print_density_matrix(dm_)
    endif
    !
  endif
  !
end subroutine ed_get_density_matrix_single



subroutine ed_get_density_matrix_lattice(dm_,dm_eig_,dm_rot_)
  implicit none
  !passed
  complex(8),allocatable,intent(out)           :: dm_(:,:,:)
  complex(8),allocatable,intent(out),optional  :: dm_rot_(:,:,:)
  real(8),allocatable,intent(out),optional     :: dm_eig_(:,:)
  !internal
  logical                                      :: diag_dm
  integer                                      :: ilat,Nlat
  real(8),allocatable                          :: dm_eig_site(:)
  complex(8),allocatable                       :: dm_site(:,:),dm_rot_site(:,:)
  !
  Nlat=size(imp_density_matrix_ii,1)
  !
  if (.not.allocated(imp_density_matrix)) then
     write(LOGfile,"(A)") "imp_density_matrix is not allocated"
     stop
  endif
  !
  diag_dm=present(dm_eig_).and.present(dm_rot_)
  !
  if(allocated(dm_))deallocate(dm_);allocate(dm_(Nlat,Nspin*Norb,Nspin*Norb));dm_ = zero
  allocate(dm_site(Nspin*Norb,Nspin*Norb));dm_site = zero
  !
  if (diag_dm) then
     if(allocated(dm_eig_))deallocate(dm_eig_);allocate(dm_eig_(Nlat,Nspin*Norb))           ; dm_eig_ = 0.0d0
     if(allocated(dm_rot_))deallocate(dm_rot_);allocate(dm_rot_(Nlat,Nspin*Norb,Nspin*Norb)); dm_rot_ = zero
     allocate(dm_eig_site(Nspin*Norb));dm_eig_site = 0.0d0
     allocate(dm_rot_site(Nspin*Norb,Nspin*Norb));dm_rot_site = zero
  endif
  !
  do ilat=1,Nlat
     !
     dm_site = nn2so_reshape(imp_density_matrix_ii(ilat,:,:,:,:),Nspin,Norb)
     !
     if(bath_type=="normal")then
        !
        if(diag_dm)then
           dm_rot_site=eye(Nspin*Norb)
           dm_eig_site=diagonal(dm_site)
           call ed_print_density_matrix(dm_site,dm_rot_site,dm_eig_site,suff=str(ilat))
        else
           call ed_print_density_matrix(dm_site,suff=str(ilat))
        endif
        !
     else
        !
        if(diag_dm)then
           dm_rot_site=dm_site
           call eigh(dm_rot_site,dm_eig_site,'V','U')
           call ed_print_density_matrix(dm_site,dm_rot_site,dm_eig_site,suff=str(ilat))
       else
          call ed_print_density_matrix(dm_site,suff=str(ilat))
       endif
       !
     endif
     !
     dm_(ilat,:,:)=dm_site
     dm_rot_(ilat,:,:)=dm_rot_site
     dm_eig_(ilat,:)=dm_eig_site
     !
  enddo
  !
end subroutine ed_get_density_matrix_lattice


subroutine ed_print_density_matrix(dm_,dm_rot_,dm_eig_,suff)
  implicit none

  complex(8),allocatable,intent(in)            :: dm_(:,:)
  complex(8),allocatable,intent(in),optional   :: dm_rot_(:,:)
  real(8),allocatable   ,intent(in),optional   :: dm_eig_(:)
  character(len=*)      ,intent(in),optional   :: suff
  !internal
  integer                                      :: unit
  character(len=24)                            :: suffix
  integer                                      :: iorb,jorb,ispin,jspin,io,jo
  !
  suffix="imp_density_matrix.dat"
  if(present(suff))suffix="imp_density_matrix_"//reg(suff)//".dat"
  !
  unit = free_unit()
  open(unit,file=suffix,action="write",position="rewind",status='unknown')
  !
  write(unit,"(A90)")"# density matrix in the impurity problem basis REAL part:"
  do io=1,Nspin*Norb
     write(unit,"(90(F15.9,1X))") (real(dm_(io,jo)),jo=1,Nspin*Norb)
  enddo
  write(unit,*)
  !
  write(unit,"(A90)")"# density matrix in the impurity problem basis IMAGINARY part:"
  do io=1,Nspin*Norb
     write(unit,"(90(F15.9,1X))") (aimag(dm_(io,jo)),jo=1,Nspin*Norb)
  enddo
  write(unit,*)
  !
  if(present(dm_eig_).and.present(dm_rot_))then
     write(unit,"(A90)")"# eigenvalues of density matrix"
     write(unit,'(10F22.12)') dm_eig_
     write(unit,*)
     !
     write(unit,"(A90)")"# density matrix eigenvector matrix REAL part:"
     do io=1,Nspin*Norb
        write(unit,"(90(F15.9,1X))") (real(dm_rot_(io,jo)),jo=1,Nspin*Norb)
     enddo
     write(unit,*)
     !
     write(unit,"(A90)")"# density matrix eigenvector matrix IMAGINARY part:"
     do io=1,Nspin*Norb
        write(unit,"(90(F15.9,1X))") (aimag(dm_rot_(io,jo)),jo=1,Nspin*Norb)
     enddo
     write(unit,*)
  endif
  !
  close(unit)
  !
end subroutine ed_print_density_matrix
