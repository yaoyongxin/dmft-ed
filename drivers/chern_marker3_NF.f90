!------------------------------------------------------------------
! This routine computes the local chern number in real space
! according to the paper by bianco & resta
! This routine is modified for a nanoflake, i.e a fully real space system
! Input are the projectors p and q of the occupied bands .
! get_local_chern is the unsymmetric version requiring less memory
! c(r) = -4 pi (2/v_uc) im tr (x_p y_q)
!------------------------------------------------------------------
subroutine get_local_chern(Hij,Smats)
  USE SCIFOR
  USE DMFT_TOOLS
  integer                                 :: Nlat,Nspin,Norb,Lmats,Nso
  real(8)                                 :: beta,xmu,lambda
  complex(8), intent(in)                  :: Hij(:,:,:)         ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nk==1]
  complex(8), intent(inout)               :: Smats(:,:,:,:,:,:) ![Nlat*Nspin*Nspin,Norb,Norb,Lreal][Nlat]
  complex(8), allocatable                 :: full_p_projector(:,:,:,:)
  complex(8), allocatable                 :: full_q_projector(:,:,:,:)
  complex(8), allocatable                 :: p_projector(:,:,:,:)
  complex(8), allocatable                 :: q_projector(:,:,:,:)
  !
  real(8), allocatable                    :: chern_marker(:,:)
  real(8), allocatable                    :: chern_marker_(:,:,:,:)
  real(8)                                 :: totalchern
  !
  real(8),allocatable                     :: Smats_rep(:,:)
  complex(8),allocatable                  :: MatTmp(:,:)
  real(8),allocatable                     :: RhoDiag(:)
  complex(8),allocatable                  :: RhoMat(:,:)
  real(8),allocatable                     :: dens_rho(:),dens_gf(:)
  integer                                 :: jx, jy, jy1, jy2, jr, js
  integer                                 :: ii,jj,iso1,iso2,iso
  real(8)                                 :: rx1,rx2,ry1,ry2
  real(8),allocatable,dimension(:)        :: r_i1,r_i2,ilat_vec,inequiv_vec
  integer                                 :: unitREAD,unitC,unit
  !
  ! These are for the projectors
  complex(8), allocatable, dimension(:,:) :: Hij_top, Evec, U,U_t
  complex(8), allocatable                 :: Gmats_rep(:,:,:)
  real(8), allocatable                    :: Eval(:),chern_loc(:)
  !
  integer                                 :: ilat, jlat, ispin, top_index, n_index
  integer                                 :: l_lat,k_lat, k_index,l_index
  complex(8), allocatable                 :: Gmats(:,:,:,:,:,:)
  complex(8), allocatable                 :: Emat(:,:,:,:)

  interface 
     function nn_reshape(MatNN,Nso,Nlat) result(Kmat)
       complex(8),dimension(Nso,Nso,Nlat,Nlat) :: MatNN
       integer                                 :: Nso,Nlat
       complex(8),dimension(Nso*Nlat,Nso*Nlat) :: Kmat
     end function nn_reshape
  end interface
  !
  interface
     function inv_nn_reshape(Kmat,Nso,Nlat) result(MatNN)
       complex(8),dimension(Nso*Nlat,Nso*Nlat) :: Kmat
       integer                                 :: Nso,Nlat
       complex(8),dimension(Nso,Nso,Nlat,Nlat) :: MatNN
     end function inv_nn_reshape
  end interface
  !
  !...> now we construct the matrix t (puts in --> full_q_projector)
  !
  Nlat  = size(Smats,1)
  Nspin = size(Smats,3)
  Norb  = size(Smats,5) ; if(Norb/=1)stop "This version is only for Norb=1"
  Lmats = size(Smats,6)
  Nso   = Nspin*Norb
  !
  call get_ctrl_var(beta,"BETA")
  call get_ctrl_var(xmu,"XMU")
  !
  !allocations
  allocate(full_p_projector(Nso,Nso,Nlat,Nlat))
  allocate(full_q_projector(Nso,Nso,Nlat,Nlat))
  allocate(p_projector(Nso,Nso,Nlat,Nlat))
  allocate(q_projector(Nso,Nso,Nlat,Nlat))
  allocate(chern_marker_(Nso,Nso,Nlat,Nlat),chern_marker(Nlat,Nlat))  
  allocate(MatTmp(Nso*Nlat,Nso*Nlat))
  allocate(r_i1(Nlat),r_i2(Nlat),ilat_vec(Nlat),inequiv_vec(Nlat))
  allocate(Hij_top(Nso*Nlat,Nso*Nlat))
  allocate(Evec(Nso*Nlat,Nso*Nlat))
  allocate(U(Nso*Nlat,Nso*Nlat),U_t(Nso*Nlat,Nso*Nlat))
  allocate(Smats_rep(Nso*Nlat,Nso*Nlat))
  allocate(Gmats_rep(Nso*Nlat,Nso*Nlat,Lmats))
  allocate(Eval(Nso*Nlat),chern_loc(Nlat))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats)) 
  allocate(Emat(Nspin,Nspin,Nlat,Nlat))
  allocate(RhoDiag(Nso*Nlat))
  allocate(RhoMat(Nso*Nlat,Nso*Nlat))
  allocate(dens_rho(Nso*Nlat))
  allocate(dens_gf(Nso*Nlat))

  ! -----------------------------------------------------------------------------------------------

  !In this part of the code, we compute the projectors, that are needed for the
  !chern marker
  p_projector = zero
  q_projector = zero
  !Here reshape
  !Smats(Nlat,Nspin,Nspin,Norb,Norb,1) => Smats_rep(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
  !This should be a block diagonal matrix where each block corresponds to a spin
  do ispin=1,Nspin
     do ilat=1,Nlat
        top_index =  ilat + (ispin-1)*Nlat
        Smats_rep(top_index,top_index) = dreal(Smats(ilat,ispin,ispin,Norb,Norb,1))
     enddo
  enddo
  !
  Hij_top = zero
  Eval = zero
  Evec = zero
  !
  !>TO BE REMOVED BEFORE FLIGHT
  Smats=zero
  !<
  Smats_rep = zero  
  Gmats_rep = zero
  dens_mat = 0d0
  f_mat = 0d0
  !  
  !Get the "Topological" Hamiltonian: i.e. the non-interacting H_ij plus static renormalization
  Hij_top = Hij(:,:,1) + Smats_rep
  Evec    = Hij_top
  !
  !Digonalize the Hij_top Hamiltonian:
  call eigh(Evec,Eval)
  !
  RhoDiag = fermi(Eval,beta)
  RhoMat  = matmul(Evec, matmul(diag(RhoDiag), transpose(conjg(Evec))) )
  !
  open(newunit=unit,file="ndens_rho.dat")
  do n_index = 1, Nspin*Nlat
     dens_rho(n_index) = dreal(RhoMat(n_index,n_index))
     write(unit,*)n_index,dens_rho(n_index)
  enddo
  close(unit)
  !
  !
  !
  ! This part is now to check if the diagonalization is correct  
  Gmats = zero
  call dmft_gloc_matsubara(Hij,[1d0],Gmats,Smats)
  !
  !Get densities
  open(newunit=unit,file="ndens_gf.dat")
  do ispin=1,Nspin
     do ilat = 1,Nlat
        n_index =  ilat + (ispin-1)*Nlat
        dens_gf(n_index) = fft_get_density(Gmats(ilat,ispin,ispin,Norb,Norb,:),beta)
        write(unit,*)n_index,dens_gf(n_index)
        write(*,*) "n_index, n_gf, n_rho:",n_index,dens_gf(n_index), dens_rho(n_index)                
     enddo
  enddo
  close(unit)
  !
  !Here we build up the P projector
  do ispin = 1, Nspin
     do ilat = 1,Nlat
        top_index =  ilat + (ispin-1)*Nlat
        lambda = real(Eval(top_index))
        if ( lambda  <  0d0 ) then
           Emat = zero
           do  l_lat = 1, Nlat
              do k_lat = 1, Nlat
                 l_index = l_lat + (ispin-1)*Nlat
                 k_index = k_lat * (ispin-1)*Nlat
                 Emat(ispin,ispin,k_lat,l_lat) =  Evec(l_index,top_index)*Evec(k_index,top_index)
              enddo
           enddo

           p_projector(ispin,ispin,:,:) = p_projector(ispin,ispin,:,:) + Emat(ispin,ispin,:,:)
           q_projector = q_projector

        elseif ( lambda  >=  0d0 ) then
           Emat = zero
           do  l_lat = 1, Nlat
              do k_lat = 1, Nlat
                 l_index = l_lat + (ispin-1)*Nlat
                 k_index = k_lat * (ispin-1)*Nlat
                 Emat(ispin,ispin,k_lat,l_lat) = Evec(l_index,top_index)*Evec(k_index,top_index)
              enddo
           enddo

           q_projector(ispin,ispin,:,:) = q_projector(ispin,ispin,:,:) + Emat(ispin,ispin,:,:)
           p_projector = p_projector
        endif
     enddo
  enddo



  ! .................................................................................................

  ! _________________________________________________________________________________________________

  ! In the following we compute the real space chern marker using the projecter
  ! computed above 

  open(free_unit(unitC),file="chern_marker.in")
  open(free_unit(unitREAD),file="RS_ave.in")

  do ilat = 1,Nlat
     read(unitREAD,*) r_i1(ilat), r_i2(ilat), ilat_vec(ilat), inequiv_vec(ilat)
  enddo

  write(*,*) "up to read its fine"
  full_p_projector = p_projector
  full_q_projector = q_projector
  rx1=0.0d0
  rx2=0.0d0
  ry1=0.0d0
  ry2=0.0d0


  do ilat = 1,Nlat
     rx1 = r_i1(ilat)
     ry1 = r_i2(ilat)
     do jlat =1,Nlat
        rx2 = r_i1(jlat)
        ry2 = r_i2(jlat)
        full_q_projector(:,:,ilat,jlat) = rx1*ry2*full_q_projector(:,:,ilat,jlat) !X1*Y2!
     enddo
  enddo
  !....>  only diagonal elements are needed
  !Ok this is very ugly, we may improve it
  full_q_projector =  inv_nn_reshape( &
       matmul (nn_reshape(full_p_projector,Nspin*Norb,Nlat) , &
       matmul( nn_reshape(full_q_projector,Nspin*Norb,Nlat) , nn_reshape(full_p_projector,Nspin*Norb,Nlat) ) ) &
       ,Nspin*Norb,Nlat)
  chern_marker=0d0
  chern_marker_=0d0
  chern_loc =0d0
  do jx = 1,Nlat
     do jy = 1,Nlat
        chern_marker_(:,:,jx,jy) = chern_marker_(:,:,jx,jy) - dimag( full_q_projector(:,:,jx,jx) )*2*pi2 !/ unit_cell_area
        chern_marker(jx,jy) = trace(chern_marker_(:,:,jx,jy))
     enddo
     !chern_loc(jx) = chern_loc(jx) + chern_marker(jx,jy)
     write(unitC,"(100F21.12)") r_i1(jx),r_i2(jy), r_i1(jy),r_i2(jx), chern_marker(jx,jy)
  enddo
  totalchern = zero
  do jr =1,Nlat
     totalchern = totalchern + dimag(trace(full_q_projector(:,:,jr,jr)))
  enddo
  ! there is an extra 1/n_sites=1/(Nlso*Nkx) in the totalchern
  totalchern = -totalchern*2*pi2/Nlat/Nspin*Norb
  write(*,*) "totalcherni=", totalchern 
  return
end subroutine get_local_chern

! _________________________________________________________________________________________________



! here are functions that are used for the computation of the chern marker

function nn_reshape(MatNN,Nso,Nlat) result(Kmat)
  complex(8),dimension(Nso,Nso,Nlat,Nlat) :: MatNN
  integer                                 :: Nso,Nlat
  complex(8),dimension(Nso*Nlat,Nso*Nlat) :: Kmat
  integer                                 :: iso,jso,i,j,ii,jj
  do i=1,Nlat
     do j=1,Nlat
        do iso=1,Nso
           do jso=1,Nso
              ii = iso + (i-1)*Nso
              jj = jso + (j-1)*Nso
              Kmat(ii,jj) = MatNN(iso,jso,i,j)
           enddo
        enddo
     enddo
  enddo
end function nn_reshape

function inv_nn_reshape(Kmat,Nso,Nlat) result(MatNN)
  complex(8),dimension(Nso*Nlat,Nso*Nlat) :: Kmat
  integer                                 :: Nso,Nlat
  complex(8),dimension(Nso,Nso,Nlat,Nlat) :: MatNN
  integer                                 :: iso,jso,i,j,ii,jj
  do i=1,Nlat
     do j=1,Nlat
        do iso=1,Nso
           do jso=1,Nso
              ii = iso + (i-1)*Nso
              jj = jso + (j-1)*Nso
              MatNN(iso,jso,i,j)  =  Kmat(ii,jj)
           enddo
        enddo
     enddo
  enddo
end function inv_nn_reshape




