!----------------------------------------------------------------------------------------!
  ! purpose: evaluate the non-local bare static spin susceptibility
  ! given the non-local Green's function and the fermi distribution on
  ! the real axis. 
  ! chi0_ij = 1/pi Im \int_{-infty}^{infty} G_ij(w) G_ji(w) f(w) dw
  !----------------------------------------------------------------------------------------!
subroutine ed_get_chi0ij(Gijmats,Nvec,e1,e2,e3)
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  complex(8),intent(in)                     :: Gijmats(:,:,:,:,:,:,:)![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
  integer,dimension(:),intent(in)           :: Nvec  ![Nx{,Ny,Nz}]
  real(8),dimension(size(Nvec)),optional    :: e1,e2,e3

  real(8),dimension(3)                      :: e1_=[1d0,0d0,0d0]
  real(8),dimension(3)                      :: e2_=[0d0,1d0,0d0]
  real(8),dimension(3)                      :: e3_=[0d0,0d0,1d0]
  !
  integer                                   :: Nlat, Lmats, Nspin, Norb, Ndim, Nk
  real(8)                                   :: beta, k_times_R, kq_times_R
  integer                                   :: Nso,ik,iq,iorb,ilat,jlat,iw,unitUP,unitDW
  real(8),allocatable,dimension(:,:)        :: Rgrid, kgrid, qgrid
  complex(8)                                :: c_n1, c_n2
  real(8)                                   :: kpoint(size(Nvec)),qpoint(size(Nvec)),kqpoint(size(Nvec)), r_i(size(Nvec)), r_j(size(Nvec))
  complex(8),allocatable,dimension(:,:,:,:) :: G_up_kq,G_up_k ![Nk][Norb][Norb][Lmats]
  complex(8),allocatable,dimension(:,:,:,:) :: G_dw_kq,G_dw_k ![Nk][Norb][Norb][Lmats]
  complex(8),allocatable                    :: chi0_up(:,:,:)
  complex(8),allocatable                    :: chi0_dw(:,:,:)
  !
  Nlat  = size(Gijmats,1)
  Nspin = size(Gijmats,3)
  Norb  = size(Gijmats,5) 
  Lmats = size(Gijmats,7)
  Ndim  = size(Nvec)    
  call get_ctrl_var(beta,"BETA")
  !
  Nso  = Nspin*Norb
  Nk = Nlat
  !
  call assert_shape(Gijmats,[Nlat,product(Nvec),Nspin,Nspin,Norb,Norb,Lmats],"ed_get_chi0ij","GijMats")
  !
  write(*,*) "Nlat=", Nlat, "Nvec=", Nvec
  !
  !> Set the direct lattice basis to the passed vectors:
  if(present(e1))e1_(:Ndim)=e1
  if(present(e2))e2_(:Ndim)=e2
  if(present(e3))e3_(:Ndim)=e3
  call TB_set_ei(e1_,e2_,e3_)
  call TB_print_ei()

  !> Compute the reciprocal basis, given the direct basis:
  call TB_build_bk(verbose=.true.)


  write(6,*) "now we are at the subroutin for chi0"



  allocate(kgrid(Nlat,Ndim))
  allocate(Rgrid(Nlat,Ndim))
  !> Build the reciprocal and direct lattice grids:
  call TB_build_kgrid(Nvec,kgrid)
  call TB_build_rgrid(Nvec,Rgrid)

  !open(free_unit(unit),file="chi0_latt_nano.ed")
  allocate(G_dw_k(Nk,Norb,Norb,Lmats))
  allocate(G_up_k(Nk,Norb,Norb,Lmats))
  allocate(G_dw_kq(Nk,Norb,Norb,Lmats))
  allocate(G_up_kq(Nk,Norb,Norb,Lmats))
  allocate(chi0_dw(Norb,Norb,Nk))
  allocate(chi0_up(Norb,Norb,Nk))

  G_dw_k  = zero
  G_up_k  = zero
  G_dw_kq = zero
  G_up_kq = zero
  chi0_up = zero
  chi0_dw = zero

  ! do ilat = 1,Nlat
  !    write(*,*) "ilat=", ilat
  !    write(*,*) "here comes momentum"
  !    write(*,*) kgrid(ilat,1), kgrid(ilat,2)
  !    write(*,*) "here comes position"
  !    write(*,*) Rgrid(ilat,1), Rgrid(ilat,2)
  ! enddo



  open(free_Unit(unitUP),file="Chi0_s1_q.dat")
  open(free_Unit(unitDW),file="Chi0_s2_q.dat")
  call start_timer()
  iq_loop: do iq = 1,Nlat
     qpoint = kgrid(iq,:)
     write(999,*) "iq=", iq     
     write(999,*) qpoint(1),qpoint(2)
     !
     ik_loop: do ik = 1,Nk
        kpoint = kgrid(ik,:)
        write(998,*) "ik=", ik        
        write(998,*) kpoint(1),kpoint(2)
        kqpoint = kpoint + qpoint
        IF (kqpoint(1) > 2*pi) kqpoint(1) = kqpoint(1) - 2*pi
        IF (kqpoint(2) > 2*pi) kqpoint(2) = kqpoint(2) - 2*pi
        !
        !Fourier transform of Gij to G_q
        ilat_loop: do ilat=1,Nlat
           R_i = Rgrid(ilat,:)
           !
           jlat_loop: do jlat=1,Nlat
              R_j = Rgrid(jlat,:)
              !
              k_times_R = dot_product(kpoint , R_i-R_j)
              c_n1 =exp(-xi*k_times_R)!exp(-i_cmplx*(kpoint(1)*(r_i(1)-r_j(1))+kpoint(2)*(r_i(2)-r_j(2))))
              G_up_k(ik,:,:,:) = G_up_k(ik,:,:,:) + Gijmats(ilat,jlat,1,1,:,:,:)*c_n1
              G_dw_k(ik,:,:,:) = G_dw_k(ik,:,:,:) + Gijmats(ilat,jlat,Nspin,Nspin,:,:,:)*c_n1
              !
              kq_times_R = dot_product(kqpoint, R_j - R_i)
              c_n2=exp(xi*kq_times_R)!exp(i_cmplx*(kqpoint(1)*(r_j(1)-r_i(1))+kqpoint(2)*(r_j(2)-r_i(2))))
              G_up_kq(ik,:,:,:) = G_up_kq(ik,:,:,:) + Gijmats(jlat,ilat,1,1,:,:,:)*c_n2
              G_dw_kq(ik,:,:,:) = G_dw_kq(ik,:,:,:) + Gijmats(jlat,ilat,Nspin,Nspin,:,:,:)*c_n2
           enddo jlat_loop
        enddo ilat_loop
        !
        G_up_k(ik,:,:,:)  = G_up_k(ik,:,:,:)/dble(Nlat)
        G_up_kq(ik,:,:,:) = G_up_kq(ik,:,:,:)/dble(Nlat)
        G_dw_k(ik,:,:,:)  = G_dw_k(ik,:,:,:)/dble(Nlat)
        G_dw_kq(ik,:,:,:) = G_dw_kq(ik,:,:,:)/dble(Nlat)
        !
        do iw = 1,Lmats
           chi0_up(:,:,iq) = chi0_up(:,:,iq) + matmul( G_up_k(ik,:,:,iw) , G_up_kq(ik,:,:,iw) )  !g_up_k(ik,iw)*G_up_kq(ik,iw)
           chi0_dw(:,:,iq) = chi0_dw(:,:,iq) + matmul( G_dw_k(ik,:,:,iw) , G_dw_kq(ik,:,:,iw) ) 
        enddo
     enddo ik_loop
     !
     chi0_up(:,:,iq) = -2d0*chi0_up(:,:,iq)/beta/dble(Nk) !the factor 2 comes from negative Matsubara frequencies
     chi0_dw(:,:,iq) = -2d0*chi0_dw(:,:,iq)/beta/dble(Nk)
     
     
     write(unitUP,"(100F21.12)")qpoint(1),qpoint(2),(dreal(chi0_up(iorb,iorb,iq)),dimag(chi0_up(iorb,iorb,iq)),iorb=1,Norb) !+dreal(chi0_dw(iq))
     write(unitDW,"(100F21.12)")qpoint(1),qpoint(2),(dreal(chi0_dw(iorb,iorb,iq)),dimag(chi0_dw(iorb,iorb,iq)),iorb=1,Norb) !+dreal(chi0_dw(iq))
     call eta(iq,Nk)
  enddo iq_loop
  call stop_timer()
  close(unitUP)
  close(unitDW)
end subroutine ed_get_chi0ij



