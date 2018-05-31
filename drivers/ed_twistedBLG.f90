  program ed_twistedBLG
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                                       :: iloop,Lk,Nso,Nlso,Nlat
  logical                                       :: converged
  integer                                       :: ispin,ilat!,i,j

  !Bath:
  integer                                       :: Nb
  real(8),allocatable,dimension(:,:)            :: Bath,Bath_prev

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Greal

  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk
  complex(8),allocatable,dimension(:,:)         :: graphHloc,hkk
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  real(8),allocatable,dimension(:)              :: Wtk

  integer,allocatable,dimension(:)              :: ik2ix,ik2iy
  real(8),dimension(2)                          :: disp,disp2
  real(8),dimension(2)                          :: a1,a2,a1z,a2z,b1,b2,RR1,RR2,GG1,GG2
  real(8),dimension(2)                          :: bk1,bk2,bklen
  real(8),dimension(2)                          :: pointK1,pointK2

  !variables for the model:
  integer                                       :: Nk,Nkpath,m0,r
  real(8)                                       :: ang0,ang0r,ts,t0,t3,tsp,phi,delta,Mh,wmixing,cut,Vsi0,Vpi0,alat,r0,Lsc,Lmo
  character(len=32)                             :: finput
  character(len=32)                             :: hkfile
  logical                                       :: spinsym
  !
  real(8),dimension(2)                          :: Eout
  real(8),allocatable,dimension(:)              :: dens
  !
  real(8),dimension(:,:),allocatable            :: Zmats,Ucell
  complex(8),dimension(:,:,:),allocatable       :: Zfoo
  complex(8),allocatable,dimension(:,:,:,:,:)   :: S0


  call parse_cmd_variable(finput,"FINPUT",default='inputGRAPHENE.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS","inputGRAPHENE.conf",default=1d0)
  call parse_input_variable(t0,"T0","inputGRAPHENE.conf",default=0.142d0*ts)
  call parse_input_variable(t3,"T3","inputGRAPHENE.conf",default=0.107d0*ts)
  call parse_input_variable(mh,"MH","inputGRAPHENE.conf",default=0d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(m0,"m0",finput,default=1)
  call parse_input_variable(r,"r",finput,default=1)
  call parse_input_variable(cut,'cutoff',finput,default=sqrt(3.d0))
  call parse_input_variable(alat,"alat",finput,default=2.46d0)
  call parse_input_variable(r0,"r0",finput,default=0.184d0)
  call parse_input_variable(Vpi0,"Vpi0",finput,default=-2.7d0)
  call parse_input_variable(Vsi0,"Vsi0",finput,default=0.48d0)
  !
  !
  call ed_read_input(trim(finput))
 
   Nlat=4*(3*m0**2+3*m0*r+r**2)   ! Number of atoms in the unit cell

  if(Norb/=1)stop "Wrong setup from input file: Norb=1"



!:::::::::: SUPERLATTICE/LATTICE VECTORS GENERATION (see: PRB 92,075402 (2015) Sboychakov et. al.) ::::::::::::::::::::::::

   ang0=(dfloat(3*m0**2+3*m0*r)+(dfloat(r**2)/2.d0))/dfloat(3*m0**2+3*m0*r+r**2)  !twisting angle calculation:
   ang0r=acos(ang0)

    a1=alat*[sqrt(3.d0)/2.d0,-1.d0/2.d0]
    a2=alat*[sqrt(3.d0)/2.d0,1.d0/2.d0]
    a1z=a1*(cos(ang0r)-sin(ang0r)/sqrt(3.d0))+a2*(2.d0*sin(ang0r))/sqrt(3.d0)
    a2z=a2*(cos(ang0r)+sin(ang0r)/sqrt(3.d0))-a1*(2.d0*sin(ang0r))/sqrt(3.d0)
    disp2=[cos(ang0r)/sqrt(3.d0),sin(ang0r)/sqrt(3.d0)]*alat
    b1=(2.d0*pi/alat)*[1.d0/sqrt(3.d0),-1.d0] 
    b2=(2.d0*pi/alat)*[1.d0/sqrt(3.d0),1.d0]
    disp=(a1+a2)/3.d0
    cutoff=cut*alat
    Lsc=alat*sqrt(dfloat(Nlat))/2.d0
    Lmo=alat/(2.d0*sin(ang0r/2.d0))

    RR1=dfloat(m0)*a1+dfloat(m0+r)*a2
    RR2=-dfloat(m0+r)*a1+dfloat(2*m0+r)*a2
    GG1=(dfloat(2*m0+r)*b1+dfloat(m0+r)*b2)/(dfloat(3*m0**2+3*m0*r+r**2))
    GG2=(-dfloat(m0+r)*b1+dfloat(m0)*b2)/(dfloat(3*m0**2+3*m0*r+r**2))
    bk1=GG1
    bk2=GG2
    pointK1=(GG1+2.d0*GG2)/3.d0
    pointK2=(2.d0*GG1+GG2)/3.d0   


   call TB_set_bk(bk1,bk2)
   
   allocate(Ucell(Nlat,3))

  call build_uni_cell(alat,m0,r,Nlat,Ucell)

  Nso=Nspin*Norb
  Nlso=Nlat*Nso
  
 


  !Allocate Weiss Field:
  !allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  !allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  !allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  !allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  !allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
  !allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb));Hloc=zero
  !allocate(S0(Nlat,Nspin,Nspin,Norb,Norb));S0=zero
  !allocate(Zmats(Nlso,Nlso));Zmats=eye(Nlso)
  !allocate(Zfoo(Nlat,Nso,Nso));Zfoo=0d0
   

  !Buil the Hamiltonian on a grid or on  path
  call build_hk(trim(hkfile))
  Hloc = lso2nnn_reshape(graphHloc,Nlat,Nspin,Norb)


  stop






  !!Setup solver
  !Nb=get_bath_dimension()
  !allocate(Bath(Nlat,Nb))
  !allocate(Bath_prev(Nlat,Nb))
  !call ed_init_solver(Bath,Hloc)


  !!DMFT loop
  !iloop=0;converged=.false.
  !do while(.not.converged.AND.iloop<nloop)
  !   iloop=iloop+1
  !   call start_loop(iloop,nloop,"DMFT-loop")

  !   !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
  !   call ed_solve(Bath,Hloc)

  !   call ed_get_sigma_matsubara(Smats,Nlat)

  !   ! compute the local gf:
  !   call dmft_gloc_matsubara(Hk,Wtk,Gmats,Smats)

  !   ! compute the Weiss field (only the Nineq ones)
  !   if(cg_scheme=='weiss')then
  !      call dmft_weiss(Gmats,Smats,Weiss,Hloc)
  !   else
  !      call dmft_delta(Gmats,Smats,Weiss,Hloc)
  !   endif

  !   !Fit the new bath, starting from the old bath + the supplied Weiss
  !   call ed_chi2_fitgf(Bath,Weiss,Hloc,ispin=1)

  !   !MIXING:
  !   if(iloop>1)Bath=wmixing*Bath + (1.d0-wmixing)*Bath_prev
  !   Bath_prev=Bath

  !   converged = check_convergence(Weiss(:,1,1,1,1,:),dmft_error,nsuccess,nloop)

  !   call end_loop
  !enddo

  !call dmft_print_gf_matsubara(Smats,"Smats",iprint=1)
 
  !! extract and print retarded self-energy and Green's function 
  !call ed_get_sigma_real(Sreal,Nlat)
  !call dmft_print_gf_realaxis(Sreal,"Sreal",iprint=1)
  !call dmft_gloc_realaxis(Hk,Wtk,Greal,Sreal)
  !call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)


contains



  ! twisted BLG
  function hk_twistedBLG_model(kvec,Nlso) result(Hk)

    real(8),dimension(:)                :: kvec
    integer                             :: Nlso,ii,jj,i,j
    complex(8),dimension(Nlso,Nlso)     :: Hk
    real(8)                             :: ddd,Vpi,Vsi,ttt,dz,a0,cutoff,pi,r00
    real(8),dimension(3)                :: RR,dist,ez,kpoint
    real(8),dimension(:,:),allocatable           :: Ucellb
    complex(8)                          :: uni

    uni=(0.d0,1.d0)     !imaginary unit!
    pi=4.d0*datan(1.d0)   !pi
    dz=3.345d0*(alat/2.46d0) !displacement along the z axis of teh two planes
    a0=alat/sqrt(3.d0)  ! nearest neighbour distance
    r00=r0*alat ! decay lenght of the transfer integral (see: PRB 96,075311 (2017) Nam et. al.)  
    cutoff=cut*alat
    
    ez=[0.d0,0.d0,1.d0]
    kpoint(1)=kvec(1)
    kpoint(2)=kvec(2)
    kpoint(3)=0.d0

    allocate(Ucellb(Nlat,3))    


    do i=1,Nlat
      do j=1,Nlat
        Hk(i,j)=(0.d0,0.d0)
      end do
    end do



     do ii=-1,1
       do jj=-1,1
                
         RR=[RR1(1)*dfloat(ii)+RR2(1)*dfloat(jj),RR1(2)*dfloat(ii)+RR2(2)*dfloat(jj),0.d0]

          do j=1,3
            do i=1,Nlat

              Ucellb(i,j)=Ucell(i,j)+RR(j)

             end do
           end do

           do i=1,Nlat
             do j=i+1,Nlat

                dist(1)=Ucell(i,1)-Ucellb(j,1)
                dist(2)=Ucell(i,2)-Ucellb(j,2)
                dist(3)=Ucell(i,3)-Ucellb(j,3)

                ddd=sqrt(dist(1)**2.d0+dist(2)**2.d0+dist(3)**2.d0)
                              
                if (ddd.lt.(cutoff)) then

                  Vpi=Vpi0*exp(-(ddd-a0)/r00)
                  Vsi=Vsi0*exp(-(ddd-dz)/r00)

                  ttt=Vpi*(1.d0-(dot_product(dist,ez)/ddd)**2.d0)+Vsi*(dot_product(dist,ez)/ddd)**2.d0

                  Hk(i,j)=Hk(i,j)+ttt*exp(-uni*dot_product(kpoint,dist)+uni*dot_product(kpoint,RR))
                 end if
                 Hk(j,i)=conjg(Hk(i,j))
               end do                     
             end do
 
                  
      end do
    end do

    deallocate(Ucellb)

  end function hk_twistedBLG_model





  subroutine build_uni_cell(alat,m0,r,Nlat,Ucell)
    integer                            :: m0,r  !define number of atoms in unit cell
    integer                            :: Nlat !unit cell dimension
    integer                            :: i,j,k,m,n,ii,jj,iii,jjj,ab,zzz,ccc
    real(8),dimension(:,:),allocatable :: Ucell ![Nlat,3]
    real(8),dimension(:),allocatable :: ez,a1p,a2p,a1pz,a2pz ! z direction versor
    real(8)                             :: a0,r0,dz,ang0,ang0r,theta  ! parameters that describe the geometry of the lattice
    real(8)                          :: A,B,cut,alat  ! some other parameters... see below
    complex(8)                     ::uni !imaginary unit

    uni=(0.d0,1.d0)     !imaginary unit!
    dz=3.345d0*(alat/2.46d0) !displacement along the z axis of teh two planes
    a0=alat/sqrt(3.d0)  ! nearest neighbour distance
    r0=r0*alat ! decay lenght of the transfer integral (see: PRB 96,075311 (2017) Nam et. al.)  
    
    ang0=(dfloat(3*m0**2+3*m0*r)+(dfloat(r**2)/2.d0))/dfloat(3*m0**2+3*m0*r+r**2)  !twisting angle calculation:
    ang0r=acos(ang0)
    theta=(ang0r*180.d0/pi)
    
      write(*,*) 'the twisting angle is:', theta
      write(*,*) 'the Hamiltonian is a', Nlat,'x',Nlat,'matrix'


     allocate(a1p(2),a2p(2),a1pz(2),a2pz(2),ez(3))

     ez=[0.d0,0.d0,1.d0]



!          :::::::::::::::::::   UNIT CELL GENERATION     :::::::::::::::::


!     FIRST LAYER (zzz=0)

     zzz=0
     ccc=0

     do i=-Nlat,Nlat
       do j=-Nlat,Nlat
         do jj=1,2

              m=i
              n=j
              ab=jj-1

              if (jj.eq.2) then
                a1p=a1*dfloat(m)+disp
                a2p=a2*dfloat(n)
                A=dot_product(GG1,a1p)+dot_product(GG1,a2p)
                B=dot_product(GG2,a1p)+dot_product(GG2,a2p)
              else
                A=dot_product(GG1,a1)*dfloat(m)+dot_product(GG1,a2)*dfloat(n)
                B=dot_product(GG2,a1)*dfloat(m)+dot_product(GG2,a2)*dfloat(n)
              end if

              if(A.ge.0.d0.and.A.lt.(2.d0*pi-0.0000001d0)) then       ! 2*pi minus a small number----> to avoid to count the first atom of the next unit cell
                if (B.ge.0.d0.and.B.lt.(2.d0*pi-0.0000001d0)) then

                   ccc=ccc+1

                   Ucell(ccc,1)=dfloat(m)*a1(1)+dfloat(n)*a2(1)+disp(1)*dfloat(ab)
                   Ucell(ccc,2)=dfloat(m)*a1(2)+dfloat(n)*a2(2)+disp(2)*dfloat(ab)
                   Ucell(ccc,3)=dfloat(zzz)*dz

                end if
              end if

         end do 
       end do
     end do


!     SECOND LAYER (zzz=1)

    zzz=1
           
    do i=-Nlat,Nlat
      do j=-Nlat,Nlat
        do jj=1,2

           m=i
           n=j
           ab=jj-1

           if (jj.eq.2) then
             a1pz=a1z*dfloat(m)-disp2
             a2pz=a2z*dfloat(n)
             !
             A=dot_product(GG1,a1pz)+dot_product(GG1,a2pz)
             B=dot_product(GG2,a1pz)+dot_product(GG2,a2pz)
           else
             A=dot_product(GG1,a1z)*dfloat(m)+dot_product(GG1,a2z)*dfloat(n)
             B=dot_product(GG2,a1z)*dfloat(m)+dot_product(GG2,a2z)*dfloat(n)
           end if

           if(A.ge.0.d0.and.A.lt.(2.d0*pi-0.0000001d0)) then
             if (B.ge.0.d0.and.B.lt.(2.d0*pi-0.0000001d0)) then
               ccc=ccc+1
               if (ccc.gt.Nlat) then
                 write(*,*) 'ERROR: wrong unit cell generation!!!'
               end if
              
               Ucell(ccc,1)=dfloat(m)*a1z(1)+dfloat(n)*a2z(1)-disp2(1)*dfloat(ab)
               Ucell(ccc,2)=dfloat(m)*a1z(2)+dfloat(n)*a2z(2)-disp2(2)*dfloat(ab)
               Ucell(ccc,3)=dfloat(zzz)*dz

              end if
            end if

        end do 
      end do
    end do

    open(1000,file='Unit_cell.dat')
    do i=1,Nlat
      write(1000,*) Ucell(i,:)
    end do
    close(1000)

    if (ccc.lt.Nlat) then
       write(*,*) 'ERROR:  the unit cell is not completed', Nlat-ccc,'atoms are missing!!!'
       stop
    end if

  end subroutine








  !---------------------------------------------------------------------
  !PURPOSE: Get BL Graphene Model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional                              :: file
    integer                                                :: i,j,ik
    integer                                                :: ix,iy
    real(8)                                                :: kx,ky  
    integer                                                :: iorb,jorb
    integer                                                :: isporb,jsporb
    integer                                                :: ispin,jspin
    integer                                                :: unit
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: Gmats,fooSmats
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: Greal,fooSreal
    real(8),dimension(2)                                   :: kvec
    real(8)                                                :: blen,area_hex,area_rect,points_in,points_tot
    real(8),allocatable,dimension(:)                       :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable                     :: KPath
    type(rgb_color),dimension(:),allocatable               :: colors
!    complex(8),dimension(Nlso,Nlso),allocatable            :: hkk

    Lk= Nk*Nk

    write(LOGfile,*)"Build H(k) of twisted Bilayer Graphene:",Lk
    write(LOGfile,*)"# of SO-bands     :",Nlso

    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nlso,Nlso,Lk));Hk=zero
    allocate(wtk(Lk));Wtk=0d0
    allocate(kxgrid(Nk),kygrid(Nk))

    ik=0
    do iy=1,Nk
       ky = dble(iy-1)/Nk
       do ix=1,Nk
          ik=ik+1
          kx=dble(ix-1)/Nk
          kvec = kx*bk1 + ky*bk2
          kxgrid(ix) = kvec(1)
          kygrid(iy) = kvec(2)
          Hk(:,:,ik)=hk_twistedBLG_model(kvec,Nlso)
     enddo
    enddo
    Wtk = 1d0/Lk


    if(present(file))then
       call TB_write_hk(Hk,"Hkrfile_BLG_twisted.data",&
            No=Nlso,&
            Nd=Norb,&
            Np=0,&
            Nineq=1,&
            Nkvec=[Nk,Nk])
    endif
    !
    allocate(graphHloc(Nlso,Nlso))
    graphHloc = sum(Hk(:,:,:),dim=3)/Lk
    where(abs(dreal(graphHloc))<1.d-4)graphHloc=0d0
    call TB_write_Hloc(graphHloc)
    call TB_write_Hloc(graphHloc,'Hloc.txt')
    !
    !


    write(*,*) 'Nlso=',Nlso
    allocate(Kpath(4,2))
    KPath(1,:)=[0,0]
    KPath(2,:)=pointK1
    KPath(3,:)=pointK2
    KPath(4,:)=[0,0]

    allocate(colors(Nlso))
    colors = gray88

    call TB_Solve_model(hk_twistedBLG_model,Nlso,KPath,Nkpath,&
         colors_name=colors,&
         points_name=[character(len=10) :: "G","K","Kp","G"],&
         file="Eigenbands.nint")



    !Build the local GF:
    !Gmats=zero
    !Greal=zero
    !fooSmats =zero
    !fooSreal =zero
    !call add_ctrl_var(beta,"BETA")
    !call add_ctrl_var(xmu,"xmu")
    !call add_ctrl_var(wini,"wini")
    !call add_ctrl_var(wfin,"wfin")
    !call add_ctrl_var(eps,"eps")
    !call dmft_gloc_matsubara(Hk,Wtk,Gmats,fooSmats)
    !call dmft_print_gf_matsubara(Gmats,"LG0",iprint=1)
    !call dmft_gloc_realaxis(Hk,Wtk,Greal,fooSreal)
    !call dmft_print_gf_realaxis(Greal,"LG0",iprint=1)
    !
  end subroutine build_hk








  end program ed_twistedBLG



