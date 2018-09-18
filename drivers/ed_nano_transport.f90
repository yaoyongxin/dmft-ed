program trans
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                                         :: Nk,Nlso,Nineq,Nlat
  integer                                         :: ilat,ineq,ispin,iorb
  !
  integer,dimension(:),allocatable                :: lat2ineq,ineq2lat
  !
  integer,dimension(:),allocatable                :: sb_field_sign
  !
  ! Green's functions
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Greal           ![Nlat][Nspin]{Nspin][Norb][Norb][Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Sreal           ![Nlat][Nspin]{Nspin][Norb][Norb][Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Sreal_ineq      ![Nineq][Nspin]{Nspin][Norb][Norb][Lreal]
  !
  !hamiltonian input:
  complex(8),allocatable                          :: Hij(:,:,:)      ![Nlso][Nlso][Nk=1]
  complex(8),allocatable                          :: nanoHloc(:,:)   ![Nlso][Nlso]
  !
  !input files:
  character(len=32)                               :: finput
  character(len=32)                               :: nfile,hijfile
  !
  logical                                         :: jbias,read_Hij


  call parse_cmd_variable(finput,"FINPUT",default='inputED_NANO.conf')
  call parse_input_variable(nfile,"NFILE",finput,default="nano.in")
  call parse_input_variable(hijfile,"HIJFILE",finput,default="hij.in")
  call parse_input_variable(jbias,"jbias",finput,default=.false.)
  call parse_input_variable(read_Hij,"readHij",finput,default=.false.)

  ! read input
  call ed_read_input(trim(finput))

  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  ! set input structure hamiltonian: allocates and sets Hij & nanoHloc
  !
  ! instead of building Hij, implement also a readin option
  ! format: wannier90-like (also compatible with build_Hij)
  ! in that case the input file would also be much "slimmer"
  ! note that many integers (e.g., Nlat, Nlso, ...) are set within build_Hij!
  if(read_Hij)then
     write(*,*)"error: readin option for Hij not implemented yet!"
  else
     call build_Hij([nfile,hijfile])
  endif

  ! allocate self-energy
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  !
  ! allocate Green's function
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))

  ! read converged self-energy
  call read_sigma_real(Sreal_ineq)
  do ilat=1,Nlat
     ineq = lat2ineq(ilat)
     Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
  enddo
  ! 
  ! get local Green's function as a sanity check
  call dmft_gloc_realaxis(Hij,[1d0],Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"LG",iprint=1)
  !
  ! evaluate the linear response (zero-bias) transmission function 
  ! if jbias=T evaluate the corresponding bias-driven current
  call ed_transport(Hij,Sreal)



contains



  !----------------------------------------------------------------------------------------!
  ! purpose: build real-space Hamiltonian for a nanostructure of size [Nlat*Nspin*Norb]**2
  !----------------------------------------------------------------------------------------!
  subroutine build_Hij(file)
    character(len=*)     :: file(2)
    integer              :: ilat,jlat,iorb,jorb,is,js,ispin,ie
    integer              :: i,isite,iineq,iineq0,isign
    integer              :: EOF
    character, parameter :: tab = achar ( 9 )
    integer              :: unit,ineq_count
    integer              :: Ns,Ne,Nb!,Nk         ! #atoms, #inequivalent, #bands
    real(8)              :: ret,imt
    logical              :: blank_at_right
    character(len=1)     :: next,prev
    character(len=6)     :: site,sign
    write(LOGfile,*)"Build H(R_i,R_j) for a NANO object:"
    ! readin generic input
    ! allocate & fill inequivalent list
    unit = free_unit()
    open(unit,file=trim(file(1)),status='old')
    read(unit,*)Ns,Ne,Nb
    !Checks:
    if(Nb/=Norb)stop "build_Hij error: Nb read from file != Norb in input.conf"
    Nk   = 1
    Nb   = Norb
    Nlat = Ns
    Nineq= Ne
    Nlso = Nlat*Nspin*Norb
    allocate(lat2ineq(Nlat),ineq2lat(Nineq))
    read(unit,"(A1)",advance='no',IOSTAT=EOF)next
    site  = next
    isite = 0
    i     = 0
    do 
       prev=next
       read(unit,"(A1)",advance='no',IOSTAT=EOF)next
       blank_at_right = ((prev/=' '.AND.prev/=tab).AND.(next==' '.OR.next==tab))
       if(.not.blank_at_right)then
          site=trim(site)//next
       else
          read(site,"(I6)")isite
          site=""
          i=i+1
          if(i>Nlat)stop "build_Hij error: lattice index > Nlat read from file"
          lat2ineq(i)=isite+1
       endif
       if(EOF<0)exit
    enddo
    if(i<Nlat)stop "build_Hij error: lattice index < Nlat read from file"
    write(*,*)"# of sites      :",Nlat
    write(*,*)"# of ineq sites :",Nineq
    write(*,*)"# of bands      :",Norb
    !
    ineq_count=1
    iineq=lat2ineq(Nlat)
    do i=Nlat,2,-1
       iineq0=lat2ineq(i-1)!iineq
       iineq =lat2ineq(i)
       if(iineq/=iineq0)then
          ineq2lat(iineq)=i
          ineq_count=ineq_count+1
       endif
       !if(ineq_count==Nineq)exit
    enddo
    iineq=lat2ineq(1)
    ineq2lat(1)=iineq
    !close(unit) ! do not close unit if readin info below
    !
    ! allocate & fill sign list of symmetry-breaking field
    allocate(sb_field_sign(Nineq))
    sign  = next
    isign = 0
    i     = 0
    do 
       prev=next
       read(unit,"(A1)",advance='no',IOSTAT=EOF)next
       blank_at_right = ((prev/=' '.AND.prev/=tab).AND.(next==' '.OR.next==tab))
       if(.not.blank_at_right)then
          sign=trim(sign)//next
       else
          read(sign,"(I6)")isign
          sign=""
          i=i+1
          if(i>Nineq)stop "build_Hij error: lattice index > Nineq read from file"
          sb_field_sign(i)=isign
       endif
       if(EOF<0)exit
    enddo
    close(unit)
    !
    ! allocate and initialize H(r_i,r_j)
    allocate(Hij(Nlso,Nlso,Nk))
    Hij = zero 
    unit = free_unit()
    open(unit,file=trim(file(2)),status='old')
    do !while(EOF>=0)
       read(unit,*,IOSTAT=EOF)ilat,iorb,jlat,jorb,ret,imt
       ilat=ilat+1
       iorb=iorb+1
       jlat=jlat+1
       jorb=jorb+1
       if(EOF<0)exit
       do ispin=1,Nspin
          is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
          js = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb
          ! symmetric hopping
          Hij(is,js,1)=dcmplx(ret,imt) 
          Hij(js,is,1)=dcmplx(ret,imt) ! symmetrize hopping
       enddo
    enddo
    close(unit)
    call assert_shape(Hij,[Nlso,Nlso,Nk],'build_Hij',"Hij")
    !
    ! basis vectors must be defined
    call TB_set_bk([1d0,0d0,0d0],[0d0,1d0,0d0],[0d0,0d0,1d0])
    call TB_write_hk(Hk=Hij,file="Hij_nano.data",&
         No=Nlso,&
         Nd=Norb,&
         Np=0,&
         Nineq=Nineq,&
         Nkvec=[1,1,1])
    !
    allocate(nanoHloc(Nlso,Nlso))
    nanoHloc = extract_Hloc(Hij,Nlat,Nspin,Norb)
    !
    !save lat2ineq,ineq2lat arrays
    unit=free_unit()
    open(unit,file="lat2ineq.ed")
    do ilat=1,Nlat
       write(unit,*)ilat,lat2ineq(ilat)
    enddo
    close(unit)
    unit=free_unit()
    open(unit,file="ineq2lat.ed")
    do i=1,Nineq
       write(unit,*)i,ineq2lat(i)
    enddo
    close(unit)
  end subroutine build_Hij



  !----------------------------------------------------------------------------------------!
  ! purpose: evaluate 
  !  - conductance (without vertex corrections) 
  !  - bias-driven current
  ! for a nanostructure on the real axis, given the non-local Green's function 
  ! and the L/R hybridization matrix, of size [Nlat*Nspin*Norb**2*Lreal]
  !----------------------------------------------------------------------------------------!
  subroutine ed_transport(Hij,Sreal)
    ! inputs: Hamiltonian and retarded self-energy
    complex(8),intent(in)                         :: Hij(:,:,:)           ![Nlso][Nlso][Nk=1]
    complex(8),intent(in)                         :: Sreal(:,:,:,:,:,:)   ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    ! retarded Green's function
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Greal(:,:,:,:,:,:)   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    ! hybridization function to environment
    complex(8),dimension(:,:),allocatable         :: self_lead            ![Nlso][Nlso]
    !
    complex(8),dimension(:,:),allocatable         :: GR,HR,GA,HL,Re,Le,Te ![Nlat*Norb][Nlat*Norb]
    integer,dimension(:,:),allocatable            :: rmask,lmask          ![Nlat*Norb][Nlat*Norb]
    !
    real(8),dimension(:),allocatable              :: wr
    !
    complex(8),dimension(:,:),allocatable         :: transe               ![Nspin][Lreal]
    real(8),dimension(:),allocatable              :: jcurr                ![Nspin]
    !
    real(8)                                       :: lbias,rbias
    !
    integer                                       :: ilat,jlat,ispin,jspin,iorb,jorb,io,jo,is,js,i,Nlso,Nlo
    integer                                       :: unit,unit_in,unit_out,eof,lfile
    character(len=30)                             :: suffix
    !
    Nlso=Nlat*Nspin*Norb
    Nlo=Nlat*Norb
    !
    allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)

    ! allocate variables for matrix-matrix multiplication
    allocate(GR(Nlo,Nlo));GR=zero
    allocate(HR(Nlo,Nlo));HR=zero
    allocate(GA(Nlo,Nlo));GA=zero
    allocate(HL(Nlo,Nlo));HL=zero
    allocate(Re(Nlo,Nlo));Re=zero
    allocate(Le(Nlo,Nlo));Le=zero
    allocate(Te(Nlo,Nlo));Te=zero

    ! -----------------------------------------------------------
    ! temporary patch:
    ! to be back-compatible with masks from Norb=1 calculations
    ! the following loop does not read the orbital indexes
    if(Norb==1)then
       ! set masks in latiice indexes
       allocate(lmask(Nlo,Nlo),rmask(Nlo,Nlo))
       lmask(:,:)=0
       rmask(:,:)=0
       lfile = get_file_length("lmask.in")
       unit = free_unit()
       open(unit,file='lmask.in',status='old')
       do i=1,lfile
          read(unit,*) ilat,jlat !does not read iorb & jorb
          ilat=ilat+1
          jlat=jlat+1
          lmask(ilat,jlat)=1
          write(6,*) ilat,jlat,lmask(ilat,jlat)
       enddo
       lfile = get_file_length("rmask.in")
       unit = free_unit()
       open(unit,file='rmask.in',status='old')
       do i=1,lfile
          read(unit,*) ilat,jlat
          ilat=ilat+1
          jlat=jlat+1
          rmask(ilat,jlat)=1
          write(6,*) ilat,jlat,rmask(ilat,jlat)
       enddo
    else
       ! set masks in latiice-orbital indexes
       allocate(lmask(Nlo,Nlo),rmask(Nlo,Nlo))
       lmask(:,:)=0
       rmask(:,:)=0
       lfile = get_file_length("lmask.in")
       unit = free_unit()
       open(unit,file='lmask.in',status='old')
       do i=1,lfile
          read(unit,*) ilat,iorb,jlat,jorb
          ilat=ilat+1
          iorb=iorb+1
          jlat=jlat+1
          jorb=jorb+1
          io = iorb +  (ilat-1)*Norb
          jo = jorb +  (jlat-1)*Norb
          lmask(io,jo)=1
          write(6,*) ilat,iorb,jlat,jorb,lmask(io,jo)
       enddo
       lfile = get_file_length("rmask.in")
       unit = free_unit()
       open(unit,file='rmask.in',status='old')
       do i=1,lfile
          read(unit,*) ilat,iorb,jlat,jorb
          ilat=ilat+1
          iorb=iorb+1
          jlat=jlat+1
          jorb=jorb+1
          io = iorb +  (ilat-1)*Norb
          jo = jorb +  (jlat-1)*Norb
          rmask(io,jo)=1
          write(6,*) ilat,iorb,jlat,jorb,rmask(io,jo)
       enddo
    endif

    ! allocate spin-resolved transmission coefficient
    allocate(transe(Nspin,Lreal))

    ! allocate self-energy of the leads for a give frequency 
    allocate(self_lead(Nlso,Nlso))
 
    ! allocate non-local Green's function for a given frequency
    allocate(Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb))


    do i=1,Lreal
       ! evaluate hybridization function for wr(i)
       call get_self_lead(wr(i),self_lead)
       ! evaluate non-local Green's function for wr(i)
       call dmft_gij_realaxis_wr(Hij(:,:,1),Greal,Sreal(:,:,:,:,:,i),wr(i),self_lead)
       !
       do ispin=1,Nspin
          ! fill auxiliary matrix [Nlso][Nlso]
          do ilat=1,Nlat
             do jlat=1,Nlat
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb +  (ilat-1)*Norb
                      jo = jorb +  (jlat-1)*Norb
                      is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb !== ilat
                      js = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb !== jlat
                      !
                      ! auxiliary retarded Green's function
                      GR(io,jo)=Greal(ilat,jlat,ispin,ispin,iorb,jorb)
                      !
                      ! set \Gamma matrix for L/R according to masks to select L-subset OR R-subset
                      ! R-subset
                      HR(io,jo)=zero
                      if(rmask(io,jo)==1) HR(io,jo) = cmplx(2.d0*dimag(self_lead(is,js)),0d0)
                      ! L-subset
                      HL(io,jo)=zero
                      if(lmask(io,jo)==1) HL(io,jo) = cmplx(2.d0*dimag(self_lead(is,js)),0d0)
                   enddo
                enddo
             enddo
          enddo
          ! advanced Green's function
          GA=conjg(transpose(GR))
          !
          ! get transmission function as T(ispin,i)=Tr[Gadvc*Hybl*Gret*Hybr]
          Re = matmul(GR,HR)
          Le = matmul(GA,HL)
          Te = matmul(Le,Re)
          transe(ispin,i) = trace_matrix(Te,Nlo)
       enddo
    enddo
    !
    ! write transport coefficient of disk
    do ispin=1,Nspin
       suffix="_s"//reg(txtfy(ispin))//"_realw.ed"
       call splot("Te"//trim(suffix),wr,transe(ispin,:))
    enddo
    
    deallocate(GR,HR,GA,HL)
    deallocate(rmask,lmask)
    deallocate(Re,Le)

    if(jbias)then
       !
       ! evaluate spin-resolved current as:
       ! J = \int_{-\infty}^{\infty} de T(e)*(f_L(e)-f_R(e))
       allocate(jcurr(Nspin));jcurr=0.d0
       !
       unit_in = free_unit()
       open(unit_in,file='jbias.in',status='old')
       unit_out= free_unit()
       open(unit_out,file="jbias.ed")
       do
          read(unit_in,*,IOSTAT=EOF)lbias,rbias
          if(EOF<0)exit
          !
          ! write L/R bias voltages
          write(unit_out,'(2f16.9)',advance='no')lbias,rbias
          !
          jcurr=0.d0
          do ispin=1,Nspin
              do i=1,Lreal
                 jcurr(ispin) = jcurr(ispin) + transe(ispin,i)* &
                                (fermi(wr(i)-lbias,beta)-fermi(wr(i)-rbias,beta))* &
                                abs(wfin-wini)/Lreal
              enddo
              !
              ! write spin-resolved current on disk
              write(unit_out,'(1f16.9)',advance='no')jcurr(ispin)
          enddo
          write(unit_out,*) ! newline
       enddo
       close(unit_in)
       close(unit_out)
       !
       deallocate(jcurr)
       !
    endif

    deallocate(Te)

  end subroutine ed_transport


  !----------------------------------------------------------------------------------------!
  ! purpose: define the hybridization matrix of size [Nlat][Nlat][Nspin][Norb][Norb][Lreal] 
  ! reading the parameters from an input file
  !----------------------------------------------------------------------------------------!
  subroutine get_self_lead(wr,self_lead)
    real(8),intent(in)                      :: wr
    complex(8),dimension(:,:),intent(inout) :: self_lead(:,:) ![Nlso][Nlso]
    complex(8),dimension(:,:),allocatable   :: lead_real ![Nlead][Nspin]
    integer                                 :: ilat,jlat,ispin,jspin,iorb,jorb,io,jo,i,Nlso
    integer                                 :: unit,l,lfile
    integer                                 :: ikind,ilead,Nlead
    real(8)                                 :: D,mu,V,epsk
    integer                                 :: k,kmax
    complex(8)                              :: ksum
    character(50)                           :: suffix
    !
    Nlso = size(self_lead,1)
    !
    ! initialize embedding hybridization function
    self_lead=zero

    ! determine Nleads & allocate lead matrix
    lfile = get_file_length("lead.in")
    unit = free_unit()
    open(unit,file='lead.in',status='old')
    read(unit,*)Nlead
    !
    allocate(lead_real(Nlead,Nspin))
    lead_real(:,:)=zero
     
    ! lead file setup lead by kind, half-bandwitdh (D) and chemical potential (mu)
    ! *** note: reading the file for each wr is slow and inefficient 
    do l=1,lfile-1 ! because Nlead was read separately above
       read(unit,*) ilead, ispin, D, mu, ikind
       ilead=ilead+1
       ispin=ispin+1
       if(ilead>Nlead)stop "set_hyb error: in input file 'lead.in' ilead > Nlead"
       if(ispin>Nspin)stop "set_hyb error: in input file 'lead.in' ispin > Nspin"
       !
       ! set the lead's Green's function, depending on ikind
       if(ikind==0)then
          ! flat DOS (analytic)
          !write(*,*) "flat DOS (analytic)"
          lead_real(ilead,ispin)=dcmplx( log(abs((D+wr+mu)/(D-wr-mu))) , -pi*heaviside(D-abs(wr+mu)) )/(2d0*D)
       elseif(ikind==1)then
          ! flat DOS (k-sum)
          !write(*,*) "flat DOS (k-sum)"
          kmax=10000
          ksum=zero
          do k=1,kmax
             epsk = -D + 2*D/kmax*(k-1)
             ksum = ksum + 1d0/( wr+xi*eps+mu - epsk)
          enddo
          lead_real(ilead,ispin)=ksum/kmax
       elseif(ikind==2)then
          ! broad-band limit
          !write(*,*) "broad-band limit (analytic)" 
          lead_real(ilead,ispin)=dcmplx(0d0,-1.d0*pi) ! to ensure DOS normalization
       elseif(ikind==3)then
          ! semicircular DOS (k-sum) 
          !write(*,*) "semicircular DOS (k-sum)"
          ksum=zero
          do k=1,kmax
             epsk = -D + 2*D/kmax*(k-1)
             ksum = ksum + (4d0/(pi*kmax))*sqrt(1d0-(epsk/D)**2)/( wr+xi*eps+mu - epsk)
          enddo
          lead_real(ilead,ispin)=ksum
       elseif(ikind==4)then
          ! readin hk DOS
          write(*,*) "readin hk DOS to be implemented and benchmarked w/ w2dynamics"
          stop
       else
          write(*,*) "set_hyb error: in input file 'lead.in' invalid ikind"
          stop
       endif
       !*** broken if wr is not an array: fix it!
       !! store lead(s) DOS on disk
       !!suffix="_ilead"//reg(txtfy(ilead))//"_s"//reg(txtfy(ispin))//"_realw.ed"
       !!call splot("lead"//trim(suffix),wr,lead_real(ilead,ispin,:))
    enddo
    close(unit)
    
    ! hybridization file determine lead-site connections 
    ! *** note: reading the file for each wr is slow and inefficient 
    lfile = get_file_length("vij.in")
    unit = free_unit()
    open(unit,file='vij.in',status='old')
    do i=1,lfile
       read(unit,*) ilat, iorb, jlat, jorb, ilead, V
       ilat=ilat+1
       iorb=iorb+1
       jlat=jlat+1
       jorb=jorb+1
       ilead=ilead+1
       if((iorb>Norb).or.(jorb>Norb))stop "set_hyb error: in input file 'vij.in' i/jorb > Norb"
       if((ilat>Nlat).or.(jlat>Nlat))stop "set_hyb error: in input file 'vij.in' i/jlat > Nlat"
       if(ilead>Nlead)stop "set_hyb error: in input file 'vij.in' ilead > Nlead"
       do ispin=1,Nspin
          ! get stride and set matrix element: no symmetrization
          io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb !== ilat
          jo = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb !== jlat
          self_lead(io,jo) = self_lead(io,jo) + lead_real(ilead,ispin)*V**2
          !*** broken if wr is not an array: fix it!
          !! store self-energy of the lead(s) on disk
          !!suffix="_i"//reg(txtfy(ilat))//"_j"//reg(txtfy(jlat))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          !!call splot("Hyb"//trim(suffix),wr,Hyb_real(io,jo,:))
       enddo
    enddo
    close(unit)
    deallocate(lead_real)
  end subroutine get_self_lead


  !----------------------------------------------------------------------------------------!
  ! purpose: evaluate the non-local Green's function for a given frequency wr
  !----------------------------------------------------------------------------------------!
  subroutine dmft_gij_realaxis_wr(Hij,Greal,Sreal,wr,self_lead)
    complex(8),dimension(:,:),intent(in)            :: Hij       ![Nlso][Nlso]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Sreal     !      [Nlat][Nspin][Nspin][Norb][Norb]
    complex(8),dimension(:,:),intent(in),optional   :: self_lead ![Nlso][Nlso]
    complex(8),dimension(:,:,:),allocatable         :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb]
    real(8),intent(in)                              :: wr
    real(8)                                         :: xmu,eps
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso
    !
    !Retrieve parameters:
    call get_ctrl_var(xmu,"XMU")
    call get_ctrl_var(eps,"EPS")
    !
    Nlat  = size(Sreal,1)
    Nspin = size(Sreal,2)
    Norb  = size(Sreal,4)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    call assert_shape(Hij,[Nlso,Nlso],"dmft_gij_realaxis_wr","Hij")
    call assert_shape(Sreal,[Nlat,Nspin,Nspin,Norb,Norb],"dmft_gij_realaxis_wr","Sreal")
    call assert_shape(Greal,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"dmft_gij_realaxis_wr","Greal")
    if(present(self_lead)) call assert_shape(self_lead,[Nlso,Nlso],"dmft_gij_realaxis_wr","self_lead")
    !
    allocate(zeta_real(Nlat,Nso,Nso));zeta_real=zero
    !
    do ilat=1,Nlat
       zeta_real(ilat,:,:) = (wr+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:),NSpin,Norb)
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    !call start_timer
    Greal=zero
    if(present(self_lead))then
       call invert_gij(zeta_real,Hij,Greal,self_lead)
    else
       call invert_gij(zeta_real,Hij,Greal)
    endif
    !call stop_timer
    !call dmft_print_gf_realaxis(Greal,"Gij",iprint=1)
  end subroutine dmft_gij_realaxis_wr


  !----------------------------------------------------------------------------------------!
  ! purpose: embed self_lead into Gij
  !----------------------------------------------------------------------------------------!
  subroutine invert_gij(zeta,Hk,Gij,self_lead)
    complex(8),dimension(:,:,:),intent(in)          :: zeta      ![Nlat][Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:),intent(in)            :: Hk        ![Nlso][Nlso]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gij       ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    complex(8),dimension(:,:),intent(in),optional   :: self_lead ![Nlso][Nlso]
    !allocatable arrays
    complex(8),dimension(:,:),allocatable           :: Gmatrix   ![Nlso][Nlso]
    integer                                         :: Nlat,Nspin,Norb,Nso,Nlso
    integer                                         :: ialt,jlat,iorb,jorb,ispin,jspin,io,jo
    !
    Nlat  = size(zeta,1)
    Nspin = size(Gij,3)
    Norb  = size(Gij,5)
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    call assert_shape(zeta,[Nlat,Nso,Nso],"invert_gij","zeta")
    call assert_shape(Hk,[Nlso,Nlso],"invert_gij","Hk")
    call assert_shape(Gij,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"invert_gij","Gij")
    if(present(self_lead)) call assert_shape(self_lead,[Nlso,Nlso],"invert_gij","self_lead")
    !
    allocate(Gmatrix(Nlso,Nlso))
    Gij=zero
    Gmatrix  = blocks_to_matrix(zeta(:,:,:),Nlat,Nso) - Hk
    if(present(self_lead)) Gmatrix = Gmatrix - self_lead
    call inv(Gmatrix) 
    !store the diagonal blocks directly into the tmp output 
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                      jo = jorb + (jspin-1)*Norb + (jlat-1)*Nspin*Norb
                      Gij(ilat,jlat,ispin,jspin,iorb,jorb) = Gmatrix(io,jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine invert_gij


  !----------------------------------------------------------------------------------------!
  ! purpose: read the real local self-energy from disk
  !----------------------------------------------------------------------------------------!
  subroutine read_sigma_real(Sreal)
    complex(8),intent(inout)         :: Sreal(:,:,:,:,:,:)
    character(len=30)                :: suffix
    integer                          :: ilat,ispin,iorb
    real(8),dimension(:),allocatable :: wm,wr
    !
    if(size(Sreal,2)/=Nspin) stop "save_sigma: error in dim 2. Nspin"
    if(size(Sreal,3)/=Nspin) stop "save_sigma: error in dim 3. Nspin"
    if(size(Sreal,4)/=Norb) stop "save_sigma: error in dim 4. Norb"
    if(size(Sreal,5)/=Norb) stop "save_sigma: error in dim 5. Norb"
    !
    allocate(wr(Lreal))
    !
    wr = linspace(wini,wfin,Lreal)
    write(LOGfile,*)"write spin-orbital diagonal elements:"
    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          call sread("LSigma"//trim(suffix),wr,Sreal(:,ispin,ispin,iorb,iorb,:))
       enddo
    enddo
  end subroutine read_sigma_real


  function blocks_to_matrix(Vblocks,Nlat,Nso) result(Matrix)
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: i,j,ip
    Matrix=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Matrix(i:j,i:j) =  Vblocks(ip,:,:)
    enddo
  end function blocks_to_matrix


  function extract_Hloc(Hk,Nlat,Nspin,Norb) result(Hloc)
    complex(8),dimension(:,:,:)                 :: Hk
    integer                                     :: Nlat,Nspin,Norb
    complex(8),dimension(size(Hk,1),size(Hk,2)) :: Hloc
    !
    integer                                     :: iorb,ispin,ilat,is
    integer                                     :: jorb,jspin,js
    Hloc = zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hloc(is,js) = sum(Hk(is,js,:))/size(Hk,3)
                enddo
             enddo
          enddo
       enddo
    enddo
    where(abs(dreal(Hloc))<1.d-9)Hloc=0d0
  end function extract_Hloc


  function trace_matrix(M,dim) result(tr)
    integer                       :: dim
    complex(8),dimension(dim,dim) :: M
    complex(8) :: tr
    integer                       :: i
    tr=dcmplx(0d0,0d0)
    do i=1,dim
       tr=tr+M(i,i)
    enddo
  end function trace_matrix


  function get_file_length(file) result(lines)
    integer           :: lines
    character(len=*)  :: file
    integer           :: ifile,ierr,pos
    logical           :: IOfile,bool,bool1,bool2
    character(len=256)::buffer
    inquire(file=reg(file),exist=IOfile)
    if(.not.IOfile)then
       inquire(file=reg(file)//".gz",exist=IOfile)
       if(IOfile)call file_gunzip(reg(file))
    endif
    lines=0
    if(.not.IOfile)then
       write(*,*) 'Cannot read +'//reg(file)//'. Skip file_size'
       return
    endif
    open(99,file=reg(file))
    ierr=0
    do while(ierr==0)
       lines=lines+1
       read(99,*,iostat=ierr)buffer
       bool1=scan(buffer,"#").ne.0
       bool2=len_trim(buffer).eq.0       
       if(bool1 .OR. bool2)lines=lines-1
    enddo
    lines=lines-1
    !write(*,'(A,I9,A)') 'there are', lines,' lines in +'//reg(file)
    rewind(99)
    close(99)
  end function get_file_length



end program trans
