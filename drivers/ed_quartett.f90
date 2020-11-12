program ed_quartett
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none


  !#########   VARIABLES DECLARATION   #########
  character(len=60)                              :: finput
  integer(8),dimension(:,:),allocatable          :: Xstride
  integer(8),dimension(:,:,:),allocatable        :: Vstride
  integer(8),dimension(:),allocatable            :: Neigh
  integer(8),dimension(:,:),allocatable          :: lat2vec,vec2lat
  real(8),dimension(:,:),allocatable             :: Radius
  real(8),dimension(:),allocatable               :: distprint
  real(8),dimension(:),allocatable               :: Vnn_used
  real(8)                                        :: bottom,Xc,Yc
  !Mpi:
  integer                                        :: comm,rank,ier
  logical                                        :: master
  !Stuff that I need by construction:
  integer                                        :: Nb
  real(8),dimension(:),allocatable               :: Bath
  complex(8),dimension(:,:,:,:),allocatable      :: Hloc_nn
  !auxiliary lattices:
  integer,dimension(:,:),allocatable             :: lattice_irred
  integer,dimension(:,:),allocatable             :: lattice_tiled
  !
  integer                                        :: Maxdist,tiling
  integer                                        :: tiling_x,tiling_y,itiling
  integer                                        :: ilatt,jlatt,ineig,distance
  integer                                        :: row,col
  logical                                        :: stripe


  !#########   MPI INITIALIZATION   #########
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  !rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  !#########    VARIABLE PARSING    #########
  call parse_cmd_variable(finput,"FINPUT",default='inputED_quartet.in')
  call ed_read_input(trim(finput),comm)
  call parse_input_variable(tiling,"TILING",finput,default=5)
  call parse_input_variable(stripe,"STRIPE",finput,default=.false.)
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nbath,"NBATH")
  call add_ctrl_var(Vnn,"VNN")
  call add_ctrl_var(HardCoreBoson,"HARDCOREBOSON")
  call add_ctrl_var(filling,"FILLING")
  call add_ctrl_var(Ust,"UST")


  !######### VECTOR-LATTICE WRAPPIN #########
  tiling_x=tiling
  tiling_y=tiling
  Maxdist=21 !tiling_x*Nbath-2
  if(stripe)then
     Maxdist=Norb*Nbath-2
     tiling_y=1
  elseif(Nbath.eq.1)then
     Maxdist=Norb-1
     tiling_y=1
  endif
  call create_lattice(Maxdist)


  !######### INTERACION REMODULATON #########
  if(HardCoreBoson.ne.0)then
     if(Nbath.ne.1)then
        allocate(Vnn_used(8));Vnn_used=0d0
        Vnn_used(1) = 0d0
        Vnn_used(2) = 0d0
        Vnn_used(3) = 2d0*(Vnn(1)+Vnn(2))+4d0*(Vnn(3)+Vnn(4))
        Vnn_used(4) = Vnn(1)+2d0*Vnn(2)+2d0*Vnn(3)+4d0*Vnn(4)
        Vnn_used(5) = Vnn(2)+4d0*Vnn(4)
        Vnn_used(6) = 2d0*Vnn(3)
        Vnn_used(7) = Vnn(3)+2d0*Vnn(4)
        Vnn_used(8) = Vnn(4)
     elseif(Nbath.eq.1)then
        allocate(Vnn_used(3));Vnn_used=0d0
        Vnn_used(1) = 0d0
        Vnn_used(2) = 2d0*(Vnn(1)+Vnn(2))+4d0*(Vnn(3)+Vnn(4))
        Vnn_used(3) = 2d0*Vnn(3)
     endif
  else
     allocate(Vnn_used(size(Vnn)));Vnn_used=0d0
     Vnn_used=Vnn
  endif
  !
  call ed_set_Vstride(Comm,0.5d0*Vnn_used,Xstride,Vstride,Neigh,lat2vec,vec2lat,Radius,distprint)
  call MPI_Barrier(Comm,ier)


  !############ PRINT SOME INFOs ############
  if(master)then
     !
     write(LOGfile,*)
     write(LOGfile,*)
     !
     do row=1,Nbath*3!tiling_y
       write(LOGfile,"(1000I4)") (lattice_tiled(row,col),col=1,Norb*3)!tiling_x)
     enddo
     write(LOGfile,*)
     write(LOGfile,*)
     !
     do row=1,Nbath
       write(LOGfile,"(1000F6.2)") (Radius(row,col),col=1,Norb)
     enddo
     write(LOGfile,*)
     write(LOGfile,*)
     !
     do ilatt=1,Norb*Nbath
       write(LOGfile,'(1I5,2A,100I3)')ilatt," Xstr "," - ",(Xstride(ilatt,ineig),ineig=1,4)
     enddo
     write(LOGfile,*)
     write(LOGfile,*)
     !
     write(LOGfile,'(A)')"Stride referred to the upper right corner:"
     do distance=1,size(Vnn_used)
       write(LOGfile,'(1I5,2F12.2,A,100I3)')distance,distprint(distance),Vnn_used(distance)," -",(Vstride(Norb,distance,ineig),ineig=1,Neigh(distance))
     enddo
     do distance=size(Vnn_used)+1,Maxdist
       write(LOGfile,'(1I5,2F12.2,A,100I3)')distance,distprint(distance),0d0," -",(Vstride(Norb,distance,ineig),ineig=1,Neigh(distance))
     enddo
     write(LOGfile,*)
     write(LOGfile,*)
     write(LOGfile,*)
     !
  endif


  !######### SOLVER INITIALIZATION  #########
  bath_type = "normal"
  ed_mode = "nonsu2"
  Nb=get_bath_dimension()
  allocate(Bath(Nb));Bath=0.0d0
  allocate(Hloc_nn(Nspin,Nspin,Norb,Norb)); Hloc_nn=zero
  !
  call ed_init_solver(Comm,Bath,Hloc_nn)


  !#########    QUARTETT SOLUTION   #########
  if (master) call start_loop(1,1,"DMFT-loop")
  call ed_solve(comm,Bath)
  !
  if(master)call end_loop



contains



   subroutine create_lattice(MaxNeig)
      implicit none
      integer,intent(in)          :: MaxNeig
      !
      logical                     :: insert
      integer                     :: Rx,Ry,idist,jdist
      integer                     :: row_i,col_i,row_j,col_j
      real(8)                     :: dist
      real(8),allocatable         :: distances(:),distances_tmp(:)
      integer,allocatable         :: order(:),doneNeigh(:)
      real                        :: shift
      !
      !
      allocate(Neigh(MaxNeig));Neigh=0
      allocate(Xstride(Nbath*Norb,4)); Xstride=0
      allocate(Vstride(Nbath*Norb,MaxNeig,50)); Vstride=0
      allocate(lattice_irred(Nbath,Norb)); lattice_irred=0
      allocate(lattice_tiled(tiling_y*Nbath,tiling_x*Norb)); lattice_tiled=0
      allocate(lat2vec(Nbath,Norb));lat2vec=0
      allocate(vec2lat(Nbath*Norb,2));vec2lat=0
      allocate(Radius(Nbath,Norb));Radius=0d0
      !
      !
      shift=1d0
      !
      !
      Xc = dble(Norb)/2  + 0.5*shift
      Yc = dble(Nbath)/2 + 0.5*shift
      ilatt=0
      do row=1,Nbath
        do col=1,Norb
            ilatt=ilatt+1
            lattice_irred(row,col) = ilatt
            lat2vec(row,col) = ilatt
            vec2lat(ilatt,1) = row
            vec2lat(ilatt,2) = col
            Radius(row,col) = sqrt( (Yc - dble(row))**2 + (Xc - dble(col))**2 )
         enddo
      enddo
      bottom=minval(Radius)
      Radius=Radius-bottom
      !
      !
      if((Norb.eq.8).and.(Nbath.eq.8).and.(Ust.ne.0d0))then
         Radius=0d0
         if(HardCoreBoson.ne.0)then
            !
            Radius(3,3)=-1d0
            !
         else
            !
            Radius(4:5,4:5)=-1d0
            !
         endif
      endif
      !
      !
      do ilatt=1,tiling_y
         do jlatt=1,tiling_x
            lattice_tiled(1+(ilatt-1)*Nbath:Nbath+(ilatt-1)*Nbath,1+(jlatt-1)*Norb:Norb+(jlatt-1)*Norb) = lattice_irred
         enddo
      enddo
      !
      !
      allocate(distances(MaxNeig**2),order(MaxNeig**2))
      distances=100000.d0
      idist=0
      do Ry=0,MaxNeig
         if((Nbath.eq.1).and.(Ry.ne.0))cycle
         do Rx=Ry,MaxNeig
            if(Rx.eq.0 .and. Ry.eq.0)cycle
            if(stripe  .and. Ry.gt.1)cycle
            dist = sqrt( dble(Rx)**2 + dble(Ry)**2 )
            !
            insert = .true.
            do jdist=1,size(distances)
               insert = insert .and. (abs(distances(jdist)-dist).gt.1e-6)
            enddo
            !
            if(insert)then
               idist=idist+1
               distances(idist) = dist
               !write(*,*)"insert ",dist," at ",idist
            endif
            !
         enddo
      enddo
      call sort_quicksort(distances,order)
      !
      allocate(distprint(Maxdist))
      distprint=distances(1:MaxNeig)
      !write(*,*)
      !write(*,*)
      !do idist=1,MaxNeig
      !   write(*,*)idist,distances(idist)
      !enddo
      !
      ! find the number of Neighbors
      row_i = vec2lat(Norb,1) + floor(tiling_y/2.)*Nbath
      col_i = vec2lat(Norb,2) + floor(tiling_x/2.)*Norb
      do row_j=1,Nbath*tiling_y
         do col_j=1,Norb*tiling_x
            !
            dist = sqrt( dble(col_i-col_j)**2 + dble(row_i-row_j)**2 )
            !
            do idist=1,MaxNeig
               if(abs(distances(idist)-dist).lt.1e-6)then
                  Neigh(idist)=Neigh(idist)+1
               endif
            enddo
            !
         enddo
      enddo
      !
      ! find the coordinates of each couple of sites
      allocate(doneNeigh(MaxNeig))
      do ilatt=1,Nbath*Norb
         !
         row_i = vec2lat(ilatt,1) + floor(tiling_y/2.)*Nbath
         col_i = vec2lat(ilatt,2) + floor(tiling_x/2.)*Norb
         !
         doneNeigh=0
         do row_j=1,Nbath*tiling_y
            do col_j=1,Norb*tiling_x
               !
               dist = sqrt( dble(col_i-col_j)**2 + dble(row_i-row_j)**2 )
               !
               do idist=1,MaxNeig
                  if(abs(distances(idist)-dist).lt.1e-6)then
                     doneNeigh(idist)=doneNeigh(idist)+1
                     Vstride(ilatt,idist,doneNeigh(idist)) = lattice_tiled(row_j,col_j)
                  endif
               enddo
               !
            enddo
         enddo
      enddo
      deallocate(doneNeigh,order,distances)
      !
      !
      if(stripe)then
         do ilatt=1,Nbath*Norb
            row_i = vec2lat(ilatt,1) + floor(tiling_y/2.)*Nbath
            col_i = vec2lat(ilatt,2) + floor(tiling_x/2.)*Norb
            if(row_i.eq.1)then
               Xstride(ilatt,1) = lattice_tiled(row_i,col_i)
               Xstride(ilatt,2) = lattice_tiled(row_i,col_i+1)
               Xstride(ilatt,3) = lattice_tiled(row_i+1,col_i+1)
               Xstride(ilatt,4) = lattice_tiled(row_i+1,col_i)
            elseif(row_i.eq.2)then
               Xstride(ilatt,1) = lattice_tiled(row_i,col_i)
               Xstride(ilatt,2) = lattice_tiled(row_i,col_i+1)
               Xstride(ilatt,3) = lattice_tiled(row_i-1,col_i+1)
               Xstride(ilatt,4) = lattice_tiled(row_i-1,col_i)
            endif
         enddo
      elseif(Nbath.eq.1)then
         do ilatt=1,Nbath*Norb
            row_i = vec2lat(ilatt,1) + floor(tiling_y/2.)*Nbath
            col_i = vec2lat(ilatt,2) + floor(tiling_x/2.)*Norb
            Xstride(ilatt,1) = lattice_tiled(row_i,col_i)
            Xstride(ilatt,2) = lattice_tiled(row_i,col_i+1)
            Xstride(ilatt,3) = lattice_tiled(row_i,col_i-1)
         enddo
      else
         do ilatt=1,Nbath*Norb
            row_i = vec2lat(ilatt,1) + floor(tiling_y/2.)*Nbath
            col_i = vec2lat(ilatt,2) + floor(tiling_x/2.)*Norb
            Xstride(ilatt,1) = lattice_tiled(row_i,col_i)
            Xstride(ilatt,2) = lattice_tiled(row_i,col_i+1)
            Xstride(ilatt,3) = lattice_tiled(row_i+1,col_i+1)
            Xstride(ilatt,4) = lattice_tiled(row_i+1,col_i)
         enddo
      endif
      !
      !
   end subroutine create_lattice


end program ed_quartett



!subroutine create_lattice_old(MaxNeig)
!   implicit none
!   integer,intent(in) :: MaxNeig
!   !
!   allocate(Neigh(MaxNeig));Neigh=zero
!   Neigh(1:3)=4
!   Neigh(4)=8
!   Neigh(5:6)=4
!   !
!   !
!   allocate(Vstride(Nbath*Norb,MaxNeig,8)); Vstride=zero
!   allocate(lattice_irred(Nbath,Norb)); lattice_irred=zero
!   allocate(lattice_tiled(tiling*Nbath,tiling*Norb)); lattice_tiled=zero
!   allocate(lat2vec(Nbath,Norb));lat2vec=0
!   allocate(vec2lat(Nbath*Norb,2));vec2lat=0
!   allocate(Radius(Nbath,Norb));Radius=0d0
!   !
!   !
!   Xc = dble(Norb)/2  + 0.5 ! 3.5 ! dble(Norb)/2  + 0.5
!   Yc = dble(Nbath)/2 + 0.5 ! 4.5 ! dble(Nbath)/2 + 0.5
!   ilatt=0
!   do row=1,Nbath
!     do col=1,Norb
!         ilatt=ilatt+1
!         lattice_irred(row,col) = ilatt
!         lat2vec(row,col) = ilatt
!         vec2lat(ilatt,1) = row
!         vec2lat(ilatt,2) = col
!         Radius(row,col) = sqrt( (Yc - dble(row))**2 + (Xc - dble(col))**2 )
!      enddo
!   enddo
!   bottom=minval(Radius)
!   Radius=Radius-bottom
!   do ilatt=1,tiling
!      do jlatt=1,tiling
!         lattice_tiled(1+(ilatt-1)*Nbath:Nbath+(ilatt-1)*Nbath,1+(jlatt-1)*Norb:Norb+(jlatt-1)*Norb) = lattice_irred
!      enddo
!   enddo
!   !
!   !
!   do ilatt=1,Nbath*Norb
!      !
!      row = vec2lat(ilatt,1) + floor(tiling/2.)*Nbath
!      col = vec2lat(ilatt,2) + floor(tiling/2.)*Norb
!      !
!      !nn-1
!      Vstride(ilatt,1,1) = lattice_tiled(row-1,col)
!      Vstride(ilatt,1,2) = lattice_tiled(row,col+1)
!      Vstride(ilatt,1,3) = lattice_tiled(row+1,col)
!      Vstride(ilatt,1,4) = lattice_tiled(row,col-1)
!      !
!      !nn-2
!      Vstride(ilatt,2,1) = lattice_tiled(row-1,col+1)
!      Vstride(ilatt,2,2) = lattice_tiled(row+1,col+1)
!      Vstride(ilatt,2,3) = lattice_tiled(row+1,col-1)
!      Vstride(ilatt,2,4) = lattice_tiled(row-1,col-1)
!      !
!      !nn-3
!      Vstride(ilatt,3,1) = lattice_tiled(row-2,col)
!      Vstride(ilatt,3,2) = lattice_tiled(row,col+2)
!      Vstride(ilatt,3,3) = lattice_tiled(row+2,col)
!      Vstride(ilatt,3,4) = lattice_tiled(row,col-2)
!      !
!      !nn-4
!      Vstride(ilatt,4,1) = lattice_tiled(row-2,col+1)
!      Vstride(ilatt,4,2) = lattice_tiled(row-1,col+2)
!      Vstride(ilatt,4,3) = lattice_tiled(row+1,col+2)
!      Vstride(ilatt,4,4) = lattice_tiled(row+2,col+1)
!      Vstride(ilatt,4,5) = lattice_tiled(row+2,col-1)
!      Vstride(ilatt,4,6) = lattice_tiled(row+1,col-2)
!      Vstride(ilatt,4,7) = lattice_tiled(row-1,col-2)
!      Vstride(ilatt,4,8) = lattice_tiled(row-2,col-1)
!      !
!      !nn-5
!      Vstride(ilatt,5,1) = lattice_tiled(row-2,col+2)
!      Vstride(ilatt,5,2) = lattice_tiled(row+2,col+2)
!      Vstride(ilatt,5,3) = lattice_tiled(row+2,col-2)
!      Vstride(ilatt,5,4) = lattice_tiled(row-2,col-2)
!      !
!      !nn-6
!      Vstride(ilatt,6,1) = lattice_tiled(row-3,col)
!      Vstride(ilatt,6,2) = lattice_tiled(row,col+3)
!      Vstride(ilatt,6,3) = lattice_tiled(row+3,col)
!      Vstride(ilatt,6,4) = lattice_tiled(row,col-3)
!      !
!   enddo
!   !
!end subroutine create_lattice_old
