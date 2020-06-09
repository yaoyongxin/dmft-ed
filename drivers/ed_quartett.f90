program ed_quartett
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none


  !#########   VARIABLES DECLARATION   #########
  character(len=60)                              :: finput
  integer(8),dimension(:,:,:),allocatable        :: Vstride
  integer(8),dimension(:),allocatable            :: Neigh
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
  integer,parameter                              :: tiling = 5
  integer                                        :: ilatt,jlatt
  integer                                        :: row,col


  !#########   MPI INITIALIZATION   #########
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  !rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  !#########    VARIABLE PARSING    #########
  call parse_cmd_variable(finput,      "FINPUT",  default='inputED_S.in')
  call ed_read_input(trim(finput),comm)
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,     "NORB")
  call add_ctrl_var(Nbath,    "NBATH")


  !######### VECTOR-LATTICE WRAPPIN #########
  allocate(Neigh(6)); Neigh=zero !Max is 7 as the max Vnn
  Neigh(1:3)=4
  Neigh(4)=8
  Neigh(5:6)=4
  !
  !
  allocate(Vstride(Nbath*Norb,6,8)); Vstride=zero
  allocate(lattice_irred(Nbath,Norb)); lattice_irred=zero
  allocate(lattice_tiled(tiling*Nbath,tiling*Norb)); lattice_tiled=zero
  !
  !
  ilatt=0
  do col=1,Norb
     do row=1,Nbath
        ilatt=ilatt+1
        lattice_irred(row,col) = ilatt
     enddo
  enddo
  do ilatt=1,tiling
     do jlatt=1,tiling
        lattice_tiled(1+(ilatt-1)*Nbath:Nbath+(ilatt-1)*Nbath,1+(jlatt-1)*Norb:Norb+(jlatt-1)*Norb) = lattice_irred
     enddo
  enddo
  !
  !
  do ilatt=1,Nbath*Norb
     !
     row = floor(tiling/2.)*Nbath + floor(dble(ilatt)/Norb)+1
     col = floor(tiling/2.)*Norb  + mod(ilatt,Norb)
     !
     !nn-1
     Vstride(ilatt,1,1) = lattice_tiled(row-1,col)
     Vstride(ilatt,1,2) = lattice_tiled(row,col+1)
     Vstride(ilatt,1,3) = lattice_tiled(row+1,col)
     Vstride(ilatt,1,4) = lattice_tiled(row,col-1)
     !
     !nn-2
     Vstride(ilatt,2,1) = lattice_tiled(row-1,col+1)
     Vstride(ilatt,2,2) = lattice_tiled(row+1,col+1)
     Vstride(ilatt,2,3) = lattice_tiled(row+1,col-1)
     Vstride(ilatt,2,4) = lattice_tiled(row-1,col-1)
     !
     !nn-3
     Vstride(ilatt,3,1) = lattice_tiled(row-2,col)
     Vstride(ilatt,3,2) = lattice_tiled(row,col+2)
     Vstride(ilatt,3,3) = lattice_tiled(row+2,col)
     Vstride(ilatt,3,4) = lattice_tiled(row,col-2)
     !
     !nn-4
     Vstride(ilatt,4,1) = lattice_tiled(row-2,col+1)
     Vstride(ilatt,4,2) = lattice_tiled(row-1,col+2)
     Vstride(ilatt,4,3) = lattice_tiled(row+1,col+2)
     Vstride(ilatt,4,4) = lattice_tiled(row+2,col+1)
     Vstride(ilatt,4,5) = lattice_tiled(row+2,col-1)
     Vstride(ilatt,4,6) = lattice_tiled(row+1,col-2)
     Vstride(ilatt,4,7) = lattice_tiled(row-1,col-2)
     Vstride(ilatt,4,8) = lattice_tiled(row-2,col-1)
     !
     !nn-5
     Vstride(ilatt,5,1) = lattice_tiled(row-2,col+2)
     Vstride(ilatt,5,2) = lattice_tiled(row+2,col+2)
     Vstride(ilatt,5,3) = lattice_tiled(row+2,col-2)
     Vstride(ilatt,5,4) = lattice_tiled(row-2,col-2)
     !
     !nn-1
     Vstride(ilatt,6,1) = lattice_tiled(row-3,col)
     Vstride(ilatt,6,2) = lattice_tiled(row,col+3)
     Vstride(ilatt,6,3) = lattice_tiled(row+3,col)
     Vstride(ilatt,6,4) = lattice_tiled(row,col-3)
     !
  enddo

  call ed_set_Vstride(Comm,Vstride,Neigh)
  call MPI_Barrier(Comm,ier)


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


end program ed_quartett
