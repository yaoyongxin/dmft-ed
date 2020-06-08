program ed_quartett
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none


  !#########   VARIABLES DECLARATION   #########
  character(len=60)                              :: finput
  integer,dimension(:,:,:),allocatable           :: Vstride
  integer,dimension(:),allocatable               :: Neigh
  !Mpi:
  integer                                        :: comm,rank,ier
  logical                                        :: master
  !Stuff that I need by construction:
  integer                                        :: Nb
  real(8),dimension(:),allocatable               :: Bath
  complex(8),dimension(:,:,:,:),allocatable      :: Hloc_nn



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


  !#########   FIELDS ALLOCATION    #########
  allocate(Vstride(Nbath*Norb,4,8)); Vstride=zero
  allocate(Neigh(4)); Neigh=zero


  !######### SOLVER INITIALIZATION  #########
  bath_type = "normal"
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
