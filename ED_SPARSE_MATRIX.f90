MODULE ED_SPARSE_MATRIX  !THIS VERSION CONTAINS ONLY COMPLX ELEMENT: (HERMITIAN MATRIX) 
  USE SF_CONSTANTS, only:zero
  USE SF_IOTOOLS, only: str,free_unit
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private

  type sparse_element
     complex(8)                            :: cval !value of the entry: double complex
     integer                               :: col  !col connected to this compress value
     type(sparse_element),pointer          :: next !link to next entry in the row
  end type sparse_element
  public :: sparse_element

  type sparse_row
     integer                               :: size=0    !size of the list
     integer                               :: max_column=0
     type(sparse_element),pointer          :: root    !head/root of the list\== list itself
  end type sparse_row

  type sparse_matrix
     integer                               :: Nrow
     integer                               :: Ncol
     logical                               :: status=.false.
     logical                               :: dryrun=.false.
     type(sparse_row),dimension(:),pointer :: row
#ifdef _MPI
     integer                               :: istart=0 !global start index for MPI storage
     integer                               :: iend=0
     integer                               :: ishift=0
     type(sparse_row),dimension(:),pointer :: loc
#endif
  end type sparse_matrix






  !INIT SPARSE MATRICES (LL)
  interface sp_init_matrix
     module procedure :: sp_init_matrix_ll
#ifdef _MPI
     module procedure :: mpi_sp_init_matrix_ll
#endif
  end interface sp_init_matrix


  !DELETE SPARSE MATRIX (LL) OR ONE OF ITS ELEMENTS (LL)
  interface sp_delete_matrix
     module procedure :: sp_delete_matrix_ll
#ifdef _MPI
     module procedure :: mpi_sp_delete_matrix_ll
#endif
  end interface sp_delete_matrix


  !INSERT ELEMENTS (D,C) IN LL-SPARSE MATRIX
  interface sp_insert_element
     module procedure :: sp_insert_element_c
  end interface sp_insert_element


  !BUILD THE COLUMNS RANGE MIN-MAX
  interface sp_columns_range_matrix
     module procedure :: sp_columns_range_matrix_ll
#ifdef _MPI
     module procedure :: mpi_sp_columns_range_matrix
#endif
  end interface sp_columns_range_matrix



  !LOAD STANDARD MATRIX INTO SPARSE MATRICES
  interface sp_load_matrix
     module procedure :: sp_load_matrix_c
#ifdef _MPI
     module procedure :: mpi_sp_load_matrix_c
#endif
  end interface sp_load_matrix


  !DUMP SPARSE MATRIX INTO STANDARD MATRIX
  interface sp_dump_matrix
     module procedure :: sp_dump_matrix_c
#ifdef _MPI
     module procedure :: mpi_sp_dump_matrix_c
#endif
  end interface sp_dump_matrix


  !PRETTY PRINTING
  interface sp_print_matrix
     module procedure sp_print_matrix_ll
  end interface sp_print_matrix


  !SPY PRINT SPARSE MATRIX
  interface sp_spy_matrix
     module procedure :: sp_spy_matrix_ll
#ifdef _MPI
     module procedure :: mpi_sp_spy_matrix_ll
#endif
  end interface sp_spy_matrix




  public :: sparse_matrix
  !
  public :: sp_init_matrix      !init the sparse matrix   !checked
  public :: sp_delete_matrix    !delete the sparse matrix !checked
  public :: sp_insert_element   !insert an element        !checked
  public :: sp_load_matrix      !create sparse from array !checked
  public :: sp_dump_matrix      !dump sparse into array   !checked
  public :: sp_print_matrix     !print sparse             !checked
  public :: sp_spy_matrix
  public :: sp_size_matrix
  public :: sp_columns_list_matrix
  public :: sp_columns_range_matrix


  integer :: MpiRank=0
  integer :: MpiSize=1
  integer :: MpiIerr
  logical :: MpiMaster=.true.
  integer :: MpiQ
  integer :: MpiR
  integer :: MpiChunk
  integer :: MpiIstart,MpiIend,MpiIshift




contains       




  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  subroutine sp_init_matrix_ll(sparse,N,N1,dryrun)
    type(sparse_matrix),intent(inout) :: sparse
    integer                           :: N
    integer,optional                  :: N1
    logical,optional                  :: dryrun
    integer                           :: i
    !put here a delete statement to avoid problems
    !
    if(sparse%status)call sp_delete_matrix_ll(sparse)
    !
    sparse%status=.true.
    !
    sparse%dryrun=.false.
    if(present(dryrun))sparse%dryrun=dryrun
    !
    select case (sparse%dryrun)
    case (.false.)
       sparse%Nrow=N       
    case (.true.)
       sparse%Nrow=1
    end select
    !
    sparse%Ncol=N
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%row(sparse%Nrow))
    do i=1,sparse%Nrow
       allocate(sparse%row(i)%root)
       sparse%row(i)%root%next => null()
       sparse%row(i)%size=0
       sparse%row(i)%max_column=0
    end do
#ifdef _MPI
    sparse%istart = 0
    sparse%iend   = 0
    sparse%ishift = 0
#endif
    !
  end subroutine sp_init_matrix_ll

#ifdef _MPI
  subroutine mpi_sp_init_matrix_ll(MpiComm,sparse,N,N1,dryrun)
    integer                           :: MpiComm
    type(sparse_matrix),intent(inout) :: sparse
    integer                           :: N
    integer,optional                  :: N1
    logical,optional                  :: dryrun
    logical                           :: dryrun_
    integer                           :: i,Ncol
    !
    if(sparse%status)call mpi_sp_delete_matrix_ll(MpiComm,sparse)
    !
    dryrun_=.false. ; if(present(dryrun))dryrun_=dryrun
    !
    call sp_MPI_setup(MpiComm,N)
    !
    Ncol = N
    if(present(N1))Ncol=N1
    !
    call sp_init_matrix_ll(sparse,MpiChunk,Ncol,dryrun_)
    !
    select case (sparse%dryrun)
    case (.false.)
       allocate(sparse%loc(sparse%Nrow))
    case (.true.)
       allocate(sparse%loc(1))
    end select
    do i=1,sparse%Nrow
       allocate(sparse%loc(i)%root)
       sparse%loc(i)%root%next => null()
       sparse%loc(i)%size=0
       sparse%loc(i)%max_column=0
    end do
    sparse%istart = MpiIstart
    sparse%iend   = MpiIend
    sparse%ishift = MpiIshift
  end subroutine mpi_sp_init_matrix_ll
#endif








  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_matrix_ll(sparse)    
    type(sparse_matrix),intent(inout) :: sparse
    integer                           :: i
    type(sparse_row),pointer          :: row
    type(sparse_element),pointer      :: p,c
    if(.not.sparse%status)stop "Warning SPARSE/sp_delete_matrix: sparse not allocated already."
    do i=1,sparse%Nrow
       row=>sparse%row(i)
       do
          p => row%root
          c => p%next
          if(.not.associated(c))exit  !empty list
          p%next => c%next !
          c%next=>null()
          deallocate(c)
       end do
       sparse%row(i)%size=0
       sparse%row(i)%max_column=0
       deallocate(sparse%row(i)%root)
    enddo
    deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%dryrun=.false.
    sparse%status=.false.
  end subroutine sp_delete_matrix_ll
  !
#ifdef _MPI
  subroutine mpi_sp_delete_matrix_ll(MpiComm,sparse)
    integer                           :: MpiComm
    type(sparse_matrix),intent(inout) :: sparse
    integer                           :: i
    type(sparse_row),pointer          :: row
    type(sparse_element),pointer      :: p,c
    if(.not.sparse%status)stop "Warning SPARSE/mpi_sp_delete_matrix: sparse not allocated already."
    do i=1,sparse%Nrow
       row=>sparse%row(i)
       do
          p => row%root
          c => p%next
          if(.not.associated(c))exit  !empty list
          p%next => c%next !
          c%next=>null()
          deallocate(c)
       end do
       sparse%row(i)%size=0
       sparse%row(i)%max_column=0
       deallocate(sparse%row(i)%root)
    enddo
    deallocate(sparse%row)
    !
    do i=1,sparse%Nrow
       row=>sparse%loc(i)
       do
          p => row%root
          c => p%next
          if(.not.associated(c))exit  !empty list
          p%next => c%next !
          c%next=>null()
          deallocate(c)
       end do
       sparse%loc(i)%size=0
       sparse%loc(i)%max_column=0
       deallocate(sparse%loc(i)%root)
    end do
    deallocate(sparse%loc)
    sparse%istart=0
    sparse%iend=0
    sparse%ishift=0
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%dryrun=.false.
    sparse%status=.false.
  end subroutine mpi_sp_delete_matrix_ll
#endif    




  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_element_c(sparse,value,i,j)
    type(sparse_matrix),intent(inout) :: sparse
    complex(8),intent(in)             :: value
    integer,intent(in)                :: i,j
    select case(sparse%dryrun)
    case (.false.)
       call sp_insert_element_main(sparse,value,i,j)
    case (.true.)
       call sp_insert_element_dryrun(sparse,j)
    end select
  end subroutine sp_insert_element_c


  subroutine sp_insert_element_main(sparse,value,i,j)
    type(sparse_matrix),intent(inout) :: sparse
    complex(8),intent(in)             :: value
    integer,intent(in)                :: i,j
    type(sparse_row),pointer          :: row
    integer                           :: column
    type(sparse_element),pointer      :: p,c
    logical                           :: iadd
    !
    column = j
    !
#ifdef _MPI
    if(column>=sparse%Istart.AND.column<=sparse%Iend)then
       row => sparse%loc(i)
    else
       row => sparse%row(i)
    endif
#else
    row => sparse%row(i)
#endif
    p => row%root
    c => p%next
    iadd = .false.                !check if column already exist
    do                            !traverse the list
       if(.not.associated(c))exit !empty list or end of the list
       if(c%col == column)then
          iadd=.true.
          exit
       endif
       !if(c%col > column)exit
       if(column <= c%col)exit
       p => c
       c => c%next
    end do
    if(iadd)then
       c%cval=c%cval + value
    else
       allocate(p%next)                !Create a new element in the list
       p%next%cval= value
       p%next%col = column
       row%size   = row%size+1
       if(column>row%max_column)row%max_column=column
       if(.not.associated(c))then !end of the list special case (current=>current%next)
          p%next%next  => null()
       else
          p%next%next  => c      !the %next of the new node come to current
       end if
    endif
  end subroutine sp_insert_element_main

  subroutine sp_insert_element_dryrun(sparse,j)
    type(sparse_matrix),intent(inout) :: sparse
    integer,intent(in)                :: j
    type(sparse_row),pointer          :: row
    integer                           :: column
    type(sparse_element),pointer      :: p,c
    logical                           :: iexist
    !
    column = j
    !
    row => sparse%row(1)
    p => row%root
    c => p%next
    iexist = .false.              !check if column already exist
    search_loop: do                            !traverse the list
       if(.not.associated(c))exit !empty list or end of the list
       if(c%col == column)then
          iexist=.true.
          exit search_loop
       endif
       if(column <= c%col)exit
       p => c
       c => c%next
    end do search_loop

    if(iexist)then
       !do nothin: this column is already present we do not need to store it
       return
    else
       allocate(p%next)                !Create a new element in the list
       p%next%col = column
       row%size   = row%size+1
       if(column>row%max_column)row%max_column=column
       if(.not.associated(c))then !end of the list special case (current=>current%next)
          p%next%next  => null()
       else
          p%next%next  => c      !the %next of the new node come to current
       end if
    endif
  end subroutine sp_insert_element_dryrun









  !+------------------------------------------------------------------+
  !PURPOSE: Get column list
  !+------------------------------------------------------------------+
  function sp_size_matrix(sparse) result(Nnz)
    type(sparse_matrix),intent(in) :: sparse
    integer                        :: i,Nnz
    type(sparse_row),pointer       :: row
    type(sparse_element),pointer   :: p,c
    if(.not.sparse%status)stop "Warning SPARSE/sp_size_matrix: sparse not allocated already."
    Nnz = 0
    do i=1,sparse%Nrow
       Nnz = Nnz + sparse%row(i)%size
    enddo
    if(Nnz==0)Nnz=1
  end function sp_size_matrix

  subroutine sp_columns_list_matrix(sparse,columns_vec)    
    type(sparse_matrix),intent(in)    :: sparse
    integer,dimension(:)              :: columns_vec
    integer                           :: i,Nsize,count
    type(sparse_row),pointer          :: row
    type(sparse_element),pointer      :: c
    if(.not.sparse%status)stop "Warning SPARSE/sp_columns_list_matrix: sparse not allocated already."
    Nsize=sp_size_matrix(sparse)
    if(size(columns_vec)/=Nsize)stop "ERROR SPARSE/sp_columns_list_matrix: size(columns_vec) incorrect"
    count=0
    columns_vec=0
    do i=1,sparse%Nrow
       row => sparse%row(i)
       c => row%root%next   !assume is associated,ie list exists
       do
          if(.not.associated(c))exit  !empty list
          count=count+1
          columns_vec(count)=c%col          
          c => c%next !          
       end do
    end do
  end subroutine sp_columns_list_matrix








  subroutine sp_columns_range_matrix_ll(sparse,columns_range)
    type(sparse_matrix),intent(in)   :: sparse
    integer,dimension(2)             :: columns_range
    integer,dimension(:),allocatable :: Vin,Vout,Vtmp
    integer                          :: Nnz
    integer                          :: k                   ! The number of unique elements
    integer                          :: i
    !
    Columns_Range = 0
    Nnz = sp_size_matrix(sparse)
    !
    allocate(Vin(Nnz))
    call sp_columns_list_matrix(sparse,Vin)
    !
    k = 1
    !
    allocate(Vtmp(Nnz))
    Vtmp = 0
    !
    Vtmp(1) = Vin(1)
    do i=2,size(Vin)
       ! if the number already exist in res check next
       if (any( Vtmp == Vin(i) )) cycle
       ! No match found so add it to the output
       k = k + 1
       Vtmp(k) = Vin(i)
    end do
    allocate(Vout(k))
    Vout = Vtmp(1:k)
    do i=0,MpiSize-1
       if(MpiRank==i)then
          write(*,"(A,2(I0,1x))")"Column list: ",minval(vout),maxval(vout)
       endif
    enddo
    !
    columns_range(1) = minval(Vout)
    columns_range(2) = maxval(Vout)
    !
  end subroutine sp_columns_range_matrix_ll

#ifdef _MPI
  subroutine mpi_sp_columns_range_matrix(MpiComm,sparse,columns_range)
    integer                            :: MpiComm
    type(sparse_matrix),intent(in)     :: sparse
    integer,dimension(0:,:)            :: columns_range
    integer,dimension(:,:),allocatable :: range_tmp
    integer,dimension(:),allocatable   :: Vin,Vout,Vtmp
    integer                            :: Nnz
    integer                            :: k,Ns !The number of unique elements
    integer                            :: i,Imin,Imax,Irank
    !
    Ns=0
    call Allreduce_MPI(MpiComm,sparse%Nrow,Ns)
    call sp_MPI_setup(MpiComm,Ns)
    !
    if( any(shape(columns_range) /= [MpiSize,2]) )stop "mpi_sp_columns_range_matrix: shape(columns_range)!= [MpiSize,2]"
    !
    Columns_Range = 0
    !
    Nnz = sp_size_matrix(sparse)
    !
    allocate(Vin(Nnz))
    call sp_columns_list_matrix(sparse,Vin)
    !
    k = 1
    !
    allocate(Vtmp(Nnz))
    Vtmp = 0
    !
    Vtmp(1) = Vin(1)
    do i=2,size(Vin)
       ! if the number already exist in res check next
       if (any( Vtmp == Vin(i) )) cycle
       ! No match found so add it to the output
       k = k + 1
       Vtmp(k) = Vin(i)
    end do
    allocate(Vout(k))
    Vout = Vtmp(1:k)
    !
    Imin = minval(Vout)/MpiQ
    Imax = maxval(Vout)/MpiQ ; if(Imax>=MpiSize-1)Imax=MpiSize-1
    !
    Columns_Range(MpiRank,1) = Imin
    Columns_Range(MpiRank,2) = Imax
    !
    ! do irank=0,MpiSize-1
    !    if(MpiRank==irank)then
    !       write(*,"(A,2(I0,1x))")"Column list     > ",minval(vout),maxval(vout)          
    !       write(*,"(A,2(I0,2x))")"Col owner Range > ",columns_range(irank,:)
    !       call MPI_Barrier(MpiComm,mpiierr)
    !    endif
    !    call MPI_Barrier(MpiComm,mpiierr)
    ! enddo
    ! call MPI_Barrier(MpiComm,mpiierr)
  end subroutine mpi_sp_columns_range_matrix
#endif








  !+------------------------------------------------------------------+
  !PURPOSE: load a regular matrix (2dim array) into a sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_load_matrix_c(matrix,sparse)
    complex(8),dimension(:,:),intent(in) :: matrix
    type(sparse_matrix),intent(inout)    :: sparse    
    integer                              :: i,j,Ndim1,Ndim2
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)   
    !
    if(sparse%status)call sp_delete_matrix_ll(sparse)
    call sp_init_matrix_ll(sparse,Ndim1,Ndim2)
    !
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=0.d0)call sp_insert_element_c(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine sp_load_matrix_c
  !
#ifdef _MPI
  subroutine mpi_sp_load_matrix_c(MpiComm,matrix,sparse)
    integer                              :: MpiComm
    complex(8),dimension(:,:),intent(in) :: matrix
    type(sparse_matrix),intent(inout)    :: sparse    
    integer                              :: i,j,Ndim1,Ndim2
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    call sp_MPI_setup(MpiComm,Ndim1) !split MPI ownership range along the rows of the matrix
    !
    if(sparse%status)call sp_delete_matrix_ll(sparse)
    call mpi_sp_init_matrix_ll(MpiComm,sparse,Ndim1,Ndim2)
    !
    do i=MpiIstart,MpiIend
       do j=1,Ndim2
          if(matrix(i,j)/=0.d0)call sp_insert_element_c(sparse,matrix(i,j),i-MpiIshift,j)
       enddo
    enddo
  end subroutine mpi_sp_load_matrix_c
#endif












  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sp_dump_matrix_c(sparse,matrix)
    type(sparse_matrix),intent(in)          :: sparse
    complex(8),dimension(:,:),intent(inout) :: matrix
    type(sparse_element),pointer            :: c
    integer                                 :: i,Ndim1,Ndim2
    !
    if(sparse%dryrun)stop "sp_dump_matrix_c ERROR: dumping with sparse%dryrun=T"
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_matrix: dimensions error"
    !
    matrix=0.d0
    do i=1,Ndim1
       c => sparse%row(i)%root%next
       do while(associated(c))
          matrix(i,c%col) = matrix(i,c%col) + c%cval
          c => c%next  !traverse list
       enddo
    enddo
  end subroutine sp_dump_matrix_c
  !
#ifdef _MPI
  subroutine mpi_sp_dump_matrix_c(MpiComm,sparse,matrix)
    integer                                 :: MpiComm
    type(sparse_matrix),intent(in)          :: sparse
    complex(8),dimension(:,:),intent(inout) :: matrix
    complex(8),dimension(:,:),allocatable   :: Mredux
    type(sparse_element),pointer            :: c
    integer                                 :: i,N1_,N2_,Ndim1,Ndim2,Nrow,Ncol
    !
    if(sparse%dryrun)stop "sp_dump_matrix_c ERROR: dumping with sparse%dryrun=T"
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)    
    !
    N1_=sparse%Nrow
    N2_=sparse%Ncol
    call MPI_AllReduce(N1_,Nrow,1,MPI_Integer,MPI_SUM,MpiComm,MpiIerr)
    call MPI_AllReduce(N2_,Ncol,1,MPI_Integer,MPI_MAX,MpiComm,MpiIerr)
    !
    if(Nrow>Ndim1 .OR. Ncol>Ndim2)stop "Warning SPARSE/mpi_dump_matrix: dimensions error"
    !
    call sp_MPI_setup(MpiComm,Ndim1)
    !
    allocate(Mredux(Ndim1,Ndim2))
    Mredux=dcmplx(0d0,0d0)
    do i=MpiIstart,MpiIend
       c => sparse%loc(i-MpiIshift)%root%next
       do while(associated(c))
          Mredux(i,c%col) = Mredux(i,c%col) + c%cval
          c => c%next  !traverse list
       enddo
       !
       c => sparse%row(i-MpiIshift)%root%next
       do while(associated(c))
          Mredux(i,c%col) = Mredux(i,c%col) + c%cval
          c => c%next  !traverse list
       enddo
    enddo
    matrix=dcmplx(0d0,0d0)
    call MPI_AllReduce(Mredux,Matrix,Ndim1*Ndim2,MPI_Double_Complex,MPI_Sum,MpiComm,MpiIerr)
    deallocate(Mredux)
  end subroutine mpi_sp_dump_matrix_c
#endif






















  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+
  subroutine sp_print_matrix_ll(sparse,unit,fmt)
    type(sparse_matrix)          :: sparse
    integer,optional             :: unit
    integer                      :: unit_
    integer                      :: i,j,Ns
    character(len=*),optional    :: fmt
    character(len=64)            :: fmt_
    type(sparse_row),pointer     :: row
    type(sparse_element),pointer :: c
    integer                      :: count=0
    unit_=6;if(present(unit))unit_=unit
    fmt_='2F15.9';if(present(fmt))fmt_=fmt
    write(*,*)"Print sparse matrix (compact mode) ->",unit_
    do i=1,sparse%Nrow
       row => sparse%row(i)
       c => row%root%next   !assume is associated,ie list exists
       do
          if(.not.associated(c))exit
          count=count+1
          write(unit_,"("//trim(fmt_)//",A1,I0,3X)",advance='no')dreal(c%cval),dimag(c%cval),',',c%col
          c => c%next  !traverse list
       end do
       write(unit_,*)
#ifdef _MPI
       row => sparse%loc(i)
       c => row%root%next   !assume is associated,ie list exists
       do
          if(.not.associated(c))exit
          count=count+1
          write(unit_,"("//trim(fmt_)//",A1,I0,3X)",advance='no')c%cval,',',c%col
          c => c%next  !traverse list
       end do
       write(unit_,*)
#endif
    enddo
    write(unit_,*)
  end subroutine sp_print_matrix_ll

  subroutine sp_spy_matrix_ll(sparse,header)
    type(sparse_matrix)          :: sparse
    character ( len = * )        :: header
    integer                      :: N1,N2
    type(sparse_element),pointer :: c
    character ( len = 255 )      :: command_filename
    integer                      :: command_unit
    character ( len = 255 )      :: data_filename
    integer                      :: data_unit
    integer                      :: i, j
    character ( len = 6 )        :: n1_s,n2_s,n1_i,n2_i
    integer                      :: nz_num
    character ( len = 255 )      :: png_filename
    !
    !  Create data file.
    !
    !
    N1 = sparse%Nrow
    N2 = sparse%Ncol
    data_filename = trim ( header ) // '_data.dat'
    open (unit=free_unit(data_unit), file = data_filename, status = 'replace' )
    nz_num = 0
    do i=1,N1
       c => sparse%row(i)%root%next
       do while(associated(c))
          write(data_unit,'(2x,i6,2x,i6)') c%col,i
          nz_num = nz_num + 1
          c => c%next  !traverse list
       enddo
    enddo
    close(data_unit)
    !
    !  Create command file.
    !
    command_filename = "plot_"//str(header)//'_commands.gp'
    open(unit = free_unit(command_unit), file = command_filename, status = 'replace' )
    write(command_unit,'(a)') '#unset key'
    write(command_unit,'(a)') 'set terminal postscript eps enhanced color font "Times-Roman,16"'
    write(command_unit,'(a)') 'set output "|ps2pdf -sEPSCrop - '//str(header)//".pdf"//'"'
    write(command_unit,'(a)') 'set size ratio -1'
    write(command_unit,'(a)') 'set xlabel "<--- J --->"'
    write(command_unit,'(a)') 'set ylabel "<--- I --->"'
    write(command_unit,'(a,i6,a)')'set title "',nz_num,' nonzeros for '//str(header)//'"'
    write(command_unit,'(a)') 'set timestamp'
    write(command_unit,'(a)' )'plot [x=1:'//str(N1)//'] [y='//str(N2)//':1] "'//&
         str(data_filename)//'" w p pt 5 ps 0.4 lc rgb "red"'
    close ( unit = command_unit )
    return
  end subroutine sp_spy_matrix_ll
#ifdef _MPI
  subroutine mpi_sp_spy_matrix_ll(MpiComm,sparse,header)
    integer                      :: MpiComm
    type(sparse_matrix)          :: sparse
    character ( len = * )        :: header
    integer                      :: N1,N2,N1_,N2_
    type(sparse_element),pointer :: c
    character ( len = 255 )      :: command_filename
    integer                      :: command_unit
    character ( len = 255 )      :: data_filename(2)
    integer                      :: data_unit
    integer                      :: i, j
    character ( len = 6 )        :: n1_s,n2_s,n1_i,n2_i
    integer                      :: nz_num
    character ( len = 255 )      :: png_filename
    !
    !  Create data file.
    !
    N1_=sparse%Nrow
    N2_=sparse%Ncol
    call MPI_AllReduce(N1_,N1,1,MPI_Integer,MPI_SUM,MpiComm,MpiIerr)
    call MPI_AllReduce(N2_,N2,1,MPI_Integer,MPI_MAX,MpiComm,MpiIerr)
    !
    call sp_MPI_setup(MpiComm,N1)
    !
    nz_num = 0
    !
    data_filename(1) = trim(header)//"_rank"//str(MpiRank,4)//'_matrix.dat'
    open(unit=free_unit(data_unit),file=data_filename(1), status = 'replace' )
    do i=1,MpiChunk
       c => sparse%row(i)%root%next
       do while(associated(c))
          write(data_unit,'(2x,i6,2x,i6)') c%col,i+MpiIshift
          nz_num = nz_num + 1
          c => c%next  !traverse list
       enddo
    enddo
    write(data_unit,'(2x,i6,2x,i6)')
    close(data_unit)
    !
    data_filename(2) = trim(header)//"_rank"//str(MpiRank,4)//'_local.dat'
    open(unit=free_unit(data_unit),file=data_filename(2), status = 'replace' )
    do i=1,MpiChunk
       c => sparse%loc(i)%root%next
       do while(associated(c))
          write(data_unit,'(2x,i6,2x,i6)') c%col,i+MpiIshift
          nz_num = nz_num + 1
          c => c%next  !traverse list
       enddo
    enddo
    write(data_unit,'(2x,i6,2x,i6)')
    close(data_unit)
    !
    call MPI_Barrier(MpiComm,MpiIerr)
    !
    !  Create command file.
    !
    command_filename = "plot_"//trim(header)//"_rank"//str(MpiRank,4)//'_commands.gp'
    open(unit = free_unit(command_unit), file = command_filename, status = 'replace' )
    write(command_unit,'(a)') '#unset key'
    write(command_unit,'(a)') 'set terminal postscript eps enhanced color font "Times-Roman,16"'
    write(command_unit,'(a)') 'set output "|ps2pdf -sEPSCrop - '//str(header)//"_rank"//str(MpiRank,4)//".pdf"//'"'
    write(command_unit,'(a)') 'set size ratio -1'
    write(command_unit,'(a)') 'set xlabel "<--- J --->"'
    write(command_unit,'(a)') 'set ylabel "<--- I --->"'
    write(command_unit,'(a,i6,a)' ) &
         'set title "',nz_num,' nonzeros for '//str(header)//"_rank"//str(MpiRank,4)//'"'
    write(command_unit,'(a)') 'set timestamp'
    write(command_unit,'(a)' )'plot [x=1:'//str(N1)//'] [y='//str(N2)//':1] "'//&
         str(data_filename(1))//'" w p pt 5 ps 0.4 lc rgb "red" title "Non-Local", "'//&
         str(data_filename(2))//'" w p pt 5 ps 0.4 lc rgb "blue" title "Local"'
    close ( unit = command_unit )
    return
  end subroutine mpi_sp_spy_matrix_ll
#endif










#ifdef _MPI
  subroutine sp_MPI_setup(MpiComm,N)
    integer :: MpiComm,N,irank
    !
    MpiRank   = get_rank_MPI(MpiComm)
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    MpiQ      = get_Q_MPI(MpiComm,N)
    MpiR      = get_R_MPI(MpiComm,N)
    MpiChunk  = MpiQ+MpiR
    MpiIstart = 1+MpiRank*MpiQ
    MpiIend   = (MpiRank+1)*MpiQ+MpiR
    MpiIshift = MpiRank*MpiQ
  end subroutine sp_MPI_setup
#endif







end module ED_SPARSE_MATRIX










!   function sp_size_local(sparse) result(Nnz)
!     type(sparse_matrix),intent(in) :: sparse
!     integer                        :: i,Nnz
!     type(sparse_row),pointer       :: row
!     type(sparse_element),pointer   :: p,c
!     if(.not.sparse%status)stop "Warning SPARSE/sp_size_local: sparse not allocated already."
!     Nnz = 0
! #ifdef _MPI
!     do i=1,sparse%Nrow
!        Nnz = Nnz + sparse%loc(i)%size
!     enddo
!     if(Nnz==0)Nnz=1
! #endif
!   end function sp_size_local

!   subroutine sp_columns_list_local(sparse,columns_vec)    
!     type(sparse_matrix),intent(in)    :: sparse
!     integer,dimension(:)              :: columns_vec
!     integer                           :: i,Nsize,count
!     type(sparse_row),pointer          :: row
!     type(sparse_element),pointer      :: c
!     if(.not.sparse%status)stop "Warning SPARSE/sp_columns_list_local: sparse not allocated already."
!     Nsize=sp_size_local(sparse)
!     if(size(columns_vec)/=Nsize)stop "ERROR SPARSE/sp_columns_list_local: size(columns_vec) incorrect"
!     columns_vec=0
! #ifdef _MPI
!     count=0
!     do i=1,sparse%Nrow
!        row => sparse%loc(i)
!        c => row%root%next   !assume is associated,ie list exists
!        do
!           if(.not.associated(c))exit  !empty list
!           count=count+1
!           columns_vec(count)=c%col          
!           c => c%next !          
!        end do
!     end do
! #endif
!   end subroutine sp_columns_list_local




! !GET NUMBER OF NON-ZERO ELEMENTS
! interface sp_get_nnz
!    module procedure sp_get_nnz_ll
! end interface sp_get_nnz
! public :: sp_get_nnz



! !INSERT DIAGONAL ENTRY IN LL-SPARSE MATRIX
! interface sp_insert_diag
!    module procedure sp_insert_diag_c
! end interface sp_insert_diag
! public :: sp_insert_diag      !insert a vector at diag  !checked



! !GET ELEMENTS ALONG THE DIAGONAL
! interface sp_get_diagonal
!    module procedure sp_get_diagonal_c
! end interface sp_get_diagonal
! public :: sp_get_diagonal     !get diagonal elements    !checked
















! !+------------------------------------------------------------------+
! !PURPOSE: delete a single element at (i,j) from the sparse matrix
! !+------------------------------------------------------------------+
! subroutine sp_delete_element(matrix,i,j)
!   type(sparse_matrix),intent(inout) :: matrix
!   integer,intent(in)                :: i,j
!   logical :: delete
!   delete = delete_element_from_row(matrix%row(i),col=j)
!   if(.not.delete)write(*,"(A,I3,I3)")"sp_delete_element: can not delete element in",i,j
! end subroutine sp_delete_element

! !+------------------------------------------------------------------+
! !PURPOSE: delete an entire row from the sparse matrix (private)
! !+------------------------------------------------------------------+
! subroutine delete_row(row)
!   type(sparse_row),intent(inout) :: row
!   type(sparse_element),pointer   :: p,c
!   do
!      p => row%root
!      c => p%next
!      if(.not.associated(c))exit  !empty list
!      p%next => c%next !
!      c%next=>null()
!      deallocate(c)
!   end do
! end subroutine delete_row



! !This shoud be better tested!
! !+------------------------------------------------------------------+
! !PURPOSE: delete a given element from a row of the sparse matrix (private)
! !+------------------------------------------------------------------+
! function delete_element_from_row(row,n,col) result(found)
!   type(sparse_row),intent(inout)    :: row
!   integer,optional                  :: n
!   integer,optional                  :: col
!   integer                           :: i,pos
!   type(sparse_element),pointer      :: p,c
!   logical                           :: found
!   pos= row%size ; if(present(n))pos=n
!   p => row%root
!   c => p%next
!   found = .false.
!   if(present(col))then
!      do 
!         if(found .OR. .not.associated(c))return
!         if(col == c%col)then
!            found=.true.
!            exit
!         else
!            p => c
!            c => c%next
!         endif
!      end do
!      if(found)then
!         p%next => c%next !reallocate skipping the deleted link
!         deallocate(c)           !free link
!         row%size=row%size-1
!      endif
!   else
!      do i=1,pos 
!         if(.not.associated(c))return !empty list
!         p => c
!         c => c%next
!      end do
!      found=.true.
!      p%next => c%next !reallocate skipping the deleted link
!      deallocate(c)           !free link
!      row%size=row%size-1
!   endif
! end function delete_element_from_row










! !+------------------------------------------------------------------+
! !PURPOSE:  return total number of non-zero elements stored in sparse
! !+------------------------------------------------------------------+
! function sp_get_nnz_ll(sparse) result(Nnz)
!   type(sparse_matrix) :: sparse
!   integer             :: i
!   integer             :: Nnz
!   Nnz=0
!   do i=1,sparse%Nrow
!      Nnz=Nnz+sparse%row(i)%size
!   enddo
! end function sp_get_nnz_ll












! !+------------------------------------------------------------------+
! !PURPOSE: insert a vector of elements at the diagonal of the sparse matrix
! !+------------------------------------------------------------------+
! ! subroutine sp_insert_diag_d(sparse,diag)
! !   type(sparse_matrix),intent(inout)  :: sparse
! !   real(8),intent(in),dimension(:)    :: diag
! !   integer                            :: i
! !   if(size(diag)/=sparse%Nrow)stop "sp_insert_diag: error in dimensions"
! !   do i=1,size(diag)
! !      call insert_element_in_row_d(sparse%row(i),diag(i),i)
! !   enddo
! ! end subroutine sp_insert_diag_d
! subroutine sp_insert_diag_c(sparse,diag)
!   type(sparse_matrix),intent(inout)  :: sparse
!   complex(8),intent(in),dimension(:) :: diag
!   integer                            :: i
!   if(size(diag)/=sparse%Nrow)stop "sp_insert_diag: error in dimensions"
!   do i=1,size(diag)
!      call insert_element_in_row_c(sparse%row(i),diag(i),i)
!   enddo
! end subroutine sp_insert_diag_c




! !+------------------------------------------------------------------+
! !PURPOSE: insert an element in a given row (private) 
! !+------------------------------------------------------------------+
! subroutine insert_element_in_row_c(row,value,column)
!   type(sparse_row),intent(inout)    :: row
!   complex(8) ,intent(in)            :: value
!   integer, intent(in)               :: column
!   type(sparse_element),pointer      :: p,c
!   logical :: iadd
!   p => row%root
!   c => p%next
!   iadd = .false.                !check if column already exist
!   do                            !traverse the list
!      if(.not.associated(c))exit !empty list or end of the list
!      if(c%col == column)then
!         iadd=.true.
!         exit
!      endif
!      !if(c%col > column)exit
!      if(column <= c%col)exit
!      p => c
!      c => c%next
!   end do
!   if(iadd)then
!      c%cval=c%cval + value
!   else
!      allocate(p%next)                !Create a new element in the list
!      p%next%cval= value
!      p%next%col = column
!      row%size   = row%size+1
!      if(.not.associated(c))then !end of the list special case (current=>current%next)
!         p%next%next  => null()
!      else
!         p%next%next  => c      !the %next of the new node come to current
!      end if
!   endif
! end subroutine insert_element_in_row_c






! !+------------------------------------------------------------------+
! !PURPOSE: get the diagonal elements of the sparse matrix
! !+------------------------------------------------------------------+
! ! subroutine sp_get_diagonal_d(sparse,diag)
! !   type(sparse_matrix),intent(inout) :: sparse
! !   real(8),dimension(:)              :: diag
! !   integer                           :: Ndim,i
! !   Ndim=size(diag);if(Ndim/=sparse%Nrow)stop "sp_get_diagonal: error in diag dimension." 
! !   do i=1,Ndim
! !      call get_element_from_row_d(sparse%row(i),diag(i),i)
! !   enddo
! ! end subroutine  sp_get_diagonal_d
! subroutine sp_get_diagonal_c(sparse,diag)
!   type(sparse_matrix),intent(inout) :: sparse
!   complex(8),dimension(:)           :: diag
!   integer                           :: Ndim,i
!   Ndim=size(diag);if(Ndim/=sparse%Nrow)stop "sp_get_diagonal: error in diag dimension." 
!   do i=1,Ndim
!      call get_element_from_row_c(sparse%row(i),diag(i),i)
!   enddo
! end subroutine sp_get_diagonal_c











! !+------------------------------------------------------------------+
! !PURPOSE: get an element from position (i,j) of the sparse matrix
! !+------------------------------------------------------------------+
! ! function sp_get_element_d(sparse,i,j) result(value)
! !   type(sparse_matrix),intent(inout) :: sparse    
! !   integer,intent(in)                :: i,j
! !   real(8)                           :: value
! !   call get_element_from_row_d(sparse%row(i),value,j)
! ! end function sp_get_element_d
! function sp_get_element_c(sparse,i,j) result(value)
!   type(sparse_matrix),intent(inout) :: sparse    
!   integer,intent(in)                :: i,j
!   complex(8)                        :: value
!   call get_element_from_row_c(sparse%row(i),value,j)
! end function sp_get_element_c



! !+------------------------------------------------------------------+
! !PURPOSE: get an element from a given row of the matrix (private)
! !+------------------------------------------------------------------+
! subroutine get_element_from_row_c(row,value,column)
!   type(sparse_row),intent(inout)    :: row
!   complex(8)                        :: value
!   integer, intent(in)               :: column
!   type(sparse_element),pointer      :: c
!   c => row%root%next
!   value=cmplx(0.d0,0.d0,8)
!   do                            !traverse the list
!      if(.not.associated(c))return !empty list or end of the list
!      if(c%col == column)exit
!      c => c%next
!   end do
!   !
!   value = c%cval
! end subroutine get_element_from_row_c












! !+------------------------------------------------------------------+
! !PURPOSE: check if a given element exists
! !+------------------------------------------------------------------+
! function sp_inquire_element(sparse,i,j) result(exist)
!   type(sparse_matrix),intent(inout) :: sparse    
!   integer,intent(in)                :: i,j
!   logical                           :: exist
!   exist = inquire_element_from_row(sparse%row(i),j)
! end function sp_inquire_element

! !+------------------------------------------------------------------+
! !PURPOSE: check if an element in a given row of the matrix exist (private)
! !+------------------------------------------------------------------+
! function inquire_element_from_row(row,column) result(exist)
!   type(sparse_row),intent(inout)    :: row
!   logical                           :: exist
!   integer, intent(in)               :: column
!   type(sparse_element),pointer      :: c
!   c => row%root%next
!   exist=.false.
!   do                            !traverse the list
!      if(.not.associated(c))return !empty list or end of the list
!      if(c%col == column)exit
!      c => c%next
!   end do
!   exist=.true.
! end function inquire_element_from_row

























! subroutine sp_print_matrix_ll(sparse,unit,fmt,full)
!   type(sparse_matrix)            :: sparse
!   integer,optional               :: unit
!   integer                        :: i,j,unit_,Ns
!   character(len=*),optional      :: fmt
!   character(len=64)              :: fmt_
!   logical,optional               :: full
!   logical                        :: full_
!   unit_=6;if(present(unit))unit_=unit
!   fmt_='F6.2';if(present(fmt))fmt_=fmt
!   full_=.false.;if(present(full))full_=full
!   Ns=sparse%Nrow
!   if(full_)then
!      write(*,*)"Print sparse matrix (full mode < 100) ->",unit_
!      do i=1,Ns
!         write(unit_,"(100("//trim(fmt_)//",A1,"//trim(fmt_)//",2X))")(&
!              real(sp_get_element_c(sparse,i,j)),",",imag(sp_get_element_c(sparse,i,j)),j=1,Ns)
!      enddo
!   else
!      write(*,*)"Print sparse matrix (compact mode) ->",unit_
!      do i=1,Ns
!         call print_row_c(sparse%row(i),unit_,fmt_)
!      enddo
!   endif
!   write(unit_,*)
! end subroutine sp_print_matrix_ll









! !+------------------------------------------------------------------+
! !PURPOSE: print an entire row of the sparse matrix (private)
! !+------------------------------------------------------------------+
! subroutine print_row_c(row,unit,fmt)
!   type(sparse_row),intent(in)   :: row
!   type(sparse_element),pointer  :: c
!   integer                       :: count=0
!   integer,optional :: unit
!   integer          :: unit_
!   character(len=*),optional :: fmt
!   character(len=64)         :: fmt_
!   unit_=6;if(present(unit))unit_=unit
!   fmt_='F15.9';if(present(fmt))fmt_=fmt
!   c => row%root%next   !assume is associated,ie list exists
!   do
!      if(.not.associated(c))exit
!      count=count+1
!      write(unit_,"(2"//trim(fmt_)//",A1,I3,3X)",advance='no')c%cval,',',c%col
!      c => c%next  !traverse list
!   end do
!   write(unit_,*)
! end subroutine print_row_c






! !+-----------------------------------------------------------------------------+!
! !PURPOSE:  test if a sparse matrix is symmetric 
! !+-----------------------------------------------------------------------------+
! subroutine sp_test_symmetric(sparse)
!   type(sparse_matrix)                   :: sparse
!   logical                               :: is_symmetric
!   complex(8),dimension(:,:),allocatable :: cM
!   integer                               :: Nrow,Ncol
!   Nrow=sparse%Nrow
!   Ncol=Nrow
!   is_symmetric=.false.
!   allocate(cM(Nrow,Ncol))
!   call sp_dump_matrix(sparse,cM)
!   if( maxval(abs(cM-conjg(transpose(cM))) ) < 1.d-12)is_symmetric=.true.
!   if(is_symmetric)then
!      write(*,"(A)")"Matrix IS Hermitian"
!   else
!      write(*,"(A)")"Matrix IS NOT Hermitian"
!   endif
! end subroutine sp_test_symmetric



