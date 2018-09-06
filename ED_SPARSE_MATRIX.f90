MODULE ED_SPARSE_MATRIX  !THIS VERSION CONTAINS ONLY COMPLX ELEMENT: (HERMITIAN MATRIX)
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private


  !SPARSE MATRIX: LINKED LIST FORMAT:
  type sparse_element_ll
     complex(8)                      :: cval !value of the entry: double complex
     integer                         :: col  !col connected to this compress value
     type(sparse_element_ll),pointer :: next !link to next entry in the row
  end type sparse_element_ll

  type sparse_row_ll
     integer                         :: size=0    !size of the list
     integer                         :: min_column= huge(1)
     integer                         :: max_column=-huge(1)
     type(sparse_element_ll),pointer :: root    !head/root of the list\== list itself
  end type sparse_row_ll

  type sparse_matrix_ll
     integer                                  :: Nrow
     integer                                  :: Ncol
     logical                                  :: status=.false.
     type(sparse_row_ll),dimension(:),pointer :: row
#ifdef _MPI
     logical                                  :: mpi=.false.
     integer                                  :: istart=0 !global start index for MPI storage
     integer                                  :: iend=0
     integer                                  :: ishift=0
#endif
  end type sparse_matrix_ll




  !INIT SPARSE MATRICES 
  interface sp_init_matrix
     module procedure :: sp_init_matrix_ll
#ifdef _MPI
     module procedure :: mpi_sp_init_matrix_ll
#endif
  end interface sp_init_matrix



  !DELETE SPARSE MATRIX 
  interface sp_delete_matrix
     module procedure :: sp_delete_matrix_ll
#ifdef _MPI
     module procedure :: mpi_sp_delete_matrix_ll
#endif
  end interface sp_delete_matrix


  !INSERT ELEMENTS
  interface sp_insert_element
     module procedure :: sp_insert_element_ll
#ifdef _MPI
     module procedure :: mpi_sp_insert_element_ll
#endif
  end interface sp_insert_element



  !DUMP SPARSE MATRIX INTO STANDARD MATRIX
  interface sp_dump_matrix
     module procedure :: sp_dump_matrix_ll
#ifdef _MPI
     module procedure :: mpi_sp_dump_matrix_ll
#endif
  end interface sp_dump_matrix




  !Linked-List Sparse Matrix
  public :: sparse_element_ll
  public :: sparse_matrix_ll

  public :: sp_init_matrix      !init the sparse matrix   !checked
  public :: sp_delete_matrix    !delete the sparse matrix !checked
  public :: sp_insert_element   !insert an element        !checked
  public :: sp_dump_matrix      !dump sparse into array   !checked
#ifdef _MPI
  public :: sp_set_mpi_ll
#endif



  integer :: MpiRank=0
  integer :: MpiSize=1
  integer :: MpiIerr
  logical :: MpiMaster=.true.






contains       


  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  subroutine sp_init_matrix_ll(sparse,N,N1)
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                           :: N
    integer,optional                  :: N1
    integer                           :: i
    !put here a delete statement to avoid problems
    if(sparse%status)stop "sp_init_matrix LL: alreay allocate can not init"
    sparse%Nrow=N
    sparse%Ncol=N
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%row(N))
    do i=1,N
       allocate(sparse%row(i)%root)
       sparse%row(i)%root%next => null()
       sparse%row(i)%size=0
       sparse%row(i)%min_column=huge(1)
       sparse%row(i)%max_column=-huge(1)
    end do
    !
    sparse%status=.true.
    !
  end subroutine sp_init_matrix_ll



#ifdef _MPI
  subroutine mpi_sp_init_matrix_ll(MpiComm,sparse,N,N1)
    integer                              :: MpiComm
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                              :: N
    integer,optional                     :: N1
    integer                              :: i,Ncol,Nloc
    !
    call sp_test_matrix_mpi(MpiComm,sparse,"mpi_sp_init_matrix_ll")
    !
    Ncol = N
    if(present(N1))Ncol=N1
    Nloc = sparse%iend-sparse%istart+1
    call sp_init_matrix_ll(sparse,Nloc,Ncol)
    !
  end subroutine mpi_sp_init_matrix_ll
#endif







  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_matrix_ll(sparse)    
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                           :: i
    type(sparse_row_ll),pointer          :: row
    type(sparse_element_ll),pointer      :: p,c
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
       sparse%row(i)%min_column=huge(1)
       sparse%row(i)%max_column=-huge(1)
       deallocate(sparse%row(i)%root)
    end do
    deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%status=.false.
#ifdef _MPI
    sparse%istart = 0
    sparse%iend   = 0
    sparse%ishift = 0
    sparse%mpi    = .false.
#endif 
  end subroutine sp_delete_matrix_ll


#ifdef _MPI
  subroutine mpi_sp_delete_matrix_ll(MpiComm,sparse)
    integer                              :: MpiComm
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                              :: i
    type(sparse_row_ll),pointer          :: row
    type(sparse_element_ll),pointer      :: p,c
    if(.not.sparse%status)stop "Error SPARSE/mpi_sp_delete_matrix: sparse is not allocated."
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
       sparse%row(i)%min_column=huge(1)
       sparse%row(i)%max_column=-huge(1)
       deallocate(sparse%row(i)%root)
       end do
    deallocate(sparse%row)
    !
    sparse%istart=0
    sparse%iend=0
    sparse%ishift=0
    sparse%mpi=.false.
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%status=.false.
  end subroutine mpi_sp_delete_matrix_ll
#endif    










  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_element_ll(sparse,value,i,j)
    type(sparse_matrix_ll),intent(inout) :: sparse
    complex(8),intent(in)                :: value
    integer,intent(in)                   :: i,j
    type(sparse_row_ll),pointer          :: row
    integer                              :: column
    type(sparse_element_ll),pointer      :: p,c
    logical                              :: iadd
    !
    column = j
    !
    row => sparse%row(i)
    !
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
       if(column<row%min_column)row%min_column=column
       if(column>row%max_column)row%max_column=column
       if(.not.associated(c))then !end of the list special case (current=>current%next)
          p%next%next  => null()
       else
          p%next%next  => c      !the %next of the new node come to current
       end if
    endif
  end subroutine sp_insert_element_ll



#ifdef _MPI
  subroutine mpi_sp_insert_element_ll(MpiComm,sparse,value,i,j)
    integer                              :: MpiComm
    type(sparse_matrix_ll),intent(inout) :: sparse
    complex(8),intent(in)                :: value
    integer,intent(in)                   :: i,j
    type(sparse_row_ll),pointer          :: row
    integer                              :: column
    type(sparse_element_ll),pointer      :: p,c
    logical                              :: iadd
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_insert_element_ll")
    !
    column = j
    !
    row => sparse%row(i-sparse%Ishift)
    !
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
       if(column<row%min_column)row%min_column=column
       if(column>row%max_column)row%max_column=column
       if(.not.associated(c))then !end of the list special case (current=>current%next)
          p%next%next  => null()
       else
          p%next%next  => c      !the %next of the new node come to current
       end if
    endif
  end subroutine mpi_sp_insert_element_ll
#endif

    !
    !
    !












  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sp_dump_matrix_ll(sparse,matrix)
    type(sparse_matrix_ll),intent(in)       :: sparse
    complex(8),dimension(:,:),intent(inout) :: matrix
    type(sparse_element_ll),pointer         :: c
    integer                                 :: i,Ndim1,Ndim2
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
  end subroutine sp_dump_matrix_ll


  !
#ifdef _MPI
  subroutine mpi_sp_dump_matrix_ll(MpiComm,sparse,matrix)
    integer                                 :: MpiComm
    type(sparse_matrix_ll),intent(in)       :: sparse
    complex(8),dimension(:,:),intent(inout) :: matrix
    type(sparse_element_ll),pointer         :: c
    integer                                 :: i,N1_,N2_,Ndim1,Ndim2,Nrow,Ncol
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_dump_matrix_ll")
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    N1_=sparse%Nrow
    N2_=sparse%Ncol
    Nrow=0
    Ncol=0
    call MPI_AllReduce(N1_,Nrow,1,MPI_Integer,MPI_SUM,MpiComm,MpiIerr)
    call MPI_AllReduce(N2_,Ncol,1,MPI_Integer,MPI_MAX,MpiComm,MpiIerr)
    !
    if(Nrow>Ndim1 .OR. Ncol>Ndim2)stop "Warning SPARSE/mpi_dump_matrix: dimensions error"
    ! !
    matrix=0d0
    do i=sparse%Istart,sparse%Iend
       c => sparse%row(i-sparse%Ishift)%root%next
       do while(associated(c))
          matrix(i,c%col) = matrix(i,c%col) + c%cval
          c => c%next  !traverse list
       enddo
    enddo
  end subroutine mpi_sp_dump_matrix_ll
#endif



















#ifdef _MPI
  subroutine sp_set_mpi_ll(MpiComm,sparse,istart,iend,ishift)
    integer                              :: MpiComm
    type(sparse_matrix_ll),intent(inout) :: sparse
    integer                              :: istart,iend,ishift
    sparse%istart = istart
    sparse%iend   = iend
    sparse%ishift = ishift
    sparse%mpi    = .true.
  end subroutine sp_set_mpi_ll

  subroutine sp_test_matrix_mpi(MpiComm,sparse,text)
    integer                              :: MpiComm
    type(sparse_matrix_ll),intent(in)    :: sparse
    character(len=*)                     :: text
    integer                              :: MpiRank
    MpiRank = get_Rank_MPI(MpiComm)
    if(.not.sparse%mpi)then
       print*,"Rank, Error in "//trim(text)//": mpi no set"
       stop
    endif
  end subroutine sp_test_matrix_mpi
#endif


end module ED_SPARSE_MATRIX
