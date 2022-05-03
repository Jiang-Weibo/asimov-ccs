submodule (partitioning) partitioning_parhip

  use kinds, only: ccs_int, ccs_real
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

contains

  !v Partition the mesh
  !
  ! Use ParHIP library to compute a k-way vertex separator given a k-way partition of the graph.
  ! The graph can be weighted or unweighted.
  module subroutine partition_kway(par_env, topo)

    ! use iso_c_binding

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition

    ! Local variables
    real(ccs_real) :: imbalance
    integer(ccs_int) :: seed
    integer(ccs_int) :: mode
    integer(ccs_int) :: suppress
    integer(ccs_int) :: edgecut
    ! type(c_ptr) :: vwgt     ! NULL pointer to be used if no weights are attached to vertices
    ! type(c_ptr) :: adjwgt   ! NULL pointer to be used if no weights are attached to edges

    ! Values hardcoded for now
    imbalance = 0.03  ! Desired balance - 0.03 = 3% 
    seed = 2022       ! "Random" seed
    mode = 3          ! FASTSOCIAL
    suppress = 0      ! Do not suppress the output
    edgecut = -1      ! XXX: silence unused variable warning

    ! ParHIP needs 0-indexing - shift array contents by -1
    topo%vtxdist = topo%vtxdist - 1
    topo%xadj = topo%xadj - 1
    topo%adjncy = topo%adjncy - 1

    ! Set weights to 1 - replace with NULL pointers later
    topo%adjwgt = 1
    topo%vwgt = 1

    ! NULL pointers
    ! vwgt = c_null_ptr
    ! adjwgt = c_null_ptr
  
    ! Partitioning an unweighted graph
    select type(par_env)
    type is (parallel_environment_mpi)

      ! call partition_parhipkway(topo%vtxdist, topo%xadj, topo%adjncy, topo%vwgt, topo%adjwgt, & 
      !                             par_env%num_procs, imbalance, suppress, seed, &
      !                             mode, edgecut, topo%part, par_env%comm)

    class default
      print*, "ERROR: Unknown parallel environment!"
    end select    

    ! Return to 1-indexing by adding 1. For safety only right now, later on 
    ! these arrays will be deallocated as they are not used beyond this point
    topo%vtxdist = topo%vtxdist + 1
    topo%xadj = topo%xadj + 1
    topo%adjncy = topo%adjncy + 1

    print*,"Number of edgecuts:",edgecut
    ! Make sure partition vector starts from 1
    topo%part = topo%part + 1

  end subroutine partition_kway

  !v Compute the input arrays for the partitioner
  !
  ! Using the topology object, compute the input arrays for the ParHIP partitioner
  ! Input arrays for the partitioner are: vtxdist, xadj and adjncy
  module subroutine compute_partitioner_input(par_env, topo)

    class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
    type(topology), target, intent(inout) :: topo                           !< The topology for which to compute the parition

    ! Local variables
    integer(ccs_int), dimension(:,:), allocatable :: wkint2d ! Temporary 2D integer array

    integer(ccs_int) :: i, j, k
    integer(ccs_int) :: irank ! MPI rank ID
    integer(ccs_int) :: isize ! Size of MPI world
    integer(ccs_int) :: start_index 
    integer(ccs_int) :: end_index  
    integer(ccs_int) :: face_nb1  
    integer(ccs_int) :: face_nb2
    integer(ccs_int) :: local_index
    integer(ccs_int) :: num_connections
 
    irank = par_env%proc_id
    isize = par_env%num_procs
  
    ! Create and populate the vtxdist array based on the total number of cells
    ! and the total number of ranks in the parallel environment
    allocate(topo%vtxdist(isize + 1)) ! vtxdist array is of size num_procs + 1 on all ranks

    topo%vtxdist(1) = 1                                 ! First element is 1
    topo%vtxdist(isize + 1) = topo%global_num_cells + 1 ! Last element is total number of cells + 1

    ! Divide the total number of cells by the world size to
    ! compute the chunk sizes
    k = int(real(topo%global_num_cells) / isize)
    j = 1

    do i = 1, isize
      topo%vtxdist(i) = j
      j = j + k
    enddo

    start_index = topo%vtxdist(irank + 1)
    end_index = topo%vtxdist(irank + 2) - 1
  
    ! Allocate adjacency array xadj based on vtxdist
    allocate(topo%xadj(topo%vtxdist(irank + 2) - topo%vtxdist(irank + 1) + 1)) 

    ! Allocated temporary 2D integer work array and initialise to 0
    allocate(wkint2d(topo%vtxdist(irank + 2) - topo%vtxdist(irank + 1), topo%max_faces + 1))
    wkint2d = 0

    ! Allocate global boundaries array
    allocate(topo%global_boundaries(topo%global_num_cells))

    ! All ranks loop over all the faces
    do i=1,topo%global_num_faces

      face_nb1 = topo%face_edge_end1(i)
      face_nb2 = topo%face_edge_end2(i)

      ! If face neighbour 1 is local to the current rank
      ! and face neighbour 2 is not 0
      if(face_nb1 .ge. start_index .and. face_nb1 .le. end_index .and. face_nb2 .ne. 0) then
         local_index = face_nb1 - start_index + 1                 ! Local cell index
         k = wkint2d(local_index, topo%max_faces + 1) + 1  ! Increment number of faces for this cell
         wkint2d(local_index, k) = face_nb2                       ! Store global index of neighbour cell
         wkint2d(local_index, topo%max_faces + 1) = k      ! Store number of faces for this cell
      endif

      ! If face neighbour 2 is local to the current rank
      ! and face neighbour 1 is not 0
      if(face_nb2 .ge. start_index .and. face_nb2 .le. end_index .and. face_nb1  .ne. 0) then
         local_index = face_nb2 - start_index + 1                 ! Local cell index
         k = wkint2d(local_index, topo%max_faces + 1) + 1  ! Increment number of faces for this cell
         wkint2d(local_index, k) = face_nb1                       ! Store global index of neighbour cell
         wkint2d(local_index, topo%max_faces + 1) = k      ! Store number of faces for this cell
      endif

      ! If face neighbour 2 is 0 we have a boundary face
      if(face_nb2 .eq. 0) then
        topo%global_boundaries(face_nb1) = topo%global_boundaries(face_nb1) + 1
      end if

    enddo

    num_connections = sum(wkint2d(:, topo%max_faces+1))
    print*, "num_connections on rank",irank,":", num_connections

    ! Allocate adjncy array based on the number of computed connections
    allocate(topo%adjncy(num_connections))
    ! Allocate partition array
    allocate(topo%part(topo%vtxdist(irank+2)-topo%vtxdist(irank+1)))

    local_index = 1

    do i=1,topo%vtxdist(irank+2)-topo%vtxdist(irank+1)  ! Loop over local cells
      
      topo%xadj(i) = local_index                          ! Pointer to start of current
       
      do j=1,wkint2d(i, topo%max_faces+1)               ! Loop over number of faces
        topo%adjncy(local_index + j - 1) = wkint2d(i,j) ! Store global IDs of neighbour cells
      end do

       local_index = local_index + wkint2d(i,topo%max_faces+1)
       topo%xadj(i+1) = local_index

    end do

    ! Allocate weight arrays
    allocate(topo%adjwgt(num_connections))
    ! Allocate partition array
    allocate(topo%vwgt(topo%vtxdist(irank+2)-topo%vtxdist(irank+1)))    

    deallocate(wkint2d)

  end subroutine compute_partitioner_input

end submodule

! allocate(part(vtxdist(pid+1)-vtxdist(pid)),stat=istat)
! idxlocal=1
! do i=1,vtxdist(pid+1)-vtxdist(pid)       ! Loop over local cells
!    xadj(i)=idxlocal                      ! Pointer to start of current
!    do j=1,wkint2d(i,icf+1)               ! Loop over number of faces
!       adjncy(idxlocal+j-1)=wkint2d(i,j)  ! Store global IDs of neighbour cells
!    enddo
!    idxlocal=idxlocal+wkint2d(i,icf+1)
!    xadj(i+1)=idxlocal
! enddo