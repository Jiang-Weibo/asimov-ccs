!> @brief Program file PETSc ex3
!>
!> @details Port of PETSc ksp/tutorial/ex3.c to ASiMoV-CCS style code - this is to help
!>          determine how to interface our design with PETSc.

program ex3

  !! External uses
  use MPI
  use petsc, only : PETSC_COMM_WORLD

  !! ASiMoV-CCS uses
  use accs_kinds, only : accs_real, accs_int, accs_err
  use accs_types, only : vector_init_data, vector, matrix_init_data, matrix
  use accsvec, only : create_vector, axpy, norm
  use accsmat, only : create_matrix
  use accs_utils, only : accs_init, accs_finalise, update

  use accs_petsctypes, only : matrix_petsc
  
  implicit none

  class(vector), allocatable :: u, b, ustar
  type(vector_init_data) :: vec_sizes
  class(matrix), allocatable :: M
  type(matrix_init_data) :: mat_sizes = matrix_init_data(rglob=-1, cglob=-1, rloc=-1, cloc=-1, comm=MPI_COMM_NULL)

  integer(accs_int), parameter :: eps = 100 ! Elements per side
                                            ! XXX: temporary parameter - this should be read from input
  
  integer(accs_int) :: istart, iend
  real(accs_real) :: h
  integer(accs_int), parameter :: npe = 4 ! Points per element

  real(accs_real) :: err_norm

  real(accs_real), dimension(16) :: K ! Element stiffness matrix
  
  integer(accs_err) :: ierr
  integer :: nrank
  
  !! Initialise program
  call accs_init()
  call MPI_Comm_rank(PETSC_COMM_WORLD, nrank, ierr)
  istart = 0; iend = 0
  
  !! Create stiffness matrix
  mat_sizes%rglob = (eps+1)**2
  mat_sizes%cglob = mat_sizes%rglob
  mat_sizes%comm = PETSC_COMM_WORLD
  call create_matrix(mat_sizes, M)

  call form_element_stiffness(h**2, K)
  call assemble_matrix(K, M)
  call update(M) ! Performs the parallel assembly

  !! Create right-hand-side and solution vectors
  vec_sizes%nloc = -1
  vec_sizes%nglob = (eps+1)**2
  vec_sizes%comm = PETSC_COMM_WORLD
  call create_vector(vec_sizes, u)
  call create_vector(vec_sizes, b)
  call update(u) ! Performs the parallel assembly
  
  !! Evaluate right-hand-side vector
  call eval_rhs(b)
  call update(b) ! Performs the parallel assembly
  
  !! Modify matrix and right-hand-side vector to apply Dirichlet boundary conditions
  !! Create linear solver & set options
  !! Solve linear system
  !! Check solution
  call create_vector(vec_sizes, ustar)
  call update(ustar) ! Performs the parallel assembly
  call axpy(-1.0_accs_real, u, ustar)
  err_norm = norm(ustar, 2)
  if (nrank == 0) then
     print *, "Norm of error = ", err_norm
  end if
  
  !! Clean up
  deallocate(u)
  deallocate(b)
  deallocate(ustar)
  deallocate(M)
  
  call accs_finalise()

contains

  subroutine eval_rhs(b)

    use accs_constants, only : add_mode
    use accs_types, only : vector_values
    use accs_utils, only : set_values
    
    class(vector), intent(inout) :: b

    integer(accs_int) :: i
    real(accs_real) :: x, y

    type(vector_values) :: val_dat

    val_dat%mode = add_mode
    allocate(val_dat%idx(npe))
    allocate(val_dat%val(npe))
    
    do i = istart, iend
       x = h * modulo(i, eps)
       y = h * (i / eps)

       call element_indices(i, val_dat%idx)
       call eval_element_rhs(x, y, h**2, val_dat%val)
       call set_values(val_dat, b)
    end do

    deallocate(val_dat%idx)
    deallocate(val_dat%val)
    
  end subroutine eval_rhs

  pure subroutine element_indices (i, idx)

    integer(accs_int), intent(in) :: i
    integer(accs_int), dimension(npe), intent(out) :: idx

    idx(1) = (eps + 1) * (i / eps) + modulo(i, eps)
    idx(2) = idx(1) + 1
    idx(3) = idx(2) + (eps + 1)
    idx(4) = idx(3) - 1
    
  end subroutine element_indices

  !> @brief Apply forcing function
  pure subroutine eval_element_rhs (x, y, H, r)
    
    real(accs_real), intent(in) :: x, y, H
    real(accs_real), dimension(npe), intent(out) :: r
    
    r(:) = 0.0 &
         + 0 * (x + y + H) ! Silence unused dummy argument error
    
  end subroutine eval_element_rhs

  pure subroutine form_element_stiffness(H, K)

    real(accs_real), intent(in) :: H
    real(accs_real), dimension(16), intent(inout) :: K

    K(1)  = H/6.0;   K(2)  = -.125*H;  K(3)  = H/12.0;  K(4)  = -.125*H;
    K(5)  = -.125*H; K(6)  = H/6.0;    K(7)  = -.125*H; K(8)  = H/12.0;
    K(9)  = H/12.0;  K(10)  = -.125*H; K(11) = H/6.0;   K(12) = -.125*H;
    K(13) = -.125*H; K(14) = H/12.0;   K(15) = -.125*H; K(16) = H/6.0;
    
  end subroutine form_element_stiffness

  subroutine assemble_matrix(K, M)

    use accs_constants, only : add_mode
    use accs_types, only : matrix_values
    use accs_utils, only : set_values
    
    real(accs_real), dimension(16), intent(in) :: K
    class(matrix), intent(inout) :: M

    type(matrix_values) :: mat_coeffs
    integer(accs_int) :: i

    allocate(mat_coeffs%rglob(4))
    allocate(mat_coeffs%cglob(4))
    allocate(mat_coeffs%val(16))
    mat_coeffs%val(:) = K(:)
    mat_coeffs%mode = add_mode
    
    do i = istart, iend
       call element_indices(i, mat_coeffs%rglob)
       mat_coeffs%cglob = mat_coeffs%rglob
       call set_values(mat_coeffs, M)
    end do

    deallocate(mat_coeffs%rglob)
    deallocate(mat_coeffs%cglob)
    deallocate(mat_coeffs%val)
  end subroutine assemble_matrix
end program ex3
