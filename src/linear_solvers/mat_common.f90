submodule (mat) mat_common

  implicit none

contains

  !> @brief Constructor for default matrix values
  !
  !> param[in/out] matrix_descriptor - the initialised matrix values
  module subroutine initialise_matrix(matrix_descriptor)
    type(matrix_init_data), intent(inout) :: matrix_descriptor
    matrix_descriptor%par_env => null()
    matrix_descriptor%par_env => null()
  end subroutine initialise_matrix

  !> @brief Setter for global matrix size
  !
  !> param[in] par_env                - the parallel environment where 
  !!                                    the matrix resides
  !> param[in] geometry               - the mesh object
  !> param[in/out] matrix_descriptor  - the matrix data object
  module subroutine set_matrix_size(par_env, geometry, matrix_descriptor)
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    class(ccs_mesh), target, intent(in) :: geometry
    type(matrix_init_data), intent(inout) :: matrix_descriptor

    matrix_descriptor%mesh => geometry
    matrix_descriptor%par_env => par_env
  end subroutine set_matrix_size

  !> @brief Setter for matrix number of non-zeros
  !
  !> param[in] nnz                   - the number of non-zeros
  !> param[in/out] matrix_descriptor - the matrix data object
  module subroutine set_nnz(nnz, matrix_descriptor)
    integer(ccs_int), intent(in) :: nnz
    type(matrix_init_data), intent(inout) :: matrix_descriptor

    matrix_descriptor%nnz = nnz
  end subroutine  set_nnz

  module procedure create_matrix_values
    allocate(val_dat%rglob(nrows))
    allocate(val_dat%cglob(nrows))
    allocate(val_dat%val(nrows))
  end procedure create_matrix_values

  module procedure set_matrix_values_mode
    val_dat%mode = mode
  end procedure set_matrix_values_mode
  
  module subroutine set_matrix_values_entry(val, val_dat)

    use constants, only : add_mode, insert_mode

    real(ccs_real), intent(in) :: val
    type(matrix_values), intent(inout) :: val_dat

    associate(x => val_dat%val(val_dat%current_entry), &
         mode => val_dat%mode)
      if (mode == insert_mode) then
        x = val
      else if (mode == add_mode) then
        x = x + val
      else
        print *, "ERROR: Unrecognised entry mode ", mode
        stop
      end if

    end associate
    
  end subroutine set_matrix_values_entry
  
end submodule 
