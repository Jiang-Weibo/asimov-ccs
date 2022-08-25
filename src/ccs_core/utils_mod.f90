!v Module file utils.mod
!
!  Provides utility functions for ASiMoV-CCS, these should be polymorphic on their input
!  and call type-specific implementations of the interface in other modules.

module utils

  use iso_c_binding

  use vec, only: set_vector_values, update_vector, begin_update_vector, end_update_vector, &
                 initialise_vector, set_vector_size, &
                 set_vector_values_mode, set_vector_values_row, set_vector_values_entry, &
                 clear_vector_values_entries, &
                 mult_vec_vec, scale_vec, zero_vector
  use mat, only: set_matrix_values, update_matrix, begin_update_matrix, end_update_matrix, &
                 initialise_matrix, finalise_matrix, set_matrix_size, &
                 set_matrix_values_mode, set_matrix_values_row, set_matrix_values_col, set_matrix_values_entry, &
                 clear_matrix_values_entries, zero_matrix
  use solver, only: initialise_equation_system
  use kinds, only: ccs_int, ccs_real

  implicit none

  private

  public :: set_values
  public :: set_entry
  public :: clear_entries
  public :: begin_update
  public :: end_update
  public :: update
  public :: finalise
  public :: initialise
  public :: set_size
  public :: mult
  public :: zero
  public :: set_mode
  public :: set_row
  public :: set_col
  public :: str
  public :: debug_print
  public :: exit_print

  !> Generic interface to set values on an object.
  interface set_values
    module procedure set_vector_values
    module procedure set_matrix_values
  end interface set_values

  !> Generic interface to store a value for later setting.
  interface set_entry
    module procedure set_vector_values_entry
    module procedure set_matrix_values_entry
  end interface set_entry

  !> Generic interface to set the storage mode (ADD/INSERT) of an object.
  interface set_mode
    module procedure set_vector_values_mode
    module procedure set_matrix_values_mode
  end interface set_mode

  !> Generic interface to set the working row.
  interface set_row
    module procedure set_vector_values_row
    module procedure set_matrix_values_row
  end interface set_row
  !> Generic interface to set the working column.
  interface set_col
    module procedure set_matrix_values_col
  end interface set_col

  !> Generic interface to indicate an object is ready for use.
  interface finalise
    module procedure finalise_matrix
  end interface finalise

  !> Generic interface to clear stored entries.
  interface clear_entries
    module procedure clear_vector_values_entries
    module procedure clear_matrix_values_entries
  end interface clear_entries

  !> Generic interface to perform parallel update of an object.
  interface update
    module procedure update_vector
    module procedure update_matrix
  end interface update

  !v Generic interface to begin parallel update of an object.
  !
  !  This is to allow overlapping comms and compute.
  interface begin_update
    module procedure begin_update_vector
    module procedure begin_update_matrix
  end interface begin_update

  !v Generic interface to end parallel update of an object.
  !
  !  This is to allow overlapping comms and compute.
  interface end_update
    module procedure end_update_vector
    module procedure end_update_matrix
  end interface end_update

  !> Generic interface to initialse vectors, matrices and linear systems
  interface initialise
    module procedure initialise_vector
    module procedure initialise_matrix
    module procedure initialise_equation_system
  end interface initialise

  !> Generic interface to set vector and matrix sizes
  interface set_size
    module procedure set_vector_size
    module procedure set_matrix_size
  end interface set_size

  !>  Generic interface to perform multiplications
  interface mult
    module procedure mult_vec_vec
  end interface mult

  !> Generic interface to zero an object
  interface zero
    module procedure zero_vector
    module procedure zero_matrix
  end interface zero

  !> Generic interface to converting numbers to strings
  interface str
    module procedure int2str
    module procedure real2str
  end interface str

  !> Generic interface to debug printer
  interface debug_print
    module procedure debug_print_actual
    module procedure noop
  end interface debug_print

contains

  !> Print a message, along with with its location.
  subroutine debug_print_actual(msg, filepath, line)
    use mpi

    character(*), intent(in) :: msg      !< text to be printed
    character(*), intent(in) :: filepath !< calling file, added by __FILE__ macro
    integer, intent(in) :: line          !< line number in calling file, added by __LINE__ macro

    character(len=:), allocatable :: filename
    integer :: slash_position
    integer :: rank, ierror
    logical :: init_flag

    slash_position = scan(trim(filepath), "/", back=.true.)
    if (slash_position .gt. 0) then
      filename = filepath(slash_position + 1:len_trim(filepath))
    else
      filename = filepath
    end if

    call mpi_initialized(init_flag, ierror)
    if (init_flag) then
      call mpi_comm_rank(MPI_COMM_WORLD, rank, ierror)
      print *, trim(filename), "(", int2str(line), ")[", int2str(rank), "] : ", msg
    else
      print *, trim(filename), "(", int2str(line), ") : ", msg
    end if
  end subroutine

  !> No-op routine, does nothing.
  subroutine noop()
  end subroutine

  !> Convert integer to string.
  function int2str(in_int, format_str) result(out_string)
    integer(ccs_int), intent(in) :: in_int           !< integer to convert
    character(*), optional, intent(in) :: format_str !< format string to use
    character(:), allocatable :: out_string          !< formatted string from input integer

    character(32) :: tmp_string

    if (present(format_str)) then
      write (tmp_string, format_str) in_int
    else
      write (tmp_string, *) in_int
    end if
    out_string = trim(adjustl(tmp_string))
  end function

  !> Convert real to string.
  function real2str(in_real, format_str) result(out_string)
    real(ccs_real), intent(in) :: in_real            !< real number to convert
    character(*), optional, intent(in) :: format_str !< format string to use
    character(:), allocatable :: out_string          !< formatted string from input real

    character(32) :: tmp_string

    if (present(format_str)) then
      write (tmp_string, format_str) in_real
    else
      write (tmp_string, *) in_real
    end if
    out_string = trim(adjustl(tmp_string))
  end function

  !> Print a message and stop program execution.
  subroutine exit_print(msg, filepath, line)
    character(*), intent(in) :: msg      !< text to be printed
    character(*), intent(in) :: filepath !< calling file, added by __FILE__ macro
    integer, intent(in) :: line          !< line number in calling file, added by __LINE__ macro

    call debug_print(msg, filepath, line)
    stop 1
  end subroutine exit_print

end module utils
