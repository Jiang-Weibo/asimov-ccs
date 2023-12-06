!v Module file utils.mod
!
!  Provides utility functions for ASiMoV-CCS, these should be polymorphic on their input
!  and call type-specific implementations of the interface in other modules.

module utils
#include "ccs_macros.inc"

  use iso_c_binding

  use vec, only: set_vector_values, update_vector, begin_update_vector, end_update_vector, &
                 initialise_vector, set_vector_size, &
                 set_vector_values_mode, set_vector_values_row, set_vector_values_entry, &
                 clear_vector_values_entries, &
                 mult_vec_vec, scale_vec, zero_vector, &
                 get_natural_data_vec
  use mat, only: set_matrix_values, update_matrix, begin_update_matrix, end_update_matrix, &
                 initialise_matrix, finalise_matrix, set_matrix_size, &
                 set_matrix_values_mode, set_matrix_values_row, set_matrix_values_col, set_matrix_values_entry, &
                 clear_matrix_values_entries, zero_matrix
  use solver, only: initialise_equation_system
  use kinds, only: ccs_int, ccs_real
  use types, only: field, fluid, fluid_solver_selector, field_ptr
  use constants, only: field_u, field_v, field_w, field_p, field_p_prime, field_mf, field_viscosity, &
                       cell_centred_central, cell_centred_upwind

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
  public :: calc_kinetic_energy
  public :: calc_enstrophy
  public :: add_field_to_outputlist
  public :: reset_outputlist_counter
  public :: get_field
  public :: get_fluid_solver_selector
  public :: set_field
  public :: set_fluid_solver_selector
  public :: allocate_fluid_fields
  public :: dealloc_fluid_fields
  public :: get_natural_data
  public :: get_scheme_name
  public :: get_scheme_id

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

  !> Generic interface to perform multiplications
  interface mult
    module procedure mult_vec_vec
  end interface mult

  !> Generic interface to zero an object
  interface zero
    module procedure zero_vector
    module procedure zero_matrix
  end interface zero

  !> Generic interface to converting numbers and bools to strings
  interface str
    module procedure int2str
    module procedure real2str
    module procedure bool2str
  end interface str

  !> Generic interface to debug printer
  interface debug_print
    module procedure debug_print_actual
    module procedure noop
  end interface debug_print

  !> Generic interface to get data in natural ordering
  interface get_natural_data
    module procedure get_natural_data_vec
  end interface get_natural_data

  integer(ccs_int), save :: outputlist_counter = 0

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

  !> Convert bool to string.
  function bool2str(in_bool) result(out_string)
    logical, intent(in) :: in_bool          !< bool to convert
    character(:), allocatable :: out_string !< string from input bool

    character(32) :: tmp_string

    write (tmp_string, *) in_bool

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

  !TODO: move this subroutine to more appropriate module
  !> Calculate kinetic energy over density
  subroutine calc_kinetic_energy(par_env, u, v, w)

    use mpi

    use constants, only: ndim, ccs_string_len
    use types, only: field, ccs_mesh, cell_locator
    use vec, only: get_vector_data, restore_vector_data
    use parallel, only: allreduce, error_handling
    use parallel_types_mpi, only: parallel_environment_mpi
    use parallel_types, only: parallel_environment
    use meshing, only: get_local_num_cells, create_cell_locator, get_volume
    use timestepping, only: get_current_step

    class(parallel_environment), allocatable, intent(in) :: par_env !< parallel environment
    class(field), intent(inout) :: u !< solve x velocity field
    class(field), intent(inout) :: v !< solve y velocity field
    class(field), intent(inout) :: w !< solve z velocity field

    integer(ccs_int) :: local_num_cells
    real(ccs_real) :: ek_local, ek_global, volume_local, volume_global
    real(ccs_real), dimension(:), pointer :: u_data, v_data, w_data
    real(ccs_real) :: rho
    integer(ccs_int) :: index_p
    character(len=ccs_string_len) :: fmt
    integer(ccs_int) :: ierr
    integer(ccs_int) :: step !< timestep

    logical, save :: first_time = .true.
    integer :: io_unit

    type(cell_locator) :: loc_p
    real(ccs_real) :: volume

    call get_current_step(step)
    rho = 1.0_ccs_real

    ek_local = 0.0_ccs_real
    ek_global = 0.0_ccs_real
    volume_local = 0.0_ccs_real
    volume_global = 0.0_ccs_real

    call get_vector_data(u%values, u_data)
    call get_vector_data(v%values, v_data)
    call get_vector_data(w%values, w_data)

    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells
      call create_cell_locator(index_p, loc_p)
      call get_volume(loc_p, volume)

      ek_local = ek_local + 0.5 * rho * volume * &
                 (u_data(index_p)**2 + v_data(index_p)**2 + w_data(index_p)**2)

      volume_local = volume_local + volume
    end do

    call restore_vector_data(u%values, u_data)
    call restore_vector_data(v%values, v_data)
    call restore_vector_data(w%values, w_data)

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_AllReduce(ek_local, ek_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, par_env%comm, ierr)
      call error_handling(ierr, "mpi", par_env)
      call MPI_AllReduce(volume_local, volume_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, par_env%comm, ierr)
      call error_handling(ierr, "mpi", par_env)
    class default
      call error_abort("ERROR: Unknown type")
    end select

    ek_global = ek_global / volume_global

    if (par_env%proc_id == par_env%root) then
      if (first_time) then
        first_time = .false.
        open (newunit=io_unit, file="tgv2d-ek.log", status="replace", form="formatted")
      else
        open (newunit=io_unit, file="tgv2d-ek.log", status="old", form="formatted", position="append")
      end if
      fmt = '(I0,1(1x,e12.4))'
      write (io_unit, fmt) step, ek_global
      close (io_unit)
    end if

  end subroutine calc_kinetic_energy

  !TODO: move this subroutine to more appropriate module
  !> Calculate enstrophy
  subroutine calc_enstrophy(par_env, u, v, w)

    use mpi

    use constants, only: ndim, ccs_string_len
    use types, only: field, ccs_mesh
    use vec, only: get_vector_data, restore_vector_data
    use parallel, only: allreduce, error_handling
    use parallel_types_mpi, only: parallel_environment_mpi
    use parallel_types, only: parallel_environment
    use meshing, only: get_local_num_cells
    use timestepping, only: get_current_step

    class(parallel_environment), allocatable, intent(in) :: par_env !< parallel environment
    class(field), intent(inout) :: u !< solve x velocity field
    class(field), intent(inout) :: v !< solve y velocity field
    class(field), intent(inout) :: w !< solve z velocity field

    integer(ccs_int) :: local_num_cells
    real(ccs_real) :: ens_local, ens_global
    real(ccs_real), dimension(:), pointer :: dudy, dudz, dvdx, dvdz, dwdx, dwdy
    integer(ccs_int) :: index_p
    character(len=ccs_string_len) :: fmt
    integer(ccs_int) :: ierr
    integer(ccs_int) :: step !< timestep

    logical, save :: first_time = .true.
    integer :: io_unit

    call get_current_step(step)
    ens_local = 0.0_ccs_real
    ens_global = 0.0_ccs_real

    call get_vector_data(u%y_gradients, dudy)
    call get_vector_data(u%z_gradients, dudz)
    call get_vector_data(v%x_gradients, dvdx)
    call get_vector_data(v%z_gradients, dvdz)
    call get_vector_data(w%x_gradients, dwdx)
    call get_vector_data(w%y_gradients, dwdy)

    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells

      ens_local = ens_local + (dwdy(index_p) - dvdz(index_p))**2 + &
                  (dudz(index_p) - dwdx(index_p))**2 + &
                  (dvdx(index_p) - dudy(index_p))**2

    end do
    ens_local = 0.5 * ens_local

    call restore_vector_data(u%y_gradients, dudy)
    call restore_vector_data(u%z_gradients, dudz)
    call restore_vector_data(v%x_gradients, dvdx)
    call restore_vector_data(v%z_gradients, dvdz)
    call restore_vector_data(w%x_gradients, dwdx)
    call restore_vector_data(w%y_gradients, dwdy)

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_AllReduce(ens_local, ens_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, par_env%comm, ierr)
      call error_handling(ierr, "mpi", par_env)
    class default
      call error_abort("ERROR: Unknown type")
    end select

    if (par_env%proc_id == par_env%root) then
      if (first_time) then
        first_time = .false.
        open (newunit=io_unit, file="tgv2d-ens.log", status="replace", form="formatted")
      else
        open (newunit=io_unit, file="tgv2d-ens.log", status="old", form="formatted", position="append")
      end if
      fmt = '(I0,1(1x,e12.4))'
      write (io_unit, fmt) step, ens_global
      close (io_unit)
    end if

  end subroutine calc_enstrophy

  !v Append new field to output list.
  !
  ! @note does not check for duplicates
  subroutine add_field_to_outputlist(var, name, list)

    use types, only: field, field_ptr

    ! Arguments
    class(field), pointer, intent(in) :: var !< The field to be added
    character(len=*), intent(in) :: name     !< The name of the field
    type(field_ptr), dimension(:), allocatable, intent(inout) :: list !< The output list

    ! Local variables
    type(field_ptr) :: new_element

    new_element%ptr => var
    new_element%name = name
    
    if (allocated(list)) then
      if (size(list) > outputlist_counter) then
        list(outputlist_counter) = new_element
      else
        list = [list, new_element]
      end if
    else
      list = [new_element]
    end if
    
    outputlist_counter = outputlist_counter + 1
    
  end subroutine add_field_to_outputlist

  subroutine reset_outputlist_counter()

    outputlist_counter = 0

  end subroutine

  !> Gets the field from the fluid structure specified by field_name
  subroutine get_field(flow, field_name, flow_field)
    !(flow, "u", u)
    type(fluid), intent(in) :: flow                   !< the structure containing all the fluid fields
    character(len=*), intent(in) :: field_name
    class(field), pointer, intent(out) :: flow_field  !< the field of interest

    integer(ccs_int), dimension(1) :: field_index
    integer(ccs_int) :: i
    character(len=:), allocatable :: msg                   !< Constructed message

    field_index(1) = 0
    do i = 1, size(flow%field_names)
      if (trim(flow%field_names(i)) == field_name) then
        field_index(1) = i
        exit
      end if           
    end do

    if (field_index(1) == 0) then
      msg = "Field " // field_name // " not found"
      call error_abort(msg)
    end if

    !field_index = findloc(flow%field_names , field_name)
    flow_field => flow%fields(field_index(1))%ptr
  end subroutine get_field

  !< Sets the pointer to the field and the corresponding field name in the fluid structure
  !subroutine set_field(field_index, flow_field, flow)
  subroutine set_field(flow_field, flow)
    !integer(ccs_int), intent(in) :: field_index     !< index of arrays at which to set the field pointer and name
    integer(ccs_int) :: field_index     !< index of arrays at which to set the field pointer and name
    class(field), target, intent(in) :: flow_field  !< the field
    type(fluid), intent(inout) :: flow              !< the fluid structure
    type(field_ptr) :: tmp_field_ptr
    logical, save :: first_call = .true.
    
    if (first_call) then
      print*, "inside first call"
      allocate (flow%fields(1))
      print*, "a"
      allocate (flow%field_names(1))
      print*, "b"
      flow%fields(1)%ptr => flow_field
      print*, "c"
      flow%field_names(1) = flow_field%name
      print*, "d"
      first_call = .false.
      print*, "e"

    else
      print*, "inside 2nd call"
      tmp_field_ptr%ptr => flow_field
      flow%fields = [ flow%fields, tmp_field_ptr]
      flow%field_names = [flow%field_names, flow_field%name]
      field_index = size(flow%fields) !creating issue
      print*, "field index=",field_index

      flow%fields(field_index)%ptr => flow_field
    end if 
    
  end subroutine set_field

  !> Gets the solver selector for a specified field
  subroutine get_fluid_solver_selector(solver_selector, field_name, selector)
    type(fluid_solver_selector), intent(in) :: solver_selector  !< Structure containing all of the solver selectors
    integer(ccs_int), intent(in) :: field_name                  !< name of field
    logical, intent(out) :: selector                            !< flag indicating whether to solve for the given field

    select case (field_name)
    case (field_u)
      selector = solver_selector%u
    case (field_v)
      selector = solver_selector%v
    case (field_w)
      selector = solver_selector%w
    case (field_p)
      selector = solver_selector%p
    case default
      call error_abort("Unrecognised field index.")
    end select
  end subroutine get_fluid_solver_selector

  !> Sets the solver selector for a specified field
  subroutine set_fluid_solver_selector(field_name, selector, solver_selector)
    integer(ccs_int), intent(in) :: field_name                      !< name of field
    logical, intent(in) :: selector                                 !< flag indicating whether to solve for the given field
    type(fluid_solver_selector), intent(inout) :: solver_selector   !< Structure containing all of the solver selectors

    select case (field_name)
    case (field_u)
      solver_selector%u = selector
    case (field_v)
      solver_selector%v = selector
    case (field_w)
      solver_selector%w = selector
    case (field_p)
      solver_selector%p = selector
    case default
      call error_abort("Unrecognised field index.")
    end select
  end subroutine set_fluid_solver_selector

  ! Allocates arrays in fluid field structure to specified size
  subroutine allocate_fluid_fields(n_fields, flow)
    integer(ccs_int), intent(in) :: n_fields  !< Size of arrays in fluid structure
    type(fluid), intent(out) :: flow          !< the fluid structure

    allocate (flow%fields(n_fields))
    allocate (flow%field_names(n_fields))
  end subroutine allocate_fluid_fields

  ! Deallocates fluid arrays
  subroutine dealloc_fluid_fields(flow)
    type(fluid), intent(inout) :: flow  !< The fluid structure to deallocate

    deallocate (flow%fields)
    deallocate (flow%field_names)
  end subroutine dealloc_fluid_fields

  !> Convert advection scheme name -> ID.
  integer(ccs_int) function get_scheme_id(scheme_name)

    character(len=*), intent(in) :: scheme_name
    integer(ccs_int) :: id
    character(len=:), allocatable :: scheme

    scheme = trim(scheme_name)
    if (scheme == "central") then
       id = cell_centred_central
    else if (scheme == "upwind") then
       id = cell_centred_upwind
    else
       call error_abort("Uknown scheme "//scheme)
    end if

    get_scheme_id = id
  end function get_scheme_id

  !> Convert advection scheme ID -> name.
  function get_scheme_name(scheme_id) result(scheme_name)

    integer(ccs_int), intent(in) :: scheme_id
    character(len=:), allocatable :: scheme_name

    if (scheme_id == cell_centred_central) then
       scheme_name = "central"
    else if (scheme_id == cell_centred_upwind) then
       scheme_name = "upwind"
    else
       call error_abort("Uknown scheme ID")
    end if

  end function get_scheme_name
    
end module utils
