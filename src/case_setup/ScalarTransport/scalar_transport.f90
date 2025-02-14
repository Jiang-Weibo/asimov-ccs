!v Program file for ScalarTransport case
!
!  @build mpi+petsc

program scalar_transport
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use ccs_base, only: mesh
  use case_config, only: num_steps, num_iters, cps, domain_size, case_name, &
                         res_target, write_gradients, dt, write_frequency, restart, unsteady
  use constants, only: cell, face, ccsconfig, ccs_string_len, &
                       face_centred, cell_centred_central, cell_centred_upwind, &
                       ccs_split_type_low_high
  use kinds, only: ccs_real, ccs_int
  use types, only: field, field_spec, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector, field_ptr, fluid
  use fields, only: create_field, set_field_config_file, set_field_n_boundaries, set_field_name, &
                    set_field_type, set_field_vector_properties, set_field_store_residuals
  use fortran_yaml_c_interface, only: parse
  use parallel, only: initialise_parallel_environment, create_new_par_env, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, &
                      is_root
  use parallel_types, only: parallel_environment
  use mesh_utils, only: build_mesh, write_mesh
  use meshing, only: get_global_num_cells, get_centre, count_neighbours, &
                     create_cell_locator, create_face_locator, create_neighbour_locator, &
                     get_local_index, get_boundary_status, get_face_normal, set_mesh_object, nullify_mesh_object, &
                     get_local_num_cells
  use vec, only: create_vector, set_vector_location, get_vector_data, restore_vector_data
  use scalars, only: update_scalars
  use utils, only: set_size, initialise, update, exit_print, add_field_to_outputlist, &
                   get_field, &
                   allocate_fluid_fields, dealloc_fluid_fields, &
                   get_scheme_name
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use read_config, only: get_boundary_count, get_boundary_names, get_store_residuals, get_variables, get_variable_types
  use io_visualisation, only: write_solution, read_solution
  use timestepping, only: set_timestep, activate_timestepping, initialise_old_values, finalise_timestep

  implicit none

  class(parallel_environment), allocatable :: par_env
  class(parallel_environment), allocatable, target :: shared_env
  character(len=:), allocatable :: input_path  ! Path to input directory
  character(len=:), allocatable :: case_path  ! Path to input directory with case name appended
  character(len=:), allocatable :: ccs_config_file ! Config file for CCS
  character(len=ccs_string_len), dimension(:), allocatable :: variable_names ! variable names for BC reading
  integer(ccs_int), dimension(:), allocatable :: variable_types              ! cell centred upwind, central, etc.

  type(vector_spec) :: vec_properties

  type(field_spec) :: field_properties
  class(field), pointer :: phi ! Generic field to work with
  
  integer(ccs_int) :: n_boundaries

  integer(ccs_int) :: it_start, it_end
  integer(ccs_int) :: irank ! MPI rank ID
  integer(ccs_int) :: isize ! Size of MPI world

  double precision :: start_time
  double precision :: init_time
  double precision :: end_time

  logical :: store_residuals
  logical :: use_mpi_splitting

  type(fluid) :: flow_fields

  real(ccs_real) :: L
  integer(ccs_int) :: t

  integer :: i

  character(len=128), dimension(:), allocatable :: bnd_names
  
  ! Launch MPI
  call initialise_parallel_environment(par_env)
  use_mpi_splitting = .false.
  call create_new_par_env(par_env, ccs_split_type_low_high, use_mpi_splitting, shared_env)

  irank = par_env%proc_id
  isize = par_env%num_procs

  call read_command_line_arguments(par_env, cps, case_name=case_name, in_dir=input_path)

  if (allocated(input_path)) then
    case_path = input_path // "/" // case_name
  else
    case_path = case_name
  end if

  ccs_config_file = case_path // ccsconfig

  call timer(start_time)

  ! Read case name from configuration file
  call read_configuration(ccs_config_file)

  if (irank == par_env%root) print *, "Starting ", case_name, " case!"

  ! Set start and end iteration numbers (read from input file)
  it_start = 1
  it_end = num_iters

  ! Create a mesh
  if (irank == par_env%root) print *, "Building mesh"
  L = 1.0_ccs_real
  call get_boundary_names(ccs_config_file, bnd_names)
  mesh = build_mesh(par_env, shared_env, cps, cps, cps, L, bnd_names)   ! 3-D mesh
  call set_mesh_object(mesh)

  ! Initialise fields
  if (irank == par_env%root) print *, "Initialise fields"

  ! Write gradients to solution file
  write_gradients = .false.

  ! Create and initialise field vectors
  call get_boundary_count(ccs_config_file, n_boundaries)
  call get_store_residuals(ccs_config_file, store_residuals)

  call initialise(vec_properties)

  call set_vector_location(cell, vec_properties)
  call set_size(par_env, mesh, vec_properties)

  call set_field_config_file(ccs_config_file, field_properties)
  call set_field_n_boundaries(n_boundaries, field_properties)
  call set_field_store_residuals(store_residuals, field_properties)

  ! Build viscosity and density explicitly (for time being)
  call set_field_vector_properties(vec_properties, field_properties)
  call set_field_type(cell_centred_central, field_properties)
  call set_field_name("viscosity", field_properties)
  call create_field(par_env, field_properties, flow_fields) 
  call set_field_type(cell_centred_central, field_properties)
  call set_field_name("density", field_properties)
  call create_field(par_env, field_properties, flow_fields) 

  if (is_root(par_env)) then
    print *, "Build field list"
  end if
  do i = 1, size(variable_names)
    if (is_root(par_env)) then
      print *, "Creating field ", trim(variable_names(i))
    end if
    call set_field_type(variable_types(i), field_properties)
    call set_field_name(variable_names(i), field_properties)
    call create_field(par_env, field_properties, flow_fields)
  end do
  if (is_root(par_env)) then
    print *, "Built ", size(variable_names), " dynamically-defined fields"
  end if
  
  !! Create mass flux as a special variable
  call set_vector_location(face, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call set_field_vector_properties(vec_properties, field_properties)
  call set_field_type(face_centred, field_properties)
  call set_field_name("mf", field_properties)
  call create_field(par_env, field_properties, flow_fields)

  do i = 1, size(flow_fields%fields)
    call get_field(flow_fields, i, phi)
    if ((phi%name == "whisky") .or. (phi%name == "water")) then
      call add_field_to_outputlist(phi)
    end if
  end do
  nullify(phi)

  ! Initialise velocity field
  if (irank == par_env%root) print *, "Initialise flow field"
  call initialise_case(flow_fields) 

  ! XXX: Opportunity for utility subroutine: update(fluid)
  do i = 1, size(flow_fields%fields)
    call get_field(flow_fields, i, phi)
    call update(phi%values)
    nullify(phi)
  end do

  if(restart) then
    if (is_root(par_env)) then
      print*, "restart capability activated"
    end if
    call read_solution(par_env, case_path, mesh, flow_fields)
  end if 

  if (irank == par_env%root) then
    call print_configuration()
  end if

  call activate_timestepping()
  call set_timestep(dt)

  call timer(init_time)

  if(.not.unsteady) then
    num_steps = 1
    print*, "steady-state activated"
  else
    print*, "unsteady-state activated"
  end if

  ! Solve using SIMPLE algorithm
  if (irank == par_env%root) print *, "Start scalar solver"
  do t = 1, num_steps
    call write_mesh(par_env, case_path, mesh)
    call update_scalars(par_env, mesh, eval_sources, flow_fields)
    if (par_env%proc_id == par_env%root) then
      print *, "TIME = ", t, " / ", num_steps
    end if

    if ((t == 1) .or. (t == num_steps) .or. (mod(t, write_frequency) == 0)) then
      if(.not. unsteady) then
        call write_solution(par_env, case_path, mesh, flow_fields)
      else
        call write_solution(par_env, case_path, mesh, flow_fields, t, num_steps, dt)
      end if 
    end if

    call finalise_timestep()
  end do

  ! Clean-up
  call dealloc_fluid_fields(flow_fields)

  call timer(end_time)

  if (irank == par_env%root) then
    print *, "Init time: ", init_time - start_time
    print *, "Elapsed time: ", end_time - start_time
  end if

  ! Finalise MPI
  call nullify_mesh_object()
  call cleanup_parallel_environment(par_env)

contains

  ! Read YAML configuration file
  subroutine read_configuration(config_filename)

    use read_config, only: get_value

    character(len=*), intent(in) :: config_filename

    class(*), pointer :: config_file  !< Pointer to CCS config file
    character(:), allocatable :: error
    
    config_file => parse(config_filename, error)
    if (allocated(error)) then
      call error_abort(trim(error))
    end if

    call get_variables(config_file, variable_names)
    if (size(variable_names) == 0) then
       call error_abort("No variables were specified.")
    end if
    call get_variable_types(config_file, variable_types)
    if (size(variable_types) /= size(variable_names)) then
       call error_abort("The number of variable types does not match the number of named variables")
    end if

    call get_value(config_file, 'restart', restart)
    
    call get_value(config_file, 'unsteady', unsteady)

    call get_value(config_file, 'iterations', num_iters) ! steady-state
    if (num_iters == huge(0)) then
      call error_abort("No value assigned to num_iters.")
    end if

    if(unsteady) then
      call get_value(config_file, 'steps', num_steps)
      if (num_steps == huge(0)) then
        call error_abort("No value assigned to num_steps.")
      end if

      call get_value(config_file, 'dt', dt)
      if (dt == huge(0.0)) then
        call error_abort("No value assigned to dt.")
      end if

      call get_value(config_file, 'write_frequency', write_frequency)
      if (write_frequency == huge(0.0)) then
        call error_abort("No value assigned to write_frequency.")
      end if
    end if 
    
    if (cps == huge(0)) then ! cps was not set on the command line
      call get_value(config_file, 'cps', cps)
      if (cps == huge(0)) then
        call error_abort("No value assigned to cps.")
      end if
    end if

    call get_value(config_file, 'L', domain_size)
    if (domain_size == huge(0.0)) then
      call error_abort("No value assigned to domain_size.")
    end if

    call get_value(config_file, 'target_residual', res_target)
    if (res_target == huge(0.0)) then
      call error_abort("No value assigned to target residual.")
    end if

  end subroutine

  ! Print test case configuration
  subroutine print_configuration()

    integer(ccs_int) :: global_num_cells

    call get_global_num_cells(global_num_cells)

    ! XXX: this should eventually be replaced by something nicely formatted that uses "write"
    print *, " "
    print *, "******************************************************************************"
    print *, "* Solving the ", case_name, " case"
    print *, "******************************************************************************"
    print *, " "
    print *, "******************************************************************************"
    print *, "* SIMULATION LENGTH"
    if (unsteady) then
      print *, "* Running for ", num_steps, "timesteps and ", num_iters, "iterations"
      write (*, '(1x, a, e10.3)') "* Time step size: ", dt
    else
      print *, "* Running for ", num_iters, "iterations"
    end if 
    print *, "******************************************************************************"
    print *, "* MESH SIZE"
    print *, "* Cells per side: ", cps
    write (*, '(1x,a,e10.3)') "* Domain size: ", domain_size
    print *, "* Global number of cells is ", global_num_cells
    print *, "******************************************************************************"
    print *, "* RELAXATION FACTORS"
    print *, "******************************************************************************"

  end subroutine

  subroutine initialise_case(flow)

    use constants, only: insert_mode
    use types, only: vector_values, cell_locator, neighbour_locator, face_locator
    use meshing, only: create_cell_locator, get_global_index, get_local_num_cells
    use fv, only: calc_cell_coords
    use utils, only: clear_entries, set_mode, set_row, set_entry, set_values
    use vec, only: get_vector_data, restore_vector_data, create_vector_values

    ! Arguments
    type(fluid), intent(inout) :: flow

    ! Local variables
    class(field), pointer :: mf, viscosity, density, phi
    integer(ccs_int) :: index_p, global_index_p, n_local
    real(ccs_real) :: whisky_val, water_val
    type(cell_locator) :: loc_p
    type(face_locator) :: loc_f
    type(neighbour_locator) :: loc_nb
    type(vector_values) :: whisky_vals, water_vals
    real(ccs_real), dimension(:), pointer :: mf_data, viscosity_data, density_data

    integer(ccs_int) :: j, nnb
    integer(ccs_int) :: index_nb, index_f
    logical :: is_boundary
    
    real(ccs_real), dimension(3) :: x, r, c, face_normal
    real(ccs_real), dimension(3) :: v
    real(ccs_real) :: theta, rmag

    call get_field(flow, "mf", mf)
    
    ! Set alias
    call get_local_num_cells(n_local)

    call create_vector_values(n_local, whisky_vals)
    call create_vector_values(n_local, water_vals)
    call set_mode(insert_mode, whisky_vals)
    call set_mode(insert_mode, water_vals)

    call get_vector_data(mf%values, mf_data)
    
    ! Coordinates of domain centre
    c = L / 2
    
    ! Set initial values for velocity fields
    do index_p = 1, n_local
      call create_cell_locator(index_p, loc_p)
      call get_global_index(loc_p, global_index_p)

      ! Get domain centre offset
      call get_centre(loc_p, x)
      r = x - c

      if (any(r <= 0.0_ccs_real)) then
        whisky_val = 1.0_ccs_real
      else
        whisky_val = 0.0_ccs_real
      end if
      water_val = 1.0_ccs_real - whisky_val

      call set_row(global_index_p, whisky_vals)
      call set_entry(whisky_val, whisky_vals)
      call set_row(global_index_p, water_vals)
      call set_entry(water_val, water_vals)

      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
         call create_face_locator(index_p, j, loc_f)
         call get_local_index(loc_f, index_f)
         call get_boundary_status(loc_f, is_boundary)

         if (.not. is_boundary) then
            call create_neighbour_locator(loc_p, j, loc_nb)
            call get_local_index(loc_nb, index_nb)

            if (index_nb > index_p) then
               call create_face_locator(index_p, j, loc_f)
               call get_face_normal(loc_f, face_normal)
               call get_centre(loc_f, x)

               ! Create a rotating field that decays to zero on the boundaries
               r = x - c
               rmag = sqrt(sum(r(1:2)**2))
               theta = asin(abs(r(2)) / rmag)
               if (r(1) >= 0) then
                  if (r(2) >= 0) then
                     theta = theta + 0 * 3.1415926 / 2
                  else
                     theta = 4 * 3.1415926 / 2 - theta
                  end if
               else
                  if (r(2) >= 0) then
                     theta = 2 * 3.1415926 / 2 - theta
                  else
                     theta = theta + 2 * 3.1415926 / 2
                  endif
               end if

               v(1) = -sin(theta)
               v(2) = cos(theta)
               v(3) = 0.0_ccs_real
               
               mf_data(index_f) = sum(v * face_normal)
            end if
         else
            mf_data(index_f) = 0.0_ccs_real
         end if

      end do

    end do

    call get_field(flow, "whisky", phi)
    call set_values(whisky_vals, phi%values)
    nullify(phi)

    call get_field(flow, "water", phi)
    call set_values(whisky_vals, phi%values)
    nullify(phi)

    deallocate (whisky_vals%global_indices)
    deallocate (water_vals%global_indices)
    deallocate (whisky_vals%values)
    deallocate (water_vals%values)

    call restore_vector_data(mf%values, mf_data)
    nullify(mf)

    call get_field(flow, "viscosity", viscosity)
    call get_vector_data(viscosity%values, viscosity_data)
    viscosity_data(:) =  1.e-2_ccs_real
    call restore_vector_data(viscosity%values, viscosity_data)
    nullify(viscosity)

    call get_field(flow, "density", density)
    call get_vector_data(density%values, density_data)
    density_data(:) = 1.0_ccs_real
    call restore_vector_data(density%values, density_data)
    nullify(density)
    
  end subroutine initialise_case

  !> Case-specific source terms
  subroutine eval_sources(flow, phi, R, S)
    use types, only: fluid, field, ccs_vector
    use fv, only: zero_sources

    type(fluid), intent(in) :: flow !< Provides access to full flow field
    class(field), intent(in) :: phi !< Field being transported
    class(ccs_vector), intent(inout) :: R !< Work vector (for evaluating linear/implicit sources)
    class(ccs_vector), intent(inout) :: S !< Work vector (for evaluating fixed/explicit sources)
    
    ! Dummy implementation - just zeros the sources, see sero_sources for example implementation
    call zero_sources(flow, phi, R, S)
    
  end subroutine eval_sources

end program scalar_transport
