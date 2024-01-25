!v Test to for finding particle function
!
!  @build mpi+petsc

program test_square_mesh_face_centres
#include "ccs_macros.inc"

  use petscvec
  use petscsys

  use testing_lib
  use constants, only: cell, face, ccsconfig, ccs_string_len
  use kinds, only: ccs_real, ccs_int
  use types, only: field, upwind_field, central_field, face_field, ccs_mesh, &
                   vector_spec, ccs_vector
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, timer, &
                      read_command_line_arguments, sync, &
                      create_shared_array
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use mesh_utils, only: build_square_mesh
  use meshing, only: set_mesh_object, nullify_mesh_object, get_local_index, get_natural_index, &
                     get_centre, create_cell_locator, get_global_index, create_neighbour_locator, &
                     create_face_locator, count_neighbours, get_face_normal, get_face_area
  use vec, only: create_vector, set_vector_location, &
                 get_vector_data, restore_vector_data
  use petsctypes, only: vector_petsc
  use pv_coupling, only: solve_nonlinear
  use utils, only: set_size, initialise, update, exit_print, debug_print, str
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays

  implicit none

  type(vector_spec) :: vec_properties
  type(cell_locator) :: loc_p
  type(neighbour_locator) :: loc_nb
  type(face_locator) :: loc_f

  integer(ccs_int) :: cps = 20 ! Default value for cells per side
  integer(ccs_int) :: part_cell_index
  integer(ccs_int) :: part_index
  integer(ccs_int) :: n_parts_global
  integer(ccs_int) :: n_parts_local
  integer(ccs_int) :: index_p, index_nb
  integer(ccs_int) :: i
  integer(ccs_int) :: n_particle_ranks
  integer(ccs_int) :: split_rank
  integer(ccs_int) :: nnb

  real(ccs_real), dimension(3) :: x_p, x_f, x_nb, norm, calc_norm
  real(ccs_real) :: dx, x_1, x_2
  real(ccs_real) :: domain_size
  real(ccs_real) :: h
  logical :: use_mpi_splitting
  
  ! Launch MPI
  call init()
  
  use_mpi_splitting = .false.
  call create_new_par_env(par_env, ccs_split_type_low_high, use_mpi_splitting, shared_env)

  domain_size = 1.0_ccs_real
  dx = domain_size/cps

  ! Create a square mesh
  mesh = build_square_mesh(par_env, shared_env, cps, domain_size)

  ! Use a random cell in interior to check whether faces and neighbours have correct coordinates
  index_p = 35
  call create_cell_locator(index_p, loc_p)
  call get_centre(loc_p, x_p)
  call count_neighbours(loc_p, nnb)
  h = mesh%geo%h
  do i = 1, nnb
    call create_neighbour_locator(loc_p, i, loc_nb)
    call create_face_locator(index_p, i, loc_f)

    call get_centre(loc_nb, x_nb)

    if (.not. (abs(abs(x_nb(1) - x_p(1)) - h) < eps .or. abs(abs(x_nb(2) - x_p(2)) - h) < eps)) then
      call error_abort("neighbour centre incorrect. x_p " // str(x_p(1)) // " " // str(x_p(2)) // " x_nb " // str(x_nb(1)) // " " // str(x_nb(2)) // " h " // str(h) // " index_p " // str(index_p) // " offset " // str(mesh%topo%shared_array_local_offset) )
    end if
  end do

  call dprint("done")

  ! Finalise MPI
  call nullify_mesh_object()
  call fin()

end program test_square_mesh_face_centres
