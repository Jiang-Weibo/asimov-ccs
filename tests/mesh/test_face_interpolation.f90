!> @brief Test the face interpolation with a 2 cells mesh
!
!> @description Generates a basic mesh geoemtry and topology for a 2 cells mesh and makes sure the face interpolation is
!!              properly computed and accessed (with get_face_interpolation)
program test_face_interpolation

  use testing_lib

  use constants, only: ndim
  use mesh_utils, only: compute_face_interpolation
  use meshing, only: get_face_interpolation, set_face_location, get_max_faces

  implicit none

  type(ccs_mesh), target :: mesh
  type(face_locator) :: loc_f
  integer(ccs_int) :: index_p, j
  real(ccs_real) :: interpol_factor
  real(ccs_real) :: face_coordinate

  call init()

  face_coordinate = 0.3_ccs_real
  mesh = generate_mesh(face_coordinate)

  call compute_face_interpolation(mesh)

  ! Test cell1 interpolation factor
  index_p = 1_ccs_int
  j = 1_ccs_int

  call set_face_location(mesh, index_p, j, loc_f)
  call get_face_interpolation(loc_f, interpol_factor)

  call assert_eq(interpol_factor, 1.0_ccs_real - face_coordinate, "FAIL: wrong interpolation for cell1")

  ! Test cell2 interpolation factor
  index_p = 2_ccs_int
  j = 1_ccs_int

  call set_face_location(mesh, index_p, j, loc_f)
  call get_face_interpolation(loc_f, interpol_factor)

  call assert_eq(interpol_factor, face_coordinate, "FAIL: wrong interpolation for cell2")

  call fin()

contains

  !v Generates a 2 cells mesh with cell1 at (0,0,0), cell2 at (1,0,0) and face1(face_coordinate, 0, 0)
  function generate_mesh(face_coordinate) result(mesh)

    use meshing, only: get_global_num_cells, set_global_num_cells, &
                       get_local_num_cells, set_local_num_cells, &
                       set_halo_num_cells, &
                       set_global_num_faces, &
                       get_total_num_cells, set_total_num_cells, &
                       set_num_faces, set_max_faces, &
                       set_cell_location, set_neighbour_location, &
                       set_local_index
    use types, only: cell_locator, neighbour_locator
    
    real(ccs_real), intent(in) :: face_coordinate
    type(ccs_mesh) :: mesh

    integer(ccs_int) :: global_num_cells
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: total_num_cells
    integer(ccs_int) :: max_faces

    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    
    ! Build 2 cells mesh topology
    call set_local_num_cells(2, mesh)
    call set_global_num_cells(2, mesh)
    call set_halo_num_cells(0, mesh)
    call set_total_num_cells(2, mesh)
    call set_global_num_faces(1, mesh)
    call set_num_faces(1, mesh)
    call set_max_faces(1, mesh)
    
    call get_local_num_cells(mesh, local_num_cells)
    call get_total_num_cells(mesh, total_num_cells)
    call get_global_num_cells(mesh, global_num_cells)
    call get_max_faces(mesh, max_faces)
    
    allocate (mesh%topo%global_indices(global_num_cells))
    mesh%topo%global_indices(1) = 1
    mesh%topo%global_indices(2) = 2

    allocate (mesh%topo%global_face_indices(max_faces, global_num_cells))
    mesh%topo%global_face_indices(:, :) = 1

    allocate (mesh%topo%face_indices(max_faces, local_num_cells))
    mesh%topo%face_indices(:, :) = 1

    allocate (mesh%topo%nb_indices(max_faces, local_num_cells))
    call set_cell_location(mesh, 1, loc_p)        ! Select cell 1
    call set_neighbour_location(loc_p, 1, loc_nb) ! Select 1st neighbour
    call set_local_index(2, loc_nb)               ! Set neighbour index = 2
    call set_cell_location(mesh, 2, loc_p)        ! Select cell 2
    call set_neighbour_location(loc_p, 1, loc_nb) ! Select 1st neighbour
    call set_local_index(1, loc_nb)               ! Set neighbour index = 1

    allocate (mesh%topo%num_nb(local_num_cells))
    mesh%topo%num_nb(:) = 1

    ! Build 2 cells mesh geometry
    allocate (mesh%geo%x_p(ndim, total_num_cells))
    mesh%geo%x_p(:, 1) = (/0.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real/)
    mesh%geo%x_p(:, 2) = (/1.0_ccs_real, 0.0_ccs_real, 0.0_ccs_real/)

    allocate (mesh%geo%x_f(ndim, max_faces, total_num_cells))
    mesh%geo%x_f(:, 1, 1) = (/face_coordinate, 0.0_ccs_real, 0.0_ccs_real/)
    mesh%geo%x_f(:, 1, 2) = (/face_coordinate, 0.0_ccs_real, 0.0_ccs_real/)

  end function

end program test_face_interpolation
