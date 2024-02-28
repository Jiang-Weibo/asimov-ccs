!> Test face values
program test_face_values

  use testing_lib
  use ccs_base, only: bnd_names_default
  use kinds, only: ccs_int, ccs_real
  use constants, only: face
  use types, only: vector_spec, ccs_mesh, field, face_field
  use mesh_utils, only: build_square_mesh
  use vec, only: create_vector, set_vector_location
  use meshing, only: create_neighbour_locator, &
                     get_global_index, get_local_index, get_face_area, get_face_normal
  use meshing, only: set_mesh_object, nullify_mesh_object
  use utils, only: initialise, set_size

  implicit none

  class(field), allocatable :: mf

  type(vector_spec) :: vec_properties

  ! integer(ccs_int) :: nfaces
  integer(ccs_int) :: cps = 5 !< Cells per side of the mesh

  call init()

  ! Create a square mesh
  mesh = build_square_mesh(par_env, shared_env, cps, 1.0_ccs_real, &
       bnd_names_default(1:4))
  call set_mesh_object(mesh)

  allocate (face_field :: mf)

  call initialise(vec_properties)

  ! Setup vector size to store face-centred values (rather than cell-centred values)
  call set_vector_location(face, vec_properties)

  call set_size(par_env, mesh, vec_properties)
  call create_vector(vec_properties, mf%values)

  call nullify_mesh_object()
  call fin()

end program test_face_values
