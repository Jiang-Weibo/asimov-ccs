!> @brief Test finding the faces on a specific boundary using the boundary name.
program find_named_faces

  use testing_lib
  use meshing, only: find_entities, set_mesh_object
  use mesh_utils, only: build_square_mesh
  
  implicit none

  integer, parameter :: cps = 4
  real(ccs_real), parameter :: l = 1.0
  character(len=128), dimension(4) :: bnd_names
  
  type(face_locator), dimension(:), allocatable :: faces

  integer :: i
  
  call init()

  ! Check we are running on one process, this just makes the counting easier rather than anything
  ! else.
  if (par_env%num_procs > 1) then
    call stop_test("This test should be run in serial")
  end if
  
  bnd_names(1) = "alpha"
  bnd_names(2) = "beta"
  bnd_names(3) = "gamma"
  bnd_names(4) = "delta"
  
  mesh = build_square_mesh(par_env, shared_env, cps, l, bnd_names)
  call set_mesh_object(mesh)
  
  do i = 1, size(bnd_names)
    call find_entities(mesh, bnd_names(1), faces)
    if (size(faces) /= cps) then
      call stop_test("Incorrect boundary face count")
    end if
  end do
  
  call fin()
  
end program find_named_faces
