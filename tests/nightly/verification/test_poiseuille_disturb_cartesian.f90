!> @brief Test grid convergence of the 2d poiseuille test case that includes a prism layer
!
!> @description Generates several square meshes varying cps and get the error to compute their convergence order
program test_poiseuille_disturb_cartesian
#include "ccs_macros.inc"

  use ccs_base, only: bnd_names_default
  use testing_lib
  use error_analysis, only: get_order, print_error_summary, disturb_cartesian
  use mesh_utils, only: build_square_mesh
  use poiseuille_core, only: run_poiseuille, domain_size
  use mesh_utils, only: compute_face_interpolation
  use meshing, only: get_total_num_cells, set_mesh_object, nullify_mesh_object

  implicit none

  integer(ccs_int), parameter :: num_cps = 3
  integer(ccs_int), parameter :: nvar = 3
  real(ccs_real), dimension(nvar, num_cps) :: error_L2
  real(ccs_real), dimension(nvar, num_cps) :: error_Linf
  real(ccs_real), dimension(:), allocatable :: orders_L2
  real(ccs_real), dimension(:), allocatable :: orders_Linf
  real(ccs_real), dimension(num_cps) :: refinements
  integer(ccs_int), dimension(num_cps) :: cps_list
  integer(ccs_int) :: cps

  character(len=12), dimension(nvar) :: variable_labels

  integer(ccs_int) :: i

  call init()

  variable_labels = (/"U", "V", "P"/)
  domain_size = 1.0_ccs_real

  cps_list = (/32, 64, 128/) !, 256/)
  refinements = real(maxval(cps_list(:))) / real(cps_list(:))

  error_L2(:, :) = 0.0_ccs_real
  error_Linf(:, :) = 0.0_ccs_real

  do i = 1, num_cps
    cps = cps_list(i)
    mesh = build_square_mesh(par_env, shared_env, cps, domain_size, &
         bnd_names_default(1:4))

    call set_mesh_object(mesh)
    call disturb_cartesian(cps, domain_size, mesh)
    call nullify_mesh_object()
    
    call run_poiseuille(par_env, shared_env, error_L2(:, i), error_Linf(:, i), input_mesh=mesh)
  end do

  if (par_env%proc_id == par_env%root) then

    call print_error_summary(variable_labels, refinements, error_L2, error_Linf)

    call get_order(refinements, error_L2, orders_L2)
    call get_order(refinements, error_Linf, orders_Linf)

    call assert_gt(orders_L2(1), 1.9_ccs_real, "U not converging in 2nd order ")
    call assert_gt(orders_L2(2), 1.9_ccs_real, "V not converging in 2nd order ")
    call assert_gt(orders_L2(3), 1.9_ccs_real, "P not converging in 2nd order ")

    call assert_gt(orders_Linf(1), 1.8_ccs_real, "U not converging in 2nd order ")
    call assert_gt(orders_Linf(2), 1.8_ccs_real, "V not converging in 2nd order ")
    call assert_gt(orders_Linf(3), 1.1_ccs_real, "P not converging in 2nd order ")

  end if

  call fin()


end program test_poiseuille_disturb_cartesian
