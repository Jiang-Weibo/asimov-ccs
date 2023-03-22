!v Module file pv_coupling.mod
!
!  An interface to pressure-velocity coupling methods (SIMPLE, etc)

module pv_coupling

  use kinds, only: ccs_int, ccs_real
  use types, only: field, ccs_mesh, fluid, fluid_solver_selector
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: solve_nonlinear

  interface

    module subroutine solve_nonlinear(par_env, mesh, it_start, it_end, res_target, &
                                      flow_solver_selector, flow, step)
      class(parallel_environment), allocatable, intent(in) :: par_env
      type(ccs_mesh), intent(in) :: mesh
      integer(ccs_int), intent(in) :: it_start, it_end
      real(ccs_real), intent(in) :: res_target
      type(fluid_solver_selector), intent(in) :: flow_solver_selector
      type(fluid), intent(inout) :: flow
      integer(ccs_int), optional, intent(in) :: step
    end subroutine solve_nonlinear

  end interface

end module pv_coupling
