!v Module file fv.mod
!
!  An interface to finite volume implementations (CDS, UDS, etc.)

module fv

  use kinds, only: ccs_real, ccs_int
  use types, only: ccs_matrix, ccs_vector, ccs_mesh, field, upwind_field, central_field, gamma_field, linear_upwind_field, bc_config, &
                   face_locator, cell_locator, neighbour_locator, bc_profile
  use constants, only: ndim

  implicit none

  private

  public :: compute_fluxes
  public :: calc_advection_coeff
  public :: calc_diffusion_coeff
  public :: calc_mass_flux
  public :: calc_cell_coords
  public :: update_gradient
  public :: compute_boundary_values
  public :: compute_boundary_coeffs
  public :: get_value_from_bc_profile
  public :: add_fixed_source
  public :: add_linear_source
  
  interface calc_advection_coeff
    module procedure calc_advection_coeff_cds
    module procedure calc_advection_coeff_uds
    module procedure calc_advection_coeff_gamma
    module procedure calc_advection_coeff_luds
  end interface calc_advection_coeff

  interface calc_mass_flux
    module procedure calc_mass_flux_uvw
    module procedure calc_mass_flux_no_uvw
  end interface calc_mass_flux

  interface

    !> Calculates advection coefficient for neighbouring cell using CDS discretisation
    module subroutine calc_advection_coeff_cds(phi, loc_f, mf, bc, coeffaP, coeffaF)
      type(central_field), intent(in) :: phi  !< scalar (central) field
      type(face_locator), intent(in) :: loc_f !< face locator
      real(ccs_real), intent(in) :: mf        !< mass flux at the face
      integer(ccs_int), intent(in) :: bc      !< flag indicating whether cell is on boundary
      real(ccs_real), intent(out) :: coeffaP  !< advection coefficient for current cell
      real(ccs_real), intent(out) :: coeffaF  !< advection coefficient for neighbour cell
    end subroutine calc_advection_coeff_cds

    !> Calculates advection coefficient for neighbouring cell using UDS discretisation
    module subroutine calc_advection_coeff_uds(phi, loc_f, mf, bc, coeffaP, coeffaF)
      type(upwind_field), intent(in) :: phi   !< scalar (upwind) field
      type(face_locator), intent(in) :: loc_f !< face locator
      real(ccs_real), intent(in) :: mf        !< mass flux at the face
      integer(ccs_int), intent(in) :: bc      !< flag indicating whether cell is on boundary
      real(ccs_real), intent(out) :: coeffaP  !< advection coefficient for current cell
      real(ccs_real), intent(out) :: coeffaF  !< advection coefficient for neighbour cell
    end subroutine calc_advection_coeff_uds

    !> Calculates advection coefficient for neighbouring cell using gamma discretisation
    module subroutine calc_advection_coeff_gamma(phi, loc_f, mf, bc, loc_p, loc_nb, coeffaP, coeffaF)
      type(gamma_field), intent(inout) :: phi       !< scalar (gamma) field
      type(face_locator), intent(in) :: loc_f       !< face locator
      real(ccs_real), intent(in) :: mf              !< mass flux at the face
      integer(ccs_int), intent(in) :: bc            !< flag indicating whether cell is on boundary
      type(cell_locator), intent(in) :: loc_p       !< current cell locator
      type(neighbour_locator), intent(in) :: loc_nb !< neighbour cell locator
      real(ccs_real), intent(out) :: coeffaP        !< advection coefficient for current cell
      real(ccs_real), intent(out) :: coeffaF        !< advection coefficient for neighbour cell
    end subroutine calc_advection_coeff_gamma

    !> Calculates advection coefficient for neighbouring cell using LUDS discretisation
    module subroutine calc_advection_coeff_luds(phi, loc_f, mf, bc, loc_p, loc_nb, coeffaP, coeffaF)
      type(linear_upwind_field), intent(inout) :: phi !< scalar (gamma) field
      type(face_locator), intent(in) :: loc_f         !< face locator
      real(ccs_real), intent(in) :: mf                !< mass flux at the face
      integer(ccs_int), intent(in) :: bc              !< flag indicating whether cell is on boundary
      type(cell_locator), intent(in) :: loc_p         !< current cell locator
      type(neighbour_locator), intent(in) :: loc_nb   !< neighbour cell locator
      real(ccs_real), intent(out) :: coeffaP          !< advection coefficient for current cell
      real(ccs_real), intent(out) :: coeffaF          !< advection coefficient for neighbour cell
    end subroutine calc_advection_coeff_luds

    !> Sets the diffusion coefficient
    ! XXX: why is this a function when the equivalent advection ones are subroutines?
    module function calc_diffusion_coeff(index_p, index_nb, mesh, enable_cell_corrections) result(coeff)
      integer(ccs_int), intent(in) :: index_p  !< the local cell index
      integer(ccs_int), intent(in) :: index_nb !< the local neigbouring cell index
      type(ccs_mesh), intent(in) :: mesh       !< the mesh structure
      logical, intent(in) :: enable_cell_corrections !< whether or not cell corrections shouls be used
      real(ccs_real) :: coeff                  !< the diffusion coefficient
    end function calc_diffusion_coeff

    !> Computes fluxes and assign to matrix and RHS
    module subroutine compute_fluxes(phi, mf, mesh, component, M, vec)
      class(field), intent(inout) :: phi             !< scalar field structure
      class(field), intent(inout) :: mf              !< mass flux field structure (defined at faces)
      type(ccs_mesh), intent(in) :: mesh          !< the mesh being used
      integer(ccs_int), intent(in) :: component   !< integer indicating direction of velocity field component
      class(ccs_matrix), intent(inout) :: M       !< Data structure containing matrix to be filled
      class(ccs_vector), intent(inout) :: vec     !< Data structure containing RHS vector to be filled
    end subroutine

    !> Calculates mass flux across given face. Note: assumes rho = 1 and uniform grid
    module function calc_mass_flux_uvw(u_field, v_field, w_field, p, dpdx, dpdy, dpdz, invAu, invAv, invAw, loc_f, enable_cell_corrections) result(flux)
      class(field), intent(inout) :: u_field               !< x velocities field
      class(field), intent(inout) :: v_field               !< y velocities field
      class(field), intent(inout) :: w_field               !< z velocities field
      real(ccs_real), dimension(:), intent(in) :: p     !< array containing pressure
      real(ccs_real), dimension(:), intent(in) :: dpdx  !< pressure gradients in x
      real(ccs_real), dimension(:), intent(in) :: dpdy  !< pressure gradients in y
      real(ccs_real), dimension(:), intent(in) :: dpdz  !< pressure gradients in z
      real(ccs_real), dimension(:), intent(in) :: invAu !< inverse momentum diagonal in x
      real(ccs_real), dimension(:), intent(in) :: invAv !< inverse momentum diagonal in y
      real(ccs_real), dimension(:), intent(in) :: invAw !< inverse momentum diagonal in z
      type(face_locator), intent(in) :: loc_f           !< face locator
      logical, intent(in) :: enable_cell_corrections    !< whether or not cell shape corrections are to be used
      real(ccs_real) :: flux                            !< the flux across the boundary
    end function calc_mass_flux_uvw

    !> Computes Rhie-Chow correction
    module function calc_mass_flux_no_uvw(p, dpdx, dpdy, dpdz, invAu, invAv, invAw, loc_f, enable_cell_corrections) result(flux)
      real(ccs_real), dimension(:), intent(in) :: p     !< array containing pressure
      real(ccs_real), dimension(:), intent(in) :: dpdx  !< pressure gradients in x
      real(ccs_real), dimension(:), intent(in) :: dpdy  !< pressure gradients in y
      real(ccs_real), dimension(:), intent(in) :: dpdz  !< pressure gradients in z
      real(ccs_real), dimension(:), intent(in) :: invAu !< inverse momentum diagonal in x
      real(ccs_real), dimension(:), intent(in) :: invAv !< inverse momentum diagonal in y
      real(ccs_real), dimension(:), intent(in) :: invAw !< inverse momentum diagonal in z
      type(face_locator), intent(in) :: loc_f           !< face locator
      logical, intent(in) :: enable_cell_corrections    !< whether or not cell shape corrections are to be used
      real(ccs_real) :: flux                            !< the flux across the boundary
    end function calc_mass_flux_no_uvw

    !> Calculates the row and column indices from flattened vector index. Assumes square mesh
    module subroutine calc_cell_coords(index, cps, row, col)
      integer(ccs_int), intent(in) :: index !< cell index
      integer(ccs_int), intent(in) :: cps   !< number of cells per side
      integer(ccs_int), intent(out) :: row  !< cell row within mesh
      integer(ccs_int), intent(out) :: col  !< cell column within mesh
    end subroutine calc_cell_coords

    !v Performs an update of the gradients of a field.
    !  @note This will perform a parallel update of the gradient fields to ensure halo cells are
    !  correctly updated on other PEs. @endnote
    module subroutine update_gradient(mesh, phi)
      type(ccs_mesh), intent(in) :: mesh !< the mesh
      class(field), intent(inout) :: phi !< the field whose gradients we want to update
    end subroutine update_gradient

    !> Computes the value of the scalar field on the boundary
    module subroutine compute_boundary_values(phi, component, loc_p, loc_f, normal, bc_value)
      class(field), intent(inout) :: phi                         !< the field for which boundary values are being computed
      integer(ccs_int), intent(in) :: component               !< integer indicating direction of velocity field component
      type(cell_locator), intent(in) :: loc_p                 !< location of cell
      type(face_locator), intent(in) :: loc_f                 !< location of face
      real(ccs_real), dimension(ndim), intent(in) :: normal   !< boundary face normal direction
      real(ccs_real), intent(out) :: bc_value                 !< the value of the scalar field at the specified boundary
    end subroutine
  
    module subroutine compute_boundary_coeffs(phi, component, loc_p, loc_f, normal, a, b)
      class(field), intent(inout) :: phi                      !< the field for which boundary values are being computed
      integer(ccs_int), intent(in) :: component               !< integer indicating direction of velocity field component
      type(cell_locator), intent(in) :: loc_p                 !< location of cell
      type(face_locator), intent(in) :: loc_f                 !< location of face
      real(ccs_real), dimension(ndim), intent(in) :: normal   !< boundary face normal direction
      real(ccs_real), intent(out) :: a                        !< The diagonal coeff (implicit)
      real(ccs_real), intent(out) :: b                        !< The RHS entry (explicit)
    end subroutine

    !> Linear interpolate of BC profile 
    module subroutine get_value_from_bc_profile(x, profile, bc_value)
        real(ccs_real), dimension(:), intent(in) :: x !< Location of the interpolation
        type(bc_profile), intent(in) :: profile       !< boundary condition profile
        real(ccs_real), intent(out) :: bc_value       !< Interpolated value
    end subroutine get_value_from_bc_profile


    !> Adds a fixed source term to the righthand side of the equation
    module subroutine add_fixed_source(mesh, S, rhs)
      type(ccs_mesh), intent(in) :: mesh     !< The mesh
      class(ccs_vector), intent(inout) :: S   !< The source field
      class(ccs_vector), intent(inout) :: rhs !< The righthand side vector
    end subroutine add_fixed_source

    !> Adds a linear source term to the system matrix
    module subroutine add_linear_source(mesh, S, M)
      type(ccs_mesh), intent(in) :: mesh    !< The mesh
      class(ccs_vector), intent(inout) :: S !< The source field
      class(ccs_matrix), intent(inout) :: M !< The system
    end subroutine add_linear_source
    
  end interface

end module fv
