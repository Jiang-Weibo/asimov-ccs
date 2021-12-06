module mesh_utils

  use kinds, only: accs_int, accs_real, accs_err
  use types, only: mesh
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  
  implicit none

  private
  public :: build_square_mesh
  public :: face_normal
  
contains

  function build_square_mesh(nps, l, par_env) result(square_mesh)

    class(parallel_environment) :: par_env
    integer(accs_int), intent(in) :: nps
    real(accs_real), intent(in) :: l

    integer(accs_int) :: istart, iend
    integer(accs_int) :: i, ii, ictr
    integer(accs_int) :: fctr
    integer(accs_int) :: comm_rank, comm_size

    type(mesh) :: square_mesh

    select type(par_env)
      type is (parallel_environment_mpi)

        square_mesh%n = nps**2            !> Number of cells
        square_mesh%h = l / nps

        associate(n=>square_mesh%n, &
                  h=>square_mesh%h)
          
          !! Setup ownership range
          comm_rank = par_env%proc_id
          comm_size = par_env%num_procs
          istart = comm_rank * (n / comm_size)
          if (modulo(square_mesh%n, comm_size) < comm_rank) then
            istart = istart + modulo(n, comm_size)
          else
            istart = istart + comm_rank
          end if
          istart = istart + 1 ! Fortran - 1 indexed
      
          iend = istart + n / comm_size
          if (modulo(square_mesh%n, comm_size) > comm_rank) then
            iend = iend + 1
          end if
          iend = iend - 1

          square_mesh%nlocal = (iend - istart) + 1

          allocate(square_mesh%idx_global(square_mesh%nlocal))
          allocate(square_mesh%nnb(square_mesh%nlocal))
          allocate(square_mesh%nbidx(4, square_mesh%nlocal))
          allocate(square_mesh%xc(2, square_mesh%nlocal))    !> @note Currently hardcoded as a 2D mesh!
          allocate(square_mesh%xf(2, 4, square_mesh%nlocal)) !> @note Currently hardcoded as a 2D mesh!
          allocate(square_mesh%vol(square_mesh%nlocal))
          allocate(square_mesh%Af(4, square_mesh%nlocal))    !> @note Currently hardcoded as a 2D mesh!
          allocate(square_mesh%nf(2, 4, square_mesh%nlocal)) !> @note Currently hardcoded as a 2D mesh!
          
          square_mesh%nnb(:) = 4 ! All cells have 4 neighbours (possibly ghost/boundary cells)
          square_mesh%vol(:) = square_mesh%h**2 !> @note Mesh is square and 2D
          square_mesh%Af(:, :) = square_mesh%h  !> @note Mesh is square and 2D
          
          !! Get neighbour indices
          !! XXX: These are global indices and thus may be off-process
          ictr = 1
          do i = istart, iend
            square_mesh%idx_global(ictr) = i
            ii = i - 1

            associate(xc => square_mesh%xc(:, ictr), &
                 xf => square_mesh%xf(:, :, ictr), &
                 nrm => square_mesh%nf(:, :, ictr))
              !! Set cell centre
              xc(1) = (modulo(ii, nps) + 0.5_accs_real) * h
              xc(2) = (ii / nps + 0.5_accs_real) * h
              
              !! Left neighbour
              fctr = 1
              if (modulo(ii, nps) == 0) then
                square_mesh%nbidx(fctr, ictr) = -1
              else
                square_mesh%nbidx(fctr, ictr) = i - 1
                if (square_mesh%nbidx(fctr, ictr) < 1) then
                  print *, "ERROR: interior neighbour idx < 1!"
                  stop
                end if
              end if
              xf(1, fctr) = xc(1) - 0.5_accs_real * h
              xf(2, fctr) = xc(2)
              nrm(1, fctr) = -1.0_accs_real
              nrm(2, fctr) = 0.0_accs_real

              !! Right neighbour
              fctr = 2
              if (modulo(ii, nps) == (nps - 1)) then
                square_mesh%nbidx(fctr, ictr) = -2
              else
                square_mesh%nbidx(fctr, ictr) = i + 1
                if (square_mesh%nbidx(fctr, ictr) > n) then
                  print *, "ERROR: interior neibour idx > N!"
                  stop
                end if
              end if
              xf(1, fctr) = xc(1) + 0.5_accs_real * h
              xf(2, fctr) = xc(2)
              nrm(1, fctr) = 1.0_accs_real
              nrm(2, fctr) = 0.0_accs_real

              !! Down neighbour
              fctr = 3
              if ((ii / nps) == 0) then
                square_mesh%nbidx(fctr, ictr) = -3
              else
                square_mesh%nbidx(fctr, ictr) = i - nps
                if (square_mesh%nbidx(fctr, ictr) < 1) then
                  print *, "ERROR: interior neighbour idx < 1!"
                  stop
                end if
              end if
              xf(1, fctr) = xc(1)
              xf(2, fctr) = xc(2) - 0.5_accs_real * h
              nrm(1, fctr) = 0.0_accs_real
              nrm(2, fctr) = -1.0_accs_real

              !! Up neighbour
              fctr = 4
              if ((ii / nps) == (nps - 1)) then
                square_mesh%nbidx(fctr, ictr) = -4
              else
                square_mesh%nbidx(fctr, ictr) = i + nps
                if (square_mesh%nbidx(fctr, ictr) > n) then
                  print *, "ERROR: interior neibour idx > N!"
                  stop
                end if
              end if
              xf(1, fctr) = xc(1)
              xf(2, fctr) = xc(2) + 0.5_accs_real * h
              nrm(1, fctr) = 0.0_accs_real
              nrm(2, fctr) = 1.0_accs_real
            end associate

            ictr = ictr + 1
          end do
        end associate

      class default
        print *, "Unknown parallel environment type!"
        stop

    end select    
  end function build_square_mesh

  function face_normal(mesh_obj, cell, face)

    type(mesh), intent(in) :: mesh_obj
    integer(accs_int), intent(in) :: cell
    integer(accs_int), intent(in) :: face

    real(accs_real), dimension(3) :: face_normal

    face_normal(1:2) = mesh_obj%nf(1:2, face, cell)
    
  end function face_normal
  
end module mesh_utils
