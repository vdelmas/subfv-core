module example2_module
  use subfv_precision_module
  implicit none
  real(kind=DOUBLE), parameter :: gamma = 1.4_DOUBLE

contains
  subroutine riemann_solver(ul, ur, n, f, amax)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in):: ul, ur
    real(kind=DOUBLE), dimension(3), intent(in) :: n
    real(kind=DOUBLE), dimension(5), intent(inout) :: f
    real(kind=DOUBLE), intent(inout) :: amax

    real(kind=DOUBLE) :: al, ar, vl, vr
    real(kind=DOUBLE), dimension(5) :: wl, wr, fnl, fnr

    wl = cons_to_prim(ul)
    al = sqrt(gamma*wl(5)/wl(1))
    vl = dot_product(wl(2:4), n)
    fnl = vl*ul
    fnl(2:4) = fnl(2:4) + wl(5)*n
    fnl(5) = fnl(5) + wl(5)*vl

    wr = cons_to_prim(ur)
    ar = sqrt(gamma*wr(5)/wr(1))
    vr = dot_product(wr(2:4), n)
    fnr = vr*ur
    fnr(2:4) = fnr(2:4) + wr(5)*n
    fnr(5) = fnr(5) + wr(5)*vr

    amax = max(abs(vl)+al, abs(vr)+ar)
    f = 0.5_DOUBLE*(fnl+fnr) - 0.5_DOUBLE*amax*(ur-ul)
  end subroutine riemann_solver

  function prim_to_cons(w) result(u)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), dimension(5) :: u

    u(1) = w(1)
    u(2:4) = w(1)*w(2:4)
    u(5) = w(1)*(w(5)/(w(1)*(gamma-1.0_DOUBLE)) &
      + 0.5_DOUBLE*dot_product(w(2:4),w(2:4)))
  end function prim_to_cons

  function cons_to_prim(u) result(w)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), dimension(5) :: w

    w(1) = u(1)
    w(2:4) = u(2:4)/u(1)
    w(5) = u(1)*(gamma-1.0_DOUBLE)*(u(5)/u(1) &
      - 0.5_DOUBLE*dot_product(w(2:4),w(2:4)))
  end function cons_to_prim

  subroutine write_sol(filename, mesh, sol)
    use subfv_mesh_module
    use subfv_io_module
    implicit none

    character(len=*), intent(in) :: filename
    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol

    integer(kind=ENTIER) :: i
    integer(kind=ENTIER) :: fn_vtu, fn_pvtu

    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(:), allocatable :: rho
    real(kind=DOUBLE), dimension(:,:), allocatable :: vel
    real(kind=DOUBLE), dimension(:), allocatable :: p, ie

    call open_file_vtu(mesh, trim(adjustl(filename)), fn_vtu, fn_pvtu)
    call write_file_vtu_start_cell_data(mesh, trim(adjustl(filename)), &
      fn_vtu, fn_pvtu)

    allocate(rho(mesh%n_elems))
    do i=1, mesh%n_elems
      rho(i) = sol(1, i)
    end do
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(filename)), &
      fn_vtu, fn_pvtu, rho, "Density")
    deallocate(rho)

    allocate(vel(3, mesh%n_elems))
    do i=1, mesh%n_elems
      vel(:, i) = sol(2:4, i)/sol(1, i)
    end do
    call write_file_vtu_cell_vector(mesh, trim(adjustl(filename)), &
      fn_vtu, fn_pvtu, vel, "Velocity")
    deallocate(vel)

    allocate(p(mesh%n_elems))
    do i=1, mesh%n_elems
      w = cons_to_prim(sol(:, i))
      p(i) = w(5)
    end do
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(filename)), &
      fn_vtu, fn_pvtu, p, "Pressure")
    deallocate(p)

    allocate(ie(mesh%n_elems))
    do i=1, mesh%n_elems
      w = cons_to_prim(sol(:, i))
      ie(i) = w(5)/(w(1)*(gamma-1.0_DOUBLE))
    end do
    call write_file_vtu_cell_scalar(mesh, trim(adjustl(filename)), &
      fn_vtu, fn_pvtu, ie, "Internal energy")
    deallocate(ie)

    call write_file_vtu_end_cell_data(mesh, trim(adjustl(filename)), &
      fn_vtu, fn_pvtu)
    call close_file_vtu(mesh, trim(adjustl(filename)), fn_vtu, fn_pvtu)
  end subroutine write_sol
end module example2_module
