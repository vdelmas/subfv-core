program main
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_geometry_module
  use mesh_connectivity_module
  use io_module
  use example2_module
  implicit none

  integer(kind=ENTIER) :: me, num_procs, mpi_ierr
  logical :: boundary_2d = .true.

  type(mesh_type) :: mesh
  type(mpi_send_recv_type) :: mpi_send_recv
  character(len=255) :: meshfile, meshfile_path

  integer(kind=ENTIER) :: n_bc
  character(len=255), dimension(:), allocatable :: bc_name

  integer(kind=ENTIER) :: i, j, k, le, re
  integer(kind=ENTIER) :: id_face, id_vert, id_sub_face, id_sub_elem
  real(kind=DOUBLE), dimension(3) :: norm

  real(kind=DOUBLE), parameter :: cfl = 0.95
  real(kind=DOUBLE) :: t, tmax, dt, a, amax
  real(kind=DOUBLE), dimension(5) :: w, sol_inflow, solr, soll
  real(kind=DOUBLE), dimension(5) :: f
  real(kind=DOUBLE), dimension(:, :), allocatable :: sol, rhs

  integer(kind=ENTIER) :: fn_vtu, fn_pvtu
  character(len=255) :: filename

  meshfile_path=""
  meshfile="example2.msh"

  n_bc = 2
  allocate(bc_name(n_bc))
  bc_name(1) = "inflow_surf"
  bc_name(2) = "outflow_surf"

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  call read_mesh_msh(mesh, meshfile_path, meshfile, &
    n_bc, bc_name, me, num_procs, mpi_send_recv)
  call build_mesh(mesh, num_procs, mpi_send_recv, .true., boundary_2d)
  call compute_geometry_mesh(mesh, .true., boundary_2d)

  allocate(sol(5, mesh%n_elems))
  allocate(rhs(5, mesh%n_elems))

  do i=1, mesh%n_elems
    w(1) = 1.4_DOUBLE
    w(2:4) = (/0.0_DOUBLE, -20.0_DOUBLE, 0.0_DOUBLE/)
    w(5) = 1.0_DOUBLE
    sol(:, i) = prim_to_cons(w)
  end do

  w(1) = 1.4_DOUBLE
  w(2:4) = (/0.0_DOUBLE, -20.0_DOUBLE, 0.0_DOUBLE/)
  w(5) = 1.0_DOUBLE
  sol_inflow = prim_to_cons(w)

  t = 0.0_DOUBLE
  tmax = 1.0_DOUBLE
  do while ( t < tmax ) 
    dt = 1e10_DOUBLE
    rhs = 0.0_DOUBLE
    do i=1, mesh%n_faces
      le = mesh%face(i)%left_neigh
      re = mesh%face(i)%right_neigh
      norm = mesh%face(i)%norm

      soll = sol(:, le)
      if( re > 0 ) then
        solr = sol(:, re)
      else
        if( re == 0 ) then
          solr = soll
          solr(2:4) = solr(2:4) &
            - 2.0_DOUBLE*dot_product(solr(2:4), norm)*norm
        else
          if(trim(adjustl(bc_name(-re))) == "inflow_surf") then
            solr = sol_inflow
          else if(trim(adjustl(bc_name(-re))) == "outflow_surf") then
            solr = soll
          end if
        end if
      end if

      call riemann_solver(soll, solr, norm, f, amax)

      rhs(:, le) = rhs(:, le) &
        - mesh%face(i)%area/mesh%elem(le)%volume*f
      dt = min(dt, cfl*mesh%elem(le)%volume/(mesh%face(i)%area*amax))
      if( re > 0 ) then
        rhs(:, re) = rhs(:, re) &
          + mesh%face(i)%area/mesh%elem(re)%volume*f
        dt = min(dt, cfl*mesh%elem(re)%volume/(mesh%face(i)%area*amax))
      end if
    end do

    sol = sol + dt * rhs
    t = t + dt
  end do

  call write_sol("final", mesh, sol)
  call MPI_FINALIZE(mpi_ierr)
end program main
