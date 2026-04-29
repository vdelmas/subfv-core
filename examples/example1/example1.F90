program example1
  use mpi
  use precision_module
  use mpi_module
  use mesh_module
  use mesh_reading_module
  use mesh_geometry_module
  use mesh_connectivity_module
  implicit none

  integer(kind=ENTIER) :: me, num_procs, mpi_ierr
  logical :: boundary_2d = .false.

  type(mesh_type) :: mesh
  type(mpi_send_recv_type) :: mpi_send_recv
  character(len=255) :: meshfile, meshfile_path

  integer(kind=ENTIER) :: n_bc
  character(len=255), dimension(:), allocatable :: bc_name

  integer(kind=ENTIER) :: i, j, k, le, re
  integer(kind=ENTIER) :: id_face, id_vert, id_sub_face, id_sub_elem
  integer(kind=ENTIER) :: n_faces_interior, n_faces_bc_unassigned
  integer(kind=ENTIER), dimension(:), allocatable :: n_faces_bc
  real(kind=DOUBLE), dimension(3) :: maxsumnorm, sumnorm, norm

  meshfile_path=""
  meshfile="example1.msh"

  n_bc = 2
  allocate(bc_name(n_bc))
  bc_name(1) = "left_surf"
  bc_name(2) = "right_surf"

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  call read_mesh_msh(mesh, meshfile_path, meshfile, &
    n_bc, bc_name, me, num_procs, mpi_send_recv)
  call build_mesh(mesh, num_procs, mpi_send_recv, .true., boundary_2d)
  call compute_geometry_mesh(mesh, .true., boundary_2d)

  !Test mesh using faces
  maxsumnorm = 0.0_DOUBLE
  do i=1, mesh%n_elems
    sumnorm = 0.0_DOUBLE
    do j=1, mesh%elem(i)%n_faces
      id_face = mesh%elem(i)%face(j)
      if( mesh%face(id_face)%left_neigh == i ) then
        norm = mesh%face(id_face)%norm
      else
        norm = -mesh%face(id_face)%norm
      end if
      sumnorm = sumnorm + mesh%face(id_face)%area*norm
    end do
    if( norm2(sumnorm) > norm2(maxsumnorm) ) maxsumnorm = sumnorm
  end do
  print*, "Maximum error using faces", norm2(maxsumnorm)

  !Test mesh using subfaces
  maxsumnorm = 0.0_DOUBLE
  do i=1, mesh%n_elems
    sumnorm = 0.0_DOUBLE
    do j=1, mesh%elem(i)%n_sub_elems
      id_sub_elem = mesh%elem(i)%sub_elem(j)
      id_vert = mesh%sub_elem(id_sub_elem)%mesh_vert
      do k=1, mesh%sub_elem(id_sub_elem)%n_sub_faces
        id_sub_face = mesh%sub_elem(id_sub_elem)%sub_face(k)
        id_face = mesh%sub_face(id_sub_face)%mesh_face
        if( mesh%face(id_face)%left_neigh == i ) then
          norm = mesh%face(id_face)%norm
        else
          norm = -mesh%face(id_face)%norm
        end if
        sumnorm = sumnorm + mesh%sub_face(id_sub_face)%area*norm
      end do
    end do
    if( norm2(sumnorm) > norm2(maxsumnorm) ) maxsumnorm = sumnorm
  end do
  print*, "Maximum error using subfaces", norm2(maxsumnorm)

  !Count faces belonging to each boundary
  n_faces_interior = 0
  allocate(n_faces_bc(n_bc))
  n_faces_bc = 0
  n_faces_bc_unassigned = 0
  do i=1, mesh%n_faces
    le = mesh%face(i)%left_neigh
    re = mesh%face(i)%right_neigh
    if( re > 0 ) then
      n_faces_interior = n_faces_interior + 1
    else
      if( re == 0 ) then
        n_faces_bc_unassigned = n_faces_bc_unassigned + 1
      else
        n_faces_bc(-re) = n_faces_bc(-re) + 1
      end if
    end if
  end do

  print*, "Total number of faces", mesh%n_faces
  print*, "Interior faces", n_faces_interior
  print*, "Unassigned boundary faces", n_faces_bc_unassigned
  do i=1, n_bc
    print*, "Faces belonging to bc ", trim(adjustl(bc_name(i))), n_faces_bc(i)
  end do

  call MPI_FINALIZE(mpi_ierr)
end program example1
