module io_module
  implicit none
contains
  subroutine open_file_vtu(mesh, filename, fn_vtu, fn_pvtu)
    use mpi
    use mpi_module
    use mesh_module
    implicit none

    type(mesh_type) :: mesh
    character(len=255), intent(in) :: filename
    integer(kind=ENTIER), intent(inout) :: fn_vtu, fn_pvtu

    integer(kind=ENTIER) :: num_procs, me, mpi_ierr
    character(len=255) :: filename_vtu, filename_pvtu, me_str

    integer(kind=ENTIER) :: i, j, k, s, id_face
    integer(kind=ENTIER) :: n_interior_elems, n_interior_vert
    real(kind=DOUBLE) :: dmin, dmax

    integer(kind=ENTIER), dimension(:), allocatable :: id_vert_no_ghost
    integer(kind=ENTIER), dimension(:), allocatable :: local_id_vert_no_ghost

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      !PVTU header
      filename_pvtu = trim(adjustl(filename))//".pvtu"
      open(newunit=fn_pvtu, file=trim(adjustl(filename_pvtu)))

      write(fn_pvtu, *) "<VTKFile type='PUnstructuredGrid' version='0.1' byte_order='LittleEndian' header_type='UInt64'>"
      write(fn_pvtu, *) "<PUnstructuredGrid GhostLevel='0'>"
      write(fn_pvtu, *) "<PPoints>"
      write(fn_pvtu, *) "<PDataArray type='Float64' Name='Points' NumberOfComponents='3'/>"
      write(fn_pvtu, *) "</PPoints>"
      write(fn_pvtu, *) "<PCells>"
      write(fn_pvtu, *) "<PDataArray type='Int64' Name='connectivity' NumberOfComponents='1'/>"
      write(fn_pvtu, *) "<PDataArray type='Int64' Name='offsets'      NumberOfComponents='1'/>"
      write(fn_pvtu, *) "<PDataArray type='UInt8' Name='types'        NumberOfComponents='1'/>"
      write(fn_pvtu, *) "</PCells>"
    end if

    !VTU
    write(me_str, *) me
    filename_vtu = trim(adjustl(me_str))//"_"//trim(adjustl(filename))//".vtu"
    open(newunit=fn_vtu, file=trim(adjustl(filename_vtu)))

    write(fn_vtu, *) "<VTKFile type='UnstructuredGrid' version='1.0' &
      &byte_order='LittleEndian' header_type='UInt64'>"
    write(fn_vtu, *) "<UnstructuredGrid>\n<Piece NumberOfPoints='", n_interior_vert, &
      &"' NumberOfCells='", n_interior_elems, "'>"

    dmin = 1e100_DOUBLE
    dmax = -1e100_DOUBLE
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        do j=1, 3
          if( mesh%vert(i)%coord(j) > dmax ) dmax = mesh%vert(i)%coord(j)
          if( mesh%vert(i)%coord(j) < dmin ) dmin = mesh%vert(i)%coord(j)
        end do
      end if
    end do

    write(fn_vtu, *) "<Points>"
    write(fn_vtu, *) " <DataArray type='Float64' Name='Points' NumberOfComponents='3' &
      & format='ascii' RangeMin='", dmin, "' RangeMax='", dmax, "'>"

    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        write(fn_vtu, *) mesh%vert(i)%coord(:)
      end if
    end do

    write(fn_vtu, *) "</DataArray>"
    write(fn_vtu, *) "</Points>"
    write(fn_vtu, *) "<Cells>"
    write(fn_vtu, *) "<DataArray type='Int64' Name='connectivity' format='ascii' &
      &RangeMin='0' RangeMax='", mesh%n_vert-1, "'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn_vtu, *) id_vert_no_ghost(mesh%elem(i)%vert) - 1
      end if
    end do
    write(fn_vtu, *) "</DataArray>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        s = s + mesh%elem(i)%n_vert
      end if
    end do

    write(fn_vtu, *) "<DataArray type='Int64' Name='offsets' format='ascii' &
      &RangeMin='", 0, "' RangeMax='", s, "'>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        s = s + mesh%elem(i)%n_vert
        write(fn_vtu, *) s
      end if
    end do

    write(fn_vtu, *) "</DataArray>"

    write(fn_vtu, *) "<DataArray type='UInt8' Name='types' format='ascii' &
      &RangeMin='42' RangeMax='42'>"

    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn_vtu, *) 42 !Type for polyhedra
      end if
    end do

    write(fn_vtu, *) "</DataArray>"
    write(fn_vtu, *) "<DataArray type='Int64' Name='faces' format='ascii' &
      &RangeMin='0' RangeMax='", mesh%n_vert-1, "'>"

    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn_vtu, *) mesh%elem(i)%n_faces
        do j=1, mesh%elem(i)%n_faces
          id_face = mesh%elem(i)%face(j)
          allocate(local_id_vert_no_ghost(mesh%face(id_face)%n_vert))
          do k=1, mesh%face(id_face)%n_vert
            local_id_vert_no_ghost(k) = id_vert_no_ghost(mesh%face(id_face)%vert(k))
          end do
          !write(fn_vtu, *) mesh%face(id_face)%n_vert, mesh%face(id_face)%vert(:) - 1
          write(fn_vtu, *) mesh%face(id_face)%n_vert, local_id_vert_no_ghost - 1
          deallocate(local_id_vert_no_ghost)
        end do
      end if
    end do

  end subroutine open_file_vtu

  subroutine close_file_vtu(mesh, filename, fn_vtu, fn_pvtu)
    use mpi
    use mpi_module
    use mesh_module
    implicit none

    type(mesh_type) :: mesh
    character(len=255), intent(in) :: filename
    integer(kind=ENTIER), intent(inout) :: fn_vtu, fn_pvtu

    integer(kind=ENTIER) :: me, num_procs, mpi_ierr

    integer(kind=ENTIER) :: i
    character(len=255) :: i_char
    
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 1 ) then
      !PVTU header
      do i=0, num_procs-1
        write(i_char, *) i
        write(fn_pvtu, *) "<Piece Source='"//trim(adjustl(i_char))//&
          trim(adjustl(filename))//".vtu'/>"
      end do
      write(fn_pvtu, *) "</PUnstructuredGrid>"
      write(fn_pvtu, *) "</VTKFile>"
      close(fn_pvtu)
    end if

    !VTU
    write(fn_vtu, *) "</Piece>"
    write(fn_vtu, *) "</UnstructuredGrid>"
    write(fn_vtu, *) "</VTKFile>"
    close(fn_vtu)
  end subroutine close_file_vtu

!  subroutine write_sol_vtu(mesh, filename, sol, grad, all_nodal_grad, vp)
!    use mpi
!    implicit none
!
!    type(mesh_type), intent(in) :: mesh
!    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
!    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
!    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(inout) :: all_nodal_grad
!    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: vp
!
!    character(len=*), intent(in) :: filename
!
!
!    integer(kind=ENTIER) :: me, num_procs, mpi_ierr
!    integer(kind=ENTIER) :: n_interior_elems, n_interior_vert
!    integer(kind=ENTIER) :: fn, i, j, s, id_face, k, id_sub_elem
!    integer(kind=ENTIER) :: id_vert
!    real(kind=DOUBLE) :: dmin, dmax, a, b, div
!    real(kind=DOUBLE), dimension(3) :: rot
!
!    real(kind=DOUBLE), dimension(5) :: w
!    real(kind=DOUBLE), dimension(3, 5) :: nodal_grad
!    integer(kind=ENTIER), dimension(:), allocatable :: id_vert_no_ghost
!    integer(kind=ENTIER), dimension(:), allocatable :: local_id_vert_no_ghost
!
!    real(kind=DOUBLE), dimension(mesh%n_vert) :: cell_size
!
!
!    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
!    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)
!
!    allocate(id_vert_no_ghost(mesh%n_vert))
!    id_vert_no_ghost = 0
!    n_interior_vert = 0
!    do i=1, mesh%n_vert
!      if ( .not. mesh%vert(i)%is_ghost ) then
!        n_interior_vert = n_interior_vert + 1
!        id_vert_no_ghost(i) = n_interior_vert
!      end if
!    end do
!
!    n_interior_elems = 0
!    do i=1, mesh%n_elems
!      if ( .not. mesh%elem(i)%is_ghost ) n_interior_elems = n_interior_elems + 1
!    end do
!
!    open(newunit=fn, file=trim(adjustl(filename)))
!
!    write(fn, *) "<VTKFile type='UnstructuredGrid' version='1.0' &
!      &byte_order='LittleEndian' header_type='UInt64'>"
!    write(fn, *) "<UnstructuredGrid>\n<Piece NumberOfPoints='", n_interior_vert, &
!      &"' NumberOfCells='", n_interior_elems, "'>"
!
!    dmin = 1e100_DOUBLE
!    dmax = -1e100_DOUBLE
!    do i=1, mesh%n_vert
!      if( .not. mesh%vert(i)%is_ghost ) then
!        do j=1, 3
!          if( mesh%vert(i)%coord(j) > dmax ) dmax = mesh%vert(i)%coord(j)
!          if( mesh%vert(i)%coord(j) < dmin ) dmin = mesh%vert(i)%coord(j)
!        end do
!      end if
!    end do
!
!    write(fn, *) "<Points>"
!    write(fn, *) " <DataArray type='Float64' Name='Points' NumberOfComponents='3' &
!      & format='ascii' RangeMin='", dmin, "' RangeMax='", dmax, "'>"
!
!    do i=1, mesh%n_vert
!      if( .not. mesh%vert(i)%is_ghost ) then
!        write(fn, *) mesh%vert(i)%coord(:)
!      end if
!    end do
!
!    write(fn, *) "</DataArray>"
!    write(fn, *) "</Points>"
!    write(fn, *) "<Cells>"
!    write(fn, *) "<DataArray type='Int64' Name='connectivity' format='ascii' &
!      &RangeMin='0' RangeMax='", mesh%n_vert-1, "'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write(fn, *) id_vert_no_ghost(mesh%elem(i)%vert) - 1
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    s = 0
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        s = s + mesh%elem(i)%n_vert
!      end if
!    end do
!
!    write(fn, *) "<DataArray type='Int64' Name='offsets' format='ascii' &
!      &RangeMin='", 0, "' RangeMax='", s, "'>"
!
!    s = 0
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        s = s + mesh%elem(i)%n_vert
!        write(fn, *) s
!      end if
!    end do
!
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='UInt8' Name='types' format='ascii' &
!      &RangeMin='42' RangeMax='42'>"
!
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write(fn, *) 42 !Type for polyhedra
!      end if
!    end do
!
!    write(fn, *) "</DataArray>"
!    write(fn, *) "<DataArray type='Int64' Name='faces' format='ascii' &
!      &RangeMin='0' RangeMax='", mesh%n_vert-1, "'>"
!
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write(fn, *) mesh%elem(i)%n_faces
!        do j=1, mesh%elem(i)%n_faces
!          id_face = mesh%elem(i)%face(j)
!          allocate(local_id_vert_no_ghost(mesh%face(id_face)%n_vert))
!          do k=1, mesh%face(id_face)%n_vert
!            local_id_vert_no_ghost(k) = id_vert_no_ghost(mesh%face(id_face)%vert(k))
!          end do
!          !write(fn, *) mesh%face(id_face)%n_vert, mesh%face(id_face)%vert(:) - 1
!          write(fn, *) mesh%face(id_face)%n_vert, local_id_vert_no_ghost - 1
!          deallocate(local_id_vert_no_ghost)
!        end do
!      end if
!    end do
!
!    write(fn, *) "</DataArray>"
!
!    s = 0
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        do j=1, mesh%elem(i)%n_faces
!          id_face = mesh%elem(i)%face(j)
!          s = s + mesh%face(id_face)%n_vert + 1
!        end do
!        s = s + 1
!      end if
!    end do
!
!    write(fn, *) "<DataArray type='Int64' Name='faceoffsets' format='ascii' &
!      &RangeMin='", 0, "' RangeMax='", s, "'>"
!
!    s = 0
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        do j=1, mesh%elem(i)%n_faces
!          id_face = mesh%elem(i)%face(j)
!          s = s + mesh%face(id_face)%n_vert + 1
!        end do
!        s = s + 1
!        write(fn, *) s
!      end if
!    end do
!
!    write(fn, *) "</DataArray>"
!    write(fn, *) "</Cells>"
!
!    !Point data
!    write(fn, *) "<PointData>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Nodal_Velocity' format='ascii' NumberOfComponents='3'>"
!    do i=1, mesh%n_vert
!      if( .not. mesh%vert(i)%is_ghost ) then
!        write (fn, *) vp(:, i)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Density' format='ascii' NumberOfComponents='3'>"
!    do i=1, mesh%n_vert
!      if( .not. mesh%vert(i)%is_ghost ) then
!        call compute_nodal_grad(mesh, i, sol, nodal_grad)
!        write (fn, *) nodal_grad(:, 1)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Velocity_X' format='ascii' NumberOfComponents='3'>"
!    do i=1, mesh%n_vert
!      if( .not. mesh%vert(i)%is_ghost ) then
!        call compute_nodal_grad(mesh, i, sol, nodal_grad)
!        write (fn, *) nodal_grad(:, 2)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Velocity_Y' format='ascii' NumberOfComponents='3'>"
!    do i=1, mesh%n_vert
!      if( .not. mesh%vert(i)%is_ghost ) then
!        call compute_nodal_grad(mesh, i, sol, nodal_grad)
!        write (fn, *) nodal_grad(:, 3)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Velocity_Z' format='ascii' NumberOfComponents='3'>"
!    do i=1, mesh%n_vert
!      if( .not. mesh%vert(i)%is_ghost ) then
!        call compute_nodal_grad(mesh, i, sol, nodal_grad)
!        write (fn, *) nodal_grad(:, 4)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Pressure' format='ascii' NumberOfComponents='3'>"
!    do i=1, mesh%n_vert
!      if( .not. mesh%vert(i)%is_ghost ) then
!        call compute_nodal_grad(mesh, i, sol, nodal_grad)
!        write (fn, *) nodal_grad(:, 5)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    call compute_cell_size(mesh, sol, grad, all_nodal_grad, cell_size)
!
!    write(fn, *) "<DataArray type='Float64' Name='cell_size' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_vert
!      if( .not. mesh%vert(i)%is_ghost ) then
!        write (fn, *) cell_size(i)/mesh%vert(i)%volume**(1.0_DOUBLE/3.0_DOUBLE)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "</PointData>"
!
!    !Cell data
!    write(fn, *) "<CellData>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Centroid' format='ascii' NumberOfComponents='3'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write(fn, *) mesh%elem(i)%coord
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Int32' Name='Tag' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write (fn, *) mesh%elem(i)%tag
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Int32' Name='MPI_Color' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write(fn, *) me
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Density' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write(fn, *) sol(1, i)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Momentum' format='ascii' NumberOfComponents='3'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write(fn, *) sol(2:4, i)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Energy' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write(fn, *) sol(5, i)/sol(1, i)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Velocity' format='ascii' NumberOfComponents='3'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write(fn, *) sol(2:4, i)/sol(1, i)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Pressure' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        w = conserv_to_primit(sol(:, i))
!        write(fn, *) w(5)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Internal_energy' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write (fn, *) sol(5, i)/sol(1, i) - 0.5_DOUBLE*norm2(sol(2:4, i)/sol(1, i))**2
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='H' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        w = conserv_to_primit(sol(:, i))
!        write (fn, *) sol(5, i)/sol(1, i) + w(5)/sol(1, i)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Mach' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        w = conserv_to_primit(sol(:, i))
!        write (fn, *) norm2(sol(2:4, i)/sol(1, i))/(sqrt(gamma*w(5)/w(1)))
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Temperature' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write (fn, *) (sol(5, i)/sol(1, i) - 0.5_DOUBLE*norm2(sol(2:4, i)/sol(1, i))**2)/Cv_p
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='mu' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        if( use_sutherland ) then
!          write (fn, *) mu(sol(:, i))
!        else
!          write (fn, *) mu_p
!        end if
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Ducros' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write (fn, *) omega_ducros(mesh, sol, i)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='r' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write (fn, *) norm2(mesh%elem(i)%coord)
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='rxy' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write (fn, *) norm2(mesh%elem(i)%coord(1:2))
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    write(fn, *) "<DataArray type='Float64' Name='Theta3' format='ascii' NumberOfComponents='1'>"
!    do i=1, mesh%n_elems
!      if( .not. mesh%elem(i)%is_ghost ) then
!        write (fn, *) atan2(-mesh%elem(i)%coord(2), -mesh%elem(i)%coord(1))
!      end if
!    end do
!    write(fn, *) "</DataArray>"
!
!    if( second_order ) then
!      write(fn, *) "<DataArray type='Float64' Name='Density_Grad' format='ascii' NumberOfComponents='3'>"
!      do i=1, mesh%n_elems
!        if( .not. mesh%elem(i)%is_ghost ) then
!          write (fn, *) grad(:, 1, i)
!        end if
!      end do
!      write(fn, *) "</DataArray>"
!
!      write(fn, *) "<DataArray type='Float64' Name='Velocity_X_Grad' format='ascii' NumberOfComponents='3'>"
!      do i=1, mesh%n_elems
!        if( .not. mesh%elem(i)%is_ghost ) then
!          write (fn, *) grad(:, 2, i)
!        end if
!      end do
!      write(fn, *) "</DataArray>"
!
!      write(fn, *) "<DataArray type='Float64' Name='Velocity_Y_Grad' format='ascii' NumberOfComponents='3'>"
!      do i=1, mesh%n_elems
!        if( .not. mesh%elem(i)%is_ghost ) then
!          write (fn, *) grad(:, 3, i)
!        end if
!      end do
!      write(fn, *) "</DataArray>"
!
!      write(fn, *) "<DataArray type='Float64' Name='Velocity_Z_Grad' format='ascii' NumberOfComponents='3'>"
!      do i=1, mesh%n_elems
!        if( .not. mesh%elem(i)%is_ghost ) then
!          write (fn, *) grad(:, 4, i)
!        end if
!      end do
!      write(fn, *) "</DataArray>"
!
!      write(fn, *) "<DataArray type='Float64' Name='Pressure_Grad' format='ascii' NumberOfComponents='3'>"
!      do i=1, mesh%n_elems
!        if( .not. mesh%elem(i)%is_ghost ) then
!          write (fn, *) grad(:, 5, i)
!        end if
!      end do
!      write(fn, *) "</DataArray>"
!    end if
!
!    write(fn, *) "</CellData>"
!
!    write(fn, *) "</Piece>"
!    write(fn, *) "</UnstructuredGrid>"
!    write(fn, *) "</VTKFile>"
!
!    close(fn)
!  end subroutine write_sol_vtu
!
!  subroutine write_sol_meta_pvtu(filename, iaff_char)
!    use mpi
!    use mpi_module
!    use global_data_module, only: second_order
!    implicit none
!
!    character(len=*), intent(in) :: filename, iaff_char
!
!    integer(kind=ENTIER) :: fn, i, num_procs, me, mpi_ierr
!    character(len=255) :: i_char
!
!    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
!    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)
!
!    if( me == 0 ) then
!      open(newunit=fn, file=trim(adjustl(filename)))
!
!      write(fn, *) "<VTKFile type='PUnstructuredGrid' version='0.1' byte_order='LittleEndian' header_type='UInt64'>"
!      write(fn, *) "<PUnstructuredGrid GhostLevel='0'>"
!      write(fn, *) "<PPoints>"
!      write(fn, *) "<PDataArray type='Float64' Name='Points' NumberOfComponents='3'/>"
!      write(fn, *) "</PPoints>"
!      write(fn, *) "<PCells>"
!      write(fn, *) "<PDataArray type='Int64' Name='connectivity' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Int64' Name='offsets'      NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='UInt8' Name='types'        NumberOfComponents='1'/>"
!      write(fn, *) "</PCells>"
!      write(fn, *) "<PPointData>"
!      write(fn, *) "<DataArray type='Float64' Name='Nodal_Velocity' format='ascii' NumberOfComponents='3'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Density' NumberOfComponents='3'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Velocity_X' NumberOfComponents='3'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Velocity_Y' NumberOfComponents='3'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Velocity_Z' NumberOfComponents='3'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Pressure' NumberOfComponents='3'/>"
!      write(fn, *) "<DataArray type='Float64' Name='cell_size' format='ascii' NumberOfComponents='1'/>"
!      write(fn, *) "</PPointData>"
!      write(fn, *) "<PCellData>"
!      write(fn, *) "<PDataArray type='Float64' Name='Centroid' NumberOfComponents='3'/>"
!      write(fn, *) "<PDataArray type='Int32' Name='Tag' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Int32' Name='MPI_Color' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Density' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Momentum' NumberOfComponents='3'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Energy' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Velocity' NumberOfComponents='3'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Pressure' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Internal_energy' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='H' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Mach' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Temperature' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='mu' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Ducros' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='r' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='rxy' NumberOfComponents='1'/>"
!      write(fn, *) "<PDataArray type='Float64' Name='Theta3' NumberOfComponents='1'/>"
!      if( second_order ) then
!        write(fn, *) "<PDataArray type='Float64' Name='Density_Grad' NumberOfComponents='3'/>"
!        write(fn, *) "<PDataArray type='Float64' Name='Velocity_X_Grad' NumberOfComponents='3'/>"
!        write(fn, *) "<PDataArray type='Float64' Name='Velocity_Y_Grad' NumberOfComponents='3'/>"
!        write(fn, *) "<PDataArray type='Float64' Name='Velocity_Z_Grad' NumberOfComponents='3'/>"
!        write(fn, *) "<PDataArray type='Float64' Name='Pressure_Grad' NumberOfComponents='3'/>"
!      end if
!      write(fn, *) "</PCellData>"
!
!      do i=0, num_procs-1
!        write(i_char, *) i
!        write(fn, *) "<Piece Source='"//trim(adjustl(i_char))//"_output_"//&
!          trim(adjustl(iaff_char))//".vtu'/>"
!      end do
!
!      write(fn, *) "</PUnstructuredGrid>"
!      write(fn, *) "</VTKFile>"
!
!      close(fn)
!    end if
!
!    call mpi_barrier(mpi_comm_world, mpi_ierr)
!  end subroutine write_sol_meta_pvtu

end module io_module
