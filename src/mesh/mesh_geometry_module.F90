module mesh_geometry_module
  use precision_module
  use mesh_module
  use mpi_module
  implicit none

  private

  public :: compute_geometry_mesh
  public :: map_cylinder
  public :: odd_even
  public :: project_sol
  public :: project_sol_box
  public :: mpi_project_sol_box
contains

  subroutine compute_geometry_mesh(mesh, use_sub_entities, b2d)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    logical, intent(in) :: use_sub_entities, b2d

    integer(kind=ENTIER) :: i

    do i=1, mesh%n_faces
      call compute_face_centroid(mesh, i)
    end do

    do i=1, mesh%n_elems
      call compute_elem_centroid(mesh, i)
    end do

    call build_face_centroid_norm(mesh)
    if (use_sub_entities) call build_sub_face_area_norm_and_sub_elem_vol(mesh)
    if (use_sub_entities) call correct_sub_face_area(mesh)
    !if (use_sub_entities) call check_mesh(mesh)
    call find_if_vert_is_bound(mesh, b2d)
    call compute_vert_volume(mesh)
  end subroutine compute_geometry_mesh

  subroutine compute_face_centroid(mesh, id_face)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    integer(kind=ENTIER), intent(in) :: id_face

    integer(kind=ENTIER) :: idv1, idv2, k
    real(kind=DOUBLE), dimension(3) :: c1, c2, cgeom
    real(kind=DOUBLE), dimension(3) :: pdv, cp

    cgeom = 0.0_DOUBLE
    do k = 1, mesh%face(id_face)%n_vert
      idv1 = mesh%face(id_face)%vert(k)
      c1 = mesh%vert(idv1)%coord
      cgeom = cgeom + c1
    end do
    cgeom = cgeom / mesh%face(id_face)%n_vert

    pdv = 0.0_DOUBLE
    do k = 1, mesh%face(id_face)%n_vert
      idv1 = mesh%face(id_face)%vert(k)
      idv2 = mesh%face(id_face)%vert(1+mod((k-1)+1,mesh%face(id_face)%n_vert))
      c1 = mesh%vert(idv1)%coord
      c2 = mesh%vert(idv2)%coord
      pdv = pdv + cross_product(c1-cgeom,c2-cgeom)
    end do
    mesh%face(id_face)%area = 0.5_DOUBLE*norm2(pdv)
    mesh%face(id_face)%norm = pdv/norm2(pdv)

    mesh%face(id_face)%coord = 0.0_DOUBLE
    do k = 1, mesh%face(id_face)%n_vert
      idv1 = mesh%face(id_face)%vert(k)
      idv2 = mesh%face(id_face)%vert(1+mod((k-1)+1,mesh%face(id_face)%n_vert))
      c1 = mesh%vert(idv1)%coord
      c2 = mesh%vert(idv2)%coord
      mesh%face(id_face)%coord = mesh%face(id_face)%coord &
        + (c1+c2+cgeom) * dot_product(cross_product(c1, c2), &
        mesh%face(id_face)%norm)
    end do
    mesh%face(id_face)%coord = mesh%face(id_face)%coord &
      / (6*mesh%face(id_face)%area)
  end subroutine compute_face_centroid

  subroutine compute_elem_centroid(mesh, i)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    integer(kind=ENTIER), intent(in) :: i

    integer(kind=ENTIER) :: j, k
    integer(kind=ENTIER) :: id_face
    integer(kind=ENTIER) :: idv1, idv2
    real(kind=DOUBLE), dimension(3) :: c1, c2
    real(kind=DOUBLE) :: vol, sgn

    mesh%elem(i)%volume = 0.0_DOUBLE
    mesh%elem(i)%coord = 0.0_DOUBLE
    do j=1, mesh%elem(i)%n_faces
      id_face = mesh%elem(i)%face(j)
      if( mesh%face(id_face)%right_neigh == i ) then
        sgn = -1.0_DOUBLE
      else
        sgn = 1.0_DOUBLE
      end if
      do k = 1, mesh%face(id_face)%n_vert
        idv1 = mesh%face(id_face)%vert(k)
        idv2 = mesh%face(id_face)%vert(1+mod((k-1)+1,mesh%face(id_face)%n_vert))
        c1 = mesh%vert(idv1)%coord
        c2 = mesh%vert(idv2)%coord
        vol = sgn*dot_product(mesh%face(id_face)%coord, &
          cross_product(c1, c2))/6.0_DOUBLE
        mesh%elem(i)%volume = mesh%elem(i)%volume + vol
        mesh%elem(i)%coord = mesh%elem(i)%coord &
          + vol*(c1+c2+mesh%face(id_face)%coord)/4.0_DOUBLE
      end do
    end do
    mesh%elem(i)%coord = mesh%elem(i)%coord/mesh%elem(i)%volume
  end subroutine compute_elem_centroid

  !subroutine compute_elem_face_centroid(mesh)
  !  implicit none

  !  type(mesh_type), intent(inout) :: mesh

  !  integer(kind=ENTIER) :: i, j, k
  !  integer(kind=ENTIER) :: id_vert, id_face, idv1, idv2
  !  real(kind=DOUBLE) :: vol, vol_tot, area, area_tot
  !  real(kind=DOUBLE), dimension(3) :: cgeom, cfgeom
  !  real(kind=DOUBLE), dimension(3) :: centroid_face, centroid_elem
  !  real(kind=DOUBLE), dimension(3) :: v1, v2, v3, pdv, c1, c2

  !  logical :: orient_chosen
  !  real(kind=DOUBLE), dimension(3) :: orient_vect

  !  !$OMP PARALLEL DO PRIVATE(j, k, id_vert, id_face, idv1, idv2, vol, vol_tot, area, area_tot,&
  !  !$OMP& cgeom, cfgeom, centroid_face, centroid_elem, v1, v2, v3, pdv, c1, c2, orient_chosen, orient_vect)
  !  do i = 1, mesh%n_elems

  !    centroid_elem = 0.0_DOUBLE
  !    cgeom = 0.0_DOUBLE
  !    do j = 1, mesh%elem(i)%n_vert
  !      id_vert = mesh%elem(i)%vert(j)
  !      cgeom = cgeom + mesh%vert(id_vert)%coord
  !    end do
  !    cgeom = cgeom/mesh%elem(i)%n_vert

  !    vol_tot = 0.0_DOUBLE
  !    do j = 1, mesh%elem(i)%n_faces
  !      id_face = mesh%elem(i)%face(j)
  !      orient_chosen = .false.

  !      cfgeom = 0.0_DOUBLE
  !      do k = 1, mesh%face(id_face)%n_vert
  !        id_vert = mesh%face(id_face)%vert(k)
  !        cfgeom = cfgeom + mesh%vert(id_vert)%coord
  !      end do
  !      cfgeom = cfgeom/mesh%face(id_face)%n_vert

  !      centroid_face = 0.0_DOUBLE
  !      area_tot = 0.0_DOUBLE
  !      vol = 0.0_DOUBLE
  !      do k = 1, mesh%face(id_face)%n_vert
  !        idv1 = mesh%face(id_face)%vert(k)
  !        idv2 = mesh%face(id_face)%vert(1 + mod((k - 1) + 1, mesh%face(id_face)%n_vert))

  !        c1 = mesh%vert(idv1)%coord
  !        c2 = mesh%vert(idv2)%coord

  !        v1 = c1 - cfgeom
  !        v2 = c2 - cfgeom
  !        v3 = cgeom - cfgeom

  !        pdv(1) = v1(2)*v2(3) - v1(3)*v2(2)
  !        pdv(2) = -(v1(1)*v2(3) - v1(3)*v2(1))
  !        pdv(3) = v1(1)*v2(2) - v1(2)*v2(1)

  !        if( .not. orient_chosen ) then
  !          orient_vect = pdv/norm2(pdv)
  !          orient_chosen = .true.
  !        end if

  !        area = 0.5_DOUBLE*norm2(pdv)*sign(1.0_DOUBLE, dot_product(pdv, orient_vect))
  !        area_tot = area_tot + area

  !        vol = vol + (1.0_DOUBLE/6.0_DOUBLE)*abs(pdv(1)*v3(1) &
  !          + pdv(2)*v3(2) + pdv(3)*v3(3))*sign(1.0_DOUBLE, dot_product(pdv, orient_vect))

  !        centroid_face = centroid_face + area*(1.0_DOUBLE/3.0_DOUBLE)*(c1 + c2 + cfgeom)
  !      end do
  !      centroid_face = centroid_face/area_tot
  !      mesh%face(id_face)%coord = centroid_face

  !      centroid_elem = centroid_elem + &
  !        vol*(centroid_face + 0.25_DOUBLE*(cgeom - centroid_face))
  !      vol_tot = vol_tot + vol
  !    end do

  !    centroid_elem = centroid_elem/vol_tot
  !    mesh%elem(i)%coord = centroid_elem
  !  end do
  !end subroutine compute_elem_face_centroid

  subroutine compute_vert_volume(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, j, id_sub_elem

    !$OMP PARALLEL DO PRIVATE(i, j, id_sub_elem)
    do i=1, mesh%n_vert
      mesh%vert(i)%volume = 0.0_DOUBLE
      do j=1, mesh%vert(i)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(i)%sub_elem_neigh(j)
        mesh%vert(i)%volume = mesh%vert(i)%volume &
          + mesh%sub_elem(id_sub_elem)%volume
      end do
    end do
  end subroutine compute_vert_volume

  subroutine build_face_centroid_norm(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, j, jp, id_left_elem, id_right_elem
    integer(kind=ENTIER) :: idv1, idv2
    real(kind=DOUBLE) :: area, tot_area, vol_left, vol_right
    real(kind=DOUBLE), dimension(3) :: centroid_face, cf, cgeom_face
    real(kind=DOUBLE), dimension(3) :: v1, v2, norm, n1, vl, vr

    !$OMP PARALLEL DO
    do i = 1, mesh%n_elems
      mesh%elem(i)%volume = 0.0_DOUBLE
    end do

    !$OMP PARALLEL DO PRIVATE(j, jp, id_left_elem, id_right_elem, &
    !$OMP& idv1, idv2, &
    !$OMP& area, tot_area, vl, vr, &
    !$OMP& cf, cgeom_face, centroid_face, v1, v2, norm, n1)
    do i = 1, mesh%n_faces
      id_left_elem = mesh%face(i)%left_neigh
      id_right_elem = mesh%face(i)%right_neigh
      cgeom_face = 0.0_DOUBLE
      do j = 1, mesh%face(i)%n_vert
        idv1 = mesh%face(i)%vert(j)
        cgeom_face = cgeom_face + mesh%vert(idv1)%coord
      end do
      cgeom_face = cgeom_face/mesh%face(i)%n_vert

      centroid_face = 0.0_DOUBLE
      tot_area = 0.0_DOUBLE
      norm = 0.0_DOUBLE
      do j = 1, mesh%face(i)%n_vert
        jp = 1 + mod((j - 1) + 1 + mesh%face(i)%n_vert, mesh%face(i)%n_vert)

        idv1 = mesh%face(i)%vert(j)
        idv2 = mesh%face(i)%vert(jp)

        cf = 1.0_DOUBLE/3.0_DOUBLE*(mesh%vert(idv1)%coord + mesh%vert(idv2)%coord + cgeom_face)
        v1 = mesh%vert(idv1)%coord - cgeom_face
        v2 = mesh%vert(idv2)%coord - cgeom_face

        n1 = cross_product(v1, v2)
        area = 0.5_DOUBLE*norm2(n1)
        n1 = 0.5_DOUBLE*n1/area

        !if (dot_product(n1, cgeom_face - mesh%elem(id_left_elem)%coord) < 0.0_DOUBLE) then
        !  n1 = -n1
        !end if

        norm = norm + area*n1
        centroid_face = centroid_face + area*cf
        tot_area = tot_area + area

        vl = mesh%elem(id_left_elem)%coord - cgeom_face
        vol_left = (1.0_DOUBLE/6.0_DOUBLE)*dot_product(cross_product(v1, v2), vl)
        !$OMP ATOMIC UPDATE
        mesh%elem(id_left_elem)%volume = mesh%elem(id_left_elem)%volume + vol_left

        if (id_right_elem > 0) then
          vr = mesh%elem(id_right_elem)%coord - cgeom_face
          vol_right = (1.0_DOUBLE/6.0_DOUBLE)*dot_product(cross_product(v1, v2), vr)
          !$OMP ATOMIC UPDATE
          mesh%elem(id_right_elem)%volume = mesh%elem(id_right_elem)%volume + vol_right
        end if

      end do
      mesh%face(i)%area = norm2(norm)
      mesh%face(i)%coord = centroid_face/tot_area
      mesh%face(i)%norm = norm/norm2(norm)
    end do
  end subroutine build_face_centroid_norm

  subroutine build_sub_face_area_norm_and_sub_elem_vol(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, j1, j3, j
    integer(kind=ENTIER) :: id_vert, id_face
    integer(kind=ENTIER) :: id1, id2, id3, id_left_elem, id_right_elem, id_sub_elem
    integer(kind=ENTIER) :: id_left_sub_elem, id_right_sub_elem
    real(kind=DOUBLE) :: area1, area2, vol_left, vol_right
    real(kind=DOUBLE), dimension(3) :: v1, v2, v3, vl, vr
    real(kind=DOUBLE), dimension(3) :: p1, p2, p3
    real(kind=DOUBLE), dimension(3) :: n1, n2, norm

    !$OMP PARALLEL DO
    do i = 1, mesh%n_sub_elems
      mesh%sub_elem(i)%volume = 0.0_DOUBLE
    end do

    !$OMP PARALLEL DO PRIVATE(j1, j3, j,&
    !$OMP& id_vert, id_face,&
    !$OMP& id1, id2, id3, id_left_elem, id_right_elem,&
    !$OMP& id_left_sub_elem, id_right_sub_elem,&
    !$OMP& area1, area2, vol_left, vol_right,&
    !$OMP& v1, v2, v3, vl, vr, p1, p2, p3, n1, n2, norm)
    do i = 1, mesh%n_sub_faces
      id_face = mesh%sub_face(i)%mesh_face
      id_vert = mesh%sub_face(i)%mesh_vert

      id_left_sub_elem = mesh%sub_face(i)%left_sub_elem_neigh
      id_right_sub_elem = mesh%sub_face(i)%right_sub_elem_neigh

      id_left_elem = mesh%face(id_face)%left_neigh
      id_right_elem = mesh%face(id_face)%right_neigh

      do j = 1, mesh%face(id_face)%n_vert
        id2 = mesh%face(id_face)%vert(j)

        if (id2 == id_vert) then
          j1 = 1 + mod((j - 1) - 1 + mesh%face(id_face)%n_vert, mesh%face(id_face)%n_vert)
          j3 = 1 + mod((j - 1) + 1 + mesh%face(id_face)%n_vert, mesh%face(id_face)%n_vert)

          id1 = mesh%face(id_face)%vert(j1)
          id3 = mesh%face(id_face)%vert(j3)

          p1 = 0.5_DOUBLE*(mesh%vert(id1)%coord + mesh%vert(id2)%coord)
          p2 = mesh%vert(id2)%coord
          p3 = 0.5_DOUBLE*(mesh%vert(id2)%coord + mesh%vert(id3)%coord)

          v1 = p1 - mesh%face(id_face)%coord
          v2 = p2 - mesh%face(id_face)%coord
          v3 = p3 - mesh%face(id_face)%coord

          n1 = cross_product(v1, v2)
          area1 = 0.5_DOUBLE*norm2(n1)
          n1 = n1/norm2(n1)

          n2 = cross_product(v2, v3)
          area2 = 0.5_DOUBLE*norm2(n2)
          n2 = n2/norm2(n2)

          norm = (area1*n1 + area2*n2)
          if (dot_product(norm, mesh%face(id_face)%coord - mesh%elem(id_left_elem)%coord) < 0.0_DOUBLE) then
            norm = -norm
          end if

          mesh%sub_face(i)%area = norm2(norm)
          ! mesh%sub_face(i)%norm = norm/norm2(norm)

          vl = mesh%elem(id_left_elem)%coord - p2
          vol_left = (1.0_DOUBLE/6.0_DOUBLE)*abs(dot_product(cross_product(v1, v2), vl)) &
            + (1.0_DOUBLE/6.0_DOUBLE)*abs(dot_product(cross_product(v2, v3), vl))

          !$OMP ATOMIC UPDATE
          mesh%sub_elem(id_left_sub_elem)%volume = mesh%sub_elem(id_left_sub_elem)%volume + vol_left

          if (id_right_elem > 0) then
            vr = mesh%elem(id_right_elem)%coord - p2
            vol_right = (1.0_DOUBLE/6.0_DOUBLE)*abs(dot_product(cross_product(v1, v2), vr)) &
              + (1.0_DOUBLE/6.0_DOUBLE)*abs(dot_product(cross_product(v2, v3), vr))
            !$OMP ATOMIC UPDATE
            mesh%sub_elem(id_right_sub_elem)%volume = mesh%sub_elem(id_right_sub_elem)%volume + vol_right
          end if

          exit
        end if
      end do
    end do

    !$OMP PARALLEL DO PRIVATE(id_sub_elem)
    do i = 1, mesh%n_elems
      mesh%elem(i)%volume = 0.0_DOUBLE
      do j = 1, mesh%elem(i)%n_sub_elems
        id_sub_elem = mesh%elem(i)%sub_elem(j)
        mesh%elem(i)%volume = mesh%elem(i)%volume + mesh%sub_elem(id_sub_elem)%volume
      end do
    end do
  end subroutine build_sub_face_area_norm_and_sub_elem_vol

  subroutine correct_sub_face_area(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i, j, id_sub_face
    real(kind=DOUBLE) :: tot_area

    !$OMP PARALLEL DO PRIVATE(j, id_sub_face, tot_area)
    do i = 1, mesh%n_faces
      tot_area = 0.0_DOUBLE
      do j = 1, mesh%face(i)%n_sub_faces
        id_sub_face = mesh%face(i)%sub_face(j)
        tot_area = tot_area + mesh%sub_face(id_sub_face)%area
      end do

      do j = 1, mesh%face(i)%n_sub_faces
        id_sub_face = mesh%face(i)%sub_face(j)
        mesh%sub_face(id_sub_face)%area = mesh%face(i)%area* &
          (mesh%sub_face(id_sub_face)%area/tot_area)
      end do
    end do
  end subroutine correct_sub_face_area

  subroutine check_mesh(mesh)
    implicit none

    type(mesh_type), intent(in) :: mesh

    logical :: mesh_is_ok

    integer(kind=ENTIER) :: i, j, id_sub_elem, k, id_sub_face, id_face
    real(kind=DOUBLE) :: test_circulation, tot_area
    real(kind=DOUBLE), dimension(3) :: test_norm, test_tot_norm

    mesh_is_ok = .TRUE.
    test_circulation = 0.
    test_tot_norm(:) = 0.

    do i = 1, mesh%n_elems

      test_norm(:) = 0.
      do j = 1, mesh%elem(i)%n_faces
        id_face = mesh%elem(i)%face(j)
        if (mesh%face(id_face)%left_neigh == i) then
          test_norm = test_norm + mesh%face(id_face)%norm*mesh%face(id_face)%area
          print*, mesh%face(id_face)%norm, mesh%face(id_face)%area
        else
          test_norm = test_norm - mesh%face(id_face)%norm*mesh%face(id_face)%area
          print*, -mesh%face(id_face)%norm, mesh%face(id_face)%area
        end if
      end do
      test_tot_norm(:) = test_tot_norm(:) + test_norm(:)
      test_circulation = test_circulation + norm2(test_norm)
      if (norm2(test_norm) > 1.e-8_DOUBLE) then
        print *, achar(27)//"[31m[-] norm face"//achar(27)//"[0m", test_norm, i, mesh%elem(i)%elem_kind
        mesh_is_ok = .FALSE.
      end if

      if( .not. mesh_is_ok ) error stop

      test_norm = 0.0_DOUBLE
      do j = 1, mesh%elem(i)%n_sub_elems
        id_sub_elem = mesh%elem(i)%sub_elem(j)
        do k = 1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          id_sub_face = mesh%sub_elem(id_sub_elem)%sub_face(k)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          if (mesh%sub_face(id_sub_face)%left_sub_elem_neigh == id_sub_elem) then
            test_norm = test_norm &
              ! + mesh%sub_face(id_sub_face)%area*mesh%sub_face(id_sub_face)%norm
            + mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
          else
            test_norm = test_norm &
              ! - mesh%sub_face(id_sub_face)%area*mesh%sub_face(id_sub_face)%norm
            - mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
          end if
        end do
      end do

      if (norm2(test_norm) > 1e-8_DOUBLE) then
        mesh_is_ok = .FALSE.
        print *, achar(27)//"[31m[-] norm subface"//achar(27)//"[0m", norm2(test_norm), i, mesh%elem(i)%elem_kind
      end if
    end do

    !if (.not. mesh_is_ok) error stop

    !do i = 1, mesh%n_faces
    !  tot_area = 0.0_DOUBLE
    !  do j = 1, mesh%face(i)%n_sub_faces
    !    id_sub_face = mesh%face(i)%sub_face(j)
    !    tot_area = tot_area + mesh%sub_face(id_sub_face)%area
    !  end do

    !  if (abs(mesh%face(i)%area - tot_area)/tot_area > 1e-8_DOUBLE) then
    !    print *, achar(27)//"[31m[-] Bad area of face computed from subfaces "//achar(27)//"[0m", &
    !      mesh%face(i)%area, tot_area
    !    mesh_is_ok = .FALSE.
    !  end if
    !end do

    !if (mesh_is_ok) then
    !  print *, achar(27)//"[32m[+] Circulation test passed"//achar(27)//"[0m", &
    !    sqrt(test_tot_norm(1)**2 + test_tot_norm(2)**2 + test_tot_norm(3)**2), &
    !    achar(27)//"[0m"
    !else
    !  print *, achar(27)//"[31m[-] Issue with the mesh, circulation test failed ! "//achar(27)//"[0m"
    !  print *, achar(27)//"[31m[-] Circulation over every element's face : "//achar(27)//"[0m", &
    !    sqrt(test_tot_norm(1)**2 + test_tot_norm(2)**2 + test_tot_norm(3)**2)
    !  error stop
    !end if

    print*,"OKOKOKOK"
    error stop
  end subroutine check_mesh

  pure subroutine map_cylinder(mesh, dimen)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    integer, intent(in) :: dimen

    integer(kind=ENTIER) :: i

    do i = 1, mesh%n_vert
      mesh%vert(i)%coord = cylinder_mapping(mesh%vert(i)%coord, dimen)
    end do
  end subroutine map_cylinder

  subroutine odd_even(mesh)
    implicit none

    type(mesh_type), intent(inout) :: mesh

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE), parameter :: PI = 4.0_DOUBLE*datan(1.0_DOUBLE)
    real(kind=DOUBLE) :: l, theta
    ! l=1e-9
    l=1e-9_DOUBLE

    do i = 1, mesh%n_vert
      if( abs(mesh%vert(i)%coord(2)) < 1e-2_DOUBLE .and. &
        abs(mesh%vert(i)%coord(3)) < 1e-2_DOUBLE) then
        theta = mesh%vert(i)%coord(1)*pi/2.0_DOUBLE
        ! print*,mesh%vert(i)%coord
        mesh%vert(i)%coord(2) = mesh%vert(i)%coord(2) + l*cos(theta)
        mesh%vert(i)%coord(3) = mesh%vert(i)%coord(3) + l*sin(theta)
        ! print*, l*cos(theta), l*sin(theta)
        ! print*,mesh%vert(i)%coord
        ! print*,"-------"
      end if
    end do
  end subroutine odd_even

  pure function cylinder_mapping(coord, dimen) result(new_coord)
    implicit none

    integer, intent(in) :: dimen
    real(kind=DOUBLE), dimension(3), intent(in) :: coord
    real(kind=DOUBLE), dimension(3) :: new_coord

    real(kind=DOUBLE) :: r, new_r, theta, d
    real(kind=DOUBLE) :: r1, new_r1, r2, new_r2

    if (dimen == 2) then
      r = norm2(coord(:2))
      theta = atan2(coord(2), -coord(1))

      r1 = 1.0
      new_r1 = r1

      r2 = 1.3851721918740656_DOUBLE
      new_r2 = r2 &
        + 0.22082080765324227_DOUBLE*theta**2 &
        + 0.07247276314520341_DOUBLE*theta**4 &
        - 0.01622444200319442_DOUBLE*theta**6 &
        + 0.00780545523147566_DOUBLE*theta**8

      new_r = (r - r1)*(new_r2 - new_r1)/(r2 - r1) + new_r1

      new_coord(1) = -new_r*cos(theta)
      new_coord(2) = new_r*sin(theta)
      new_coord(3) = coord(3)
    else if (dimen == 21) then !For cylinder with radius 0.5
      r = norm2(2.0_DOUBLE*coord(:2))
      theta = atan2(2.0_DOUBLE*coord(2), -2.0_DOUBLE*coord(1))

      r1 = 1.0
      new_r1 = r1

      r2 = 1.3851721918740656_DOUBLE
      new_r2 = r2 &
        + 0.22082080765324227_DOUBLE*theta**2 &
        + 0.07247276314520341_DOUBLE*theta**4 &
        - 0.01622444200319442_DOUBLE*theta**6 &
        + 0.00780545523147566_DOUBLE*theta**8

      new_r = (r - r1)*(new_r2 - new_r1)/(r2 - r1) + new_r1

      new_coord(1) = -new_r*cos(theta)/2.0_DOUBLE
      new_coord(2) = new_r*sin(theta)/2.0_DOUBLE
      new_coord(3) = coord(3)
    else if (dimen == 3) then
      r = norm2(coord)
      d = norm2(coord(2:))
      theta = atan2(d, -coord(1))

      r1 = 1.0
      new_r1 = r1

      r2 = 1.3851721918740656_DOUBLE
      new_r2 = r2 &
        + 0.22082080765324227_DOUBLE*theta**2 &
        + 0.07247276314520341_DOUBLE*theta**4 &
        - 0.01622444200319442_DOUBLE*theta**6 &
        + 0.00780545523147566_DOUBLE*theta**8

      new_r = (r - r1)*(new_r2 - new_r1)/(r2 - r1) + new_r1
      new_coord = coord + (new_r - r)*coord
    end if

  end function cylinder_mapping

  pure function tensor_product(a, b)
    implicit none

    real(kind=DOUBLE), intent(in) :: a(:), b(:)
    real(kind=DOUBLE), dimension(size(a), size(a)) :: tensor_product

    integer(kind=ENTIER) :: i, j

    do i = 1, size(a)
      do j = 1, size(a)
        tensor_product(i, j) = a(i)*b(j)
      end do
    end do
  end function tensor_product

  pure function cross_product(a, b)
    implicit none

    real(kind=DOUBLE), intent(in) :: a(3), b(3)
    real(kind=DOUBLE), dimension(3) :: cross_product

    cross_product(1) =   a(2)*b(3) - a(3)*b(2)
    cross_product(2) = -(a(1)*b(3) - a(3)*b(1))
    cross_product(3) =   a(1)*b(2) - a(2)*b(1)
  end function cross_product

  subroutine project_sol_box(sol_size, mesh1, sol1, mesh2, sol2)
    implicit none

    integer(kind=ENTIER) :: sol_size
    type(mesh_type), intent(in) :: mesh1, mesh2
    real(kind=DOUBLE), dimension(sol_size, mesh1%n_elems), intent(in) :: sol1
    real(kind=DOUBLE), dimension(sol_size, mesh2%n_elems), intent(out) :: sol2

    type :: box_type
      integer(kind=ENTIER), dimension(:), allocatable :: elem
    end type box_type

    real(kind=DOUBLE), parameter :: r_per_box=3.0_DOUBLE

    integer(kind=ENTIER) :: i, nx, ny, nz, ix, iy, iz, iclosest, j
    integer(kind=ENTIER) :: kx, ky, kz, idbx, idby, idbz
    real(kind=DOUBLE) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(kind=DOUBLE) :: r, tot_vol, dx, dy, dz, d, d1
    integer(kind=ENTIER), dimension(:,:,:), allocatable :: n_box
    type(box_type), dimension(:,:,:), allocatable :: box

    xmin = mesh1%elem(1)%coord(1)
    xmax = mesh1%elem(1)%coord(1)
    ymin = mesh1%elem(1)%coord(2)
    ymax = mesh1%elem(1)%coord(2)
    zmin = mesh1%elem(1)%coord(3)
    zmax = mesh1%elem(1)%coord(3)

    tot_vol = 0.0_DOUBLE
    do i = 1, mesh1%n_elems
      if(mesh1%elem(i)%coord(1) < xmin) xmin = mesh1%elem(i)%coord(1)
      if(mesh1%elem(i)%coord(1) > xmax) xmax = mesh1%elem(i)%coord(1)
      if(mesh1%elem(i)%coord(2) < ymin) ymin = mesh1%elem(i)%coord(2)
      if(mesh1%elem(i)%coord(2) > ymax) ymax = mesh1%elem(i)%coord(2)
      if(mesh1%elem(i)%coord(3) < zmin) zmin = mesh1%elem(i)%coord(3)
      if(mesh1%elem(i)%coord(3) > zmax) zmax = mesh1%elem(i)%coord(3)
      tot_vol = tot_vol + mesh1%elem(i)%volume
    end do

    r = (tot_vol/mesh1%n_elems)**(1.0_DOUBLE/3.0_DOUBLE)

    xmin = xmin - 1e-1_DOUBLE*r
    xmax = xmax + 1e-1_DOUBLE*r
    ymin = ymin - 1e-1_DOUBLE*r
    ymax = ymax + 1e-1_DOUBLE*r
    zmin = zmin - 1e-1_DOUBLE*r
    zmax = zmax + 1e-1_DOUBLE*r

    if( xmax - xmin < r ) then
      dx = xmax - xmin
    else
      dx = r_per_box*r
    end if
    nx = ceiling( (xmax - xmin)/dx )

    if( ymax - ymin < r ) then
      dy = ymax - ymin
    else
      dy = r_per_box*r
    end if
    ny = ceiling( (ymax - ymin)/dy )

    if( zmax - zmin < r ) then
      dz = zmax - zmin
    else
      dz = r_per_box*r
    end if
    nz = ceiling( (zmax - zmin)/dz )

    allocate(n_box(nx,ny,nz))
    n_box = 0
    do i=1,mesh1%n_elems
      ix = ceiling( (mesh1%elem(i)%coord(1) - xmin)/dx )
      iy = ceiling( (mesh1%elem(i)%coord(2) - ymin)/dy )
      iz = ceiling( (mesh1%elem(i)%coord(3) - zmin)/dz )
      n_box(ix,iy,iz) = n_box(ix,iy,iz) + 1
    end do

    allocate(box(nx,ny,nz))
    do ix=1,nx
      do iy=1,ny
        do iz=1,nz
          if(n_box(ix,iy,iz) > 0) then
            allocate(box(ix,iy,iz)%elem(n_box(ix,iy,iz)))
            box(ix,iy,iz)%elem = 0
          end if
        end do
      end do
    end do

    n_box = 0
    do i=1,mesh1%n_elems
      ix = ceiling( (mesh1%elem(i)%coord(1) - xmin)/dx )
      iy = ceiling( (mesh1%elem(i)%coord(2) - ymin)/dy )
      iz = ceiling( (mesh1%elem(i)%coord(3) - zmin)/dz )
      n_box(ix,iy,iz) = n_box(ix,iy,iz) + 1
      box(ix,iy,iz)%elem(n_box(ix,iy,iz)) = i
    enddo

    do i=1,mesh2%n_elems
      ix = ceiling( (mesh2%elem(i)%coord(1) - xmin)/dx )
      iy = ceiling( (mesh2%elem(i)%coord(2) - ymin)/dy )
      iz = ceiling( (mesh2%elem(i)%coord(3) - zmin)/dz )

      d = 10.0_DOUBLE*r
      iclosest = 0
      do kx=1,3
        do ky=1,3
          do kz=1,3
            idbx = min(max(1,ix+kx-2),nx)
            idby = min(max(1,iy+ky-2),ny)
            idbz = min(max(1,iz+kz-2),nz)
            do j=1,n_box(idbx, idby, idbz)
              d1 = norm2(mesh2%elem(i)%coord&
                -mesh1%elem(box(idbx,idby,idbz)%elem(j))%coord)
              if( d1 < d ) then
                iclosest = box(idbx,idby,idbz)%elem(j)
                d = d1
              end if
            end do
          end do
        end do
      end do

      if(iclosest == 0) then
        print*, "Issue to find closet elems in box"
        error stop
      else
        sol2(:, i) = sol1(:, iclosest)
      end if
    end do
  end subroutine project_sol_box

  subroutine project_sol(sol_size, mesh1, sol1, mesh2, sol2)
    implicit none

    integer(kind=ENTIER) :: sol_size
    type(mesh_type), intent(in) :: mesh1, mesh2
    real(kind=DOUBLE), dimension(sol_size, mesh1%n_elems), intent(in) :: sol1
    real(kind=DOUBLE), dimension(sol_size, mesh2%n_elems), intent(out) :: sol2

    integer(kind=ENTIER) :: i1, i2, iclosest, itest, k

    i1 = 1
    do i2 = 1, mesh2%n_elems
      !Find i1 in mesh1 closest to i2 in mesh2
      do while (.TRUE.)
        iclosest = i1
        do k = 1, mesh1%elem(i1)%n_neigh_by_vert
          itest = mesh1%elem(i1)%neigh_by_vert(k)
          if (norm2(mesh1%elem(itest)%coord - mesh2%elem(i2)%coord) < &
            (norm2(mesh1%elem(iclosest)%coord - mesh2%elem(i2)%coord))) then
            iclosest = itest
          end if
        end do

        if (i1 == iclosest) then
          exit
        else
          i1 = iclosest
        end if
      end do

      !Interpolate solution
      sol2(:, i2) = sol1(:, i1)
    end do
  end subroutine project_sol

  subroutine find_if_vert_is_bound(mesh, b2d)
    implicit none

    type(mesh_type), intent(inout) :: mesh
    logical, intent(in) :: b2d

    integer(kind=DOUBLE) :: i, j, id_face

    do i = 1, mesh%n_vert
      mesh%vert(i)%is_bound = .FALSE.
      do j = 1, mesh%vert(i)%n_faces_neigh
        id_face = mesh%vert(i)%face_neigh(j)
        if( b2d ) then
          if (mesh%face(id_face)%right_neigh <= 0 .and. &
            abs(mesh%face(id_face)%norm(3)) < 1e-12_DOUBLE ) then
            mesh%vert(i)%is_bound = .TRUE.
          end if
        else
          if (mesh%face(id_face)%right_neigh <= 0) then
            mesh%vert(i)%is_bound = .TRUE.
          end if
        end if
      end do
    end do
  end subroutine find_if_vert_is_bound

  subroutine mpi_project_sol_box(sol_size, mesh1, sol1, mesh2, sol2)
    use mpi
    implicit none

    integer(kind=ENTIER) :: sol_size
    type(mesh_type), intent(in) :: mesh1, mesh2
    real(kind=DOUBLE), dimension(sol_size, mesh1%n_elems), intent(in) :: sol1
    real(kind=DOUBLE), dimension(sol_size, mesh2%n_elems), intent(inout) :: sol2

    type :: box_type
      integer(kind=ENTIER), dimension(:), allocatable :: elem
      real(kind=DOUBLE), dimension(:,:), allocatable :: coord
      real(kind=DOUBLE), dimension(:,:), allocatable :: sol
    end type box_type

    real(kind=DOUBLE), parameter :: r_per_box = 5.0_DOUBLE

    integer(kind=ENTIER) :: me, num_procs, p, mpi_ierr
    integer(kind=ENTIER) :: i, nx, ny, nz, ix, iy, iz, iclosest, j
    integer(kind=ENTIER) :: kx, ky, kz, idbx, idby, idbz
    real(kind=DOUBLE) :: xmin, xmax, ymin, ymax, zmin, zmax
    real(kind=DOUBLE) :: r, tot_vol, dx, dy, dz, d1
    real(kind=DOUBLE) :: dx_current, dy_current, dz_current
    integer(kind=ENTIER), dimension(:, :, :), allocatable :: n_box
    type(box_type), dimension(:, :, :), allocatable :: box
    real(kind=double), dimension(mesh2%n_elems) :: d_closest

    real(kind=DOUBLE) :: t1, t2
    real(kind=DOUBLE) :: xmin_current, xmax_current, &
      ymin_current, ymax_current, &
      zmin_current, zmax_current
    integer(kind=ENTIER) :: nx_current, ny_current, nz_current
    integer(kind=ENTIER), dimension(:, :, :), allocatable :: n_box_current
    type(box_type), dimension(:, :, :), allocatable :: box_current

    real(kind=double), dimension(3, 5, mesh2%n_elems) :: grad
    real(kind=double), dimension(mesh2%n_elems) :: residu
    integer(kind=ENTIER), dimension(mesh2%n_elems) :: color
    character(len=255) :: me_str

    logical :: found

    sol2 = -1.
    grad = 0.
    residu = 0.
    color = 0

    call cpu_time(t1)

    xmin = mesh1%elem(1)%coord(1)
    xmax = mesh1%elem(1)%coord(1)
    ymin = mesh1%elem(1)%coord(2)
    ymax = mesh1%elem(1)%coord(2)
    zmin = mesh1%elem(1)%coord(3)
    zmax = mesh1%elem(1)%coord(3)

    tot_vol = 0.0_DOUBLE
    do i = 1, mesh1%n_elems
      if (mesh1%elem(i)%coord(1) < xmin) xmin = mesh1%elem(i)%coord(1)
      if (mesh1%elem(i)%coord(1) > xmax) xmax = mesh1%elem(i)%coord(1)
      if (mesh1%elem(i)%coord(2) < ymin) ymin = mesh1%elem(i)%coord(2)
      if (mesh1%elem(i)%coord(2) > ymax) ymax = mesh1%elem(i)%coord(2)
      if (mesh1%elem(i)%coord(3) < zmin) zmin = mesh1%elem(i)%coord(3)
      if (mesh1%elem(i)%coord(3) > zmax) zmax = mesh1%elem(i)%coord(3)
      tot_vol = tot_vol + mesh1%elem(i)%volume
    end do

    r = (tot_vol/mesh1%n_elems)**(1.0_DOUBLE/3.0_DOUBLE)

    xmin = xmin - r
    xmax = xmax + r
    ymin = ymin - r
    ymax = ymax + r
    zmin = zmin - r
    zmax = zmax + r

    if (xmax - xmin < r) then
      dx = xmax - xmin
    else
      dx = r_per_box*r
    end if
    nx = ceiling((xmax - xmin)/dx)+1

    if (ymax - ymin < r) then
      dy = ymax - ymin
    else
      dy = r_per_box*r
    end if
    ny = ceiling((ymax - ymin)/dy)+1

    if (zmax - zmin < r) then
      dz = zmax - zmin
    else
      dz = r_per_box*r
    end if
    nz = ceiling((zmax - zmin)/dz)+1

    allocate (n_box(nx, ny, nz))
    n_box = 0
    do i = 1, mesh1%n_elems
      ix = ceiling((mesh1%elem(i)%coord(1) - xmin)/dx)
      iy = ceiling((mesh1%elem(i)%coord(2) - ymin)/dy)
      iz = ceiling((mesh1%elem(i)%coord(3) - zmin)/dz)
      n_box(ix, iy, iz) = n_box(ix, iy, iz) + 1
    end do

    allocate (box(nx, ny, nz))
    do ix = 1, nx
      do iy = 1, ny
        do iz = 1, nz
          if (n_box(ix, iy, iz) > 0) then
            allocate (box(ix, iy, iz)%elem(n_box(ix, iy, iz)))
            box(ix, iy, iz)%elem = 0
            allocate (box(ix, iy, iz)%sol(sol_size, n_box(ix, iy, iz)))
            box(ix, iy, iz)%sol = 0
            allocate (box(ix, iy, iz)%coord(3, n_box(ix, iy, iz)))
            box(ix, iy, iz)%coord = 0
          end if
        end do
      end do
    end do

    n_box = 0
    do i = 1, mesh1%n_elems
      ix = ceiling((mesh1%elem(i)%coord(1) - xmin)/dx)
      iy = ceiling((mesh1%elem(i)%coord(2) - ymin)/dy)
      iz = ceiling((mesh1%elem(i)%coord(3) - zmin)/dz)
      n_box(ix, iy, iz) = n_box(ix, iy, iz) + 1
      box(ix, iy, iz)%elem(n_box(ix, iy, iz)) = i
      box(ix, iy, iz)%sol(:, n_box(ix, iy, iz)) = sol1(:, i)
      box(ix, iy, iz)%coord(:, n_box(ix, iy, iz)) = mesh1%elem(i)%coord
    end do

    !MPI
    sol2 = 0.0_DOUBLE
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    d_closest = 1e10

    do p=0, num_procs-1
      if( p == me ) then
        print*, "Proj from", me
        xmin_current = xmin
        xmax_current = xmax
        ymin_current = ymin
        ymax_current = ymax
        zmin_current = zmin
        zmax_current = zmax
        nx_current = nx
        ny_current = ny
        nz_current = nz
        dx_current = dx
        dy_current = dy
        dz_current = dz
      end if

      call mpi_bcast(xmin_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(xmax_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(ymin_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(ymax_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(zmin_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(zmax_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(dx_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(dy_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(dz_current, 1, MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(nx_current, 1, MPI_INT, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(ny_current, 1, MPI_INT, p, MPI_COMM_WORLD, mpi_ierr)
      call mpi_bcast(nz_current, 1, MPI_INT, p, MPI_COMM_WORLD, mpi_ierr)
      allocate (n_box_current(nx_current, ny_current, nz_current))
      if( me == p ) n_box_current = n_box
      call mpi_bcast(n_box_current, nx_current*ny_current*nz_current, &
        MPI_INT, p, MPI_COMM_WORLD, mpi_ierr)

      allocate (box_current(nx_current, ny_current, nz_current))
      do ix = 1, nx_current
        do iy = 1, ny_current
          do iz = 1, nz_current
            if (n_box_current(ix, iy, iz) > 0) then
              allocate (box_current(ix, iy, iz)%elem(n_box_current(ix, iy, iz)))
              if( me == p ) box_current(ix, iy, iz)%elem = box(ix, iy, iz)%elem
              call mpi_bcast(box_current(ix,iy,iz)%elem, n_box_current(ix,iy,iz), &
                MPI_INT, p, MPI_COMM_WORLD, mpi_ierr)
              allocate (box_current(ix, iy, iz)%sol(sol_size, n_box_current(ix, iy, iz)))
              if( me == p ) box_current(ix, iy, iz)%sol = box(ix, iy, iz)%sol
              call mpi_bcast(box_current(ix,iy,iz)%sol, sol_size*n_box_current(ix,iy,iz), &
                MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
              allocate (box_current(ix, iy, iz)%coord(3, n_box_current(ix, iy, iz)))
              if( me == p ) box_current(ix, iy, iz)%coord = box(ix, iy, iz)%coord
              call mpi_bcast(box_current(ix,iy,iz)%coord, 3*n_box_current(ix,iy,iz), &
                MPI_DOUBLE, p, MPI_COMM_WORLD, mpi_ierr)
            end if
          end do
        end do
      end do

      do i = 1, mesh2%n_elems
        if( mesh2%elem(i)%coord(1) > xmin_current .and. &
          mesh2%elem(i)%coord(1) < xmax_current .and.&
          mesh2%elem(i)%coord(2) > ymin_current .and. &
          mesh2%elem(i)%coord(2) < ymax_current .and.&
          mesh2%elem(i)%coord(3) > zmin_current .and. &
          mesh2%elem(i)%coord(3) < zmax_current ) then

          ix = ceiling((mesh2%elem(i)%coord(1) - xmin_current)/dx_current)
          iy = ceiling((mesh2%elem(i)%coord(2) - ymin_current)/dy_current)
          iz = ceiling((mesh2%elem(i)%coord(3) - zmin_current)/dz_current)

          found = .false.
          iclosest = 0
          do kx = 1, 5
            do ky = 1, 5
              do kz = 1, 5
                idbx = min(max(1, ix + kx - 3), nx_current)
                idby = min(max(1, iy + ky - 3), ny_current)
                idbz = min(max(1, iz + kz - 3), nz_current)
                do j = 1, n_box_current(idbx, idby, idbz)
                  d1 = norm2(mesh2%elem(i)%coord &
                    - box_current(idbx, idby, idbz)%coord(:, j))
                  if (d1 < d_closest(i)) then
                    iclosest = box_current(idbx, idby, idbz)%elem(j)
                    sol2(:, i) = box_current(idbx, idby, idbz)%sol(:, j)
                    d_closest(i) = d1
                    found = .true.
                  end if
                end do
              end do
            end do
          end do

        end if !bounding box
      end do
      deallocate (n_box_current, box_current)
      call mpi_barrier(mpi_comm_world, mpi_ierr)
    end do

    call mpi_barrier(mpi_comm_world, mpi_ierr)
    call cpu_time(t2)
    if( me == 0 ) then
      print *, ""//achar(27)//"[33m[*] Time for projection onto next mesh :"//achar(27)//"[0m", t2 - t1
    end if
  end subroutine mpi_project_sol_box
end module mesh_geometry_module
